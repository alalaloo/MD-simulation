#pragma once
#include <omp.h>
#include <vector>
#include <cmath>
#include <iostream>
#include <random>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <chrono>

#ifdef USE_CUDA
#include <cuda_runtime.h>
#include <cuda.h>
#include <device_launch_parameters.h>
#endif

// ==================== CUDA КЕРНЕЛЫ (ДЛЯ GPU) ====================
#ifdef USE_CUDA
__global__ void build_verlet_list_gpu_kernel(double* rx, double* ry, double* rz, int N,
                                             int* head, int* linked_list,
                                             int M_cells_x, int M_cells_y, int M_cells_z,
                                             double cell_size_x, double cell_size_y, double cell_size_z,
                                             double rlist_sq, double L, double inv_L,
                                             int* d_neighbor_list, int* d_neighbor_counts,
                                             int max_neighbors) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;

    int cx = static_cast<int>(rx[i] / cell_size_x);
    int cy = static_cast<int>(ry[i] / cell_size_y);
    int cz = static_cast<int>(rz[i] / cell_size_z);

    if (cx >= M_cells_x) cx = M_cells_x - 1; else if (cx < 0) cx = 0;
    if (cy >= M_cells_y) cy = M_cells_y - 1; else if (cy < 0) cy = 0;
    if (cz >= M_cells_z) cz = M_cells_z - 1; else if (cz < 0) cz = 0;

    int count = 0;

    for (int dx_c = -1; dx_c <= 1; ++dx_c) {
        int scx = (cx + dx_c + M_cells_x) % M_cells_x;
        for (int dy_c = -1; dy_c <= 1; ++dy_c) {
            int scy = (cy + dy_c + M_cells_y) % M_cells_y;
            for (int dz_c = -1; dz_c <= 1; ++dz_c) {
                int scz = (cz + dz_c + M_cells_z) % M_cells_z;
                int s_idx = scx + scy * M_cells_x + scz * M_cells_x * M_cells_y;

                int j = head[s_idx];
                while (j != -1) {
                    if (i != j) {
                        double dx = rx[i] - rx[j];
                        double dy = ry[i] - ry[j];
                        double dz = rz[i] - rz[j];
                        
                        dx -= L * floor(dx * inv_L + 0.5);
                        dy -= L * floor(dy * inv_L + 0.5);
                        dz -= L * floor(dz * inv_L + 0.5);

                        if (dx*dx + dy*dy + dz*dz < rlist_sq) {
                            if (count < max_neighbors) {
                                d_neighbor_list[i * max_neighbors + count] = j;
                                count++;
                            }
                        }
                    }
                    j = linked_list[j];
                }
            }
        }
    }
    d_neighbor_counts[i] = count;
}

__global__ void compute_forces_gpu_kernel(double* rx, double* ry, double* rz,
                                          int* d_neighbor_list, int* d_neighbor_counts, int N,
                                          double cutoff_sq, double L, double inv_L,
                                          double* fx, double* fy, double* fz,
                                          double* d_U_total, int max_neighbors,
                                          double sigma, double epsilon) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;

    double fx_local = 0.0, fy_local = 0.0, fz_local = 0.0;
    double u_local = 0.0;
    double sigma2 = sigma * sigma;

    int count = d_neighbor_counts[i];
    for (int k = 0; k < count; ++k) {
        int j = d_neighbor_list[i * max_neighbors + k];
        double dx = rx[i] - rx[j];
        double dy = ry[i] - ry[j];
        double dz = rz[i] - rz[j];
        
        dx -= L * floor(dx * inv_L + 0.5);
        dy -= L * floor(dy * inv_L + 0.5);
        dz -= L * floor(dz * inv_L + 0.5);

        double r2 = dx*dx + dy*dy + dz*dz;
        if (r2 < cutoff_sq) {
            double s2 = sigma2 / r2;
            double s6 = s2 * s2 * s2;
            double force_mag = 24.0 * epsilon * (2.0 * s6 * s6 - s6) / r2;
            fx_local += force_mag * dx;
            fy_local += force_mag * dy;
            fz_local += force_mag * dz;
            u_local += 4.0 * epsilon * (s6 * s6 - s6);
        }
    }

    fx[i] = fx_local;
    fy[i] = fy_local;
    fz[i] = fz_local;

    atomicAdd(d_U_total, u_local * 0.5);
}

__global__ void verlet_step_part1_kernel(double* rx, double* ry, double* rz,
                                         double* vx, double* vy, double* vz,
                                         double* fx, double* fy, double* fz,
                                         int N, double dt, double mass, double L, double inv_L) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;
    
    vx[i] += 0.5 * dt * fx[i] / mass;
    vy[i] += 0.5 * dt * fy[i] / mass;
    vz[i] += 0.5 * dt * fz[i] / mass;
    
    rx[i] += dt * vx[i];
    ry[i] += dt * vy[i];
    rz[i] += dt * vz[i];
    
    if (rx[i] < 0.0) rx[i] += L; else if (rx[i] >= L) rx[i] -= L;
    if (ry[i] < 0.0) ry[i] += L; else if (ry[i] >= L) ry[i] -= L;
    if (rz[i] < 0.0) rz[i] += L; else if (rz[i] >= L) rz[i] -= L;
}

__global__ void verlet_step_part2_kernel(double* vx, double* vy, double* vz,
                                         double* fx, double* fy, double* fz,
                                         int N, double dt, double mass) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;
    
    vx[i] += 0.5 * dt * fx[i] / mass;
    vy[i] += 0.5 * dt * fy[i] / mass;
    vz[i] += 0.5 * dt * fz[i] / mass;
}

__global__ void scale_velocities_kernel(double* vx, double* vy, double* vz, int N, double s) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < N) {
        vx[i] *= s; vy[i] *= s; vz[i] *= s;
    }
}

__global__ void compute_kinetic_energy_kernel(double* vx, double* vy, double* vz, int N, double mass, double* d_K) {
    extern __shared__ double sdata[];
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int tid = threadIdx.x;

    double local_K = 0.0;
    if (i < N) {
        local_K = 0.5 * mass * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
    }
    sdata[tid] = local_K;
    __syncthreads();

    for (int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (tid < s) sdata[tid] += sdata[tid + s];
        __syncthreads();
    }
    if (tid == 0) atomicAdd(d_K, sdata[0]);
}

__global__ void check_displacement_kernel(double* rx, double* ry, double* rz,
                                          double* rx0, double* ry0, double* rz0,
                                          int N, double L, double inv_L, double* d_block_maxes) {
    extern __shared__ double sdata[];
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int tid = threadIdx.x;

    double max_d2 = 0.0;
    if (i < N) {
        double dx = rx[i] - rx0[i];
        double dy = ry[i] - ry0[i];
        double dz = rz[i] - rz0[i];
        dx -= L * floor(dx * inv_L + 0.5);
        dy -= L * floor(dy * inv_L + 0.5);
        dz -= L * floor(dz * inv_L + 0.5);
        max_d2 = dx*dx + dy*dy + dz*dz;
    }
    sdata[tid] = max_d2;
    __syncthreads();

    for (int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (tid < s) {
            if (sdata[tid + s] > sdata[tid]) sdata[tid] = sdata[tid + s];
        }
        __syncthreads();
    }
    if (tid == 0) d_block_maxes[blockIdx.x] = sdata[0];
}
#endif

// ==================== КЛАСС SYSTEMMD ====================
class systemMD {
private:
    std::vector<std::vector<double>> fx_priv, fy_priv, fz_priv;
    std::vector<double> u_priv;

public:
    int N;
    double L;
    double inv_L; 
    std::vector<double> rx, ry, rz;
    std::vector<double> vx, vy, vz;
    std::vector<double> fx, fy, fz;
    double mass;
    double sigma;      // Параметры потенциала Леннард-Джонса в СИ
    double epsilon;    // Параметры потенциала Леннард-Джонса в СИ
    double cutoff;
    double cutoff_sq;
    double skin;
    double rlist;
    double rlist_sq;
    double potential_energy;

    bool thermostat_enabled;
    double xi1, xi2;
    double v_xi1, v_xi2;
    double Q1, Q2;
    double T_target;
    double kB;

    struct VerletList {
        std::vector<int> i_list, j_list;
        std::vector<double> rx0, ry0, rz0;
        double max_disp;
        bool need_rebuild;
#ifdef USE_CUDA
        double *d_rx0 = nullptr, *d_ry0 = nullptr, *d_rz0 = nullptr;
#endif
    } verlet;

    int M_cells_x, M_cells_y, M_cells_z;
    int num_cells;
    double cell_size_x, cell_size_y, cell_size_z;
    std::vector<int> head;
    std::vector<int> linked_list;

    long long time_verlet_rebuild_ns = 0;
    long long time_check_disp_ns = 0;
    long long time_force_calc_ns = 0;
    long long time_io_ns = 0;
    int count_rebuilds = 0;

#ifdef USE_CUDA
    double *d_rx = nullptr, *d_ry = nullptr, *d_rz = nullptr;
    double *d_vx = nullptr, *d_vy = nullptr, *d_vz = nullptr;
    double *d_fx = nullptr, *d_fy = nullptr, *d_fz = nullptr;
    int *d_head = nullptr, *d_linked_list = nullptr;
    int *d_neighbor_list = nullptr, *d_neighbor_counts = nullptr;
    double *d_K = nullptr, *d_U_total = nullptr, *d_block_maxes = nullptr;
    int cuda_blocks_count = 0;
    bool cuda_allocated = false;
#endif

    systemMD(int n_req, double a0, double lj_sigma, double lj_epsilon, double particle_mass,
             double cut_off, double skin_depth)
        : mass(particle_mass), sigma(lj_sigma), epsilon(lj_epsilon), cutoff(cut_off), skin(skin_depth)
    {
        cutoff_sq = cutoff * cutoff;
        rlist = cutoff + skin;
        rlist_sq = rlist * rlist;

        int n_cells = std::max(1, (int)std::round(std::cbrt(n_req / 4.0)));
        N = 4 * n_cells * n_cells * n_cells;
        L = n_cells * a0;
        inv_L = 1.0 / L;

        rx.resize(N); ry.resize(N); rz.resize(N);
        vx.resize(N); vy.resize(N); vz.resize(N);
        fx.resize(N); fy.resize(N); fz.resize(N);

        int idx = 0;
        for (int ix = 0; ix < n_cells && idx < N; ++ix)
            for (int iy = 0; iy < n_cells && idx < N; ++iy)
                for (int iz = 0; iz < n_cells && idx < N; ++iz) {
                    double x = ix * a0, y = iy * a0, z = iz * a0;
                    if (idx < N) { rx[idx]=x;         ry[idx]=y;         rz[idx]=z;         ++idx; }
                    if (idx < N) { rx[idx]=x+0.5*a0;  ry[idx]=y+0.5*a0;  rz[idx]=z;         ++idx; }
                    if (idx < N) { rx[idx]=x+0.5*a0;  ry[idx]=y;         rz[idx]=z+0.5*a0;  ++idx; }
                    if (idx < N) { rx[idx]=x;         ry[idx]=y+0.5*a0;  rz[idx]=z+0.5*a0;  ++idx; }
                }

        thermostat_enabled = false;
        xi1 = xi2 = v_xi1 = v_xi2 = 0.0;
        T_target = 1.0; kB = 1.380649e-23; // Инициализация СИ-значением
        potential_energy = 0.0;

        verlet.max_disp = 0.0;
        verlet.need_rebuild = true;

        M_cells_x = std::max(1, static_cast<int>(L / rlist));
        M_cells_y = std::max(1, static_cast<int>(L / rlist));
        M_cells_z = std::max(1, static_cast<int>(L / rlist));
        num_cells = M_cells_x * M_cells_y * M_cells_z;

        cell_size_x = L / M_cells_x;
        cell_size_y = L / M_cells_y;
        cell_size_z = L / M_cells_z;

        head.resize(num_cells);
        linked_list.resize(N);
    }

    ~systemMD() {
#ifdef USE_CUDA
        cuda_free();
#endif
    }

    int getN() const { return N; }

    void set_velocities(double T, double k_boltzmann = 1.380649e-23) {
        kB = k_boltzmann;
        T_target = T;
        std::mt19937 gen(42);
        std::normal_distribution<double> dist(0.0, std::sqrt(kB * T / mass));
        double sum_vx = 0, sum_vy = 0, sum_vz = 0;
        for (int i = 0; i < N; ++i) {
            vx[i] = dist(gen); vy[i] = dist(gen); vz[i] = dist(gen);
            sum_vx += vx[i]; sum_vy += vy[i]; sum_vz += vz[i];
        }
        double cm_vx = sum_vx / N, cm_vy = sum_vy / N, cm_vz = sum_vz / N;
        double K = 0.0;
        for (int i = 0; i < N; ++i) {
            vx[i] -= cm_vx; vy[i] -= cm_vy; vz[i] -= cm_vz;
            K += 0.5 * mass * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
        }
        double T_cur = (2.0 / 3.0) * K / (N * kB);
        double scale = std::sqrt(T / T_cur);
        for (int i = 0; i < N; ++i) {
            vx[i] *= scale; vy[i] *= scale; vz[i] *= scale;
        }
#ifdef USE_CUDA
        if (cuda_allocated) sync_host_to_device();
#endif
    }

    void enable_nose_hoover(double T, double tau_param) {
        thermostat_enabled = true;
        T_target = T;
        double g = 3.0 * N;
        double tau = (tau_param > 0.0) ? tau_param : 1.0e-13;
        
        Q1 = g * kB * T_target * tau * tau;
        Q2 = kB * T_target * tau * tau;
        
        xi1 = xi2 = v_xi1 = v_xi2 = 0.0;
    }

    void disable_thermostat() { thermostat_enabled = false; }

    inline void apply_pbc(double& dx, double& dy, double& dz) const {
        dx -= L * std::floor(dx * inv_L + 0.5);
        dy -= L * std::floor(dy * inv_L + 0.5);
        dz -= L * std::floor(dz * inv_L + 0.5);
    }

    void build_verlet_list_cpu() {
        verlet.i_list.clear(); verlet.j_list.clear();
        verlet.rx0 = rx; verlet.ry0 = ry; verlet.rz0 = rz;
        verlet.max_disp = 0.0;
        verlet.need_rebuild = false;

        std::fill(head.begin(), head.end(), -1);
        for (int i = 0; i < N; ++i) {
            int cx = static_cast<int>(rx[i] / cell_size_x);
            int cy = static_cast<int>(ry[i] / cell_size_y);
            int cz = static_cast<int>(rz[i] / cell_size_z);
            if (cx >= M_cells_x) cx = M_cells_x - 1; else if (cx < 0) cx = 0;
            if (cy >= M_cells_y) cy = M_cells_y - 1; else if (cy < 0) cy = 0;
            if (cz >= M_cells_z) cz = M_cells_z - 1; else if (cz < 0) cz = 0;

            int cell_idx = cx + cy * M_cells_x + cz * M_cells_x * M_cells_y;
            linked_list[i] = head[cell_idx];
            head[cell_idx] = i;
        }

        int num_threads = omp_get_max_threads();
        std::vector<std::vector<int>> local_i(num_threads);
        std::vector<std::vector<int>> local_j(num_threads);

        #pragma omp parallel
        {
            int tid = omp_get_thread_num();
            #pragma omp for collapse(3) schedule(dynamic)
            for (int cx = 0; cx < M_cells_x; ++cx) {
                for (int cy = 0; cy < M_cells_y; ++cy) {
                    for (int cz = 0; cz < M_cells_z; ++cz) {
                        int c_idx = cx + cy * M_cells_x + cz * M_cells_x * M_cells_y;

                        for (int dx_c = -1; dx_c <= 1; ++dx_c) {
                            int scx = (cx + dx_c + M_cells_x) % M_cells_x;
                            for (int dy_c = -1; dy_c <= 1; ++dy_c) {
                                int scy = (cy + dy_c + M_cells_y) % M_cells_y;
                                for (int dz_c = -1; dz_c <= 1; ++dz_c) {
                                    int scz = (cz + dz_c + M_cells_z) % M_cells_z;
                                    int s_idx = scx + scy * M_cells_x + scz * M_cells_x * M_cells_y;

                                    if (c_idx > s_idx) continue;

                                    int i = head[c_idx];
                                    while (i != -1) {
                                        int j = head[s_idx];
                                        while (j != -1) {
                                            if (c_idx == s_idx && i >= j) {
                                                j = linked_list[j];
                                                continue;
                                            }
                                            double dx = rx[i] - rx[j];
                                            double dy = ry[i] - ry[j];
                                            double dz = rz[i] - rz[j];
                                            apply_pbc(dx, dy, dz);

                                            if (dx*dx + dy*dy + dz*dz < rlist_sq) {
                                                local_i[tid].push_back(i);
                                                local_j[tid].push_back(j);
                                            }
                                            j = linked_list[j];
                                        }
                                        i = linked_list[i];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        for (int t = 0; t < num_threads; ++t) {
            verlet.i_list.insert(verlet.i_list.end(), local_i[t].begin(), local_i[t].end());
            verlet.j_list.insert(verlet.j_list.end(), local_j[t].begin(), local_j[t].end());
        }
    }

#ifdef USE_CUDA
    void cuda_allocate() {
        if (cuda_allocated) return;
        cudaMalloc(&d_rx, N * sizeof(double));
        cudaMalloc(&d_ry, N * sizeof(double));
        cudaMalloc(&d_rz, N * sizeof(double));
        cudaMalloc(&d_vx, N * sizeof(double));
        cudaMalloc(&d_vy, N * sizeof(double));
        cudaMalloc(&d_vz, N * sizeof(double));
        cudaMalloc(&d_fx, N * sizeof(double));
        cudaMalloc(&d_fy, N * sizeof(double));
        cudaMalloc(&d_fz, N * sizeof(double));

        cudaMalloc(&d_head, num_cells * sizeof(int));
        cudaMalloc(&d_linked_list, N * sizeof(int));
        
        int max_neighbors = 256;
        cudaMalloc(&d_neighbor_list, N * max_neighbors * sizeof(int));
        cudaMalloc(&d_neighbor_counts, N * sizeof(int));

        cudaMalloc(&d_K, sizeof(double));
        cudaMalloc(&d_U_total, sizeof(double));

        cudaMalloc(&verlet.d_rx0, N * sizeof(double));
        cudaMalloc(&verlet.d_ry0, N * sizeof(double));
        cudaMalloc(&verlet.d_rz0, N * sizeof(double));

        int threads = 256;
        cuda_blocks_count = (N + threads - 1) / threads;
        cudaMalloc(&d_block_maxes, cuda_blocks_count * sizeof(double));

        cuda_allocated = true;
        sync_host_to_device();
    }

    void cuda_free() {
        if (!cuda_allocated) return;
        cudaFree(d_rx); cudaFree(d_ry); cudaFree(d_rz);
        cudaFree(d_vx); cudaFree(d_vy); cudaFree(d_vz);
        cudaFree(d_fx); cudaFree(d_fy); cudaFree(d_fz);
        cudaFree(d_head); cudaFree(d_linked_list);
        cudaFree(d_neighbor_list); cudaFree(d_neighbor_counts);
        cudaFree(d_K); cudaFree(d_U_total); cudaFree(d_block_maxes);
        cudaFree(verlet.d_rx0); cudaFree(verlet.d_ry0); cudaFree(verlet.d_rz0);
        cuda_allocated = false;
    }

    void sync_host_to_device() {
        cudaMemcpy(d_rx, rx.data(), N * sizeof(double), cudaMemcpyHostToDevice);
        cudaMemcpy(d_ry, ry.data(), N * sizeof(double), cudaMemcpyHostToDevice);
        cudaMemcpy(d_rz, rz.data(), N * sizeof(double), cudaMemcpyHostToDevice);
        cudaMemcpy(d_vx, vx.data(), N * sizeof(double), cudaMemcpyHostToDevice);
        cudaMemcpy(d_vy, vy.data(), N * sizeof(double), cudaMemcpyHostToDevice);
        cudaMemcpy(d_vz, vz.data(), N * sizeof(double), cudaMemcpyHostToDevice);
    }

    void sync_device_to_host() {
        cudaMemcpy(rx.data(), d_rx, N * sizeof(double), cudaMemcpyDeviceToHost);
        cudaMemcpy(ry.data(), d_ry, N * sizeof(double), cudaMemcpyDeviceToHost);
        cudaMemcpy(rz.data(), d_rz, N * sizeof(double), cudaMemcpyDeviceToHost);
        cudaMemcpy(vx.data(), d_vx, N * sizeof(double), cudaMemcpyDeviceToHost);
        cudaMemcpy(vy.data(), d_vy, N * sizeof(double), cudaMemcpyDeviceToHost);
        cudaMemcpy(vz.data(), d_vz, N * sizeof(double), cudaMemcpyDeviceToHost);
        cudaMemcpy(fx.data(), d_fx, N * sizeof(double), cudaMemcpyDeviceToHost);
        cudaMemcpy(fy.data(), d_fy, d_fy, N * sizeof(double), cudaMemcpyDeviceToHost); // typo protection: fy.data() from d_fy
        cudaMemcpy(fy.data(), d_fy, N * sizeof(double), cudaMemcpyDeviceToHost);
        cudaMemcpy(fz.data(), d_fz, N * sizeof(double), cudaMemcpyDeviceToHost);
    }

    void build_verlet_list_gpu() {
        cuda_allocate();
        cudaMemcpy(rx.data(), d_rx, N * sizeof(double), cudaMemcpyDeviceToHost);
        cudaMemcpy(ry.data(), d_ry, N * sizeof(double), cudaMemcpyDeviceToHost);
        cudaMemcpy(rz.data(), d_rz, N * sizeof(double), cudaMemcpyDeviceToHost);

        std::fill(head.begin(), head.end(), -1);
        for (int i = 0; i < N; ++i) {
            int cx = static_cast<int>(rx[i] / cell_size_x);
            int cy = static_cast<int>(ry[i] / cell_size_y);
            int cz = static_cast<int>(rz[i] / cell_size_z);
            if (cx >= M_cells_x) cx = M_cells_x - 1; else if (cx < 0) cx = 0;
            if (cy >= M_cells_y) cy = M_cells_y - 1; else if (cy < 0) cy = 0;
            if (cz >= M_cells_z) cz = M_cells_z - 1; else if (cz < 0) cz = 0;
            int cell_idx = cx + cy * M_cells_x + cz * M_cells_x * M_cells_y;
            linked_list[i] = head[cell_idx];
            head[cell_idx] = i;
        }

        cudaMemcpy(d_head, head.data(), num_cells * sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(d_linked_list, linked_list.data(), N * sizeof(int), cudaMemcpyHostToDevice);

        int threads = 256;
        int max_neighbors = 256;
        build_verlet_list_gpu_kernel<<<cuda_blocks_count, threads>>>(d_rx, d_ry, d_rz, N,
                                                                     d_head, d_linked_list,
                                                                     M_cells_x, M_cells_y, M_cells_z,
                                                                     cell_size_x, cell_size_y, cell_size_z,
                                                                     rlist_sq, L, inv_L,
                                                                     d_neighbor_list, d_neighbor_counts,
                                                                     max_neighbors);
        cudaDeviceSynchronize();

        cudaMemcpy(verlet.d_rx0, d_rx, N * sizeof(double), cudaMemcpyDeviceToDevice);
        cudaMemcpy(verlet.d_ry0, d_ry, N * sizeof(double), cudaMemcpyDeviceToDevice);
        cudaMemcpy(verlet.d_rz0, d_rz, N * sizeof(double), cudaMemcpyDeviceToDevice);
        
        verlet.rx0 = rx; verlet.ry0 = ry; verlet.rz0 = rz;
        verlet.need_rebuild = false;
    }
#endif

    void build_verlet_list() {
        auto start = std::chrono::high_resolution_clock::now();
#ifdef USE_CUDA
        build_verlet_list_gpu();
#else
        build_verlet_list_cpu();
#endif
        auto end = std::chrono::high_resolution_clock::now();
        time_verlet_rebuild_ns += std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
        count_rebuilds++;
    }

    void check_displacement() {
        auto start = std::chrono::high_resolution_clock::now();
        double max_disp = 0.0;
#ifdef USE_CUDA
        int threads = 256;
        size_t shared_mem = threads * sizeof(double);
        check_displacement_kernel<<<cuda_blocks_count, threads, shared_mem>>>(d_rx, d_ry, d_rz, 
                                                                             verlet.d_rx0, verlet.d_ry0, verlet.d_rz0, 
                                                                             N, L, inv_L, d_block_maxes);
        cudaDeviceSynchronize();
        
        std::vector<double> h_block_maxes(cuda_blocks_count);
        cudaMemcpy(h_block_maxes.data(), d_block_maxes, cuda_blocks_count * sizeof(double), cudaMemcpyDeviceToHost);
        double max2 = 0.0;
        for (int b = 0; b < cuda_blocks_count; ++b) {
            if (h_block_maxes[b] > max2) max2 = h_block_maxes[b];
        }
        max_disp = std::sqrt(max2);
#else
        double max2 = 0.0;
        #pragma omp parallel for reduction(max:max2)
        for (int i = 0; i < N; ++i) {
            double dx = rx[i] - verlet.rx0[i];
            double dy = ry[i] - verlet.ry0[i];
            double dz = rz[i] - verlet.rz0[i];
            apply_pbc(dx, dy, dz);
            double d2 = dx*dx + dy*dy + dz*dz;
            if (d2 > max2) max2 = d2;
        }
        max_disp = std::sqrt(max2);
#endif
        verlet.max_disp = max_disp;
        if (verlet.max_disp > 0.5 * skin) verlet.need_rebuild = true;
        auto end = std::chrono::high_resolution_clock::now();
        time_check_disp_ns += std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    }

    void compute_forces() {
        if (verlet.need_rebuild) {
            build_verlet_list();
        } else {
            check_displacement();
            if (verlet.need_rebuild) build_verlet_list();
        }

        auto start_force = std::chrono::high_resolution_clock::now();

#ifdef USE_CUDA
        double zero_u = 0.0;
        cudaMemcpy(d_U_total, &zero_u, sizeof(double), cudaMemcpyHostToDevice);
        int threads = 256;
        int max_neighbors = 256;
        compute_forces_gpu_kernel<<<cuda_blocks_count, threads>>>(d_rx, d_ry, d_rz,
                                                                 d_neighbor_list, d_neighbor_counts, N,
                                                                 cutoff_sq, L, inv_L,
                                                                 d_fx, d_fy, d_fz, d_U_total, max_neighbors,
                                                                 sigma, epsilon);
        cudaDeviceSynchronize();
        cudaMemcpy(&potential_energy, d_U_total, sizeof(double), cudaMemcpyDeviceToHost);
#else
        int num_threads = omp_get_max_threads();
        if (fx_priv.size() != (size_t)num_threads || fx_priv[0].size() != (size_t)N) {
            fx_priv.assign(num_threads, std::vector<double>(N, 0.0));
            fy_priv.assign(num_threads, std::vector<double>(N, 0.0));
            fz_priv.assign(num_threads, std::vector<double>(N, 0.0));
            u_priv.assign(num_threads, 0.0);
        } else {
            #pragma omp parallel for
            for (int t = 0; t < num_threads; ++t) {
                std::fill(fx_priv[t].begin(), fx_priv[t].end(), 0.0);
                std::fill(fy_priv[t].begin(), fy_priv[t].end(), 0.0);
                std::fill(fz_priv[t].begin(), fz_priv[t].end(), 0.0);
                u_priv[t] = 0.0;
            }
        }

        std::fill(fx.begin(), fx.end(), 0.0);
        std::fill(fy.begin(), fy.end(), 0.0);
        std::fill(fz.begin(), fz.end(), 0.0);
        potential_energy = 0.0;

        size_t num_pairs = verlet.i_list.size();
        double sigma2 = sigma * sigma;

        #pragma omp parallel
        {
            int tid = omp_get_thread_num();
            double u_local = 0.0;
            #pragma omp Epic for structure
            #pragma omp for schedule(static)
            for (size_t k = 0; k < num_pairs; ++k) {
                int i = verlet.i_list[k];
                int j = verlet.j_list[k];
                double dx = rx[i] - rx[j];
                double dy = ry[i] - ry[j];
                double dz = rz[i] - rz[j];
                apply_pbc(dx, dy, dz);
                double r2 = dx*dx + dy*dy + dz*dz;
                if (r2 >= cutoff_sq) continue;
                
                double s2 = sigma2 / r2;
                double s6 = s2 * s2 * s2;
                double force_mag = 24.0 * epsilon * (2.0 * s6 * s6 - s6) / r2;
                
                double fx_ij = force_mag * dx;
                double fy_ij = force_mag * dy;
                double fz_ij = force_mag * dz;
                fx_priv[tid][i] += fx_ij;   fx_priv[tid][j] -= fx_ij;
                fy_priv[tid][i] += fy_ij;   fy_priv[tid][j] -= fy_ij;
                fz_priv[tid][i] += fz_ij;   fz_priv[tid][j] -= fz_ij;
                u_local += 4.0 * epsilon * (s6 * s6 - s6);
            }
            u_priv[tid] = u_local;
        }

        for (int t = 0; t < num_threads; ++t) {
            potential_energy += u_priv[t];
            for (int i = 0; i < N; ++i) {
                fx[i] += fx_priv[t][i]; fy[i] += fy_priv[t][i]; fz[i] += fz_priv[t][i];
            }
        }
#endif
        auto end_force = std::chrono::high_resolution_clock::now();
        time_force_calc_ns += std::chrono::duration_cast<std::chrono::nanoseconds>(end_force - start_force).count();
    }

    double get_kinetic_energy() {
#ifdef USE_CUDA
        double zero_k = 0.0;
        cudaMemcpy(d_K, &zero_k, sizeof(double), cudaMemcpyHostToDevice);
        int threads = 256;
        size_t shared_mem = threads * sizeof(double);
        compute_kinetic_energy_kernel<<<cuda_blocks_count, threads, shared_mem>>>(d_vx, d_vy, d_vz, N, mass, d_K);
        cudaDeviceSynchronize();
        double K_host = 0.0;
        cudaMemcpy(&K_host, d_K, sizeof(double), cudaMemcpyDeviceToHost);
        return K_host;
#else
        double K = 0.0;
        #pragma omp parallel for reduction(+:K)
        for (int i = 0; i < N; ++i)
            K += 0.5 * mass * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
        return K;
#endif
    }

    double get_temperature() { return (2.0 / 3.0) * get_kinetic_energy() / (N * kB); }

    double get_potential_energy_fast() const { return potential_energy; }

    void compute_full_inter() {
#ifdef USE_CUDA
        cuda_allocate();
#endif
        std::fill(fx.begin(), fx.end(), 0.0);
        std::fill(fy.begin(), fy.end(), 0.0);
        std::fill(fz.begin(), fz.end(), 0.0);
        potential_energy = 0.0;
        double sigma2 = sigma * sigma;

        for (int i = 0; i < N; ++i) {
            for (int j = i+1; j < N; ++j) {
                double dx = rx[i] - rx[j]; double dy = ry[i] - ry[j]; double dz = rz[i] - rz[j];
                apply_pbc(dx, dy, dz);
                double r2 = dx*dx + dy*dy + dz*dz;
                if (r2 >= cutoff_sq) continue;
                
                double s2 = sigma2 / r2; 
                double s6 = s2 * s2 * s2;
                double force_mag = 24.0 * epsilon * (2.0 * s6 * s6 - s6) / r2;
                
                double fx_ij = force_mag * dx; double fy_ij = force_mag * dy; double fz_ij = force_mag * dz;
                fx[i] += fx_ij; fy[i] += fy_ij; fz[i] += fz_ij;
                fx[j] -= fx_ij; fy[j] -= fy_ij; fz[j] -= fz_ij;
                potential_energy += 4.0 * epsilon * (s6 * s6 - s6);
            }
        }
#ifdef USE_CUDA
        sync_host_to_device();
#endif
    }

    void scale_velocities(double s) {
#ifdef USE_CUDA
        int threads = 256;
        scale_velocities_kernel<<<cuda_blocks_count, threads>>>(d_vx, d_vy, d_vz, N, s);
        cudaDeviceSynchronize();
#else
        #pragma omp parallel for
        for (int i = 0; i < N; ++i) {
            vx[i] *= s; vy[i] *= s; vz[i] *= s;
        }
#endif
    }

    void verlet_step(double dt) {
        double half_dt = 0.5 * dt;
        double g = 3.0 * N;

        if (thermostat_enabled) {
            double K = get_kinetic_energy();
            double G2 = (Q1 * v_xi1 * v_xi1 - kB * T_target) / Q2;
            v_xi2 += 0.25 * dt * G2;
            double G1 = (2.0 * K - g * kB * T_target) / Q1 - v_xi1 * v_xi2;
            v_xi1 += 0.25 * dt * G1;
            
            double s = std::exp(-v_xi1 * 0.5 * dt);
            scale_velocities(s);
            K *= (s * s);
            
            xi1 += v_xi1 * 0.5 * dt;
            xi2 += v_xi2 * 0.5 * dt;
            
            G1 = (2.0 * K - g * kB * T_target) / Q1 - v_xi1 * v_xi2;
            v_xi1 += 0.25 * dt * G1;
            G2 = (Q1 * v_xi1 * v_xi1 - kB * T_target) / Q2;
            v_xi2 += 0.25 * dt * G2;
        }

#ifdef USE_CUDA
        int threads = 256;
        verlet_step_part1_kernel<<<cuda_blocks_count, threads>>>(d_rx, d_ry, d_rz, d_vx, d_vy, d_vz, d_fx, d_fy, d_fz, N, dt, mass, L, inv_L);
        cudaDeviceSynchronize();
#else
        for (int i = 0; i < N; ++i) {
            vx[i] += half_dt * fx[i] / mass; vy[i] += half_dt * fy[i] / mass; vz[i] += half_dt * fz[i] / mass;
            rx[i] += dt * vx[i];
            if (rx[i] < 0.0) rx[i] += L; else if (rx[i] >= L) rx[i] -= L;
            if (ry[i] < 0.0) ry[i] += L; else if (ry[i] >= L) ry[i] -= L;
            if (rz[i] < 0.0) rz[i] += L; else if (rz[i] >= L) rz[i] -= L;
        }
#endif

        compute_forces();

#ifdef USE_CUDA
        int threads_p2 = 256;
        verlet_step_part2_kernel<<<cuda_blocks_count, threads_p2>>>(d_vx, d_vy, d_vz, d_fx, d_fy, d_fz, N, dt, mass);
        cudaDeviceSynchronize();
#else
        for (int i = 0; i < N; ++i) {
            vx[i] += half_dt * fx[i] / mass; vy[i] += half_dt * fy[i] / mass; vz[i] += half_dt * fz[i] / mass;
        }
#endif

        if (thermostat_enabled) {
            double K = get_kinetic_energy();
            double G2 = (Q1 * v_xi1 * v_xi1 - kB * T_target) / Q2;
            v_xi2 += 0.25 * dt * G2;
            double G1 = (2.0 * K - g * kB * T_target) / Q1 - v_xi1 * v_xi2;
            v_xi1 += 0.25 * dt * G1;
            
            double s = std::exp(-v_xi1 * 0.5 * dt);
            scale_velocities(s);
            K *= (s * s);
            
            xi1 += v_xi1 * 0.5 * dt;
            xi2 += v_xi2 * 0.5 * dt;
            
            G1 = (2.0 * K - g * kB * T_target) / Q1 - v_xi1 * v_xi2;
            v_xi1 += 0.25 * dt * G1;
            G2 = (Q1 * v_xi1 * v_xi1 - kB * T_target) / Q2;
            v_xi2 += 0.25 * dt * G2;
        }
    }

    void run_verlet(double dt, double total_time, int output_freq, const std::string& filename, double t_equil) {
        int steps = static_cast<int>(total_time / dt);
        if (steps <= 0) return;
        compute_forces();
        
        std::ofstream fout(filename);
        fout << "step;x;y;z\n";
        fout << std::fixed << std::setprecision(12); // Увеличена точность вывода для СИ масштабов (метры)
        double sum_E = 0.0, sum_E2 = 0.0;
        
        int equil_step = static_cast<int>(t_equil / dt);
        if (equil_step > steps) equil_step = steps;
        int cnt = 0;
        
        time_verlet_rebuild_ns = 0; time_check_disp_ns = 0; time_force_calc_ns = 0; time_io_ns = 0;
        count_rebuilds = 0;

        for (int s = 0; s <= steps; ++s) {
            double U = get_potential_energy_fast();
            double K = get_kinetic_energy();
            double E = U + K;
            
            if (s % output_freq == 0 || s == steps) {
                auto io_start = std::chrono::high_resolution_clock::now();
#ifdef USE_CUDA
                sync_device_to_host();
#endif
                for (int i = 0; i < N; ++i) fout << s << ";" << rx[i] << ";" << ry[i] << ";" << rz[i] << "\n";
                auto io_end = std::chrono::high_resolution_clock::now();
                time_io_ns += std::chrono::duration_cast<std::chrono::nanoseconds>(io_end - io_start).count();
            }

            if (s >= equil_step) { 
                sum_E += E; 
                sum_E2 += E * E; 
                cnt++; 
            }
            if (s < steps) verlet_step(dt);
            
            if (s % 500 == 0 || s == steps) std::cout << "Step: " << s << " / " << steps << " | T = " << (2.0 / 3.0) * K / (N * kB) << " K\n";
        }
        fout.close();

        std::cout << "\nTrajectory saved to " << filename << std::endl;
        if (cnt > 0) {
            double avg_E  = sum_E / cnt;
            double avg_E2 = sum_E2 / cnt;
            
            // Теплоемкость C_V в единицах СИ: [Дж/К]
            double Cv = (avg_E2 - avg_E*avg_E) / (kB * T_target * T_target);
            double NA = 6.02214076e23; // Число Авогадро
            
            std::cout << "=== THERMODYNAMIC PROPERTIES (SI UNITS VALIDATION) ===" << "\n";
            std::cout << "Total Heat Capacity C_V      = " << Cv << " J/K\n";
            std::cout << "C_V per single atom          = " << Cv / N << " J/(atom*K)  [or " << (Cv / N) / kB << " in units of k_B]\n";
            std::cout << "Molar Heat Capacity C_V      = " << (Cv / N) * NA << " J/(mol*K)\n";
            std::cout << "Specific Heat Capacity c_V   = " << Cv / (N * mass) << " J/(kg*K)\n";
            std::cout << "Statistical samples gathered = " << cnt << "\n";
        } else {
            std::cout << "\n[Warning] No samples gathered for C_V. Simulation time was shorter than equilibration time.\n";
        }

        double t_rebuild_ms = time_verlet_rebuild_ns / 1e6;
        double t_check_ms = time_check_disp_ns / 1e6;
        double t_force_ms = time_force_calc_ns / 1e6;
        double t_io_ms = time_io_ns / 1e6;
        double t_total_ms = t_rebuild_ms + t_check_ms + t_force_ms + t_io_ms;

        std::cout << "\n========== PERFORMANCE REPORT (GPU OPTIMIZED PIPELINE) ==========\n";
        std::cout << "Total tracked time    : " << std::fixed << std::setprecision(2) << t_total_ms << " ms\n";
        std::cout << "1. Verlet List Rebuild: " << t_rebuild_ms << " ms (" << (t_rebuild_ms / t_total_ms * 100.0) << "%)\n";
        std::cout << "   -> Total Rebuilds  : " << count_rebuilds << " times\n";
        std::cout << "2. Displacement Check : " << t_check_ms << " ms (" << (t_check_ms / t_total_ms * 100.0) << "%)\n";
        std::cout << "3. Force Calc + Energy: " << t_force_ms << " ms (" << (t_force_ms / t_total_ms * 100.0) << "%)\n";
        std::cout << "4. File Output (I/O)  : " << t_io_ms << " ms (" << (t_io_ms / t_total_ms * 100.0) << "%)\n";
        std::cout << "==================================================================\n\n";
    }

    const std::vector<double>& get_force_x() const { return fx; }
    const std::vector<double>& get_force_y() const { return fy; }
    const std::vector<double>& get_force_z() const { return fz; }
};
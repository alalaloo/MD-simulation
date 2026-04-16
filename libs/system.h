#pragma once
#include <vector>
#include <cmath>
#include <iostream>
#include <random>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <limits>
#include <array>

using namespace std;

class point3D {
public:
    vector<double> pos;
    vector<double> vel;
    vector<double> force;
    
    point3D() : pos(3, 0.0), vel(3, 0.0), force(3, 0.0) {}
    point3D(double x, double y, double z) : pos{x, y, z}, vel(3, 0.0), force(3, 0.0) {}
    
    void set_pos(double x, double y, double z) { pos = {x, y, z}; }
    void set_vel(double vx, double vy, double vz) { vel = {vx, vy, vz}; }
    void set_force(double fx, double fy, double fz) { force = {fx, fy, fz}; }
    void add_force(double fx, double fy, double fz) {
        force[0] += fx; force[1] += fy; force[2] += fz;
    }
};

// cell-list
struct Cell {
    vector<int> particle_indices;
    void clear() { particle_indices.clear(); }
};

class systemMD {
private:
    int N;                       // число частиц
    vector<point3D> particles;   // 1d вектор с частицами
    double mass;                 // масса частиц
    double cutoff;		 // радиус обрезания потенциала
    double cutoff_sq;		 
    double boxX, boxY, boxZ;	 // размер системы
    
    // параметры для потенциала Леннард-Джонса
    double lj_epsilon = 1.0;	 
    double lj_sigma   = 1.0;
    
    // Cell list параметры
    double cell_size;
    int gridX, gridY, gridZ;
    vector<Cell> cells;
    
    // периодические граничные условия
    void apply_pbc(double& dx, double& dy, double& dz) const {
        dx -= boxX * round(dx / boxX);
        dy -= boxY * round(dy / boxY);
        dz -= boxZ * round(dz / boxZ);
    }
    
    // преобразование координат в индекс
    int getCellIndex(double x, double y, double z) const {
        int ix = static_cast<int>(floor(x / cell_size));
        int iy = static_cast<int>(floor(y / cell_size));
        int iz = static_cast<int>(floor(z / cell_size));
        ix = (ix % gridX + gridX) % gridX;
        iy = (iy % gridY + gridY) % gridY;
        iz = (iz % gridZ + gridZ) % gridZ;
        return (ix * gridY + iy) * gridZ + iz;
    }
    
    void rebuildCellList() {
        for (auto& cell : cells) cell.clear();
        for (int i = 0; i < N; ++i) {
            const auto& p = particles[i];
            int cid = getCellIndex(p.pos[0], p.pos[1], p.pos[2]);
            cells[cid].particle_indices.push_back(i);
        }
    }

    bool lj_interaction(int i, int j, double dx, double dy, double dz,
                        bool add_energy, double* energy_acc) {
        double r2 = dx*dx + dy*dy + dz*dz;
        if (r2 >= cutoff_sq || r2 < 1e-12) return false;
        
        double inv2 = 1.0 / r2;
        double inv6 = inv2 * inv2 * inv2;
        double force_mag = 24.0 * lj_epsilon * (2.0 * inv6 * inv6 - inv6) * inv2;
        double fx = force_mag * dx;
        double fy = force_mag * dy;
        double fz = force_mag * dz;
        
        particles[i].add_force( fx,  fy,  fz);
        particles[j].add_force(-fx, -fy, -fz);
        
        if (add_energy && energy_acc) {
            *energy_acc += 4.0 * lj_epsilon * (inv6 * inv6 - inv6);
        }
        return true;
    }
    
    double lj_energy_pair(int i, int j, double dx, double dy, double dz) const {
        double r2 = dx*dx + dy*dy + dz*dz;
        if (r2 >= cutoff_sq || r2 < 1e-12) return 0.0;
        double inv2 = 1.0 / r2;
        double inv6 = inv2 * inv2 * inv2;
        return 4.0 * lj_epsilon * (inv6 * inv6 - inv6);
    }
    
public:
    systemMD(int n_particles, double lattice_constant, double particle_mass = 1.0, double cut_off = 3.0)
        : mass(particle_mass), cutoff(cut_off)
    {
        cutoff_sq = cutoff * cutoff;
        
        int n_per_side = static_cast<int>(round(cbrt(n_particles)));
        int nx = n_per_side;
        int ny = n_per_side;
        int nz = n_per_side;
        N = nx * ny * nz;
        particles.resize(N);
        
        int idx = 0;
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                for (int k = 0; k < nz; ++k) {
                    double x = i * lattice_constant;
                    double y = j * lattice_constant;
                    double z = k * lattice_constant;
                    particles[idx++].set_pos(x, y, z);
                }
            }
        }
        
        boxX = nx * lattice_constant;
        boxY = ny * lattice_constant;
        boxZ = nz * lattice_constant;
        
	// разбиение пространства на ячейки
        cell_size = cutoff;
        gridX = max(3, static_cast<int>(ceil(boxX / cell_size)));
        gridY = max(3, static_cast<int>(ceil(boxY / cell_size)));
        gridZ = max(3, static_cast<int>(ceil(boxZ / cell_size)));
        cells.resize(gridX * gridY * gridZ);
        rebuildCellList();
    }
    
    int getN() const { return N; }
    
    // распределение Максвелла
    void set_velocities(double T_target, double k_boltzmann = 1.0) {
        mt19937 gen(42); 
        normal_distribution<double> dist(0.0, sqrt(k_boltzmann * T_target / mass));
        
        double vx_cm = 0.0, vy_cm = 0.0, vz_cm = 0.0;
        for (auto& p : particles) {
            double vx = dist(gen);
            double vy = dist(gen);
            double vz = dist(gen);
            p.set_vel(vx, vy, vz);
            vx_cm += vx; vy_cm += vy; vz_cm += vz;
        }
        vx_cm /= N; vy_cm /= N; vz_cm /= N;
        
        double kinetic = 0.0;
        for (auto& p : particles) {
            double vx = p.vel[0] - vx_cm;
            double vy = p.vel[1] - vy_cm;
            double vz = p.vel[2] - vz_cm;
            p.set_vel(vx, vy, vz);
            kinetic += 0.5 * mass * (vx*vx + vy*vy + vz*vz);
        }
        double T_current = (2.0 / 3.0) * kinetic / (N * k_boltzmann);
        double scale = sqrt(T_target / T_current);
        for (auto& p : particles) {
            p.vel[0] *= scale;
            p.vel[1] *= scale;
            p.vel[2] *= scale;
        }
    }

    // вычисление силы полным перебором
    void compute_full_inter() {
        
        for (auto& p : particles) p.set_force(0.0, 0.0, 0.0);
        
        for (int i = 0; i < N; ++i) {
            for (int j = i+1; j < N; ++j) {
                double dx = particles[j].pos[0] - particles[i].pos[0];
                double dy = particles[j].pos[1] - particles[i].pos[1];
                double dz = particles[j].pos[2] - particles[i].pos[2];
                apply_pbc(dx, dy, dz);
                
                lj_interaction(i, j, dx, dy, dz, false, nullptr);
            }
        }
    }
    
    // вычисление силы методом cell-list
    void compute_forces() {
        rebuildCellList();
        for (auto& p : particles) p.set_force(0.0, 0.0, 0.0);
       
        for (int cx = 0; cx < gridX; ++cx) {
            for (int cy = 0; cy < gridY; ++cy) {
                for (int cz = 0; cz < gridZ; ++cz) {
                    int cell_idx = (cx * gridY + cy) * gridZ + cz;
                    const auto& cellA = cells[cell_idx].particle_indices;
                    
                    for (size_t a = 0; a < cellA.size(); ++a) {
                        for (size_t b = a+1; b < cellA.size(); ++b) {
                            int i = cellA[a];
                            int j = cellA[b];
                            double dx = particles[j].pos[0] - particles[i].pos[0];
                            double dy = particles[j].pos[1] - particles[i].pos[1];
                            double dz = particles[j].pos[2] - particles[i].pos[2];
                            apply_pbc(dx, dy, dz);
                            lj_interaction(i, j, dx, dy, dz, false, nullptr);
                        }
                    }
                    
                    for (int dcx = -1; dcx <= 1; ++dcx) {
                        for (int dcy = -1; dcy <= 1; ++dcy) {
                            for (int dcz = -1; dcz <= 1; ++dcz) {
                                if (dcx == 0 && dcy == 0 && dcz == 0) continue;
                                int nx = (cx + dcx + gridX) % gridX;
                                int ny = (cy + dcy + gridY) % gridY;
                                int nz = (cz + dcz + gridZ) % gridZ;
                                int neighbor_idx = (nx * gridY + ny) * gridZ + nz;
                                if (cell_idx < neighbor_idx) {
                                    const auto& cellB = cells[neighbor_idx].particle_indices;
                                    for (int i : cellA) {
                                        for (int j : cellB) {
                                            double dx = particles[j].pos[0] - particles[i].pos[0];
                                            double dy = particles[j].pos[1] - particles[i].pos[1];
                                            double dz = particles[j].pos[2] - particles[i].pos[2];
                                            apply_pbc(dx, dy, dz);
                                            lj_interaction(i, j, dx, dy, dz, false, nullptr);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    
    // вычисление энергии полным перебором
    double get_potential_energy_slow() const {
        double energy = 0.0;
        for (int i = 0; i < N; ++i) {
            for (int j = i+1; j < N; ++j) {
                double dx = particles[j].pos[0] - particles[i].pos[0];
                double dy = particles[j].pos[1] - particles[i].pos[1];
                double dz = particles[j].pos[2] - particles[i].pos[2];
                apply_pbc(dx, dy, dz);
                energy += lj_energy_pair(i, j, dx, dy, dz);
            }
        }
        return energy;
    }
    
    // вычисление энергии методом cell-list
    double get_potential_energy_fast() const {
        vector<Cell> local_cells(gridX * gridY * gridZ);
        for (int i = 0; i < N; ++i) {
            const auto& p = particles[i];
            int cid = getCellIndex(p.pos[0], p.pos[1], p.pos[2]);
            local_cells[cid].particle_indices.push_back(i);
        }
        
        double energy = 0.0;
        for (int cx = 0; cx < gridX; ++cx) {
            for (int cy = 0; cy < gridY; ++cy) {
                for (int cz = 0; cz < gridZ; ++cz) {
                    int cell_idx = (cx * gridY + cy) * gridZ + cz;
                    const auto& cellA = local_cells[cell_idx].particle_indices;
                    
                    for (size_t a = 0; a < cellA.size(); ++a) {
                        for (size_t b = a+1; b < cellA.size(); ++b) {
                            int i = cellA[a];
                            int j = cellA[b];
                            double dx = particles[j].pos[0] - particles[i].pos[0];
                            double dy = particles[j].pos[1] - particles[i].pos[1];
                            double dz = particles[j].pos[2] - particles[i].pos[2];
                            apply_pbc(dx, dy, dz);
                            energy += lj_energy_pair(i, j, dx, dy, dz);
                        }
                    }
                    
                    for (int dcx = -1; dcx <= 1; ++dcx) {
                        for (int dcy = -1; dcy <= 1; ++dcy) {
                            for (int dcz = -1; dcz <= 1; ++dcz) {
                                if (dcx == 0 && dcy == 0 && dcz == 0) continue;
                                int nx = (cx + dcx + gridX) % gridX;
                                int ny = (cy + dcy + gridY) % gridY;
                                int nz = (cz + dcz + gridZ) % gridZ;
                                int neighbor_idx = (nx * gridY + ny) * gridZ + nz;
                                if (cell_idx < neighbor_idx) {
                                    const auto& cellB = local_cells[neighbor_idx].particle_indices;
                                    for (int i : cellA) {
                                        for (int j : cellB) {
                                            double dx = particles[j].pos[0] - particles[i].pos[0];
                                            double dy = particles[j].pos[1] - particles[i].pos[1];
                                            double dz = particles[j].pos[2] - particles[i].pos[2];
                                            apply_pbc(dx, dy, dz);
                                            energy += lj_energy_pair(i, j, dx, dy, dz);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        return energy;
    }
    
    double get_kinetic_energy() const {
        double ke = 0.0;
        for (const auto& p : particles) {
            ke += 0.5 * mass * (p.vel[0]*p.vel[0] + p.vel[1]*p.vel[1] + p.vel[2]*p.vel[2]);
        }
        return ke;
    }
    
    double get_temperature(double k_boltzmann = 1.0) const {
        return (2.0 / 3.0) * get_kinetic_energy() / (N * k_boltzmann);
    }

    // один шаг скоростного Верле
    void verlet_step(double dt) {
        // gолушаг скоростей: v += 0.5 * dt * F / m
        for (auto& p : particles) {
            p.vel[0] += 0.5 * dt * p.force[0] / mass;
            p.vel[1] += 0.5 * dt * p.force[1] / mass;
            p.vel[2] += 0.5 * dt * p.force[2] / mass;
        }

        // обновление координат 
        for (auto& p : particles) {
            p.pos[0] += dt * p.vel[0];
            p.pos[1] += dt * p.vel[1];
            p.pos[2] += dt * p.vel[2];

            // периодические граничные условия
            p.pos[0] = fmod(p.pos[0], boxX);
            if (p.pos[0] < 0) p.pos[0] += boxX;
            p.pos[1] = fmod(p.pos[1], boxY);
            if (p.pos[1] < 0) p.pos[1] += boxY;
            p.pos[2] = fmod(p.pos[2], boxZ);
            if (p.pos[2] < 0) p.pos[2] += boxZ;
        }

        // новые силы F(t+dt)
        compute_forces();

        // завершающий полушаг скоростей
        for (auto& p : particles) {
            p.vel[0] += 0.5 * dt * p.force[0] / mass;
            p.vel[1] += 0.5 * dt * p.force[1] / mass;
            p.vel[2] += 0.5 * dt * p.force[2] / mass;
        }
    }

    // симуляция методов списков Верле
    void run_verlet(double dt, double total_time, int output_freq, const string& traj_filename) {
        int steps = static_cast<int>(total_time / dt);
        if (steps <= 0) return;

        ofstream traj_file(traj_filename);
        traj_file << "step,x,y,z,vx,vy,vz\n";
        traj_file << fixed << setprecision(6);

        // шаг 0
        for (int i = 0; i < N; ++i) {
            const auto& p = particles[i];
            traj_file << 0 << ","
                      << p.pos[0] << "," << p.pos[1] << "," << p.pos[2] << ","
                      << p.vel[0] << "," << p.vel[1] << "," << p.vel[2] << "\n";
        }

        for (int step = 1; step <= steps; ++step) {
            verlet_step(dt);

            if (step % output_freq == 0 || step == steps) {
                for (int i = 0; i < N; ++i) {
                    const auto& p = particles[i];
                    traj_file << step << ","
                              << p.pos[0] << "," << p.pos[1] << "," << p.pos[2] << ","
                              << p.vel[0] << "," << p.vel[1] << "," << p.vel[2] << "\n";
                }
            }
        }
        traj_file.close();
        cout << "Trajectory saved to: " << traj_filename << endl;
    }
    
    const point3D& getParticle(int idx) const { return particles[idx]; }
    
    void saveToFile(const string& filename) const {
        ofstream file(filename);
        if (!file.is_open()) {
            cerr << "Error opening file: " << filename << endl;
            return;
        }
        file << "index,x,y,z,vx,vy,vz,mass,fx,fy,fz\n";
        for (int i = 0; i < N; ++i) {
            const auto& p = particles[i];
            file << i << ","
                 << p.pos[0] << "," << p.pos[1] << "," << p.pos[2] << ","
                 << p.vel[0] << "," << p.vel[1] << "," << p.vel[2] << ","
                 << mass << ","
                 << p.force[0] << "," << p.force[1] << "," << p.force[2] << "\n";
        }
        file.close();
        cout << "Data saved to: " << filename << endl;
    }
};

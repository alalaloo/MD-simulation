#pragma once
#include <vector>
#include <cmath>
#include <iostream>
#include <random>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <chrono>

class systemMD {
public:
    int N;                         // число частиц
    double L;                      // размер кубической ячейки
    std::vector<double> rx, ry, rz; // координаты
    std::vector<double> vx, vy, vz; // скорости
    std::vector<double> fx, fy, fz; // силы
    double mass;                   // масса частицы
    double cutoff;                 // радиус обрезки
    double cutoff_sq;              // квадрат радиуса обрезки
    double skin;                   // "запас" для списка Верле
    double rlist;                  // радиус построения списка Верле = cutoff + skin
    double rlist_sq;

    // Параметры термостата Нозе–Хувера
    bool thermostat_enabled;
    double xi;                     // переменная термостата ξ
    double Q;                      // тепловая инерция
    double T_target;               // целевая температура
    double kB;                     // постоянная Больцмана

    struct VerletList {
        std::vector<int> i_list, j_list;       // пары (i,j)
        std::vector<double> rx0, ry0, rz0;     // позиции при последнем построении
        double max_disp;                       // максимальное смещение
        bool need_rebuild;
    } verlet;

    systemMD(int n_req, double a0, double particle_mass = 1.0,
             double cut_off = 3.0, double skin_depth = 0.3)
        : mass(particle_mass), cutoff(cut_off), skin(skin_depth)
    {
        cutoff_sq = cutoff * cutoff;
        rlist = cutoff + skin;
        rlist_sq = rlist * rlist;

        // ГЦК: 4 атома на кубическую ячейку
        int n_cells = std::max(1, (int)std::round(std::cbrt(n_req / 4.0)));
        N = 4 * n_cells * n_cells * n_cells;
        L = n_cells * a0;

        rx.resize(N); ry.resize(N); rz.resize(N);
        vx.resize(N); vy.resize(N); vz.resize(N);
        fx.resize(N); fy.resize(N); fz.resize(N);

        // Инициализация ГЦК решётки
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
        xi = 0.0;
        T_target = 1.0;
        kB = 1.0;

        verlet.max_disp = 0.0;
        verlet.need_rebuild = true;
    }

    int getN() const { return N; }

    // Инициализация скоростей
    void set_velocities(double T, double k_boltzmann = 1.0) {
        kB = k_boltzmann;
        T_target = T;
        std::mt19937 gen(42);
        std::normal_distribution<double> dist(0.0, std::sqrt(kB * T / mass));
        double sum_vx = 0, sum_vy = 0, sum_vz = 0;
        for (int i = 0; i < N; ++i) {
            vx[i] = dist(gen);
            vy[i] = dist(gen);
            vz[i] = dist(gen);
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
    }

    // термостат Нозе–Хувера
    void enable_nose_hoover(double T, double Q_param) {
        thermostat_enabled = true;
        T_target = T;
        Q = Q_param;
        xi = 0.0;
    }

    void disable_thermostat() { thermostat_enabled = false; }

    // Периодические граничные условия для разности координат
    void apply_pbc(double& dx, double& dy, double& dz) const {
        dx -= L * std::round(dx / L);
        dy -= L * std::round(dy / L);
        dz -= L * std::round(dz / L);
    }

    // Полный перебор сил (O(N²))
    void compute_full_inter() {
        std::fill(fx.begin(), fx.end(), 0.0);
        std::fill(fy.begin(), fy.end(), 0.0);
        std::fill(fz.begin(), fz.end(), 0.0);
        for (int i = 0; i < N; ++i) {
            for (int j = i+1; j < N; ++j) {
                double dx = rx[i] - rx[j];
                double dy = ry[i] - ry[j];
                double dz = rz[i] - rz[j];
                apply_pbc(dx, dy, dz);
                double r2 = dx*dx + dy*dy + dz*dz;
                if (r2 >= cutoff_sq) continue;
                double inv2 = 1.0 / r2;
                double inv6 = inv2 * inv2 * inv2;
                double force_mag = 24.0 * (2.0 * inv6*inv6 - inv6) * inv2;
                double fx_ij = force_mag * dx;
                double fy_ij = force_mag * dy;
                double fz_ij = force_mag * dz;
                fx[i] += fx_ij; fy[i] += fy_ij; fz[i] += fz_ij;
                fx[j] -= fx_ij; fy[j] -= fy_ij; fz[j] -= fz_ij;
            }
        }
    }

    // Построение Verlet-списка
    void build_verlet_list() {
        verlet.i_list.clear(); verlet.j_list.clear();
        verlet.rx0 = rx; verlet.ry0 = ry; verlet.rz0 = rz;
        verlet.max_disp = 0.0;
        verlet.need_rebuild = false;
        for (int i = 0; i < N; ++i) {
            for (int j = i+1; j < N; ++j) {
                double dx = rx[i] - rx[j];
                double dy = ry[i] - ry[j];
                double dz = rz[i] - rz[j];
                apply_pbc(dx, dy, dz);
                if (dx*dx + dy*dy + dz*dz < rlist_sq) {
                    verlet.i_list.push_back(i);
                    verlet.j_list.push_back(j);
                }
            }
        }
    }

    void check_displacement() {
        double max2 = 0.0;
        for (int i = 0; i < N; ++i) {
            double dx = rx[i] - verlet.rx0[i];
            double dy = ry[i] - verlet.ry0[i];
            double dz = rz[i] - verlet.rz0[i];
            apply_pbc(dx, dy, dz);
            double d2 = dx*dx + dy*dy + dz*dz;
            if (d2 > max2) max2 = d2;
        }
        verlet.max_disp = std::sqrt(max2);
        if (verlet.max_disp > 0.5 * skin) verlet.need_rebuild = true;
    }

    // Вычисление сил с использованием Verlet-списка
    void compute_forces() {
        if (verlet.need_rebuild) build_verlet_list();
        else {
            check_displacement();
            if (verlet.need_rebuild) build_verlet_list();
        }
        std::fill(fx.begin(), fx.end(), 0.0);
        std::fill(fy.begin(), fy.end(), 0.0);
        std::fill(fz.begin(), fz.end(), 0.0);
        for (size_t k = 0; k < verlet.i_list.size(); ++k) {
            int i = verlet.i_list[k];
            int j = verlet.j_list[k];
            double dx = rx[i] - rx[j];
            double dy = ry[i] - ry[j];
            double dz = rz[i] - rz[j];
            apply_pbc(dx, dy, dz);
            double r2 = dx*dx + dy*dy + dz*dz;
            if (r2 >= cutoff_sq) continue;
            double inv2 = 1.0 / r2;
            double inv6 = inv2 * inv2 * inv2;
            double force_mag = 24.0 * (2.0 * inv6*inv6 - inv6) * inv2;
            double fx_ij = force_mag * dx;
            double fy_ij = force_mag * dy;
            double fz_ij = force_mag * dz;
            fx[i] += fx_ij; fy[i] += fy_ij; fz[i] += fz_ij;
            fx[j] -= fx_ij; fy[j] -= fy_ij; fz[j] -= fz_ij;
        }
    }

    // Потенциальная энергия (полный перебор)
    double get_potential_energy_slow() const {
        double U = 0.0;
        for (int i = 0; i < N; ++i) {
            for (int j = i+1; j < N; ++j) {
                double dx = rx[i] - rx[j];
                double dy = ry[i] - ry[j];
                double dz = rz[i] - rz[j];
                apply_pbc(dx, dy, dz);
                double r2 = dx*dx + dy*dy + dz*dz;
                if (r2 >= cutoff_sq) continue;
                double inv2 = 1.0 / r2;
                double inv6 = inv2 * inv2 * inv2;
                U += 4.0 * (inv6*inv6 - inv6);
            }
        }
        return U;
    }

    // Потенциальная энергия (по Verlet-списку)
    double get_potential_energy_fast() const {
        double U = 0.0;
        for (size_t k = 0; k < verlet.i_list.size(); ++k) {
            int i = verlet.i_list[k];
            int j = verlet.j_list[k];
            double dx = rx[i] - rx[j];
            double dy = ry[i] - ry[j];
            double dz = rz[i] - rz[j];
            apply_pbc(dx, dy, dz);
            double r2 = dx*dx + dy*dy + dz*dz;
            if (r2 >= cutoff_sq) continue;
            double inv2 = 1.0 / r2;
            double inv6 = inv2 * inv2 * inv2;
            U += 4.0 * (inv6*inv6 - inv6);
        }
        return U;
    }

    double get_kinetic_energy() const {
        double K = 0.0;
        for (int i = 0; i < N; ++i)
            K += 0.5 * mass * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
        return K;
    }

    double get_temperature() const {
        return (2.0 / 3.0) * get_kinetic_energy() / (N * kB);
    }

    void verlet_step(double dt) {
        double half_dt = 0.5 * dt;

        // Первый полушаг скоростей (без термостата)
        for (int i = 0; i < N; ++i) {
            vx[i] += half_dt * fx[i] / mass;
            vy[i] += half_dt * fy[i] / mass;
            vz[i] += half_dt * fz[i] / mass;
        }

        if (thermostat_enabled) {
            // Термостат: первый полушаг
            double K = get_kinetic_energy();
            double g = 3.0 * N;
            double G1 = (2.0 * K - g * kB * T_target) / Q;
            xi += half_dt * G1;
            double s = std::exp(-xi * half_dt);
            for (int i = 0; i < N; ++i) {
                vx[i] *= s; vy[i] *= s; vz[i] *= s;
            }
        }

        // Обновление координат и периодические границы
        for (int i = 0; i < N; ++i) {
            rx[i] += dt * vx[i];
            if (rx[i] < 0) rx[i] += L; if (rx[i] >= L) rx[i] -= L;
            ry[i] += dt * vy[i];
            if (ry[i] < 0) ry[i] += L; if (ry[i] >= L) ry[i] -= L;
            rz[i] += dt * vz[i];
            if (rz[i] < 0) rz[i] += L; if (rz[i] >= L) rz[i] -= L;
        }

        compute_forces();

        // Второй полушаг скоростей
        for (int i = 0; i < N; ++i) {
            vx[i] += half_dt * fx[i] / mass;
            vy[i] += half_dt * fy[i] / mass;
            vz[i] += half_dt * fz[i] / mass;
        }

        if (thermostat_enabled) {
            // Термостат: второй полушаг
            double K = get_kinetic_energy();
            double g = 3.0 * N;
            double G2 = (2.0 * K - g * kB * T_target) / Q;
            xi += half_dt * G2;
            double s = std::exp(-xi * half_dt);
            for (int i = 0; i < N; ++i) {
                vx[i] *= s; vy[i] *= s; vz[i] *= s;
            }
        }
    }

    void run_verlet(double dt, double total_time, int output_freq,
                    const std::string& filename) {
        int steps = static_cast<int>(total_time / dt);
        if (steps <= 0) return;

        compute_forces();

        std::ofstream fout(filename);
        fout << "step;x;y;z;vx;vy;vz;U;K;E;T\n";
        fout << std::fixed << std::setprecision(6);

        double sum_E = 0.0, sum_E2 = 0.0;
        int equil = steps / 5;  
        int cnt = 0;

        for (int s = 0; s <= steps; ++s) {
            double U = get_potential_energy_fast();
            double K = get_kinetic_energy();
            double E = U + K;
            double T = get_temperature();
            if (s % output_freq == 0 || s == steps) {
                for (int i = 0; i < N; ++i) {
                    fout << s << ";"
                         << rx[i] << ";" << ry[i] << ";" << rz[i] << ";"
                         << vx[i] << ";" << vy[i] << ";" << vz[i] << ";"
                         << U << ";" << K << ";" << E << ";" << T << "\n";
                }
            }
            if (s >= equil) {
                sum_E  += E;
                sum_E2 += E * E;
                cnt++;
            }
            if (s < steps) verlet_step(dt);
        }
        fout.close();
        std::cout << "Trajectory saved to " << filename << std::endl;
        if (cnt > 0) {
            double avg_E  = sum_E / cnt;
            double avg_E2 = sum_E2 / cnt;
            double Cv = (avg_E2 - avg_E*avg_E) / (kB * T_target * T_target);
            std::cout << "C_V = " << Cv / N << " (units of k_B), samples = " << cnt << "\n";
        }
    }

    const std::vector<double>& get_force_x() const { return fx; }
    const std::vector<double>& get_force_y() const { return fy; }
    const std::vector<double>& get_force_z() const { return fz; }
};

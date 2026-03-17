#include <vector>
#include <cmath>
#include <iostream>
#include <random>
#include <fstream>
#include <iomanip>
#include <algorithm> 
#include <limits>   

using namespace std;

class point3D {
private:
    vector<double> pos;   // координата
    vector<double> force; // сила
    vector<double> vel;   // скорость

public:
    point3D() : pos(3, 0), force(3, 0), vel(3, 0) {}
    point3D(double x, double y, double z) : pos({x, y, z}), force(3, 0), vel(3, 0) {}

    void set_pos(double x, double y, double z) {
        pos = {x, y, z};
    }
    void set_force(double fx, double fy, double fz) {
        force = {fx, fy, fz};
    }
    void set_vel(double vx, double vy, double vz) {
        vel = {vx, vy, vz};
    }
    
    void add_force(double fx, double fy, double fz) {
        force[0] += fx;
        force[1] += fy;
        force[2] += fz;
    }

    vector<double> get_pos() const { return pos; }
    vector<double> get_force() const { return force; }
    vector<double> get_vel() const { return vel; }
};

// Блоки с частицами
struct Cell {
    vector<int> particle_indices;

    void clear() { 
        particle_indices.clear(); 
        particle_indices.shrink_to_fit();
    }
};

class systemMD {
private:
    int N;                          // число частиц
    int sizeX, sizeY, sizeZ;        // размеры решетки
    point3D*** lattice;             // 3D массив частиц
    double avg_mass;

    double cutoff;                  // радиус обрезания потенциала
    double boxX, boxY, boxZ;        // размеры симуляционной коробки
    double lj_sigma;                // параметр sigma Леннарда-Джонса
    double lj_epsilon;              // параметр epsilon Леннарда-Джонса
    
    double cell_size;               // размер ячейки 
    int gridX, gridY, gridZ;        // количество ячеек по осям
    vector<Cell> cells;             // вектор ячеек

    int getLinearIndex(int i, int j, int k) const {
        return i * (sizeY * sizeZ) + j * sizeZ + k;
    }

    int getCellIndex(double x, double y, double z) const {
        int ix = static_cast<int>(floor(x / cell_size));
        int iy = static_cast<int>(floor(y / cell_size));
        int iz = static_cast<int>(floor(z / cell_size));
        
        // периодические ГУ
        if (ix < 0) ix += gridX; if (ix >= gridX) ix -= gridX;
        if (iy < 0) iy += gridY; if (iy >= gridY) iy -= gridY;
        if (iz < 0) iz += gridZ; if (iz >= gridZ) iz -= gridZ;
        
        return ix * (gridY * gridZ) + iy * gridZ + iz;
    }

    void rebuildCellList() {
        for (auto& cell : cells) cell.clear();
        
        for (int i = 0; i < sizeX; ++i) {
            for (int j = 0; j < sizeY; ++j) {
                for (int k = 0; k < sizeZ; ++k) {
                    vector<double> pos = lattice[i][j][k].get_pos();
                    int cell_idx = getCellIndex(pos[0], pos[1], pos[2]);
                    cells[cell_idx].particle_indices.push_back(getLinearIndex(i, j, k));
                }
            }
        }
    }

    void apply_pbc(double& dx, double& dy, double& dz) const {
        if (dx > boxX * 0.5) dx -= boxX;
        if (dx < -boxX * 0.5) dx += boxX;
        if (dy > boxY * 0.5) dy -= boxY;
        if (dy < -boxY * 0.5) dy += boxY;
        if (dz > boxZ * 0.5) dz -= boxZ;
        if (dz < -boxZ * 0.5) dz += boxZ;
    }

    void compute_lj_pair(int idx1, int idx2, double cutoff_sq) {

	int i1 = idx1 / (sizeY * sizeZ);
        int j1 = (idx1 / sizeZ) % sizeY;
        int k1 = idx1 % sizeZ;

        int i2 = idx2 / (sizeY * sizeZ);
        int j2 = (idx2 / sizeZ) % sizeY;
        int k2 = idx2 % sizeZ;

        vector<double> pos1 = lattice[i1][j1][k1].get_pos();
        vector<double> pos2 = lattice[i2][j2][k2].get_pos();

        double dx = pos2[0] - pos1[0];
        double dy = pos2[1] - pos1[1];
        double dz = pos2[2] - pos1[2];

        apply_pbc(dx, dy, dz);

        double r2 = dx*dx + dy*dy + dz*dz;

        if ((r2 < cutoff_sq) && (r2 > 1e-12)) {
            double r = sqrt(r2);
            double sig_r = lj_sigma / r;
            double sig_r6 = sig_r * sig_r * sig_r * sig_r * sig_r * sig_r;
            double sig_r12 = sig_r6 * sig_r6;

            double force_mag = 24.0 * lj_epsilon * (2.0 * sig_r12 - sig_r6) / r2;

            double fx = force_mag * dx;
            double fy = force_mag * dy;
            double fz = force_mag * dz;

            lattice[i1][j1][k1].add_force(fx, fy, fz);
            lattice[i2][j2][k2].add_force(-fx, -fy, -fz);
        }
    }

    void set_null_force() {
	    for (int i = 0; i < sizeX; ++i)
            	for (int j = 0; j < sizeY; ++j)
                	for (int k = 0; k < sizeZ; ++k)
                    		lattice[i][j][k].set_force(0, 0, 0);
    }

public:
    systemMD() : N(0), sizeX(0), sizeY(0), sizeZ(0), avg_mass(1.0), 
                 cutoff(2.5), lj_sigma(1.0), lj_epsilon(1.0) {
        lattice = nullptr;
        boxX = boxY = boxZ = 0;
        gridX = gridY = gridZ = 0;
    }

    systemMD(int n, double distance, double mass = 1.0, double cut_off = 2.5) 
        : avg_mass(mass), cutoff(cut_off), lj_sigma(1.0), lj_epsilon(1.0) {
        
        if (n >= 0) {
            sizeX = sizeY = sizeZ = static_cast<int>(ceil(cbrt(n)));
            N = sizeX * sizeY * sizeZ;
            
            // Выделение памяти под решетку
            lattice = new point3D**[sizeX];
            for (int i = 0; i < sizeX; ++i) {
                lattice[i] = new point3D*[sizeY];
                for (int j = 0; j < sizeY; ++j) {
                    lattice[i][j] = new point3D[sizeZ];
                }
            }
            
            // Заполнение частицами
            for (int i = 0; i < sizeX; ++i) {
                for (int j = 0; j < sizeY; ++j) {
                    for (int k = 0; k < sizeZ; ++k) {
                        double x = i * distance;
                        double y = j * distance;
                        double z = k * distance;
                        lattice[i][j][k] = point3D(x, y, z);
                    }
                }
            }
            
            // Размеры коробки
            boxX = sizeX * distance;
            boxY = sizeY * distance;
            boxZ = sizeZ * distance;
            
            // Инициализация Cell List
            cell_size = max(cutoff, distance);
            gridX = static_cast<int>(ceil(boxX / cell_size));
            gridY = static_cast<int>(ceil(boxY / cell_size));
            gridZ = static_cast<int>(ceil(boxZ / cell_size));
            cells.resize(gridX * gridY * gridZ);
            
            rebuildCellList();
        } else {
            sizeX = sizeY = sizeZ = 0;
            N = 0;
            lattice = nullptr;
            cout << "(N < 0) =>> lattice = nullptr!!!" << endl;
        }
    }

    ~systemMD() {
        if (lattice != nullptr) {
            for (int i = 0; i < sizeX; ++i) {
                for (int j = 0; j < sizeY; ++j) {
                    delete[] lattice[i][j];
                }
                delete[] lattice[i];
            }
            delete[] lattice;
        }
    }

    int getSizeX() const { return sizeX; }
    int getSizeY() const { return sizeY; }
    int getSizeZ() const { return sizeZ; }
    int getN() const { return N; }
    
    point3D& getPoint(int i, int j, int k) { return lattice[i][j][k]; }
    const point3D& getPoint(int i, int j, int k) const { return lattice[i][j][k]; }

    void set_velocities(double t_target = 298.0, double k_boltzmann = 1.0) {
        random_device rd;
        mt19937 gen(rd());
        uniform_real_distribution<double> dist(0.0, 1.0);
        double sigma = sqrt(k_boltzmann * t_target / avg_mass);
        
        vector<double> all_vx, all_vy, all_vz;
        all_vx.reserve(N); all_vy.reserve(N); all_vz.reserve(N);
        
        for (int i = 0; i < sizeX; ++i) {
            for (int j = 0; j < sizeY; ++j) {
                for (int k = 0; k < sizeZ; ++k) {
                    double u1 = max(dist(gen), numeric_limits<double>::min());
                    double u2 = dist(gen);
                    double vx = sigma * sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
                    
                    u1 = max(dist(gen), numeric_limits<double>::min());
                    u2 = dist(gen);
                    double vy = sigma * sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
                    
                    u1 = max(dist(gen), numeric_limits<double>::min());
                    u2 = dist(gen);
                    double vz = sigma * sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
                    
                    lattice[i][j][k].set_vel(vx, vy, vz);
                    all_vx.push_back(vx); all_vy.push_back(vy); all_vz.push_back(vz);
                }
            }
        }
        
        // Удаление дрейфа ЦМ
        double vx_cm = 0, vy_cm = 0, vz_cm = 0;
        for (int idx = 0; idx < N; ++idx) {
            vx_cm += all_vx[idx]; vy_cm += all_vy[idx]; vz_cm += all_vz[idx];
        }
        vx_cm /= N; vy_cm /= N; vz_cm /= N;
        
        int idx = 0;
        for (int i = 0; i < sizeX; ++i) {
            for (int j = 0; j < sizeY; ++j) {
                for (int k = 0; k < sizeZ; ++k) {
                    vector<double> vel = lattice[i][j][k].get_vel();
                    lattice[i][j][k].set_vel(vel[0] - vx_cm, vel[1] - vy_cm, vel[2] - vz_cm);
                    all_vx[idx] = vel[0] - vx_cm;
                    all_vy[idx] = vel[1] - vy_cm;
                    all_vz[idx] = vel[2] - vz_cm;
                    idx++;
                }
            }
        }
        
        // Масштабирование до t_target
        double kinetic_energy = 0;
        for (int idx = 0; idx < N; ++idx) {
            kinetic_energy += 0.5 * avg_mass * (all_vx[idx]*all_vx[idx] + all_vy[idx]*all_vy[idx] + all_vz[idx]*all_vz[idx]);
        }
        double t_current = (2.0 / 3.0) * kinetic_energy / (N * k_boltzmann);
        double scale_factor = sqrt(t_target / t_current);
        
        for (int i = 0; i < sizeX; ++i) {
            for (int j = 0; j < sizeY; ++j) {
                for (int k = 0; k < sizeZ; ++k) {
                    vector<double> vel = lattice[i][j][k].get_vel();
                    lattice[i][j][k].set_vel(vel[0] * scale_factor, vel[1] * scale_factor, vel[2] * scale_factor);
                }
            }
        }
    }

    // полный перебор
    void compute_full_inter() {
	
	set_null_force();
        
        double cutoff_sq = cutoff * cutoff;
        
        for (int i = 0; i < sizeX; ++i) {
            for (int j = 0; j < sizeY; ++j) {
                for (int k = 0; k < sizeZ; ++k) {
                    int idx1 = getLinearIndex(i, j, k);
                    vector<double> pos1 = lattice[i][j][k].get_pos();
                    
                    for (int ii = i; ii < sizeX; ++ii) {
                        int jj_start = (ii == i) ? j + 1 : 0;
                        for (int jj = jj_start; jj < sizeY; ++jj) {
                            int kk_start = (ii == i && jj == j) ? k + 1 : 0;
                            for (int kk = kk_start; kk < sizeZ; ++kk) {
                                int idx2 = getLinearIndex(ii, jj, kk);
                                vector<double> pos2 = lattice[ii][jj][kk].get_pos();
                                
                                double dx = pos2[0] - pos1[0];
                                double dy = pos2[1] - pos1[1];
                                double dz = pos2[2] - pos1[2];
                                
                                apply_pbc(dx, dy, dz);
                                
                                double r2 = dx*dx + dy*dy + dz*dz;
                                
                                if (r2 < cutoff_sq && r2 > 1e-12) {
                                    double r = sqrt(r2);
                                    double sig_r = lj_sigma / r;
                                    double sig_r6 = pow(sig_r, 6);
                                    double sig_r12 = sig_r6 * sig_r6;
                                    double force_mag = 24.0 * lj_epsilon * (2.0 * sig_r12 - sig_r6) / r2;
                                    
                                    lattice[i][j][k].add_force(force_mag * dx, force_mag * dy, force_mag * dz);
                                    lattice[ii][jj][kk].add_force(-force_mag * dx, -force_mag * dy, -force_mag * dz);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // перебор в соседних ячейках
    void compute_forces() {

        rebuildCellList(); 
        
        set_null_force();
        
        double cutoff_sq = cutoff * cutoff;
        
        for (int ix = 0; ix < gridX; ++ix) {
            for (int iy = 0; iy < gridY; ++iy) {
                for (int iz = 0; iz < gridZ; ++iz) {
                    int current_idx = ix * (gridY * gridZ) + iy * gridZ + iz;
                    const vector<int>& current_particles = cells[current_idx].particle_indices;
                    
                    if (current_particles.empty()) continue;
                    
		    // перебор соседних ячеек
                    for (int dx = -1; dx <= 1; ++dx) {
                        for (int dy = -1; dy <= 1; ++dy) {
                            for (int dz = -1; dz <= 1; ++dz) {
                                int nx = ix + dx;
                                int ny = iy + dy;
                                int nz = iz + dz;
                                
                                if (nx < 0) nx += gridX; if (nx >= gridX) nx -= gridX;
                                if (ny < 0) ny += gridY; if (ny >= gridY) ny -= gridY;
                                if (nz < 0) nz += gridZ; if (nz >= gridZ) nz -= gridZ;
                                
                                int neighbor_idx = nx * (gridY * gridZ) + ny * gridZ + nz;
                                const vector<int>& neighbor_particles = cells[neighbor_idx].particle_indices;
                                
                                if (neighbor_particles.empty()) continue;
                                
                                for (int idx_i : current_particles) {
                                    for (int idx_j : neighbor_particles) {
                                        if (idx_i >= idx_j && dx == 0 && dy == 0 && dz == 0) continue;
                                        compute_lj_pair(idx_i, idx_j, cutoff_sq);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    double get_potential_energy() {
        double energy = 0;
        double cutoff_sq = cutoff * cutoff;

        for (int i = 0; i < sizeX; ++i) {
            for (int j = 0; j < sizeY; ++j) {
                for (int k = 0; k < sizeZ; ++k) {

                    int idx1 = getLinearIndex(i, j, k);
                    vector<double> pos1 = lattice[i][j][k].get_pos();

                    for (int ii = i; ii < sizeX; ++ii) {
                        int jj_start = (ii == i) ? j + 1 : 0;
                        for (int jj = jj_start; jj < sizeY; ++jj) {
                            int kk_start = (ii == i && jj == j) ? k + 1 : 0;
                            for (int kk = kk_start; kk < sizeZ; ++kk) {

                                vector<double> pos2 = lattice[ii][jj][kk].get_pos();
                                double dx = pos2[0] - pos1[0];
                                double dy = pos2[1] - pos1[1];
                                double dz = pos2[2] - pos1[2];
                                apply_pbc(dx, dy, dz);
                                double r2 = dx*dx + dy*dy + dz*dz;
                                if (r2 < cutoff_sq && r2 > 1e-12) {
                                    double r = sqrt(r2);
                                    double sig_r = lj_sigma / r;
                                    double sig_r6 = pow(sig_r, 6);
                                    energy += 4.0 * lj_epsilon * (sig_r6 * sig_r6 - sig_r6);
                                }
                            }
                        }
                    }
                }
            }
        }
        return energy;
    }

    double get_kinetic_energy() {
        double energy = 0;
        for (int i = 0; i < sizeX; ++i)
            for (int j = 0; j < sizeY; ++j)
                for (int k = 0; k < sizeZ; ++k) {

                    vector<double> vel = lattice[i][j][k].get_vel();
                    energy += 0.5 * avg_mass * (vel[0]*vel[0] + vel[1]*vel[1] + vel[2]*vel[2]);
                }
        return energy;
    }

    double get_temperature(double k_boltzmann = 1.0) {
        return (2.0 / 3.0) * get_kinetic_energy() / (N * k_boltzmann);
    }

    void saveToFile(const string& filename) const {
        ofstream file(filename);
        if (!file.is_open()) { cerr << "FILE OPENING ERROR: " << filename << endl; return; }
        file << "index;x;y;z;vx;vy;vz;mass;forceX;forceY;forceZ\n";
        int index = 0;
        for (int i = 0; i < sizeX; ++i) {
            for (int j = 0; j < sizeY; ++j) {
                for (int k = 0; k < sizeZ; ++k) {
                    const point3D& p = lattice[i][j][k];
                    vector<double> pos = p.get_pos();
                    vector<double> vel = p.get_vel();
		    vector<double> force = p.get_force();
                    file << index++ << ";" << pos[0] << ";" << pos[1] << ";" << pos[2] << ";"
                         << vel[0] << ";" << vel[1] << ";" << vel[2] << ";" << avg_mass << ";" << force[0] << ";" << force[1] << ";" << force[2] << "\n";
                }
            }
        }
        file.close();
        cout << "The data saved to: " << filename << endl;
    }

    vector<double> getAllVx() const {
        vector<double> all_vx; all_vx.reserve(N);
        for (int i = 0; i < sizeX; ++i)
            for (int j = 0; j < sizeY; ++j)
                for (int k = 0; k < sizeZ; ++k)
                    all_vx.push_back(lattice[i][j][k].get_vel()[0]);
        return all_vx;
    }
    vector<double> getAllVy() const {
        vector<double> all_vy; all_vy.reserve(N);
        for (int i = 0; i < sizeX; ++i)
            for (int j = 0; j < sizeY; ++j)
                for (int k = 0; k < sizeZ; ++k)
                    all_vy.push_back(lattice[i][j][k].get_vel()[1]);
        return all_vy;
    }
    vector<double> getAllVz() const {
        vector<double> all_vz; all_vz.reserve(N);
        for (int i = 0; i < sizeX; ++i)
            for (int j = 0; j < sizeY; ++j)
                for (int k = 0; k < sizeZ; ++k)
                    all_vz.push_back(lattice[i][j][k].get_vel()[2]);
        return all_vz;
    }
    vector<double> getAllSpeed() const {
        vector<double> speeds; speeds.reserve(N);
        for (int i = 0; i < sizeX; ++i)
            for (int j = 0; j < sizeY; ++j)
                for (int k = 0; k < sizeZ; ++k) {
                    vector<double> vel = lattice[i][j][k].get_vel();
                    speeds.push_back(sqrt(vel[0]*vel[0] + vel[1]*vel[1] + vel[2]*vel[2]));
                }
        return speeds;
    }
};


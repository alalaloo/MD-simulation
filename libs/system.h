#include <vector>
#include <cmath>
#include <iostream>
#include <random>
#include <fstream>
#include <iomanip>

using namespace std;

class point3D {
	private:
		double mass;		//масса
		vector<double> pos;	//координата
		vector<double> force;	//сила
		vector<double> vel;	//скорость
	
	public:
		point3D() : mass(0), pos(3, 0), force(3, 0), vel(3, 0) {}

		point3D(double m, double x, double y, double z): mass(m), pos({x, y, z}), force(3, 0), vel(3, 0) {}

		
		void set_pos(double x, double y, double z) {
			pos = {x, y, z};
		}

		void set_force(double fx, double fy, double fz) {
			force = {fx, fy, fz};
		}

		void set_vel(double vx, double vy, double vz) {
			vel = {vx, vy, vz};
		}

		void set_mass(double m) {
			mass = m;
		}

		
		vector<double> get_pos() {
			return pos;
		}

		vector<double> get_force() {
			return force;
		}

		vector<double> get_vel() {
			return vel;
		}

		double get_mass() {
			return mass;
		}
};

class system {
	private:
		int N;				//число частиц
		int sizeX, sizeY, sizeZ;	//размеры решетки
		point3D*** lattice;		//инициализация кристаллической решетки
		double avg_mass;
		
	public:
		system() : N(0), sizeX(0), sizeY(0), sizeZ(0), avg_mass(1.0) {
			lattice = nullptr;
		}

		system(int n, double distance, double mass = 1.0) : avg_mass(mass) {
			if (n >= 0) {
				sizeX = sizeY = sizeZ = static_cast<int>(ceil(cbrt(n)));	//достраеваем до куба (N >= n)
				N = sizeX * sizeY * sizeZ;
				
				lattice = new point3D**[sizeX];
				for (int i = 0; i < sizeX; ++i) {
					lattice[i] = new point3D*[sizeY];
					for (int j = 0; j < sizeY; ++j) {
						lattice[i][j] = new point3D[sizeZ];
					} 
				}
				
				for (int i = 0; i < sizeX; ++i) {
					for (int j = 0; j < sizeY; ++j) {
						for (int k = 0; k < sizeZ; ++k) {
							double x = i * distance;
							double y = j * distance;
							double z = k * distance;
							
							lattice[i][j][k] = point3D(mass, x, y, z);
						}
					}
				}	
			} else {
				sizeX = sizeY = sizeZ = 0;
				N = 0;
				lattice = nullptr;
				cout << "(N < 0) =>> lattice = nullptr!!!" << endl;
			}
		}

		~system() {
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
		
		point3D& getPoint(int i, int j, int k) { 
			return lattice[i][j][k]; 
		}
		
		const point3D& getPoint(int i, int j, int k) const { 
			return lattice[i][j][k]; 
		}

		void set_velocities(double t_target = 298.0, double k_boltzmann = 1.0) {
    			random_device rd; 
    			mt19937 gen(rd()); // генератор Мерсенна Твистера
    			uniform_real_distribution<double> dist(0.0, 1.0);

			double sigma = sqrt(k_boltzmann * t_target / avg_mass);
			
			vector<double> all_vx, all_vy, all_vz;
			all_vx.reserve(N);
			all_vy.reserve(N);
			all_vz.reserve(N);
			
			// метод Бокса-Мюллера
			for (int i = 0; i < sizeX; ++i) {
				for (int j = 0; j < sizeY; ++j) {
			    		for (int k = 0; k < sizeZ; ++k) {
						
						// x
						double u1 = dist(gen);
						double u2 = dist(gen);
						u1 = max(u1, numeric_limits<double>::min());
						double vx = sigma * sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
				
						// y
						u1 = dist(gen);
						u2 = dist(gen);
						u1 = max(u1, numeric_limits<double>::min());
						double vy = sigma * sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);

						// z
						u1 = dist(gen);
						u2 = dist(gen);
						u1 = max(u1, numeric_limits<double>::min());
						double vz = sigma * sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
				
						lattice[i][j][k].set_vel(vx, vy, vz);
				
						all_vx.push_back(vx);
						all_vy.push_back(vy);
						all_vz.push_back(vz);
			    		}
				}
			}

			// удаление дрейфа цм
			double vx_cm = 0, vy_cm = 0, vz_cm = 0;

			for (int idx = 0; idx < N; ++idx) {
				vx_cm += all_vx[idx];
				vy_cm += all_vy[idx];
				vz_cm += all_vz[idx];
			}

			vx_cm /= N;
			vy_cm /= N;
			vz_cm /= N;
			
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
			
			// масштабирование до t_target
			double kinetic_energy = 0;
			for (int idx = 0; idx < N; ++idx) {
				kinetic_energy += 0.5 * avg_mass * (all_vx[idx] * all_vx[idx] + all_vy[idx] * all_vy[idx] + all_vz[idx] * all_vz[idx]);
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
		
		void saveToFile(const string& filename) const {
			ofstream file(filename);
			if (!file.is_open()) {
				cerr << "FILE OPENING ERROR: " << filename << endl;
				return;
			}
			
			// Заголовок таблицы
			file << "# index\tx\ty\tz\tvx\tvy\tvz\tmass\n";
			file << fixed << setprecision(6);
			
			int index = 0;
			for (int i = 0; i < sizeX; ++i) {
				for (int j = 0; j < sizeY; ++j) {
					for (int k = 0; k < sizeZ; ++k) {
						const point3D& p = lattice[i][j][k];
						vector<double> pos = p.get_pos();
						vector<double> vel = p.get_vel();
						
						file << index++ << "\t"
							 << pos[0] << "\t" << pos[1] << "\t" << pos[2] << "\t"
							 << vel[0] << "\t" << vel[1] << "\t" << vel[2] << "\t"
							 << p.get_mass() << "\n";
					}
				}
			}
			
			file.close();
			cout << "The data saved to: " << filename << endl;
		}
};

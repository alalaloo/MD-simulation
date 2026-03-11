#include <vector>
#include <cmath>
#include <iostream>
#include <random>

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
				
				for (int i = 0; i < sizeX; ++i) {
					for (int j = 0; j < sizeY; ++j) {
						for (int k = 0; k < sizeZ; ++k) {

								double x = i * distance;
								double y = j * distance;
								double z = k * distance;

								lattice[i][j][k].set_pos(x, y, z);
							}
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

		void set_velocities(double t_target = 298.0, double k_boltzmann = 1.0) {
    			random_device rd; 
    			mt19937 gen(rd()); // генератор Мерсенна Твистера
    			uniform_real_distribution<double> dist(0.0, 1.0);

			double sigma = sqrt(k_boltzmann * t_target / avg_mass);
			
			//vector<double> all_vx, all_vy, all_vz;

			//all_vx.reserve(N);
			//all_vy.reserve(N);
			//all_vz.reserve(N);
			
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
				
						//all_vx.push_back(vx);
						//all_vy.push_back(vy);
						//all_vz.push_back(vz);
			    		}
				}
			}
		}	
				
};


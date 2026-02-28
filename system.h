#include <vector>
#include <cmath>
#include <iostream>

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

	public:
		system() : N(0), sizeX(0), sizeY(0), sizeZ(0) {
			lattice = nullptr;
		}

		system(int n, double distance) {
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

};


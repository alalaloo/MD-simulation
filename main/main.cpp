#include "system.h"
#include <iostream>
#include <iomanip>
#include <chrono>

using namespace std;
using namespace chrono;

int main() {

	const int atoms_count = 10000;
	const double lattice_constant = 1.5874;
	const double mass = 1.0;
	const double temperature = 0.5;
	const double k_boltzmann = 1.0;

	cout << "=== Molecular Dynamics Simulation ===" << endl;
	cout << "Particles: " << atoms_count << endl;
	cout << "Lattice constant: " << lattice_constant << endl;
	cout << "Temperature: " << temperature << endl;
	cout << endl;

	auto total_start = high_resolution_clock::now();

	auto start = high_resolution_clock::now();

	systemMD crystal(atoms_count, lattice_constant, mass, 3.0);

	auto end = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(end - start);
	cout << "[Timing] System initialization: " << duration.count() << " ms" << endl;

	start = high_resolution_clock::now();

	crystal.set_velocities(temperature, k_boltzmann);

	end = high_resolution_clock::now();
	duration = duration_cast<milliseconds>(end - start);
	cout << "[Timing] Velocity initialization: " << duration.count() << " ms" << endl;

	start = high_resolution_clock::now();

	//crystal.compute_full_inter();

	crystal.compute_forces();
	end = high_resolution_clock::now();
	duration = duration_cast<milliseconds>(end - start);
	cout << "[Timing] Force calculation: " << duration.count() << " ms" << endl;
	
	start = high_resolution_clock::now();
	double potential_slow = crystal.get_potential_energy_slow();
	end = high_resolution_clock::now();
	duration = duration_cast<milliseconds>(end - start);
	cout << "[Timing] Potential energy (slow): " << duration.count() << " ms" << endl;

	start = high_resolution_clock::now();
	double potential_fast = crystal.get_potential_energy_fast();
	end = high_resolution_clock::now();
	duration = duration_cast<milliseconds>(end - start);
	cout << "[Timing] Potential energy (fast): " << duration.count() << " ms" << endl;

	// Проверка корректности (разница должна быть близка к нулю)
	cout << "\n=== Energy Comparison ===" << endl;
	cout << "Slow potential: " << setw(15) << potential_slow << endl;
	cout << "Fast potential: " << setw(15) << potential_fast << endl;
	cout << "Difference:     " << setw(15) << (potential_slow - potential_fast) << endl;
    
	start = high_resolution_clock::now();

	double kinetic = crystal.get_kinetic_energy();
	end = high_resolution_clock::now();
	duration = duration_cast<milliseconds>(end - start);
	cout << "[Timing] Kinetic Energy calculation: " << duration.count() << " ms" << endl;

	double potential = potential_slow;
	double total_energy = kinetic + potential;
	
	cout << "\n=== Energy Report ===" << endl;
	cout << "Kinetic energy:   " << setw(15) << kinetic << " (reduced units)" << endl;
	cout << "Potential energy: " << setw(15) << potential << " (reduced units)" << endl;
	cout << "Total energy:     " << setw(15) << total_energy << " (reduced units)" << endl;
	cout << "Temperature:      " << setw(15) << crystal.get_temperature(k_boltzmann) << " (reduced units)" << endl;
	cout << "===================" << endl;

	start = high_resolution_clock::now();
	crystal.saveToFile("lattice_data.csv");
	end = high_resolution_clock::now();
	duration = duration_cast<milliseconds>(end - start);
	cout << "[Timing] File saving: " << duration.count() << " ms" << endl;

	auto total_end = high_resolution_clock::now();
	auto total_duration = duration_cast<milliseconds>(total_end - total_start);
	cout << "\n[Timing] Total execution time: " << total_duration.count() << " ms" << endl;

	return 0;
}

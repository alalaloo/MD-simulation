#include "system.h"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <cstring>
#include <cstdlib>

using namespace std;
using namespace chrono;

void print_usage(const char* prog_name) {
    cout << "Usage: " << prog_name << " [options]\n"
         << "Options:\n"
         << "  -N <int>         Number of particles (approximated to nearest cube) [default: 1000]\n"
         << "  -l <double>      Lattice constant [default: 1.5874]\n"
         << "  -m <double>      Particle mass [default: 1.0]\n"
         << "  -T <double>      Temperature [default: 0.5]\n"
         << "  -k <double>      Boltzmann constant [default: 1.0]\n"
         << "  -c <double>      Cutoff radius [default: 3.0]\n"
         << "  -h               Show this help\n";
}

int main(int argc, char* argv[]) {
    int atoms_count = 1000;
    double lattice_constant = 1.5874;
    double mass = 1.0;
    double temperature = 0.5;
    double k_boltzmann = 1.0;
    double cutoff = 3.0;

    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            print_usage(argv[0]);
            return 0;
        } else if (strcmp(argv[i], "-N") == 0 && i + 1 < argc) {
            atoms_count = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-l") == 0 && i + 1 < argc) {
            lattice_constant = atof(argv[++i]);
        } else if (strcmp(argv[i], "-m") == 0 && i + 1 < argc) {
            mass = atof(argv[++i]);
        } else if (strcmp(argv[i], "-T") == 0 && i + 1 < argc) {
            temperature = atof(argv[++i]);
        } else if (strcmp(argv[i], "-k") == 0 && i + 1 < argc) {
            k_boltzmann = atof(argv[++i]);
        } else if (strcmp(argv[i], "-c") == 0 && i + 1 < argc) {
            cutoff = atof(argv[++i]);
        } else {
            cerr << "Unknown option: " << argv[i] << "\n";
            print_usage(argv[0]);
            return 1;
        }
    }

    cout << "=== Molecular Dynamics Simulation ===" << endl;
    cout << "Requested particles: " << atoms_count << endl;
    cout << "Lattice constant: " << lattice_constant << endl;
    cout << "Mass: " << mass << endl;
    cout << "Temperature: " << temperature << endl;
    cout << "Boltzmann constant: " << k_boltzmann << endl;
    cout << "Cutoff radius: " << cutoff << endl;
    cout << endl;

    auto total_start = high_resolution_clock::now();

    auto start = high_resolution_clock::now();
    systemMD crystal(atoms_count, lattice_constant, mass, cutoff);
    auto end = high_resolution_clock::now();
    cout << "[Timing] System initialization: "
         << duration_cast<milliseconds>(end - start).count() << " ms" << endl;
    cout << "Actual number of particles: " << crystal.getN() << endl;

    start = high_resolution_clock::now();
    crystal.set_velocities(temperature, k_boltzmann);
    end = high_resolution_clock::now();
    cout << "[Timing] Velocity initialization: "
         << duration_cast<milliseconds>(end - start).count() << " ms" << endl;

    // полный перебор (О(N^2))
    start = high_resolution_clock::now();
    crystal.compute_full_inter();
    end = high_resolution_clock::now();
    auto full_force_time = duration_cast<milliseconds>(end - start).count();
    cout << "[Timing] Force calculation (full O(N^2)): " << full_force_time << " ms" << endl;

    vector<double> full_force_0 = crystal.getParticle(0).force;
    vector<double> full_force_1 = crystal.getParticle(1).force;
    vector<double> full_force_2 = crystal.getParticle(2).force;

    // быстрый алгоритм с cell-list
    start = high_resolution_clock::now();
    crystal.compute_forces();
    end = high_resolution_clock::now();
    auto cell_force_time = duration_cast<milliseconds>(end - start).count();
    cout << "[Timing] Force calculation (cell list): " << cell_force_time << " ms" << endl;

    vector<double> cell_force_0 = crystal.getParticle(0).force;
    vector<double> cell_force_1 = crystal.getParticle(1).force;
    vector<double> cell_force_2 = crystal.getParticle(2).force;

    cout << "\n=== Force Comparison (first three particles) ===" << endl;
    cout << "Particle 0 full:   " << full_force_0[0] << " " << full_force_0[1] << " " << full_force_0[2] << endl;
    cout << "Particle 0 cell:   " << cell_force_0[0] << " " << cell_force_0[1] << " " << cell_force_0[2] << endl;
    cout << "Particle 1 full:   " << full_force_1[0] << " " << full_force_1[1] << " " << full_force_1[2] << endl;
    cout << "Particle 1 cell:   " << cell_force_1[0] << " " << cell_force_1[1] << " " << cell_force_1[2] << endl;
    cout << "Particle 2 full:   " << full_force_2[0] << " " << full_force_2[1] << " " << full_force_2[2] << endl;
    cout << "Particle 2 cell:   " << cell_force_2[0] << " " << cell_force_2[1] << " " << cell_force_2[2] << endl;

    // сравнение алгоритмов вычисления энергии
    cout << "\n=== Potential Energy Comparison ===" << endl;
    start = high_resolution_clock::now();
    double pe_slow = crystal.get_potential_energy_slow();
    end = high_resolution_clock::now();
    auto pe_slow_time = duration_cast<milliseconds>(end - start).count();
    cout << "[Timing] Potential energy (full O(N^2)): " << pe_slow_time << " ms" << endl;

    start = high_resolution_clock::now();
    double pe_fast = crystal.get_potential_energy_fast();
    end = high_resolution_clock::now();
    auto pe_fast_time = duration_cast<milliseconds>(end - start).count();
    cout << "[Timing] Potential energy (cell list): " << pe_fast_time << " ms" << endl;

    cout << "Slow (full)   : " << pe_slow << " (reduced units)" << endl;
    cout << "Fast (cell)   : " << pe_fast << " (reduced units)" << endl;
    cout << "Difference    : " << pe_slow - pe_fast << endl;
    if (pe_slow != 0.0)
        cout << "Relative error: " << 100.0 * (pe_slow - pe_fast) / pe_slow << " %" << endl;

    double kinetic = crystal.get_kinetic_energy();
    cout << "\n=== Total Energy ===" << endl;
    cout << "Kinetic energy:   " << setw(15) << kinetic << " (reduced units)" << endl;
    cout << "Potential energy: " << setw(15) << pe_slow << " (reduced units)" << endl;
    cout << "Total energy:     " << setw(15) << kinetic + pe_slow << " (reduced units)" << endl;
    cout << "Temperature:      " << setw(15) << crystal.get_temperature(k_boltzmann) << " (reduced units)" << endl;

    start = high_resolution_clock::now();
    crystal.saveToFile("lattice_data.csv");
    end = high_resolution_clock::now();
    cout << "\n[Timing] File saving: "
         << duration_cast<milliseconds>(end - start).count() << " ms" << endl;

    auto total_end = high_resolution_clock::now();
    cout << "[Timing] Total execution time: "
         << duration_cast<milliseconds>(total_end - total_start).count() << " ms" << endl;

    return 0;
}

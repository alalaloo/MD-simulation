// ──────────────────────────────── main.cpp ────────────────────────────────
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
         << "  -N <int>        Approx. number of particles [default: 4000]\n"
         << "  -l <double>     Lattice constant (FCC unit cell size) [default: 1.5874]\n"
         << "  -m <double>     Particle mass [default: 1.0]\n"
         << "  -T <double>     Temperature [default: 1.0]\n"
         << "  -rc <double>    Cutoff radius [default: 3.0]\n"
         << "  -skin <double>  Verlet list skin [default: 0.3]\n"
         << "  -dt <double>    Time step [default: 0.001]\n"
         << "  -t <double>     Simulation time (if ≤0, only test) [default: -1]\n"
         << "  -out <int>      Output frequency [default: 100]\n"
         << "  -Noz            Enable Nosé-Hoover thermostat\n"
         << "  -Q <double>     Thermostat inertia [default: auto]\n"
         << "  -h              Show this help\n";
}

int main(int argc, char* argv[]) {
    int N_req = 4000;
    double a0 = 1.5874;          // постоянная решётки
    double mass = 1.0;
    double T_req = 1.0;
    double rc = 3.0;
    double skin = 0.3;
    double dt = 0.001;
    double tmax = -1.0;
    int out_freq = 100;
    bool use_nh = false;
    double Q_nh = -1.0;           // отрицательное – автоматический выбор

    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "-h")==0 || strcmp(argv[i], "--help")==0) {
            print_usage(argv[0]); return 0;
        } else if (strcmp(argv[i], "-N")==0 && i+1<argc) {
            N_req = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-l")==0 && i+1<argc) {
            a0 = atof(argv[++i]);
        } else if (strcmp(argv[i], "-m")==0 && i+1<argc) {
            mass = atof(argv[++i]);
        } else if (strcmp(argv[i], "-T")==0 && i+1<argc) {
            T_req = atof(argv[++i]);
        } else if (strcmp(argv[i], "-rc")==0 && i+1<argc) {
            rc = atof(argv[++i]);
        } else if (strcmp(argv[i], "-skin")==0 && i+1<argc) {
            skin = atof(argv[++i]);
        } else if (strcmp(argv[i], "-dt")==0 && i+1<argc) {
            dt = atof(argv[++i]);
        } else if (strcmp(argv[i], "-t")==0 && i+1<argc) {
            tmax = atof(argv[++i]);
        } else if (strcmp(argv[i], "-out")==0 && i+1<argc) {
            out_freq = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-Noz")==0) {
            use_nh = true;
        } else if (strcmp(argv[i], "-Q")==0 && i+1<argc) {
            Q_nh = atof(argv[++i]);
        } else {
            cerr << "Unknown option: " << argv[i] << "\n";
            print_usage(argv[0]); return 1;
        }
    }

    cout << "=== MD heat capacity (Verlet list) ===\n";
    cout << "N ~ " << N_req << ", a0 = " << a0 << ", T = " << T_req
         << ", rc = " << rc << ", skin = " << skin << endl;

    systemMD sys(N_req, a0, mass, rc, skin);
    cout << "Actual number of particles: " << sys.getN() << endl;
    sys.set_velocities(T_req, 1.0);

    // Полный перебор сил
    auto t1 = high_resolution_clock::now();
    sys.compute_full_inter();
    auto t2 = high_resolution_clock::now();
    auto full_fx = sys.get_force_x(), full_fy = sys.get_force_y(), full_fz = sys.get_force_z();

    // Verlet-список
    sys.build_verlet_list();
    auto t3 = high_resolution_clock::now();
    sys.compute_forces();
    auto t4 = high_resolution_clock::now();
    auto fast_fx = sys.get_force_x(), fast_fy = sys.get_force_y(), fast_fz = sys.get_force_z();

    cout << "Force check (particle 0): full (" << full_fx[0] << "," << full_fy[0] << "," << full_fz[0] << ")\n";
    cout << "Force check (particle 0): fast (" << fast_fx[0] << "," << fast_fy[0] << "," << fast_fz[0] << ")\n";

    double U_full = sys.get_potential_energy_slow();
    double U_fast = sys.get_potential_energy_fast();
    cout << "Potential energy: full = " << U_full << ", fast = " << U_fast
         << ", diff = " << U_full - U_fast << endl;

    cout << "\nTiming:\n";
    cout << "Full O(N^2) force: " << duration_cast<milliseconds>(t2 - t1).count() << " ms\n";
    cout << "Verlet list force: " << duration_cast<milliseconds>(t4 - t3).count() << " ms\n";

    if (tmax > 0.0) {
        if (use_nh) {
            if (Q_nh <= 0.0) {
                double g = 3.0 * sys.getN();
                Q_nh = g * T_req * (0.1 * 0.1);  // τ_T ≈ 0.1
            }
            sys.enable_nose_hoover(T_req, Q_nh);
        }
        cout << "\nStarting Verlet integration...\n";
        auto t0 = high_resolution_clock::now();
        sys.run_verlet(dt, tmax, out_freq, "trajectory.csv");
        auto t5 = high_resolution_clock::now();
        cout << "Integration time: " << duration_cast<milliseconds>(t5 - t0).count() << " ms\n";
    }

    return 0;
}

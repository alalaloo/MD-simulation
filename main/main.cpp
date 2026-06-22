#include "../libs/system.h"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <cstring>
#include <cstdlib>

using namespace std;
using namespace chrono;

void print_usage(const char* prog_name) {
    cout << "Usage: " << prog_name << " [options]\n"
         << "Options (All parameters are strictly in SI units where applicable):\n"
         << "  -N <int>        Approx. number of particles [default: 4000]\n"
         << "  -l <double>     Lattice constant a0 (m) [default: 5.256e-10 (Argon)]\n"
         << "  -sigma <double> LJ sigma parameter (m) [default: 3.405e-10]\n"
         << "  -eps <double>   LJ epsilon parameter (J) [default: 1.654e-21]\n"
         << "  -m <double>     Particle mass (kg) [default: 6.634e-26]\n"
         << "  -T <double>     Target Temperature (K) [default: 87.3]\n"
         << "  -rc <double>    Cutoff radius multiplier of sigma (rc = rc_mult * sigma) [default: 1.9435 (3*sigma)]\n"
         << "  -skin <double>  Verlet list skin fraction of sigma (skin = skin_mult * sigma) [default: 0.19435 (0.3*sigma)]\n"
         << "  -dt <double>    Time step (s) [default: 2.0e-15 (2 fs)]\n"
         << "  -steps <int>    Total simulation steps (if <=0, only test) [default: -1]\n"
         << "  -out <int>      Output frequency [default: 100]\n"
         << "  -Noz            Enable Nosé-Hoover Chains thermostat (NVT ensemble)\n"
         << "  -tau <double>   Thermostat relaxation time multiplier of dt (tau = tau_mult * dt) [default: 50.0]\n"
         << "  -tequil <double> Equilibration time multiplier of dt (tequil = tequil_mult * dt) [default: 1500.0]\n"
         << "  -h              Show this help\n";
}

int main(int argc, char* argv[]) {
    // Дефолтные значения физических параметров в СИ для Аргона (Ar)
    int N_req = 4000;
    double a0 = 5.256e-10;          
    double sigma = 3.405e-10;
    double epsilon = 1.654e-21;
    double mass = 6.634e-26;
    double T_req = 87.3; // Температура кипения Аргона
    double dt = 2.0e-15; // 2 фс
    
    int num_steps = -1; // Вместо tmax задаем количество шагов
    int out_freq = 100;
    bool use_nh = false;
    
    // Множители по умолчанию (соответствуют старым дефолтным значениям)
    double rc_factor = 1.9435;      // rc / a0 = (3 * 3.405e-10) / 5.256e-10
    double skin_factor = 0.19435;   // skin / a0 = (0.3 * 3.405e-10) / 5.256e-10
    double tau_factor = 50.0;       // tau / dt = 1.0e-13 / 2.0e-15
    double tequil_factor = 1500.0;  // tequil / dt = 3.0e-12 / 2.0e-15

    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            print_usage(argv[0]); return 0;
        } else if (strcmp(argv[i], "-N") == 0 && i + 1 < argc) {
            N_req = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-l") == 0 && i + 1 < argc) {
            a0 = atof(argv[++i]);
        } else if (strcmp(argv[i], "-sigma") == 0 && i + 1 < argc) {
            sigma = atof(argv[++i]);
        } else if (strcmp(argv[i], "-eps") == 0 && i + 1 < argc) {
            epsilon = atof(argv[++i]);
        } else if (strcmp(argv[i], "-m") == 0 && i + 1 < argc) {
            mass = atof(argv[++i]);
        } else if (strcmp(argv[i], "-T") == 0 && i + 1 < argc) {
            T_req = atof(argv[++i]);
        } else if (strcmp(argv[i], "-rc") == 0 && i + 1 < argc) {
            rc_factor = atof(argv[++i]); // Теперь считываем как множитель
        } else if (strcmp(argv[i], "-skin") == 0 && i + 1 < argc) {
            skin_factor = atof(argv[++i]); // Теперь считываем как долю/множитель
        } else if (strcmp(argv[i], "-dt") == 0 && i + 1 < argc) {
            dt = atof(argv[++i]);
        } else if (strcmp(argv[i], "-steps") == 0 && i + 1 < argc) {
            num_steps = atoi(argv[++i]); // Считываем количество шагов
        } else if (strcmp(argv[i], "-out") == 0 && i + 1 < argc) {
            out_freq = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-Noz") == 0) {
            use_nh = true;
        } else if (strcmp(argv[i], "-tau") == 0 && i + 1 < argc) {
            tau_factor = atof(argv[++i]); // Считываем как множитель для dt
        } else if (strcmp(argv[i], "-tequil") == 0 && i + 1 < argc) {
            tequil_factor = atof(argv[++i]); // Считываем как множитель для dt
        } else {
            cerr << "Unknown option: " << argv[i] << "\n";
            print_usage(argv[0]); return 1;
        }
    }

    // Вычисляем абсолютные физические величины на основе множителей
    double rc = rc_factor * sigma;
    double skin = skin_factor * sigma;
    double tmax = num_steps * dt;
    double tau_nh = tau_factor * dt;
    double t_equil = tequil_factor * dt;

    cout << "=== MD Simulation (Strict SI Units Pipeline) ===\n";
    cout << "Parameters Set:\n"
         << "  N_req   = " << N_req << "\n"
         << "  a0      = " << a0 << " m\n"
         << "  sigma   = " << sigma << " m\n"
         << "  epsilon = " << epsilon << " J\n"
         << "  mass    = " << mass << " kg\n"
         << "  T_req   = " << T_req << " K\n"
         << "  rc      = " << rc << " m (multiplier of sigma: " << rc_factor << ")\n"
         << "  skin    = " << skin << " m (fraction of sigma: " << skin_factor << ")\n"
         << "  dt      = " << dt << " s\n"
         << "  steps   = " << num_steps << " (calculated total time: " << tmax << " s)\n";

    systemMD sys(N_req, a0, sigma, epsilon, mass, rc, skin);
    cout << "Actual number of particles (FCC): " << sys.getN() << "\n";
    
    // Передаем реальную константу Больцмана СИ: 1.380649e-23 Дж/К

    sys.set_velocities(T_req, 1.380649e-23);

    auto t1 = high_resolution_clock::now();
    sys.compute_full_inter();
    auto t2 = high_resolution_clock::now();
    auto full_fx = sys.get_force_x(), full_fy = sys.get_force_y(), full_fz = sys.get_force_z();

    sys.build_verlet_list();
    auto t3 = high_resolution_clock::now();
    sys.compute_forces(); 
    auto t4 = high_resolution_clock::now();
    
#ifdef USE_CUDA
    sys.sync_device_to_host();
#endif
    auto fast_fx = sys.get_force_x(), fast_fy = sys.get_force_y(), fast_fz = sys.get_force_z();

    cout << "Force validation (particle 0) [N]:\n";
    cout << "  -> Full O(N^2) : (" << full_fx[0] << ", " << full_fy[0] << ", " << full_fz[0] << ")\n";
    cout << "  -> Optimized   : (" << fast_fx[0] << ", " << fast_fy[0] << ", " << fast_fz[0] << ")\n";

    cout << "Time spent full O(N^2): " << duration_cast<milliseconds>(t2 - t1).count() << " ms\n";
    cout << "Time spent optimized list: " << duration_cast<milliseconds>(t4 - t3).count() << " ms\n";

    if (num_steps > 0) {
        if (use_nh) {
            sys.enable_nose_hoover(T_req, tau_nh);
            cout << "Nosé-Hoover Chains Thermostat ENABLED (tau = " << tau_nh << " s, multiplier = " << tau_factor << ").\n";
        }
        cout << "Equilibration time set to: " << t_equil << " s (multiplier = " << tequil_factor << ")\n";
        cout << "\nStarting production integration pipeline...\n";
        
        auto t0 = high_resolution_clock::now();
        sys.run_verlet(dt, tmax, out_freq, "trajectory.csv", t_equil);
        auto t5 = high_resolution_clock::now();
        
        cout << "Integration completely finished. Total pipeline time: " << duration_cast<milliseconds>(t5 - t0).count() << " ms\n";
    }

    return 0;
}

#include "system.h"

int main() {

    const int atoms_count = 125;      
    const double lattice_constant = 1.5874;  
    const double mass = 1.0;          
    const double temperature = 0.5;    
    const double k_boltzmann = 1.0;    
    
    systemMD crystal(atoms_count, lattice_constant, mass);
    
    crystal.set_velocities(temperature, k_boltzmann);

    crystal.saveToFile("lattice_data.csv");
    
    return 0;
}

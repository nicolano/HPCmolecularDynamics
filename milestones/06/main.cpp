//
// Created by Nicolas von Trott on 25.05.22.
//

#include "neighbors.h"
#include "set_cubic_lattice.h"
#include "berendsen_thermostat.h"
#include "kin_energy.h"
#include "lj.h"
#include "verlet.h"
#include "xyz.h"
#include <iostream>
#include <cmath>
#include <ctime>

int main() {
    // Declare simulation parameters.
    const double epsilon = 1.0;
    const double sigma = 1.0;
    const double mass = 1.0;

    // Declare physical constants.
    const double k_b = 1.380649e-23;

    // Set the cutoff range for the neighbor list.
    const double cutoff_range = 1.8 * sigma;

    // Target temperature for the berendsen thermostat.
    const double temperature = 10;
    // Set relaxation time for the berendsen thermostat.
    const double relaxation_time = 0.1;

    // Set the maximal simulation time.
    const double t_max = 1000 * sqrt(mass * sigma * sigma / epsilon);
    // Set the time step.
    const double delta_t = 0.02 * sqrt(mass * sigma * sigma / epsilon);
    // Calculate the number of time steps to reach the max simulation time.
    const int nb_steps = int (t_max / delta_t);

    // Set number of atoms for which the simulation should be performed.
    const int nb_atoms_arr[13] = {5,25,50,75,100,125,150,175,200,225,250,275,300};
    // Set the lattice constant for the cubic atom lattice.
    const double lattice_constant = 1 * sigma;

    for (int nb_atoms : nb_atoms_arr) {
        // Get starting time.
        clock_t t_start = clock();

        // Create Atoms structure for the cluster.
        Atoms atoms(nb_atoms);
        // Place atoms in a cubic grid.
        set_cubic_lattice(atoms, lattice_constant);
        // Create NeighborList structure.
        NeighborList neighbor_list(cutoff_range);

        std::ofstream traj("/Users/nicolasvontrott/Desktop/molDyn/milestones/06/result.xyz");
        write_xyz(traj, atoms);

        // Perform the simulation for the cluster.
        for (int i = 0; i < nb_steps; ++i) {
            // Update neighbor list.
            neighbor_list.update(atoms);

            // Perform first verlet step
            verlet_step1( atoms, delta_t, mass);

            // Calculate LJ potential with cutoff and get potential energy.
            double e_pot{lj(atoms, neighbor_list, cutoff_range, epsilon, sigma)};
            // Calculate kinetic energy from atoms velocities and masses.
            double e_kin{kin_energy(atoms, mass)};
            // Calculate total energy
            double e_tot{e_kin + e_pot};

            // Perform second verlet step.
            verlet_step2( atoms, delta_t, mass);

            if (i > 0){
                // Apply the berendsen thermostat.
                berendsen_thermostat(atoms, e_kin, temperature, delta_t, relaxation_time, k_b);
            }

            // Write positions to file after a number of  timesteps.
            if (i % 100 == 0) {
                // Calculate current time.
                //double t{i * delta_t};

                // Print current time.
                //std::cout << t << " of " << t_max << ": ";

                // Print current energies.
                //std::cout << "e_pot = " << e_pot << ", ";
                //std::cout << "e_kin = " << e_kin << ", ";
                //std::cout << "e_tot = " << e_tot << '\n';

                // Write trajectory to file.
                write_xyz(traj, atoms);
            }
        }

        // Get finish time.
        clock_t t_stop = clock();
        // Calculate elapsed time between start and finish.
        double t_elapsed = (double) (t_stop - t_start) / CLOCKS_PER_SEC;
        // Print result.
        std::cout << "nb_atoms_arr = " << nb_atoms << ", ";
        std::cout << "t_elapsed time = " << t_elapsed << '\n';

        // Close file streams.
        traj.close();
    }

    return 0;
}

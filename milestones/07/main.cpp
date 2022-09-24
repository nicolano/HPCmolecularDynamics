//
// Created by Nicolas von Trott on 25.05.22.
//

#include "gupta.h"
#include "neighbors.h"
#include "berendsen_thermostat.h"
#include "kin_energy.h"
#include "verlet.h"
#include "xyz.h"
#include <iostream>
#include <cmath>

int main() {
    // Declare physical constants.
    const double mass = 196.96657 / 0.009649,
                 k_b = 8.617333262 * pow(10,-5),
                 A = 0.2061,
                 xi = 1.790,
                 p = 10.229,
                 q = 4.036,
                 re = 4.079 / sqrt(2);

    // Set the cutoff range for the neighbor list.
    const double cutoff_range = 10;

    // Set relaxation time for the berendsen thermostat.
    const double relaxation_time = 1000;

    // Set the maximal simulation time.
    const double t_max = 2 * relaxation_time;
    // Set the time step.
    const double delta_t = 1.;

    // Read names, positions and velocities from file.
    auto [names, positions, velocities]{read_xyz_with_velocities("/Users/nicolasvontrott/Desktop/molDyn/milestones/07/cluster_923.xyz")};
    // Store names, positions and velocities in atoms structure.
    Atoms atoms(positions);

    // Create NeighborList structure.
    NeighborList neighbor_list(cutoff_range);

    // Open stream for the trajectory.
    std::ofstream traj("/Users/nicolasvontrott/Desktop/molDyn/milestones/07/traj_923.xyz");
    // Open stream for the energy.
    std::ofstream ener("/Users/nicolasvontrott/Desktop/molDyn/milestones/07/epot_923.xyz");

    // Write initial positions to file.
    write_xyz(traj, atoms);

    // Heat the cluster up by raising the temperature multiple times.
    for (int j = 0; j < 40; j++) {
        // Set the target temperature for the current run.
        double temperature = j * 50;

        // Declare variables for mean potential and kinetic energy.
        double mean_pot = 0,
               mean_kin = 0;

        // Perform the simulation for the current target temperature.
        for (int t = 0; t <= t_max; t += int(delta_t)) {
            // Update neighbor list
            neighbor_list.update(atoms);

            // Perform first verlet step.
            verlet_step1(atoms, delta_t, mass);

            // Calculate gupta potential and get potential energy.
            double e_pot{gupta(atoms, neighbor_list, cutoff_range, A, xi, p, q, re)};
            // Calculate kinetic energy from atoms velocities and masses.
            double e_kin{kin_energy(atoms, mass)};

            // Perform second verlet step.
            verlet_step2(atoms, delta_t, mass);

            if (t > 0) {
                // Apply the berendsen thermostat.
                berendsen_thermostat(atoms, e_kin, temperature, delta_t, relaxation_time, k_b);
            }

            // Write positions to file after x timesteps.
            if (t % 100 == 0) {
                write_xyz(traj, atoms);
                //std::cout << t<< "\t\t";
                //std::cout << "e_pot =";
                //std::cout << e_pot << "\t\t";
                //std::cout << "e_kin =";
                //std::cout << e_kin << "\t\t";
                //std::cout << "e_tot =";
                //std::cout << e_pot + e_kin << '\n';
            }

            // Wait after the relaxation time is over.
            if (t >= relaxation_time) {
                // Sum the kinetic and potential energy nb_steps / 2 steps.
                mean_kin += e_kin;
                mean_pot += e_pot;
            }
        }

        // Calculate the mean potential and kinetic energy.
        mean_pot = mean_pot / (relaxation_time / 2);
        mean_kin = mean_kin / (relaxation_time / 2);
        // Calculate the mean total energy.
        double mean_tot{mean_pot + mean_kin};

        // Write mean total energy and respective temperature to file.
        write_energy(ener, mean_tot, temperature);

        // Print results.
        std::cout << temperature;
        std::cout << "\t\t";
        std::cout << "e_pot =";
        std::cout << mean_pot << "\t\t";
        std::cout << "e_kin =";
        std::cout << mean_kin << "\t\t";
        std::cout << "e_tot =";
        std::cout << mean_pot + mean_kin << '\n';
    }

    // Close streams.
    traj.close();
    ener.close();

    return 0;
}

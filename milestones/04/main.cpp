#include "kin_energy.h"
#include "lj_direct_summation.h"
#include "verlet.h"
#include "xyz.h"
#include <iostream>
#include <cmath>

int main() {
    // Read names, positions and velocities from file.
    auto [names, positions, velocities]{read_xyz_with_velocities("/Users/nicolasvontrott/Desktop/molDyn/milestones/04/lj54.xyz")};
    // Store names, positions and velocities in atoms structure
    Atoms atoms(positions, velocities);

    // Declare simulation parameters
    constexpr double epsilon = 1.0;
    constexpr double sigma = 1.0;
    constexpr double mass = 1.0;

    // Declare Variable to store the total energy.
    double e_tot = 0;

    // Set the maximal simulation time.
    double t_max = 100 * sqrt(mass * sigma * sigma / epsilon);
    // Set the time step.
    double delta_t = 0.02 * sqrt(mass * sigma * sigma / epsilon);
    // Calculate the number of timesteps to reach the max simulation time.
    int nb_steps = int (t_max / delta_t);

    // Open stream for writing the trajectory of the atoms to a file.
    std::ofstream traj("/Users/nicolasvontrott/Desktop/molDyn/milestones/04/result.xyz");
    // Write initial positions to file.
    write_xyz(traj, atoms);

    // Open stream for writing the total energy e_tot and time t to a file.
    std::ofstream energy("/Users/nicolasvontrott/Desktop/molDyn/milestones/04/t001.txt");

    // Perform simulation.
    for (int i = 0; i < nb_steps; ++i) {
        // Perform first verlet step.
        verlet_step1( atoms, delta_t, mass);

        // Do LJ direct summation and get potential energy.
        double e_pot{lj_direct_summation(atoms, epsilon, sigma)};
        // Calculate kinetic Energy from atoms velocities and masses.
        double e_kin{kin_energy(atoms, mass)};
        // Calculate total energy
        e_tot = e_kin + e_pot;

        // Perform second verlet step.
        verlet_step2( atoms, delta_t, mass);

        // Write positions to file after 10 timesteps.
        if (i % 10 == 0) {
            // Calculate current time.
            double t{i * delta_t};

            // Print current time.
            std::cout << t << " of " << t_max << ": ";

            // Print current energies.
            std::cout << "e_pot = " << e_pot << ", ";
            std::cout << "e_kin = " << e_kin << ", ";
            std::cout << "e_tot = " << e_tot << '\n';

            // Write trajectory to file.
            write_xyz(traj, atoms);

            // Write time and energy to file.
            write_energy(energy, t, e_tot);
        }
    }

    // Close file streams.
    traj.close();
    energy.close();

    return 0;
}

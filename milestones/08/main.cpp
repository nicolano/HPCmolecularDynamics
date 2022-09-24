//
// Created by Nicolas von Trott on 25.05.22.
//

#include "domain.h"
#include "mpi.h"
#include "gupta.h"
#include "neighbors.h"
#include "berendsen_thermostat.h"
#include "kin_energy.h"
#include "verlet.h"
#include "xyz.h"
#include <iostream>
#include <cmath>

int main(int argc, char *argv[]) {
    // Initialize Communicator.
    MPI_Init(&argc, &argv);

    // Get rank local rank and cluster size from Communicator.
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Set up domain.
    Domain domain(
        MPI_COMM_WORLD,
        {35.0, 35.0, 35.0},
        {1, 1, size },
        {1, 1, 1});

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

    // Set the border width.
    const double border_width = 2 * cutoff_range;

    // Set relaxation time for the berendsen thermostat.
    const double relaxation_time = 1000;

    // Set the maximal simulation time.
    const double t_max = 2 * relaxation_time;
    // Set the time step.
    const double delta_t = 1.;

    // Read names, positions and velocities from file.
    auto [names, positions, velocities]{read_xyz_with_velocities("/Users/nicolasvontrott/Desktop/molDyn/milestones/08/cluster_923.xyz")};
    // Store names, positions and velocities in atoms structure.
    Atoms atoms(positions);
    // Set the atoms masses.
    atoms.masses.setConstant(mass);

    // Create NeighborList structure.
    NeighborList neighbor_list(cutoff_range);

    // Open stream for the trajectory.
    std::ofstream traj("/Users/nicolasvontrott/Desktop/molDyn/milestones/08/traj_n2_923.xyz");
    // Open stream for the energy.
    std::ofstream ener("/Users/nicolasvontrott/Desktop/molDyn/milestones/08/epot_n2_923.xyz");

    // Write initial positions to file.
    write_xyz(traj, atoms);

    // Switch into the decomposed state.
    domain.enable(atoms);

    // Perform the simulation for the current target temperature.
    for (int t = 0; t <= t_max; t += int(delta_t)) {
        // Perform first verlet step.
        verlet_step1(atoms, delta_t, mass);

        // Move the current subdomain leaving atoms to the neighboring domain.
        domain.exchange_atoms(atoms);

        // Populate the Atoms object with ghost atoms from the neighboring subdomain
        domain.update_ghosts(atoms, border_width);

        // Update neighbor list.
        neighbor_list.update(atoms);

        // Calculate gupta potential for local domain and gather results for the potential energy.
        double epot_global{MPI::allreduce(gupta_sub(atoms, neighbor_list, domain.nb_local(),cutoff_range, A, xi, p, q, re),
                           MPI_SUM,
                           MPI_COMM_WORLD)};
        // Calculate kinetic energy for local domain and gather results from atoms velocities and masses.
        double ekin_global{MPI::allreduce(kin_energy_sub(atoms, mass, domain.nb_local()),
                           MPI_SUM,
                           MPI_COMM_WORLD)};

        // Perform second verlet step.
        verlet_step2(atoms, delta_t, mass);

        // Write positions to file after x timesteps.
        if (t % 10 == 0) {
            // Switch back to the replicated state for writing trajectories
            domain.disable(atoms);

            if (rank == 0) {
                double e_tot_global{epot_global + ekin_global};
                double time{double (t)};
                write_xyz(traj, atoms);
                write_energy(ener, e_tot_global, time);

                std::cout << t << "\t\t";
                std::cout << "e_pot =";
                std::cout << epot_global << "\t\t";
                std::cout << "e_kin =";
                std::cout << ekin_global << "\t\t";
                std::cout << "e_tot =";
                std::cout << epot_global + ekin_global << "\n";
                // std::cout << domain.nb_local() << "\n";
                // std::cout << atoms.nb_atoms() - domain.nb_local() << "\n";
            }

            // Switch back to the decomposed state
            domain.enable(atoms);
        }
    }

    // Close file streams.
    traj.close();
    ener.close();

    // Shut down cluster.
    MPI_Finalize();

    return 0;
}

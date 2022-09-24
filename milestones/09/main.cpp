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
#include "pressure.h"
#include "force_on_whisker.h"
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
    Eigen::Array3d domain_length = {40, 40, 144.0};
    Domain domain(
        MPI_COMM_WORLD,
        {domain_length(0), domain_length(1), domain_length(2)},
        {1, 1, size },
        {0, 0, 1});

    // Declare physical constants.
    const double mass = 196.96657 / 0.009649,
                 k_b = 8.617333262 * pow(10,-5),
                 A = 0.2061,
                 xi = 1.790,
                 p = 10.229,
                 q = 4.036,
                 re = 4.079 / sqrt(2);

    // Set the target temperature for the simulation.
    double temperature = 0;

    // Set the cutoff range for the neighbor list.
    const double cutoff_range = 10;

    // Set the border width.
    const double border_width = 2 * cutoff_range;

    // Set relaxation time for the berendsen thermostat.
    const double relaxation_time = 1000;

    // Set the maximal simulation time.
    const double t_max = 6000 * relaxation_time;
    // Set the time step.
    const double delta_t = 1;

    // Read names, positions and velocities from file.
    auto [names, positions, velocities]{read_xyz_with_velocities("/Users/nicolasvontrott/Desktop/molDyn/milestones/09/whisker_small.xyz")};
    // Store names, positions and velocities in atoms structure.
    Atoms atoms(positions);
    // Set the atoms masses.
    atoms.masses.setConstant(mass);
    int numberOfAtoms = atoms.nb_atoms();

    double length_0 = domain.domain_length(2);

    // Create NeighborList structure.
    NeighborList neighbor_list(cutoff_range);

    // Open stream for the trajectory.
    std::ofstream traj("/Users/nicolasvontrott/Desktop/molDyn/milestones/09/traj_small_s000014.xyz");
    // Open stream for the force.
    std::ofstream force_l("/Users/nicolasvontrott/Desktop/molDyn/milestones/09/force_small_whisker_s000014.xyz");
    // Open stream for the force.
    std::ofstream force_r("/Users/nicolasvontrott/Desktop/molDyn/milestones/09/force_small_whisker_s000014.xyz");

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

         if (t > 0) {
            // Apply the berendsen thermostat.
            berendsen_thermostat(atoms, ekin_global, temperature, delta_t, relaxation_time, k_b);
         }

         if ((t > relaxation_time) && (t % 1000 == 0) || (t == 0)) {
             // Calculate the force on the left side of the whisker
            if (rank == 0) {
                double strain{(domain.domain_length(2) - length_0 ) / length_0};
                double force_left{force_on_whisker(atoms, domain.nb_local(),0, true)};
                write_energy(force_l, strain, force_left);
                std::cout << "Strain: ";
                std::cout << strain << "\t\t";

                std::cout << "Force left: ";
                std::cout << force_left << "\t\t";
            } else if (rank == size - 1) {
                // Calculate the force on the right side of the whisker
                double strain{(domain.domain_length(2) - length_0 ) / length_0};
                double force_right{force_on_whisker(atoms, domain.nb_local(),domain.domain_length(2), true)};
                write_energy(force_r, strain, force_right);
                std::cout << "Force right: ";
                std::cout << force_right << "\n";
            }

            // Switch back to the replicated state for writing trajectories
            domain.disable(atoms);

            if (rank == 0) {
                write_xyz(traj, atoms);
            }

            // Switch back to the decomposed state
            domain.enable(atoms);
        }

        // Stretch the whisker.
        if (t > relaxation_time) {
            domain_length(2) += 0.000014;
            domain.scale(atoms, domain_length);
        }
    }

    // Close file streams.
    traj.close();
    force_l.close();
    force_r.close();

    // Shut down cluster.
    MPI_Finalize();

    return 0;
}

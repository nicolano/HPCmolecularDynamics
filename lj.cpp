//
// Created by Nicolas von Trott on 01.06.22.
//

#include "lj.h"
#include "neighbors.h"

double lj(Atoms &atoms, const NeighborList& neighbor_list, double cutoff_range, double epsilon, double sigma) {
    // Variable to store the sum of the potential energy of all atoms.
    double e_pot = 0.0;
    // Calc sigma to the power of 6 and 12
    double sigma6 = pow(sigma, 6);
    double sigma12 = pow(sigma6, 2);
    // Calc cutoff to the power of 6, 8, 12, and 14.
    double cutoff6 = pow(cutoff_range, 6);
    double cutoff8 = pow(cutoff_range, 8);
    double cutoff12 = pow(cutoff_range, 12);
    double cutoff14 = pow(cutoff_range, 14);

    // Set all forces to zero
    atoms.forces.setZero();

    // Loop over all atom pairs in the neighbor list
    // (e.g. atoms within the cutoff range).
    for (auto[i, j] : neighbor_list) {
        // Make sure to count atom pairs only one time.
        if (i < j) {
            // calculate cartesian distance vector
            Eigen::Vector3d d_vec_ij {
                atoms.positions.col(i) - atoms.positions.col(j)
            };

            // Calculate vector norm resp. distance between the two atoms.
            double r {d_vec_ij.norm()};
            // Calculate the vector norm to the power of 2, 6, and 12.
            double r2 {pow(r, 2)};
            double r6 {pow(r2, 3)};
            double r12 {pow(r6, 2)};

            // Calculate potential energy from the lj potential.
            e_pot += 4 * epsilon * (sigma12 / r12 - sigma6 / r6);
            // Calculate energy at cutoff distance.
            double offset_e{4 * epsilon * (sigma12 / cutoff12 - sigma6 / cutoff6)};
            // Shift energy to zero at cutoff distance.
            e_pot -= offset_e;
                
            // Calculate the absolute force acting between atom i and j.
            double f_abs_ik {4 * epsilon * (12 * sigma12 / (r12 * r2) - 6 * sigma6 / (r6 * r2))};
            // Calculate force at cutoff distance.
            double offset_f {4 * epsilon * (12 * sigma12 / cutoff14 - 6 * sigma6 / cutoff8)};
            // Calculate the force vector between atom i and j.
            // Shift force to zero at cutoff distance.
            Eigen::Array3d force_vec {(f_abs_ik - offset_f) * d_vec_ij};

            // Add calculated force vector two force of atom i.
            atoms.forces.col(i) += force_vec;
            // Distract calculated force vector from forces of atom j.
            atoms.forces.col(j) -= force_vec;
        }
    }

    // Return the sum of the potential energy between all atoms in the neighbor list.
    return e_pot;
}
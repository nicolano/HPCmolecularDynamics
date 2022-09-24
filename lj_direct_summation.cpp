#include "lj_direct_summation.h"

double lj_direct_summation(Atoms &atoms, double epsilon, double sigma) {
    // Variable to store the sum of the potential energy of all atoms.
    double e_pot = 0.0;
    // Get the number of atoms in the Atom structure.
    long nb_atoms = atoms.nb_atoms();
    // Calc the value of sigma to the power of 6 and 12.
    double sigma6 = pow(sigma, 6);
    double sigma12 = pow(sigma6, 2);

    // Set all forces to zero.
    atoms.forces.setZero();

    // Loop over all pairs of atoms in the structure.
    for (long i = 0; i < nb_atoms - 1; i++) {
        for (long k = i + 1; k < nb_atoms; k++) {
            // Calculate cartesian distance vector between atom i and j.
            Eigen::Vector3d d_vec_ik{
                atoms.positions.col(i) - atoms.positions.col(k)
            };

            // Calculate vector norm resp. distance between the two atoms.
            double r {d_vec_ik.norm()};
            // Calculate the vector norm to the power of 2, 6, and 12.
            double r2 {pow(r, 2)};
            double r6 {pow(r2, 3)};
            double r12 {pow(r6, 2)};

            // Calculate potential energy from the lj potential.
            e_pot += 4 * epsilon * (sigma12 / r12 - sigma6 / r6);

            // Calculate the absolute force acting between atom i and j.
            double f_abs_ik{4 * epsilon * (12 * sigma12 / (r12 * r2) - 6 * sigma6 / (r6 * r2))};
            // Calculate the force vector between atom i and j.
            Eigen::Array3d f_vec_ik{f_abs_ik * d_vec_ik};

            // Add calculated force vector two force of atom i.
            atoms.forces.col(i) += f_vec_ik;
            // Distract calculated force vector from forces of atom j.
            atoms.forces.col(k) -= f_vec_ik;
        }
    }

    // Return the sum of the potential energy between all atoms in the atom structure.
    return e_pot;
}

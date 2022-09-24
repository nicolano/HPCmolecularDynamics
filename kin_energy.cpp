#include "kin_energy.h"

double kin_energy(Atoms &atoms, double mass) {
    // Declare variable to store the kinetic energy of all atoms.
    double e_kin = 0;

    // Loop over all atoms in the Atoms structure.
    for (int i = 0; i < atoms.nb_atoms(); ++i) {
        // Get the velocity in vector form from the atom with index i.
        Eigen::Vector3d velocity {atoms.velocities.col(i)};
        // Calculate the kinetic energy of atom i.
        double e_kin_i{0.5 * mass * velocity.squaredNorm()};
        // Sum the kinetic energy of atom i to the total kinetic energy.
        e_kin += e_kin_i;
    }

    // Return the sum of the kinetic energy.
    return e_kin;
}

double kin_energy_sub(Atoms &atoms, double mass, long nb_local) {
    // Declare variable to store the kinetic energy of all atoms in the subset.
    double e_kin = 0;

    // Loop over all atoms in the Atoms subset.
    for (int i = 0; i < nb_local; ++i) {
        // Get the velocity in vector form from the atom with index i.
        Eigen::Vector3d velocity {atoms.velocities.col(i)};
        // Calculate the kinetic energy of atom i.
        double e_kin_i{0.5 * mass * velocity.squaredNorm()};
        // Sum the kinetic energy of atom i to the total kinetic energy.
        e_kin += e_kin_i;
    }

    // Return the sum of the kinetic energy.
    return e_kin;
}
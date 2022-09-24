//
// Created by Nicolas von Trott on 01.06.22.
//
#include "set_cubic_lattice.h"

void set_cubic_lattice(Atoms &atoms, double lattice_constant) {
    // Get number of atoms in the Atoms structure.
    long nb_atoms = atoms.nb_atoms();
    // Get the size of the next smallest cubic lattice in which all atoms will fit.
    // Not every grid point must be occupied by an atom.
    int grid_size = ceil(cbrt(nb_atoms) / lattice_constant);

    // Variable to count the placed atoms.
    int i = 0;
    // Variable for breaking the loops when all atoms are positioned.
    int flag = 1;
    // Loop over the cartesian coordinates x, y, and z.
    for (int x = 0; x < grid_size; ++x) {
        for (int y = 0; y < grid_size; ++y) {
            for (int z = 0; z < grid_size; ++z) {

                // Break loop if all atoms are positioned.
                if (i == nb_atoms) {
                    flag = 0;
                    break;
                }

                // Get and set the x, y, and z coordinates of the atom.
                atoms.positions(0, i) = x * lattice_constant;
                atoms.positions(1, i) = y * lattice_constant;
                atoms.positions(2, i) = z * lattice_constant;

                // Count the placed atom.
                i++;
            }

            // Break loop if all atoms are positioned.
            if (flag == 0) {
                break;
            }
        }

        // Break loop if all atoms are positioned.
        if (flag == 0) {
            break;
        }
    }
}

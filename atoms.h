#ifndef ATOMS_H
#define ATOMS_H

#include "types.h"
#include <cfloat>
#include <iostream>
class Atoms { 
public: 
    Positions_t positions; 
    Velocities_t velocities; 
    Forces_t forces;
    Masses_t masses;
 
    Atoms(const int &nb_atoms) : 
            positions{3, nb_atoms}, velocities{3, nb_atoms}, forces{3, nb_atoms}, masses{nb_atoms}{
        positions.setZero();
        velocities.setZero(); 
        forces.setZero();
        masses.setOnes();
    }

    Atoms(const Positions_t &p) : 
            positions{p}, velocities{3, p.cols()}, forces{3, p.cols()}, masses{p.cols()} {
        velocities.setZero(); 
        forces.setZero();
        masses.setOnes();
    } 
 
    Atoms(const Positions_t &p, const Velocities_t &v) : 
            positions{p}, velocities{v}, forces{3, p.cols()}, masses{p.cols()} {
        assert(p.cols() == v.cols());
        forces.setZero();
        masses.setOnes();
    }

    Eigen::Index nb_atoms() const {
        return positions.cols(); 
    }

    // Returns the length of the Atom structure in z direction.
    double length() {
        double smallest_z = DBL_MAX;
        double largest_z = 0;

        // Loop over all atoms.
        for (long i = 0; i < positions.cols() - 1; i++) {
            // Check is atom is positioned most left.
            if (positions(2, i) < smallest_z) {
                smallest_z = positions(2, i);
            }

            // Check is atom is positioned most right.
            if (positions(2, i) > largest_z) {
                largest_z = positions(2, i);
            }
        }

        // Return the length.
        return largest_z - smallest_z;
    }

    // Resizes size of the Atom structure to resize_factor.
    void resize(int resize_factor) {
        // Resize positions, velocities, forces and masses array of the Atom structure.
        positions.conservativeResize(3, resize_factor);
        velocities.conservativeResize(3, resize_factor);
        forces.conservativeResize(3, resize_factor);
        masses.conservativeResize(resize_factor);
    }
};

#endif
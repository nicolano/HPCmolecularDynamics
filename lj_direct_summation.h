#ifndef LJ_DIRECT_SUMMATION_H
#define LJ_DIRECT_SUMMATION_H

#include "atoms.h"

/*
 * Implementation of the Lennard-Jones potential with direct summation.
 * Returns the sum of the potential energy between all atoms in the Atoms structure.
 */
double lj_direct_summation(Atoms &atoms, double epsilon = 1.0, double sigma = 1.0);

#endif
#ifndef MOLDYN_KIN_ENERGY_H
#define MOLDYN_KIN_ENERGY_H

#include "atoms.h"

/*
 * Calculates and returns the kinetic energy of all atoms in an Atoms structure.
 */
double kin_energy(Atoms &atoms, double mass);

/*
 * Calculates and returns the kinetic energy for a subset of atoms in an Atoms
 * structure. The subset contains the atoms with indices i = 0 to i = nb_local.
 */
double kin_energy_sub(Atoms &atoms, double mass, long nb_local);

#endif // MOLDYN_KIN_ENERGY_H

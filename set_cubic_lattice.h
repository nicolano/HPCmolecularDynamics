#ifndef MOLDYN_SET_CUBIC_LATTICE_H
#define MOLDYN_SET_CUBIC_LATTICE_H

#include "atoms.h"

/*
 * Places all Atoms of an Atoms structure in a cubic lattice.
 */
void set_cubic_lattice(Atoms &atoms, double lattice_constant = 1.0);

#endif // MOLDYN_SET_CUBIC_LATTICE_H

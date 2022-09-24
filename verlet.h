#ifndef __VERLET_H
#define __VERLET_H

#include <Eigen/Dense>
#include "atoms.h"

/*
 * Perform the first step of the verlet integration.
 */
void verlet_step1(Atoms &atoms, double delta_t, double mass);

/*
 * Perform second step of the verlet integration.
 */
void verlet_step2(Atoms &atoms, double delta_t, double mass);

#endif  // __VERLET_H
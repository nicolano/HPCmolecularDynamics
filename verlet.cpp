#include "verlet.h"

void verlet_step1(Atoms &atoms, double delta_t, double mass) {
    // Calculate the new velocity from the forces acting on the atoms.
    atoms.velocities += 0.5 * atoms.forces * delta_t / mass;
    // Calculate the new velocity from the forces acting on the atoms.
    atoms.positions += atoms.velocities * delta_t;
}

void verlet_step2(Atoms &atoms, double delta_t, double mass) {
    // Calculate the new velocity from the forces acting on the atoms.
    atoms.velocities += 0.5 * atoms.forces * delta_t / mass;
}

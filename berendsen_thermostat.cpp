//
// Created by Nicolas von Trott on 25.05.22.
//
#include "berendsen_thermostat.h"

void berendsen_thermostat(Atoms &atoms, double e_kin, double temperature, double delta_t, double tau, double k_b) {
    // Calculate the current temperature in the cluster.
    double current_temperature {2 * e_kin / (3 * atoms.nb_atoms() * k_b)};
    // Calculate the rescale factor.
    double lambda {sqrt(1 + ((temperature / current_temperature) - 1) * delta_t / tau)};
    // Rescale velocities with the rescale factor.
    atoms.velocities *= lambda;
}

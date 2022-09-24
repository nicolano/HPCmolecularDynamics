//
// Created by Nicolas von Trott on 25.05.22.
//

#ifndef MOLDYN_BERENDSEN_THERMOSTAT_H
#define MOLDYN_BERENDSEN_THERMOSTAT_H

#include "atoms.h"

/*
 * Rescales the temperature in the cluster with a berendsen thermostat.
 */
void berendsen_thermostat(Atoms &atoms, double e_kin, double temperature, double delta_t, double tau, double k_b);

#endif // MOLDYN_BERENDSEN_THERMOSTAT_H

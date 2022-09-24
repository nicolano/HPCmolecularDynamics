//
// Created by Nicolas von Trott on 01.06.22.
//

#include <gtest/gtest.h>

#include "berendsen_thermostat.h"
#include "set_cubic_lattice.h"
#include "verlet.h"
#include "kin_energy.h"
#include "lj_direct_summation.h"
#include "atoms.h"


TEST(BerendsenThermostatTest, Energy) {
    //constexpr double BOLTZMANN_CONST = 1.380649e-23;
    // Declare variables
    const double BOLTZMANN_CONST = 1.0;
    const double epsilon = 1.0;
    const double sigma = 1.0;
    const double mass = 1.0;
    const double temperature = 0.5;
    const double relaxation_time = 0.1 ;
    const double t_max = 100 * sqrt(mass * sigma * sigma / epsilon);
    const double timestep = 0.001 * sqrt(mass * sigma * sigma / epsilon);
    const int nb_steps = t_max/timestep;
    const int nb_atoms = 27;
    const double lattice_constant = 1 * sigma;

    double e_kin_theo = 3.0 / 2.0 * nb_atoms * temperature * BOLTZMANN_CONST;
    double e_kin_avg = 0;
    double temperature_t = 0;

    // Store names, positions and velocities in atoms structure
    Atoms atoms(nb_atoms);
    // Setup atom grid
    set_cubic_lattice(atoms, lattice_constant);

    double start_temperature = 2.0 * kin_energy(atoms, mass) / ( 3.0 * atoms.nb_atoms() * BOLTZMANN_CONST );

    // Simulation
    for (int i = 0; i < nb_steps; ++i) {
        // Do verlet step 1
        verlet_step1( atoms, timestep, 1.0);

        double e_kin{kin_energy(atoms, mass)};
        // Calculate forces
        if (i >= 99000 && i < 100000) {
            e_kin_avg += e_kin;
        }
        lj_direct_summation(atoms, epsilon, sigma);

        // Do verlet step 2
        verlet_step2( atoms, timestep, 1.0);

        double current_temperature{2.0 * e_kin / ( 3.0 * atoms.nb_atoms() * BOLTZMANN_CONST )};
        double current_temperature_theo{temperature + ( start_temperature - temperature) * exp(- i / relaxation_time)};
        EXPECT_NEAR(current_temperature, current_temperature_theo, 1);

        // Do berendsen rescalling
        if (i > 0){
            berendsen_thermostat(atoms, e_kin, temperature, timestep, relaxation_time, BOLTZMANN_CONST);
        }
    }

    e_kin_avg = e_kin_avg / 1000;

    EXPECT_NEAR(e_kin_theo, e_kin_avg, 1e-5);
}
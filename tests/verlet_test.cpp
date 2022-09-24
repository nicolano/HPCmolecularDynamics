#include "verlet.h"
#include "types.h"
#include "bwgl.h"

#include <gtest/gtest.h>
#include <math.h>
#include <iostream>


TEST(VerletTest, positions)
{
    int nb_atoms = 1;
    int nb_steps = 10;
    int delta_t = 1;

    Atoms atoms(nb_atoms);

    Positions_t pos_0(3, nb_atoms);
    Positions_t pos_analytical(3, nb_atoms);
    Velocities_t vel_0(3, nb_atoms);

    atoms.positions.col(0) << 0.0, 5.0, 0.0;
    pos_0.col(0) << 0.0, 5.0, 0.0;
    pos_analytical.col(0) << 0.0, 5.0, 0.0;

    atoms.velocities.col(0) << 1.0, 0.0, 0.0;
    vel_0.col(0) << 1.0, 0.0, 0.0;

    atoms.forces.col(0) << 1.0, 0.5, 2.0;
    
    for (int i = 0; i < nb_steps; ++i) {
        // Do verlet step integration
        verlet_step1( atoms, delta_t, 1.0);
        verlet_step2( atoms, delta_t, 1.0);

        // get analytical solution for current time=(i+1)*delta_t
        bwgl(pos_analytical, vel_0, atoms.forces, (i+1)*delta_t);
        
        // loop through each position element 
        for (int n = 0; n < nb_atoms; ++n) {
            for (int m = 0; m < 3; ++m) {
                // Test analytical with numerical solution
                ASSERT_NEAR( atoms.positions(m,n), pos_analytical(m,n), 0.01);
            }
        }      

        // reset pos_analytical to start position
        pos_analytical = pos_0;
    }
}



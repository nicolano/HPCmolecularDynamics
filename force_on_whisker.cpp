//
// Created by Nicolas von Trott on 01.09.22.
//
#include "force_on_whisker.h"

double force_on_whisker(Atoms &atoms, long nb_local, double z_limit, bool isLeft) {
    double force = 0;

    for (long i = 0; i < atoms.nb_atoms() - 1; i++) {
        double z = atoms.positions.coeff(2,i);
        Eigen::Vector3d f(atoms.forces.col(i));
        if (((isLeft) && (z < z_limit))) {
            // sum up.
            force += f.norm();
        } else if ((!isLeft) && (z > z_limit)) {
            // sum up.
            force += f.norm();
        }
    }
    return force;
}


//Eigen::Vector3d force_on_whisker(Atoms &atoms, long nb_local, double z_length, bool isLeft) {
//    Eigen::Vector3d force(0,0,0);
//
//    std::cout << "----------------------------------" << "\n";
//
//    for (long i = 0; i < atoms.nb_atoms() - 1; i++) {
//        double z = atoms.positions.coeff(2,i);
//        if (((isLeft) && (z < 0)) || ((!isLeft) && (z > z_length))) {
//            Eigen::Vector3d f(atoms.forces.col(i));
//            // sum up.
//            force += f;
//        }
//    }
//
//    return force;
//}
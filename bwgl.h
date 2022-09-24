#ifndef __BWGL_H
#define __BWGL_H

#include <Eigen/Dense>

void bwgl(Eigen::Array3Xd &positions, Eigen::Array3Xd &velocities,
                    const Eigen::Array3Xd &forces, double t);

#endif  // __VERLET_H
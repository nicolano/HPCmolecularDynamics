#include "bwgl.h"

void bwgl(Eigen::Array3Xd &positions, Eigen::Array3Xd &velocities,
                    const Eigen::Array3Xd &forces, double t) {

    positions += 0.5 * forces * t * t + velocities * t;
}
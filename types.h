#ifndef TYPES_H
#define TYPES_H

#include <Eigen/Dense>

using Names_t = Eigen::RowVectorX<std::string>;
using Positions_t = Eigen::Array3Xd;
using Velocities_t = Eigen::Array3Xd;
using Forces_t = Eigen::Array3Xd;
using Masses_t = Eigen::RowVectorX<std::double_t>;

#endif
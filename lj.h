//
// Created by Nicolas von Trott on 01.06.22.
//

#ifndef MOLDYN_LJ_H
#define MOLDYN_LJ_H

#include "neighbors.h"

double lj(Atoms &atoms, const NeighborList& neighbor_list, double cutoff_range, double epsilon = 1.0, double sigma = 1.0);

#endif // MOLDYN_LJ_H

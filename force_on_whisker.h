//
// Created by Nicolas von Trott on 01.09.22.
//

#ifndef MOLDYN_FORCE_ON_WHISKER_H
#define MOLDYN_FORCE_ON_WHISKER_H

#include "Atoms.h"

/*
 * Calculates the force on a whisker.
 */
double force_on_whisker(Atoms &atoms, long nb_local, double z_limit, bool isLeft);

#endif // MOLDYN_FORCE_ON_WHISKER_H

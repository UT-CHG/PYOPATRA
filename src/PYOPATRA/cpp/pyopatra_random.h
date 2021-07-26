//
// Created by georgia on 7/26/21.
//

#ifndef PYOPATRA_PYOPATRA_RANDOM_H
#define PYOPATRA_PYOPATRA_RANDOM_H

#include <random>

static std::normal_distribution<double> normal(0, 1);
static std::uniform_real_distribution<double> unif_pi(0, M_PI);
static std::uniform_real_distribution<double> unif(0, 1);
static std::default_random_engine generator;

#endif //PYOPATRA_PYOPATRA_RANDOM_H

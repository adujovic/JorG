#ifndef _SOLVER_H
#define _SOLVER_H

#include <iostream>
#include <iomanip>
#include <fstream>

#include <vector>
#include <array>
#include <unordered_set>
#include <numeric>
#include <random>
#include <regex>

#include <cmath>
#include <cstring>
#include <string>

#include "../asa.h"
#include "../ising.h"
#include "../arithmeticvector.h"

#include "aux.h"

extern "C"
int solver(char _basis[],char _supercell[],char _flippable[], size_t reference, size_t unique_flips, size_t ansatz=0U);

#endif

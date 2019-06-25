#ifndef _D_ISING_H
#define _D_ISING_H

#include "ising.h"
#include <memory>

namespace ising{

std::unique_ptr< ising::AbstractIsingModel > factoryIsing(const size_t& N);

} //end of namespace ising
#endif

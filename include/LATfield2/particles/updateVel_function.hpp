#ifndef LATFIELD2_UPDATEVEL_FUNCTION_HPP
#define LATFIELD2_UPDATEVEL_FUNCTION_HPP

#include "LATfield2/macros.hpp"
#include "LATfield2/Field.hpp"
#include "LATfield2/Site.hpp"
#include "LATfield2/Imag.hpp"
#include "LATfield2/particles/simple.hpp"

namespace LATfield2
{

/**
 * \addtogroup prartClass
 * @{
 */

Real updateVel_simple(double /*dtau*/,
                          double /*lat_resolution*/,
                          part_simple * part,
                          double * /*ref_dist*/,
                          part_simple_info /*partInfo*/,
                          Field<Real> ** /*fields*/,
                          Site * /*sites*/,
                          int /*nfield*/,
                          double * /*params*/,
                          double * /*outputs*/,
                          int /*noutputs*/);

/**@}*/
}
#endif

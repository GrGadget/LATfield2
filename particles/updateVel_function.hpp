#ifndef LATFIELD2_UPDATEVEL_FUNCTION_HPP
#define LATFIELD2_UPDATEVEL_FUNCTION_HPP


#include "LATfield2_macros.hpp"
#include "LATfield2_Field.hpp"
#include "LATfield2_Site.hpp"
#include "Imag.hpp"
#include "particles/LATfield2_particle_simple.hpp"

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

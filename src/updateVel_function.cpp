#include "LATfield2/particles/updateVel_function.hpp"
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
                          int /*noutputs*/)
{

    double v2{};

    
    for(int i=0;i<3;i++)
    {
        (*part).vel[i] = (*part).vel[i];
        v2 += (*part).vel[i] * (*part).vel[i];
    }
    
    
    return v2;

}

/**@}*/
}

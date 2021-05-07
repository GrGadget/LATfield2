#include "config.h"
#include "particles/LATfield2_particle_rk4.hpp"

#include <iostream>

namespace LATfield2
{
std::ostream& operator<<(std::ostream& os, const part_rk4& p)
{
    os << "ID: "<<p.ID<<", mass: "<< p.mass<<" , Pos: ("<< p.pos[0]<<","<< p.pos[1]<<","<< p.pos[2]<<") , Vel: (" << p.vel[0]<<","<< p.vel[1]<<","<< p.vel[2]<<")";
    return os;
}
}

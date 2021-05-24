#include "LATfield2/particles/tools.hpp"
#include "LATfield2/parallel2d.hpp"

#include <iostream>

namespace LATfield2{
using std::endl;
using std::cout;

Real get_lattice_resolution(const int npts[3],const LATfield2::Real boxSize[3])
{
  LATfield2::Real latRes[3];

  for(int i =0;i<3;i++)latRes[i]=boxSize[i]/ (LATfield2::Real)npts[i];

  if(latRes[0]==latRes[1] && latRes[0]==latRes[2])return latRes[0];
  else{
    COUT<< "wrong physical box size and lattice size, relosution must be same in each dimensions"<<endl;
    COUT<< " Exiting... "<<endl;
    std::exit(222);
      return -1;
  }
}
}

#ifndef LATFIELD2_HPP
#define LATFIELD2_HPP
#include "config.h"

/*! \file LATfield2.hpp
 \brief LATfield2 header
 \author David Daverio,Neil Bevis

 */

#include "LATfield2_macros.hpp"

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <typeinfo>
#include <list>
#include <vector>


#include "LATfield2_Lattice.hpp"

#ifdef EXTERNAL_IO
#include "LATfield2_IO_server.hpp"
namespace LATfield2
{
    IOserver ioserver;
}
#endif

#include "LATfield2_parallel2d.hpp"
namespace LATfield2
{
    extern Parallel2d parallel;
}

#include "LATfield2_SettingsFile.hpp"



#ifdef H5_HAVE_PIXIE
#include "LATfield2_save_hdf5_pixie.h"
#else
#include "LATfield2_save_hdf5.hpp"
#endif


#include "Imag.hpp"
#include "int2string.hpp"
#include "LATfield2_Lattice.hpp"
#include "LATfield2_Site.hpp"
#include "LATfield2_Field.hpp"
#ifdef FFT3D
    #include "LATfield2_PlanFFT.hpp"
#endif
#include "particles/LATfield2_Particles.hpp"
#ifdef CATALAT
    #include "LATfield2_catalyst.hpp"
#endif

#endif

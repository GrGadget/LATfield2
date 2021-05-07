#ifndef LATFIELD2_HPP
#define LATFIELD2_HPP
#include "LATfield2/config.h"

/*! \file LATfield2.hpp
 \brief LATfield2 header
 \author David Daverio,Neil Bevis

 */


#ifdef EXTERNAL_IO
#include "LATfield2/IO_server.hpp"
namespace LATfield2
{
    IOserver ioserver;
}
#endif

#include "LATfield2/macros.hpp"
#include "LATfield2/Lattice.hpp"
#include "LATfield2/SettingsFile.hpp"
#include "LATfield2/parallel2d.hpp"
namespace LATfield2
{
    extern Parallel2d parallel;
}

#ifdef H5_HAVE_PIXIE
#include "LATfield2/save_hdf5_pixie.h"
#else
#include "LATfield2/save_hdf5.hpp"
#endif


#include "LATfield2/Imag.hpp"
#include "LATfield2/int2string.hpp"
#include "LATfield2/Lattice.hpp"
#include "LATfield2/Site.hpp"
#include "LATfield2/Field.hpp"
#ifdef FFT3D
    #include "LATfield2/PlanFFT.hpp"
#endif
#include "LATfield2/particles/Particles.hpp"
#ifdef CATALAT
    #include "LATfield2/catalyst.hpp"
#endif

#endif

#ifndef PARTICLES_TOOLS_HPP
#define PARTICLES_TOOLS_HPP
#include "config.h"

#include "LATfield2_macros.hpp"
#include "Imag.hpp"
#include <iostream>
namespace LATfield2
{

#ifndef DOXYGEN_SHOULD_SKIP_THIS
Real get_lattice_resolution(int npts[3],Real boxSize[3]);

#define CREATE_MEMBER_DETECTOR(X)            \
  template<typename T> struct has_##X {      \
    struct Fallback {int X; };               \
    struct Derived : T, Fallback { };        \
    template<typename C, C> struct ChT;      \
    template<typename C> static char (&f(ChT<int Fallback::*, &C::X>*))[1]; \
    template<typename C> static char (&f(...))[2]; \
    static bool const value = sizeof(f<Derived>(0)) == 2; \
  };


#define CREATE_MEMBER_DETECTOR_MAXI(X)          \
    template<typename T> struct has_maxi_##X {      \
    struct Fallback {int X; };               \
    struct Derived : T, Fallback { };        \
    template<typename C, C> struct ChT;      \
    template<typename C,typename CC> int f(ChT<int Fallback::*, &C::X>*){return -1;} \
    template<typename C,typename CC> int f(...){return offsetof(T,X);} \
    int gos(){return f<Derived,T>(0);} \
};

#endif

}

#endif

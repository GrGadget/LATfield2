#ifndef LATFIELD_MACROS
#define LATFIELD_MACROS 

#include "config.h"

#ifdef SINGLE
#define REAL_TYPE H5T_NATIVE_FLOAT
#else
#define REAL_TYPE H5T_NATIVE_DOUBLE
#endif

#ifndef RealC
#ifdef SINGLE
#define RealC float
#define MPI_RealC MPI_FLOAT
#else
#define RealC double
#define MPI_RealC MPI_DOUBLE
#endif
#endif

#ifdef BIG_ENDIAN_ORDER
  #define DATA_ORDER H5T_ORDER_BE
  #define LONG_TYPE_FILE H5T_STD_I64BE
  #define INT_TYPE_FILE H5T_STD_I32BE
  #define FLOAT_TYPE_FILE H5T_IEEE_F32BE
  #define DOUBLE_TYPE_FILE H5T_IEEE_F64BE
  #define BOOL_TYPE_FILE H5T_STD_I8BE
  #ifdef SINGLE
    #define REAL_TYPE_FILE H5T_IEEE_F32BE
  #else
    #define REAL_TYPE_FILE H5T_IEEE_F64BE
  #endif
#else
  #define  DATA_ORDER H5T_ORDER_LE
  #define LONG_TYPE_FILE H5T_STD_I64LE
  #define INT_TYPE_FILE H5T_STD_I32LE
  #define FLOAT_TYPE_FILE H5T_IEEE_F32LE
  #define DOUBLE_TYPE_FILE H5T_IEEE_F64LE
  #define BOOL_TYPE_FILE H5T_STD_I8LE
  #ifdef SINGLE
    #define REAL_TYPE_FILE H5T_IEEE_F32LE
  #else
    #define REAL_TYPE_FILE H5T_IEEE_F64LE
  #endif
#endif



/*! \file looping_macro.hpp
 \brief looping over the lattice site macro
 \author David Daverio

 */

#define forallboundary_start(latttt,xxxx) \
Site xxxx(latttt); \
for(xxxx.haloFirst();xxxx.haloTest();xxxx.haloNext()) \
{ \
bool inBoundary = false; \
for(int i=0; i<latttt.dim(); i++)if(xxxx.coord(i)>=latttt.size(i) || xxxx.coord(i)<0)inBoundary=true; \
if(inBoundary==true) \
{



#define forallboundary_stop }};

#endif

#ifndef LATFIELD_MACROS
#define LATFIELD_MACROS 


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

#define FFT3D


#endif

#ifndef LATFIELD2_LATTICE_DECL_HPP
#define LATFIELD2_LATTICE_DECL_HPP

/*! \class Lattice
 
 \brief The Lattice class describe a cartesian mesh (with 2 or more dimensions). The updateHalo method of the Field class generate the periodicity.
 
 
 It store the global and local geometry of the mesh. The last 2 dimension of the lattice are scattered into the MPI processes grid.
 
 
 */

#include "LATfield2/macros.hpp"
#include <string>

namespace LATfield2
{
 
 using std::string;
 
class Lattice
{
public:
    //! Constructor.
    Lattice();
    
    Lattice(const Lattice&) = delete;
    Lattice& operator=(const Lattice&) = delete;
    
    /*!
     Constructor with initialization
     \sa initialize(int dim, const int* size, int halo);
     \param dim : number of dimension
     \param size : array containing the size of each dimension.
     \param halo : size of the halo (ghost cells, same for each dimension)
     */
    Lattice(int dim, const int* size, int halo);
    
    /*!
     Constructor with initialization
     \sa initialize(int dim, const int size, int halo);
     \param dim : number of dimension
     \param size : size of each dimension (same for each dimension)
     \param halo : size of the halo (same for each dimension)
     */
    Lattice(int dim, const int size, int halo);
    
    //! Destructor.
    ~Lattice();
    
    /*!
     Initialization of a dim-dimensional lattice, the size of each dimension is set by the second parameter: int *size. The ghost cell number (halo) is the same for each dimension.
     
     \param dim : number of dimension
     \param size : array containing the size of each dimension.
     \param halo : size of the halo (same for each dimension)
     */
    void initialize(int dim, const int* size, int halo);
    
    /*!
     Initialization of a dim-dimensional lattice, each dimension have the same size. The ghost cell number (halo) is the same for each dimension.
     
     
     \param dim : number of dimension
     \param size : size of each dimension (same for each dimension)
     \param halo : size of the halo (same for each dimension)
     */
    void initialize(int dim, const int size, int halo);
    
    
#ifdef FFT3D
    /*!
     Initialization of a lattice for Fourier space in case of real to complex transform. The Fourier space lattice size is defined according to the real space one. The fourier space lattice have "halo" ghost cells in each dimension (which can be different than the halo of the real space lattice).
     \param lat_real : pointer to a real space lattice.
     \param halo : size of the halo (same for each dimension)
     */
    void initializeRealFFT(Lattice & lat_real, int halo);
    
    /*!
     Initialization of a lattice for Fourier space in case of complex to complex transform. The Fourier space lattice size is defined according to the real space one.. The fourier space lattice have "halo" ghost cells in each dimension (which can be different than the halo of the real space lattice).
     \param lat_real : pointer to a real space lattice.
     \param halo : size of the halo (same for each dimension)
     */
    void initializeComplexFFT(Lattice & lat_real, int halo);
#endif
    
    
    
    /*!
     \return int. Number of dimensions of the lattice.
     */
    int  dim()const {return dim_;}
    /*!
     \return int. Size of the halo (ghost cells).
     */
    int  halo() const { return halo_; }
    
    /*!
     \return int*. Pointer to the array of the size of each dimension of the lattice.
     */
    const int* size() const{ return size_; }
    
    /*!
     Function which returns the size of a given dimension of the lattice.
     \param direction : asked dimension.
     \return int. Global size of the lattice in the given dimension.
     */
    int  size(int i) const { return size_[i]; }
    
    /*!
     \return int*. Pointer to the array of the size of each dimension of the sublattice stored in this MPI process.
     */
    const int * sizeLocal()const { return sizeLocal_; }
    
    /*!
     Function which returns the size of a given dimension of the sublattice stored in this MPI process.
     \param direction : asked dimension.
     \return int. Global size of the sublattice (of this MPI process) in the given dimension.
     */
    int  sizeLocal(int i)const { return sizeLocal_[i]; }
    
    
    /*!
     \return long. Number of sites on the lattice (excluding halo sites).
     */
    long  sites() const { return sites_; }
    /*!
     \return long. Number of sites on the lattice (including halo sites).
     */
    long  sitesGross() const { return sitesGross_; }
    
    /*!
     \return long. Number of sites (excluding halo sites) of the sublattice stored in this MPI process.
     */
    long  sitesLocal() const { return sitesLocal_; }
    
    /*!
     \return long. Number of sites (including halo sites) of the sublattice stored in this MPI process.
     */
    long  sitesLocalGross()const { return sitesLocalGross_; }
    
    /*!
     \return long. Array index of the first site which is not within the halo.
     */
    long siteFirst() const { return siteFirst_; }
    
    /*!
     \return long. Array index of the last site which is not within the halo.
     */
    long siteLast() const { return siteLast_; }
    
    
    
    /*!
     Function which return the number of data_ array elements to jump to move to the next site in the given direction. (does not take into account the number of component of the fields, therefor should be multiplied by Field.components().) Should not be used by user.
     \param direction : asked direction.
     \return long. Number of array elements to jump.
     */
    long  jump(int i) const { return jump_[i]; }
    
    /*!
     \return long. Number of sites before first local site in lattice. Should not be used by users.
     */
    long  sitesSkip()const { return sitesSkip_; }
    
    /*!
     \return long. Number of sites before first local site in lattice. Should not be used by users.
     */
    long  sitesSkip2d() const { return sitesSkip2d_; }
    
    /*!
     \return *long. Pointer to an array which store the last 2 dimensions coordinate of the first local(in this MPI process) sites. Index 0 is for dim-1, index 1 is for dim-2/
     */
    const long*  coordSkip()const { return coordSkip_; }
    
    /*!
     Function which save in serial and in ASCII the global and local description of the Lattice. Usefull to read a file writen by fast_save or fast_write methods of the Field class.
     \param filename : filename of the architectur file.
     */
    void save_arch(const string filename);
    
    /*!
     \return return true if the description of the lattice has been written on disk.
     
     \sa save_arch(const string filename)
     */
    bool is_arch_saved() const {return arch_saved_;}
    
    int getRank(int* coord) ; //return the world rank of the process who get the lattices site "coord"
    int getRankDim0(int coord) ;
    int getRankDim1(int coord) ;
    
    const int * sizeLocalAllProcDim0()const{ return sizeLocalAllProcDim0_; }
    const int * sizeLocalAllProcDim1()const{ return sizeLocalAllProcDim1_; }
    
    
private:
    int        status_;
    static int initialized;
    
    //Global variables==============
    int  dim_;              //ok//Number of dimensions
    int* size_;             //ok//Number of lattice sites in each direction
    long  sites_;            //ok//Number of sites in lattice
    long  sitesGross_;       //ok//Number of sites in lattice plus halos
    int  halo_;             //ok//Number of sites extra in each direction
    
    //Local variables===============
    int* sizeLocal_;       //ok//Number of local lattice sites in each direction
    int* sizeLocalAllProcDim0_;
    int* sizeLocalAllProcDim1_;
    long  sitesLocal_;      //ok//Number of local sites in lattice
    long  sitesLocalGross_; //ok//Number of local sites in lattice plus halo
    long* jump_;            //ok//Jumps needed to move up one in each direction
    
    
    long  siteFirst_;       //ok//Index of first local site in lattice
    
    long  siteLast_;        //ok//Index of last local site in lattice
    long  sitesSkip_;      //Number of global lattice sites before first local site (say in a file)
    long  sitesSkip2d_;
    long  coordSkip_[2];       //Number to add to coord[dim_-1] and coord[dim_-2] to get global value
    
    
    //save variable for fast save
    int arch_saved_;
    
};

}
#endif

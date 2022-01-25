#pragma once

/*! \class Lattice
 
 \brief The Lattice class describe a cartesian mesh (with 2 or more dimensions). The updateHalo method of the Field class generate the periodicity.
 
 
 It store the global and local geometry of the mesh. The last 2 dimension of the lattice are scattered into the MPI processes grid.
 */

#include "LATfield2/exception.hpp"
#include <boost/mpi/cartesian_communicator.hpp>

namespace LATfield2
{

// TODO: API for fourier space lattice
// TODO: restrict to 3dimensions
 
class Lattice 
{
    using communicator = ::boost::mpi::communicator;
    using cartesian_topology = ::boost::mpi::cartesian_topology;
    using cartesian_communicator = ::boost::mpi::cartesian_communicator;
    
    struct domain_decomposition
    {
        using int_type = long long;
        
        static int_type local_size(
            int_type total_size, int_type N_domains, int_type my_domain)
        {
            const int_type d = total_size / N_domains;
            const int_type n1 = total_size - d*N_domains; // have (d+1) size
            const int_type n0 = N_domains - n1; // have (d) size
            
            return my_domain < n0 ? d : d+1;
        }
        static int_type local_cumulative_size(
            int_type total_size, int_type N_domains, int_type my_domain)
        {
            const int_type d = total_size / N_domains;
            const int_type n1 = total_size - d*N_domains; // have (d+1) size
            const int_type n0 = N_domains - n1; // have (d) size
            
            return my_domain < n0 ? d*my_domain : n0*d + (my_domain-n0)*(d+1);
        }
        static int_type domain_coordinate(
            int_type total_size, int_type N_domains, int_type mesh_position)
        {
            const int_type d = total_size / N_domains;
            const int_type n1 = total_size - d*N_domains; // have (d+1) size
            const int_type n0 = N_domains - n1; // have (d) size
            
            const int_type len0 = n0*d;
            
            return mesh_position < len0 ? 
                mesh_position/d :
                n0 + (mesh_position-len0)/(d+1) ;
        }
    };
    
    struct check_comm_size
    {
        check_comm_size(
            const cartesian_topology& mesh_top, 
            const cartesian_topology& comm_top, 
            const communicator& comm
            )
        {
            if(mesh_top.stl().size() != 3)
                throw bad_dimensions(
                    "Lattice::Lattice: mesh topology must have 3 dimensions");
            if(comm_top.stl().size() != 2)
                throw bad_dimensions(
                    "Lattice::Lattice: communicator topology "
                    "must have 2 dimensions");
                    
            long long tot_size = 1;
            for(auto dim : comm_top.stl())
            {
                tot_size *= dim.size;
            }
            if(tot_size!=comm.size())
                throw bad_dimensions(
                    "Lattice::check_comm_size: "
                    "the processor grid size is different from the communicator size");
            
            for(auto dim : comm_top.stl())
            if(dim.periodic==false)
            {
                throw bad_dimensions(
                    "Lattice::check_comm_size: "
                    "the communicator topology must be periodic");
            }
            for(auto dim : mesh_top.stl())
            if(dim.periodic==false)
            {
                throw bad_dimensions(
                    "Lattice::check_comm_size: "
                    "the mesh topology must be periodic");
            }
        }
        
    } mesh_check_comm_size;

    
    
    // entire mesh data
    cartesian_topology mesh_topology;
    cartesian_communicator mesh_communicator;
    const int mesh_halo;
    std::vector<int> mesh_size; //Number of Global lattice sites in each direction
    std::vector<int> mesh_commSize;
    
    // this process cached data
    std::vector<int> my_size; // Number of local lattice sites in each direction
    std::vector<int> my_offset; // Starting position in the global lattice
    
public:
    //! Constructor.
    Lattice(
        const cartesian_topology& mesh_top, 
        const cartesian_topology& comm_top, 
        const communicator& comm,
        int halo):
            mesh_check_comm_size(mesh_top,comm_top,comm),
            mesh_topology{mesh_top},
            mesh_communicator(comm,comm_top),
            mesh_halo{halo},
            
            mesh_size(mesh_communicator.ndims()),
            mesh_commSize(mesh_communicator.ndims()),
            
            my_size(mesh_communicator.ndims()),
            my_offset(mesh_communicator.ndims())
    {
        
        const cartesian_topology domain_topology = mesh_communicator.topology();
         
        const auto domain_topvector = domain_topology.stl();
        const auto mesh_topvector = mesh_topology.stl();
        const auto coordinates = mesh_communicator.coordinates(mesh_communicator.rank());
        const int ndim = domain_topvector.size();
        for(int i = 0 ;i<ndim; ++i) 
        {
            mesh_commSize[i] = domain_topvector[i].size;
            
            mesh_size[i] = mesh_topvector[i].size;
            
            my_size[i] =
                domain_decomposition::local_size(
                    mesh_topvector[i].size,
                    domain_topvector[i].size,
                    coordinates[i]);
                    
            my_offset[i]=
                domain_decomposition::local_cumulative_size(
                    mesh_topvector[i].size,
                    domain_topvector[i].size,
                    coordinates[i]);
        }
    }
    
    // int mesh_offset(int dom_coordinate, int dim) const
    // {
    //     return domain_decomposition::local_cumulative_size(
    //         my_sizeGlobal[dim],
    //         my_commSize[dim],
    //         dom_coordinate
    //         );
    // }
    // int domain_coordinate(int mesh_position, int dim) const
    // {
    //     return domain_decomposition::domain_coordinate(
    //         my_sizeGlobal[dim],
    //         my_commSize[dim],
    //         mesh_position
    //     );
    // }
    
    
    ~Lattice(){}
    
    /*!
     \return int. Number of dimensions of the lattice.
     */
    int  dim()const {return mesh_communicator.ndims();}
    /*!
     \return int. Size of the halo (ghost cells).
     */
    int  halo() const { return mesh_halo; }
    
    /*!
     \return int*. Pointer to the array of the size of each dimension of the lattice.
     */
    const auto& size() const{ return mesh_size; }
    
    /*!
     Function which returns the size of a given dimension of the lattice.
     \param direction : asked dimension.
     \return int. Global size of the lattice in the given dimension.
     */
    int  size(int i) const { return mesh_size[i]; }
    
    /*!
     \return int*. Pointer to the array of the size of each dimension of the sublattice stored in this MPI process.
     */
    const auto& sizeLocal()const { return my_size; }
    
    /*!
     Function which returns the size of a given dimension of the sublattice stored in this MPI process.
     \param direction : asked dimension.
     \return int. Global size of the sublattice (of this MPI process) in the given dimension.
     */
    int sizeLocal(int i)const { return my_size[i]; }
    
    
    // /*!
    //  \return long. Number of sites on the lattice (excluding halo sites).
    //  */
    // long  sites() const { return sites_; }
    // /*!
    //  \return long. Number of sites on the lattice (including halo sites).
    //  */
    // long  sitesGross() const { return sitesGross_; }
    // 
    // /*!
    //  \return long. Number of sites (excluding halo sites) of the sublattice stored in this MPI process.
    //  */
    // long  sitesLocal() const { return sitesLocal_; }
    // 
    // /*!
    //  \return long. Number of sites (including halo sites) of the sublattice stored in this MPI process.
    //  */
    // long  sitesLocalGross()const { return sitesLocalGross_; }
    // 
    // /*!
    //  \return long. Array index of the first site which is not within the halo.
    //  */
    // long siteFirst() const { return siteFirst_; }
    // 
    // /*!
    //  \return long. Array index of the last site which is not within the halo.
    //  */
    // long siteLast() const { return siteLast_; }
    // 
    // 
    // 
    // /*!
    //  Function which return the number of data_ array elements to jump to move to the next site in the given direction. (does not take into account the number of component of the fields, therefor should be multiplied by Field.components().) Should not be used by user.
    //  \param direction : asked direction.
    //  \return long. Number of array elements to jump.
    //  */
    // long  jump(int i) const { return jump_[i]; }
    // 
    // /*!
    //  \return long. Number of sites before first local site in lattice. Should not be used by users.
    //  */
    // long  sitesSkip()const { return sitesSkip_; }
    // 
    // /*!
    //  \return long. Number of sites before first local site in lattice. Should not be used by users.
    //  */
    // long  sitesSkip2d() const { return sitesSkip2d_; }
    // 
    // /*!
    //  \return *long. Pointer to an array which store the last 2 dimensions coordinate of the first local(in this MPI process) sites. Index 0 is for dim-1, index 1 is for dim-2/
    //  */
    // const long*  coordSkip()const { return coordSkip_; }
    
private:
    
    // //Global variables==============
    // long  sites_;            //ok//Number of sites in lattice
    // long  sitesGross_;       //ok//Number of sites in lattice plus halos
    // 
    // //Local variables===============
    // long  sitesLocalGross_; //ok//Number of local sites in lattice plus halo
    // std::vector<long> jump_;            //ok//Jumps needed to move up one in each direction
    // 
    // 
    // long  siteFirst_;       //ok//Index of first local site in lattice
    // 
    // long  siteLast_;        //ok//Index of last local site in lattice
    // long  sitesSkip_;      //Number of global lattice sites before first local site (say in a file)
    // long  sitesSkip2d_;
    // long  coordSkip_[2];       //Number to add to coord[dim_-1] and coord[dim_-2] to get global value
};

}

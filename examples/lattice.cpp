#include <LATfield2.hpp>
#include <iostream>

namespace mpi = ::boost::mpi;
int main()
{
    mpi::environment env; // required, it sets the MPI context, like MPI_Init
    mpi::communicator comm_world;
    
    std::vector < mpi::cartesian_dimension > comm_size{{2,true},{3,true}};
    std::vector < mpi::cartesian_dimension > mesh_size{{100,true},{100,true},{100,true}};
    
    // notice that std::vector< cartesian_dimension > is convertible to
    // cartesian_topology 
    try{
        LATfield2::Lattice lat(mesh_size,comm_size,comm_world,0); 
        
        for(int p=0;p<comm_world.size();++p)
        {
            if(comm_world.rank()==p)
            {
                std::cout << "rank " << p << "\n";
                std::cout << "  local size = ";
                for(size_t i=0;i<mesh_size.size();++i)
                    std::cout << " " << lat.sizeLocal(i) ;
                std::cout << "\n";
                std::cout << "  coordinates in the processor grid:";
                for(size_t i=0;i<comm_size.size();++i)
                    std::cout << " " << lat.coordinates(i) ;
                std::cout << '\n';
            }
        }
        
    }catch(std::exception & e)
    {
        std::cerr <<  "std::exception: " << e.what() << "\n";
    }catch(...)
    {
        std::cerr << "Unexpected exception.\n";
    }
    return 0;
}

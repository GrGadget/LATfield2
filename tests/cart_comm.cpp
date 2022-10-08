#define BOOST_TEST_MODULE Test CART_COMM
#define FFT3D
#include <boost/mpi.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/test/unit_test.hpp>

#include "LATfield2.hpp"

using namespace LATfield2;
namespace ut = boost::unit_test;
namespace mpi = boost::mpi;

/*
    From Stroustrup's The C++ Programming Language 4ed
    section 13.3.1
*/
template <class F>
struct Final_action
{
    F clean;
    Final_action(F f) : clean{f} {}
    ~Final_action() { clean(); }
};
template <class F>
Final_action<F> finally(F f)
{
    return Final_action<F>(f);
}

struct fixture
{
    fixture() {}

    ~fixture() {}

    static mpi::environment env;
    static mpi::communicator world;
};

mpi::environment fixture::env;
mpi::communicator fixture::world;

BOOST_TEST_GLOBAL_FIXTURE(fixture);

BOOST_AUTO_TEST_CASE(parallel_initialize)
{   
    // parameters
    const int n=4,m=3;
    
    // latfield setup
    parallel.initialize(fixture::world,n,m);
    
    const auto cart_comm = parallel.cartesian_communicator();
    
    BOOST_CHECK_EQUAL(parallel.rank(),cart_comm.rank());
    BOOST_CHECK_EQUAL(parallel.size(),cart_comm.size());
    
    BOOST_CHECK_EQUAL(parallel.world_rank(),cart_comm.rank());
    BOOST_CHECK_EQUAL(parallel.world_size(),cart_comm.size());
    
    const auto coord = parallel.coordinates();
    const auto top = parallel.grid_topology();
   
    for(int d=0;d<2;++d)
    {
        BOOST_CHECK_EQUAL(parallel.grid_size()[d],top[d]);
        BOOST_CHECK_EQUAL(parallel.grid_rank()[d],coord[d]);
    }
    
    BOOST_CHECK_EQUAL(parallel.grid2world(coord[0],coord[1]),cart_comm.rank());
    BOOST_CHECK_EQUAL(parallel.rank(coord),cart_comm.rank());
    
    // TODO: uncomment this line when all the bullshit has been removed
    // BOOST_CHECK_EQUAL(cart_comm.rank(coord),cart_comm.rank());
}

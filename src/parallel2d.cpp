#include "LATfield2/parallel2d.hpp"
#include <iostream>

namespace LATfield2
{
    using std::cerr;
    using std::endl;
    using std::cout;

/////////// add verif for the 2*2 process and dim>2.



/*! \file LATfield2_parallel2d.hpp
 \brief LATfield2_parallel2d.hpp contains the class Parallel2d implementation.
 \author David Daverio, edited by Wessel Valkenburg
 */


Parallel2d::Parallel2d() : neverFinalizeMPI(false)
{
}

void Parallel2d::initialize(MPI_Comm com,int proc_size0, int proc_size1)
{
    this->initialize(com,proc_size0, proc_size1,0,0);
}


void Parallel2d::initialize(MPI_Comm com, int proc_size0, int proc_size1,int, int)
// initializes the state of the Parallel2d variable once we know the MPI communicator
// and the processor grid.

// Example:
// proc_size0 = 4
// proc_size1 = 3
//
//   Process grid will look like this:
//
//     8 9 10 11
//     4 5 6  7
//     0 1 2  3
//
// and they will be split into several MPI groups
// 
// dim0_group_[0] = {0,1,2,3}
// dim0_group_[1] = {4,5,6,7}
// dim0_group_[2] = {8,9,10,11}
//
// dim1_group_[0] = {0,4,8}
// dim1_group_[1] = {1,5,9}
// dim1_group_[2] = {2,6,10}
// dim1_group_[3] = {3,7,11}

{
    world_comm_ = com;
    MPI_Comm_rank( world_comm_, &world_rank_ );
    MPI_Comm_size( world_comm_, &world_size_ );
    
    /*
        For Boost::MPI the fastest running index is the last,
        while for LATfield is the first.
        To introduce a LATfield2-compatible process topology
        I had to invert the order of proc_sizeX.
    */
    my_topology = ::boost::mpi::cartesian_topology{
        {proc_size1,/* periodicity = */ true},
        {proc_size0,/* periodicity = */ true}};
    
    my_communicator = ::boost::mpi::communicator(
                        world_comm_,
                        ::boost::mpi::comm_duplicate);
    
    // the root process has special powers,
    // eg. root writes the snapshot headers and error messages
    root_=0;
    isRoot_ = (root_ == my_communicator.rank());
                        
    my_cartesian_communicator.reset(new ::boost::mpi::cartesian_communicator(
                    my_communicator,
                    my_topology
                    ) );
    
    const auto coord = my_cartesian_communicator->coordinates(my_communicator.rank());
    
    // communicator for processes that share the same coord[0]
    my_axis_communicator[0] = my_communicator.split(coord[0],my_communicator.rank());
    
    // communicator for processes that share the same coord[1]
    my_axis_communicator[1] = my_communicator.split(coord[1],my_communicator.rank());
                    
                    
    #ifndef EXTERNAL_IO
    lat_world_comm_ = com;
    MPI_Comm_rank( lat_world_comm_, &lat_world_rank_ );
    MPI_Comm_size( lat_world_comm_, &lat_world_size_ );
    #endif

    grid_size_[0]=proc_size0;
	grid_size_[1]=proc_size1;
    
	dim0_comm_ = (MPI_Comm *)malloc(grid_size_[1]*sizeof(MPI_Comm));
	dim1_comm_ = (MPI_Comm *)malloc(grid_size_[0]*sizeof(MPI_Comm));

	dim0_group_ = (MPI_Group *)malloc(grid_size_[1]*sizeof(MPI_Group));
	dim1_group_ = (MPI_Group *)malloc(grid_size_[0]*sizeof(MPI_Group));

    // to create MPI_Group
    // one range is a triplet: (first rank, last rank, stride)
    int range[3];
    
    int comm_rank;



#ifdef EXTERNAL_IO

    if(world_rank_==0)
	{
		if(proc_size0*proc_size1+IO_total_size!=world_size_)
		{
			cerr<<"Latfield::Parallel2d::initialization - wrong number of process"<<endl;
			cerr<<"Latfield::Parallel2d::initialization - Number of total process must be equal to proc_size0*proc_size1+IO_total_size"<<endl;
			cerr<<"Latfield::Parallel2d::initialization - Within the call : Parallel2d(int proc_size0, int proc_size1, int IO_total_size)"<<endl;
			this->abortForce();
		}



	}

    MPI_Comm_group(world_comm_,&world_group_);

    range[0]=0;
    range[1]=proc_size0*proc_size1-1;
    range[2]=1;

    MPI_Group_range_incl(world_group_,1,&range,&lat_world_group_);
    MPI_Comm_create(world_comm_,lat_world_group_ , &lat_world_comm_);

    range[0]=proc_size0*proc_size1;
    range[1]=proc_size0*proc_size1 + IO_total_size - 1;
    range[2]=1;

    MPI_Group_range_incl(world_group_,1,&range,&IO_group_);
    MPI_Comm_create(world_comm_,IO_group_ , &IO_comm_);


    MPI_Group_rank(lat_world_group_, &comm_rank);
    if(comm_rank!=MPI_UNDEFINED)
    {
        lat_world_rank_=comm_rank;
        MPI_Comm_size( lat_world_comm_, &lat_world_size_ );

        range[2]=1;
        for(int j=0;j<grid_size_[1];j++)
        {
            range[0] = j * grid_size_[0];
            range[1] = grid_size_[0] - 1 + j*grid_size_[0];
            MPI_Group_range_incl(lat_world_group_,1,&range,&dim0_group_[j]);
            MPI_Comm_create(lat_world_comm_, dim0_group_[j], &dim0_comm_[j]);
        }


        range[2]=grid_size_[0];
        for(int i=0;i<grid_size_[0];i++)
        {
            range[0]=i;
            range[1]=i+(grid_size_[1]-1)*grid_size_[0];
            MPI_Group_range_incl(lat_world_group_,1,&range,&dim1_group_[i]);
            MPI_Comm_create(lat_world_comm_, dim1_group_[i], &dim1_comm_[i]);
        }


        for(int i=0;i<grid_size_[0];i++)
        {
            MPI_Group_rank(dim1_group_[i], &comm_rank);
            if(comm_rank!=MPI_UNDEFINED)grid_rank_[1]=comm_rank;
        }

        for(int j=0;j<grid_size_[1];j++)
        {
            MPI_Group_rank(dim0_group_[j], &comm_rank);
            if(comm_rank!=MPI_UNDEFINED)grid_rank_[0]=comm_rank;
        }


        isIO_=false;
    }
    else
    {
        lat_world_rank_=-1;
        grid_rank_[1]=-1;
        grid_rank_[0]=-1;
        isIO_=true;
    }

    if(grid_rank_[0]==grid_size_[0]-1)last_proc_[0]=true;
    else last_proc_[0]=false;
    if(grid_rank_[1]==grid_size_[1]-1)last_proc_[1]=true;
    else last_proc_[1]=false;




    ioserver.initialize(proc_size0,proc_size1,IO_total_size,IO_node_size);

#else

	if(lat_world_rank_==0)
	{
		if(proc_size0*proc_size1!=lat_world_size_)
		{
			cerr<<"Latfield::Parallel2d::initialization - wrong number of process"<<endl;
			cerr<<"Latfield::Parallel2d::initialization - Number of total process must be equal to proc_size0*proc_size1"<<endl;
			cerr<<"Latfield::Parallel2d::initialization - Within the call : Parallel2d(int proc_size0, int proc_size1)"<<endl;
			this->abortForce();
		}
	}
    
    // create the global group
    MPI_Comm_group(lat_world_comm_,&lat_world_group_);


    // create the groups and communicators along the X direction
	range[2]=1;
	for(int j=0;j<grid_size_[1];j++)
	{
		range[0] = j * grid_size_[0];
		range[1] = grid_size_[0] - 1 + j*grid_size_[0];
		MPI_Group_range_incl(lat_world_group_,1,&range,&dim0_group_[j]);
		MPI_Comm_create(lat_world_comm_, dim0_group_[j], &dim0_comm_[j]);
	}


    // create the groups and communicators along the Y direction
	range[2]=grid_size_[0];
	for(int i=0;i<grid_size_[0];i++)
	{
		range[0]=i;
		range[1]=i+(grid_size_[1]-1)*grid_size_[0];
		MPI_Group_range_incl(lat_world_group_,1,&range,&dim1_group_[i]);
		MPI_Comm_create(lat_world_comm_, dim1_group_[i], &dim1_comm_[i]);
	}



	for(int i=0;i<grid_size_[0];i++)
	{
		MPI_Group_rank(dim1_group_[i], &comm_rank);
		if(comm_rank!=MPI_UNDEFINED)grid_rank_[1]=comm_rank;
	}

	for(int j=0;j<grid_size_[1];j++)
	{
		MPI_Group_rank(dim0_group_[j], &comm_rank);
		if(comm_rank!=MPI_UNDEFINED)grid_rank_[0]=comm_rank;
	}



    if(grid_rank_[0]==grid_size_[0]-1)last_proc_[0]=true;
	else last_proc_[0]=false;
	if(grid_rank_[1]==grid_size_[1]-1)last_proc_[1]=true;
	else last_proc_[1]=false;



#endif

   
   
    /*
        The use of Boost::MPI cartesian_communicator defeats the
        purpose of grid_rank_[], dim0_comm_[], dim1_comm_[],
        dim0_group_[], dim1_group_[].
        I will try to remove those in the future.
        For the moment they will coexist with boost
        but we must make sure they're compatible.
    */
    
    assert(my_axis_communicator[0].rank()==grid_rank_[0]);
    assert(my_axis_communicator[1].rank()==grid_rank_[1]);
    
    assert(my_axis_communicator[0].size()==grid_size_[0]);
    assert(my_axis_communicator[1].size()==grid_size_[1]);
   
    {
        const auto rank = my_cartesian_communicator->rank();
        assert(rank==world_rank_);
        
        const auto coord = my_cartesian_communicator->coordinates(rank);
        assert(coord[1]==grid_rank_[0]);
        assert(coord[0]==grid_rank_[1]);
    }
}


Parallel2d::~Parallel2d()
{
	int finalized;
  MPI_Finalized(&finalized);
  if((!finalized) && (!neverFinalizeMPI)) { 
    // MPI_Finalize();
  }
}

//ABORT AND BARRIER===============================

void Parallel2d::abortForce()
{
	MPI_Abort( world_comm_, EXIT_FAILURE);
}

void Parallel2d::abortRequest()
{
	char failure;
	if(isRoot())
	{
		failure=char(1);
		broadcast(failure, root_);
		exit(EXIT_FAILURE);
	}
	else
	{
		cout<<"Parallel::abortRequest() called from non-Root process."<<endl;
		cout<<"Parallel::abortRequest() calling abortForce()..."<<endl;
		abortForce();
	}
}

void Parallel2d::barrier()
{

    MPI_Barrier( lat_world_comm_ );
}



//specializations: using MPI_ALLREDUCE...
template<> void Parallel2d::sum<double>(double* array, int len)
{
	MPI_Allreduce(MPI_IN_PLACE,array,len, MPI_DOUBLE,MPI_SUM,lat_world_comm_);
}
template<> void Parallel2d::sum<float>(float* array, int len)
{
	MPI_Allreduce(MPI_IN_PLACE,array,len, MPI_FLOAT,MPI_SUM,lat_world_comm_);
}
template<> void Parallel2d::sum<int>(int* array, int len)
{
	MPI_Allreduce(MPI_IN_PLACE,array,len, MPI_INT,MPI_SUM,lat_world_comm_);
}
template<> void Parallel2d::sum<long>(long* array, int len)
{
	MPI_Allreduce(MPI_IN_PLACE,array,len, MPI_LONG,MPI_SUM,lat_world_comm_);
}



//specializations: using MPI_ALLREDUCE...
template<> void Parallel2d::sum_dim0<double>(double* array, int len)
{
	MPI_Allreduce(MPI_IN_PLACE,array,len, MPI_DOUBLE,MPI_SUM,dim0_comm_[grid_rank_[1]]);
}
template<> void Parallel2d::sum_dim0<float>(float* array, int len)
{
	MPI_Allreduce(MPI_IN_PLACE,array,len, MPI_FLOAT,MPI_SUM,dim0_comm_[grid_rank_[1]]);
}
template<> void Parallel2d::sum_dim0<int>(int* array, int len)
{
	MPI_Allreduce(MPI_IN_PLACE,array,len, MPI_INT,MPI_SUM,dim0_comm_[grid_rank_[1]]);
}
template<> void Parallel2d::sum_dim0<long>(long* array, int len)
{
	MPI_Allreduce(MPI_IN_PLACE,array,len, MPI_LONG,MPI_SUM,dim0_comm_[grid_rank_[1]]);
}



//specializations: using MPI_ALLREDUCE...
template<> void Parallel2d::sum_dim1<double>(double* array, int len)
{
	MPI_Allreduce(MPI_IN_PLACE,array,len, MPI_DOUBLE,MPI_SUM,dim1_comm_[grid_rank_[0]]);
}
template<> void Parallel2d::sum_dim1<float>(float* array, int len)
{
	MPI_Allreduce(MPI_IN_PLACE,array,len, MPI_FLOAT,MPI_SUM,dim1_comm_[grid_rank_[0]]);
}
template<> void Parallel2d::sum_dim1<int>(int* array, int len)
{
	MPI_Allreduce(MPI_IN_PLACE,array,len, MPI_INT,MPI_SUM,dim1_comm_[grid_rank_[0]]);
}
template<> void Parallel2d::sum_dim1<long>(long* array, int len)
{
	MPI_Allreduce(MPI_IN_PLACE,array,len, MPI_LONG,MPI_SUM,dim1_comm_[grid_rank_[0]]);
}



//specializations using MPI_REDUCE:
template<>
void Parallel2d::sum_to<double>(double* array, int len, int dest)
{
	if(rank() == dest)
		MPI_Reduce(MPI_IN_PLACE,(void*)array,len,MPI_DOUBLE,MPI_SUM,dest,lat_world_comm_);
	else
		MPI_Reduce((void*)array,(void*)array,len,MPI_DOUBLE,MPI_SUM,dest,lat_world_comm_);
}
template<>
void Parallel2d::sum_to<float>(float* array, int len, int dest)
{
	if(rank() == dest)
		MPI_Reduce(MPI_IN_PLACE,(void*)array,len,MPI_FLOAT,MPI_SUM,dest,lat_world_comm_);
	else
		MPI_Reduce((void*)array,(void*)array,len,MPI_FLOAT,MPI_SUM,dest,lat_world_comm_);
}
template<>
void Parallel2d::sum_to<int>(int* array, int len, int dest)
{
	if(rank() == dest)
		MPI_Reduce(MPI_IN_PLACE,(void*)array,len,MPI_INT,MPI_SUM,dest,lat_world_comm_);
	else
		MPI_Reduce((void*)array,(void*)array,len,MPI_INT,MPI_SUM,dest,lat_world_comm_);
}
template<>
void Parallel2d::sum_to<long>(long* array, int len, int dest)
{
	if(rank() == dest)
		MPI_Reduce(MPI_IN_PLACE,(void*)array,len,MPI_LONG,MPI_SUM,dest,lat_world_comm_);
	else
		MPI_Reduce((void*)array,(void*)array,len,MPI_LONG,MPI_SUM,dest,lat_world_comm_);
}


//specializations using MPI_REDUCE:
template<>
void Parallel2d::sum_dim0_to<double>(double* array, int len, int dest)
{
	if(grid_rank_[0] == dest)
		MPI_Reduce(MPI_IN_PLACE,(void*)array,len,MPI_DOUBLE,MPI_SUM,dest,dim0_comm_[grid_rank_[1]]);
	else
		MPI_Reduce((void*)array,(void*)array,len,MPI_DOUBLE,MPI_SUM,dest,dim0_comm_[grid_rank_[1]]);
}
template<>
void Parallel2d::sum_dim0_to<float>(float* array, int len, int dest)
{
	if(grid_rank_[0] == dest)
		MPI_Reduce(MPI_IN_PLACE,(void*)array,len,MPI_FLOAT,MPI_SUM,dest,dim0_comm_[grid_rank_[1]]);
	else
		MPI_Reduce((void*)array,(void*)array,len,MPI_FLOAT,MPI_SUM,dest,dim0_comm_[grid_rank_[1]]);
}
template<>
void Parallel2d::sum_dim0_to<int>(int* array, int len, int dest)
{
	if(grid_rank_[0] == dest)
		MPI_Reduce(MPI_IN_PLACE,(void*)array,len,MPI_INT,MPI_SUM,dest,dim0_comm_[grid_rank_[1]]);
	else
		MPI_Reduce((void*)array,(void*)array,len,MPI_INT,MPI_SUM,dest,dim0_comm_[grid_rank_[1]]);
}
template<>
void Parallel2d::sum_dim0_to<long>(long* array, int len, int dest)
{
	if(grid_rank_[0] == dest)
		MPI_Reduce(MPI_IN_PLACE,(void*)array,len,MPI_LONG,MPI_SUM,dest,dim0_comm_[grid_rank_[1]]);
	else
		MPI_Reduce((void*)array,(void*)array,len,MPI_LONG,MPI_SUM,dest,dim0_comm_[grid_rank_[1]]);
}


//specializations using MPI_REDUCE:
template<>
void Parallel2d::sum_dim1_to<double>(double* array, int len, int dest)
{
	if(grid_rank_[1] == dest)
		MPI_Reduce(MPI_IN_PLACE,(void*)array,len,MPI_DOUBLE,MPI_SUM,dest,dim1_comm_[grid_rank_[0]]);
	else
		MPI_Reduce((void*)array,(void*)array,len,MPI_DOUBLE,MPI_SUM,dest,dim1_comm_[grid_rank_[0]]);
}
template<>
void Parallel2d::sum_dim1_to<float>(float* array, int len, int dest)
{
	if(grid_rank_[1] == dest)
		MPI_Reduce(MPI_IN_PLACE,(void*)array,len,MPI_FLOAT,MPI_SUM,dest,dim1_comm_[grid_rank_[0]]);
	else
		MPI_Reduce((void*)array,(void*)array,len,MPI_FLOAT,MPI_SUM,dest,dim1_comm_[grid_rank_[0]]);
}
template<>
void Parallel2d::sum_dim1_to<int>(int* array, int len, int dest)
{
	if(grid_rank_[1] == dest)
		MPI_Reduce(MPI_IN_PLACE,(void*)array,len,MPI_INT,MPI_SUM,dest,dim1_comm_[grid_rank_[0]]);
	else
		MPI_Reduce((void*)array,(void*)array,len,MPI_INT,MPI_SUM,dest,dim1_comm_[grid_rank_[0]]);
}
template<>
void Parallel2d::sum_dim1_to<long>(long* array, int len, int dest)
{
	if(grid_rank_[1] == dest)
		MPI_Reduce(MPI_IN_PLACE,(void*)array,len,MPI_LONG,MPI_SUM,dest,dim1_comm_[grid_rank_[0]]);
	else
		MPI_Reduce((void*)array,(void*)array,len,MPI_LONG,MPI_SUM,dest,dim1_comm_[grid_rank_[0]]);
}


template<> void Parallel2d::max<double>(double* array, int len)
{
	MPI_Allreduce(MPI_IN_PLACE,array,len, MPI_DOUBLE,MPI_MAX,lat_world_comm_);
}
template<> void Parallel2d::max<float>(float* array, int len)
{
	MPI_Allreduce(MPI_IN_PLACE,array,len, MPI_FLOAT,MPI_MAX,lat_world_comm_);
}
template<> void Parallel2d::max<int>(int* array, int len)
{
	MPI_Allreduce(MPI_IN_PLACE,array,len, MPI_INT,MPI_MAX,lat_world_comm_);
}
template<> void Parallel2d::max<long>(long* array, int len)
{
	MPI_Allreduce(MPI_IN_PLACE,array,len, MPI_LONG,MPI_MAX,lat_world_comm_);
}


//specializations: using MPI_ALLREDUCE...
template<> void Parallel2d::max_dim0<double>(double* array, int len)
{
	MPI_Allreduce(MPI_IN_PLACE,array,len, MPI_DOUBLE,MPI_MAX,dim0_comm_[grid_rank_[1]]);
}
template<> void Parallel2d::max_dim0<float>(float* array, int len)
{
	MPI_Allreduce(MPI_IN_PLACE,array,len, MPI_FLOAT,MPI_MAX,dim0_comm_[grid_rank_[1]]);
}
template<> void Parallel2d::max_dim0<int>(int* array, int len)
{
	MPI_Allreduce(MPI_IN_PLACE,array,len, MPI_INT,MPI_MAX,dim0_comm_[grid_rank_[1]]);
}
template<> void Parallel2d::max_dim0<long>(long* array, int len)
{
	MPI_Allreduce(MPI_IN_PLACE,array,len, MPI_LONG,MPI_MAX,dim0_comm_[grid_rank_[1]]);
}


//specializations: using MPI_ALLREDUCE...
template<> void Parallel2d::max_dim1<double>(double* array, int len)
{
	MPI_Allreduce(MPI_IN_PLACE,array,len, MPI_DOUBLE,MPI_MAX,dim1_comm_[grid_rank_[0]]);
}
template<> void Parallel2d::max_dim1<float>(float* array, int len)
{
	MPI_Allreduce(MPI_IN_PLACE,array,len, MPI_FLOAT,MPI_MAX,dim1_comm_[grid_rank_[0]]);
}
template<> void Parallel2d::max_dim1<int>(int* array, int len)
{
	MPI_Allreduce(MPI_IN_PLACE,array,len, MPI_INT,MPI_MAX,dim1_comm_[grid_rank_[0]]);
}
template<> void Parallel2d::max_dim1<long>(long* array, int len)
{
	MPI_Allreduce(MPI_IN_PLACE,array,len, MPI_LONG,MPI_MAX,dim1_comm_[grid_rank_[0]]);
}



//specializations using MPI_REDUCE:
template<>
void Parallel2d::max_to<double>(double* array, int len, int dest)
{
	if(rank() == dest)
		MPI_Reduce(MPI_IN_PLACE,(void*)array,len,MPI_DOUBLE,MPI_MAX,dest,lat_world_comm_);
	else
		MPI_Reduce((void*)array,(void*)array,len,MPI_DOUBLE,MPI_MAX,dest,lat_world_comm_);
}
template<>
void Parallel2d::max_to<float>(float* array, int len, int dest)
{
	if(rank() == dest)
		MPI_Reduce(MPI_IN_PLACE,(void*)array,len,MPI_FLOAT,MPI_MAX,dest,lat_world_comm_);
	else
		MPI_Reduce((void*)array,(void*)array,len,MPI_FLOAT,MPI_MAX,dest,lat_world_comm_);
}
template<>
void Parallel2d::max_to<int>(int* array, int len, int dest)
{
	if(rank() == dest)
		MPI_Reduce(MPI_IN_PLACE,(void*)array,len,MPI_INT,MPI_MAX,dest,lat_world_comm_);
	else
		MPI_Reduce((void*)array,(void*)array,len,MPI_INT,MPI_MAX,dest,lat_world_comm_);
}
template<>
void Parallel2d::max_to<long>(long* array, int len, int dest)
{
	if(rank() == dest)
		MPI_Reduce(MPI_IN_PLACE,(void*)array,len,MPI_LONG,MPI_MAX,dest,lat_world_comm_);
	else
		MPI_Reduce((void*)array,(void*)array,len,MPI_LONG,MPI_MAX,dest,lat_world_comm_);
}



//specializations using MPI_REDUCE:
template<>
void Parallel2d::max_dim0_to<double>(double* array, int len, int dest)
{
	if(grid_rank_[0] == dest)
		MPI_Reduce(MPI_IN_PLACE,(void*)array,len,MPI_DOUBLE,MPI_MAX,dest,dim0_comm_[grid_rank_[1]]);
	else
		MPI_Reduce((void*)array,(void*)array,len,MPI_DOUBLE,MPI_MAX,dest,dim0_comm_[grid_rank_[1]]);
}
template<>
void Parallel2d::max_dim0_to<float>(float* array, int len, int dest)
{
	if(grid_rank_[0] == dest)
		MPI_Reduce(MPI_IN_PLACE,(void*)array,len,MPI_FLOAT,MPI_MAX,dest,dim0_comm_[grid_rank_[1]]);
	else
		MPI_Reduce((void*)array,(void*)array,len,MPI_FLOAT,MPI_MAX,dest,dim0_comm_[grid_rank_[1]]);
}
template<>
void Parallel2d::max_dim0_to<int>(int* array, int len, int dest)
{
	if(grid_rank_[0] == dest)
		MPI_Reduce(MPI_IN_PLACE,(void*)array,len,MPI_INT,MPI_MAX,dest,dim0_comm_[grid_rank_[1]]);
	else
		MPI_Reduce((void*)array,(void*)array,len,MPI_INT,MPI_MAX,dest,dim0_comm_[grid_rank_[1]]);
}
template<>
void Parallel2d::max_dim0_to<long>(long* array, int len, int dest)
{
	if(grid_rank_[0] == dest)
		MPI_Reduce(MPI_IN_PLACE,(void*)array,len,MPI_LONG,MPI_MAX,dest,dim0_comm_[grid_rank_[1]]);
	else
		MPI_Reduce((void*)array,(void*)array,len,MPI_LONG,MPI_MAX,dest,dim0_comm_[grid_rank_[1]]);
}



//specializations using MPI_REDUCE:
template<>
void Parallel2d::max_dim1_to<double>(double* array, int len, int dest)
{
	if(grid_rank_[1] == dest)
		MPI_Reduce(MPI_IN_PLACE,(void*)array,len,MPI_DOUBLE,MPI_MAX,dest,dim1_comm_[grid_rank_[0]]);
	else
		MPI_Reduce((void*)array,(void*)array,len,MPI_DOUBLE,MPI_MAX,dest,dim1_comm_[grid_rank_[0]]);
}
template<>
void Parallel2d::max_dim1_to<float>(float* array, int len, int dest)
{
	if(grid_rank_[1] == dest)
		MPI_Reduce(MPI_IN_PLACE,(void*)array,len,MPI_FLOAT,MPI_MAX,dest,dim1_comm_[grid_rank_[0]]);
	else
		MPI_Reduce((void*)array,(void*)array,len,MPI_FLOAT,MPI_MAX,dest,dim1_comm_[grid_rank_[0]]);
}
template<>
void Parallel2d::max_dim1_to<int>(int* array, int len, int dest)
{
	if(grid_rank_[1] == dest)
		MPI_Reduce(MPI_IN_PLACE,(void*)array,len,MPI_INT,MPI_MAX,dest,dim1_comm_[grid_rank_[0]]);
	else
		MPI_Reduce((void*)array,(void*)array,len,MPI_INT,MPI_MAX,dest,dim1_comm_[grid_rank_[0]]);
}
template<>
void Parallel2d::max_dim1_to<long>(long* array, int len, int dest)
{
	if(grid_rank_[1] == dest)
		MPI_Reduce(MPI_IN_PLACE,(void*)array,len,MPI_LONG,MPI_MAX,dest,dim1_comm_[grid_rank_[0]]);
	else
		MPI_Reduce((void*)array,(void*)array,len,MPI_LONG,MPI_MAX,dest,dim1_comm_[grid_rank_[0]]);
}


template<> void Parallel2d::min<double>(double* array, int len)
{
	MPI_Allreduce(MPI_IN_PLACE,array,len, MPI_DOUBLE,MPI_MIN,lat_world_comm_);
}
template<> void Parallel2d::min<float>(float* array, int len)
{
	MPI_Allreduce(MPI_IN_PLACE,array,len, MPI_FLOAT,MPI_MIN,lat_world_comm_);
}
template<> void Parallel2d::min<int>(int* array, int len)
{
	MPI_Allreduce(MPI_IN_PLACE,array,len, MPI_INT,MPI_MIN,lat_world_comm_);
}
template<> void Parallel2d::min<long>(long* array, int len)
{
	MPI_Allreduce(MPI_IN_PLACE,array,len, MPI_LONG,MPI_MIN,lat_world_comm_);
}


//specializations: using MPI_ALLREDUCE...
template<> void Parallel2d::min_dim0<double>(double* array, int len)
{
	MPI_Allreduce(MPI_IN_PLACE,array,len, MPI_DOUBLE,MPI_MIN,dim0_comm_[grid_rank_[1]]);
}
template<> void Parallel2d::min_dim0<float>(float* array, int len)
{
	MPI_Allreduce(MPI_IN_PLACE,array,len, MPI_FLOAT,MPI_MIN,dim0_comm_[grid_rank_[1]]);
}
template<> void Parallel2d::min_dim0<int>(int* array, int len)
{
	MPI_Allreduce(MPI_IN_PLACE,array,len, MPI_INT,MPI_MIN,dim0_comm_[grid_rank_[1]]);
}
template<> void Parallel2d::min_dim0<long>(long* array, int len)
{
	MPI_Allreduce(MPI_IN_PLACE,array,len, MPI_LONG,MPI_MIN,dim0_comm_[grid_rank_[1]]);
}


//specializations: using MPI_ALLREDUCE...
template<> void Parallel2d::min_dim1<double>(double* array, int len)
{
	MPI_Allreduce(MPI_IN_PLACE,array,len, MPI_DOUBLE,MPI_MIN,dim1_comm_[grid_rank_[0]]);
}
template<> void Parallel2d::min_dim1<float>(float* array, int len)
{
	MPI_Allreduce(MPI_IN_PLACE,array,len, MPI_FLOAT,MPI_MIN,dim1_comm_[grid_rank_[0]]);
}
template<> void Parallel2d::min_dim1<int>(int* array, int len)
{
	MPI_Allreduce(MPI_IN_PLACE,array,len, MPI_INT,MPI_MIN,dim1_comm_[grid_rank_[0]]);
}
template<> void Parallel2d::min_dim1<long>(long* array, int len)
{
	MPI_Allreduce(MPI_IN_PLACE,array,len, MPI_LONG,MPI_MIN,dim1_comm_[grid_rank_[0]]);
}


//specializations using MPI_REDUCE:
template<>
void Parallel2d::min_to<double>(double* array, int len, int dest)
{
	if(rank() == dest)
		MPI_Reduce(MPI_IN_PLACE,(void*)array,len,MPI_DOUBLE,MPI_MIN,dest,lat_world_comm_);
	else
		MPI_Reduce((void*)array,(void*)array,len,MPI_DOUBLE,MPI_MIN,dest,lat_world_comm_);
}
template<>
void Parallel2d::min_to<float>(float* array, int len, int dest)
{
	if(rank() == dest)
		MPI_Reduce(MPI_IN_PLACE,(void*)array,len,MPI_FLOAT,MPI_MIN,dest,lat_world_comm_);
	else
		MPI_Reduce((void*)array,(void*)array,len,MPI_FLOAT,MPI_MIN,dest,lat_world_comm_);
}
template<>
void Parallel2d::min_to<int>(int* array, int len, int dest)
{
	if(rank() == dest)
		MPI_Reduce(MPI_IN_PLACE,(void*)array,len,MPI_INT,MPI_MIN,dest,lat_world_comm_);
	else
		MPI_Reduce((void*)array,(void*)array,len,MPI_INT,MPI_MIN,dest,lat_world_comm_);
}
template<>
void Parallel2d::min_to<long>(long* array, int len, int dest)
{
	if(rank() == dest)
		MPI_Reduce(MPI_IN_PLACE,(void*)array,len,MPI_LONG,MPI_MIN,dest,lat_world_comm_);
	else
		MPI_Reduce((void*)array,(void*)array,len,MPI_LONG,MPI_MIN,dest,lat_world_comm_);
}



//specializations using MPI_REDUCE:
template<>
void Parallel2d::min_dim0_to<double>(double* array, int len, int dest)
{
	if(grid_rank_[0] == dest)
		MPI_Reduce(MPI_IN_PLACE,(void*)array,len,MPI_DOUBLE,MPI_MIN,dest,dim0_comm_[grid_rank_[1]]);
	else
		MPI_Reduce((void*)array,(void*)array,len,MPI_DOUBLE,MPI_MIN,dest,dim0_comm_[grid_rank_[1]]);
}
template<>
void Parallel2d::min_dim0_to<float>(float* array, int len, int dest)
{
	if(grid_rank_[0] == dest)
		MPI_Reduce(MPI_IN_PLACE,(void*)array,len,MPI_FLOAT,MPI_MIN,dest,dim0_comm_[grid_rank_[1]]);
	else
		MPI_Reduce((void*)array,(void*)array,len,MPI_FLOAT,MPI_MIN,dest,dim0_comm_[grid_rank_[1]]);
}
template<>
void Parallel2d::min_dim0_to<int>(int* array, int len, int dest)
{
	if(grid_rank_[0] == dest)
		MPI_Reduce(MPI_IN_PLACE,(void*)array,len,MPI_INT,MPI_MIN,dest,dim0_comm_[grid_rank_[1]]);
	else
		MPI_Reduce((void*)array,(void*)array,len,MPI_INT,MPI_MIN,dest,dim0_comm_[grid_rank_[1]]);
}
template<>
void Parallel2d::min_dim0_to<long>(long* array, int len, int dest)
{
	if(grid_rank_[0] == dest)
		MPI_Reduce(MPI_IN_PLACE,(void*)array,len,MPI_LONG,MPI_MIN,dest,dim0_comm_[grid_rank_[1]]);
	else
		MPI_Reduce((void*)array,(void*)array,len,MPI_LONG,MPI_MIN,dest,dim0_comm_[grid_rank_[1]]);
}


//specializations using MPI_REDUCE:
template<>
void Parallel2d::min_dim1_to<double>(double* array, int len, int dest)
{
	if(grid_rank_[1] == dest)
		MPI_Reduce(MPI_IN_PLACE,(void*)array,len,MPI_DOUBLE,MPI_MIN,dest,dim1_comm_[grid_rank_[0]]);
	else
		MPI_Reduce((void*)array,(void*)array,len,MPI_DOUBLE,MPI_MIN,dest,dim1_comm_[grid_rank_[0]]);
}
template<>
void Parallel2d::min_dim1_to<float>(float* array, int len, int dest)
{
	if(grid_rank_[1] == dest)
		MPI_Reduce(MPI_IN_PLACE,(void*)array,len,MPI_FLOAT,MPI_MIN,dest,dim1_comm_[grid_rank_[0]]);
	else
		MPI_Reduce((void*)array,(void*)array,len,MPI_FLOAT,MPI_MIN,dest,dim1_comm_[grid_rank_[0]]);
}
template<>
void Parallel2d::min_dim1_to<int>(int* array, int len, int dest)
{
	if(grid_rank_[1] == dest)
		MPI_Reduce(MPI_IN_PLACE,(void*)array,len,MPI_INT,MPI_MIN,dest,dim1_comm_[grid_rank_[0]]);
	else
		MPI_Reduce((void*)array,(void*)array,len,MPI_INT,MPI_MIN,dest,dim1_comm_[grid_rank_[0]]);
}
template<>
void Parallel2d::min_dim1_to<long>(long* array, int len, int dest)
{
	if(grid_rank_[1] == dest)
		MPI_Reduce(MPI_IN_PLACE,(void*)array,len,MPI_LONG,MPI_MIN,dest,dim1_comm_[grid_rank_[0]]);
	else
		MPI_Reduce((void*)array,(void*)array,len,MPI_LONG,MPI_MIN,dest,dim1_comm_[grid_rank_[0]]);
}


}

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
{
    #ifndef EXTERNAL_IO
    lat_world_comm_ = com;
    world_comm_ = com;
    MPI_Comm_rank( lat_world_comm_, &lat_world_rank_ );
    MPI_Comm_size( lat_world_comm_, &lat_world_size_ );
    MPI_Comm_rank( world_comm_, &world_rank_ );
    MPI_Comm_size( world_comm_, &world_size_ );
    #else
    world_comm_ = com;
    MPI_Comm_rank( world_comm_, &world_rank_ );
    MPI_Comm_size( world_comm_, &world_size_ );
    #endif

    grid_size_[0]=proc_size0;
	grid_size_[1]=proc_size1;

	dim0_comm_ = (MPI_Comm *)malloc(grid_size_[1]*sizeof(MPI_Comm));
	dim1_comm_ = (MPI_Comm *)malloc(grid_size_[0]*sizeof(MPI_Comm));

	dim0_group_ = (MPI_Group *)malloc(grid_size_[1]*sizeof(MPI_Group));
	dim1_group_ = (MPI_Group *)malloc(grid_size_[0]*sizeof(MPI_Group));

	int rang[3],i,j,comm_rank;



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

    rang[0]=0;
    rang[1]=proc_size0*proc_size1-1;
    rang[2]=1;

    MPI_Group_range_incl(world_group_,1,&rang,&lat_world_group_);
    MPI_Comm_create(world_comm_,lat_world_group_ , &lat_world_comm_);

    rang[0]=proc_size0*proc_size1;
    rang[1]=proc_size0*proc_size1 + IO_total_size - 1;
    rang[2]=1;

    MPI_Group_range_incl(world_group_,1,&rang,&IO_group_);
    MPI_Comm_create(world_comm_,IO_group_ , &IO_comm_);


    MPI_Group_rank(lat_world_group_, &comm_rank);
    if(comm_rank!=MPI_UNDEFINED)
    {
        lat_world_rank_=comm_rank;
        MPI_Comm_size( lat_world_comm_, &lat_world_size_ );

        rang[2]=1;
        for(j=0;j<grid_size_[1];j++)
        {
            rang[0] = j * grid_size_[0];
            rang[1] = grid_size_[0] - 1 + j*grid_size_[0];
            MPI_Group_range_incl(lat_world_group_,1,&rang,&dim0_group_[j]);
            MPI_Comm_create(lat_world_comm_, dim0_group_[j], &dim0_comm_[j]);
        }


        rang[2]=grid_size_[0];
        for(i=0;i<grid_size_[0];i++)
        {
            rang[0]=i;
            rang[1]=i+(grid_size_[1]-1)*grid_size_[0];
            MPI_Group_range_incl(lat_world_group_,1,&rang,&dim1_group_[i]);
            MPI_Comm_create(lat_world_comm_, dim1_group_[i], &dim1_comm_[i]);
        }


        for(i=0;i<grid_size_[0];i++)
        {
            MPI_Group_rank(dim1_group_[i], &comm_rank);
            if(comm_rank!=MPI_UNDEFINED)grid_rank_[1]=comm_rank;
        }

        for(j=0;j<grid_size_[1];j++)
        {
            MPI_Group_rank(dim0_group_[j], &comm_rank);
            if(comm_rank!=MPI_UNDEFINED)grid_rank_[0]=comm_rank;
        }


        root_=0;
        isIO_=false;
    }
    else
    {
        lat_world_rank_=-1;
        root_=0;
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

    MPI_Comm_group(lat_world_comm_,&lat_world_group_);




	rang[2]=1;
	for(j=0;j<grid_size_[1];j++)
	{
		rang[0] = j * grid_size_[0];
		rang[1] = grid_size_[0] - 1 + j*grid_size_[0];
		MPI_Group_range_incl(lat_world_group_,1,&rang,&dim0_group_[j]);
		MPI_Comm_create(lat_world_comm_, dim0_group_[j], &dim0_comm_[j]);
	}


	rang[2]=grid_size_[0];
	for(i=0;i<grid_size_[0];i++)
	{
		rang[0]=i;
		rang[1]=i+(grid_size_[1]-1)*grid_size_[0];
		MPI_Group_range_incl(lat_world_group_,1,&rang,&dim1_group_[i]);
		MPI_Comm_create(lat_world_comm_, dim1_group_[i], &dim1_comm_[i]);
	}



	for(i=0;i<grid_size_[0];i++)
	{
		MPI_Group_rank(dim1_group_[i], &comm_rank);
		if(comm_rank!=MPI_UNDEFINED)grid_rank_[1]=comm_rank;
	}

	for(j=0;j<grid_size_[1];j++)
	{
		MPI_Group_rank(dim0_group_[j], &comm_rank);
		if(comm_rank!=MPI_UNDEFINED)grid_rank_[0]=comm_rank;
	}


	root_=0;

    if(grid_rank_[0]==grid_size_[0]-1)last_proc_[0]=true;
	else last_proc_[0]=false;
	if(grid_rank_[1]==grid_size_[1]-1)last_proc_[1]=true;
	else last_proc_[1]=false;



#endif

	if(root_ == lat_world_rank_)isRoot_=true;
    else isRoot_=false;


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

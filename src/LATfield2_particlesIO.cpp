#include "config.h"
#ifdef HDF5
#include "particles/LATfield2_particlesIO.hpp"

namespace LATfield2
{


void get_fileDsc_global(std::string filename,fileDsc &fd)
{

  #ifdef H5_HAVE_PARALLEL

  hid_t plist_id,file_id,attr_id;//root_id;


  MPI_Info info  = MPI_INFO_NULL;


  plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id,parallel.lat_world_comm(),info);
  file_id = H5Fopen(filename.c_str(),H5F_ACC_RDONLY,plist_id);

  if(file_id<0)std::cout<< "get_fileDsc_global: cant open file: "<<filename<<std::endl;

  if(parallel.rank()==0)
  {
  	attr_id = H5Aopen_name(file_id, "fileNumber");
  	H5Aread(attr_id, H5T_NATIVE_INT, &(fd.fileNumber));
  	H5Aclose(attr_id);

          attr_id = H5Aopen_name(file_id, "numProcPerFile");
  	H5Aread(attr_id, H5T_NATIVE_INT, &(fd.numProcPerFile));
  	H5Aclose(attr_id);

  	attr_id = H5Aopen_name(file_id, "world_size");
  	H5Aread(attr_id, H5T_NATIVE_INT, &(fd.world_size));
  	H5Aclose(attr_id);

  	attr_id = H5Aopen_name(file_id, "grid_size");
  	H5Aread(attr_id, H5T_NATIVE_INT, fd.grid_size);
  	H5Aclose(attr_id);

  	attr_id = H5Aopen_name(file_id, "boxSize");
  	H5Aread(attr_id, REAL_TYPE, fd.boxSize);
  	H5Aclose(attr_id);



  	attr_id = H5Aopen_name(file_id, "fileBoxSize");
  	H5Aread(attr_id, REAL_TYPE, &(fd.fileBoxSize));
  	H5Aclose(attr_id);

  	attr_id = H5Aopen_name(file_id, "fileBoxOffset");
  	H5Aread(attr_id, REAL_TYPE, &(fd.fileBoxOffset));
  	H5Aclose(attr_id);
  }
  H5Pclose(plist_id);
	H5Fclose(file_id);


  #else

    if(parallel.rank()==0)
    {
        hid_t plist_id,file_id,attr_id,root_id;

        plist_id = H5Pcreate(H5P_FILE_ACCESS);
        file_id = H5Fopen(filename.c_str(),H5F_ACC_RDONLY,plist_id);




	attr_id = H5Aopen_name(file_id, "fileNumber");
	H5Aread(attr_id, H5T_NATIVE_INT, &(fd.fileNumber));
	H5Aclose(attr_id);

        attr_id = H5Aopen_name(file_id, "numProcPerFile");
	H5Aread(attr_id, H5T_NATIVE_INT, &(fd.numProcPerFile));
	H5Aclose(attr_id);

	attr_id = H5Aopen_name(file_id, "world_size");
	H5Aread(attr_id, H5T_NATIVE_INT, &(fd.world_size));
	H5Aclose(attr_id);

	attr_id = H5Aopen_name(file_id, "grid_size");
	H5Aread(attr_id, H5T_NATIVE_INT, fd.grid_size);
	H5Aclose(attr_id);

	attr_id = H5Aopen_name(file_id, "boxSize");
	H5Aread(attr_id, REAL_TYPE, fd.boxSize);
	H5Aclose(attr_id);



	attr_id = H5Aopen_name(file_id, "fileBoxSize");
	H5Aread(attr_id, REAL_TYPE, &(fd.fileBoxSize));
	H5Aclose(attr_id);

	attr_id = H5Aopen_name(file_id, "fileBoxOffset");
	H5Aread(attr_id, REAL_TYPE, &(fd.fileBoxOffset));
	H5Aclose(attr_id);
	//cout<<"bosize loaded is "<< fd.boxSize[0] <<endl;
	H5Pclose(plist_id);
	H5Fclose(file_id);
    }

#endif
  parallel.broadcast<fileDsc>(fd,0);
  //cout<<"bosize loaded is "<< fd.boxSize[0] <<endl;


}
void get_fileDsc_local(std::string filename,long * numParts, RealC * localBoxOffset, RealC * localBoxSize, int numProcPerfile)
{
  #ifdef H5_HAVE_PARALLEL

  hid_t plist_id,file_id,dataset_id;//attr_id,root_id;


  MPI_Info info  = MPI_INFO_NULL;


  plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id,parallel.lat_world_comm(),info);
  file_id = H5Fopen(filename.c_str(),H5F_ACC_RDONLY,plist_id);
  if(file_id<0)std::cout<< "get_fileDsc_local: cant open file:  "<<filename<<std::endl;

  if(parallel.rank()==0)
  {
    dataset_id = H5Dopen(file_id, "/numParts", H5P_DEFAULT);
  	H5Dread(dataset_id,H5T_NATIVE_LONG , H5S_ALL, H5S_ALL, H5P_DEFAULT,numParts);
  	H5Dclose(dataset_id);
  	dataset_id = H5Dopen(file_id, "/localBoxOffset", H5P_DEFAULT);
  	H5Dread(dataset_id,REAL_TYPE , H5S_ALL, H5S_ALL, H5P_DEFAULT,localBoxOffset);
  	H5Dclose(dataset_id);
  	dataset_id = H5Dopen(file_id, "/localBoxSize", H5P_DEFAULT);
  	H5Dread(dataset_id,REAL_TYPE , H5S_ALL, H5S_ALL, H5P_DEFAULT,localBoxSize);
  	H5Dclose(dataset_id);

  }
  H5Pclose(plist_id);
  H5Fclose(file_id);
  #else
    if(parallel.rank()==0)
    {
	hid_t plist_id,file_id,dataset_id;

	plist_id = H5Pcreate(H5P_FILE_ACCESS);
	file_id = H5Fopen(filename.c_str(),H5F_ACC_RDONLY,plist_id);

  if(file_id<0)std::cout<< "get_fileDsc_local: cant open file:  "<<filename<<std::endl;

	dataset_id = H5Dopen(file_id, "/numParts", H5P_DEFAULT);
	H5Dread(dataset_id,H5T_NATIVE_LONG , H5S_ALL, H5S_ALL, H5P_DEFAULT,numParts);
	H5Dclose(dataset_id);
	dataset_id = H5Dopen(file_id, "/localBoxOffset", H5P_DEFAULT);
	H5Dread(dataset_id,REAL_TYPE , H5S_ALL, H5S_ALL, H5P_DEFAULT,localBoxOffset);
	H5Dclose(dataset_id);
	dataset_id = H5Dopen(file_id, "/localBoxSize", H5P_DEFAULT);
	H5Dread(dataset_id,REAL_TYPE , H5S_ALL, H5S_ALL, H5P_DEFAULT,localBoxSize);
	H5Dclose(dataset_id);
	H5Pclose(plist_id);
	H5Fclose(file_id);
    }

    #endif
    parallel.broadcast(numParts,numProcPerfile,0);
    parallel.broadcast(localBoxOffset,3*numProcPerfile,0);
    parallel.broadcast(localBoxSize,3*numProcPerfile,0);

}



}
#endif

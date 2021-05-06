#pragma once
/*! \file LATfield2_save_hdf5.h
 \brief LATfield2_save_hdf5.h contains the definition of the function used for hdf5 i/o.
 \author David Daverio
 */

#include "LATfield2_macros.hpp"
#include <cmath>
#include <hdf5.h>
#include <cstdlib>
#include <cstring>
#include <string>



namespace LATfield2
{

   int save_hdf5_externC(char *data,long file_offset[2],int *size,int *
   sizeLocal,int halo, int lat_dim,int comp,hid_t array_type,int
   array_size,std::string  filename_str, std::string dataset_name_str);


	int load_hdf5_externC(char *data,long file_offset[2],int *size,int *
    sizeLocal,int halo, int lat_dim,std::string  filename_str, std::string
    dataset_name_str);

template<class fieldType>
int save_hdf5(fieldType *data,hid_t type_id,int array_size,long
file_offset[2],int *size,int * sizeLocal,int halo, int lat_dim,int
comp,std::string  filename_str, std::string dataset_name_str)
{

	return save_hdf5_externC((char*)data, file_offset, size, sizeLocal, halo, lat_dim, comp, type_id, array_size, filename_str, dataset_name_str);
}
template<class fieldType>
int load_hdf5(fieldType *data,long file_offset[2],int *size,int * sizeLocal,int
halo, int lat_dim,int comp,std::string  filename_str, std::string dataset_name_str)
{
    return load_hdf5_externC( (char*) data, file_offset, size, sizeLocal, halo, lat_dim,  filename_str, dataset_name_str);
}

}

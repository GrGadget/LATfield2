#include "config.h"
#include "particles/projections.hpp"

namespace LATfield2
{
void VecVecProjectionCIC_comm(Field<Real> * Tij)
{


    if(Tij->lattice().halo() == 0)
    {
        cout<< "LATfield2::scalarProjectionCIC_proj: the field has to have at least a halo of 1" <<endl;
        cout<< "LATfield2::scalarProjectionCIC_proj: aborting" <<endl;
        exit(-1);
    }

    Real *bufferSendUp;
    Real *bufferSendDown;

    Real *bufferRecUp;
    Real *bufferRecDown;


    long bufferSizeY;
    long bufferSizeZ;


    int sizeLocal[3];
    long sizeLocalGross[3];
    int sizeLocalOne[3];
    int halo = Tij->lattice().halo();

    for(int i=0;i<3;i++)
    {
        sizeLocal[i]=Tij->lattice().sizeLocal(i);
        sizeLocalGross[i] = sizeLocal[i] + 2 * halo;
        sizeLocalOne[i]=sizeLocal[i]+2;
    }

    int distHaloOne = halo - 1;








    int imax;

    int iref1=sizeLocalGross[0]-halo;
    int iref2=sizeLocalGross[0]-halo-1;

    for(int k=distHaloOne;k<sizeLocalOne[2]+distHaloOne;k++)
    {
        for(int j=distHaloOne;j<sizeLocalOne[1]+distHaloOne;j++)
        {
            for(int comp=0;comp<6;comp++)(*Tij)(setIndex(sizeLocalGross,halo,j,k),comp) += (*Tij)(setIndex(sizeLocalGross,iref1,j,k),comp);

            (*Tij)(setIndex(sizeLocalGross,iref2,j,k),0,1) += (*Tij)(setIndex(sizeLocalGross,distHaloOne,j,k),0,1);
            (*Tij)(setIndex(sizeLocalGross,iref2,j,k),0,2) += (*Tij)(setIndex(sizeLocalGross,distHaloOne,j,k),0,2);
            (*Tij)(setIndex(sizeLocalGross,iref2,j,k),1,2) += (*Tij)(setIndex(sizeLocalGross,distHaloOne,j,k),1,2);
        }
    }



    //build buffer size

    bufferSizeY = (long)(sizeLocalOne[2]) * (long)sizeLocal[0];
    bufferSizeZ = ((long)sizeLocal[0]) * ((long)sizeLocal[1]);



    if(bufferSizeY >= bufferSizeZ)
    {
        bufferSendUp = (Real*)malloc(6*bufferSizeY * sizeof(Real));
        bufferRecUp = (Real*)malloc(6*bufferSizeY * sizeof(Real));
        bufferSendDown = (Real*)malloc(3*bufferSizeY * sizeof(Real));
        bufferRecDown = (Real*)malloc(3*bufferSizeY * sizeof(Real));
    }
    else //if(bufferSizeZ> bufferSizeY )
    {
        bufferSendUp = (Real*)malloc(6*bufferSizeZ * sizeof(Real));
        bufferRecUp = (Real*)malloc(6*bufferSizeZ * sizeof(Real));
        bufferSendDown = (Real*)malloc(3*bufferSizeZ * sizeof(Real));
        bufferRecDown = (Real*)malloc(3*bufferSizeZ * sizeof(Real));
    }


    //working on Y dimension;


    ///verif:



    imax=sizeLocal[0];
    iref1 = sizeLocalGross[1]- halo;

    //COUT<<(*Tij)(setIndex(sizeLocalGross,6+halo,iref1,6+halo),0)<<endl;

    for(int k=0;k<sizeLocalOne[2];k++)
    {
        for(int i=0;i<imax;i++)
        {
            for(int comp =0;comp<6;comp++)bufferSendUp[comp + 6l * (i+k*imax)]=(*Tij)(setIndex(sizeLocalGross,i+halo,iref1,k+distHaloOne),comp);
        }

        for(int i=0;i<imax;i++)
        {
            bufferSendDown[    3l * (i+k*imax)]=(*Tij)(setIndex(sizeLocalGross,i+halo,distHaloOne,k+distHaloOne),0,1);
            bufferSendDown[1l + 3l * (i+k*imax)]=(*Tij)(setIndex(sizeLocalGross,i+halo,distHaloOne,k+distHaloOne),0,2);
            bufferSendDown[2l + 3l * (i+k*imax)]=(*Tij)(setIndex(sizeLocalGross,i+halo,distHaloOne,k+distHaloOne),1,2);
        }

    }

    //COUT<<bufferSendUp[0l + 6l * (6+7*imax)]<<endl;

    parallel.sendUpDown_dim1(bufferSendUp,bufferRecUp,bufferSizeY * 6l,bufferSendDown,bufferRecDown,bufferSizeY * 3l);

    //if(parallel.rank()==2)cout<<bufferRecUp[0l + 6l * (6+7*imax)]<<endl;

    //unpack data
    iref1 = sizeLocalGross[1]- halo - 1;


    for(int k=0;k<sizeLocalOne[2];k++)
    {
        for(int i=0;i<imax;i++)
        {
            for(int comp =0;comp<6;comp++)(*Tij)(setIndex(sizeLocalGross,i+halo,halo,k+distHaloOne),comp) += bufferRecUp[comp + 6l * (i+k*imax)];
        }

        for(int i=0;i<imax;i++)
        {
            (*Tij)(setIndex(sizeLocalGross,i+halo,iref1,k+distHaloOne),0,1) += bufferRecDown[     3l * (i+k*imax)];
            (*Tij)(setIndex(sizeLocalGross,i+halo,iref1,k+distHaloOne),0,2) += bufferRecDown[1l + 3l * (i+k*imax)];
            (*Tij)(setIndex(sizeLocalGross,i+halo,iref1,k+distHaloOne),1,2) += bufferRecDown[2l + 3l * (i+k*imax)];
        }

    }

    //if(parallel.rank()==2)cout<<(*Tij)(setIndex(sizeLocalGross,6+halo,halo ,6+halo),0)<<endl;

    //working on Z dimension;

    //pack data

    iref1=sizeLocalGross[2]- halo;

    for(int j=0;j<(sizeLocalOne[1]-2);j++)
    {
        for(int i=0;i<imax;i++)
        {
            for(int comp =0;comp<6;comp++)bufferSendUp[comp + 6l * (i+j*imax)]=(*Tij)(setIndex(sizeLocalGross,i+halo,j+halo,iref1),comp);
        }
        for(int i=0;i<imax;i++)
        {
            bufferSendDown[     3l * (i+j*imax)]=(*Tij)(setIndex(sizeLocalGross,i+halo,j+halo,distHaloOne),0,1);
            bufferSendDown[1l + 3l * (i+j*imax)]=(*Tij)(setIndex(sizeLocalGross,i+halo,j+halo,distHaloOne),0,2);
            bufferSendDown[2l + 3l * (i+j*imax)]=(*Tij)(setIndex(sizeLocalGross,i+halo,j+halo,distHaloOne),1,2);
        }

    }

    parallel.sendUpDown_dim0(bufferSendUp,bufferRecUp,bufferSizeZ * 6l,bufferSendDown,bufferRecDown,bufferSizeZ * 3l);

    //unpack data

    iref1 = sizeLocalGross[2]- halo - 1;

    for(int j=0;j<(sizeLocalOne[1]-2);j++)
    {
        for(int i=0;i<imax;i++)
        {
            for(int comp =0;comp<6;comp++)(*Tij)(setIndex(sizeLocalGross,i+halo,j+halo,halo),comp) += bufferRecUp[comp + 6l * (i+j*imax)];
        }

        for(int i=0;i<imax;i++)
        {

            (*Tij)(setIndex(sizeLocalGross,i+halo,j+halo,iref1),0,1) += bufferRecDown[     3l * (i+j*imax)];
            (*Tij)(setIndex(sizeLocalGross,i+halo,j+halo,iref1),0,2) += bufferRecDown[1l + 3l * (i+j*imax)];
            (*Tij)(setIndex(sizeLocalGross,i+halo,j+halo,iref1),1,2) += bufferRecDown[2l + 3l * (i+j*imax)];
        }


    }

    free(bufferSendUp);
    free(bufferSendDown);
    free(bufferRecUp);
    free(bufferRecDown);
}
void vectorProjectionCIC_comm(Field<Real> * vel)
{


    if(vel->lattice().halo() == 0)
    {
        cout<< "LATfield2::scalarProjectionCIC_proj: the field has to have at least a halo of 1" <<endl;
        cout<< "LATfield2::scalarProjectionCIC_proj: aborting" <<endl;
        exit(-1);
    }

    Real *bufferSendUp;
    Real *bufferSendDown;

    Real *bufferRecUp;
    Real *bufferRecDown;


    long bufferSizeY;
    long bufferSizeZ;


    int sizeLocal[3];
    long sizeLocalGross[3];
    int sizeLocalOne[3];
    int halo = vel->lattice().halo();

    for(int i=0;i<3;i++)
    {
        sizeLocal[i]=vel->lattice().sizeLocal(i);
        sizeLocalGross[i] = sizeLocal[i] + 2 * halo;
        sizeLocalOne[i]=sizeLocal[i]+2;
    }

    int distHaloOne = halo - 1;

    int imax;
    //update from halo data in X;

    int iref1=sizeLocalGross[0]-halo;
    int iref2=sizeLocalGross[0]-halo-1;

    for(int k=distHaloOne;k<sizeLocalOne[2]+distHaloOne;k++)
    {
        for(int j=distHaloOne;j<sizeLocalOne[1]+distHaloOne;j++)
        {
            for(int comp=0;comp<3;comp++)(*vel)(setIndex(sizeLocalGross,halo,j,k),comp) += (*vel)(setIndex(sizeLocalGross,iref1,j,k),comp);

            for(int comp=0;comp<3;comp++)(*vel)(setIndex(sizeLocalGross,iref2,j,k),comp) += (*vel)(setIndex(sizeLocalGross,distHaloOne,j,k),comp);
        }
    }

    //communication



    //build buffer size

    bufferSizeY = (long)(sizeLocalOne[2]) * sizeLocal[0];
    bufferSizeZ = sizeLocal[0] * sizeLocal[1];

    if(bufferSizeY>bufferSizeZ)
    {
        bufferSendUp = (Real*)malloc(3*bufferSizeY * sizeof(Real));
        bufferRecUp = (Real*)malloc(3*bufferSizeY * sizeof(Real));
        bufferSendDown = (Real*)malloc(3*bufferSizeY * sizeof(Real));
        bufferRecDown = (Real*)malloc(3*bufferSizeY * sizeof(Real));
    }
    else
    {
        bufferSendUp = (Real*)malloc(3*bufferSizeZ * sizeof(Real));
        bufferRecUp = (Real*)malloc(3*bufferSizeZ * sizeof(Real));
        bufferSendDown = (Real*)malloc(3*bufferSizeZ * sizeof(Real));
        bufferRecDown = (Real*)malloc(3*bufferSizeZ * sizeof(Real));
    }


    //working on Y dimension


    imax=sizeLocal[0];
    iref1 = sizeLocalGross[1]- halo;


    for(int k=0;k<sizeLocalOne[2];k++)
    {
        for(int i=0;i<imax;i++)
        {
            for(int comp =0;comp<3;comp++)bufferSendUp[comp + 3l * (i+k*imax)]=(*vel)(setIndex(sizeLocalGross,i+halo,iref1,k+distHaloOne),comp);
        }
        for(int i=0;i<imax;i++)
        {
            for(int comp =0;comp<3;comp++)bufferSendDown[comp + 3l * (i+k*imax)]=(*vel)(setIndex(sizeLocalGross,i+halo,distHaloOne,k+distHaloOne),comp);
        }
    }

    parallel.sendUpDown_dim1(bufferSendUp,bufferRecUp,bufferSizeY * 3l,bufferSendDown,bufferRecDown,bufferSizeY * 3l);

    //unpack data
    iref1 = sizeLocalGross[1]- halo - 1;

    for(int k=0;k<sizeLocalOne[2];k++)
    {
        for(int i=0;i<imax;i++)
        {
            for(int comp =0;comp<3;comp++)(*vel)(setIndex(sizeLocalGross,i+halo,halo,k+distHaloOne),comp) += bufferRecUp[comp + 3l * (i+k*imax)];
        }
        for(int i=0;i<imax;i++)
        {
            for(int comp =0;comp<3;comp++)(*vel)(setIndex(sizeLocalGross,i+halo,iref1,k+distHaloOne),comp) += bufferRecDown[comp + 3l * (i+k*imax)];

        }
    }

    //workin on dim z

    //pack data

    iref1=sizeLocalGross[2]- halo;

    for(int j=0;j<(sizeLocalOne[1]-2);j++)
    {
        for(int i=0;i<imax;i++)
        {
            for(int comp =0;comp<3;comp++)bufferSendUp[comp + 3l * (i+j*imax)]=(*vel)(setIndex(sizeLocalGross,i+halo,j+halo,iref1),comp);
        }
        for(int i=0;i<imax;i++)
        {
            for(int comp =0;comp<3;comp++)bufferSendDown[comp + 3l * (i+j*imax)]=(*vel)(setIndex(sizeLocalGross,i+halo,j+halo,distHaloOne),comp);
        }

    }

    parallel.sendUpDown_dim0(bufferSendUp,bufferRecUp,bufferSizeY * 3l,bufferSendDown,bufferRecDown,bufferSizeY * 3l);

    //unpack data

    iref1 = sizeLocalGross[2]- halo - 1;

    for(int j=0;j<(sizeLocalOne[1]-2);j++)
    {
        for(int i=0;i<imax;i++)
        {
            for(int comp =0;comp<3;comp++)(*vel)(setIndex(sizeLocalGross,i+halo,j+halo,halo),comp) += bufferRecUp[comp + 3l * (i+j*imax)];
        }
        for(int i=0;i<imax;i++)
        {
            for(int comp =0;comp<3;comp++)(*vel)(setIndex(sizeLocalGross,i+halo,j+halo,iref1),comp) += bufferRecDown[comp + 3l * (i+j*imax)];
        }

    }


    free(bufferSendUp);
    free(bufferSendDown);
    free(bufferRecUp);
    free(bufferRecDown);


}
void symtensorProjectionCICNGP_comm(Field<Real> * Tij)
{

    if(Tij->lattice().halo() == 0)
    {
        cout<< "LATfield2::scalarProjectionCIC_proj: the field has to have at least a halo of 1" <<endl;
        cout<< "LATfield2::scalarProjectionCIC_proj: aborting" <<endl;
        exit(-1);
    }

    Real *bufferSend;
    Real *bufferRec;


    long bufferSizeY;
    long bufferSizeZ;


    int sizeLocal[3];
    long sizeLocalGross[3];
    int sizeLocalOne[3];
    int halo = Tij->lattice().halo();

    for(int i=0;i<3;i++)
    {
        sizeLocal[i]=Tij->lattice().sizeLocal(i);
        sizeLocalGross[i] = sizeLocal[i] + 2 * halo;
        sizeLocalOne[i]=sizeLocal[i]+2;
    }

    int distHaloOne = halo - 1;

    int iref;
    int imax;




    int comp=6;

    iref = sizeLocalGross[0]-halo;
    for(int k=distHaloOne;k<sizeLocalOne[2]+distHaloOne;k++)
    {
        for(int j=distHaloOne;j<sizeLocalOne[1]+distHaloOne;j++)
        {
            for(int c=0;c<comp;c++)(*Tij)(setIndex(sizeLocalGross,halo,j,k),c) += (*Tij)(setIndex(sizeLocalGross,iref,j,k),c);
        }
    }


    //send halo in direction Y
    bufferSizeY =  (long)(sizeLocalOne[2]-1)*sizeLocal[0] * comp;
    bufferSizeZ = sizeLocal[0] * sizeLocal[1] * comp;
    if(bufferSizeY>bufferSizeZ)
    {
        bufferSend = (Real*)malloc(sizeof(Real)*bufferSizeY);
        bufferRec = (Real*)malloc(sizeof(Real)*bufferSizeY);
    }
    else
    {
        bufferSend = (Real*)malloc(sizeof(Real)*bufferSizeZ);
        bufferRec = (Real*)malloc(sizeof(Real)*bufferSizeZ);
    }


    //pack data
    imax=sizeLocalGross[0]-2* halo;
    iref=sizeLocalGross[1]- halo;
    for(int k=0;k<(sizeLocalOne[2]-1);k++)
    {
        for(int i=0;i<imax;i++)
        {
            for(int c=0;c<comp;c++)bufferSend[c+comp*(i+k*imax)]=(*Tij)(setIndex(sizeLocalGross,i+halo,iref,k+halo),c);
        }
    }


    parallel.sendUp_dim1(bufferSend,bufferRec,bufferSizeY);


    //unpack data
    for(int k=0;k<(sizeLocalOne[2]-1);k++)
    {
        for(int i=0;i<imax;i++)
        {
            for(int c=0;c<comp;c++)(*Tij)(setIndex(sizeLocalGross,i+halo,halo,k+halo),c)+=bufferRec[c+comp*(i+k*imax)];
        }

    }

    //send halo in direction Z

    //pack data
    iref=sizeLocalGross[2]-halo;
    for(int j=0;j<(sizeLocalOne[1]-2);j++)
    {
        for(int i=0;i<imax;i++)
        {
            for(int c=0;c<comp;c++)bufferSend[c+comp*(i+j*imax)]=(*Tij)(setIndex(sizeLocalGross,i+halo,j+halo,iref),c);
        }
    }

    parallel.sendUp_dim0(bufferSend,bufferRec,bufferSizeZ);


    //unpack data

    for(int j=0;j<(sizeLocalOne[1]-2);j++)
    {
        for(int i=0;i<imax;i++)
        {
            for(int c=0;c<comp;c++)(*Tij)(setIndex(sizeLocalGross,i+halo,j+halo,halo),c)+=bufferRec[c+comp*(i+j*imax)];
        }
    }

    free(bufferRec);
    free(bufferSend);



}
void vectorProjectionCICNGP_comm(Field<Real> * vel)
{

    if(vel->lattice().halo() == 0)
    {
        cout<< "LATfield2::scalarProjectionCIC_proj: the field has to have at least a halo of 1" <<endl;
        cout<< "LATfield2::scalarProjectionCIC_proj: aborting" <<endl;
        exit(-1);
    }

    Real *bufferSend;
    Real *bufferRec;


    long bufferSizeY;
    long bufferSizeZ;


    int sizeLocal[3];
    long sizeLocalGross[3];
    int sizeLocalOne[3];
    int halo = vel->lattice().halo();

    for(int i=0;i<3;i++)
    {
        sizeLocal[i]=vel->lattice().sizeLocal(i);
        sizeLocalGross[i] = sizeLocal[i] + 2 * halo;
        sizeLocalOne[i]=sizeLocal[i]+2;
    }

    int distHaloOne = halo - 1;

    int iref;
    int imax;




    int comp=3;


    iref =sizeLocalGross[0] - halo ;
    for(int k=distHaloOne;k<sizeLocalOne[2]+distHaloOne;k++)
    {
        for(int j=distHaloOne;j<sizeLocalOne[1]+distHaloOne;j++)
        {
            for(int c=0;c<comp;c++)(*vel)(setIndex(sizeLocalGross,halo,j,k),c) += (*vel)(setIndex(sizeLocalGross,iref,j,k),c);
        }
    }


    //send halo in direction Y
    bufferSizeY =  (long)(sizeLocalOne[2]-1)*sizeLocal[0] * comp;
    bufferSizeZ = sizeLocal[0] * sizeLocal[1] * comp;
    if(bufferSizeY>bufferSizeZ)
    {
        bufferSend = (Real*)malloc(sizeof(Real)*bufferSizeY);
        bufferRec = (Real*)malloc(sizeof(Real)*bufferSizeY);
    }
    else
    {
        bufferSend = (Real*)malloc(sizeof(Real)*bufferSizeZ);
        bufferRec = (Real*)malloc(sizeof(Real)*bufferSizeZ);
    }


    //pack data
    imax=sizeLocalGross[0]-2* halo;
    iref=sizeLocalGross[1]- halo;
    for(int k=0;k<(sizeLocalOne[2]-1);k++)
    {
        for(int i=0;i<imax;i++)
        {
            for(int c=0;c<comp;c++)bufferSend[c+comp*(i+k*imax)]=(*vel)(setIndex(sizeLocalGross,i+halo,iref,k+halo),c);
        }
    }


    parallel.sendUp_dim1(bufferSend,bufferRec,bufferSizeY);


    //unpack data
    for(int k=0;k<(sizeLocalOne[2]-1);k++)
    {
        for(int i=0;i<imax;i++)
        {
            for(int c=0;c<comp;c++)(*vel)(setIndex(sizeLocalGross,i+halo,halo,k+halo),c)+=bufferRec[c+comp*(i+k*imax)];
        }

    }

    //send halo in direction Z

    //pack data

    //cout<<"okok"<<endl;

    iref=sizeLocalGross[2]-halo;
    for(int j=0;j<(sizeLocalOne[1]-2);j++)
    {
        for(int i=0;i<imax;i++)
        {
            for(int c=0;c<comp;c++)bufferSend[c+comp*(i+j*imax)]=(*vel)(setIndex(sizeLocalGross,i+halo,j+halo,iref),c);
        }
    }

    parallel.sendUp_dim0(bufferSend,bufferRec,bufferSizeZ);


    //unpack data

    for(int j=0;j<(sizeLocalOne[1]-2);j++)
    {
        for(int i=0;i<imax;i++)
        {
            for(int c=0;c<comp;c++)(*vel)(setIndex(sizeLocalGross,i+halo,j+halo,halo),c)+=bufferRec[c+comp*(i+j*imax)];
        }
    }

    free(bufferRec);
    free(bufferSend);


}
void vertexProjectionCIC_comm(Field<Real> * vel)
{

      if(vel->lattice().halo() == 0)
      {
          cout<< "NBCC_CIC_vel_comm: the field has to have at least a halo of 1" <<endl;
          cout<< "NBCC_CIC_vel_comm: aborting" <<endl;
          exit(-1);
      }

      Real *bufferSend;
      Real *bufferRec;


      long bufferSizeY;
      long bufferSizeZ;


      int sizeLocal[3];
      long sizeLocalGross[3];
      int sizeLocalOne[3];
      int halo = vel->lattice().halo();

      for(int i=0;i<3;i++)
      {
          sizeLocal[i]=vel->lattice().sizeLocal(i);
          sizeLocalGross[i] = sizeLocal[i] + 2 * halo;
          sizeLocalOne[i]=sizeLocal[i]+2;
      }

      int distHaloOne = halo - 1;

      int iref;
      int imax;

      bufferSizeY =  (long)(sizeLocalOne[2]-1) * (long)sizeLocal[0];
      bufferSizeZ = (long)sizeLocal[0] * (long)sizeLocal[1] ;
      if(bufferSizeY>bufferSizeZ)
      {
          bufferSend = (Real*)malloc(sizeof(Real)*bufferSizeY*3);
          bufferRec = (Real*)malloc(sizeof(Real)*bufferSizeY*3);
      }
      else
      {
          bufferSend = (Real*)malloc(sizeof(Real)*bufferSizeZ*3);
          bufferRec = (Real*)malloc(sizeof(Real)*bufferSizeZ*3);
      }

      //send halo in direction X
      iref = sizeLocalGross[0] - halo;
      for(int k=distHaloOne;k<sizeLocalOne[2]+distHaloOne;k++)
      {
          for(int j=distHaloOne;j<sizeLocalOne[1]+distHaloOne;j++)
          {
              for(int c=0;c<3;c++)(*vel)(setIndex(sizeLocalGross,halo,j,k),c) += (*vel)(setIndex(sizeLocalGross,iref,j,k),c);
          }
      }
      //send halo in direction Y
      imax = sizeLocal[0];
      iref = sizeLocalGross[1]- halo;
      for(int k=0;k<(sizeLocalOne[2]-1);k++)
      {
          for(int i=0;i<imax;i++)
          {
              for(int c=0;c<3;c++)bufferSend[c+3*(i+k*imax)]=(*vel)(setIndex(sizeLocalGross,i+halo,iref,k+halo),c);
          }
      }
      parallel.sendUp_dim1(bufferSend,bufferRec,bufferSizeY*3);
      for(int k=0;k<(sizeLocalOne[2]-1);k++)
      {
          for(int i=0;i<imax;i++)
          {
            for(int c=0;c<3;c++)(*vel)(setIndex(sizeLocalGross,i+halo,halo,k+halo),c)+=bufferRec[c+3*(i+k*imax)];
          }

      }

      //send halo in direction Z
      iref=sizeLocalGross[2]-halo;
      for(int j=0;j<(sizeLocalOne[1]-2);j++)
      {
          for(int i=0;i<imax;i++)
          {
              for(int c=0;c<3;c++)bufferSend[c+3*(i+j*imax)]=(*vel)(setIndex(sizeLocalGross,i+halo,j+halo,iref),c);
          }
      }
      parallel.sendUp_dim0(bufferSend,bufferRec,bufferSizeZ*3);
      for(int j=0;j<(sizeLocalOne[1]-2);j++)
      {
          for(int i=0;i<imax;i++)
          {
            for(int c=0;c<3;c++)(*vel)(setIndex(sizeLocalGross,i+halo,j+halo,halo),c)+=bufferRec[c+3*(i+j*imax)];
          }
      }

      free(bufferRec);
      free(bufferSend);
}
void scalarProjectionCIC_comm(Field<Real> * rho)
{

    if(rho->lattice().halo() == 0)
    {
        cout<< "LATfield2::scalarProjectionCIC_proj: the field has to have at least a halo of 1" <<endl;
        cout<< "LATfield2::scalarProjectionCIC_proj: aborting" <<endl;
        exit(-1);
    }

    Real *bufferSend;
    Real *bufferRec;


    long bufferSizeY;
    long bufferSizeZ;


    int sizeLocal[3];
    long sizeLocalGross[3];
    int sizeLocalOne[3];
    int halo = rho->lattice().halo();

    for(int i=0;i<3;i++)
    {
        sizeLocal[i]=rho->lattice().sizeLocal(i);
        sizeLocalGross[i] = sizeLocal[i] + 2 * halo;
        sizeLocalOne[i]=sizeLocal[i]+2;
    }

    int distHaloOne = halo - 1;

    int iref;
    int imax;


    iref = sizeLocalGross[0] - halo;

    for(int k=distHaloOne;k<sizeLocalOne[2]+distHaloOne;k++)
    {
        for(int j=distHaloOne;j<sizeLocalOne[1]+distHaloOne;j++)
        {
            (*rho)(setIndex(sizeLocalGross,halo,j,k)) += (*rho)(setIndex(sizeLocalGross,iref,j,k));
        }
    }


    //send halo in direction Y
    bufferSizeY =  (long)(sizeLocalOne[2]-1) * (long)sizeLocal[0];
    bufferSizeZ = (long)sizeLocal[0] * (long)sizeLocal[1] ;
    if(bufferSizeY>bufferSizeZ)
    {
        bufferSend = (Real*)malloc(sizeof(Real)*bufferSizeY);
        bufferRec = (Real*)malloc(sizeof(Real)*bufferSizeY);
    }
    else
    {
        bufferSend = (Real*)malloc(sizeof(Real)*bufferSizeZ);
        bufferRec = (Real*)malloc(sizeof(Real)*bufferSizeZ);
    }

    //pack data
    imax = sizeLocal[0];
    iref = sizeLocalGross[1]- halo;
    for(int k=0;k<(sizeLocalOne[2]-1);k++)
    {
        for(int i=0;i<imax;i++)
        {
            bufferSend[i+k*imax]=(*rho)(setIndex(sizeLocalGross,i+halo,iref,k+halo));
        }
    }
    parallel.sendUp_dim1(bufferSend,bufferRec,bufferSizeY);
    //unpack data
    for(int k=0;k<(sizeLocalOne[2]-1);k++)
    {
        for(int i=0;i<imax;i++)(*rho)(setIndex(sizeLocalGross,i+halo,halo,k+halo))+=bufferRec[i+k*imax];

    }

    //send halo in direction Z



    //pack data
    iref=sizeLocalGross[2]-halo;
    for(int j=0;j<(sizeLocalOne[1]-2);j++)
    {
        for(int i=0;i<imax;i++)
        {
            bufferSend[i+j*imax]=(*rho)(setIndex(sizeLocalGross,i+halo,j+halo,iref));
        }
    }

    parallel.sendUp_dim0(bufferSend,bufferRec,bufferSizeZ);


    //unpack data

    for(int j=0;j<(sizeLocalOne[1]-2);j++)
    {
        for(int i=0;i<imax;i++)(*rho)(setIndex(sizeLocalGross,i+halo,j+halo,halo))+=bufferRec[i+j*imax];
    }

    free(bufferRec);
    free(bufferSend);

}
}

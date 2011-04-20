#include "mpi_utils.h"

int utils::init_mpi(int *argc, char ***argv, int *mpi_rank, int *mpi_size)
{
  int retval = 0;

  retval = MPI_Init(argc,argv);
  if(retval)
    return retval;
  
  retval = MPI_Comm_rank(MPI_COMM_WORLD,mpi_rank);
  if(retval)
    return retval;
  retval = MPI_Comm_size(MPI_COMM_WORLD,mpi_size);

  return retval;
}


void utils::mpi_recombine(GLAlgorithm& gls, MPI_Comm plane_creators)
{

  double *shared_array = NULL;
  size_t rows,columns;
  int retval = 0;
  Plane<math::Complex> *lens = NULL;

  lens = gls.getDeflectionPlane();
  if(mpi_data.rank == utils::MASTER)
    {
      if(lens == NULL)
	throw DavidException("utils::mpi_recombine: Master has a NULL lens plane. Cannot recombine.",DavidException::DATA_FORMAT_ERROR);
      rows = lens->numberOfRows();
      columns = lens->numberOfColumns();
      MPI_Bcast(&rows,1,MPI_DOUBLE,mpi_data.rank,plane_creators);
      MPI_Bcast(&columns,1,MPI_DOUBLE,utils::MASTER,plane_creators);
    }
  else
    {
      if(lens == NULL)
	return;

      MPI_Bcast(&rows,1,MPI_DOUBLE,utils::MASTER,plane_creators);
      MPI_Bcast(&columns,1,MPI_DOUBLE,utils::MASTER,plane_creators);
      if(rows != lens->numberOfRows() || columns != lens->numberOfColumns())
	{
	  std::ostringstream error;
	  error << "utils::mpi_recombine: Process " << mpi_data.rank 
		<< " has an array of differing dimensions as the Master." 
		<< std::endl<<
		<< "Master Plane is " << rows << " by " << columns <<std::endl
		<< "Process " <<  mpi_data.rank << " is " << 
	    lens->numberOfRows() << " by " << != lens->numberOfColumns() ;
	  throw DavidException(error,DavidException::DATA_FORMAT_ERROR);
	}
    }
  shared_array = (double*)malloc(2*rows*columns*sizeof(double));
  
for(size_t i = 0;i<rows;i++)
  {
    for(size_t j = 0;j<columns;j++)
      {
	for(size_k = 0;k<2;k++)
	  shared_array[2*j+i*rows+k] = lens->getValue(i,j).getValue(k);
      }
  }
 retval = MPI_Allreduce(MPI_IN_PLACE,shared_array,2*rows*columns,MPI_DOUBLE,MPI_SUM,plane_creators);
    
 if(retval)
   {
     std::ostringstream error;
     char mpierror[MPI_MAX_ERROR_STRING];
     size dummy;
     error << "utils::mpi_recombine: Process " << mpi_data.rank
	   << "Could not MPI_Allreduce data. " << std::endl
	   << "Reason: " << MPI_Error_string(retval,mpierror,&size);
     throw DavidException(error,DavidException::MPI_ERROR_CODE);
   }

 for(size_t i = 0;i<rows;i++)
  {
    for(size_t j = 0;j<columns;j++)
      {
	math::Complex cnum(shared_array[2*j+i*rows],shared_array[2*j+i*rows+1]);
	lens->setValue(i,j, cnum);
      }
  }

 free(shared_array);

 MPI_Barrier(plane_creators);

}



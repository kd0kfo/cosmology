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

void utils::mpi_recombine(Plane<math::Complex> *plane, MPI_Comm plane_creators)
{
  double *shared_array = NULL;
  size_t rows,columns;
  int retval = 0;

  if(mpi_data.rank == utils::MASTER_RANK)
    {
      if(plane == NULL)
	throw DavidException("utils::mpi_recombine: Master has a NULL plane. Cannot recombine.",DavidException::DATA_FORMAT_ERROR);
      rows = plane->numberOfRows();
      columns = plane->numberOfColumns();
      MPI_Bcast(&rows,1,MPI_DOUBLE,mpi_data.rank,plane_creators);
      MPI_Bcast(&columns,1,MPI_DOUBLE,utils::MASTER_RANK,plane_creators);
    }
  else
    {
      if(plane == NULL)
	return;

      MPI_Bcast(&rows,1,MPI_DOUBLE,utils::MASTER_RANK,plane_creators);
      MPI_Bcast(&columns,1,MPI_DOUBLE,utils::MASTER_RANK,plane_creators);
      if(rows != plane->numberOfRows() || columns != plane->numberOfColumns())
	{
	  std::ostringstream error;
	  error << "utils::mpi_recombine: Process " << mpi_data.rank 
		<< " has an array of differing dimensions as the Master." 
		<< std::endl
		<< "Master Plane is " << rows << " by " << columns <<std::endl
		<< "Process " <<  mpi_data.rank << " is " << 
	    plane->numberOfRows() << " by " <<  plane->numberOfColumns() ;
	  throw DavidException(error,DavidException::DATA_FORMAT_ERROR);
	}
    }
  shared_array = (double*)malloc(2*rows*columns*sizeof(double));
  
for(size_t i = 0;i<rows;i++)
  {
    for(size_t j = 0;j<columns;j++)
      {
    	shared_array[2*j+i*rows] = plane->getValue(i,j).real();
    	shared_array[2*j+i*rows+1] = plane->getValue(i,j).imag();
      }
  }
 retval = MPI_Allreduce(MPI_IN_PLACE,shared_array,2*rows*columns,MPI_DOUBLE,MPI_SUM,plane_creators);
    
 if(retval)
   {
     std::ostringstream error;
     char mpierror[MPI_MAX_ERROR_STRING];
     int size;
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
    	plane->setValue(i,j, cnum);
      }
  }

 free(shared_array);

 MPI_Barrier(plane_creators);

}

void utils::mpi_recombine(GLAlgorithm& gls, MPI_Comm plane_creators)
{

  Plane<math::Complex> *lens  = gls.getDeflectionPlane();

  mpi_recombine(lens,plane_creators);

}

void utils::mpi_adjust_glellipsebounds(int *glellipseBounds, size_t N)
{
	int *bound_vals = NULL;
	int old_bounds[4];
	size_t X,Y,increment;
	if(glellipseBounds == NULL && N == 0)
		return;

	if(glellipseBounds == NULL)
	  {
	    old_bounds[1] = old_bounds[3] = N;
	    old_bounds[0] = old_bounds[2] = 0;
	  }
	else
	  {
	    if(glellipseBounds[0] == -1)
	      glellipseBounds[0] = 0;
	    if(glellipseBounds[1] == -1)
	      glellipseBounds[1] = N;
	    if(glellipseBounds[2] == -1)
	      glellipseBounds[2] = 0;
	    if(glellipseBounds[3] == -1)
	      glellipseBounds[3] = N;
	    memcpy(old_bounds,glellipseBounds,4*sizeof(int));
	    
	  }
	X = old_bounds[1]- old_bounds[0];
	Y = old_bounds[3] - old_bounds[2];
	N = X*Y;

	if(mpi_data.num_ranks < Y)
	  {
	    glellipseBounds[2] = old_bounds[2] + mpi_data.rank * (int)ceil(Y/mpi_data.num_ranks);
	    glellipseBounds[3] = old_bounds[2] + (mpi_data.rank+1) * (int)ceil(Y/mpi_data.num_ranks);
	  }
	else
	  {
	    float subdivision = (float)X/Y;
	    glellipseBounds[3] = old_bounds[2] + mpi_data.rank * (int)ceil(Y/mpi_data.num_ranks);
	    if(glellipseBounds[3] > old_bounds[3])
	      glellipseBounds[3] = old_bounds[3];
	    
	    if(mpi_data.rank % 2 == 0)
	      glellipseBounds[1] = old_bounds[0] + (int)ceil(subdivision);
	    else
	      {
		glellipseBounds[0] = old_bounds[0] + (int)ceil(subdivision);
	      }

	    if(mpi_data.num_ranks % 2 != 0)
	      {
		// if there are an odd number of slaves, have the last slave
		// do a whole row.
		if(mpi_data.rank = mpi_data.num_ranks - 1)
		    memcpy(&glellipseBounds[2],&old_bounds[2],2*sizeof(int));
	      }
	  }
}



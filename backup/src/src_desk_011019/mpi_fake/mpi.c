int MPI_Init(int *argc, char ***argv){

  *argc = *argc;
  ***argv = ***argv;
  return 0;
}


int MPI_Comm_rank(int comm, int *rank){

  comm = 0;
  *rank = 0;
  return 0;
}

int MPI_Comm_size(int comm, int *size){

  comm = 0;
  *size = 1;
  return 0;
}


int MPI_Barrier(int comm){

  comm = 0;

  return 0;
}


int MPI_Finalize(void){
  
  return 0;
}

int MPI_Send(const void *buf, int count, int datatype, int dest,
	     int tag, int comm){

  buf= 0;
  count = 0;
  datatype = 0;
  dest = 0;
  tag = 0;
  comm = 0;  

  return 0;
}

int MPI_Recv(void *buf, int count, int datatype, int source,
	     int tag, int comm, int status){

  buf= 0;
  count = 0;
  datatype = 0;
  source = 0;
  tag = 0;
  comm =  0;
  status = 0 ;

  return 0;
}


// NOTE: the real mpi.h is in /opt/local/include/openmpi-gcc49/ on Charles's CLASP laptop

int MPI_Init(int *argc, char ***argv);
int MPI_Comm_rank(int comm, int *rank);
int MPI_Comm_size(int comm, int *size);
int MPI_Barrier(int comm);
int MPI_Finalize(void);
int MPI_Send(const void *buf, int count, int datatype, int dest,
	     int tag, int comm);
int MPI_Recv(void *buf, int count, int datatype, int source,
	     int tag, int comm, int status);

#define MPI_COMM_WORLD 1
#define MPI_DOUBLE 1
#define MPI_STATUS_IGNORE 1
#define MPI_INT 1 

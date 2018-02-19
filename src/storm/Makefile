ifdef COMPILER
	CC=${PATH_COMPILER} -DINLINE -Wall -g -std=gnu99 -O0 -Wextra
	DEPS = ../options.h ../propagator.h ..nrlmsise-00.h prop_math.h fly_storm.h ../mpi_fake/mpi.h

	OBJ = ./fly_storm.o ../nrlmsise-00.o ../nrlmsise-00_data.o \
	../prop_math.o \
        ../mpi_fake/mpi.o \
	../propagator.o ../load_options.o \
	../generate_ephemerides.o ../initialize_constellation.o 
	MPI_H_DUMMY := $(shell ln -s ./mpi_fake/mpi.h ../mpi.h)

else
	CC=${PATH_COMPILER} -DINLINE -Wall -g -std=gnu99 -O0 -Wextra -w
	DEPS = ../options.h ../propagator.h ..nrlmsise-00.h prop_math.h fly_storm.h

	OBJ = ./fly_storm.o ../nrlmsise-00.o ../nrlmsise-00_data.o \
	../prop_math.o \
	../propagator.o ../load_options.o \
	../generate_ephemerides.o ../initialize_constellation.o \

	MPI_H_REMOVE_DUMMY := $(shell rm -f ./mpi.h)
	MPI_H_REMOVE_DUMMY_AGAIN := $(shell rm -f ../mpi.h)
endif

SPICE_DIR = ${PATH_SPICE}
GSL_DIR   = ${PATH_GSL}


CFLAGS=-g -I. -I./../ \
-I$(GSL_DIR) -I$(GSL_DIR)/include -I$(SPICE_DIR) -I$(SPICE_DIR)/include

LIBS=-lgsl -lgslcblas $(SPICE_DIR)/lib/csupport.a -lm \
$(SPICE_DIR)/lib/cspice.a 

LDFLAGS=-L$(GSL_DIR)/lib



%.o: %.c $(DEPS)
	$(CC)  -c -o $@ $<  $(CFLAGS) $(LIBS) $(LDFLAGS)

all: storm

storm: $(OBJ)
	$(CC) -o ${PATH_EXECUTABLE}/$@ $^ $(CFLAGS) $(LIBS) $(LDFLAGS)

clean: 
	rm -f ${PATH_EXECUTABLE}/storm
	rm -f ../*.o
	rm -f ./*o

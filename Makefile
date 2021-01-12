ifdef COMPILER
	CC=${PATH_COMPILER} -DINLINE -Wall -g -std=gnu99 -O0 -Wextra -fdiagnostics-color=auto
	DEPS = options.h propagator.h nrlmsise-00.h prop_math.h kalman_9state_new_with_tau.h ./src/mpi_fake/mpi.h
	OBJ = ./src/main.o ./src/nrlmsise-00.o ./src/nrlmsise-00_data.o \
	./src/prop_math.o \
	./src/kalman_9state_new_with_tau.o \
	./src/mpi_fake/mpi.o \
	./src/propagator_earth_pressure_faster.o ./src/load_options.o \
	./src/generate_ephemerides.o ./src/initialize_constellation.o 
	MPI_H_DUMMY := $(shell ln -s mpi_fake/mpi.h ./src/mpi.h)

else
	CC=${PATH_COMPILER} -DINLINE -g -std=gnu99 -O0 -fdiagnostics-color=auto 
	DEPS = options.h propagator.h nrlmsise-00.h prop_math.h kalman_9state_new_with_tau.h
	OBJ = ./src/main.o ./src/nrlmsise-00.o ./src/nrlmsise-00_data.o \
	./src/prop_math.o \
	./src/kalman_9state_new_with_tau.o \
	./src/propagator_earth_pressure_faster.o ./src/load_options.o \
	./src/generate_ephemerides.o ./src/initialize_constellation.o 
	MPI_H_REMOVE_DUMMY := $(shell rm -f ./src/mpi.h)
endif


SPICE_DIR = ${PATH_SPICE}
GSL_DIR   = ${PATH_GSL}


CFLAGS=-g -I. -I./src \
-I$(GSL_DIR) -I$(GSL_DIR)/include -I$(SPICE_DIR) -I$(SPICE_DIR)/include

LIBS=-lgsl -lgslcblas $(SPICE_DIR)/lib/csupport.a -lm \
$(SPICE_DIR)/lib/cspice.a 

LDFLAGS=-L$(GSL_DIR)/lib

%.o: %.c $(DEPS)
	$(CC)  -c -o $@ $<  $(CFLAGS) $(LIBS) $(LDFLAGS)

all: spock 

spock: $(OBJ)
	$(CC) -o ${PATH_EXECUTABLE}/$@ $^ $(CFLAGS) $(LIBS) $(LDFLAGS) 

clean:
	find  ${PATH_EXECUTABLE} -type f -maxdepth 1 -name spock -delete #rm -f  ${PATH_EXECUTABLE}/spock_dev_again # ignore the warning with maxdepth (it does work)
	rm -f ./src/*.o
	rm -f ./src/mpi_fake/*.o


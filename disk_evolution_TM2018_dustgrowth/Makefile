CPP = g++

CFLAGS = -g -c -O3 -march=native -std=c++17

# LFLAGS = -lgsl -lgslcblas -lm

LFLAGS =

DEPENDPATH = -I.

objects = src/main.o \
		  src/utils.o \
	      src/explicit_solver/explicit_solver.o \
		  src/sim_data.o \
		  src/bonner_ebert_sphere.o \
		  src/opacity.o \
		  src/reconstruction.o \
		  src/disk_evolution.o 


run: $(objects)
	$(CPP) -o run $(objects) $(LFLAGS)

src/main.o: src/main.cpp
	(cd src; $(CPP) $(CFLAGS) $(DEPENDPATH) main.cpp)

src/utils.o: src/utils.hpp src/utils.cpp 
	(cd src; $(CPP) $(CFLAGS) $(DEPENDPATH) utils.cpp)

src/explicit_solver/explicit_solver.o: src/explicit_solver/explicit_solver.hpp src/explicit_solver/explicit_solver.cpp 
	(cd src/explicit_solver; $(CPP) $(CFLAGS) $(DEPENDPATH) explicit_solver.cpp)

src/sim_data.o: src/sim_data.hpp src/sim_data.cpp
	(cd src; $(CPP) $(CFLAGS) $(DEPENDPATH) sim_data.cpp)

src/bonner_ebert_sphere.o: src/bonner_ebert_sphere.hpp src/bonner_ebert_sphere.cpp 
	(cd src; $(CPP) $(CFLAGS) $(DEPENDPATH) bonner_ebert_sphere.cpp)

src/opacity.o: src/opacity.hpp src/opacity.cpp
	(cd src; $(CPP) $(CFLAGS) $(DEPENDPATH) opacity.cpp)

src/reconstruction.o: src/reconstruction.hpp src/reconstruction.cpp 
	(cd src; $(CPP) $(CFLAGS) $(DEPENDPATH) reconstruction.cpp)

src/disk_evolution.o: src/disk_evolution.hpp src/disk_evolution.cpp 
	(cd src; $(CPP) $(CFLAGS) $(DEPENDPATH) disk_evolution.cpp)


clean:
	rm run $(objects)
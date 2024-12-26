CPP = g++

CFLAGS = -g -c -O3 -std=c++17

LFLAGS =

DEPENDPATH = -I.

objects = src/main.o \
		  src/utils.o \
		  src/sim_data.o \
		  src/bonner_ebert_sphere.o \
		  src/grid.o \
		  src/disk_evolution_driver.o 


run: $(objects)
	$(CPP) -o run $(objects) $(LFLAGS)

src/main.o: src/main.cpp
	(cd src; $(CPP) $(CFLAGS) $(DEPENDPATH) main.cpp)

src/utils.o: src/utils.hpp src/utils.cpp 
	(cd src; $(CPP) $(CFLAGS) $(DEPENDPATH) utils.cpp)

src/sim_data.o: src/sim_data.hpp src/sim_data.cpp
	(cd src; $(CPP) $(CFLAGS) $(DEPENDPATH) sim_data.cpp)

src/bonner_ebert_sphere.o: src/bonner_ebert_sphere.hpp src/bonner_ebert_sphere.cpp 
	(cd src; $(CPP) $(CFLAGS) $(DEPENDPATH) bonner_ebert_sphere.cpp)

src/grid.o: src/grid.hpp src/grid.cpp
	(cd src; $(CPP) $(CFLAGS) $(DEPENDPATH) grid.cpp)

src/disk_evolution_driver.o: src/disk_evolution_driver.hpp src/disk_evolution_driver.cpp 
	(cd src; $(CPP) $(CFLAGS) $(DEPENDPATH) disk_evolution_driver.cpp)


clean:
	rm run $(objects)
PROBDIR = TM2018
# PROBDIR = TM2018_2
# PROBDIR = TM2018_dustgrowth

CPP = g++
CFLAGS = -g -c -O3 -std=c++17
LFLAGS =
DEPENDPATH = -I.

objects = src/main.o \
		  src/utils.o \
		  src/bonner_ebert_sphere.o \
		  src/grid.o \
		  src/sim_data.o \
		  src/$(PROBDIR)/dust_input.o \
		  src/driver_base.o \
		  src/$(PROBDIR)/driver.o \
		  src/quadpack.o


run: $(objects)
	$(CPP) -o run $(objects) $(LFLAGS)

src/main.o: src/main.cpp
	(cd src; $(CPP) $(CFLAGS) $(DEPENDPATH) main.cpp)

src/utils.o: src/utils.hpp src/utils.cpp 
	(cd src; $(CPP) $(CFLAGS) $(DEPENDPATH) utils.cpp)

src/bonner_ebert_sphere.o: src/bonner_ebert_sphere.hpp src/bonner_ebert_sphere.cpp 
	(cd src; $(CPP) $(CFLAGS) $(DEPENDPATH) bonner_ebert_sphere.cpp)

src/grid.o: src/grid.hpp src/grid.cpp
	(cd src; $(CPP) $(CFLAGS) $(DEPENDPATH) grid.cpp)

src/sim_data.o: src/sim_data.hpp src/sim_data.cpp
	(cd src; $(CPP) $(CFLAGS) $(DEPENDPATH) sim_data.cpp)

src/$(PROBDIR)/dust_input.o: src/sim_data.hpp src/$(PROBDIR)/dust_input.cpp
	(cd src/$(PROBDIR); $(CPP) $(CFLAGS) $(DEPENDPATH) dust_input.cpp)

src/driver_base.o: src/driver_base.hpp src/driver_base.cpp 
	(cd src; $(CPP) $(CFLAGS) $(DEPENDPATH) driver_base.cpp)

src/$(PROBDIR)/driver.o: src/driver.hpp src/$(PROBDIR)/driver.cpp 
	(cd src/$(PROBDIR); $(CPP) $(CFLAGS) $(DEPENDPATH) driver.cpp)

src/quadpack.o: src/quadpack.hpp src/quadpack.cpp
	(cd src; $(CPP) $(CFLAGS) $(DEPENDPATH) quadpack.cpp)

clean:
	rm run $(objects)
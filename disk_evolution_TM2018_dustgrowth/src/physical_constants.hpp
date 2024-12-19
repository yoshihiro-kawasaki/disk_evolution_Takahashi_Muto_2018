#ifndef _PHYSICALCONSTANT_HPP_
#define _PHYSICALCONSTANT_HPP_

#define _USE_MATH_DEFINES
#include <cmath>

namespace physical_constnats
{
    constexpr double SPEED_OF_LIGHT = 2.99792458e10;
    constexpr double GRAVITATIONAL_CONSTANT = 6.67408e-08;
    constexpr double PLANCK_CONSTANT = 6.6260755e-27;
    constexpr double BOLTZMANN_CONSTANT = 1.38064852e-16;
    constexpr double STEFAN_BOLTZMANN_CONSTANT = 5.6705e-5;
    constexpr double ATOMIC_MASS_UNIT = 1.660539067e-24;
    constexpr double ELECTRON_MASS = 9.109383702e-28;
    constexpr double PROTON_MASS = 1.672621924e-24;
    constexpr double NEUTRON_MASS = 1.674927498e-24;
    constexpr double CHARGE_UNIT = 4.803204673e-10;
    constexpr double ELECTRON_VOLT = 1.60218e-12;
    constexpr double BOHR_RADIUS = 5.2917720859e-9;
    constexpr double AVOGADRO_CONSTANT = 6.0221e23;
    constexpr double GAS_CONSTANT_MOL = 8.3145e7;
    constexpr double RADIATION_CONSTANT = 7.5646e-15;
    constexpr double ASTRONOMICAL_UNIT = 1.495979e+13;
    constexpr double PARSEC = 3.085677e18;
    constexpr double LIGHT_YEAR = 9.460730473e17;
    constexpr double SOLAR_YEAR = 3.1556925e7;
    constexpr double SOLAR_MASS = 1.9884e33;
    constexpr double SOLAR_RADIUS = 6.955080e10;
    constexpr double SOLAR_LUMINOSITY = 3.828e33;
    constexpr double GAS_MOLECULAR_MASS = 2.34*PROTON_MASS;
    constexpr double H2_CROSS_SECTION = 2.0e-15;
}

namespace pc = physical_constnats;

#endif /* _PHYSICALCONSTANT_HPP_ */
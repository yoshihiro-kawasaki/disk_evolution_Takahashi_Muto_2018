#ifndef CONSTANTS_HPP_
#define CONSTANTS_HPP_

#define _USE_MATH_DEFINES
#include <cmath>

namespace constnats
/*
    cgs unit
*/
{
    constexpr double SPEED_OF_LIGHT            = 2.99792458e10;
    constexpr double GRAVITATIONAL_CONSTANT    = 6.67408e-08;
    constexpr double PLANCK_CONSTANT           = 6.6260755e-27;
    constexpr double BOLTZMANN_CONSTANT        = 1.38064852e-16;
    constexpr double STEFAN_BOLTZMANN_CONSTANT = 5.6705e-5;
    constexpr double ATOMIC_MASS_UNIT          = 1.660539067e-24;
    constexpr double ELECTRON_MASS             = 9.109383702e-28;
    constexpr double PROTON_MASS               = 1.672621924e-24;
    constexpr double NEUTRON_MASS              = 1.674927498e-24;
    constexpr double CHARGE_UNIT               = 4.803204673e-10;
    constexpr double ELECTRON_VOLT             = 1.60218e-12;
    constexpr double BOHR_RADIUS               = 5.2917720859e-9;
    constexpr double AVOGADRO_CONSTANT         = 6.0221e23;
    constexpr double GAS_CONSTANT_MOL          = 8.3145e7;
    constexpr double RADIATION_CONSTANT        = 7.5646e-15;
    constexpr double ASTRONOMICAL_UNIT         = 1.495979e+13;
    constexpr double PARSEC                    = 3.085677e18;
    constexpr double LIGHT_YEAR                = 9.460730473e17;
    constexpr double SOLAR_YEAR                = 3.1556925e7;
    constexpr double SOLAR_MASS                = 1.9884e33;
    constexpr double SOLAR_RADIUS              = 6.955080e10;
    constexpr double SOLAR_LUMINOSITY          = 3.828e33;
    constexpr double MEAN_MOLECULAR_WEIGHT     = 2.34;
    constexpr double GAS_MOLECULAR_MASS        = MEAN_MOLECULAR_WEIGHT*PROTON_MASS;
    constexpr double H2_CROSS_SECTION          = 2.0e-15;

    constexpr double INV_MG                    = 1.0 / GAS_MOLECULAR_MASS;
    constexpr double SQRT_PI                   = 1.7724538509055159; // sqrt(pi)
    constexpr double SQRT_2PI                  = 2.5066282746310002; // sqrt(2pi)
    constexpr double INV_3PI                   = 1.0 / (3.0 * M_PI); // 1 / (3pi)
    constexpr double INV_4PI                   = 1.0 / (4.0 * M_PI); // 1 / (4pi)
    constexpr double SQRT_8PI_3                = 2.8944050182330705; // sqrt(8pi/3)
    constexpr double PI_8                      = M_PI / 8.0;         // pi / 8
    constexpr double SQRT_8_PI                 = 1.5957691216057308; // sqrt(8/pi)
    constexpr double SQRT_2_PI                 = 0.7978845608028654; // sqrt(2/pi)
    constexpr double INV_SQRT_2PI              = 0.3989422804014327; // 1 / sqrt(2pi)

    constexpr double GAMMA_HEAT                = 1.4;
    constexpr double INV_GAMMA_HEAT_M1         = 1.0 / (GAMMA_HEAT - 1.0);
    constexpr double COE_E_TO_T                = (GAMMA_HEAT - 1.0) * GAS_MOLECULAR_MASS / BOLTZMANN_CONSTANT; 
}

namespace cst = constnats;

#endif /* CONSTANTS_HPP_ */
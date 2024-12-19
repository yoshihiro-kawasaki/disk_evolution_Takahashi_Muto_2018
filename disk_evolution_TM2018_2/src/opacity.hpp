#ifndef OPACITY_HPP_
#define OPACITY_HPP_

#include "physical_constants.hpp"
#include "define_macro.hpp"

#include <iostream>
#include <cmath>

double OpacityModelZhu2012(const double P, const double T);

double OpacityModelSuzukiPlusZhu(const double P, const double T);

#endif /* OPACITY_HPP_ */
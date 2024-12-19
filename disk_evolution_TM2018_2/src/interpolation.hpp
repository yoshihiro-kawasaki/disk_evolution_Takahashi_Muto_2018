#ifndef INTERPOLATION_HPP_
#define INTERPOLATION_HPP_

#include "array.hpp"
#include "define_macro.hpp"
#include "sim_data.hpp"

class Interpolation
{
public:

    Interpolation(SimulationData*pdata);

    double QuadraticInterpolation(const double f1, const double f2, const double f3) {
        return cp_(1)*f1 + cp_(2)*f2 + cp_(3)*f3;
    }


private:

    SimulationData *pdata_;
    Array1D<double> cp_; // 2次(外挿)補間のための係数

};

#endif /* INTERPOLATION_HPP_ */
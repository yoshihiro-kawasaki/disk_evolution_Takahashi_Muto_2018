#ifndef EXPLICIT_SOLVER
#define EXPLICIT_SOLVER

#include "../array.hpp"
#include "../define_macro.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <cstdint>
#include <fstream>
#include <cmath>

typedef void (*EXPLICIT_FUNCTION)(double x, 
                                  Array1D<double> &y, 
                                  Array1D<double> &dydx, 
                                  void *);

class ExplicitOdeSolver
{
public:

    ExplicitOdeSolver(std::string method="rk4");
    ExplicitOdeSolver(int n, std::string method="rk4");
    ~ExplicitOdeSolver();

    void ChangeMethod(std::string method);

    void SetRungeKutta();
    void SetRungeKuttaFehlberg();
    void SetRungeKuttaDormandPrince();

    void RungeKutta(EXPLICIT_FUNCTION func, 
                    double &x, 
                    double h, 
                    Array1D<double> &y, 
                    void *data);

    void RungeKuttaFehlberg(EXPLICIT_FUNCTION f, 
                            double &x, 
                            double &h, 
                            Array1D<double> &y, 
                            double xbound, 
                            int &info, 
                            double tol, 
                            void *data);

    void RungeKuttaDormandPrince(EXPLICIT_FUNCTION f, 
                                 double &x, 
                                 double &h, 
                                 Array1D<double> &y, 
                                 double xbound, 
                                 int &info, 
                                 double tol, 
                                 Array1D<double> &work0, 
                                 void *data);

private:

    int n_;
    std::string method_;
    Array2D<double> a_;
    Array1D<double> b1_, b2_, c_, rc_;
    
};

#endif /* EXPLICIT_SOLVER */
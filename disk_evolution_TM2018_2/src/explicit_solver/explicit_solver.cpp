#include "explicit_solver.hpp"
#include "../array.hpp"

ExplicitOdeSolver::ExplicitOdeSolver(std::string method)
    : method_(method)
/*
    method = rk4    : Runge Kutta (defalut)
           = drkf45 : Runge Kutta Fehlberg
           = ddp45  : Runge Kutta Dormand Prince
*/
{
    if (method == "rk4") {

        SetRungeKutta();

    } else if (method == "drkf45") {

        SetRungeKuttaFehlberg();

    } else if (method == "ddp45" || method == "dDP45") {

        SetRungeKuttaDormandPrince();

    } else {

        std::cerr << "Odeint error : no method " << method << std::endl;
        exit(1);
    }

}


ExplicitOdeSolver::ExplicitOdeSolver(int n, std::string method)
    : n_(n), method_(method)
/*
    n : number of equations
    method = rk4    : Runge Kutta (defalut)
           = drkf45 : Runge Kutta Fehlberg
           = ddp45  : Runge Kutta Dormand Prince
*/
{
    if (method == "rk4") {

        SetRungeKutta();

    } else if (method == "drkf45") {

        SetRungeKuttaFehlberg();

    } else if (method == "ddp45" || method == "dDP45") {

        SetRungeKuttaDormandPrince();

    } else {

        std::cerr << "Odeint error : no method " << method << std::endl;
        exit(1);
    }

} // end Odeint


ExplicitOdeSolver::~ExplicitOdeSolver()
{

}


void ExplicitOdeSolver::ChangeMethod(std::string method)
{
    method_ = method;
    if (method == "rk4") {

        SetRungeKutta();

    } else if (method == "drkf45") {

        SetRungeKuttaFehlberg();

    } else if (method == "ddp45" || method == "dDP45") {

        SetRungeKuttaDormandPrince();

    } else {

        std::cerr << "Odeint error : no method " << method << std::endl;
        exit(1);
    }

    return;
}

void ExplicitOdeSolver::SetRungeKutta()
{
    int s = 4;
    c_.Resize(s+1, 0.0);
    c_(1) = 0.0;
    c_(2) = 0.5;
    c_(3) = 0.5;
    c_(4) = 1.0;

    return;
}


void ExplicitOdeSolver::RungeKutta(EXPLICIT_FUNCTION func, 
                                   double &x, 
                                   double h, 
                                   Array1D<double> &y, 
                                   void *data)
{
    int i, j;
    int n = y.GetDim1(), s = 4;
    Array1D<double> tmp(n, 0.0), f(n, 0.0);
    Array2D<double> K(n, s+1, 0.0);
    double tx;

    if (method_ != "rk4") {
        std::cerr << "ExplicitOdeSolver::RungeKutta error : method = " << method_ << std::endl;
        exit(1);
    }

    for (j = 1; j <= s; ++j) {
        for (i = 0; i < n; ++i) {
            tmp(i) = y(i) + K(i, j-1)*c_(j);
        }
        tx = x + c_(j)*h;
        (*func)(tx, tmp, f, &data);
        for (i = 0; i < n; ++i) {
            K(i, j) = h * f(i);
        }
    }

    x = x + h;
    for(i = 0; i < n; ++i) {
        y(i) = y(i) + (K(i, 1) + K(i, 4))/6.0 + (K(i, 2) + K(i, 3))/3.0;
    }

    return;
}

void ExplicitOdeSolver::SetRungeKuttaFehlberg()
{
    int s = 6;

    a_.Resize(s, s, 0.0);
    b1_.Resize(s, 0.0);
    b2_.Resize(s, 0.0);
    c_.Resize(s, 0.0);
    rc_.Resize(s, 0.0);

    a_(1, 0) = 0.25;
    a_(2, 0) = 0.09375;
    a_(2, 1) = 0.28125;
    a_(3, 0) = 0.8793809740555302685480200273099681383705;
    a_(3, 1) = -3.277196176604460628129267182521620391443;
    a_(3, 2) = 3.320892125625853436504324078288575329995;
    a_(4, 0) = 2.032407407407407407407407407407407407407;
    a_(4, 1) = -8.0;
    a_(4, 2) = 7.173489278752436647173489278752436647173;
    a_(4, 3) = -0.2058966861598440545808966861598440545809;
    a_(5, 0) = -0.2962962962962962962962962962962962962963;
    a_(5, 1) = 2.0;
    a_(5, 2) = -1.381676413255360623781676413255360623782;
    a_(5, 3) = 0.4529727095516569200779727095516569200780;
    a_(5, 4) = -0.275;

    b1_(0) = 0.1157407407407407407407407407407407407407;
    b1_(2) = 0.5489278752436647173489278752436647173489;
    b1_(3) = 0.5353313840155945419103313840155945419103;
    b1_(4) = -0.2;

    b2_(0) = 0.1185185185185185185185185185185185185185;
    b2_(2) = 0.5189863547758284600389863547758284600390;
    b2_(3) = 0.5061314903420166578061314903420166578061;
    b2_(4) = -0.18;
    b2_(5) = 0.03636363636363636363636363636363636363636;

    c_(1) = 0.25;
    c_(2) = 0.375;
    c_(3) = 0.9230769230769230769230769230769230769231;
    c_(4) = 1.0;
    c_(5) = 0.50;

    rc_(0) = 0.002777777777777777777777777777777777777778;
    rc_(2) = -0.02994152046783625730994152046783625730994;
    rc_(3) = -0.02919989367357788410419989367357788410420;
    rc_(4) = 0.02;
    rc_(5) = 0.03636363636363636363636363636363636363636;

    return;
}


void ExplicitOdeSolver::RungeKuttaFehlberg(EXPLICIT_FUNCTION f, 
                                           double &x, 
                                           double &h, 
                                           Array1D<double> &y, 
                                           double xbound, 
                                           int &info, 
                                           double tol, 
                                           void *data)
{
    int n = y.GetDim1(), s = 6;
    int i, j, k, flag, key;
    double r, rr, delta, tx, sy, err;
    double hmin = 1.0e-14, hmax = 0.50;
    Array1D<double> ty(n , 0.0), tf(n, 0.0);
    Array2D<double> K(s, n, 0.0);

    if (method_ != "drkf45") {
        std::cerr << "ExplicitOdeSolver::RungeKuttaFehlberg error : method = " << method_ << std::endl;
        exit(1);
    }

    key = 0;

    if (std::abs(h) >= hmax) h = SIGN(h) * hmax;
    if (h >= std::abs(xbound-x)) h = xbound - x;

    flag = 1;
    if (std::abs(x - xbound) <= hmin) {
        info = 1;
        flag = 0;
    }

    while (flag == 1) {
        tx = x; 
        for (j = 0; j < s; ++j) {
            tx = x + c_(j)*h;
            ty = y;
            for (i = 0; i <= j-1; ++i) {
                for (k = 0; k < n; ++k) {
                    ty(k) = ty(k) + K(i, k)*a_(j, i);
                }
            }
            (*f)(tx, ty, tf, &data);
            for (k = 0; k < n; ++k) {
                K(j, k) = h * tf(k);
            }
        }

        // step4
        r = 0.0;
        for (i = 0; i < n; ++i) {
            rr = rc_(0)*K(0, i) + rc_(2)*K(2, i) + rc_(3)*K(3, i) + rc_(4)*K(4, i) 
                + rc_(5)*K(5, i);
            r = r + SQR(rr);
        }
        r = std::abs(sqrt(r)/h/double(n));

        sy = 0.0;
        for (i = 0; i < n; ++i) sy = sy + y(i)*y(i);
        sy = sqrt(sy);
        if (sy >= 1.0) {
            err = tol * sy;
        } else {
            err = tol;
        }

        // step 5
        if (r <= err || key == 1) {
            x = x + h;
            for (i = 0; i < n; ++i) {
                y(i) = y(i) + b1_(0)*K(0, i) + b1_(2)*K(2, i) + b1_(3)*K(3, i) 
                    + b1_(4)*K(4, i);
            }
            flag = 0;
        }

        // step 6
        if (r >= 1.0e-20) {
            delta = pow(err/(2.0*r), 0.25);
        } else {
            delta = 4.0;
        }

        // step 7
        if (delta <= 0.1) {
            h = 0.1 * h;
        } else if (delta >= 4.0) {
            delta = 4.0 * h;
        } else {
            h = delta * h;
        }

        // step 8
        if (std::abs(h) >= hmax) {
            h = SIGN(h) * hmax;
        } else if (std::abs(h) < hmin) {
            h = SIGN(h) * hmin;
            key = 1;
        } 

        // step 9
        if (std::abs(xbound - x) <= std::abs(h)) {
            h = xbound - x;
            if (std::abs(h) <= hmin) {
                info = 1;
                flag = 0;
            }
        }

        if (h <= 0.0 && (xbound - x) >= 0.0) {
            info = 1;
            flag = 0;
        } else if (h >= 0.0 && (xbound - x) <= 0.0) {
            info = 1;
            flag = 0;
        }

    } // end while 

    if (key == 1) {
        std::cout << "Strange point between " << (x-h) << " and " << x << std::endl;
        info = -9; 
    }

    return;
}


void ExplicitOdeSolver::SetRungeKuttaDormandPrince()
{
    int s = 7;

    a_.Resize(s, s, 0.0);
    b1_.Resize(s, 0.0);
    b2_.Resize(s, 0.0);
    c_.Resize(s, 0.0);
    rc_.Resize(s, 0.0);

    a_(1, 0) = 0.2;
    a_(2, 0) = 0.075;
    a_(2, 1) = 0.225;
    a_(3, 0) = 0.977777777777777777777777777777777777778;
    a_(3, 1) = -3.733333333333333333333333333333333333333;
    a_(3, 2) = 3.555555555555555555555555555555555555556;
    a_(4, 0) = 2.95259868922420362749580856576741350403902;
    a_(4, 1) = -11.59579332418838591678097850937357110196616;
    a_(4, 2) = 9.8228928516994360615759792714525224813290657;
    a_(4, 3) = -0.29080932784636488340192043895747599451303155;
    a_(5, 0) = 2.846275252525252525252525252525252525252525;
    a_(5, 1) = -10.757575757575757575757575757575757575757576;
    a_(5, 2) = 8.906422717743472460453592529064227177434725;
    a_(5, 3) = 0.278409090909090909090909090909090909090909;
    a_(5, 4) = -0.273531303602058319039451114922813036020583;
    a_(6, 0) = 0.09114583333333333333333333333333333333333333;
    a_(6, 2) = 0.4492362982929020664869721473495058400718778077;
    a_(6, 3) = 0.6510416666666666666666666666666666666666666667;
    a_(6, 4) = -0.3223761792452830188679245283018867924528301887;
    a_(6, 5) = 0.13095238095238095238095238095238095238095238095;

    for (int i = 0; i < s; ++i) b1_(i) = a_(6, i);

    b2_(0) = 0.0899131944444444444444444444444444444444444444;
    b2_(2) = 0.45348906858340820604971548367774782869122491764;
    b2_(3) = 0.6140625;
    b2_(4) = -0.271512382075471698113207547169811320754716981;
    b2_(5) = 0.089047619047619047619047619047619047619047619;
    b2_(6) = 0.025;

    c_(1) = 0.2;
    c_(2) = 0.3;
    c_(3) = 0.8;
    c_(4) = 0.888888888888888888888888888888888888888889;
    c_(5) = 1.0;
    c_(6) = 1.0;

    rc_(0) = 0.001232638888888888888888888888888888888888888889;
    rc_(2) = -0.004252770290506139562743336328241988619347109913;
    rc_(3) = 0.036979166666666666666666666666666666666666666;
    rc_(4) = -0.050863797169811320754716981132075471698113207547;
    rc_(5) = 0.0419047619047619047619047619047619047619047619;
    rc_(6) = -0.025;

    return;
}

void ExplicitOdeSolver::RungeKuttaDormandPrince(EXPLICIT_FUNCTION f, 
                                                double &x, 
                                                double &h, 
                                                Array1D<double> &y, 
                                                double xbound, 
                                                int &info, 
                                                double tol, 
                                                Array1D<double> &work0, 
                                                void *data)
{
    int n = y.GetDim1(), s = 7;
    int i, j, k, flag, key;
    double r, rr, delta, tx, sy, err;
    double hmin = 1.0e-12, hmax = 0.50;
    Array1D<double> tmp(n , 0.0), tf(n, 0.0), work(n, 0.0);
    Array2D<double> K(s, n, 0.0);

    if (method_ != "ddp45" && method_ != "dDP45") {
        std::cerr << "ExplicitOdeSolver::RungeKuttaDormandPrince error : method = " << method_ << std::endl;
        exit(1);
    }

    key = 0;

    if (std::abs(h) >= hmax) h = SIGN(h) * hmax;
    if (h >= std::abs(xbound-x)) h = xbound - x;

    flag = 1;
    if (std::abs(x - xbound) <= hmin) {
        info = 1;
        flag = 0;
    }

    while (flag == 1) {
        tx = x;

        if (info != -1) {
            for (i = 0; i < n; ++i) K(0, i) = h * work0(i);
        } else {
            (*f)(x, y, work, &data);
            for (i = 0; i < n; ++i) {
                K(0, i) = h * work(i);
                work0(i) = work(i);
            }
            info = 0;
        }

        for (i = 1; i < s; ++i) {
            tx = x + c_(i)*h;
            tmp = y;
            for (j = 0; j <= i-1; ++j) {
                for (k = 0; k < n; ++k) {
                    tmp(k) = tmp(k) + K(j, k)*a_(i, j);
                }
            }
            (*f)(tx, tmp, work, &data);
            for (j = 0; j < n; ++j) {
                K(i, j) = h * work(j);
            }
        }

        // step 4
        r = 0.0;
        for (i = 0; i < n; ++i) {
            rr = rc_(0)*K(0, i) + rc_(2)*K(2, i) + rc_(3)*K(3, i) + rc_(4)*K(4, i) 
                + rc_(5)*K(5, i) + rc_(6)*K(6, i);
            r = r + SQR(rr);
        }
        r = std::abs(sqrt(r)/h/double(n));

        sy = 0.0;
        for (i = 0; i < n; ++i) sy = sy + y(i)*y(i);
        sy = sqrt(sy);
        if (sy >= 1.0) {
            err = tol * sy;
        } else {
            err = tol;
        }

        // step 5
        if (r <= err || key == 1) {
            x = x + h;
            for (i = 0; i < n; ++i) {
                y(i) = y(i) + b1_(0)*K(0, i) + b1_(2)*K(2, i) + b1_(3)*K(3, i) + b1_(4)*K(4, i)
                    + b1_(5)*K(5, i);
                work0(i) = work(i);
            }
            flag = 0;
        }

        // step 6
        if (r >= 1.0e-20) {
            delta = pow(err/(2.0*r), 0.20);
        } else {
            delta = 4.0;
        }

        // step 7
        if (delta <= 0.1) {
            h = 0.1 * h;
        } else if (delta >= 4.0) {
            delta = 4.0 * h;
        } else {
            h = delta * h;
        }

        // step 8
        if (std::abs(h) >= hmax) {
            h = SIGN(h) * hmax;
        } else if (std::abs(h) < hmin) {
            h = SIGN(h) * hmin;
            key = 1;
        } 

        // step 9
        if (std::abs(xbound - x) <= std::abs(h)) {
            h = xbound - x;
            if (std::abs(h) <= hmin) {
                info = 1;
                flag = 0;
            }
        }

        if (h <= 0.0 && (xbound - x) >= 0.0) {
            info = 1;
            flag = 0;
        } else if (h >= 0.0 && (xbound - x) <= 0.0) {
            info = 1;
            flag = 0;
        }

    } // end while

    if (key == 1) {
        std::cout << "Strange point between " << (x-h) << " and " << x << std::endl;
        info = -2; 
    }
    return;
}
/**
 * @file quadpack.cpp
 * @brief 
 * @author yoshihiro kawasaki
 * @date 2025/01/14
 * @details 
 * @note
 *      https://www.netlib.org/quadpack/
 *      https://github.com/jacobwilliams/quadpack
*/

#include "quadpack.hpp"

namespace quadpack_cpp
{

/**
* @def ARRAYF
* @brief C++ array index --> Fortran array index
* @details
*/
#define ARRAYF(a, i) (a[(i)-1])

/**
 * @def MATF
 * @brief C++ array index --> Fortran 2d-array index
 * @details
 *  a[nrows][ncols] = a[nrows * ncols]
 *  a[i][j] = a[i*ncols + j]
*/
#define MATF(a, ncols,  i, j) (a[((i) - 1) * (ncols) + ((j) - 1)])

// ******************************************************************************************************************* //
// ******************************************************************************************************************* //
// ******************************************************************************************************************* //

using QUADPACK_CPP_WEIGHT_FUNCTION = double(*)(const double x, const double a, const double b, 
                                              const double c, const double d, const int i);

// ******************************************************************************************************************* //
// ******************************************************************************************************************* //
// ******************************************************************************************************************* //

const double d1mach[5] = {
    std::numeric_limits<double>::min(),
    std::numeric_limits<double>::max(),
    std::pow(static_cast<double>(std::numeric_limits<double>::radix), -std::numeric_limits<double>::digits),
    std::numeric_limits<double>::epsilon(),
    std::log10(static_cast<double>(std::numeric_limits<double>::radix))
};

// ******************************************************************************************************************* //
// ******************************************************************************************************************* //
// ******************************************************************************************************************* //

/**
 * @fn void dgtsl
 * @brief dgtsl given a general tridiagonal matrix and a right hand side will find the solution.
 *
 * @param [in] int n 
 *      is the order of the tridiagonal matrix.
 *
 * @param [in] double* c
 *      is the subdiagonal of the tridiagonal matrix.
 *      c(2) through c(n) should contain the subdiagonal.
 *      on output c is destroyed.
 *   
 * @param [in] double *d
 *      is the diagonal of the tridiagonal matrix.
 *      on output d is destroyed.
 *
 * @param [in] double *e
 *      is the superdiagonal of the tridiagonal matrix.
 *      e(1) through e(n-1) should contain the superdiagonal.
 *      on output e is destroyed.
 *
 * @param [inout] double *b
 *      is the right hand side vector.
 *
 * @param [inout] int &info
 *      info = 0 normal value.
 *      info = k if the k-th element of the diagonal becomes
 *             exactly zero.  the subroutine returns when
 *             this is detected.
 *
 * @note 
 *      linpack. this version dated 08/14/78 .
 *      jack dongarra, argonne national laboratory.
*/
void dgtsl(const int n, double *c, double *d, double *e, double *b, int &info);

void dqc25c(QUADPACK_CPP_FUNCTION f, const double a, const double b, const double c, 
            double &result, double &abserr, int &krul, int &neval, void *user_data);

void dqc25f(QUADPACK_CPP_FUNCTION f, const double a, const double b, const double omega, const int integr, 
            const int nrmom, const int maxp1, const int ksave, double &result, double &abserr, int &neval, 
            double &resabs, double &resasc, int &momcom, double *chebmo, void *user_data);

void dqc25s(QUADPACK_CPP_FUNCTION f, const double a, const double b, const double bl, const double br,
            const double alfa, const double beta, const double *ri, const double *rj, const double *rg,
            const double *rh, double &result, double &abserr, double &resasc, const int integr, int &nev,
            void *user_data);

void dqcheb(const double *x, double *fval, double *cheb12, double *cheb24);

void dqelg(int &n, double *epstab, double &result, double &abserr, double *res3la, int &nres);

void dqk15(QUADPACK_CPP_FUNCTION f, const double a, const double b, double &result, double &abserr,
           double &resabs, double &resasc, void *user_data);

void dqk15i(QUADPACK_CPP_FUNCTION f, const double boun, const int inf, const double a, const double b,
            double &result, double &abserr, double &resabs, double &resasc, void *user_data);

void dqk15w(QUADPACK_CPP_FUNCTION f, QUADPACK_CPP_WEIGHT_FUNCTION w, const double p1, const double p2,
            const double p3, const double p4, const int kp, const double a, const double b, double &result,
            double &abserr, double &resabs, double &resasc, void *user_data);

void dqk21(QUADPACK_CPP_FUNCTION f, const double a, const double b, double &result, double &abserr,
           double &resabs, double &resasc, void *user_data);

void dqk31(QUADPACK_CPP_FUNCTION f, const double a, const double b, double &result, double &abserr,
           double &resabs, double &resasc, void *user_data);

void dqk41(QUADPACK_CPP_FUNCTION f, const double a, const double b, double &result, double &abserr,
           double &resabs, double &resasc, void *user_data);

void dqk51(QUADPACK_CPP_FUNCTION f, const double a, const double b, double &result, double &abserr,
           double &resabs, double &resasc, void *user_data);

void dqk61(QUADPACK_CPP_FUNCTION f, const double a, const double b, double &result, double &abserr,
           double &resabs, double &resasc, void *user_data);

void dqmomo(const double alfa, const double beta, double *ri, double *rj, double *rg, double *rh, const int integr);

void dqpsrt(const int limit, const int last, int &maxerr, double &ermax, const double *elist, int *iord, int &nrmax);

double dqwgtc(const double x, const double c, const double p2, const double p3, const double p4, const int kp);

double dqwgtf(const double x, const double omega, const double p2, const double p3, const double p4, const int integr);

double dqwgts(const double x, const double a, const double b, const double alfa, const double beta, const int integr);

void xerror(const std::string messg, const int nerr, const int level);

// ******************************************************************************************************************* //
// ******************************************************************************************************************* //
// ******************************************************************************************************************* //

/**
 * @fn void dgag
*/
void dqag(QUADPACK_CPP_FUNCTION f, const double a, const double b, const double epsabs, const double epsrel,
          const int key, double &result, double &abserr, int &neval, int &ier, const int limit, 
          const int lenw, int &last, int *iwork, double *work, void *user_data)
{
    int lvl, l1, l2, l3;
    //
    // check validity of lenw.
    //
    // first executable statement  dqag
    ier    = 6;
    neval  = 0;
    last   = 0;
    result = 0.0;
    abserr = 0.0;
    if (limit < 1 || lenw < limit*4) goto LABEL_10;
    //
    // prepare call for dqage.
    //
    l1 = limit + 1;
    l2 = limit + l1;
    l3 = limit + l2;
    //
    dqage(f, a, b, epsabs, epsrel, key, limit, result, abserr, neval, ier,
        &ARRAYF(work, 1) ,&ARRAYF(work, l1) ,&ARRAYF(work, l2) ,&ARRAYF(work, l3), iwork, last, user_data);
    //
    // call error handler if necessary.
    //
    lvl = 0;
LABEL_10:
    if (ier == 6) lvl = 1;
    if (ier != 0) {
        std::string messg = "abnormal return from dqag";
        xerror(messg, ier,lvl);
    }
    return;
}

/**
 * @fn void dgage
*/
void dqage(QUADPACK_CPP_FUNCTION f, const double a, const double b, const double epsabs, const double epsrel,
           const int key, const int limit, double &result, double &abserr, int &neval, int &ier, 
           double *alist, double *blist, double *rlist, double *elist, int *iord, int &last, void *user_data)
{
    double area1, a1, b1, defab1, error1;
    double area2, a2, b2, defab2, error2;
    double area, area12, erro12, errsum, errmax, errbnd, defabs, resabs;
    int maxerr, iroff1, iroff2, k, keyf, nrmax;
    double epmach, uflow;
    // first executable statement  dqage
    epmach = ARRAYF(d1mach, 4);
    uflow  = ARRAYF(d1mach, 1);
    //
    // test on validity of parameters
    // ------------------------------
    //
    ier              = 0;
    neval            = 0;
    last             = 0;
    result           = 0.0;
    abserr           = 0.0;
    ARRAYF(alist, 1) = a;
    ARRAYF(blist, 1) = b;
    ARRAYF(rlist, 1) = 0.0;
    ARRAYF(elist, 1) = 0.0;
    if (epsabs <= 0.0 && epsrel < std::max(50.0*epmach, 5.0e-29)) ier = 6;
    if (ier == 6) goto LABEL_999;
    //
    // first approximation to the integral
    // -----------------------------------
    //
    keyf = key;
    if (key <= 0) keyf = 1;
    if (key >= 7) keyf = 6;
    neval = 0;
    if (keyf == 1) dqk15(f, a, b, result, abserr, defabs, resabs, user_data);
    if (keyf == 2) dqk21(f, a, b, result, abserr, defabs, resabs, user_data);
    if (keyf == 3) dqk31(f, a, b, result, abserr, defabs, resabs, user_data);
    if (keyf == 4) dqk41(f, a, b, result, abserr, defabs, resabs, user_data);
    if (keyf == 5) dqk51(f, a, b, result, abserr, defabs, resabs, user_data);
    if (keyf == 6) dqk61(f, a, b, result, abserr, defabs, resabs, user_data);
    last             = 1;
    ARRAYF(rlist, 1) = result;
    ARRAYF(elist, 1) = abserr;
    ARRAYF(iord,  1) = 1;
    //
    // test on accuracy.
    //
    errbnd = std::max(epsabs, epsrel * std::abs(result));
    if (abserr <= 50.0 * epmach * defabs && abserr > errbnd) ier = 2;
    if(limit == 1) ier = 1;
    if(ier != 0 || (abserr <= errbnd && abserr != resabs) || abserr == 0.0) goto LABEL_60;
    //
    // initialization
    // --------------
    //
    errmax = abserr;
    maxerr = 1;
    area   = result;
    errsum = abserr;
    nrmax  = 1;
    iroff1 = 0;
    iroff2 = 0;
    //
    // main for-loop
    // -------------
    //
    for (last = 2; last <= limit; ++last) {
        //
        // bisect the subinterval with the largest error estimate.
        //
        a1 = ARRAYF(alist, maxerr);
        b1 = 0.5 * (ARRAYF(alist, maxerr) + ARRAYF(blist, maxerr));
        a2 = b1;
        b2 = ARRAYF(blist, maxerr);
        if (keyf == 1) dqk15(f, a1, b1, area1, error1, resabs, defab1, user_data);
        if (keyf == 2) dqk21(f, a1, b1, area1, error1, resabs, defab1, user_data);
        if (keyf == 3) dqk31(f, a1, b1, area1, error1, resabs, defab1, user_data);
        if (keyf == 4) dqk41(f, a1, b1, area1, error1, resabs, defab1, user_data);
        if (keyf == 5) dqk51(f, a1, b1, area1, error1, resabs, defab1, user_data);
        if (keyf == 6) dqk61(f, a1, b1, area1, error1, resabs, defab1, user_data);
        if (keyf == 1) dqk15(f, a2, b2, area2, error2, resabs, defab2, user_data);
        if (keyf == 2) dqk21(f, a2, b2, area2, error2, resabs, defab2, user_data);
        if (keyf == 3) dqk31(f, a2, b2, area2, error2, resabs, defab2, user_data);
        if (keyf == 4) dqk41(f, a2, b2, area2, error2, resabs, defab2, user_data);
        if (keyf == 5) dqk51(f, a2, b2, area2, error2, resabs, defab2, user_data);
        if (keyf == 6) dqk61(f, a2, b2, area2, error2, resabs, defab2, user_data);
        //
        // improve previous approximations to integral
        // and error and test for accuracy.
        //
        neval++;
        area12 = area1  + area2;
        erro12 = error1 + error2;
        errsum = errsum + erro12 - errmax;
        area   = area   + area12 - ARRAYF(rlist, maxerr);
        if (defab1 == error1 || defab2 == error2) goto LABEL_5;
        if (std::abs(ARRAYF(rlist, maxerr) - area12) <= 1.0e-5 * std::abs(area12) && erro12 >= 0.99 * errmax) iroff1++;
        if (last > 10 && erro12 > errmax) iroff2++;
LABEL_5:
        ARRAYF(rlist, maxerr) = area1;
        ARRAYF(rlist, last)   = area2;
        errbnd                = std::max(epsabs, epsrel * std::abs(area));
        if (errsum <= errbnd) goto LABEL_8;
        //
        // test for roundoff error and eventually set error flag.
        //
        if (iroff1 >= 6 || iroff2 >= 20) ier = 2;
        //
        // set error flag in the case that the number of subintervals
        // equals limit.
        //
        if (last == limit) ier = 1;
        //
        // set error flag in the case of bad integrand behaviour
        // at a point of the integration range.
        //
        if(std::max(std::abs(a1), std::abs(b2)) <= (1.0 + 1.0e2 * epmach) * (std::abs(a2) + 1.0e3 * uflow)) ier = 3;
        //
        // append the newly-created intervals to the list.
        //
LABEL_8:
        if (error2 > error1) goto LABEL_10;
        ARRAYF(alist, last)   = a2;
        ARRAYF(blist, maxerr) = b1;
        ARRAYF(blist, last)   = b2;
        ARRAYF(elist, maxerr) = error1;
        ARRAYF(elist, last)   = error2;
        goto LABEL_20;
LABEL_10:
        ARRAYF(alist, maxerr) = a2;
        ARRAYF(alist, last)   = a1;
        ARRAYF(blist, last)   = b1;
        ARRAYF(rlist, maxerr) = area2;
        ARRAYF(rlist, last)   = area1;
        ARRAYF(elist, maxerr) = error2;
        ARRAYF(elist, last)   = error1;
        //
        // call subroutine dqpsrt to maintain the descending ordering
        // in the list of error estimates and select the subinterval
        // with the largest error estimate (to be bisected next).
        //
LABEL_20:
        dqpsrt(limit, last, maxerr, errmax, elist, iord, nrmax);
        // jump out of for-loop
        if(ier != 0.0 || errsum <= errbnd) goto LABEL_40;
    }
    //
    // compute final result.
    // ---------------------
    //
LABEL_40:
    result = 0.0;
    for (k = 1; k <= last; ++k) {
        result += ARRAYF(rlist, k);
    }
    abserr = errsum;
LABEL_60:
    if (keyf != 1) neval = (10 * keyf + 1) * (2 * neval + 1);
    if (keyf == 1) neval = 30 * neval + 15;
LABEL_999:
    return;
}

/**
 * @fn void dgagi
*/
void dqagi(QUADPACK_CPP_FUNCTION f, const double bound, const int inf, const double epsabs, const double epsrel,
           double &result, double &abserr, int &neval, int &ier, const int limit, const int lenw,
           int &last, int *iwork, double *work, void *user_data)
{
    int lvl, l1, l2, l3;
    // first executable statement  dqagi
    ier    = 6;
    neval  = 0;
    last   = 0;
    result = 0.0;
    abserr = 0.0;
    if(limit < 1 || lenw < limit*4) goto LABEL_10;
    //
    // prepare call for dqagie.
    //
    l1 = limit + 1;
    l2 = limit + l1;
    l3 = limit + l2;
    //
    dqagie(f, bound, inf, epsabs, epsrel, limit, result, abserr, neval, ier,
        &ARRAYF(work, 1) ,&ARRAYF(work, l1) ,&ARRAYF(work, l2) ,&ARRAYF(work, l3), iwork, last, user_data);
    //
    // call error handler if necessary.
    //
    lvl = 0;
LABEL_10:
    if (ier == 6) lvl = 1;
    if (ier != 0) {
        std::string messg = "abnormal return from dqagi";
        xerror(messg, ier, lvl);   
    }
    return;
}

/**
 * @fn void dgagie
*/
void dqagie(QUADPACK_CPP_FUNCTION f, const double bound, const int inf, const double epsabs, const double epsrel,
            const int limit, double &result, double &abserr, int &neval, int &ier, double *alist,
            double *blist, double *rlist, double *elist, int *iord, int &last, void *user_data)
{
    int id, ierro, iroff1, iroff2, iroff3, jupbnd, k, ksgn, ktmin, L, maxerr, nres, nrmax, numrl2;
    bool extrap, noext;
    double a1, a2, abseps, area, area1, area12, area2, b1, b2, boun, correc, defabs;
    double defab1, defab2, dres, erlarg, erlast, errbnd, errmax, error1, error2;
    double erro12, errsum, ertest, resabs, reseps, small;
    double rlist2[52], res3la[3];
    double epmach, uflow, oflow;
    // first executable statement  dqagie
    epmach = ARRAYF(d1mach, 4);
    //
    // test on validity of parameters
    // -----------------------------
    //
    ier              = 0;
    neval            = 0;
    last             = 0;
    result           = 0.0;
    abserr           = 0.0;
    ARRAYF(alist, 1) = 0.0;
    ARRAYF(blist, 1) = 1.0;
    ARRAYF(rlist, 1) = 0.0;
    ARRAYF(elist, 1) = 0.0;
    ARRAYF(iord,  1) = 0;
    if (epsabs <= 0.0 && epsrel < std::max(50.0 * epmach, 5.0e-29)) ier = 6;
    if (ier == 6) goto LABEL_999;
    //
    //
    // first approximation to the integral
    // -----------------------------------
    //
    // determine the interval to be mapped onto (0,1).
    // if inf = 2 the integral is computed as i = i1+i2, where
    // i1 = integral of f over (-infinity,0),
    // i2 = integral of f over (0,+infinity).
    //
    boun = bound;
    if (inf == 2) boun = 0.0;
    dqk15i(f, boun, inf, 0.0, 1.0, result, abserr, defabs, resabs, user_data);
    //
    // test on accuracy
    //
    last             = 1;
    ARRAYF(rlist, 1) = result;
    ARRAYF(elist, 1) = abserr;
    ARRAYF(iord,  1) = 1;
    dres             = std::abs(result);
    errbnd           = std::max(epsabs, epsrel * dres);
    if (abserr <= 1.0e2 * epmach * defabs && abserr > errbnd) ier = 2;
    if (limit == 1) ier = 1;
    if (ier != 0 || (abserr <= errbnd && abserr != resabs) || abserr == 0.0) goto LABEL_130;
    //
    // initialization
    // --------------
    //
    uflow             = ARRAYF(d1mach, 1);
    oflow             = ARRAYF(d1mach, 2);
    ARRAYF(rlist2, 1) = result;
    errmax            = abserr;
    maxerr            = 1;
    area              = result;
    errsum            = abserr;
    abserr            = oflow;
    nrmax             = 1;
    nres              = 0;
    ktmin             = 0;
    numrl2            = 2;
    extrap            = false;
    noext             = false;
    ierro             = 0;
    iroff1            = 0;
    iroff2            = 0;
    iroff3            = 0;
    ksgn              = -1;
    if(dres >= (1.0 - 50.0*epmach)*defabs) ksgn = 1;
    //
    // main for-loop
    // -------------
    //
    for (last = 2; last <= limit; ++last) {
        //
        // bisect the subinterval with the nrmax-th largest error estimate.
        //
        a1 = ARRAYF(alist, maxerr);
        b1 = 0.5 * (ARRAYF(alist, maxerr) + ARRAYF(blist, maxerr));
        a2 = b1;
        b2 = ARRAYF(blist, maxerr);
        erlast = errmax;
        dqk15i(f, boun, inf, a1, b1, area1, error1, resabs, defab1, user_data);
        dqk15i(f, boun, inf, a2, b2, area2, error2, resabs, defab2, user_data);
        // 
        // improve previous approximations to integral 
        // and error and test for accuracy.
        //
        area12 = area1  + area2;
        erro12 = error1 + error2;
        errsum = errsum + erro12 - errmax;
        area   = area   + area12 - ARRAYF(rlist, maxerr);
        if (defab1 == error1 || defab2 == error2) goto LABEL_15;
        if (std::abs(ARRAYF(rlist, maxerr) - area12) > 1.0e-5 * std::abs(area12) || erro12 < 0.99 * errmax) goto LABEL_10;
        if (extrap)  iroff2++;
        if (!extrap) iroff1++;
LABEL_10:
        if (last > 10 && erro12 > errmax) iroff3++;
LABEL_15:
        ARRAYF(rlist, maxerr) = area1;
        ARRAYF(rlist,   last) = area2;
        errbnd                = std::max(epsabs, epsrel * std::abs(area));
        // 
        // test for roundoff error and eventually set error flag.
        //
        if (iroff1 + iroff2 >= 10 || iroff3 >= 20) ier = 2;
        if (iroff2 >= 5) ierro = 3;
        //
        // set error flag in the case that the number of subintervals 
        // equals limit.
        //
        if (last == limit) ier = 1;
        //
        // set error flag in the case of bad integrand behaviour 
        // at a point of the integration range.
        //
        if (std::max(std::abs(a1), std::abs(b2)) <= (1.0 + 1.0e2 * epmach) * (std::abs(a2) + 1.0e3 * uflow)) ier = 4;
        //
        // append the newly-created intervals to the list.
        //
        if (error2 > error1) goto LABEL_20;
        ARRAYF(alist, last)   = a2;
        ARRAYF(blist, maxerr) = b1;
        ARRAYF(blist, last)   = b2;
        ARRAYF(elist, maxerr) = error1;
        ARRAYF(elist, last)   = error2;
        goto LABEL_30;
LABEL_20:
        ARRAYF(alist, maxerr) = a2;
        ARRAYF(alist, last)   = a1;
        ARRAYF(blist, last)   = b1;
        ARRAYF(rlist, maxerr) = area2;
        ARRAYF(rlist, last)   = area1;
        ARRAYF(elist, maxerr) = error2;
        ARRAYF(elist, last)   = error1;
        //
        // call subroutine dqpsrt to maintain the descending ordering
        // in the list of error estimates and select the subinterval
        // with nrmax-th largest error estimate (to be bisected next).
        //
LABEL_30:
        dqpsrt(limit, last, maxerr, errmax, elist, iord, nrmax);
        if (errsum <= errbnd) goto LABEL_115;
        if (ier != 0)  goto LABEL_100;
        if (last == 2) goto LABEL_80;
        if (noext) goto LABEL_90;
        erlarg -= erlast;
        if (std::abs(b1 - a1) > small) erlarg += erro12;
        if (extrap) goto LABEL_40;
        // 
        // test whether the interval to be bisected next is the
        // smallest interval.
        // 
        if (std::abs(ARRAYF(blist, maxerr) - ARRAYF(alist, maxerr)) > small) goto LABEL_90;
        extrap = true;
        nrmax  = 2;
LABEL_40:
        if(ierro == 3 || erlarg <= ertest) goto LABEL_60;
        //
        // the smallest interval has the largest error.
        // before bisecting decrease the sum of the errors over the
        // larger intervals (erlarg) and perform extrapolation.
        // 
        id     = nrmax;
        jupbnd = last;
        if (last > (2 + limit / 2)) jupbnd = limit + 3 - last;
        for (k = id; k <= jupbnd; ++k) {
            maxerr = ARRAYF(iord, nrmax);
            errmax = ARRAYF(elist, maxerr);
            if (std::abs(ARRAYF(blist, maxerr) - ARRAYF(alist, maxerr)) > small) goto LABEL_90;
            nrmax++;
        }
        //
        // perform extrapolation.
        // 
LABEL_60:
        numrl2++;
        ARRAYF(rlist2, numrl2) = area;
        dqelg(numrl2, rlist2, reseps, abseps, res3la, nres);
        ktmin++;
        if(ktmin > 5 && abserr < 1.0e-3 * errsum) ier = 5;
        if(abseps >= abserr) goto LABEL_70;
        ktmin  = 0;
        abserr = abseps;
        result = reseps;
        correc = erlarg;
        ertest = std::max(epsabs, epsrel * std::abs(reseps));
        if (abserr <= ertest) goto LABEL_100;
        //
        // prepare bisection of the smallest interval.
        // 
LABEL_70:
        if (numrl2 == 1) noext = true;
        if (ier == 5) goto LABEL_100;
        maxerr = ARRAYF(iord, 1);
        errmax = ARRAYF(elist, maxerr);
        nrmax  = 1;
        extrap = false;
        small  = small * 0.5;
        erlarg = errsum;
        goto LABEL_90;
LABEL_80:
        small  = 0.375;
        erlarg = errsum;
        ertest = errbnd;
        ARRAYF(rlist2, 2) = area;
LABEL_90:
        continue;
    }
    //
    // set final result and error estimate.
    //
LABEL_100:
    if (abserr == oflow) goto LABEL_115;
    if (ier + ierro == 0) goto LABEL_110;
    if (ierro == 3) abserr += correc;
    if (ier == 0) ier = 3;
    if (result != 0.0 && area != 0.0) goto LABEL_105;
    if (abserr > errsum) goto LABEL_115;
    if (area == 0.0) goto LABEL_130;
    goto LABEL_110;
LABEL_105:
    if (abserr / std::abs(result) > errsum / std::abs(area)) goto LABEL_115;
    //
    // test on divergence.
    //
LABEL_110:
    if (ksgn == -1 && std::max(std::abs(result), std::abs(area)) <= defabs * 1.0e-2) goto LABEL_130;
    if (1.0e-2 > (result / area) || (result / area) > 1.0e2 || errsum > std::abs(area)) ier = 6;
    goto LABEL_130;
    // 
    // compute global integral sum.
    //
LABEL_115:
    result = 0.0;
    for (k = 1; k <= last; ++k) {
        result += ARRAYF(rlist, k);
    }
    abserr = errsum;
LABEL_130:
    neval = 30 * last - 15;
    if (inf == 2)  neval = 2 * neval;
    if (ier > 2) ier--;
LABEL_999:
    return;
}

/**
 * @fn void dgagp
*/
void dqagp(QUADPACK_CPP_FUNCTION f, const double a, const double b, const int npts2, const double *points,
           const double epsabs, const double epsrel, double &result, double &abserr, int &neval, int &ier,
           const int leniw, const int lenw, int &last, int *iwork, double *work, void *user_data)
{
    int limit, lvl, l1, l2, l3, l4;
    // first executable statement  dqagp
    ier    = 6;
    neval  = 0;
    last   = 0;
    result = 0.0;
    abserr = 0.0;
    if (leniw < (3*npts2 - 2) || lenw < (leniw*2 - npts2) || npts2 < 2) goto LABEL_10;
    //
    // prepare call for dqagpe.
    //
    limit = (leniw - npts2) / 2;
    l1    = limit + 1;
    l2    = limit + l1;
    l3    = limit + l2;
    l4    = limit + l3;
    //
    dqagpe(f , a, b, npts2, points, epsabs, epsrel, limit, result, abserr, 
        neval, ier, &ARRAYF(work, 1), &ARRAYF(work, l1), &ARRAYF(work, l2), &ARRAYF(work, l3) ,&ARRAYF(work, l4),
        &ARRAYF(iwork, 1), &ARRAYF(iwork, l1), &ARRAYF(iwork, l2), last, user_data);
    //
    // call error handler if necessary.
    //
    lvl = 0;
LABEL_10:
    if (ier == 6) lvl = 1;
    if (ier != 0) {
        std::string messg = "abnormal return from dqagp";
        xerror(messg, ier, lvl);   
    }
    return;
}

/**
 * @fn void dgagpe
*/
void dqagpe(QUADPACK_CPP_FUNCTION f, const double a, const double b, const int npts2, const double *points,
            const double epsabs, const double epsrel, const int limit, double &result, double &abserr, int &neval,
            int &ier, double *alist, double *blist, double *rlist, double *elist, double *pts, int *iord,
            int *level, int *ndin, int &last, void *user_data)
{
    int i, ip1, id, ierro, ind1, ind2, iroff1, iroff2, iroff3, j, jlow, jupbnd;
    int k, ksgn, ktmin, L, levcur, levmax, maxerr, nint, nintp1, npts, nres, nrmax, numrl2;
    bool extrap, noext;
    double abseps, area, area1, area12, area2, a1, a2, b1, b2, correc;
    double defabs, defab1, defab2, dres, erlarg, erlast, errbnd, errmax, error1;
    double erro12, error2, errsum, ertest, resa, resabs, reseps, sign, temp;
    double res3la[3], rlist2[52];
    double epmach, uflow, oflow;
    // first executable statement  dqagpe
    epmach = ARRAYF(d1mach, 4);
    //
    // test on validity of parameters
    // -----------------------------
    //
    ier              = 0;
    neval            = 0;
    last             = 0;
    result           = 0.0;
    abserr           = 0.0;
    ARRAYF(alist, 1) = a;
    ARRAYF(blist, 1) = b;
    ARRAYF(rlist, 1) = 0.0;
    ARRAYF(elist, 1) = 0.0;
    ARRAYF(iord,  1) = 0;
    ARRAYF(level, 1) = 0;
    npts = npts2 - 2;
    if (npts2 < 2 || limit <= npts || (epsabs <= 0.0 && epsrel < std::max(50.0 * epmach, 5.0e-29))) ier = 6;
    if (ier == 6) goto LABEL_999;
    //
    // if any break points are provided, sort them into an
    // ascending sequence.
    //
    sign = 1.0;
    if (a > b) sign = -1.0;
    ARRAYF(pts, 1) = std::min(a, b);
    if (npts == 0) goto LABEL_15;
    for (i = 1; i <= npts; ++i) {
        ARRAYF(pts, i+1) = ARRAYF(points, i);
    }
LABEL_15:
    ARRAYF(pts, npts + 2) = std::max(a, b);
    nint = npts + 1;
    a1 = ARRAYF(pts, 1);
    if (npts == 0) goto LABEL_40;
    nintp1 = nint + 1;
    for (i = 1; i <= nint; ++i) {
        ip1 = i + 1;
        for (j = ip1; j <= nintp1; ++j) {
          if (ARRAYF(pts, i) <= ARRAYF(pts, j)) continue;
          temp           = ARRAYF(pts, i);
          ARRAYF(pts, i) = ARRAYF(pts, j);
          ARRAYF(pts, j) = temp;
        }
    }
    if (ARRAYF(pts, 1) != std::min(a, b) || ARRAYF(pts, nintp1) != std::max(a, b)) ier = 6;
    if (ier == 6) goto LABEL_999;
    //
    // compute first integral and error approximations.
    // ------------------------------------------------
    //
LABEL_40:
    resabs = 0.0;
    for (i = 1; i <= nint; ++i) {
        b1 = ARRAYF(pts, i+1);
        dqk21(f, a1, b1, area1, error1, defabs, resa, user_data);
        abserr          += error1;
        result          += area1;
        ARRAYF(ndin, i)  = 0;
        if(error1 == resa && error1 != 0.0) ARRAYF(ndin, i) = 1;
        resabs          += defabs;
        ARRAYF(level, i) = 0;
        ARRAYF(elist, i) = error1;
        ARRAYF(alist, i) = a1;
        ARRAYF(blist, i) = b1;
        ARRAYF(rlist, i) = area1;
        ARRAYF(iord,  i) = i;
        a1               = b1;
    }
    errsum = 0.0;
    for (i = 1; i <= nint; ++i) {
        if (ARRAYF(ndin, i) == 1) ARRAYF(elist, i) = abserr;
        errsum += ARRAYF(elist, i);
    }
    //
    // test on accuracy.
    //
    last   = nint;
    neval  = 21 * nint;
    dres   = std::abs(result);
    errbnd = std::max(epsabs, epsrel * dres);
    if (abserr <= 1.0e2 * epmach * resabs && abserr > errbnd) ier = 2;
    if (nint == 1) goto LABEL_80;
    for (i = 1; i <= npts; ++i) {
        jlow = i + 1;
        ind1 = ARRAYF(iord, i);
        for (j = jlow; j <= nint; ++j) {
            ind2 = ARRAYF(iord, j);
            if (ARRAYF(elist, ind1) > ARRAYF(elist, ind2)) continue;
            ind1 = ind2;
            k    = j;
        }
        if (ind1 == ARRAYF(iord, i)) continue;
        ARRAYF(iord, k) = ARRAYF(iord, i);
        ARRAYF(iord, i) = ind1;
    }
    if (limit < npts2) ier = 1;
LABEL_80:
    if (ier != 0 || abserr <= errbnd) goto LABEL_210;
    //
    // initialization
    // --------------
    //
    ARRAYF(rlist2, 1) = result;
    maxerr            = ARRAYF(iord, 1);
    errmax            = ARRAYF(elist, maxerr);
    area              = result;
    nrmax             = 1;
    nres              = 0;
    numrl2            = 1;
    ktmin             = 0;
    extrap            = false;
    noext             = false;
    erlarg            = errsum;
    ertest            = errbnd;
    levmax            = 1;
    iroff1            = 0;
    iroff2            = 0;
    iroff3            = 0;
    ierro             = 0;
    uflow             = ARRAYF(d1mach, 1);
    oflow             = ARRAYF(d1mach, 2);
    abserr            = oflow;
    ksgn              = -1;
    if (dres >= (1.0 - 50.0 * epmach) * resabs) ksgn = 1;
    //
    // main for-loop
    // ------------
    //
    for (last = npts2; last <= limit; ++last) {
        //
        // bisect the subinterval with the nrmax-th largest error
        // estimate.
        //
        levcur = ARRAYF(level, maxerr) + 1;
        a1     = ARRAYF(alist, maxerr);
        b1     = 0.5 * (ARRAYF(alist, maxerr) + ARRAYF(blist, maxerr));
        a2     = b1;
        b2     = ARRAYF(blist, maxerr);
        erlast = errmax;
        dqk21(f, a1, b1, area1, error1, resa, defab1, user_data);
        dqk21(f, a2, b2, area2, error2, resa, defab2, user_data);
        //
        // improve previous approximations to integral
        // and error and test for accuracy.
        //
        neval  = neval + 42;
        area12 = area1 + area2;
        erro12 = error1 + error2;
        errsum = errsum + erro12 - errmax;
        area   = area + area12 - ARRAYF(rlist, maxerr);
        if (defab1 == error1 || defab2 == error2) goto LABEL_95;
        if (std::abs(ARRAYF(rlist, maxerr) - area12) > 1.0e-5 * std::abs(area12) || erro12 < 0.99 * errmax) goto LABEL_90;
        if (extrap)  iroff2++;
        if (!extrap) iroff1++;
LABEL_90:
        if (last > 10 && erro12 > errmax) iroff3++;
LABEL_95:
        ARRAYF(level, maxerr) = levcur;
        ARRAYF(level, last)   = levcur;
        ARRAYF(rlist, maxerr) = area1;
        ARRAYF(rlist, last)   = area2;
        errbnd = std::max(epsabs, epsrel * std::abs(area));
        //
        // test for roundoff error and eventually set error flag.
        //
        if (iroff1 + iroff2 >= 10 || iroff3 >= 20) ier = 2;
        if (iroff2 >= 5) ierro = 3;
        //
        // set error flag in the case that the number of
        // subintervals equals limit.
        //
        if (last == limit) ier = 1;
        //
        // set error flag in the case of bad integrand behaviour
        // at a point of the integration range
        //
        if (std::max(std::abs(a1), std::abs(b2)) <= (1.0 + 1.0e2 * epmach) * (std::abs(a2) + 1.0e3 * uflow)) ier = 4;
        //
        // append the newly-created intervals to the list.
        //
        if(error2 > error1) goto LABEL_100;
        ARRAYF(alist, last)   = a2;
        ARRAYF(blist, maxerr) = b1;
        ARRAYF(blist, last)   = b2;
        ARRAYF(elist, maxerr) = error1;
        ARRAYF(elist, last)   = error2;
        goto LABEL_110;
LABEL_100:
        ARRAYF(alist, maxerr) = a2;
        ARRAYF(alist, last)   = a1;
        ARRAYF(blist, last)   = b1;
        ARRAYF(rlist, maxerr) = area2;
        ARRAYF(rlist, last)   = area1;
        ARRAYF(elist, maxerr) = error2;
        ARRAYF(elist, last)   = error1;
        //
        // call subroutine dqpsrt to maintain the descending ordering
        // in the list of error estimates and select the subinterval
        // with nrmax-th largest error estimate (to be bisected next).
        //
LABEL_110:
        dqpsrt(limit, last, maxerr, errmax, elist, iord, nrmax);
        // jump out of for-loop
        if (errsum <= errbnd) goto LABEL_190;
        // jump out of for-loop
        if (ier != 0) goto LABEL_170;
        if (noext) goto LABEL_160;
        erlarg -= erlast;
        if (levcur + 1 <= levmax) erlarg += erro12;
        if (extrap) goto LABEL_120;
        //
        // test whether the interval to be bisected next is the
        // smallest interval.
        //
        if (ARRAYF(level, maxerr) + 1 <= levmax) goto LABEL_160;
        extrap = true;
        nrmax  = 2;
LABEL_120:
        if (ierro == 3 || erlarg <= ertest) goto LABEL_140;
        //
        // the smallest interval has the largest error.
        // before bisecting decrease the sum of the errors over
        // the larger intervals (erlarg) and perform extrapolation.
        //
        id     = nrmax;
        jupbnd = last;
        if (last > (2 + limit / 2)) jupbnd = limit + 3 - last;
        for (k = id; k <= jupbnd; ++k) {
            maxerr = ARRAYF(iord, nrmax);
            errmax = ARRAYF(elist, maxerr);
            // jump out of do-loop
            if (ARRAYF(level, maxerr) + 1 <= levmax) goto LABEL_160;
            nrmax = nrmax+1;
        }
        //
        // perform extrapolation.
        //
LABEL_140:
        numrl2++;
        ARRAYF(rlist2, numrl2) = area;
        if (numrl2 <= 2) goto LABEL_155;
        dqelg(numrl2, rlist2, reseps, abseps, res3la, nres);
        ktmin = ktmin + 1;
        if (ktmin > 5 && abserr < 1.0e-3 * errsum) ier = 5;
        if (abseps >= abserr) goto LABEL_150;
        ktmin  = 0;
        abserr = abseps;
        result = reseps;
        correc = erlarg;
        ertest = std::max(epsabs, epsrel * std::abs(reseps));
        // jump out of do-loop
        if (abserr < ertest) goto LABEL_170;
        //
        // prepare bisection of the smallest interval.
        //
LABEL_150:
        if (numrl2 == 1) noext = true;
        if(ier >= 5) goto LABEL_170;
LABEL_155:
        maxerr = ARRAYF(iord, 1);
        errmax = ARRAYF(elist, maxerr);
        nrmax  = 1;
        extrap = false;
        levmax++;
        erlarg = errsum;
LABEL_160:
        continue;
    }
    //
    // set the final result.
    // ---------------------
    //
    //
LABEL_170:
    if (abserr == oflow) goto LABEL_190;
    if (ier + ierro == 0) goto LABEL_180;
    if (ierro == 3) abserr = abserr + correc;
    if (ier == 0) ier = 3;
    if (result != 0.0 && area != 0.0) goto LABEL_175;
    if (abserr > errsum) goto LABEL_190;
    if (area == 0.0) goto LABEL_210;
    goto LABEL_180;
LABEL_175:
    if (abserr / std::abs(result) > errsum / std::abs(area)) goto LABEL_190;
    //
    // test on divergence.
    //
LABEL_180:
    if(ksgn == -1 && std::max(std::abs(result), std::abs(area)) <= resabs * 1.0e-2) goto LABEL_210;
    if (1.0e-2 > result/area || result/area > 1.0e2 || errsum > std::abs(area)) ier = 6;
    goto LABEL_210;
    //
    // compute global integral sum.
    //
LABEL_190:
    result = 0.0;
    for (k = 1; k <= last; ++k) {
        result += ARRAYF(rlist, k);
    }
    abserr = errsum;
LABEL_210:
    if (ier > 2) ier -= 1;
    result *= sign;
LABEL_999:
    return;
}

/**
 * @fn void dgags
*/
void dqags(QUADPACK_CPP_FUNCTION f, const double a, const double b, const double epsabs, const double epsrel, 
           double &result, double &abserr, int &neval, int &ier, const int limit, const int lenw, int &last, 
           int *iwork, double *work, void *user_data)
{
    int lvl, l1, l2, l3;
    //
    // check validity of limit and lenw.
    //
    // first executable statement  dqags
    ier    = 6;
    neval  = 0;
    last   = 0;
    result = 0.0;
    abserr = 0.0;
    if (limit < 1 || lenw < 4*limit) goto LABEL_10;
    //
    // prepare call for dqagse.
    //
    l1 = limit + 1;
    l2 = limit + l1;
    l3 = limit + l2;
    //
    dqagse(f, a, b, epsabs, epsrel, limit, result, abserr, neval, ier,
        &ARRAYF(work, 1), &ARRAYF(work, l1), &ARRAYF(work, l2), &ARRAYF(work, l3), iwork, last, user_data);
    //
    // call error handler if necessary.
    //
    lvl = 0;
LABEL_10:
    if (ier == 6) lvl = 1;
    if (ier != 0) {
        std::string messg = "abnormal return from dqags";
        xerror(messg, ier, lvl);
    }
    return;
}

/**
 * @fn void dgagse
*/
void dqagse(QUADPACK_CPP_FUNCTION f, const double a, const double b, const double epsabs, const double epsrel, 
            const int limit, double &result, double &abserr, int &neval, int &ier, double *alist, double *blist,
            double *rlist, double *elist, int *iord, int &last, void *user_data)
{
    int id, ierro, iroff1, iroff2, iroff3, k, ksgn, ktmin, maxerr, nres, nrmax, numrl2;
    bool noext, extrap;
    double a1, a2, area, area1, area12, area2, b1, b2, defab1, defab2, defabs;
    double dres, erlarg, erlast, errbnd, errmax, error1, erro12, error2, errsum;
    double ertest, small, jupbnd, correc, resabs, reseps, abseps;
    double res3la[3], rlist2[52];
    double uflow, oflow, epmach;
    // first executable statement  dqagse
    epmach = ARRAYF(d1mach, 4);
    //
    // test on validity of parameters
    // ------------------------------
    //
    ier              = 0;
    neval            = 0;
    last             = 0;
    result           = 0.0;
    abserr           = 0.0;
    ARRAYF(alist, 1) = a;
    ARRAYF(blist, 1) = b;
    ARRAYF(rlist, 1) = 0.0;
    ARRAYF(elist, 1) = 0.0;
    if (epsabs <= 0.0 && epsrel < std::max(50.0 * epmach, 5.0e-29)) ier = 6;
    if (ier == 6) goto LABEL_999;
    // 
    // first approximation to the integral
    // -----------------------------------
    //
    uflow = ARRAYF(d1mach, 1);
    oflow = ARRAYF(d1mach, 2);
    ierro = 0;
    dqk21(f, a, b, result, abserr, defabs, resabs, user_data);
    // 
    // test on accuracy
    // 
    dres             = std::abs(result);
    errbnd           = std::max(epsabs, epsrel * dres);
    last             = 1;
    ARRAYF(rlist, 1) = result;
    ARRAYF(elist, 1) = abserr;
    ARRAYF(iord,  1) = 1;
    if (abserr <= 1.0e2 * epmach * defabs && abserr > errbnd) ier = 2;
    if (limit == 1) ier = 1;
    if (ier != 0 || (abserr <= errbnd && abserr != resabs) || abserr == 0.0) goto LABEL_140;
    //
    // initialization
    // --------------
    //
    ARRAYF(rlist2, 1) = result;
    errmax = abserr;
    maxerr = 1;
    area   = result;
    errsum = abserr;
    abserr = oflow;
    nrmax  = 1;
    nres   = 0;
    numrl2 = 2;
    ktmin  = 0;
    extrap = false;
    noext  = false;
    iroff1 = 0;
    iroff2 = 0;
    iroff3 = 0;
    ksgn   = -1;
    if (dres >= (1.0 - 50.0 * epmach) * defabs) ksgn = 1;
    //
    // main for-loop
    // ---------
    //
    for (last = 2; last <= limit; ++last) {
        //
        // bisect the subinterval with the nrmax-th largest error 
        // estimate.
        //
        a1 = ARRAYF(alist, maxerr);
        b1 = 0.5 * (ARRAYF(alist, maxerr) + ARRAYF(blist, maxerr));
        a2 = b1;
        b2 = ARRAYF(blist, maxerr);
        erlast = errmax;
        dqk21(f, a1, b1, area1, error1, resabs, defab1, user_data);
        dqk21(f, a2, b2, area2, error2, resabs, defab2, user_data);
        // 
        // improve previous approximations to integral 
        // and error and test for accuracy.
        //
        area12 = area1  + area2;
        erro12 = error1 + error2;
        errsum = errsum + erro12 - errmax;
        area   = area   + area12 - ARRAYF(rlist, maxerr);
        if (defab1 == error1 || defab2 == error2) goto LABEL_15;
        if (std::abs(ARRAYF(rlist, maxerr) - area12) > 1.0e-5 * std::abs(area12) || erro12 < 0.99 * errmax) goto LABEL_10;
        if (extrap)  iroff2++;
        if (!extrap) iroff1++;
LABEL_10:
        if (last > 10 && erro12 > errmax) iroff3++;
LABEL_15:
        ARRAYF(rlist, maxerr) = area1;
        ARRAYF(rlist, last)   = area2;
        errbnd                = std::max(epsabs, epsrel * std::abs(area));
        // 
        // test for roundoff error and eventually set error flag.
        //
        if (iroff1 + iroff2 >= 10 || iroff3 >= 20) ier = 2;
        if (iroff2 >= 5) ierro = 3;
        //
        // set error flag in the case that the number of subintervals 
        // equals limit.
        //
        if (last == limit) ier = 1;
        //
        // set error flag in the case of bad integrand behaviour 
        // at a point of the integration range.
        //
        if (std::max(std::abs(a1), std::abs(b2)) <= (1.0 + 1.0e2 * epmach) * (std::abs(a2) + 1.0e3 * uflow)) ier = 4;
        //
        // append the newly-created intervals to the list.
        //
        if (error2 > error1) goto LABEL_20;
        ARRAYF(alist, last)   = a2;
        ARRAYF(blist, maxerr) = b1;
        ARRAYF(blist, last)   = b2;
        ARRAYF(elist, maxerr) = error1;
        ARRAYF(elist, last)   = error2;
        goto LABEL_30;
LABEL_20:
        ARRAYF(alist, maxerr) = a2;
        ARRAYF(alist, last)   = a1;
        ARRAYF(blist, last)   = b1;
        ARRAYF(rlist, maxerr) = area2;
        ARRAYF(rlist, last)   = area1;
        ARRAYF(elist, maxerr) = error2;
        ARRAYF(elist, last)   = error1;
        //
        // call subroutine dqpsrt to maintain the descending ordering
        // in the list of error estimates and select the subinterval
        // with nrmax-th largest error estimate (to be bisected next).
        //
LABEL_30:
        dqpsrt(limit, last, maxerr, errmax, elist, iord, nrmax);
        // jump out of for-loop
        if (errsum <= errbnd) goto LABEL_115;
        // jump out of for-loop
        if (ier != 0)  goto LABEL_100;
        if (last == 2) goto LABEL_80;
        if (noext) goto LABEL_90;
        erlarg -= erlast;
        if (std::abs(b1-a1) > small) erlarg += erro12;
        if (extrap) goto LABEL_40;
        // 
        // test whether the interval to be bisected next is the
        // smallest interval.
        // 
        if(std::abs(ARRAYF(blist, maxerr) - ARRAYF(alist, maxerr)) > small) goto LABEL_90;
        extrap = true;
        nrmax  = 2;
LABEL_40:
        if(ierro == 3 || erlarg <= ertest) goto LABEL_60;
        //
        // the smallest interval has the largest error.
        // before bisecting decrease the sum of the errors over the
        // larger intervals (erlarg) and perform extrapolation.
        // 
        id     = nrmax;
        jupbnd = last;
        if (last > (2 + limit / 2)) jupbnd = limit + 3 - last;
        for (k = id; k <= jupbnd; ++k) {
            maxerr = ARRAYF(iord, nrmax);
            errmax = ARRAYF(elist, maxerr);
            // jump out of for-loop
            if (std::abs(ARRAYF(blist, maxerr) - ARRAYF(alist, maxerr)) > small) goto LABEL_90;
            nrmax += 1;
        }
        //
        // perform extrapolation.
        // 
LABEL_60:
        numrl2++;
        ARRAYF(rlist2, numrl2) = area;
        dqelg(numrl2, rlist2, reseps, abseps, res3la, nres);
        ktmin++;
        if(ktmin > 5 && abserr < 1.0e-3 * errsum) ier = 5;
        if(abseps >= abserr) goto LABEL_70;
        ktmin  = 0;
        abserr = abseps;
        result = reseps;
        correc = erlarg;
        ertest = std::max(epsabs, epsrel * std::abs(reseps));
        // jump out of for-loop
        if (abserr <= ertest) goto LABEL_100;
        //
        // prepare bisection of the smallest interval.
        // 
LABEL_70:
        if (numrl2 == 1) noext = true;
        if (ier == 5) goto LABEL_100;
        maxerr = ARRAYF(iord, 1);
        errmax = ARRAYF(elist, maxerr);
        nrmax  = 1;
        extrap = false;
        small  = small*0.5;
        erlarg = errsum;
        goto LABEL_90;
LABEL_80:
        small  = std::abs(b - a) * 0.375;
        erlarg = errsum;
        ertest = errbnd;
        ARRAYF(rlist2, 2) = area;
LABEL_90:
        continue;
    }
    //
    // set final result and error estimate.
    //
LABEL_100:
    if (abserr == oflow) goto LABEL_115;
    if (ier + ierro == 0) goto LABEL_110;
    if (ierro == 3) abserr += correc;
    if (ier == 0) ier = 3;
    if (result != 0.0 && area != 0.0) goto LABEL_105;
    if (abserr > errsum) goto LABEL_115;
    if (area == 0.0) goto LABEL_130;
    goto LABEL_110;
LABEL_105:
    if (abserr / std::abs(result) > errsum / std::abs(area)) goto LABEL_115;
    //
    // test on divergence.
    //
LABEL_110:
    if(ksgn == -1 && std::max(std::abs(result), std::abs(area)) <= defabs * 1.0e-2) goto LABEL_130;
    if(1.0e-2 > (result/area) || (result/area) > 1.0e2 || errsum > std::abs(area)) ier = 6;
    goto LABEL_130;
    // 
    // compute global integral sum.
    //
LABEL_115:
    result = 0.0;
    for (k = 1; k <= last; ++k) {
        result += ARRAYF(rlist, k);
    }
    abserr = errsum;
LABEL_130:
    if(ier > 2) ier -= 1;
LABEL_140:
    neval = 42 * last - 21;
LABEL_999:
    return;
}

void dqawc(QUADPACK_CPP_FUNCTION f, const double a, const double b, const double c, const double epsabs, const double epsrel,
           double &result, double &abserr, int &neval, int &ier, const int limit, const int lenw, int &last, int *iwork,
           double *work, void *user_data)
{
    int lvl, l1, l2, l3;
    // first executable statement  dqawc
    ier    = 6;
    neval  = 0;
    last   = 0;
    result = 0.0;
    abserr = 0.0;
    if(limit < 1 || lenw < limit*4) goto LABEL_10;
    //
    // prepare call for dqawce.
    //
    l1 = limit + 1;
    l2 = limit + l1;
    l3 = limit + l2;
    dqawce(f, a, b, c, epsabs, epsrel, limit, result, abserr, neval, ier,
        &ARRAYF(work,1 ), &ARRAYF(work,l1), &ARRAYF(work, l2), &ARRAYF(work, l3), iwork, last, user_data);
    //
    // call error handler if necessary.
    //
    lvl = 0;
LABEL_10:
    if (ier == 6) ier = 1;
    if (ier != 0) {
        std::string messg = "abnormal return from dqawc";
        xerror(messg, ier, lvl);
    }
    return;
}

void dqawce(QUADPACK_CPP_FUNCTION f, const double a, const double b, const double c, const double epsabs, const double epsrel,
            const int limit, double &result, double &abserr, int &neval, int &ier, double *alist, double *blist, double *rlist,
            double *elist, int *iord, int &last, void *user_data)
{
    int iroff1, iroff2, k, krule, L, maxerr, nev, nrmax;
    double a1, a2, aa, area, area1, area12, area2, b1, b2, bb, errbnd, errmax;
    double error1, erro12, error2, errsum;
    double epmach, uflow;
    // first executable statement  dqawce
    epmach = ARRAYF(d1mach, 4);
    uflow  = ARRAYF(d1mach, 1);
    //
    //
    // test on validity of parameters
    // ------------------------------
    //
    ier              = 6;
    neval            = 0;
    last             = 0;
    ARRAYF(alist, 1) = a;
    ARRAYF(blist, 1) = b;
    ARRAYF(rlist, 1) = 0.0;
    ARRAYF(elist, 1) = 0.0;
    ARRAYF(iord,  1) = 0;
    result           = 0.0;
    abserr          = 0.0;
    if (c == a || c == b || (epsabs <= 0.0 && epsrel < std::max(50.0 * epmach, 5.0e-29))) goto LABEL_999;
    //
    // first approximation to the integral
    // -----------------------------------
    //
    aa    = a;
    bb    = b;
    if (a <= b) goto LABEL_10;
    aa    = b;
    bb    = a;
    LABEL_10:
    ier   = 0;
    krule = 1;
    dqc25c(f, aa, bb, c, result, abserr, krule, neval, user_data);
    last  = 1;
    ARRAYF(rlist, 1) = result;
    ARRAYF(elist, 1) = abserr;
    ARRAYF(iord,  1) = 1;
    ARRAYF(alist, 1) = a;
    ARRAYF(blist, 1) = b;
    //
    // test on accuracy
    //
    errbnd = std::max(epsabs, epsrel * std::abs(result));
    if(limit == 1) ier = 1;
    if(abserr < std::min(1.0e-2 * std::abs(result), errbnd) || ier == 1) goto LABEL_70;
    //
    // initialization
    // --------------
    //
    ARRAYF(alist, 1) = aa;
    ARRAYF(blist, 1) = bb;
    ARRAYF(rlist, 1) = result;
    errmax = abserr;
    maxerr = 1;
    area   = result;
    errsum = abserr;
    nrmax  = 1;
    iroff1 = 0;
    iroff2 = 0;
    //
    // main for-loop
    // ------------
    //
    for (last = 2; last <= limit; ++last) {
        //
        // bisect the subinterval with nrmax-th largest
        // error estimate.
        //
        a1 = ARRAYF(alist, maxerr);
        b1 = 0.5 * (ARRAYF(alist, maxerr) + ARRAYF(blist, maxerr));
        b2 = ARRAYF(blist, maxerr);
        if (c <= b1 && c > a1) b1 = 0.5 * (c + b2);
        if (c > b1  && c < b2) b1 = 0.5 * (a1 + c);
        a2 = b1;
        krule = 2;
        dqc25c(f, a1, b1, c, area1, error1, krule, nev, user_data);
        neval += nev;
        dqc25c(f, a2, b2, c, area2, error2, krule, nev, user_data);
        neval += nev;
        //
        // improve previous approximations to integral
        // and error and test for accuracy.
        //
        area12 = area1  + area2;
        erro12 = error1 + error2;
        errsum = errsum + erro12 - errmax;
        area   = area   + area12 - ARRAYF(rlist, maxerr);
        if (std::abs(ARRAYF(rlist, maxerr) - area12) < 1.0e-5 * std::abs(area12) 
           && erro12 >= 0.99 * errmax && krule == 0) iroff1++;
        if (last > 10 && erro12 > errmax && krule == 0) iroff2++;
        ARRAYF(rlist, maxerr) = area1;
        ARRAYF(rlist,   last) = area2;
        errbnd = std::max(epsabs, epsrel * std::abs(area));
        if (errsum <= errbnd) goto LABEL_15;
        //
        // test for roundoff error and eventually set error flag.
        //
        if (iroff1 >= 6 && iroff2 > 20) ier = 2;
        //
        // set error flag in the case that number of interval
        // bisections exceeds limit.
        //
        if (last == limit) ier = 1;
        //
        // set error flag in the case of bad integrand behaviour
        // at a point of the integration range.
        //
        if (std::max(std::abs(a1), std::abs(b2)) <= (1.0 + 1.0e2 * epmach) * (std::abs(a2) + 1.0e3 * uflow)) ier = 3;
        //
        // append the newly-created intervals to the list.
        //
LABEL_15:
        if(error2 > error1) goto LABEL_20;
        ARRAYF(alist, last)   = a2;
        ARRAYF(blist, maxerr) = b1;
        ARRAYF(blist, last)   = b2;
        ARRAYF(elist, maxerr) = error1;
        ARRAYF(elist, last)   = error2;
        goto LABEL_30;
LABEL_20:
        ARRAYF(alist, maxerr) = a2;
        ARRAYF(alist, last)   = a1;
        ARRAYF(blist, last)   = b1;
        ARRAYF(rlist, maxerr) = area2;
        ARRAYF(rlist, last)   = area1;
        ARRAYF(elist, maxerr) = error2;
        ARRAYF(elist, last)   = error1;
        //
        // call subroutine dqpsrt to maintain the descending ordering
        // in the list of error estimates and select the subinterval
        // with nrmax-th largest error estimate (to be bisected next).
        //
LABEL_30:
        dqpsrt(limit, last, maxerr, errmax, elist, iord, nrmax);
        // jump out of for-loop
        if(ier != 0 || errsum <= errbnd) goto LABEL_50;
    }
    //
    // compute final result.
    // ---------------------
    //
LABEL_50:
    result = 0.0;
    for (k = 1; k <= last; ++k) {
        result += ARRAYF(rlist, k);
    }
    abserr = errsum;
LABEL_70:
    if (aa == b) result = -result;
LABEL_999:
    return;
}

void dqawf(QUADPACK_CPP_FUNCTION f, const double a, const double omega, const int integr, const double epsabs, double &result,
           double &abserr, int &neval, int &ier, const int limlst, int &lst, const int leniw, const int maxp1, const int lenw,
           int *iwork, double *work, void *user_data)
{
    int last, limit, ll2, lvl, l1, l2, l3, l4, l5, l6;
    // first executable statement  dqawf
    ier    = 6;
    neval  = 0;
    last   = 0;
    result = 0.0;
    abserr = 0.0;
    if (limlst < 3 || leniw < limlst + 2 || maxp1 < 1 || lenw < (leniw*2 + maxp1*25)) goto LABEL_10;
    //
    // prepare call for dqawfe
    //
    limit = (leniw - limlst) / 2;
    l1    = limlst + 1;
    l2    = limlst + l1;
    l3    = limit  + l2;
    l4    = limit  + l3;
    l5    = limit  + l4;
    l6    = limit  + l5;
    ll2   = limit  + l1;
    dqawfe(f, a, omega, integr, epsabs, limlst, limit, maxp1, result, abserr, neval, ier,
        &ARRAYF(work, 1), &ARRAYF(work, l1), &ARRAYF(iwork, 1), lst, &ARRAYF(work, l2),
        &ARRAYF(work,l3), &ARRAYF(work, l4) ,&ARRAYF(work, l5), &ARRAYF(iwork, l1), &ARRAYF(iwork, ll2), &ARRAYF(work, l6),
        user_data);
    //
    // call error handler if necessary
    //
    lvl = 0;
LABEL_10:
    if (ier == 6) ier = 1;
    if (ier != 0) {
        std::string messg = "abnormal return from dqawf";
        xerror(messg, ier, lvl);
    }
    return;
}

void dqawfe(QUADPACK_CPP_FUNCTION f, const double a, const double omega, const int integr, const double epsabs, const int limlst,
            const int limit, const int maxp1, double &result, double &abserr, int &neval, int &ier, double *rslst, double *erlst,
            int *ierlst, int &lst, double *alist, double *blist, double *rlist, double *elist, int *iord, int *nnlog, double *chebmo, 
            void *user_data)
{
    int ktmin, l, last, ll, iter, momcom, nev, nres, numrl2;
    double abseps, correc, cycle, c1, c2, dl, drl, ep, eps, epsa;
    double errsum, fact, p1, reseps;
    double psum[52], res3la[3];
    const double p = 0.9;
    const double pi = M_PI;
    double uflow;
    //
    // test on validity of parameters
    // ------------------------------
    //
    // first executable statement  dqawfe
    result = 0.0;
    abserr = 0.0;
    neval  = 0;
    lst    = 0;
    ier    = 0;
    if ((integr != 1 && integr != 2) || epsabs <= 0.0 || limlst < 3) ier = 6;
    if (ier == 6) goto LABEL_999;
    if (omega != 0.0) goto LABEL_10;
    //
    // integration by dqagie if omega is zero
    // --------------------------------------
    //
    if(integr == 1) dqagie(f, 0.0, 1, epsabs, 0.0,limit, result, abserr, neval, ier,
                        alist, blist, rlist, elist, iord, last, user_data);
    ARRAYF(rslst,  1) = result;
    ARRAYF(erlst,  1) = abserr;
    ARRAYF(ierlst, 1) = ier;
    lst               = 1;
    goto LABEL_999;
    //
    // initializations
    // ---------------
    //
LABEL_10:
    l      = static_cast<int>(std::abs(omega));
    dl     = 2 * static_cast<double>(l) + 1;
    cycle  = dl * pi / std::abs(omega);
    ier    = 0;
    ktmin  = 0;
    neval  = 0;
    numrl2 = 0;
    nres   = 0;
    c1     = a;
    c2     = cycle + a;
    p1     = 1.0 - p;
    uflow  = ARRAYF(d1mach, 1);
    eps    = epsabs;
    if(epsabs > uflow/p1) eps = epsabs * p1;
    ep     = eps;
    fact   = 1.0;
    correc = 0.0;
    abserr = 0.0;
    errsum = 0.0;
    ll     = 0;
    //
    // main for-loop
    // -------------
    //
    for (lst = 1; lst <= limlst; ++lst) {
        //
        // integrate over current subinterval.
        //
        epsa = eps * fact;
        dqawoe(f, c1, c2, omega, integr, epsa, 0.0, limit, lst, maxp1,
            ARRAYF(rslst, lst), ARRAYF(erlst, lst), nev, ARRAYF(ierlst, lst), 
            last, alist, blist, rlist, elist, iord, nnlog, momcom, chebmo, user_data);
        neval  += nev;
        fact   *= p;
        errsum += ARRAYF(erlst, lst);
        drl     = 50.0 * std::abs(ARRAYF(rslst, lst));
        //
        // test on accuracy with partial sum
        //
        if (errsum + drl <= epsabs && lst >= 6) goto LABEL_80;
        correc = std::max(correc, ARRAYF(erlst, lst));
        if (ARRAYF(ierlst, lst) != 0) eps = std::max(ep, correc * p1);
        if (ARRAYF(ierlst, lst) != 0) ier = 7;
        if (ier == 7 && errsum + drl <= 10.0 * correc && lst > 5) goto LABEL_80;
        numrl2++;
        if (lst > 1) goto LABEL_20;
        ARRAYF(psum, 1) = ARRAYF(rslst, 1);
        goto LABEL_40;
LABEL_20:
        ARRAYF(psum, numrl2) = ARRAYF(psum, ll) + ARRAYF(rslst, lst);
        if(lst == 2) goto LABEL_40;
        //
        // test on maximum number of subintervals
        //
        if (lst == limlst) ier = 1;
        //
        // perform new extrapolation
        //
        dqelg(numrl2, psum, reseps, abseps, res3la, nres);
        //
        // test whether extrapolated result is influenced by roundoff
        //
        ktmin++;
        if (ktmin >= 15 && abserr <= 1.0e-3 * (errsum + drl)) ier = 4;
        if (abseps > abserr && lst != 3) goto LABEL_30;
        abserr = abseps;
        result = reseps;
        ktmin  = 0;
        //
        // if ier is not 0, check whether direct result (partial sum)
        // or extrapolated result yields the best integral
        // approximation
        // 
        if ((abserr + 10.0 * correc) <= epsabs || (abserr <= epsabs && 10.0*correc >= epsabs)) goto LABEL_60;
LABEL_30:
        if (ier != 0 && ier != 7) goto LABEL_60;
LABEL_40:
        ll  = numrl2;
        c1  = c2;
        c2 += cycle;
    }
    //
    // set final result and error estimate
    // -----------------------------------
    //
LABEL_60:
    abserr += 10.0 * correc;
    if (ier == 0) goto LABEL_999;
    if (result != 0.0 && ARRAYF(psum, numrl2) != 0.0) goto LABEL_70;
    if (abserr > errsum) goto LABEL_80;
    if (ARRAYF(psum, numrl2) == 0.0) goto LABEL_999;
LABEL_70:
    if (abserr / std::abs(result) > (errsum + drl) / std::abs(ARRAYF(psum, numrl2))) goto LABEL_80;
    if (ier >= 1 && ier != 7) abserr += drl;
    goto LABEL_999;
LABEL_80:
    result = ARRAYF(psum, numrl2);
    abserr = errsum + drl;
LABEL_999:
    return;
}


void dqawo(QUADPACK_CPP_FUNCTION f, const double a, const double b, const double omega, const int integr, const double epsabs, 
           const double epsrel, double &result, double &abserr, int &neval, int &ier, const int leniw, const int maxp1, const int lenw,
           int &last, int *iwork, double *work, void *user_data)
{
    int limit, lvl, l1, l2, l3, l4, momcom;
    // first executable statement  dqawo
    ier    = 6;
    neval  = 0;
    last   = 0;
    result = 0.0;
    abserr = 0.0;
    if (leniw < 2 || maxp1 < 1 || lenw < (leniw*2 + maxp1*25)) goto LABEL_10;
    //
    // prepare call for dqawoe
    //
    limit = leniw / 2;
    l1    = limit + 1;
    l2    = limit + l1;
    l3    = limit + l2;
    l4    = limit + l3;
    dqawoe(f, a, b, omega, integr, epsabs, epsrel, limit, 1 ,maxp1, result, abserr, neval, ier, last,
        &ARRAYF( work, 1), &ARRAYF( work, l1), &ARRAYF(work, l2), &ARRAYF(work, l3), 
        &ARRAYF(iwork, 1), &ARRAYF(iwork, l1), momcom, &ARRAYF(work, l4), user_data);
    //
    // call error handler if necessary
    //
    lvl = 0;
LABEL_10:
    if (ier == 6) ier = 1;
    if (ier != 0) {
        std::string messg = "abnormal return from dqawo";
        xerror(messg, ier, lvl);
    }
    return;
}

void dqawoe(QUADPACK_CPP_FUNCTION f, const double a, const double b, const double omega, const int integr, const double epsabs,
            const double epsrel, const int limit, const int icall, const int maxp1, double &result, double &abserr, int  &neval,
            int &ier, int &last, double *alist, double *blist, double *rlist, double *elist, int *iord, int *nnlog, int &momcom,
            double *chebmo, void *user_data)
{
    int id, ierro, iroff1, iroff2, iroff3, jupbnd, k, ksgn;
    int ktmin, maxerr, nev, noext, nres, nrmax, nrmom, numrl2;
    bool extrap, extall;
    double abseps, area, area1, area12, area2, a1, a2, b1, b2, correc, defab1;
    double defab2, defabs, domega, dres, erlarg, erlast, errbnd, errmax, error1;
    double erro12, error2, errsum, ertest, resabs, reseps, small, width;
    double rlist2[52], res3la[3];
    double epmach, uflow, oflow;
    // first executable statement  dqawoe
    epmach = ARRAYF(d1mach, 4);
    //
    // test on validity of parameters
    // ------------------------------
    //
    ier              = 0;
    neval            = 0;
    last             = 0;
    result           = 0.0;
    abserr           = 0.0;
    ARRAYF(alist, 1) = a;
    ARRAYF(blist, 1) = b;
    ARRAYF(rlist, 1) = 0.0;
    ARRAYF(elist, 1) = 0.0;
    ARRAYF(iord,  1) = 0;
    ARRAYF(nnlog, 1) = 0;
    if((integr != 1 && integr != 2) || 
       (epsabs <= 0.0 && epsrel < std::max(50.0 * epmach, 5.0e-29)) || 
       icall < 1 || maxp1 < 1) ier = 6;
    if(ier == 6) goto LABEL_999;
    //
    // first approximation to the integral
    // -----------------------------------
    //
    domega = std::abs(omega);
    nrmom  = 0;
    if (icall > 1) goto LABEL_5;
    momcom = 0;
LABEL_5:
    dqc25f(f, a, b, domega, integr, nrmom, maxp1, 0, result, abserr,
        neval, defabs, resabs, momcom, chebmo, user_data);
    //
    // test on accuracy.
    //
    dres             = std::abs(result);
    errbnd           = std::max(epsabs, epsrel * dres);
    ARRAYF(rlist, 1) = result;
    ARRAYF(elist, 1) = abserr;
    ARRAYF(iord,  1) = 1;
    if (abserr <= 1.0e2 * epmach * defabs && abserr > errbnd) ier = 2;
    if (limit == 1) ier = 1;
    if (ier != 0 || abserr <= errbnd) goto LABEL_200;
    //
    // initializations
    // ---------------
    //
    uflow  = ARRAYF(d1mach, 1);
    oflow  = ARRAYF(d1mach, 2);
    errmax = abserr;
    maxerr = 1;
    area   = result;
    errsum = abserr;
    abserr = oflow;
    nrmax  = 1;
    extrap = false;
    noext  = false;
    ierro  = 0;
    iroff1 = 0;
    iroff2 = 0;
    iroff3 = 0;
    ktmin  = 0;
    small  = 0.75 * std::abs(b - a);
    nres   = 0;
    numrl2 = 0;
    extall = false;
    if (0.5 * std::abs(b - a) * domega > 2.0) goto LABEL_10;
    numrl2 = 1;
    extall = true;
    ARRAYF(rlist2, 1) = result;
LABEL_10:
    if (0.25 * std::abs(b - a) * domega <= 2.0) extall = true;
    ksgn = -1;
    if (dres >= (1.0 - 50.0*epmach)*defabs) ksgn = 1;
    //
    // main for-loop
    // -------------
    //
    for (last = 2; last <= limit; ++last) {
        //
        // bisect the subinterval with the nrmax-th largest
        // error estimate.
        //
        nrmom  = ARRAYF(nnlog, maxerr) + 1;
        a1     = ARRAYF(alist, maxerr);
        b1     = 0.5 * (ARRAYF(alist, maxerr) + ARRAYF(blist, maxerr));
        a2     = b1;
        b2     = ARRAYF(blist, maxerr);
        erlast = errmax;
        dqc25f(f, a1, b1, domega, integr, nrmom, maxp1, 0, area1, 
            error1, nev, resabs, defab1, momcom, chebmo, user_data);
        neval += nev;
        dqc25f(f, a2, b2, domega, integr, nrmom, maxp1, 1, area2, 
            error2, nev, resabs, defab2, momcom, chebmo, user_data);
        neval += nev;
        //
        // improve previous approximations to integral
        // and error and test for accuracy.
        //
        area12 = area1  + area2;
        erro12 = error1 + error2;
        errsum = errsum + erro12 - errmax;
        area   = area   + area12 - ARRAYF(rlist, maxerr);
        if (defab1 == error1 || defab2 == error2) goto LABEL_25;
        if (std::abs(ARRAYF(rlist, maxerr) - area12) > 1.0e-5 * std::abs(area12) || erro12 < 0.99*errmax) goto LABEL_20;
        if (extrap)  iroff2++;
        if (!extrap) iroff1++;
LABEL_20:
        if (last > 10 && erro12 > errmax) iroff3++;
LABEL_25:
        ARRAYF(rlist, maxerr) = area1;
        ARRAYF(rlist,   last) = area2;
        ARRAYF(nnlog, maxerr) = nrmom;
        ARRAYF(nnlog,   last) = nrmom;
        errbnd                = std::max(epsabs, epsrel * std::abs(area));
        //
        // test for roundoff error and eventually set error flag.
        //
        if (iroff1 + iroff2 >= 10 || iroff3 >= 20) ier = 2;
        if (iroff2 >= 5) ierro = 3;
        //
        // set error flag in the case that the number of
        // subintervals equals limit.
        //
        if (last == limit) ier = 1;
        //
        // set error flag in the case of bad integrand behaviour
        // at a point of the integration range.
        //
        if(std::max(std::abs(a1), std::abs(b2)) <= (1.0 + 1.0e2 * epmach) * (std::abs(a2) + 1.0e3 * uflow)) ier = 4;
        //
        // append the newly-created intervals to the list.
        //
        if (error2 > error1) goto LABEL_30;
        ARRAYF(alist,   last) = a2;
        ARRAYF(blist, maxerr) = b1;
        ARRAYF(blist,   last) = b2;
        ARRAYF(elist, maxerr) = error1;
        ARRAYF(elist,   last) = error2;
        goto LABEL_40;
LABEL_30:
        ARRAYF(alist, maxerr) = a2;
        ARRAYF(alist,   last) = a1;
        ARRAYF(blist,   last) = b1;
        ARRAYF(rlist, maxerr) = area2;
        ARRAYF(rlist,   last) = area1;
        ARRAYF(elist, maxerr) = error2;
        ARRAYF(elist,   last) = error1;
        //
        // call subroutine dqpsrt to maintain the descending ordering
        // in the list of error estimates and select the subinterval
        // with nrmax-th largest error estimate (to bisected next).
        //
LABEL_40:
        dqpsrt(limit, last, maxerr, errmax, elist, iord, nrmax);
        // jump out of for-loop
        if (errsum <= errbnd) goto LABEL_170;
        if (ier != 0) goto LABEL_150;
        if (last == 2 && extall) goto LABEL_120;
        if (noext) goto LABEL_140;
        if (!extall) goto LABEL_50;
        erlarg -= erlast;
        if (std::abs(b1 - a1) > small) erlarg += erro12;
        if (extrap) goto LABEL_70;
        //
        // test whether the interval to be bisected next is the
        // smallest interval.
        //
LABEL_50:
        width = std::abs(ARRAYF(blist, maxerr) - ARRAYF(alist, maxerr));
        if (width > small) goto LABEL_140;
        if (extall) goto LABEL_60;
        //
        // test whether we can start with the extrapolation procedure
        // (we do this if we integrate over the next interval with
        // use of a gauss-kronrod rule - see subroutine dqc25f).
        //
        small  = small * 0.5;
        if (0.25*width*domega > 2.0) goto LABEL_140;
        extall = true;
        goto LABEL_130;
LABEL_60:
        extrap = true;
        nrmax  = 2;
LABEL_70:
        if (ierro == 3 || erlarg <= ertest) goto LABEL_90;
        //
        // the smallest interval has the largest error.
        // before bisecting decrease the sum of the errors over
        // the larger intervals (erlarg) and perform extrapolation.
        //
        jupbnd = last;
        if (last > (limit / 2 + 2)) jupbnd = limit + 3 - last;
        id     = nrmax;
        for (k = id; k <= jupbnd; ++k) {
            maxerr = ARRAYF(iord,   nrmax);
            errmax = ARRAYF(elist, maxerr);
            if (std::abs(ARRAYF(blist, maxerr) - ARRAYF(alist, maxerr)) > small) goto LABEL_140;
            nrmax++;
        }
        // 
        // perform extrapolation.
        //
        LABEL_90:
        numrl2++;
        ARRAYF(rlist2, numrl2) = area;
        if (numrl2 < 3) goto LABEL_110;
        dqelg(numrl2, rlist2, reseps, abseps, res3la, nres);
        ktmin++;
        if (ktmin > 5 && abserr < 1.0e-3*errsum) ier = 5;
        if (abseps >= abserr) goto LABEL_100;
        ktmin  = 0;
        abserr = abseps;
        result = reseps;
        correc = erlarg;
        ertest = std::max(epsabs, epsrel * std::abs(reseps));
        // jump out of for-loop
        if (abserr <= ertest) goto LABEL_150;
        //
        // prepare bisection of the smallest interval.
        //
LABEL_100:
        if (numrl2 == 1) noext = true;
        if (ier == 5) goto LABEL_150;
LABEL_110:
        maxerr = ARRAYF(iord, 1);
        errmax = ARRAYF(elist, maxerr);
        nrmax  = 1;
        extrap = false;
        small *= 0.5;
        erlarg = errsum;
        goto LABEL_140;
LABEL_120:
        small  *= 0.5;
        numrl2++;
        ARRAYF(rlist2, numrl2) = area;
LABEL_130:
        ertest = errbnd;
        erlarg = errsum;
LABEL_140:
        continue;
    }
    //
    // set the final result.
    // ---------------------
    //
LABEL_150:
    if (abserr == oflow || nres == 0) goto LABEL_170;
    if (ier + ierro == 0) goto LABEL_165;
    if (ierro == 3) abserr += correc;
    if (ier == 0) ier = 3;
    if (result != 0.0 && area != 0.0) goto LABEL_160;
    if (abserr > errsum) goto LABEL_170;
    if (area == 0.0) goto LABEL_190;
    goto LABEL_165;
LABEL_160:
    if (abserr/std::abs(result) > errsum/std::abs(area)) goto LABEL_170;
    //
    // test on divergence.
    //
LABEL_165:
    if (ksgn == -1 && std::max(std::abs(result), std::abs(area)) <= defabs*1.0e-2) goto LABEL_190;
    if (1.0e-2 > result/area || result/area > 1.0e2 || errsum >= std::abs(area)) ier = 6;
    goto LABEL_190;
    //
    // compute global integral sum.
    //
LABEL_170:
    result = 0.0;
    for (k = 1; k <= last; ++k) {
        result += ARRAYF(rlist, k);
    }
    abserr = errsum;
LABEL_190:
    if (ier > 2) ier -= 1;
LABEL_200:
    if (integr == 2 && omega < 0.0) result = -result;
LABEL_999:
    return;
}

void dqaws(QUADPACK_CPP_FUNCTION f, const double a, const double b, const double alfa, const double beta, const int integr, const double epsabs,
           const double epsrel, double &result, double &abserr, int &neval, int &ier, const int limit, const int lenw, int &last, int *iwork,
           double *work, void *user_data)
{
    int lvl, l1, l2, l3;
    // first executable statement  dqaws
    ier    = 6;
    neval  = 0;
    last   = 0;
    result = 0.0;
    abserr = 0.0;
    if (limit < 2 || lenw < limit*4) goto LABEL_10;
    //
    // prepare call for dqawse.
    //
    l1 = limit + 1;
    l2 = limit + l1;
    l3 = limit + l2;
    //
    dqawse(f, a, b, alfa, beta, integr, epsabs, epsrel, limit, result, abserr, neval, ier, 
        &ARRAYF(work, 1), &ARRAYF(work, l1), &ARRAYF(work, l2), &ARRAYF(work, l3), iwork, last, user_data);
    //
    // call error handler if necessary
    //
    lvl = 0;
LABEL_10:
    if (ier == 6) ier = 1;
    if (ier != 0) {
        std::string messg = "abnormal return from dqaws";
        xerror(messg, ier, lvl);
    }
    return;
}

void dqawse(QUADPACK_CPP_FUNCTION f, const double a, const double b, const double alfa, const double beta, const int integr, const double epsabs,
            const double epsrel, const int limit, double &result, double &abserr, int &neval, int &ier, double *alist, double *blist, double *rlist,
            double *elist, int *iord, int &last, void *user_data)
{
    int iroff1, iroff2, k, maxerr, nev, nrmax;
    double area, area1, area12, area2, a1, a2, b1, b2, centre, errbnd, errmax;
    double error1, erro12, error2, errsum, resas1, resas2;
    double ri[25], rj[25], rh[25], rg[25];
    double epmach, uflow;
    // first executable statement  dqawse
    epmach = ARRAYF(d1mach, 4);
    uflow  = ARRAYF(d1mach, 1);
    //
    // test on validity of parameters
    // ------------------------------
    //
    ier              = 6;
    neval            = 0;
    last             = 0;
    ARRAYF(rlist, 1) = 0.0;
    ARRAYF(elist, 1) = 0.0;
    ARRAYF(iord,  1) = 0;
    result           = 0.0;
    abserr           = 0.0;
    if (b <= a || (epsabs == 0.0 && epsrel < std::max(50.0*epmach, 5.0e-29)) || 
        alfa <= -1.0 || beta <= -1.0 || integr < 1 || integr > 4 || limit < 2) goto LABEL_999;
    ier              = 0;
    //
    // compute the modified chebyshev moments.
    //
    dqmomo(alfa, beta, ri, rj, rg, rh, integr);
    //
    // integrate over the intervals (a,(a+b)/2) and ((a+b)/2,b).
    //
    centre = 0.5 * (b + a);
    dqc25s(f, a, b, a, centre, alfa, beta, ri, rj, rg, rh, area1, error1, resas1, integr, nev, user_data);
    neval  = nev;
    dqc25s(f, a, b, centre, b, alfa, beta, ri, rj, rg, rh, area2, error2, resas2, integr, nev, user_data);
    last   = 2;
    neval += nev;
    result = area1 + area2;
    abserr = error1 + error2;
    //
    // test on accuracy.
    //
    errbnd = std::max(epsabs, epsrel * std::abs(result));
    //
    // initialization
    // --------------
    //
    if (error2 > error1) goto LABEL_10;
    ARRAYF(alist, 1) = a;
    ARRAYF(alist, 2) = centre;
    ARRAYF(blist, 1) = centre;
    ARRAYF(blist, 2) = b;
    ARRAYF(rlist, 1) = area1;
    ARRAYF(rlist, 2) = area2;
    ARRAYF(elist, 1) = error1;
    ARRAYF(elist, 2) = error2;
    goto LABEL_20;
LABEL_10:
    ARRAYF(alist, 1) = centre;
    ARRAYF(alist, 2) = a;
    ARRAYF(blist, 1) = b;
    ARRAYF(blist, 2) = centre;
    ARRAYF(rlist, 1) = area2;
    ARRAYF(rlist, 2) = area1;
    ARRAYF(elist, 1) = error2;
    ARRAYF(elist, 2) = error1;
LABEL_20:
    ARRAYF(iord,  1) = 1;
    ARRAYF(iord,  2) = 2;
    if (limit == 2) ier = 1;
    if (abserr <= errbnd || ier == 1) goto LABEL_999;
    errmax           = ARRAYF(elist, 1);
    maxerr           = 1;
    nrmax            = 1;
    area             = result;
    errsum           = abserr;
    iroff1           = 0;
    iroff2           = 0;
    //
    // main for-loop
    // -------------
    //
    for (last = 3; last <= limit; ++last) {
        //
        // bisect the subinterval with largest error estimate.
        //
        a1 = ARRAYF(alist, maxerr);
        b1 = 0.5 * (ARRAYF(alist, maxerr) + ARRAYF(blist, maxerr));
        a2 = b1;
        b2 = ARRAYF(blist, maxerr);
        dqc25s(f, a, b, a1, b1, alfa, beta, ri, rj, rg, rh, area1, error1, resas1, integr, nev, user_data);
        neval += nev;
        dqc25s(f, a, b, a2, b2, alfa, beta, ri, rj, rg, rh, area2, error2, resas2, integr, nev, user_data);
        neval += nev;
        //
        // improve previous approximations integral and error
        // and test for accuracy.
        //
        area12 = area1  + area2;
        erro12 = error1 + error2;
        errsum = errsum + erro12 - errmax;
        area   = area   + area12 - ARRAYF(rlist, maxerr);
        if (a == a1 || b == b2) goto LABEL_30;
        if (resas1 == error1 || resas2 == error2) goto LABEL_30;
        //
        // test for roundoff error.
        //
        if (std::abs(ARRAYF(rlist, maxerr) - area12) < 1.0e-5 * std::abs(area12) && erro12 >= 0.99 * errmax) iroff1++;
        if (last > 10 && erro12 > errmax) iroff2++;
LABEL_30:
        ARRAYF(rlist, maxerr) = area1;
        ARRAYF(rlist,   last) = area2;
        //
        // test on accuracy.
        //
        errbnd = std::max(epsabs, epsrel * std::abs(area));
        if (errsum <= errbnd) goto LABEL_35;
        //
        // set error flag in the case that the number of interval
        // bisections exceeds limit.
        //
        if (last == limit) ier = 1;
        //
        // set error flag in the case of roundoff error.
        //
        if (iroff1 >= 6 || iroff2 >= 20) ier = 2;
        //
        // set error flag in the case of bad integrand behaviour
        // at interior points of integration range.
        //
        if(std::max(std::abs(a1), std::abs(b2)) <= (1.0 + 1.0e2*epmach) * (std::abs(a2) + 1.0e3*uflow)) ier = 3;
        //
        // append the newly-created intervals to the list.
        //
LABEL_35:
        if (error2 > error1) goto LABEL_40;
        ARRAYF(alist,   last) = a2;
        ARRAYF(blist, maxerr) = b1;
        ARRAYF(blist,   last) = b2;
        ARRAYF(elist, maxerr) = error1;
        ARRAYF(elist,   last) = error2;
        goto LABEL_50;
LABEL_40:
        ARRAYF(alist, maxerr) = a2;
        ARRAYF(alist,   last) = a1;
        ARRAYF(blist,   last) = b1;
        ARRAYF(rlist, maxerr) = area2;
        ARRAYF(rlist,   last) = area1;
        ARRAYF(elist, maxerr) = error2;
        ARRAYF(elist,   last) = error1;
        //
        // call subroutine dqpsrt to maintain the descending ordering
        // in the list of error estimates and select the subinterval
        // with largest error estimate (to be bisected next).
        //
LABEL_50:
        dqpsrt(limit, last, maxerr, errmax, elist, iord, nrmax);
        // jump out of for-loop
        if (ier != 0 || errsum <= errbnd) goto LABEL_70;
    }
    //
    // compute final result.
    // ---------------------
    //
LABEL_70:
    result = 0.0;
    for (k = 1; k <= last; ++k) {
        result += ARRAYF(rlist, k);
    }
    abserr = errsum;
LABEL_999:
    return;
}

/**
 * @fn void dgtsl
*/
void dgtsl(const int n, double *c, double *d, double *e, double *b, int &info)
{
    int k, kb, kp1, nm1, nm2;
    double t;
    // begin block permitting ...exits to 100
    info = 0;
    ARRAYF(c, 1) = ARRAYF(d, 1);
    nm1          = n - 1;
    if (nm1 < 1) goto LABEL_40;
    ARRAYF(d, 1) = ARRAYF(e, 1);
    ARRAYF(e, 1) = 0.0;
    ARRAYF(e, n) = 0.0;
    //
    for (k = 1; k <= nm1; ++k) {
        kp1 = k + 1;
        //
        // find the largest of the two rows
        //
        if (std::abs(ARRAYF(c, kp1)) < std::abs(ARRAYF(c, k))) goto LABEL_10;
        //
        // interchange row
        //
        t              = ARRAYF(c, kp1);
        ARRAYF(c, kp1) = ARRAYF(c,   k);
        ARRAYF(c,   k) = t;
        t              = ARRAYF(d, kp1);
        ARRAYF(d, kp1) = ARRAYF(d,   k);
        ARRAYF(d,   k) = t;
        t              = ARRAYF(e, kp1);
        ARRAYF(e, kp1) = ARRAYF(e,   k);
        ARRAYF(e,   k) = t;
        t              = ARRAYF(b, kp1);
        ARRAYF(b, kp1) = ARRAYF(b,   k);
        ARRAYF(b,   k) = t;
LABEL_10:
        //
        // zero elements
        //
        if (ARRAYF(c, k) != 0.0) goto LABEL_20;
        info = k;
        // ............exit
        goto LABEL_100;
LABEL_20:
        t              = -ARRAYF(c, kp1) / ARRAYF(c, k);
        ARRAYF(c, kp1) =  ARRAYF(d, kp1) + t*ARRAYF(d, k);
        ARRAYF(d, kp1) =  ARRAYF(e, kp1) + t*ARRAYF(e, k);
        ARRAYF(e, kp1) =  0.0;
        ARRAYF(b, kp1) =  ARRAYF(b, kp1) + t*ARRAYF(b, k);
    }
LABEL_40:
    if (ARRAYF(c, n) != 0.0) goto LABEL_50;
    info = n;
    goto LABEL_90;
LABEL_50:
    //
    // back solve
    //
    nm2          = n - 2;
    ARRAYF(b, n) = ARRAYF(b, n) / ARRAYF(c, n);
    if (n == 1) goto LABEL_80;
    ARRAYF(b, nm1) = (ARRAYF(b, nm1) - ARRAYF(d, nm1)*ARRAYF(b, n)) / ARRAYF(c, nm1);
    if (nm2 < 1) goto LABEL_70;
    for (kb = 1; kb <= nm2; ++kb) {
        k            = nm2 - kb + 1;
        ARRAYF(b, k) = (ARRAYF(b, k) - ARRAYF(d, k)*ARRAYF(b, k+1) - ARRAYF(e, k)*ARRAYF(b, k+2)) / ARRAYF(c, k);
    }
LABEL_70:
LABEL_80:
LABEL_90:
LABEL_100:
    //
    return;
}

void dqc25c(QUADPACK_CPP_FUNCTION f, const double a, const double b, const double c, 
            double &result, double &abserr, int &krul, int &neval, void *user_data)
{
    int i, isym, k;
    double ak22, amom0, amom1, amom2, cc, centr, hlgth, p2, p3, p4, resabs, resasc, res12, res24, u;
    double fval[25], cheb12[13], cheb24[25];
    const int kp = 0;
    //
    // the vector x contains the values cos(k*pi/24),
    // k = 1, ..., 11, to be used for the chebyshev series
    // expansion of f
    //
    const double x[11] = {
        0.991444861373810411144557526928563,
        0.965925826289068286749743199728897,
        0.923879532511286756128183189396788,
        0.866025403784438646763723170752936,
        0.793353340291235164579776961501299,
        0.707106781186547524400844362104849,
        0.608761429008720639416097542898164,
        0.500000000000000000000000000000000,
        0.382683432365089771728459984030399,
        0.258819045102520762348898837624048,
        0.130526192220051591548406227895489
    };

    // first executable statement  dqc25c
    cc = (2.0*c - b - a) / (b - a);
    if (std::abs(cc) < 1.1) goto LABEL_10;
    //
    // apply the 15-point gauss-kronrod scheme.
    //
    krul--;
    dqk15w(f, dqwgtc, c, p2, p3, p4, kp, a, b, result, abserr, resabs, resasc, user_data);
    neval = 15;
    if (resasc == abserr) krul++;
    goto LABEL_50;
    //
    // use the generalized clenshaw-curtis method.
    //
LABEL_10:
    hlgth = 0.5 * (b - a);
    centr = 0.5 * (b + a);
    neval = 25;
    ARRAYF(fval,  1) = 0.5 * f(hlgth + centr, user_data);
    ARRAYF(fval, 13) =       f(centr,         user_data);
    ARRAYF(fval, 25) = 0.5 * f(centr - hlgth, user_data);
    for (i = 2; i <= 12; ++i) {
        u    = hlgth * ARRAYF(x, i-1);
        isym = 26 - i;
        ARRAYF(fval,    i) = f(u + centr, user_data);
        ARRAYF(fval, isym) = f(centr - u, user_data);
    }
    //
    // compute the chebyshev series expansion.
    //
    dqcheb(x, fval, cheb12, cheb24);
    //
    // the modified chebyshev moments are computed by forward
    // recursion, using amom0 and amom1 as starting values.
    //
    amom0 = std::log(std::abs((1.0 - cc) / (1.0 + cc)));
    amom1 = 2.0 + cc * amom0;
    res12 = ARRAYF(cheb12, 1) * amom0 + ARRAYF(cheb12, 2) * amom1;
    res24 = ARRAYF(cheb24, 1) * amom0 + ARRAYF(cheb24, 2) * amom1;
    for (k = 3; k <= 13; ++k) {
        amom2 = 2.0 * cc * amom1 - amom0;
        ak22  = (k - 2) * (k - 2);
        if ((k / 2) * 2 == k) amom2 = amom2 - 4.0 / (ak22 - 1.0);
        res12 = res12 + ARRAYF(cheb12, k) * amom2;
        res24 = res24 + ARRAYF(cheb24, k) * amom2;
        amom0 = amom1;
        amom1 = amom2;
    }
    for (k = 14; k <= 25; ++k) {
        amom2 = 2.0 * cc * amom1 - amom0;
        ak22  = (k - 2) * (k - 2);
        if ((k / 2) * 2 == k) amom2 = amom2 - 4.0 / (ak22 - 1.0);
        res24 = res24 +  ARRAYF(cheb24, k)*amom2;
        amom0 = amom1;
        amom1 = amom2;
    }
    result = res24;
    abserr = std::abs(res24 - res12);
LABEL_50:
    return;
}

void dqc25f(QUADPACK_CPP_FUNCTION f, const double a, const double b, const double omega, const int integr, 
            const int nrmom, const int maxp1, const int ksave, double &result, double &abserr, int &neval, 
            double &resabs, double &resasc, int &momcom, double *chebmo, void *user_data)
{
    int i, iers, isym, j, k, m, noequ, noeq1;
    double ac, an, an2, as, asap, ass, centr, conc, cons, cospar, estc, ests;
    double hlgth, parint, par2, par22, resc12, resc24, ress12, ress24, sinpar;
    double p2, p3, p4;
    double cheb12[13], cheb24[25], d[25], d1[25], d2[25], fval[25], v[28];
    double oflow;
    const int chebmo_nrows = 25;
    //
    // the vector x contains the values cos(k*pi/24)
    // k = 1, ...,11, to be used for the chebyshev expansion of f
    //
    const double x[11] = {
        0.991444861373810411144557526928563,
        0.965925826289068286749743199728897,
        0.923879532511286756128183189396788,
        0.866025403784438646763723170752936,
        0.793353340291235164579776961501299,
        0.707106781186547524400844362104849,
        0.608761429008720639416097542898164,
        0.500000000000000000000000000000000,
        0.382683432365089771728459984030399,
        0.258819045102520762348898837624048,
        0.130526192220051591548406227895489
    };
    // first executable statement  dqc25f
    oflow  = ARRAYF(d1mach, 2);
    //
    centr  = 0.5 * (b + a);
    hlgth  = 0.5 * (b - a);
    parint = omega * hlgth;
    //
    // compute the integral using the 15-point gauss-kronrod
    // formula if the value of the parameter in the integrand
    // is small.
    //
    if(std::abs(parint) > 2.0) goto LABEL_10;
    dqk15w(f, dqwgtf, omega, p2, p3, p4, integr, a, b, result, abserr, resabs, resasc, user_data);
    neval = 15;
    goto LABEL_170;
    //
    // compute the integral using the generalized clenshaw-
    // curtis method.
    //
LABEL_10:
    conc   = hlgth * std::cos(centr * omega);
    cons   = hlgth * std::sin(centr * omega);
    resasc = oflow;
    neval  = 25;
    //
    // check whether the chebyshev moments for this interval
    // have already been computed.
    // 
    if (nrmom < momcom || ksave == 1) goto LABEL_120;
    //
    // compute a new set of chebyshev moments.
    // 
    m      = momcom + 1;
    par2   = parint * parint;
    par22  = par2 + 2.0;
    sinpar = std::sin(parint);
    cospar = std::cos(parint);
    //
    // compute the chebyshev moments with respect to cosine.
    //
    ARRAYF(v, 1) = 2.0 * sinpar / parint; 
    ARRAYF(v, 2) = (8.0 * cospar + (par2 + par2 - 8.0) * sinpar / parint) / par2;
    ARRAYF(v, 3) = (32.0 * (par2 - 12.0) * cospar + (2.0 * ((par2 - 80.0) * par2 + 192.0) * sinpar) / parint) / (par2 * par2);
    ac =  8.0 * cospar;
    as = 24.0 * parint * sinpar;
    if (std::abs(parint) > 24.0) goto LABEL_30;
    //
    // compute the chebyshev moments as the solutions of a
    // boundary value problem with 1 initial value (v(3)) and 1
    // end value (computed using an asymptotic formula).
    //
    noequ = 25;
    noeq1 = noequ - 1;
    an    = 6.0;
    for (k = 1; k <= noeq1; ++k) {
        an2             = an * an;
        ARRAYF(d,    k) = -2.0 * (an2 - 4.0) * (par22 - an2 - an2);
        ARRAYF(d2,   k) = (an - 1.0) * (an - 2.0) * par2;
        ARRAYF(d1, k+1) = (an + 3.0) * (an + 4.0) * par2;
        ARRAYF(v,  k+3) = as - (an2 - 4.0) * ac;
        an             += 2.0;
    }
    an2                = an * an;
    ARRAYF(d,   noequ) = -2.0 * (an2 - 4.0) * (par22 - an2 - an2);
    ARRAYF(v, noequ+3) = as - (an2 - 4.0) * ac;
    ARRAYF(v,       4) = ARRAYF(v, 4) - 56.0 * par2 * ARRAYF(v, 3);
    ass                = parint * sinpar;
    asap               = (((((210.0*par2 - 1.0)*cospar - (105.0*par2 - 63.0)*ass) / an2 - (1.0 - 15.0*par2)*cospar
                        + 15.0*ass) / an2 - cospar + 3.0*ass) / an2 - cospar) / an2;
    ARRAYF(v, noequ+3) = ARRAYF(v, noequ+3) - 2.0 * asap * par2 * (an - 1.0) * (an - 2.0);
    //
    // solve the tridiagonal system by means of gaussian
    // elimination with partial pivoting.
    //
    dgtsl(noequ, d1, d, d2, &ARRAYF(v, 4), iers);
    goto LABEL_50;
    //
    // compute the chebyshev moments by means of forward
    // recursion.
    //
LABEL_30:
    an = 4.0;
    for (i = 4; i <= 13; ++i) {
        an2          = an * an;
        ARRAYF(v, i) = ((an2 - 4.0) * (2.0 * (par22 - an2 - an2) * ARRAYF(v, i-1) - ac) 
                     + as - par2 * (an + 1.0) * (an + 2.0) * ARRAYF(v, i-2)) / (par2 * (an - 1.0) * (an - 2.0));
        an          += 2.0;
    }
LABEL_50:
    for (j = 1; j <= 13; ++j) {
        MATF(chebmo, chebmo_nrows, m, 2*j-1) = ARRAYF(v, j);
    }
    //
    // compute the chebyshev moments with respect to sine.
    //
    ARRAYF(v, 1) = 2.0 * (sinpar - parint * cospar) / par2;
    ARRAYF(v, 2) = (18.0 - 48.0 / par2) * sinpar / par2 + (-2.0 + 48.0 / par2) * cospar / parint;
    ac           = -24.0 * parint * cospar;
    as           = -8.0 * sinpar;
    if (std::abs(parint) > 24.0) goto LABEL_80;
    //
    // compute the chebyshev moments as the solutions of a boundary
    // value problem with 1 initial value (v(2)) and 1 end value
    // (computed using an asymptotic formula).
    //
    an = 5.0;
    for (k = 1; k <= noeq1; ++k) {
        an2             = an * an;
        ARRAYF(d,    k) = -2.0 * (an2 - 4.0) * (par22 - an2 - an2);
        ARRAYF(d2,   k) = (an - 1.0) * (an - 2.0) * par2;
        ARRAYF(d1, k+1) = (an + 3.0) * (an + 4.0) * par2;
        ARRAYF(v,  k+2) = ac + (an2 - 4.0) * as;
        an             += 2.0;
    }
    an2                = an * an;
    ARRAYF(d,   noequ) = -2.0 * (an2 - 4.0) * (par22 - an2 - an2);
    ARRAYF(v, noequ+2) = ac + (an2 - 4.0) * as;
    ARRAYF(v,       3) = ARRAYF(v, 3) - 42.0 * par2 * ARRAYF(v, 2);
    ass                = parint * cospar;
    asap               = (((((105.0*par2 - 63.0)*ass + (210.0*par2 - 1.0)*sinpar) / an2 + (15.0*par2 - 1.0)*sinpar - 15.0*ass) / an2 - 3.0*ass - sinpar) / an2 - sinpar) / an2;
    ARRAYF(v, noequ+2) = ARRAYF(v, noequ+2) - 2.0 * asap * par2 * (an - 1.0) *(an - 2.0);
    //
    // solve the tridiagonal system by means of gaussian
    // elimination with partial pivoting.
    //
    dgtsl(noequ, d1, d, d2, &ARRAYF(v, 3), iers);
    goto LABEL_100;
    //
    // compute the chebyshev moments by means of forward recursion.
    //
LABEL_80:
    an = 3.0;
    for (i = 3; i <= 12; ++i) {
        an2          = an * an;
        ARRAYF(v, i) = ((an2 - 4.0) * (2.0 * (par22 - an2 - an2) * ARRAYF(v, i-1) + as) + ac - par2 * (an + 1.0) * (an + 2.0) * ARRAYF(v, i-2)) / (par2 * (an - 1.0) * (an - 2.0));
        an          += 2.0;
    }
LABEL_100:
    for (j = 1; j <= 12; ++j) {
        MATF(chebmo, chebmo_nrows, m, 2*j) = ARRAYF(v, j);
    }
LABEL_120:
    if (nrmom < momcom) m = nrmom + 1;
    if (momcom < maxp1 - 1 && nrmom >= momcom) momcom += 1;
    //
    // compute the coefficients of the chebyshev expansions
    // of degrees 12 and 24 of the function f.
    //
    ARRAYF(fval,  1) = 0.5 * f(centr + hlgth, user_data);
    ARRAYF(fval, 13) =       f(centr,         user_data);
    ARRAYF(fval, 25) = 0.5 * f(centr - hlgth, user_data);
    for (i = 2; i <= 12; ++i) {
        isym = 26 - i;
        ARRAYF(fval,    i) = f(centr + hlgth * ARRAYF(x, i-1), user_data);
        ARRAYF(fval, isym) = f(centr - hlgth * ARRAYF(x, i-1), user_data);
    }
    dqcheb(x, fval, cheb12, cheb24);
    //
    // compute the integral and error estimates.
    //
    resc12 = ARRAYF(cheb12, 13) * MATF(chebmo, chebmo_nrows, m, 13);
    ress12 = 0.0;
    k      = 11;
    for (j = 1; j <= 6; ++j) {
        resc12 = resc12 + ARRAYF(cheb12,   k) * MATF(chebmo, chebmo_nrows, m,   k);
        ress12 = ress12 + ARRAYF(cheb12, k+1) * MATF(chebmo, chebmo_nrows, m, k+1);
        k     -= 2;
    }
    resc24 = ARRAYF(cheb24, 25) * MATF(chebmo, maxp1, m, 25);
    ress24 = 0.0;
    resabs = std::abs(ARRAYF(cheb24, 25));
    k      = 23;
    for (j = 1; j <= 12; ++j) {
        resc24 = resc24 + ARRAYF(cheb24,   k) * MATF(chebmo, chebmo_nrows, m,   k);
        ress24 = ress24 + ARRAYF(cheb24, k+1) * MATF(chebmo, chebmo_nrows, m, k+1);
        resabs = resabs + std::abs(ARRAYF(cheb24, k)) + std::abs(ARRAYF(cheb24, k+1));
        k -= 2;
    }
    estc   = std::abs(resc24 - resc12);
    ests   = std::abs(ress24 - ress12);
    resabs = resabs * std::abs(hlgth);
    if (integr == 2) goto LABEL_160;
    result = conc*resc24 - cons*ress24;
    abserr = std::abs(conc*estc) + std::abs(cons*ests);
    goto LABEL_170;
LABEL_160:
    result = conc*ress24 + cons*resc24;
    abserr = std::abs(conc*ests) + std::abs(cons*estc);
LABEL_170:
    return;
}

void dqc25s(QUADPACK_CPP_FUNCTION f, const double a, const double b, const double bl, const double br,
            const double alfa, const double beta, const double *ri, const double *rj, const double *rg,
            const double *rh, double &result, double &abserr, double &resasc, const int integr, int &nev,
            void *user_data)
{
    int i, isym;
    double centr, dc, factor, fix, hlgth, resabs, res12, res24, u;
    double cheb12[13], cheb24[25], fval[25];
    //
    // the vector x contains the values cos(k*pi/24)
    // k = 1, ..., 11, to be used for the computation of the
    // chebyshev series expansion of f.
    //
    const double x[11]= {
        0.991444861373810411144557526928563,
        0.965925826289068286749743199728897,
        0.923879532511286756128183189396788,
        0.866025403784438646763723170752936,
        0.793353340291235164579776961501299,
        0.707106781186547524400844362104849,
        0.608761429008720639416097542898164,
        0.500000000000000000000000000000000,
        0.382683432365089771728459984030399,
        0.258819045102520762348898837624048,
        0.130526192220051591548406227895489
    };
    // first executable statement  dqc25s
    nev = 25;
    if (bl == a && (alfa != 0.0 || integr == 2 || integr == 4)) goto LABEL_10;
    if (br == b && (beta != 0.0 || integr == 3 || integr == 4)) goto LABEL_140;
    //
    // if a.gt.bl and b.lt.br, apply the 15-point gauss-kronrod
    // scheme.
    //
    dqk15w(f, dqwgts, a, b, alfa, beta, integr, bl, br, result, abserr, resabs, resasc, user_data);
    nev = 15;
    goto LABEL_270;
    //
    // this part of the program is executed only if a = bl.
    // ----------------------------------------------------
    //
    // compute the chebyshev series expansion of the
    // following function
    // f1 = (0.5*(b+b-br-a)-0.5*(br-a)*x)**beta
    //        *f(0.5*(br-a)*x+0.5*(br+a))
    //
LABEL_10:
    hlgth            = 0.5 * (br - bl);
    centr            = 0.5 * (br + bl);
    fix              = b - centr;
    ARRAYF(fval,  1) = 0.5 * f(hlgth + centr, user_data) * std::pow(fix - hlgth, beta);
    ARRAYF(fval, 13) =       f(centr,         user_data) * std::pow(fix,         beta);
    ARRAYF(fval, 25) = 0.5 * f(centr - hlgth, user_data) * std::pow(fix + hlgth, beta);
    for (i = 2; i <= 12; ++i) {
        u                  = hlgth * ARRAYF(x, i-1);
        isym               = 26 - i;
        ARRAYF(fval,    i) = f(u + centr, user_data) * std::pow(fix - u, beta);
        ARRAYF(fval, isym) = f(centr - u, user_data) * std::pow(fix + u, beta);
    }
    factor = std::pow(hlgth, alfa + 1.0);
    result = 0.0;
    abserr = 0.0;
    res12  = 0.0;
    res24  = 0.0;
    if (integr > 2) goto LABEL_70;
    dqcheb(x, fval, cheb12, cheb24);
    //
    // integr = 1  (or 2)
    //
    for (i = 1; i <= 13; ++i) {
        res12 = res12 + ARRAYF(cheb12, i) * ARRAYF(ri, i);
        res24 = res24 + ARRAYF(cheb24, i) * ARRAYF(ri, i);
    }
    for (i = 14; i <= 25; ++i) {
        res24 = res24 + ARRAYF(cheb24, i) * ARRAYF(ri, i);
    }
    if (integr == 1) goto LABEL_130;
    //
    // integr = 2
    //
    dc     = std::log(br - bl);
    result = res24 * dc;
    abserr = std::abs((res24 - res12) * dc);
    res12  = 0.0;
    res24  = 0.0;
    for (i = 1; i <= 13; ++i) {
        res12 = res12 + ARRAYF(cheb12, i) * ARRAYF(rg, i);
        res24 = res12 + ARRAYF(cheb24, i) * ARRAYF(rg, i);
    }
    for (i = 14; i <= 25; ++i) {
        res24 = res24 + ARRAYF(cheb24, i) * ARRAYF(rg, i);
    }
    goto LABEL_130;
    //
    // compute the chebyshev series expansion of the
    // following function
    // f4 = f1*log(0.5*(b+b-br-a)-0.5*(br-a)*x)
    //
LABEL_70:
    ARRAYF(fval,  1) = ARRAYF(fval,  1) * std::log(fix - hlgth);
    ARRAYF(fval, 13) = ARRAYF(fval, 13) * std::log(fix);
    ARRAYF(fval, 25) = ARRAYF(fval, 25) * std::log(fix + hlgth);
    for (i = 2; i <= 12; ++i) {
        u                  = hlgth * ARRAYF(x, i-1);
        isym               = 26 - i;
        ARRAYF(fval,    i) = ARRAYF(fval,    i) * std::log(fix - u);
        ARRAYF(fval, isym) = ARRAYF(fval, isym) * std::log(fix + u);
    }
    dqcheb(x, fval, cheb12, cheb24);
    //
    // integr = 3  (or 4)
    //
    for (i = 1; i <= 13; ++i) {
        res12 = res12 + ARRAYF(cheb12, i) * ARRAYF(ri, i);
        res24 = res24 + ARRAYF(cheb24, i) * ARRAYF(ri, i);
    }
    for (i = 14; i <= 25; ++i) {
        res24 = res24 + ARRAYF(cheb24, i) * ARRAYF(ri, i);
    }
    if (integr == 3) goto LABEL_130;
    //
    // integr = 4
    //
    dc     = std::log(br - bl);
    result = res24 * dc;
    abserr = std::abs((res24 - res12) * dc);
    res12  = 0.0;
    res24  = 0.0;
    for (i = 1; i <= 13; ++i) {
        res12 = res12 + ARRAYF(cheb12, i) * ARRAYF(rg, i);
        res24 = res24 + ARRAYF(cheb24, i) * ARRAYF(rg, i);
    }
    for (i = 14; i <= 25; ++i) {
        res24 = res24 + ARRAYF(cheb24, i) * ARRAYF(rg, i);
    }
LABEL_130:
    result = (result + res24) * factor;
    abserr = (abserr + std::abs(res24 - res12)) * factor;
    goto LABEL_270;
    //
    // this part of the program is executed only if b = br.
    // ----------------------------------------------------
    //
    // compute the chebyshev series expansion of the
    // following function
    // f2 = (0.5*(b+bl-a-a)+0.5*(b-bl)*x)**alfa
    //      *f(0.5*(b-bl)*x+0.5*(b+bl))
    //
LABEL_140:
    hlgth            = 0.5 * (br - bl);
    centr            = 0.5 * (br + bl);
    fix              = centr - a;
    ARRAYF(fval,  1) = 0.5 * f(hlgth + centr, user_data) * std::pow(fix + hlgth, alfa);
    ARRAYF(fval, 13) =       f(centr,         user_data) * std::pow(fix,         alfa);
    ARRAYF(fval, 25) = 0.5 * f(centr - hlgth, user_data) * std::pow(fix - hlgth, alfa);
    for (i = 2; i <= 12; ++i) {
        u                  = hlgth * ARRAYF(x, i-1);
        isym               = 26 - i;
        ARRAYF(fval,    i) = f(u + centr, user_data) * std::pow(fix + u, alfa);
        ARRAYF(fval, isym) = f(centr - u, user_data) * std::pow(fix - u, alfa);
    }
    factor = std::pow(hlgth, beta + 1.0);
    result = 0.0;
    abserr = 0.0;
    res12  = 0.0;
    res24  = 0.0;
    if (integr == 2 || integr == 4) goto LABEL_200;
    //
    // integr = 1  (or 3)
    //
    dqcheb(x, fval, cheb12, cheb24);
    for (i = 1; i <= 13; ++i) {
        res12 = res12 + ARRAYF(cheb12, i) * ARRAYF(rj, i);
        res24 = res24 + ARRAYF(cheb24, i) * ARRAYF(rj, i);
    }
    for (i = 14; i <= 25; ++i) {
        res24 = res24 + ARRAYF(cheb24, i) * ARRAYF(rj, i);
    }
    if (integr == 1) goto LABEL_260;
    //
    // integr = 3
    //
    dc     = std::log(br - bl);
    result = res24 * dc;
    abserr = std::abs((res24 - res12) * dc);
    res12  = 0.0;
    res24  = 0.0;
    for (i = 1; i <= 13; ++i) {
        res12 = res12 + ARRAYF(cheb12, i) * ARRAYF(rh, i);
        res24 = res24 + ARRAYF(cheb24, i) * ARRAYF(rh, i);
    }
    for (i = 14; i <= 25; ++i) {
        res24 = res24 + ARRAYF(cheb24, i) * ARRAYF(rh, i);
    }
    goto LABEL_260;
    //
    // compute the chebyshev series expansion of the
    // following function
    // f3 = f2*log(0.5*(b-bl)*x+0.5*(b+bl-a-a))
    //
LABEL_200:
    ARRAYF(fval,  1) = ARRAYF(fval,  1) * std::log(hlgth + fix);
    ARRAYF(fval, 13) = ARRAYF(fval, 13) * std::log(fix);
    ARRAYF(fval, 25) = ARRAYF(fval, 25) * std::log(fix - hlgth);
    for (i = 2; i <= 12; ++i) {
        u                  = hlgth * ARRAYF(x, i-1);
        isym               = 26 - i;
        ARRAYF(fval,    i) = ARRAYF(fval,    i) * std::log(u + fix);
        ARRAYF(fval, isym) = ARRAYF(fval, isym) * std::log(fix - u);
    }
    dqcheb(x, fval, cheb12, cheb24);
    //
    // integr = 2  (or 4)
    //
    for (i = 1; i <= 13; ++i) {
        res12 = res12 + ARRAYF(cheb12, i) * ARRAYF(rj, i);
        res24 = res24 + ARRAYF(cheb24, i) * ARRAYF(rj, i);
    }
    for (i = 14; i <= 25; ++i) {
        res24 = res24 + ARRAYF(cheb24, i) * ARRAYF(rj, i);
    }
    if (integr == 2) goto LABEL_260;
    dc     = std::log(br - bl);
    result = res24 * dc;
    abserr = std::abs((res24 - res12) * dc);
    res12  = 0.0;
    res24  = 0.0;
    //
    // integr = 4
    //
    for (i = 1; i <= 13; ++i) {
        res12 = res12 + ARRAYF(cheb12, i) * ARRAYF(rh, i);
        res24 = res24 + ARRAYF(cheb24, i) * ARRAYF(rh, i);
    }
    for (i = 14; i <= 25; ++i) {
        res24 = res24 + ARRAYF(cheb24, i) * ARRAYF(rh, i);
    }
LABEL_260:
    result = (result + res24) * factor;
    abserr = (abserr + std::abs(res24 - res12)) * factor;
LABEL_270:
    return;
}


void dqcheb(const double *x, double *fval, double *cheb12, double *cheb24)
{
    int i, j;
    double alam, alam1, alam2, part1, part2, part3, v[12];
    // first executable statement  dqcheb
    for (i = 1; i <= 12; ++i) {
        j               = 26 - i;
        ARRAYF(v,    i) = ARRAYF(fval, i) - ARRAYF(fval, j);
        ARRAYF(fval, i) = ARRAYF(fval, i) + ARRAYF(fval, j);
    }
    alam1              = ARRAYF(v, 1) - ARRAYF(v, 9);
    alam2              = ARRAYF(x, 6) * (ARRAYF(v, 3) - ARRAYF(v, 7) - ARRAYF(v, 11));
    ARRAYF(cheb12,  4) = alam1 + alam2;
    ARRAYF(cheb12, 10) = alam1 - alam2;
    alam1              = ARRAYF(v, 2) - ARRAYF(v, 8) - ARRAYF(v, 10);
    alam2              = ARRAYF(v, 4) - ARRAYF(v, 6) - ARRAYF(v, 12);
    alam               = ARRAYF(x, 3)*alam1 + ARRAYF(x, 9)*alam2;
    ARRAYF(cheb24,  4) = ARRAYF(cheb12, 4) + alam;
    ARRAYF(cheb24, 22) = ARRAYF(cheb12, 4) - alam;
    alam               = ARRAYF(x, 9)*alam1 - ARRAYF(x, 3)*alam2;
    ARRAYF(cheb24, 10) = ARRAYF(cheb12, 10) + alam;
    ARRAYF(cheb24, 16) = ARRAYF(cheb12, 10) - alam;
    part1              = ARRAYF(x, 4) * ARRAYF(v, 5);
    part2              = ARRAYF(x, 8) * ARRAYF(v, 9);
    part3              = ARRAYF(x, 6) * ARRAYF(v, 7);
    alam1              = ARRAYF(v, 1) + part1 + part2;
    alam2              = ARRAYF(x, 2)*ARRAYF(v, 3) + part3 + ARRAYF(x, 10)*ARRAYF(v, 11);
    ARRAYF(cheb12,  2) = alam1 + alam2;
    ARRAYF(cheb12, 12) = alam1 - alam2;
    alam               = ARRAYF(x, 1)*ARRAYF(v, 2) + ARRAYF(x, 3)*ARRAYF(v,  4) + ARRAYF(x,  5)*ARRAYF(v,  6) 
                       + ARRAYF(x, 7)*ARRAYF(v, 8) + ARRAYF(x, 9)*ARRAYF(v, 10) + ARRAYF(x, 11)*ARRAYF(v, 12);
    ARRAYF(cheb24,  2) = ARRAYF(cheb12, 2) + alam;
    ARRAYF(cheb24, 24) = ARRAYF(cheb12, 2) - alam;
    alam               = ARRAYF(x, 11)*ARRAYF(v, 2) - ARRAYF(x, 9)*ARRAYF(v,  4) + ARRAYF(x, 7)*ARRAYF(v,  6)
                       - ARRAYF(x,  5)*ARRAYF(v, 8) + ARRAYF(x, 3)*ARRAYF(v, 10) - ARRAYF(x, 1)*ARRAYF(v, 12);
    ARRAYF(cheb24, 12) = ARRAYF(cheb12, 12) + alam;
    ARRAYF(cheb24, 14) = ARRAYF(cheb12, 12) - alam;
    alam1              = ARRAYF(v, 1) - part1 + part2;
    alam2              = ARRAYF(x, 10)*ARRAYF(v, 3) - part3 + ARRAYF(x, 2)*ARRAYF(v, 11);
    ARRAYF(cheb12, 6)  = alam1 + alam2;
    ARRAYF(cheb12, 8)  = alam1 - alam2;
    alam               = ARRAYF(x,  5)*ARRAYF(v, 2) - ARRAYF(x, 9)*ARRAYF(v,  4) - ARRAYF(x, 1)*ARRAYF(v,  6)
                       - ARRAYF(x, 11)*ARRAYF(v, 8) + ARRAYF(x, 3)*ARRAYF(v, 10) + ARRAYF(x, 7)*ARRAYF(v, 12);
    ARRAYF(cheb24,  6) = ARRAYF(cheb12, 6) + alam;
    ARRAYF(cheb24, 20) = ARRAYF(cheb12, 6) - alam;
    alam               = ARRAYF(x, 7)*ARRAYF(v, 2) - ARRAYF(x, 3)*ARRAYF(v,  4) - ARRAYF(x, 11)*ARRAYF(v,  6)
                       + ARRAYF(x, 1)*ARRAYF(v, 8) - ARRAYF(x, 9)*ARRAYF(v, 10) - ARRAYF(x,  5)*ARRAYF(v, 12);
    ARRAYF(cheb24,  8) = ARRAYF(cheb12, 8) + alam;
    ARRAYF(cheb24, 18) = ARRAYF(cheb12, 8) - alam;
    for (i = 1; i <= 6; ++i) {
        j               = 14 - i;
        ARRAYF(v,    i) = ARRAYF(fval, i) - ARRAYF(fval, j);
        ARRAYF(fval, i) = ARRAYF(fval, i) + ARRAYF(fval, j);
    }
    alam1              = ARRAYF(v, 1) + ARRAYF(x, 8)*ARRAYF(v, 5);
    alam2              = ARRAYF(x, 4) * ARRAYF(v, 3);
    ARRAYF(cheb12,  3) = alam1 + alam2;
    ARRAYF(cheb12, 11) = alam1 - alam2;
    ARRAYF(cheb12,  7) = ARRAYF(v, 1) - ARRAYF(v, 5);
    alam               = ARRAYF(x, 2)*ARRAYF(v, 2) + ARRAYF(x, 6)*ARRAYF(v, 4) + ARRAYF(x, 10)*ARRAYF(v, 6);
    ARRAYF(cheb24,  3) = ARRAYF(cheb12, 3) + alam;
    ARRAYF(cheb24, 23) = ARRAYF(cheb12, 3) - alam;
    alam               = ARRAYF(x, 6) * (ARRAYF(v, 2) - ARRAYF(v, 4) - ARRAYF(v, 6));
    ARRAYF(cheb24,  7) = ARRAYF(cheb12, 7) + alam;
    ARRAYF(cheb24, 19) = ARRAYF(cheb12, 7) - alam;
    alam               = ARRAYF(x, 10)*ARRAYF(v, 2) - ARRAYF(x, 6)*ARRAYF(v, 4) + ARRAYF(x, 2)*ARRAYF(v, 6);
    ARRAYF(cheb24, 11) = ARRAYF(cheb12, 11) + alam;
    ARRAYF(cheb24, 15) = ARRAYF(cheb12, 11) - alam;
    for (i = 1; i <= 3; ++i) {
        j               = 8 - i;
        ARRAYF(v,    i) = ARRAYF(fval, i) - ARRAYF(fval, j);
        ARRAYF(fval, i) = ARRAYF(fval, i) + ARRAYF(fval, j);
    }
    ARRAYF(cheb12, 5)  = ARRAYF(v,    1) + ARRAYF(x, 8)*ARRAYF(v,    3);
    ARRAYF(cheb12, 9)  = ARRAYF(fval, 1) - ARRAYF(x, 8)*ARRAYF(fval, 3);
    alam               = ARRAYF(x, 4) * ARRAYF(v, 2);
    ARRAYF(cheb24,  5) = ARRAYF(cheb12, 5) + alam;
    ARRAYF(cheb24, 21) = ARRAYF(cheb12, 5) - alam;
    alam               = ARRAYF(x, 8)*ARRAYF(fval, 2) - ARRAYF(fval, 4);
    ARRAYF(cheb24,  9) = ARRAYF(cheb12, 9) + alam;
    ARRAYF(cheb24, 17) = ARRAYF(cheb12, 9) - alam;
    ARRAYF(cheb12,  1) = ARRAYF(fval, 1) + ARRAYF(fval, 3);
    alam               = ARRAYF(fval, 2) + ARRAYF(fval, 4);
    ARRAYF(cheb24,  1) = ARRAYF(cheb12, 1) + alam;
    ARRAYF(cheb24, 25) = ARRAYF(cheb12, 1) - alam;
    ARRAYF(cheb12, 13) = ARRAYF(v, 1) - ARRAYF(v, 3);
    ARRAYF(cheb24, 13) = ARRAYF(cheb12, 13);
    alam               = 1.0 / 6.0;
    for (i = 2; i <= 12; ++i) {
        ARRAYF(cheb12, i) = ARRAYF(cheb12, i) * alam;
    }
    alam               = 0.5 * alam;
    ARRAYF(cheb12,  1) = ARRAYF(cheb12,  1) * alam;
    ARRAYF(cheb12, 13) = ARRAYF(cheb12, 13) * alam;
    for (i = 2; i <= 24; ++i) {
        ARRAYF(cheb24, i) = ARRAYF(cheb24, i) * alam;
    }
    ARRAYF(cheb24,  1) = 0.5 * alam * ARRAYF(cheb24,  1);
    ARRAYF(cheb24, 25) = 0.5 * alam * ARRAYF(cheb24, 25);
    return;
}

void dqelg(int &n, double *epstab, double &result, double &abserr, double *res3la, int &nres)
{
    double delta1, delta2, delta3, epsinf, error,
           err1, err2, err3, e0, e1, e1abs,
           e2, e3, res, ss, tol1, tol2, tol3;
    double epmach, oflow;
    int i, ib, ib2, ie, indx, k1, k2, k3, num, limexp, newelm;
    // first executable statement  dqelg
    epmach = ARRAYF(d1mach, 4);
    oflow  = ARRAYF(d1mach, 2);
    nres  += 1;
    abserr = oflow;
    result = ARRAYF(epstab, n);
    if(n < 3) goto LABEL_100;
    limexp = 50;
    ARRAYF(epstab, n+2) = ARRAYF(epstab, n);
    newelm = (n - 1) / 2;
    ARRAYF(epstab, n) = oflow;
    num    = n;
    k1     = n;
    for (i = 1; i <= newelm; ++i) {
        k2     = k1 - 1;
        k3     = k1 - 2;
        res    = ARRAYF(epstab, k1+2);
        e0     = ARRAYF(epstab,   k3);
        e1     = ARRAYF(epstab,   k2);
        e2     = res;
        e1abs  = std::abs(e1);
        delta2 = e2 - e1;
        err2   = std::abs(delta2);
        tol2   = std::max(std::abs(e2), e1abs) * epmach;
        delta3 = e1 - e0;
        err3   = std::abs(delta3);
        tol3   = std::max(e1abs, std::abs(e0)) * epmach;
        if (err2 > tol2 || err3 > tol3) goto LABEL_10;
        //
        // if e0, e1 and e2 are equal to within machine
        // accuracy, convergence is assumed.
        // result = e2
        // abserr = abs(e1-e0)+abs(e2-e1)
        //
        result = res;
        abserr = err2 + err3;
        // jump out of for-loop
        goto LABEL_100;
LABEL_10:
        e3           = ARRAYF(epstab, k1);
        ARRAYF(epstab, k1) = e1;
        delta1       = e1 - e3;
        err1         = std::abs(delta1);
        tol1         = std::max(e1abs, std::abs(e3)) * epmach;
        //
        // if two elements are very close to each other, omit
        // a part of the table by adjusting the value of n
        //
        if (err1 <= tol1 || err2 <= tol2 || err3 <= tol3) goto LABEL_20;
        ss = 1.0 / delta1 + 1.0 / delta2 - 1.0 / delta3;
        epsinf = std::abs(ss*e1);
        //
        // test to detect irregular behaviour in the table, and
        // eventually omit a part of the table adjusting the value
        // of n.
        // 
        if (epsinf > 1.0e-4) goto LABEL_30;
LABEL_20:
        n = i + i - 1;
        // jump out of for-loop
        goto LABEL_50;
        //
        // compute a new element and eventually adjust
        // the value of result.
        //
LABEL_30:
        res          = e1 + 1.0 / ss;
        ARRAYF(epstab, k1) = res;
        k1           = k1 - 2;
        error        = err2 + std::abs(res - e2) + err3;
        if (error > abserr) goto LABEL_40;
        abserr = error;
        result = res;
LABEL_40:
        continue;
    }
    //
    // shift the table.
    //
LABEL_50:
    if (n == limexp) n = 2 * (limexp / 2) - 1;
    ib = 1;
    if ((num / 2) * 2 == num) ib = 2;
    ie = newelm + 1;
    for (i = 1; i <= ie; ++i) {
        ib2 = ib + 2;
        ARRAYF(epstab, ib) = ARRAYF(epstab, ib2);
        ib = ib2;
    }
    if(num == n) goto LABEL_80;
    indx = num - n + 1;
    for (i = 1; i <= n; ++i) {
        ARRAYF(epstab, i) = ARRAYF(epstab, indx);
        indx += 1;
    }
LABEL_80:
    if(nres >= 4) goto LABEL_90;
    ARRAYF(res3la, nres) = result;
    abserr = oflow;
    goto LABEL_100;
    //
    // compute error estimate
    // 
LABEL_90:
    abserr = std::abs(result - ARRAYF(res3la, 3)) + std::abs(result - ARRAYF(res3la, 2)) + std::abs(result - ARRAYF(res3la, 1));
    ARRAYF(res3la, 1) = ARRAYF(res3la, 2);
    ARRAYF(res3la, 2) = ARRAYF(res3la, 3);
    ARRAYF(res3la, 3) = result;
LABEL_100:
    abserr = std::max(abserr, 5.0 * epmach * std::abs(result));
    return;
}

void dqk15(QUADPACK_CPP_FUNCTION f, const double a, const double b, double &result, double &abserr,
           double &resabs, double &resasc, void *user_data)
{
    int j, jtw, jtwm1;
    double absc, centr, dhlgth, fc, fsum, fval1, fval2, hlgth, resg, resk,reskh;
    double fv1[7], fv2[7];
    double uflow, epmach;
    //
    // the abscissae and weights are given for the interval (-1,1).
    // because of symmetry only the positive abscissae and their
    // corresponding weights are given.
    //
    // xgk    - abscissae of the 15-point kronrod rule
    //          xgk(2), xgk(4), ...  abscissae of the 7-point
    //          gauss rule
    //          xgk(1), xgk(3), ...  abscissae which are optimally
    //          added to the 7-point gauss rule
    // 
    // wgk    - weights of the 15-point kronrod rule
    //
    // wg     - weights of the 7-point gauss rule
    //
    // gauss quadrature weights and kronron quadrature abscissae and weights
    // as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
    // bell labs, nov. 1981.
    // 
    const double wg[4] = {
        0.129484966168869693270611432679082,
        0.279705391489276667901467771423780,
        0.381830050505118944950369775488975,
        0.417959183673469387755102040816327
    };
    const double xgk[8] = {
        0.991455371120812639206854697526329,
        0.949107912342758524526189684047851,
        0.864864423359769072789712788640926,
        0.741531185599394439863864773280788,
        0.586087235467691130294144838258730,
        0.405845151377397166906606412076961,
        0.207784955007898467600689403773245,
        0.000000000000000000000000000000000
    };
    const double wgk[8] = {
        0.022935322010529224963732008058970,
        0.063092092629978553290700663189204,
        0.104790010322250183839876322541518,
        0.140653259715525918745189590510238,
        0.169004726639267902826583426598550,
        0.190350578064785409913256402421014,
        0.204432940075298892414161999234649,
        0.209482141084727828012999174891714
    };
    // first executable statement  dqk15
    epmach = ARRAYF(d1mach, 4);
    uflow  = ARRAYF(d1mach, 1);
    //
    centr  = 0.5 * (a + b);
    hlgth  = 0.5 * (b - a);
    dhlgth = std::abs(hlgth);
    //
    // compute the 15-point kronrod approximation to
    // the integral, and estimate the absolute error.
    //
    fc     = f(centr, user_data);
    resg   = fc * ARRAYF(wg,  4);
    resk   = fc * ARRAYF(wgk, 8);
    resabs = std::abs(resk);
    for (j = 1; j <= 3; ++j) {
        jtw    = j * 2;
        absc   = hlgth * ARRAYF(xgk, jtw);
        fval1  = f(centr - absc, user_data);
        fval2  = f(centr + absc, user_data);
        ARRAYF(fv1, jtw) = fval1;
        ARRAYF(fv2, jtw) = fval2;
        fsum   = fval1 + fval2;
        resg   = resg + ARRAYF(wg,    j) * fsum;
        resk   = resk + ARRAYF(wgk, jtw) * fsum;
        resabs = resabs + ARRAYF(wgk, jtw) * (std::abs(fval1) + std::abs(fval2));
    }
    for (j = 1; j <= 4; ++j) {
        jtwm1  = 2 * j - 1;
        absc   = hlgth * ARRAYF(xgk, jtwm1);
        fval1  = f(centr - absc, user_data);
        fval2  = f(centr + absc, user_data);
        ARRAYF(fv1, jtwm1) = fval1;
        ARRAYF(fv2, jtwm1) = fval2;
        fsum   = fval1 + fval2;
        resk   = resk   + ARRAYF(wgk, jtwm1) * fsum;
        resabs = resabs + ARRAYF(wgk, jtwm1) * (std::abs(fval1) + std::abs(fval2));
    }
    reskh  = resk * 0.5;
    resasc = ARRAYF(wgk, 8) * std::abs(fc - reskh);
    for (j = 1; j <= 7; ++j) {
        resasc = resasc + ARRAYF(wgk, j) * (std::abs(ARRAYF(fv1, j) - reskh) + std::abs(ARRAYF(fv2, j) - reskh));
    }
    result = resk   * hlgth;
    resabs = resabs * dhlgth;
    resasc = resasc * dhlgth;
    abserr = std::abs((resk - resg) * hlgth);
    if (resasc != 0.0 && abserr != 0.0) abserr = resasc * std::min(1.0, std::pow(2.0e2*abserr/resasc, 1.5));
    if (resabs > uflow/(50.0*epmach)) abserr = std::max((epmach*50.0)*resabs, abserr);
    return;
}

void dqk15i(QUADPACK_CPP_FUNCTION f, const double boun, const int inf, const double a, const double b,
            double &result, double &abserr, double &resabs, double &resasc, void *user_data)
{
    int j;
    double absc, absc1, absc2, centr, dinf, fc, fsum, fval1, fval2, hlgth, resg;
    double resk, reskh, tabsc1, tabsc2;
    double fv1[7], fv2[7];
    double uflow, epmach;
    //
    // the abscissae and weights are supplied for the interval
    // (-1,1).  because of symmetry only the positive abscissae and
    // their corresponding weights are given.
    //
    // xgk    - abscissae of the 15-point kronrod rule
    //          xgk(2), xgk(4), ... abscissae of the 7-point
    //          gauss rule
    //          xgk(1), xgk(3), ...  abscissae which are optimally
    //          added to the 7-point gauss rule
    //
    // wgk    - weights of the 15-point kronrod rule
    //
    // wg     - weights of the 7-point gauss rule, corresponding
    //          to the abscissae xgk(2), xgk(4), ...
    //          wg(1), wg(3), ... are set to zero.
    //
    const double wg[8] = {
        0.000000000000000000000000000000000,
        0.129484966168869693270611432679082,
        0.000000000000000000000000000000000,
        0.279705391489276667901467771423780,
        0.000000000000000000000000000000000,
        0.381830050505118944950369775488975,
        0.000000000000000000000000000000000,
        0.417959183673469387755102040816327
    };
    const double xgk[8] = {
        0.991455371120812639206854697526329,
        0.949107912342758524526189684047851,
        0.864864423359769072789712788640926,
        0.741531185599394439863864773280788,
        0.586087235467691130294144838258730,
        0.405845151377397166906606412076961,
        0.207784955007898467600689403773245,
        0.000000000000000000000000000000000
    };
    const double wgk[8] = {
        0.022935322010529224963732008058970,
        0.063092092629978553290700663189204,
        0.104790010322250183839876322541518,
        0.140653259715525918745189590510238,
        0.169004726639267902826583426598550,
        0.190350578064785409913256402421014,
        0.204432940075298892414161999234649,
        0.209482141084727828012999174891714
    };
    // first executable statement  dqk15i
    epmach = ARRAYF(d1mach, 4);
    uflow  = ARRAYF(d1mach, 1);
    dinf   = std::max(1.0, static_cast<double>(inf));
    //
    centr  = 0.5 * (a + b);
    hlgth  = 0.5 * (b - a);
    tabsc1 = boun + dinf * (1.0 - centr) / centr;
    fval1  = f(tabsc1, user_data);
    if (inf == 2) fval1 = fval1 + f(-tabsc1, user_data);
    fc     = (fval1 / centr) / centr;
    //
    // compute the 15-point kronrod approximation to
    // the integral, and estimate the error.
    //
    resg   = ARRAYF(wg,  8) * fc;
    resk   = ARRAYF(wgk, 8) * fc;
    resabs = std::abs(resk);
    for (j = 1; j <= 7; ++j) {
        absc   = hlgth * ARRAYF(xgk, j);
        absc1  = centr - absc;
        absc2  = centr + absc;
        tabsc1 = boun + dinf * (1.0 - absc1) / absc1;
        tabsc2 = boun + dinf * (1.0 - absc2) / absc2;
        fval1  = f(tabsc1, user_data);
        fval2  = f(tabsc2, user_data);
        if(inf == 2) fval1 = fval1 + f(-tabsc1, user_data);
        if(inf == 2) fval2 = fval2 + f(-tabsc2, user_data);
        fval1  = (fval1 / absc1) / absc1;
        fval2  = (fval2 / absc2) / absc2;
        ARRAYF(fv1, j) = fval1;
        ARRAYF(fv2, j) = fval2;
        fsum   = fval1  + fval2;
        resg   = resg   + ARRAYF(wg,  j) * fsum;
        resk   = resk   + ARRAYF(wgk, j) * fsum;
        resabs = resabs + ARRAYF(wgk, j) * (std::abs(fval1) + std::abs(fval2));
    }
    reskh  = resk * 0.5;
    resasc = ARRAYF(wgk, 8) * std::abs(fc - reskh);
    for (j = 1; j <= 7; ++j) {
        resasc += ARRAYF(wgk, j) * (std::abs(ARRAYF(fv1, j) - reskh) + std::abs(ARRAYF(fv2, j) - reskh));
    }
    result = resk   * hlgth;
    resabs = resabs * hlgth;
    resasc = resasc * hlgth;
    abserr = std::abs((resk - resg) * hlgth);
    if (resasc != 0.0 && abserr != 0.0) abserr = resasc * std::min(1.0, std::pow(2.0e2 * abserr / resasc, 1.5));
    if (resabs > uflow / (50.0 * epmach)) abserr = std::max(50.0 * epmach * resabs, abserr);
    return;
}

void dqk15w(QUADPACK_CPP_FUNCTION f, QUADPACK_CPP_WEIGHT_FUNCTION w, const double p1, const double p2,
            const double p3, const double p4, const int kp, const double a, const double b, double &result,
            double &abserr, double &resabs, double &resasc, void *user_data)
{
    int j, jtw, jtwm1;
    double absc, absc1, absc2, centr, dhlgth, fc, fsum, fval1, fval2, hlgth;
    double resg, resk, reskh;
    double fv1[7], fv2[7];
    double epmach, uflow;
    //
    // the abscissae and weights are given for the interval (-1,1).
    // because of symmetry only the positive abscissae and their
    // corresponding weights are given.
    //
    // xgk    - abscissae of the 15-point gauss-kronrod rule
    //          xgk(2), xgk(4), ... abscissae of the 7-point
    //          gauss rule
    //          xgk(1), xgk(3), ... abscissae which are optimally
    //          added to the 7-point gauss rule
    //
    // wgk    - weights of the 15-point gauss-kronrod rule
    //
    // wg     - weights of the 7-point gauss rule
    //
    const double xgk[8] = {
        9.91455371120812639206854697526328516642e-1,
        9.49107912342758524526189684047851262401e-1,
        8.64864423359769072789712788640926201211e-1,
        7.41531185599394439863864773280788407074e-1,
        5.86087235467691130294144838258729598437e-1,
        4.05845151377397166906606412076961463347e-1,
        2.07784955007898467600689403773244913480e-1,
        0.00000000000000000000000000000000000000000
    };
    const double wgk[8] = {
        2.29353220105292249637320080589695919936e-2,
        6.30920926299785532907006631892042866651e-2,
        1.04790010322250183839876322541518017444e-1,
        1.40653259715525918745189590510237920400e-1,
        1.69004726639267902826583426598550284106e-1,
        1.90350578064785409913256402421013682826e-1,
        2.04432940075298892414161999234649084717e-1,
        2.09482141084727828012999174891714263698e-1
    };
    const double wg[4] = {
        1.29484966168869693270611432679082018329e-1,
        2.79705391489276667901467771423779582487e-1,
        3.81830050505118944950369775488975133878e-1,
        4.17959183673469387755102040816326530612e-1
    };
    // first executable statement  dqk15w
    epmach = ARRAYF(d1mach, 4);
    uflow  = ARRAYF(d1mach, 1);
    //
    centr  = 0.5 * (a + b);
    hlgth  = 0.5 * (b - a);
    dhlgth = std::abs(hlgth);
    //
    // compute the 15-point kronrod approximation to the
    // integral, and estimate the error.
    //
    fc   = f(centr, user_data) * w(centr, p1, p2, p3, p4, kp);
    resg = ARRAYF(wg,  4) * fc;
    resk = ARRAYF(wgk, 8) * fc;
    resabs = std::abs(resk);
    for (j = 1; j <= 3; ++j) {
        jtw    = j*2;
        absc   = hlgth * ARRAYF(xgk, jtw); 
        absc1  = centr - absc;
        absc2  = centr + absc;
        fval1  = f(absc1, user_data) * w(absc1, p1, p2, p3, p4, kp);
        fval2  = f(absc2, user_data) * w(absc2, p1, p2, p3, p4, kp);
        ARRAYF(fv1, jtw) = fval1;
        ARRAYF(fv2, jtw) = fval2;
        fsum   = fval1  + fval2;
        resg   = resg   + ARRAYF(wg,    j) * fsum;
        resk   = resk   + ARRAYF(wgk, jtw) * fsum;
        resabs = resabs + ARRAYF(wgk, jtw) * (std::abs(fval1) + std::abs(fval2));
    }
    for (j = 1; j <= 4; ++j) {
        jtwm1  = j*2 - 1;
        absc   = hlgth * ARRAYF(xgk, jtwm1);
        absc1  = centr - absc;
        absc2  = centr + absc;
        fval1  = f(absc1, user_data) * w(absc1, p1, p2, p3, p4, kp);
        fval2  = f(absc2, user_data) * w(absc2, p1, p2, p3, p4, kp);
        ARRAYF(fv1, jtwm1) = fval1;
        ARRAYF(fv2, jtwm1) = fval2;
        fsum   = fval1  + fval2;
        resk   = resk   + ARRAYF(wgk, jtwm1) * fsum;
        resabs = resabs + ARRAYF(wgk, jtwm1) * (std::abs(fval1) + std::abs(fval2));
    }
    reskh = resk * 0.5;
    resasc = ARRAYF(wgk, 8) * std::abs(fc - reskh);
    for (j = 1; j <= 7; ++j) {
        resasc = resasc + ARRAYF(wgk, j) * (std::abs(ARRAYF(fv1, j) - reskh) + std::abs(ARRAYF(fv2, j) - reskh));
    }
    result = resk   * hlgth;
    resabs = resabs * dhlgth;
    resasc = resasc * dhlgth;
    abserr = std::abs((resk - resg) * hlgth);
    if (resasc != 0.0 && abserr != 0.0)   abserr = resasc * std::min(1.0, std::pow(2.0e2 * abserr / resasc, 1.5));
    if (resabs > uflow / (50.0 * epmach)) abserr = std::max(50.0 * epmach * resabs, abserr);
    return;
}

void dqk21(QUADPACK_CPP_FUNCTION f, const double a, const double b, double &result, double &abserr,
           double &resabs, double &resasc, void *user_data)
{
    // double dhlgth, fc, fsum, fv1[10], fv2[10];
    // double uflow, epmach;
    // int j, jtw, jtwm1;
    // double centr;
    // double hlgth;
    // double absc;
    // double fval1;
    // double fval2;
    // double resg;
    // double resk;
    // double reskh;

    int j, jtw, jtwm1;
    double absc, centr, dhlgth, fc, fsum, fval1, fval2, hlgth, resg, resk, reskh;
    double fv1[10], fv2[10];
    double uflow, epmach;
    //
    // the abscissae and weights are given for the interval (-1,1).
    // because of symmetry only the positive abscissae and their
    // corresponding weights are given.
    //
    // xgk    - abscissae of the 21-point kronrod rule
    //          xgk(2), xgk(4), ...  abscissae of the 10-point
    //          gauss rule
    //          xgk(1), xgk(3), ...  abscissae which are optimally
    //          added to the 10-point gauss rule
    //
    // wgk    - weights of the 21-point kronrod rule
    //
    // wg     - weights of the 10-point gauss rule
    //
    //
    // gauss quadrature weights and kronron quadrature abscissae and weights
    // as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
    // c bell labs, nov. 1981.
    //
    const double wg[5] = {
        6.66713443086881375935688098933317928579e-2,
        1.49451349150580593145776339657697332403e-1,
        2.19086362515982043995534934228163192459e-1,
        2.69266719309996355091226921569469352860e-1,
        2.95524224714752870173892994651338329421e-1
    };
    const double xgk[11] = {
        9.95657163025808080735527280689002847921e-1,
        9.73906528517171720077964012084452053428e-1,
        9.30157491355708226001207180059508346225e-1,
        8.65063366688984510732096688423493048528e-1,
        7.80817726586416897063717578345042377163e-1,
        6.79409568299024406234327365114873575769e-1,
        5.62757134668604683339000099272694140843e-1,
        4.33395394129247190799265943165784162200e-1,
        2.94392862701460198131126603103865566163e-1,
        1.48874338981631210884826001129719984618e-1,
        0.00000000000000000000000000000000000000000
    };
    const double wgk[11] = {
        1.16946388673718742780643960621920483962e-2,
        3.25581623079647274788189724593897606174e-2,
        5.47558965743519960313813002445801763737e-2,
        7.50396748109199527670431409161900093952e-2,
        9.31254545836976055350654650833663443900e-2,
        1.09387158802297641899210590325804960272e-1,
        1.23491976262065851077958109831074159512e-1,
        1.34709217311473325928054001771706832761e-1,
        1.42775938577060080797094273138717060886e-1,
        1.47739104901338491374841515972068045524e-1,
        1.49445554002916905664936468389821203745e-1
    };
    // first executable statement  dqk21
    epmach = ARRAYF(d1mach, 4);
    uflow  = ARRAYF(d1mach, 1);
    //
    centr  = 0.5 * (a + b);
    hlgth  = 0.5 * (b - a);
    dhlgth = std::abs(hlgth);
    //
    // compute the 21-point kronrod approximation to
    // the integral, and estimate the absolute error.
    //
    resg   = 0.0;
    fc     = f(centr, user_data);
    resk   = ARRAYF(wgk, 11) * fc;
    resabs = std::abs(resk);
    for (j = 1; j <= 5; ++j) {
        jtw              = 2*j;
        absc             = hlgth * ARRAYF(xgk, jtw);
        fval1            = f(centr - absc, user_data);
        fval2            = f(centr + absc, user_data);
        ARRAYF(fv1, jtw) = fval1;
        ARRAYF(fv2, jtw) = fval2;
        fsum             = fval1 + fval2;
        resg             = resg + ARRAYF(wg,    j) * fsum;
        resk             = resk + ARRAYF(wgk, jtw) * fsum;
        resabs           = resabs + ARRAYF(wgk, jtw) * (std::abs(fval1) + std::abs(fval2));
    }
    for (j = 1; j <= 5; ++j) {
        jtwm1      = 2*j - 1;
        absc       = hlgth * ARRAYF(xgk, jtwm1);
        fval1      = f(centr - absc, user_data);
        fval2      = f(centr + absc, user_data);
        ARRAYF(fv1, jtwm1) = fval1;
        ARRAYF(fv2, jtwm1) = fval2;
        fsum       = fval1 + fval2;
        resk       = resk   + ARRAYF(wgk, jtwm1) * fsum;
        resabs     = resabs + ARRAYF(wgk, jtwm1) * (std::abs(fval1) + std::abs(fval2));
    }
    reskh  = resk * 0.5;
    resasc = ARRAYF(wgk, 11) * std::abs(fc - reskh);
    for (j = 1; j <= 10; ++j) {
        resasc += ARRAYF(wgk, j) * (std::abs(ARRAYF(fv1, j) - reskh) + std::abs(ARRAYF(fv2, j) - reskh));
    }
    result = resk * hlgth;
    resabs = resabs * dhlgth;
    resasc = resasc * dhlgth;
    abserr = std::abs((resk - resg) * hlgth);
    if (resasc != 0.0 && abserr != 0.0) abserr = resasc * std::min(1.0, std::pow(2.0e2 * abserr / resasc, 1.5));
    if (resabs > uflow / (50.0 * epmach)) abserr = std::max(50.0 * epmach * resabs, abserr);
    return;
}

void dqk31(QUADPACK_CPP_FUNCTION f, const double a, const double b, double &result, double &abserr,
           double &resabs, double &resasc, void *user_data)
{
    int j, jtw, jtwm1;
    double absc, centr, dhlgth, fc, fsum, fval1, fval2, hlgth, resg, resk, reskh;
    double fv1[15], fv2[15];
    double uflow, epmach;
    //
    // the abscissae and weights are given for the interval (-1,1).
    // because of symmetry only the positive abscissae and their
    // corresponding weights are given.
    //
    // xgk    - abscissae of the 31-point kronrod rule
    //          xgk(2), xgk(4), ...  abscissae of the 15-point
    //          gauss rule
    //          xgk(1), xgk(3), ...  abscissae which are optimally
    //          added to the 15-point gauss rule
    //
    // wgk    - weights of the 31-point kronrod rule
    //
    // wg     - weights of the 15-point gauss rule
    //
    //
    // gauss quadrature weights and kronron quadrature abscissae and weights
    // as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
    // bell labs, nov. 1981.
    //
    const double wg[8] = {
        3.07532419961172683546283935772044177217e-2,
        7.03660474881081247092674164506673384667e-2,
        1.07159220467171935011869546685869303416e-1,
        1.39570677926154314447804794511028322521e-1,
        1.66269205816993933553200860481208811131e-1,
        1.86161000015562211026800561866422824506e-1,
        1.98431485327111576456118326443839324819e-1,
        2.02578241925561272880620199967519314839e-1
    };
    const double xgk[16] = {
        9.98002298693397060285172840152271209073e-1,
        9.87992518020485428489565718586612581147e-1,
        9.67739075679139134257347978784337225283e-1,
        9.37273392400705904307758947710209471244e-1,
        8.97264532344081900882509656454495882832e-1,
        8.48206583410427216200648320774216851366e-1,
        7.90418501442465932967649294817947346862e-1,
        7.24417731360170047416186054613938009631e-1,
        6.50996741297416970533735895313274692547e-1,
        5.70972172608538847537226737253910641238e-1,
        4.85081863640239680693655740232350612866e-1,
        3.94151347077563369897207370981045468363e-1,
        2.99180007153168812166780024266388962662e-1,
        2.01194093997434522300628303394596207813e-1,
        1.01142066918717499027074231447392338787e-1,
        0.00000000000000000000000000000000000000000
    };
    const double wgk[16] = {
        5.37747987292334898779205143012764981831e-3,
        1.50079473293161225383747630758072680946e-2,
        2.54608473267153201868740010196533593973e-2,
        3.53463607913758462220379484783600481226e-2,
        4.45897513247648766082272993732796902233e-2,
        5.34815246909280872653431472394302967716e-2,
        6.20095678006706402851392309608029321904e-2,
        6.98541213187282587095200770991474757860e-2,
        7.68496807577203788944327774826590067221e-2,
        8.30805028231330210382892472861037896016e-2,
        8.85644430562117706472754436937743032123e-2,
        9.31265981708253212254868727473457185619e-2,
        9.66427269836236785051799076275893351367e-2,
        9.91735987217919593323931734846031310596e-2,
        1.00769845523875595044946662617569721916e-1,
        1.01330007014791549017374792767492546771e-1
    };
    // first executable statement  dqk31
    epmach = ARRAYF(d1mach, 4);
    uflow  = ARRAYF(d1mach, 1);
    //
    centr  = 0.5 * (a + b);
    hlgth  = 0.5 * (b - a);
    dhlgth = std::abs(hlgth);
    //
    // compute the 31-point kronrod approximation to
    // the integral, and estimate the absolute error.
    //
    fc     = f(centr, user_data);
    resg   = ARRAYF(wg,   8) * fc;
    resk   = ARRAYF(wgk, 16) * fc;
    resabs = std::abs(resk);
    for (j = 1; j <= 7; ++j) {
        jtw              = 2*j;
        absc             = hlgth * ARRAYF(xgk, jtw);
        fval1            = f(centr - absc, user_data);
        fval2            = f(centr + absc, user_data);
        ARRAYF(fv1, jtw) = fval1;
        ARRAYF(fv2, jtw) = fval2;
        fsum             = fval1 + fval2;
        resg             = resg + ARRAYF(wg,    j) * fsum;
        resk             = resk + ARRAYF(wgk, jtw) * fsum;
        resabs           = resabs + ARRAYF(wgk, jtw) * (std::abs(fval1) + std::abs(fval2));
    }
    for (j = 1; j <= 8; ++j) {
        jtwm1      = 2*j - 1;
        absc       = hlgth * ARRAYF(xgk, jtwm1);
        fval1      = f(centr - absc, user_data);
        fval2      = f(centr + absc, user_data);
        ARRAYF(fv1, jtwm1) = fval1;
        ARRAYF(fv2, jtwm1) = fval2;
        fsum       = fval1 + fval2;
        resk       = resk + ARRAYF(wgk, jtwm1) * fsum;
        resabs     = resabs + ARRAYF(wgk, jtwm1) * (std::abs(fval1) + std::abs(fval2));
    }
    reskh  = resk * 0.5;
    resasc = ARRAYF(wgk, 16) * std::abs(fc - reskh);
    for (j = 1; j <= 15; ++j) {
        resasc += ARRAYF(wgk, j) * (std::abs(ARRAYF(fv1, j) - reskh) + std::abs(ARRAYF(fv2, j) - reskh));
    }
    result = resk   * hlgth;
    resabs = resabs * dhlgth;
    resasc = resasc * dhlgth;
    abserr = std::abs((resk - resg) * hlgth);
    if (resasc != 0.0 && abserr != 0.0) abserr = resasc * std::min(1.0, std::pow(2.0e2 * abserr / resasc, 1.5));
    if (resabs > uflow / (50.0 * epmach)) abserr = std::max(50.0 * epmach * resabs, abserr);
    return;
}

void dqk41(QUADPACK_CPP_FUNCTION f, const double a, const double b, double &result, double &abserr,
           double &resabs, double &resasc, void *user_data)
{
    int j, jtw, jtwm1;
    double absc, centr, dhlgth, fc, fsum, fval1, fval2, hlgth, resg, resk, reskh;
    double fv1[20], fv2[20];
    double uflow, epmach;
    //
    // the abscissae and weights are given for the interval (-1,1).
    // because of symmetry only the positive abscissae and their
    // corresponding weights are given.
    //
    // xgk    - abscissae of the 41-point kronrod rule
    //          xgk(2), xgk(4), ...  abscissae of the 20-point
    //          gauss rule
    //          xgk(1), xgk(3), ...  abscissae which are optimally
    //          added to the 20-point gauss rule
    //
    // wgk    - weights of the 41-point kronrod rule
    //
    // wg     - weights of the 20-point gauss rule
    //
    //
    // gauss quadrature weights and kronron quadrature abscissae and weights
    // as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
    // bell labs, nov. 1981.
    //
    const double wg[10] = {
        1.76140071391521183118619623518528163621e-2,
        4.06014298003869413310399522749321098791e-2,
        6.26720483341090635695065351870416063516e-2,
        8.32767415767047487247581432220462061002e-2,
        1.01930119817240435036750135480349876167e-1,
        1.18194531961518417312377377711382287005e-1,
        1.31688638449176626898494499748163134916e-1,
        1.42096109318382051329298325067164933035e-1,
        1.49172986472603746787828737001969436693e-1,
        1.52753387130725850698084331955097593492e-1
    };
    const double xgk[21] = {
        9.98859031588277663838315576545863010000e-1,
        9.93128599185094924786122388471320278223e-1,
        9.81507877450250259193342994720216944567e-1,
        9.63971927277913791267666131197277221912e-1,
        9.40822633831754753519982722212443380274e-1,
        9.12234428251325905867752441203298113049e-1,
        8.78276811252281976077442995113078466711e-1,
        8.39116971822218823394529061701520685330e-1,
        7.95041428837551198350638833272787942959e-1,
        7.46331906460150792614305070355641590311e-1,
        6.93237656334751384805490711845931533386e-1,
        6.36053680726515025452836696226285936743e-1,
        5.75140446819710315342946036586425132814e-1,
        5.10867001950827098004364050955250998425e-1,
        4.43593175238725103199992213492640107840e-1,
        3.73706088715419560672548177024927237396e-1,
        3.01627868114913004320555356858592260615e-1,
        2.27785851141645078080496195368574624743e-1,
        1.52605465240922675505220241022677527912e-1,
        7.65265211334973337546404093988382110048e-2,
        0.00000000000000000000000000000000000000000
    };
    const double wgk[21] = {
        3.07358371852053150121829324603098748803e-3,
        8.60026985564294219866178795010234725213e-3,
        1.46261692569712529837879603088683561639e-2,
        2.03883734612665235980102314327547051228e-2,
        2.58821336049511588345050670961531429995e-2,
        3.12873067770327989585431193238007378878e-2,
        3.66001697582007980305572407072110084875e-2,
        4.16688733279736862637883059368947380440e-2,
        4.64348218674976747202318809261075168421e-2,
        5.09445739237286919327076700503449486648e-2,
        5.51951053482859947448323724197773291948e-2,
        5.91114008806395723749672206485942171364e-2,
        6.26532375547811680258701221742549805858e-2,
        6.58345971336184221115635569693979431472e-2,
        6.86486729285216193456234118853678017155e-2,
        7.10544235534440683057903617232101674129e-2,
        7.30306903327866674951894176589131127606e-2,
        7.45828754004991889865814183624875286161e-2,
        7.57044976845566746595427753766165582634e-2,
        7.63778676720807367055028350380610018008e-2,
        7.66007119179996564450499015301017408279e-2
    };
    // first executable statement  dqk31
    epmach = ARRAYF(d1mach, 4);
    uflow  = ARRAYF(d1mach, 1);
    //
    centr  = 0.5 * (a + b);
    hlgth  = 0.5 * (b - a);
    dhlgth = std::abs(hlgth);
    //
    // compute the 41-point gauss-kronrod approximation to
    // the integral, and estimate the absolute error.
    //
    resg   = 0.0;
    fc     = f(centr, user_data);
    resk   = ARRAYF(wgk, 21) * fc;
    resabs = std::abs(resk);
    for (j = 1; j <= 10; ++j) {
        jtw              = 2*j;
        absc             = hlgth * ARRAYF(xgk, jtw);
        fval1            = f(centr - absc, user_data);
        fval2            = f(centr + absc, user_data);
        ARRAYF(fv1, jtw) = fval1;
        ARRAYF(fv2, jtw) = fval2;
        fsum             = fval1 + fval2;
        resg             = resg + ARRAYF(wg,    j) * fsum;
        resk             = resk + ARRAYF(wgk, jtw) * fsum;
        resabs           = resabs + ARRAYF(wgk, jtw) * (std::abs(fval1) + std::abs(fval2));
    }
    for (j = 1; j <= 10; ++j) {
        jtwm1              = 2*j - 1;
        absc               = hlgth * ARRAYF(xgk, jtwm1);
        fval1              = f(centr - absc, user_data);
        fval2              = f(centr + absc, user_data);
        ARRAYF(fv1, jtwm1) = fval1;
        ARRAYF(fv2, jtwm1) = fval2;
        fsum               = fval1 + fval2;
        resk               = resk   + ARRAYF(wgk, jtwm1) * fsum;
        resabs             = resabs + ARRAYF(wgk, jtwm1) * (std::abs(fval1) + std::abs(fval2));
    }
    reskh  = resk * 0.5;
    resasc = ARRAYF(wgk, 21) * std::abs(fc - reskh);
    for (j = 1; j <= 20; ++j) {
        resasc += ARRAYF(wgk, j) * (std::abs(ARRAYF(fv1, j) - reskh) + std::abs(ARRAYF(fv2, j) - reskh));
    }
    result = resk   * hlgth;
    resabs = resabs * dhlgth;
    resasc = resasc * dhlgth;
    abserr = std::abs((resk - resg) * hlgth);
    if (resasc != 0.0 && abserr != 0.0) abserr = resasc * std::min(1.0, std::pow(2.0e2 * abserr / resasc, 1.5));
    if (resabs > uflow / (50.0 * epmach)) abserr = std::max(50.0 * epmach * resabs, abserr);
    return;
}

void dqk51(QUADPACK_CPP_FUNCTION f, const double a, const double b, double &result, double &abserr,
           double &resabs, double &resasc, void *user_data)
{
    int j, jtw, jtwm1;
    double absc, centr, dhlgth, fc, fsum, fval1, fval2, hlgth, resg, resk, reskh;
    double fv1[25], fv2[25];
    double uflow, epmach;
    //
    // the abscissae and weights are given for the interval (-1,1).
    // because of symmetry only the positive abscissae and their
    // corresponding weights are given.
    //
    // xgk    - abscissae of the 51-point kronrod rule
    //          xgk(2), xgk(4), ...  abscissae of the 25-point
    //          gauss rule
    //          xgk(1), xgk(3), ...  abscissae which are optimally
    //          added to the 25-point gauss rule
    //
    // wgk    - weights of the 51-point kronrod rule
    //
    // wg     - weights of the 25-point gauss rule
    //
    //
    // gauss quadrature weights and kronron quadrature abscissae and weights
    // as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
    // bell labs, nov. 1981.
    //
    const double wg[13] = {
        1.13937985010262879479029641132347736033e-2,
        2.63549866150321372619018152952991449360e-2,
        4.09391567013063126556234877116459536608e-2,
        5.49046959758351919259368915404733241601e-2,
        6.80383338123569172071871856567079685547e-2,
        8.01407003350010180132349596691113022902e-2,
        9.10282619829636498114972207028916533810e-2,
        1.00535949067050644202206890392685826988e-1,
        1.08519624474263653116093957050116619340e-1,
        1.14858259145711648339325545869555808641e-1,
        1.19455763535784772228178126512901047390e-1,
        1.22242442990310041688959518945851505835e-1,
        1.23176053726715451203902873079050142438e-1
    };
    const double xgk[26] = {
        9.99262104992609834193457486540340593705e-1,
        9.95556969790498097908784946893901617258e-1,
        9.88035794534077247637331014577406227072e-1,
        9.76663921459517511498315386479594067745e-1,
        9.61614986425842512418130033660167241692e-1,
        9.42974571228974339414011169658470531905e-1,
        9.20747115281701561746346084546330631575e-1,
        8.94991997878275368851042006782804954175e-1,
        8.65847065293275595448996969588340088203e-1,
        8.33442628760834001421021108693569569461e-1,
        7.97873797998500059410410904994306569409e-1,
        7.59259263037357630577282865204360976388e-1,
        7.17766406813084388186654079773297780598e-1,
        6.73566368473468364485120633247622175883e-1,
        6.26810099010317412788122681624517881020e-1,
        5.77662930241222967723689841612654067396e-1,
        5.26325284334719182599623778158010178037e-1,
        4.73002731445714960522182115009192041332e-1,
        4.17885382193037748851814394594572487093e-1,
        3.61172305809387837735821730127640667422e-1,
        3.03089538931107830167478909980339329200e-1,
        2.43866883720988432045190362797451586406e-1,
        1.83718939421048892015969888759528415785e-1,
        1.22864692610710396387359818808036805532e-1,
        6.15444830056850788865463923667966312817e-2,
        0.00000000000000000000000000000000000000000
    };
    const double wgk[26] = {
        1.98738389233031592650785188284340988943e-3,
        5.56193213535671375804023690106552207018e-3,
        9.47397338617415160720771052365532387165e-3,
        1.32362291955716748136564058469762380776e-2,
        1.68478177091282982315166675363363158404e-2,
        2.04353711458828354565682922359389736788e-2,
        2.40099456069532162200924891648810813929e-2,
        2.74753175878517378029484555178110786148e-2,
        3.07923001673874888911090202152285856009e-2,
        3.40021302743293378367487952295512032257e-2,
        3.71162714834155435603306253676198759960e-2,
        4.00838255040323820748392844670756464014e-2,
        4.28728450201700494768957924394951611020e-2,
        4.55029130499217889098705847526603930437e-2,
        4.79825371388367139063922557569147549836e-2,
        5.02776790807156719633252594334400844406e-2,
        5.23628858064074758643667121378727148874e-2,
        5.42511298885454901445433704598756068261e-2,
        5.59508112204123173082406863827473468203e-2,
        5.74371163615678328535826939395064719948e-2,
        5.86896800223942079619741758567877641398e-2,
        5.97203403241740599790992919325618538354e-2,
        6.05394553760458629453602675175654271623e-2,
        6.11285097170530483058590304162927119227e-2,
        6.14711898714253166615441319652641775865e-2,
        6.15808180678329350787598242400645531904e-2
    };
    // first executable statement  dqk31
    epmach = ARRAYF(d1mach, 4);
    uflow  = ARRAYF(d1mach, 1);
    //
    centr  = 0.5 * (a + b);
    hlgth  = 0.5 * (b - a);
    dhlgth = std::abs(hlgth);
    //
    // compute the 51-point kronrod approximation to
    // the integral, and estimate the absolute error.
    //
    fc     = f(centr, user_data);
    resg   = ARRAYF(wg, 13) * fc;
    resk   = ARRAYF(wgk, 26) * fc;
    resabs = std::abs(resk);
    for (j = 1; j <= 12; ++j) {
        jtw              = 2*j;
        absc             = hlgth * ARRAYF(xgk, jtw);
        fval1            = f(centr - absc, user_data);
        fval2            = f(centr + absc, user_data);
        ARRAYF(fv1, jtw) = fval1;
        ARRAYF(fv2, jtw) = fval2;
        fsum             = fval1 + fval2;
        resg             = resg + ARRAYF(wg,    j) * fsum;
        resk             = resk + ARRAYF(wgk, jtw) * fsum;
        resabs           = resabs + ARRAYF(wgk, jtw) * (std::abs(fval1) + std::abs(fval2));
    }
    for (j = 1; j <= 13; ++j) {
        jtwm1              = 2*j - 1;
        absc               = hlgth * ARRAYF(xgk, jtwm1);
        fval1              = f(centr - absc, user_data);
        fval2              = f(centr + absc, user_data);
        ARRAYF(fv1, jtwm1) = fval1;
        ARRAYF(fv2, jtwm1) = fval2;
        fsum               = fval1 + fval2;
        resk               = resk + ARRAYF(wgk, jtwm1) * fsum;
        resabs             = resabs + ARRAYF(wgk, jtwm1) * (std::abs(fval1) + std::abs(fval2));
    }
    reskh  = resk * 0.5;
    resasc = ARRAYF(wgk, 26) * std::abs(fc - reskh);
    for (j = 1; j <= 25; ++j) {
        resasc += ARRAYF(wgk, j) * (std::abs(ARRAYF(fv1, j) - reskh) + std::abs(ARRAYF(fv2, j) - reskh));
    }
    result = resk   * hlgth;
    resabs = resabs * dhlgth;
    resasc = resasc * dhlgth;
    abserr = std::abs((resk - resg) * hlgth);
    if (resasc != 0.0 && abserr != 0.0) abserr = resasc * std::min(1.0, std::pow(2.0e2 * abserr / resasc, 1.5));
    if (resabs > uflow / (50.0 * epmach)) abserr = std::max(50.0 * epmach * resabs, abserr);
    return;
}

void dqk61(QUADPACK_CPP_FUNCTION f, const double a, const double b, double &result, double &abserr,
           double &resabs, double &resasc, void *user_data)
{
    int j, jtw, jtwm1;
    double absc, centr, dhlgth, fc, fsum, fval1, fval2, hlgth, resg, resk, reskh;
    double fv1[30], fv2[30];
    double uflow, epmach;
    //
    // the abscissae and weights are given for the interval (-1,1).
    // because of symmetry only the positive abscissae and their
    // corresponding weights are given.
    //
    // xgk    - abscissae of the 61-point kronrod rule
    //          xgk(2), xgk(4), ...  abscissae of the 30-point
    //          gauss rule
    //          xgk(1), xgk(3), ...  abscissae which are optimally
    //          added to the 30-point gauss rule
    //
    // wgk    - weights of the 61-point kronrod rule
    //
    // wg     - weights of the 30-point gauss rule
    //
    //
    // gauss quadrature weights and kronron quadrature abscissae and weights
    // as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
    // bell labs, nov. 1981.
    //
    const double wg[15] = {
        7.96819249616660561546588347467362245048e-3,
        1.84664683110909591423021319120472690962e-2,
        2.87847078833233693497191796112920436396e-2,
        3.87991925696270495968019364463476920332e-2,
        4.84026728305940529029381404228075178153e-2,
        5.74931562176190664817216894020561287971e-2,
        6.59742298821804951281285151159623612374e-2,
        7.37559747377052062682438500221907341538e-2,
        8.07558952294202153546949384605297308759e-2,
        8.68997872010829798023875307151257025768e-2,
        9.21225222377861287176327070876187671969e-2,
        9.63687371746442596394686263518098650964e-2,
        9.95934205867952670627802821035694765299e-2,
        1.01762389748405504596428952168554044633e-1,
        1.02852652893558840341285636705415043868e-1
        };
    const double xgk[31] = {
        9.99484410050490637571325895705810819469e-1,
        9.96893484074649540271630050918695283341e-1,
        9.91630996870404594858628366109485724851e-1,
        9.83668123279747209970032581605662801940e-1,
        9.73116322501126268374693868423706884888e-1,
        9.60021864968307512216871025581797662930e-1,
        9.44374444748559979415831324037439121586e-1,
        9.26200047429274325879324277080474004086e-1,
        9.05573307699907798546522558925958319569e-1,
        8.82560535792052681543116462530225590057e-1,
        8.57205233546061098958658510658943856821e-1,
        8.29565762382768397442898119732501916439e-1,
        7.99727835821839083013668942322683240736e-1,
        7.67777432104826194917977340974503131695e-1,
        7.33790062453226804726171131369527645669e-1,
        6.97850494793315796932292388026640068382e-1,
        6.60061064126626961370053668149270753038e-1,
        6.20526182989242861140477556431189299207e-1,
        5.79345235826361691756024932172540495907e-1,
        5.36624148142019899264169793311072794164e-1,
        4.92480467861778574993693061207708795644e-1,
        4.47033769538089176780609900322854000162e-1,
        4.00401254830394392535476211542660633611e-1,
        3.52704725530878113471037207089373860654e-1,
        3.04073202273625077372677107199256553531e-1,
        2.54636926167889846439805129817805107883e-1,
        2.04525116682309891438957671002024709524e-1,
        1.53869913608583546963794672743255920419e-1,
        1.02806937966737030147096751318000592472e-1,
        5.14718425553176958330252131667225737491e-2,
        0.00000000000000000000000000000000000000
        };
    const double wgk[31] = {
        1.38901369867700762455159122675969968105e-3,
        3.89046112709988405126720184451550327852e-3,
        6.63070391593129217331982636975016813363e-3,
        9.27327965951776342844114689202436042127e-3,
        1.18230152534963417422328988532505928963e-2,
        1.43697295070458048124514324435800101958e-2,
        1.69208891890532726275722894203220923686e-2,
        1.94141411939423811734089510501284558514e-2,
        2.18280358216091922971674857383389934015e-2,
        2.41911620780806013656863707252320267604e-2,
        2.65099548823331016106017093350754143665e-2,
        2.87540487650412928439787853543342111447e-2,
        3.09072575623877624728842529430922726353e-2,
        3.29814470574837260318141910168539275106e-2,
        3.49793380280600241374996707314678750972e-2,
        3.68823646518212292239110656171359677370e-2,
        3.86789456247275929503486515322810502509e-2,
        4.03745389515359591119952797524681142161e-2,
        4.19698102151642461471475412859697577901e-2,
        4.34525397013560693168317281170732580746e-2,
        4.48148001331626631923555516167232437574e-2,
        4.60592382710069881162717355593735805947e-2,
        4.71855465692991539452614781810994864829e-2,
        4.81858617570871291407794922983045926058e-2,
        4.90554345550297788875281653672381736059e-2,
        4.97956834270742063578115693799423285392e-2,
        5.04059214027823468408930856535850289022e-2,
        5.08817958987496064922974730498046918534e-2,
        5.12215478492587721706562826049442082511e-2,
        5.14261285374590259338628792157812598296e-2,
        5.14947294294515675583404336470993075327e-2
        };
    epmach = ARRAYF(d1mach, 4);
    uflow  = ARRAYF(d1mach, 1);
    //
    centr  = 0.5 * (a + b);
    hlgth  = 0.5 * (b - a);
    dhlgth = std::abs(hlgth);
    //
    // compute the 51-point kronrod approximation to
    // the integral, and estimate the absolute error.
    //
    resg   = 0.0;
    fc     = f(centr, user_data);
    resk   = ARRAYF(wgk, 31) * fc;
    resabs = std::abs(resk);
    for (j = 1; j <= 15; ++j) {
        jtw              = 2*j;
        absc             = hlgth * ARRAYF(xgk, jtw);
        fval1            = f(centr - absc, user_data);
        fval2            = f(centr + absc, user_data);
        ARRAYF(fv1, jtw) = fval1;
        ARRAYF(fv2, jtw) = fval2;
        fsum             = fval1 + fval2;
        resg             = resg + ARRAYF(wg,    j) * fsum;
        resk             = resk + ARRAYF(wgk, jtw) * fsum;
        resabs           = resabs + ARRAYF(wgk, jtw) * (std::abs(fval1) + std::abs(fval2));
    }
    for (j = 1; j <= 15; ++j) {
        jtwm1              = 2*j - 1;
        absc               = hlgth * ARRAYF(xgk, jtwm1);
        fval1              = f(centr - absc, user_data);
        fval2              = f(centr + absc, user_data);
        ARRAYF(fv1, jtwm1) = fval1;
        ARRAYF(fv2, jtwm1) = fval2;
        fsum               = fval1 + fval2;
        resk               = resk + ARRAYF(wgk, jtwm1) * fsum;
        resabs             = resabs + ARRAYF(wgk, jtwm1) * (std::abs(fval1) + std::abs(fval2));
    }
    reskh  = resk * 0.5;
    resasc = ARRAYF(wgk, 31) * std::abs(fc - reskh);
    for (j = 1; j <= 30; ++j) {
        resasc += ARRAYF(wgk, j) * (std::abs(ARRAYF(fv1, j) - reskh) + std::abs(ARRAYF(fv2, j) - reskh));
    }
    result = resk   * hlgth;
    resabs = resabs * dhlgth;
    resasc = resasc * dhlgth;
    abserr = std::abs((resk - resg) * hlgth);
    if (resasc != 0.0 && abserr != 0.0) abserr = resasc * std::min(1.0, std::pow(2.0e2 * abserr / resasc, 1.5));
    if (resabs > uflow / (50.0 * epmach)) abserr = std::max(50.0 * epmach * resabs, abserr);
    return;
}

void dqmomo(const double alfa, const double beta, double *ri, double *rj, double *rg, double *rh, const int integr)
{
    double alfp1, alfp2, an, anm1, betp1, betp2, ralf, rbet;
    int i, im1;
    //
    // first executable statement  dqmomo
    alfp1 = alfa + 1.0;
    betp1 = beta + 1.0;
    alfp2 = alfa + 2.0;
    betp2 = beta + 2.0;
    ralf  = std::pow(2.0, alfp1);
    rbet  = std::pow(2.0, betp1);
    //
    // compute ri, rj using a forward recurrence relation.
    //
    ARRAYF(ri, 1) = ralf / alfp1;
    ARRAYF(rj, 1) = rbet / betp1;
    ARRAYF(ri, 2) = ARRAYF(ri, 1) * alfa / alfp2;
    ARRAYF(rj, 2) = ARRAYF(rj, 1) * beta / betp2;
    an            = 2.0;
    anm1          = 1.0;
    for (i = 3; i <= 25; ++i) {
        ARRAYF(ri, i) = -(ralf + an * (an - alfp2) * ARRAYF(ri, i-1)) / (anm1 * (an + alfp1));
        ARRAYF(rj, i) = -(rbet + an * (an - betp2) * ARRAYF(rj, i-1)) / (anm1 * (an + betp1));
        anm1          = an;
        an            = an + 1.0;
    }
    if (integr == 1) goto LABEL_70;
    if (integr == 3) goto LABEL_40;
    //
    // compute rg using a forward recurrence relation.
    //
    ARRAYF(rg, 1) = -ARRAYF(ri, 1) / alfp1;
    ARRAYF(rg, 2) = -(ralf + ralf) / (alfp2 * alfp2) - ARRAYF(rg, 1);
    an            = 2.0;
    anm1          = 1.0;
    im1           = 2;
    for (i = 3; i <= 25; ++i) {
        ARRAYF(rg, i) = -(an * (an - alfp2) * ARRAYF(rg, im1) - an * ARRAYF(ri, im1) + anm1 * ARRAYF(ri, i)) / (anm1 * (an + alfp1));
        anm1          = an;
        an           += 1.0;
        im1           = i;
    }
    if (integr == 2) goto LABEL_70;
    //
    // compute rh using a forward recurrence relation.
    //
LABEL_40:
    ARRAYF(rh, 1) = -ARRAYF(rj, 1) / betp1;
    ARRAYF(rh, 2) = -(rbet + rbet) / (betp2 * betp2) - ARRAYF(rh, 1);
    an            = 2.0;
    anm1          = 1.0;
    im1           = 2;
    for (i = 3; i <= 25; ++i) {
        ARRAYF(rh, i) = -(an * (an - betp2) * ARRAYF(rh, im1) - an * ARRAYF(rj, im1) + anm1 * ARRAYF(rj, i)) / (anm1 * (an + betp1));
        anm1          = an;
        an           += 1.0;
        im1           = i;
    }
    for (i = 2; i <= 25; i += 2) {
        ARRAYF(rh, i) = -ARRAYF(rh, i);
    }
LABEL_70:
    for (i = 2; i <= 25; i += 2) {
        ARRAYF(rj, i) = -ARRAYF(rj, i);
    }
LABEL_90:
    return;
}

void dqng(QUADPACK_CPP_FUNCTION f, const double a, const double b, const double epsabs, const double epsrel,
          double &result, double &abserr, int &neval, int &ier, void *user_data)
{
    // double dhlgth, fval1, fval2, fv1[5], fv2[5], fv3[5], fv4[5], reskh;
    // int ipx, k, l;
    // double centr;
    // double hlgth;
    // double fcentr;
    // double absc;
    // double fval;
    // double savfun[21];
    // double res10;
    // double res21;
    // double res43;
    // double res87;
    // double resabs;
    // double resasc;
    // double epmach, uflow;
    // std::string messg;
    int ipx, k, l;
    double fv1[5], fv2[5], fv3[5], fv4[5], savfun[21];
    double absc, centr, dhlgth, fcentr, fval, fval1, fval2, hlgth;
    double res10, res21, res43, res87, resabs, resasc, reskh;
    double epmach, uflow;
    std::string messg;
    //
    // the following data statements contain the
    // abscissae and weights of the integration rules used.
    //
    // x1      abscissae common to the 10-, 21-, 43- and 87-
    //         point rule
    // x2      abscissae common to the 21-, 43- and 87-point rule
    // x3      abscissae common to the 43- and 87-point rule
    // x4      abscissae of the 87-point rule
    // w10     weights of the 10-point formula
    // w21a    weights of the 21-point formula for abscissae x1
    // w21b    weights of the 21-point formula for abscissae x2
    // w43a    weights of the 43-point formula for abscissae x1, x3
    // w43b    weights of the 43-point formula for abscissae x3
    // w87a    weights of the 87-point formula for abscissae x1,
    //         x2, x3
    // w87b    weights of the 87-point formula for abscissae x4
    //
    //
    // gauss-kronrod-patterson quadrature coefficients for use in
    // quadpack routine qng.  these coefficients were calculated with
    // 101 decimal digit arithmetic by l. w. fullerton, bell labs, nov 1981.
    //
    const double x1[5] = {
        9.73906528517171720077964012084452053428e-1,
        8.65063366688984510732096688423493048528e-1,
        6.79409568299024406234327365114873575769e-1,
        4.33395394129247190799265943165784162200e-1,
        1.48874338981631210884826001129719984618e-1
    };
    const double w10[5] = {
        6.66713443086881375935688098933317928579e-2,
        1.49451349150580593145776339657697332403e-1,
        2.19086362515982043995534934228163192459e-1,
        2.69266719309996355091226921569469352860e-1,
        2.95524224714752870173892994651338329421e-1
    };
    const double x2[5] = {
        9.95657163025808080735527280689002847921e-1,
        9.30157491355708226001207180059508346225e-1,
        7.80817726586416897063717578345042377163e-1,
        5.62757134668604683339000099272694140843e-1,
        2.94392862701460198131126603103865566163e-1
    };
    const double w21a[5] = {
        3.25581623079647274788189724593897606174e-2,
        7.50396748109199527670431409161900093952e-2,
        1.09387158802297641899210590325804960272e-1,
        1.34709217311473325928054001771706832761e-1,
        1.47739104901338491374841515972068045524e-1
    };
    const double w21b[6] = {
        1.16946388673718742780643960621920483962e-2,
        5.47558965743519960313813002445801763737e-2,
        9.31254545836976055350654650833663443900e-2,
        1.23491976262065851077958109831074159512e-1,
        1.42775938577060080797094273138717060886e-1,
        1.49445554002916905664936468389821203745e-1
    };
    const double x3[11] = {
        0.999333360901932081394099323919911,
        0.987433402908088869795961478381209,
        0.954807934814266299257919200290473,
        0.900148695748328293625099494069092,
        0.825198314983114150847066732588520,
        0.732148388989304982612354848755461,
        0.622847970537725238641159120344323,
        0.499479574071056499952214885499755,
        0.364901661346580768043989548502644,
        0.222254919776601296498260928066212,
        0.074650617461383322043914435796506
    };
    const double w43a[10] = {
        0.016296734289666564924281974617663,
        0.037522876120869501461613795898115,
        0.054694902058255442147212685465005,
        0.067355414609478086075553166302174,
        0.073870199632393953432140695251367,
        0.005768556059769796184184327908655,
        0.027371890593248842081276069289151,
        0.046560826910428830743339154433824,
        0.061744995201442564496240336030883,
        0.071387267268693397768559114425516
    };
    const double w43b[12] = {
        0.001844477640212414100389106552965,
        0.010798689585891651740465406741293,
        0.021895363867795428102523123075149,
        0.032597463975345689443882222526137,
        0.042163137935191811847627924327955,
        0.050741939600184577780189020092084,
        0.058379395542619248375475369330206,
        0.064746404951445885544689259517511,
        0.069566197912356484528633315038405,
        0.072824441471833208150939535192842,
        0.074507751014175118273571813842889,
        0.074722147517403005594425168280423
    };
    const double x4[22] = {
        0.999902977262729234490529830591582,
        0.997989895986678745427496322365960,
        0.992175497860687222808523352251425,
        0.981358163572712773571916941623894,
        0.965057623858384619128284110607926,
        0.943167613133670596816416634507426,
        0.915806414685507209591826430720050,
        0.883221657771316501372117548744163,
        0.845710748462415666605902011504855,
        0.803557658035230982788739474980964,
        0.757005730685495558328942793432020,
        0.706273209787321819824094274740840,
        0.651589466501177922534422205016736,
        0.593223374057961088875273770349144,
        0.531493605970831932285268948562671,
        0.466763623042022844871966781659270,
        0.399424847859218804732101665817923,
        0.329874877106188288265053371824597,
        0.258503559202161551802280975429025,
        0.185695396568346652015917141167606,
        0.111842213179907468172398359241362,
        0.037352123394619870814998165437704
    };
    const double w87a[21] = {
        0.008148377384149172900002878448190,
        0.018761438201562822243935059003794,
        0.027347451050052286161582829741283,
        0.033677707311637930046581056957588,
        0.036935099820427907614589586742499,
        0.002884872430211530501334156248695,
        0.013685946022712701888950035273128,
        0.023280413502888311123409291030404,
        0.030872497611713358675466394126442,
        0.035693633639418770719351355457044,
        0.000915283345202241360843392549948,
        0.005399280219300471367738743391053,
        0.010947679601118931134327826856808,
        0.016298731696787335262665703223280,
        0.021081568889203835112433060188190,
        0.025370969769253827243467999831710,
        0.029189697756475752501446154084920,
        0.032373202467202789685788194889595,
        0.034783098950365142750781997949596,
        0.036412220731351787562801163687577,
        0.037253875503047708539592001191226
    };
    const double w87b[23] = {
        0.000274145563762072350016527092881,
        0.001807124155057942948341311753254,
        0.004096869282759164864458070683480,
        0.006758290051847378699816577897424,
        0.009549957672201646536053581325377,
        0.012329447652244853694626639963780,
        0.015010447346388952376697286041943,
        0.017548967986243191099665352925900,
        0.019938037786440888202278192730714,
        0.022194935961012286796332102959499,
        0.024339147126000805470360647041454,
        0.026374505414839207241503786552615,
        0.028286910788771200659968002987960,
        0.030052581128092695322521110347341,
        0.031646751371439929404586051078883,
        0.033050413419978503290785944862689,
        0.034255099704226061787082821046821,
        0.035262412660156681033782717998428,
        0.036076989622888701185500318003895,
        0.036698604498456094498018047441094,
        0.037120549269832576114119958413599,
        0.037334228751935040321235449094698,
        0.037361073762679023410321241766599
    };
    // first executable statement  dqng
    epmach = ARRAYF(d1mach, 4);
    uflow  = ARRAYF(d1mach, 1);
    //
    // test on validity of parameters
    // ------------------------------
    //
    result = 0.0;
    abserr = 0.0;
    neval  = 0;
    ier    = 6;
    if (epsabs <= 0.0 && epsrel < std::max(50.0 * epmach, 5.0e-29)) goto LABEL_80;
    hlgth  = 0.5 * (b - a);
    dhlgth = std::abs(hlgth);
    centr  = 0.5 * (b + a);
    fcentr = f(centr, user_data);
    neval  = 21;
    ier    = 1;
    //
    // compute the integral using the 10- and 21-point formula.
    //
    for (l = 1; l <= 3; ++l) {
        if (l == 1) goto LABEL_5;
        if (l == 2) goto LABEL_25;
        if (l == 3) goto LABEL_45;
LABEL_5:
        res10  = 0.0;
        res21  = ARRAYF(w21b, 6) * fcentr;
        resabs = ARRAYF(w21b, 6) * std::abs(fcentr);
        for (k = 1; k <= 5; ++k) {
            absc = hlgth * ARRAYF(x1, k);
            fval1             = f(centr + absc, user_data);
            fval2             = f(centr - absc, user_data);
            fval              = fval1 + fval2;
            res10             = res10  + ARRAYF( w10, k) * fval;
            res21             = res21  + ARRAYF(w21a, k) * fval;
            resabs            = resabs + ARRAYF(w21a, k) * (std::abs(fval1) + std::abs(fval2));
            ARRAYF(savfun, k) = fval;
            ARRAYF(   fv1, k) = fval1;
            ARRAYF(   fv2, k) = fval2;
        }
        ipx = 5;
        for (k = 1; k <= 5; ++k) {
            ipx                 = ipx + 1;
            absc                = hlgth * ARRAYF(x2, k);
            fval1               = f(centr + absc, user_data);
            fval2               = f(centr - absc, user_data);
            fval                = fval1 + fval2;
            res21               = res21 + ARRAYF(w21b, k) * fval;
            resabs              = resabs + ARRAYF(w21b, k) * (std::abs(fval1) + std::abs(fval2));
            ARRAYF(savfun, ipx) = fval;
            ARRAYF(   fv3,   k) = fval1;
            ARRAYF(   fv4,   k) = fval2;
        }
        //
        // test for convergence.
        //
        result = res21  * hlgth;
        resabs = resabs * dhlgth;
        reskh  = 0.5 * res21;
        resasc = ARRAYF(w21b, 6) * std::abs(fcentr - reskh);
        for (k = 1; k <= 5; ++k) {
            resasc = resasc + ARRAYF(w21a, k) * (std::abs(ARRAYF(fv1, k) - reskh) + std::abs(ARRAYF(fv2, k) - reskh))
                            + ARRAYF(w21b, k) * (std::abs(ARRAYF(fv3, k) - reskh) + std::abs(ARRAYF(fv4, k) - reskh));
        }
        abserr = std::abs((res21 - res10) * hlgth);
        resasc = resasc * dhlgth;
        goto LABEL_65;
        //
        // compute the integral using the 43-point formula.
        //
LABEL_25:
        res43 = ARRAYF(w43b, 12) * fcentr;
        neval = 43;
        for (k = 1; k <= 10; ++k) {
            res43 = res43 + ARRAYF(savfun, k) * ARRAYF(w43a, k);
        }
        for (k = 1; k <= 11; ++k) {
            ipx                 = ipx + 1;
            absc                = hlgth * ARRAYF(x3, k);
            fval                = f(absc + centr, user_data) + f(centr - absc, user_data);
            res43               = res43 + fval * ARRAYF(w43b, k);
            ARRAYF(savfun, ipx) = fval;
        }
        //
        // test for convergence.
        //
        result = res43 * hlgth;
        abserr = std::abs((res43 - res21) * hlgth);
        goto LABEL_65;
        //
        // compute the integral using the 87-point formula.
        //
LABEL_45:
        res87 = ARRAYF(w87b, 23) * fcentr;
        neval = 87;
        for (k = 1; k <= 21; ++k) {
            res87 = res87 + ARRAYF(savfun, k) * ARRAYF(w87a, k);
        }
        for (k = 1; k <= 22; ++k) {
            absc  = hlgth * ARRAYF(x4, k);
            res87 = res87 + ARRAYF(w87b, k) * (f(absc + centr, user_data) + f(centr - absc, user_data));
        }
        result = res87 * hlgth;
        abserr = std::abs((res87 - res43) * hlgth);
LABEL_65:
        if (resasc != 0.0 && abserr != 0.0)   abserr = resasc * std::min(1.0, std::pow(2.0e2 * abserr / resasc, 1.5));
        if (resabs > uflow / (50.0 * epmach)) abserr = std::max(50.0 * epmach * resabs, abserr);
        if (abserr <= std::max(epsabs, epsrel * std::abs(result))) ier = 0;
        // jump out of for-loop
        if (ier == 0) goto LABEL_999;
    }
LABEL_80:
    messg = "abnormal return from dqng";
    xerror(messg, ier, 0);
LABEL_999:
    return;
}

void dqpsrt(const int limit, const int last, int &maxerr, double &ermax, const double *elist, int *iord, int &nrmax)
{
    double errmax, errmin;
    int i, ibeg, ido, isucc, j, jbnd, jupbn, k;
    //
    // check whether the list contains more than
    // two error estimates.
    //
    // first executable statement  dqpsrt
    if (last > 2) goto LABEL_10;
    ARRAYF(iord, 1) = 1;
    ARRAYF(iord, 2) = 2;
    goto LABEL_90;
    //
    // this part of the routine is only executed if, due to a
    // difficult integrand, subdivision increased the error
    // estimate. in the normal case the insert procedure should
    // start after the nrmax-th largest error estimate.
    //
LABEL_10:
    errmax = ARRAYF(elist, maxerr);
    if (nrmax == 1) goto LABEL_30;
    ido = nrmax - 1;
    for (i = 1; i <= ido; ++i) {
        isucc = ARRAYF(iord, nrmax-1);
        // jump out of do-loop
        if (errmax <= ARRAYF(elist, isucc)) goto LABEL_30;
        ARRAYF(iord, nrmax) = isucc;
        nrmax--;
    }
    //
    // compute the number of elements in the list to be maintained
    // in descending order. this number depends on the number of
    // subdivisions still allowed.
    //
LABEL_30:
    jupbn = last;
    if (last > (limit / 2 + 2)) jupbn = limit + 3 - last;
    errmin = ARRAYF(elist, last);
    //
    // insert errmax by traversing the list top-down,
    // starting comparison from the element elist(iord(nrmax+1)).
    // 
    jbnd = jupbn - 1;
    ibeg = nrmax + 1;
    if (ibeg > jbnd) goto LABEL_50;
    for (i = ibeg; i <= jbnd; ++i) {
        isucc = ARRAYF(iord, i);
        // jump out of for-loop
        if(errmax >= ARRAYF(elist, isucc)) goto LABEL_60;
        ARRAYF(iord, i-1) = isucc;
    }
LABEL_50:
    ARRAYF(iord, jbnd)  = maxerr;
    ARRAYF(iord, jupbn) = last;
    goto LABEL_90;
    //
    // insert errmin by traversing the list bottom-up.
    //
LABEL_60:
    ARRAYF(iord, i-1) = maxerr;
    k = jbnd;
    for (j = i; j <= jbnd; ++j) {
        isucc = ARRAYF(iord, k);
        // jump out of for-loop
        if (errmin < ARRAYF(elist, isucc)) goto LABEL_80;
        ARRAYF(iord, k+1) = isucc;
        k--;
    }
    ARRAYF(iord, i) = last;
    goto LABEL_90;
LABEL_80:
    ARRAYF(iord, k+1) = last;
    //
    // set maxerr and ermax.
    //
LABEL_90:
    maxerr = ARRAYF(iord,   nrmax);
    ermax  = ARRAYF(elist, maxerr);
    return;
}

double dqwgtc(const double x, const double c, const double p2, const double p3, const double p4, const int kp)
{
    return 1.0 / (x - c);
}

double dqwgtf(const double x, const double omega, const double p2, const double p3, const double p4, const int integr)
{
    double omx;
    // first executable statement  dqwgtf
    omx = omega * x;
    if (integr == 1) {
        return std::cos(omx);
    } else if (integr == 2) {
        return std::sin(omx);
    } else {
        return 0.0;
    }
}

double dqwgts(const double x, const double a, const double b, const double alfa, const double beta, const int integr)
{
    double xma = x - a;
    double bmx = b - x;
    double coef_dqwgts = std::pow(xma, alfa) * std::pow(bmx, beta);

    switch (integr)
    {
        case 1:
            return coef_dqwgts;
        case 2:
            return coef_dqwgts * std::log(xma);
        case 3:
            return coef_dqwgts * std::log(bmx);
        case 4:
            return coef_dqwgts * std::log(xma) * std::log(bmx);
        default:
            return 0.0;
    }
}

void xerror(const std::string messg, const int nerr, const int level)
{
    std::cerr << nerr << " " << messg << std::endl;
    if (level == 2) std::exit(1);
    return;
}

}; // end namespace quadpack
/**
 * @file quadpack.hpp
 * @brief 
 * @author yoshihiro kawasaki
 * @date 2025/01/14
 * @details 
 * @note
 * @ref
 *      https://www.netlib.org/quadpack/
 *      https://github.com/jacobwilliams/quadpack
*/

#ifndef QUADPACK_HPP_
#define QUADPACK_HPP_

#include <iostream>
#include <string>
#include <limits>
#include <cmath>

/**
* @namespace quadpack_cpp
* @brief
* @details
*/
namespace quadpack_cpp
{

using QUADPACK_CPP_FUNCTION = double(*)(double x, void *user_data);

/**
 * @fn void dqag
 * @brief
 * 
 * @param [in] QUADPACK_CPP_FUNCTION f
 *      function subprogam defining the integrand function f(x).
 * 
 * @param [in] double a
 *      lower limit of integration.
 * 
 * @param [in] double b
 *      upper limit of integration.
 * 
 * @param [in] double epsabs
 *      absolute accoracy requested.
 * 
 * @param [in] double epsrel
 *      relative accuracy requested.
 *      if epsabs <= 0 and epsrel < max(50*rel.mach.acc. , 5.0e-29), 
 *      the routine will end with ier = 6.
 * 
 * @param [in] int key
 *      key for choice of local integration rule
 *      a gauss-kronrod pair is used with
 *       7 - 15 points if key < 2,
 *      10 - 21 points if key = 2,
 *      15 - 31 points if key = 3,
 *      20 - 41 points if key = 4,
 *      25 - 51 points if key = 5,
 *      30 - 61 points if key > 5.
 * 
 * @param [out] double result
 *      approximation to the integral.
 * 
 * @param [out] double abserr
 *      estimate of the modulus of the absolute error,
 *      which should equal or exceed abs(i-result).
 * 
 * @param [out] int neval
 *      number of integrand evaluations.
 * 
 * @param [out] int ier
 *      ier = 0 normal and reliable termination of the
 *              routine. it is assumed that the requested
 *              accuracy has been achieved.
 *      ier > 0 abnormal termination of the routine
 *              the estimates for result and error are
 *              less reliable. it is assumed that the
 *              requested accuracy has not been achieved.
 *      error messages
 *      ier = 1 maximum number of subdivisions allowed
 *              has been achieved. one can allow more
 *              subdivisions by increasing the value of
 *              limit (and taking the according dimension
 *              adjustments into account). however, if
 *              this yield no improvement it is advised
 *              to analyze the integrand in order to
 *              determine the integration difficulaties.
 *              if the position of a local difficulty can
 *              be determined (i.e.singularity,
 *              discontinuity within the interval) one
 *              will probably gain from splitting up the
 *              interval at this point and calling the
 *              integrator on the subranges. if possible,
 *              an appropriate special-purpose integrator
 *              should be used which is designed for
 *              handling the type of difficulty involved.
 *          = 2 the occurrence of roundoff error is
 *              detected, which prevents the requested
 *              tolerance from being achieved.
 *          = 3 extremely bad integrand behaviour occurs
 *              at some points of the integration
 *              interval.
 *          = 6 the input is invalid, because
 *              (epsabs <= 0 and
 *              epsrel < max(50*rel.mach.acc. , 5.0e-29))
 *              or limit < 1 or lenw < limit*4.
 *              result, abserr, neval, last are set
 *              to zero.
 *              except when lenw is invalid, iwork(1),
 *              work(limit*2+1) and work(limit*3+1) are
 *              set to zero, work(1) is set to a and
 *              work(limit+1) to b.
 * 
 * @param [in] int limit
 *      dimensioning parameter for iwork.
 *      limit determines the maximum number of subintervals
 *      in the partition of the given integration interval
 *      (a,b), limit >= 1.
 *      if limit < 1, the routine will end with ier = 6.
 * 
 * @param [in] int lenw
 *      dimensioning parameter for work.
 *      lenw must be at least limit*4.
 *      if lenw < limit*4, the routine will end with
 *      ier = 6.
 * 
 * @param [inout] int last
 *      on return, last equals the number of subintervals
 *      produced in the subdiviosion process, which
 *      determines the number of significant elements
 *      actually in the work arrays.
 * 
 * @param [inout] int *iwork
 *      vector of dimension at least limit, the first k
 *      elements of which contain pointers to the error
 *      estimates over the subintervals, such that
 *      work(limit*3+iwork(1)),... , work(limit*3+iwork(k))
 *      form a decreasing sequence with k = last if
 *      last <= (limit/2+2), and k = limit+1-last otherwise.      
 * 
 * @param [inout] double *work
 * 
 * @param [in] void *data
 * 
 * @sa dqage, xerror
 *      
*/
void dqag(QUADPACK_CPP_FUNCTION f, const double a, const double b, const double epsabs, const double epsrel,
          const int key, double &result, double &abserr, int &neval, int &ier, const int limit, 
          const int lenw, int &last, int *iwork, double *work, void *user_data);

void dqage(QUADPACK_CPP_FUNCTION f, const double a, const double b, const double epsabs, const double epsrel,
           const int key, const int limit, double &result, double &abserr, int &neval, int &ier, 
           double *alist, double *blist, double *rlist, double *elist, int *iord, int &last, void *user_data);

void dqagi(QUADPACK_CPP_FUNCTION f, const double bound, const int inf, const double epsabs, const double epsrel,
           double &result, double &abserr, int &neval, int &ier, const int limit, const int lenw,
           int &last, int *iwork, double *work, void *user_data);

void dqagie(QUADPACK_CPP_FUNCTION f, const double bound, const int inf, const double epsabs, const double epsrel,
            const int limit, double &result, double &abserr, int &neval, int &ier, double *alist,
            double *blist, double *rlist, double *elist, int *iord, int &last, void *user_data);

void dqagp(QUADPACK_CPP_FUNCTION f, const double a, const double b, const int npts2, const double *points,
           const double epsabs, const double epsrel, double &result, double &abserr, int &neval, int &ier,
           const int leniw, const int lenw, int &last, int *iwork, double *work, void *user_data);

void dqagpe(QUADPACK_CPP_FUNCTION f, const double a, const double b, const int npts2, const double *points,
            const double epsabs, const double epsrel, const int limit, double &result, double &abserr, int &neval,
            int &ier, double *alist, double *blist, double *rlist, double *elist, double *pts, int *iord,
            int *level, int *ndin, int &last, void *user_data);

void dqags(QUADPACK_CPP_FUNCTION f, const double a, const double b, const double epsabs, const double epsrel, 
           double &result, double &abserr, int &neval, int &ier, const int limit, const int lenw, int &last, 
           int *iwork, double *work, void *user_data);

void dqagse(QUADPACK_CPP_FUNCTION f, const double a, const double b, const double epsabs, const double epsrel, 
            const int limit, double &result, double &abserr, int &neval, int &ier, double *alist, double *blist,
            double *rlist, double *elist, int *iord, int &last, void *user_data);

void dqawc(QUADPACK_CPP_FUNCTION f, const double a, const double b, const double c, const double epsabs, const double epsrel,
           double &result, double &abserr, int &neval, int &ier, const int limit, const int lenw, int &last, int *iwork,
           double *work, void *user_data);

void dqawce(QUADPACK_CPP_FUNCTION f, const double a, const double b, const double c, const double epsabs, const double epsrel,
            const int limit, double &result, double &abserr, int &neval, int &ier, double *alist, double *blist, double *rlist,
            double *elist, int *iord, int &last, void *user_data);

void dqawf(QUADPACK_CPP_FUNCTION f, const double a, const double omega, const int integr, const double epsabs, double &result,
           double &abserr, int &neval, int &ier, const int limlst, int &lst, const int leniw, const int maxp1, const int lenw,
           int *iwork, double *work, void *user_data);

void dqawfe(QUADPACK_CPP_FUNCTION f, const double a, const double omega, const int integr, const double epsabs, const int limlst,
            const int limit, const int maxp1, double &result, double &abserr, int &neval, int &ier, double *rslst, double *erlst,
            int *ierlst, int &lst, double *alist, double *blist, double *rlist, double *elist, int *iord, int *nnlog, double *chebmo, 
            void *user_data);

void dqawo(QUADPACK_CPP_FUNCTION f, const double a, const double b, const double omega, const int integr, const double epsabs, 
           const double epsrel, double &result, double &abserr, int &neval, int &ier, const int leniw, const int maxp1, const int lenw,
           int &last, int *iwork, double *work, void *user_data);

void dqawoe(QUADPACK_CPP_FUNCTION f, const double a, const double b, const double omega, const int integr, const double epsabs,
            const double epsrel, const int limit, const int icall, const int maxp1, double &result, double &abserr, int  &neval,
            int &ier, int &last, double *alist, double *blist, double *rlist, double *elist, int *iord, int *nnlog, int &momcom,
            double *chebmo, void *user_data);

void dqaws(QUADPACK_CPP_FUNCTION f, const double a, const double b, const double alfa, const double beta, const int integr, const double epsabs,
           const double epsrel, double &result, double &abserr, int &neval, int &ier, const int limit, const int lenw, int &last, int *iwork,
           double *work, void *user_data);

void dqawse(QUADPACK_CPP_FUNCTION f, const double a, const double b, const double alfa, const double beta, const int integr, const double epsabs,
            const double epsrel, const int limit, double &result, double &abserr, int &neval, int &ier, double *alist, double *blist, double *rlist,
            double *elist, int *iord, int &last, void *user_data);

void dqng(QUADPACK_CPP_FUNCTION f, const double a, const double b, const double epsabs, const double epsrel,
          double &result, double &abserr, int &neval, int &ier, void *user_data);
          
};

#endif /* QUADPACK_HPP_ */
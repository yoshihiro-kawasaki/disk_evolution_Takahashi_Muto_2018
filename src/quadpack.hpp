/**
 * @file quadpack.hpp
 * @brief 
 * @author yoshihiro kawasaki
 * @date 2025/01/14
 * @details 
 * @note
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

typedef double (*QUADPACK_CPP_FUNCTION)(const double x, void *data);

typedef double (*QUADPACK_CPP_WEIGHT_FUNCTION)(const double x, 
                                               const double a, 
                                               const double b, 
                                               const double c, 
                                               const double d, 
                                               const int i);
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
void dgtsl(const int n, 
           double *c,
           double *d,
           double *e,
           double *b,
           int &info);

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
void dqag(QUADPACK_CPP_FUNCTION f, 
          const double a, 
          const double b, 
          const double epsabs, 
          const double epsrel,
          const int key, 
          double &result,
          double &abserr, 
          int &neval, 
          int &ier,
          const int limit, 
          const int lenw, 
          int &last, 
          int *iwork,
          double *work,
          void *data);

void dqage(QUADPACK_CPP_FUNCTION f, 
           const double a, 
           const double b, 
           const double epsabs, 
           const double epsrel,
           const int key, 
           const int limit, 
           double &result, 
           double &abserr, 
           int &neval,
           int &ier, 
           double *alist, 
           double *blist,
           double *rlist, 
           double *elist, 
           int *iord, 
           int &last,
           void *data);

void dqagi(QUADPACK_CPP_FUNCTION f,
           const double bound,
           const int inf,
           const double epsabs,
           const double epsrel,
           double &result,
           double &abserr,
           int &neval,
           int &ier,
           const int limit,
           const int lenw,
           int &last,
           int *iwork,
           double *work,
           void *data);

void dqagie(QUADPACK_CPP_FUNCTION f,
            const double bound,
            const int inf,
            const double epsabs,
            const double epsrel,
            const int limit,
            double &result,
            double &abserr,
            int &neval,
            int &ier,
            double *alist,
            double *blist,
            double *rlist,
            double *elist,
            int *iord,
            int &last,
            void *data);

void dqagp(QUADPACK_CPP_FUNCTION f,
           const double a,
           const double b,
           const int npts2,
           const double *points,
           const double epsabs,
           const double epsrel,
           double &result,
           double &abserr,
           int &neval,
           int &ier,
           const int leniw,
           const int lenw,
           int &last,
           int *iwork,
           double *work,
           void *data);

void dqagpe(QUADPACK_CPP_FUNCTION f,
            double a,
            double b,
            const int npts2,
            const double *points,
            const double epsabs,
            const double epsrel,
            const int limit,
            double &result,
            double &abserr,
            int &neval,
            int &ier,
            double *alist,
            double *blist,
            double *rlist,
            double *elist,
            double *pts,
            int *iord,
            int *level,
            int *ndin,
            int &last,
            void *data);

void dqags(QUADPACK_CPP_FUNCTION f, 
           const double a, 
           const double b, 
           const double epsabs, 
           const double epsrel, 
           double &result,
           double &abserr, 
           int &neval, 
           int &ier,
           const int limit, 
           const int lenw, 
           int &last, 
           int *iwork,
           double *work,
           void *data);


void dqagse(QUADPACK_CPP_FUNCTION f, 
            const double a, 
            const double b, 
            const double epsabs, 
            const double epsrel, 
            const int limit, 
            double &result, 
            double &abserr, 
            int &neval,
            int &ier, 
            double *alist, 
            double *blist,
            double *rlist, 
            double *elist, 
            int *iord, 
            int &last,
            void *data);


void dqawc(QUADPACK_CPP_FUNCTION f,
           const double a,
           const double b,
           const double c,
           const double epsabs,
           const double epsrel,
           double &result,
           double &abserr,
           int &neval,
           int &ier,
           const int limit,
           const int lenw,
           int &last,
           int *iwork,
           double *work,
           void *data);

void dqawce(QUADPACK_CPP_FUNCTION f,
            const double a,
            const double b,
            const double c,
            const double epsabs,
            const double epsrel,
            const int limit,
            double &result,
            double &abserr,
            int &neval,
            int &ier,
            double *alist,
            double *blist,
            double *rlist,
            double *elist,
            int *iord,
            int &last,
            void *data);

void dqawf(QUADPACK_CPP_FUNCTION f,
           const double a,
           const double omega,
           const int integr,
           const double epsabs,
           double &result,
           double &abserr,
           int &neval,
           int &ier,
           const int limlst,
           int &lst,
           const int leniw,
           const int maxp1,
           const int lenw,
           int *iwork,
           double *work,
           void *data);

void dqawfe(QUADPACK_CPP_FUNCTION f,
            const double a,
            const double omega,
            const int integr,
            const double epsabs,
            const int limlst,
            const int limit,
            const int maxp1,
            double &result,
            double &abserr,
            int &neval,
            int &ier,
            double *rslst,
            double *erlst,
            int *ierlst,
            int &lst,
            double *alist,
            double *blist,
            double *rlist,
            double *elist,
            int *iord,
            int *nnlog,
            double *chebmo, // **chebmo
            void *data);

void dqawo(QUADPACK_CPP_FUNCTION f,
           const double a,
           const double b,
           const double omega,
           const int integr,
           const double epsabs,
           const double epsrel,
           double &result,
           double &abserr,
           int &neval,
           int &ier,
           const int leniw,
           const int maxp1,
           const int lenw,
           int &last,
           int *iwork,
           double *work,
           void *data);

void dqawoe(QUADPACK_CPP_FUNCTION f,
            const double a,
            const double b,
            const double omega,
            const int integr,
            const double epsabs,
            const double epsrel,
            const int limit,
            const int icall,
            const int maxp1,
            double &result,
            double &abserr,
            int  &neval,
            int &ier,
            int &last,
            double *alist,
            double *blist,
            double *rlist,
            double *elist,
            int *iord,
            int *nnlog,
            int &momcom,
            double *chebmo,
            void * data);

void dqaws(QUADPACK_CPP_FUNCTION f,
           const double a,
           const double b,
           const double alfa,
           const double beta,
           const int integr,
           const double epsabs,
           const double epsrel,
           double &result,
           double &abserr,
           int &neval,
           int &ier,
           const int limit,
           const int lenw,
           int &last,
           int *iwork,
           double *work,
           void *data);

void dqawse(QUADPACK_CPP_FUNCTION f,
            const double a,
            const double b,
            const double alfa,
            const double beta,
            const int integr,
            const double epsabs,
            const double epsrel,
            const int limit,
            double &result,
            double &abserr,
            int &neval,
            int &ier,
            double *alist,
            double *blist,
            double *rlist,
            double *elist,
            int *iord,
            int &last,
            void *data);

void dqc25c(QUADPACK_CPP_FUNCTION f,
            const double a,
            const double b,
            const double c,
            double &result,
            double &abserr,
            int &krul,
            int &neval,
            void *data);

void dqc25f(QUADPACK_CPP_FUNCTION f,
            const double a,
            const double b,
            const double omega,
            const int integr,
            const int nrmom,
            const int maxp1,
            const int ksave,
            double &result,
            double &abserr,
            int &neval,
            double &resabs,
            double &resasc,
            int &momcom,
            double *chebmo, // chebmo
            void *data);

void dqc25s(QUADPACK_CPP_FUNCTION f,
            const double a,
            const double b,
            const double bl,
            const double br,
            const double alfa,
            const double beta,
            const double *ri,
            const double *rj,
            const double *rg,
            const double *rh,
            double &result,
            double &abserr,
            double &resasc,
            const int integr,
            int &nev,
            void *data);

void dqcheb(const double *x,
            double *fval,
            double *cheb12,
            double *cheb24);

void dqelg(int &n,
           double *epstab,
           double &result,
           double &abserr,
           double *res3la,
           int &nres);

void dqk15(QUADPACK_CPP_FUNCTION f,
           const double a,
           const double b,
           double &result,
           double &abserr,
           double &resabs,
           double &resasc,
           void *data);

void dqk15i(QUADPACK_CPP_FUNCTION f,
            const double boun,
            const int inf,
            const double a,
            const double b,
            double &result,
            double &abserr,
            double &resabs,
            double &resasc,
            void *data);

void dqk15w(QUADPACK_CPP_FUNCTION f,
            QUADPACK_CPP_WEIGHT_FUNCTION w,
            const double p1,
            const double p2,
            const double p3,
            const double p4,
            const int kp,
            const double a,
            const double b,
            double &result,
            double &abserr,
            double &resabs,
            double &resasc,
            void *data);

void dqk21(QUADPACK_CPP_FUNCTION f,
           const double a,
           const double b,
           double &result,
           double &abserr,
           double &resabs,
           double &resasc,
           void *data);

void dqk31(QUADPACK_CPP_FUNCTION f,
           const double a,
           const double b,
           double &result,
           double &abserr,
           double &resabs,
           double &resasc,
           void *data);

void dqk41(QUADPACK_CPP_FUNCTION f,
           const double a,
           const double b,
           double &result,
           double &abserr,
           double &resabs,
           double &resasc,
           void *data);

void dqk51(QUADPACK_CPP_FUNCTION f,
           const double a,
           const double b,
           double &result,
           double &abserr,
           double &resabs,
           double &resasc,
           void *data);

void dqk61(QUADPACK_CPP_FUNCTION f,
           const double a,
           const double b,
           double &result,
           double &abserr,
           double &resabs,
           double &resasc,
           void *data);

void dqmomo(const double alfa,
            const double beta,
            double *ri,
            double *rj,
            double *rg,
            double *rh,
            const int integr);

void dqng(QUADPACK_CPP_FUNCTION f,
          const double a,
          const double b,
          const double epsabs,
          const double epsrel,
          double &result,
          double &abserr,
          int &neval,
          int &ier,
          void *data);

void dqpsrt(const int limit,
            const int last,
            int &maxerr,
            double &ermax,
            const double *elist,
            int *iord,
            int &nrmax);

double dqwgtc(const double x,
              const double c,
              const double p2,
              const double p3,
              const double p4,
              const int kp);

double dqwgtf(const double x,
              const double omega,
              const double p2,
              const double p3,
              const double p4,
              const int integr);

double dqwgts(const double x,
              const double a,
              const double b,
              const double alfa,
              const double beta,
              const int integr);

void xerror(const std::string messg, 
           const int nerr, 
           const int level);

};

#endif /* QUADPACK_HPP_ */
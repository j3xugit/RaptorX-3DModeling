#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
#include <string>
#include <vector>
using namespace std;

//----- typedef ------//
typedef short Int2;
typedef long  Int4;
typedef long long Int8;
typedef unsigned short Uint2;
typedef unsigned long  Uint4;
typedef unsigned long long Uint8;

//--- other definition here ----//



/** Number of derivatives of log(x) to carry in gamma-related 
    computations */
#define LOGDERIV_ORDER_MAX	4  
/** Number of derivatives of polygamma(x) to carry in gamma-related 
    computations for non-integral values of x */
#define POLYGAMMA_ORDER_MAX	4

/** value of pi is only used in gamma-related computations */
#define NCBIMATH_PI	3.1415926535897932384626433832795

/** Natural log(2) */
#define NCBIMATH_LN2	0.69314718055994530941723212145818
/** Natural log(PI) */
#define NCBIMATH_LNPI	1.1447298858494001741434273513531

/** size of the next series term that indicates convergence
    in the log and polygamma functions */
#ifndef DBL_EPSILON
#define DBL_EPSILON 2.2204460492503131e-16
#endif

/** HUGE_VAL definition **/
#ifndef HUGE_VAL
#define HUGE_VAL 1e+37
#endif /*HUGE_VAL*/



/** Natural logarithm with shifted input
 *  @param x input operand (x > -1)
 *  @return log(x+1)
 */
extern double BLAST_Log1p (double x);

/** Exponentional with base e 
 *  @param x input operand
 *  @return exp(x) - 1
 */
extern double BLAST_Expm1 (double x);

/** Factorial function
 *  @param n input operand
 *  @return (double)(1 * 2 * 3 * ... * n)
 */
extern double BLAST_Factorial(int n);

/** Logarithm of the factorial 
 *  @param x input operand
 *  @return log(1 * 2 * 3 * ... * x)
 */
extern double BLAST_LnFactorial (double x);

/** log(gamma(n)), integral n 
 *  @param n input operand
 *  @return log(1 * 2 * 3 * ... (n-1))
 */
extern double BLAST_LnGammaInt (int n);

/** Romberg numerical integrator 
 *  @param f Pointer to the function to integrate; the first argument
 *               is the variable to integrate over, the second is a pointer
 *               to a list of additional arguments that f may need
 *  @param fargs Pointer to an array of extra arguments or parameters
 *               needed to compute the function to be integrated. None
 *               of the items in this list may vary over the region
 *               of integration
 *  @param p Left-hand endpoint of the integration interval
 *  @param q Right-hand endpoint of the integration interval
 *           (q is assumed > p)
 *  @param eps The relative error tolerance that indicates convergence
 *  @param epsit The number of consecutive diagonal entries in the 
 *               Romberg array whose relative difference must be less than
 *               eps before convergence is assumed. This is presently 
 *               limited to 1, 2, or 3
 *  @param itmin The minimum number of diagnonal Romberg entries that
 *               will be computed
 *  @return The computed integral of f() between p and q
 */
extern double BLAST_RombergIntegrate (double (*f) (double, void*), 
                               void* fargs, double p, double q, 
                               double eps, int epsit, int itmin);

/** Greatest common divisor 
 *  @param a First operand (any integer)
 *  @param b Second operand (any integer)
 *  @return The largest integer that evenly divides a and b
 */
extern int BLAST_Gcd (int a, int b);

/** Divide 3 numbers by their greatest common divisor
 * @param a First integer [in] [out]
 * @param b Second integer [in] [out]
 * @param c Third integer [in] [out]
 * @return The greatest common divisor
 */
extern int BLAST_Gdb3(int* a, int* b, int* c);

/** Nearest integer 
 *  @param x Input to round (rounded value must be representable
 *           as a 32-bit signed integer)
 *  @return floor(x + 0.5);
 */
extern long BLAST_Nint (double x);

/** Integral power of x 
 * @param x floating-point base of the exponential
 * @param n (integer) exponent
 * @return x multiplied by itself n times
 */
extern double BLAST_Powi (double x, int n);

/** The error function of x: the integral from 0 to x of e(-t*t) dt,
 *  scaled by 2/sqrt(pi) to fall within the range (-1,1). */
extern double BLAST_Erf (double x);

/** The complementary error function of x: 1 - erf(x), but calculated
 *  more accurately for large x (where erf(x) approaches unity). */
extern double BLAST_ErfC (double x);


//------ other utilities please put here --------//

//-> frequency normalize
extern void BLAST_FreqNormalize (double *input,int size,double factor=1.0);


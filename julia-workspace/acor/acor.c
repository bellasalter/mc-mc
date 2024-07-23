/*
 * This is a very slightly re-factored version of code by Jonathan Goodman
 * For the original code, see:
 *  http://www.math.nyu.edu/faculty/goodman/software/acor/index.html
 * Ported for use in Python with permission from the original author.
 */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include "acor.h"

#define WINMULT 5
#define MINFAC  1

int main() {
    double data[200] = {-6.364792,  -2.405985,  -3.060521,  -5.867774,  -2.761165,  3.763935,  2.274887,  -1.191277,  0.036821,  -2.356774,  -2.210410,  -3.148326,  0.832340,  2.219213,  -1.216001,  1.252149,  2.510949,  3.405035,  1.535183,  1.608925,  0.589727,  1.109422,  2.902956,  3.077903,  1.905187,  3.052903,  -0.662005,  3.009243,  6.677305,  1.745914,  3.535656,  3.447538,  -1.672488,  -2.910748,  -2.391768,  -2.194846,  2.161494,  -0.239168,  -2.427846,  1.974192,  -0.156623,  -0.053817,  3.101896,  -1.102728,  -1.733144,  1.695124,  -0.271328,  -0.494322,  5.255137,  4.295133,  4.319540,  3.570357,  1.019935,  2.463013,  1.647400,  3.041588,  1.286103,  -0.092964,  1.357805,  2.499508,  3.057822,  2.375573,  1.257744,  1.496246,  2.637574,  -0.474342,  2.831710,  5.401437,  3.464833,  4.741225,  1.565969,  1.833384,  3.424827,  -0.908102,  3.167409,  -2.203770,  -0.185427,  4.252108,  5.815527,  2.673012,  2.747388,  1.576633,  1.228272,  0.221400,  0.297982,  1.241262,  0.901301,  0.899656,  1.134218,  0.772028,  -1.457878,  -0.887193,  -0.775462,  -2.399401,  -2.850860,  -0.981666,  0.108077,  -1.338109,  0.075879,  -1.735573,  0.251565,  1.146866,  2.549492,  0.515489,  1.023993,  -0.587518,  1.320450,  -0.453406,  -0.316391,  1.700082,  -0.729463,  0.233800,  0.700348,  -1.038779,  0.961980,  2.433019,  -1.930679,  -3.455162,  -1.261099,  -2.931517,  -4.079239,  -0.404268,  2.172780,  -1.430507,  -2.775100,  0.335369,  0.742811,  2.361869,  1.478141,  2.271457,  -1.124313,  -1.882459,  -1.166330,  -2.553729,  -2.070031,  -0.319841,  -1.128319,  -0.492066,  1.750031,  -3.842028,  -0.441428,  0.318812,  3.401601,  6.805297,  3.675329,  2.442296,  2.704897,  1.519433,  -0.475357,  1.854670,  6.958490,  2.781858,  3.221526,  3.353330,  -2.505073,  -0.015047,  -0.653368,  -2.690398,  -3.157167,  -1.487320,  -0.775718,  1.003947,  2.921457,  2.266223,  2.395964,  -1.103396,  -0.227373,  0.397888,  -0.316360,  0.516011,  -0.311217,  3.579063,  1.896212,  3.149646,  -1.511802,  -1.028559,  0.107889,  6.141254,  5.519213,  3.176690,  2.844158,  5.458626,  4.122481,  3.749172,  3.514198,  2.028676,  0.897435,  1.012948,  -0.905149,  -0.703770,  0.754011,  -1.841134,  -4.829756,  -0.769792,  1.979243,  1.488422,  -2.862107,  -4.815402,  -4.644148,  -4.492501};
    double *mn = 0;
    double *sa = 0;
    double *tu = 0;
    int L = 199;
    int max_lag = 199;
    int res = acor(mn, sa, tu, data, L, max_lag);
    printf("%d\n", res);
}

int acor(double *mean, double *sigma, double *tau, double X[], int L, int maxlag)
{
    int i, s;
    double D;
    double *C;
    int iMax = L - maxlag;

    /* Declare variables for the pairwise analysis */
    int Lh = L/2, j1 = 0, j2 = 1;
    double newMean;

    /* Fail if the chain is too short */
    if (L < MINFAC*maxlag)
        return 1;

    /* Allocate memory for the auto-covariance function */
    if (!(C = (double*)malloc((maxlag+1)*sizeof(double))))
        return -1;

    /* compute the mean of X and normalize the data */
    *mean = 0.;
    for (i = 0; i < L; i++)
        *mean += X[i];
    *mean = *mean/L;
    for (i = 0; i <  L; i++ )
        X[i] -= *mean;

    /* Compute the auto-covariance function */
    for (s = 0; s <= maxlag; s++)
        C[s] = 0.;
    for (i = 0; i < iMax; i++)
        for (s = 0; s <= maxlag; s++)
            C[s] += X[i]*X[i+s];
    for (s = 0; s <= maxlag; s++)
        C[s] = C[s]/iMax; /* Normalize */

    /* Calculate the "diffusion coefficient" as the sum of the autocovariances */
    D = C[0];
    for (s = 1; s <= maxlag; s++ )
        D += 2*C[s];

    /* Calculate the standard error, if D were the complete sum. */
    /* FIXME: why is D sometimes negative?? */
    if ( D < 0 ) {
        return 2;
    }
    *sigma = sqrt( D / L );
    /* A provisional estimate, since D is only part of the complete sum. */
    *tau   = D / C[0];

    if ( *tau*WINMULT < maxlag ) {
        free(C);
        return 0; /* Stop if the D sum includes the given multiple of tau. */
    }

    /* If the provisional tau is so large that we don't think tau is accurate,
     * apply the acor procedure to the pairwise sums of X */
    for (i = 0; i < Lh; i++) {
        X[i] = X[j1] + X[j2];
        j1  += 2;
        j2  += 2;
    }

    acor( &newMean, sigma, tau, X, Lh, maxlag);
    D      = .25*(*sigma) * (*sigma) * L;
    *tau   = D/C[0];
    *sigma = sqrt( D/L );

    free(C);
    return 0;
}


int acor_fn(double *mean, double *fn, double X[], int L, int maxt)
{
    int t, n;

    /* compute the mean of X */
    *mean = 0.;
    for (n = 0; n < L; n++)
        *mean += X[n];
    *mean = *mean / L;

    /* Normalize */
    for (n = 0; n <  L; n++ )
        X[n] -= *mean;

    /* Compute the auto-covariance function */
    for (t = 0; t < maxt; t++)
        fn[t] = 0.;

    for (t = 0; t < maxt; t++)
        for (n = 0; n < L - t; n++)
            fn[t] += X[n] * X[n + t];

    for (t = 0; t < maxt; t++)
        fn[t] /= L - t;

    return 0;
}

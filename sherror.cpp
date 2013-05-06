// Copyright (c) 2013 Robert Kooima
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.

#include <getopt.h>

#include <algorithm>
#include <cassert>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <sys/time.h>

#include <image.h>

#include "sh.hpp"

//------------------------------------------------------------------------------

static int usage(const char *exe)
{
    fprintf(stderr,
            "%s [-n n]\n"
            "\t-n ... Harmonic degree (1)\n\n", exe);

    return -1;
}

template <typename real> real rmserr(Flm<real>& F)
{
    const int n = F.n;

    real t = 0;

    for     (int l =  0; l <  n; l++)
        for (int m = -l; m <= l; m++)

            t += (F(l, m, 0) - 1.0) * (F(l, m, 0) - 1.0);

    return sqrt(t / (n * n));
}

template <typename real> real maxerr(Flm<real>& F)
{
    const int n = F.n;

    real e = 0;

    for     (int l =  0; l <  n; l++)
        for (int m = -l; m <= l; m++)

            if (e < fabs(F(l, m, 0) - 1.0))
                e = fabs(F(l, m, 0) - 1.0);

    return e;
}

template <typename real> void test(int n)
{
    sht<real> T(n, 1);

    struct timeval t0;
    struct timeval t1;

    // Construct a white noise input.

    for     (int l =  0; l <  n; l++)
        for (int m = -l; m <= l; m++)
            T.F(l, m, 0) = 1.0;

    // Instance the transformer and perform a timed round trip.

    gettimeofday(&t0, 0);
    {
        T.syn();
        T.ana();
    }
    gettimeofday(&t1, 0);

    float d = (t1.tv_sec  - t0.tv_sec) +
              (t1.tv_usec - t0.tv_usec) / 1000000.0f;

    // Compute and report error bounds.

    float r = rmserr<real>(T.F);
    float m = maxerr<real>(T.F);

    printf("%4d %8.3f %20.16e %8.3f %20.16e %8.3f\n",
            n, d, r, -log2f(r), m, -log2f(m));
}

int main(int argc, char **argv)
{
    // Parse the options.

    int t = 'L';
    int n = 1;
    int o;

    while ((o = getopt(argc, argv, "n:FDLWQ")) != -1)
        switch (o)
        {
            case 'n': n = strtol(optarg, 0, 0); break;
            case 'F': t = o;                    break;
            case 'D': t = o;                    break;
            case 'L': t = o;                    break;
            case 'W': t = o;                    break;
            case 'Q': t = o;                    break;
            default: return usage(argv[0]);
        }

    // Confirm a reasonable configuration and run.

    if (n > 0)
    {
        switch (t)
        {
            case 'F': test<      float>(n); break;
            case 'D': test<     double>(n); break;
            case 'L': test<long double>(n); break;
#if 0
            case 'W': test< __float80 >(n); break;
            case 'Q': test< __float128>(n); break;
#endif
        }
    }
    else usage(argv[0]);

    return 0;
}

//------------------------------------------------------------------------------

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

#include <image.h>

#include "sh.hpp"

//------------------------------------------------------------------------------

template <typename real> void grace_cathedral(sht<real>& T)
{
    T.F(0,  0, 0) =  0.79; T.F(0,  0, 1) =  0.44; T.F(0,  0, 2) =  0.54;
    T.F(1, -1, 0) =  0.39; T.F(1, -1, 1) =  0.35; T.F(1, -1, 2) =  0.60;
    T.F(1,  0, 0) = -0.34; T.F(1,  0, 1) = -0.18; T.F(1,  0, 2) = -0.27;
    T.F(1,  1, 0) = -0.29; T.F(1,  1, 1) = -0.06; T.F(1,  1, 2) =  0.01;
    T.F(2, -2, 0) = -0.11; T.F(2, -2, 1) = -0.05; T.F(2, -2, 2) = -0.12;
    T.F(2, -1, 0) = -0.26; T.F(2, -1, 1) = -0.22; T.F(2, -1, 2) = -0.47;
    T.F(2,  0, 0) = -0.16; T.F(2,  0, 1) = -0.09; T.F(2,  0, 2) = -0.15;
    T.F(2,  1, 0) =  0.56; T.F(2,  1, 1) =  0.21; T.F(2,  1, 2) =  0.14;
    T.F(2,  2, 0) =  0.21; T.F(2,  2, 1) = -0.05; T.F(2,  2, 2) = -0.30;
}

template <typename real> void eucalyptus_grove(sht<real>& T)
{
    T.F(0,  0, 0) =  0.38; T.F(0,  0, 1) =  0.43; T.F(0,  0, 2) =  0.45;
    T.F(1, -1, 0) =  0.29; T.F(1, -1, 1) =  0.36; T.F(1, -1, 2) =  0.41;
    T.F(1,  0, 0) =  0.04; T.F(1,  0, 1) =  0.03; T.F(1,  0, 2) =  0.01;
    T.F(1,  1, 0) = -0.10; T.F(1,  1, 1) = -0.10; T.F(1,  1, 2) = -0.09;
    T.F(2, -2, 0) = -0.06; T.F(2, -2, 1) = -0.06; T.F(2, -2, 2) = -0.04;
    T.F(2, -1, 0) =  0.01; T.F(2, -1, 1) = -0.01; T.F(2, -1, 2) = -0.05;
    T.F(2,  0, 0) = -0.09; T.F(2,  0, 1) = -0.13; T.F(2,  0, 2) = -0.15;
    T.F(2,  1, 0) = -0.06; T.F(2,  1, 1) = -0.05; T.F(2,  1, 2) = -0.04;
    T.F(2,  2, 0) =  0.02; T.F(2,  2, 1) = -0.00; T.F(2,  2, 2) = -0.05;
}

template <typename real> void st_peters_basilica(sht<real>& T)
{
    T.F(0,  0, 0) =  0.36; T.F(0,  0, 1) =  0.26; T.F(0,  0, 2) =  0.23;
    T.F(1, -1, 0) =  0.18; T.F(1, -1, 1) =  0.14; T.F(1, -1, 2) =  0.13;
    T.F(1,  0, 0) = -0.02; T.F(1,  0, 1) = -0.01; T.F(1,  0, 2) = -0.00;
    T.F(1,  1, 0) =  0.03; T.F(1,  1, 1) =  0.02; T.F(1,  1, 2) =  0.01;
    T.F(2, -2, 0) =  0.02; T.F(2, -2, 1) =  0.01; T.F(2, -2, 2) =  0.00;
    T.F(2, -1, 0) = -0.05; T.F(2, -1, 1) = -0.03; T.F(2, -1, 2) = -0.01;
    T.F(2,  0, 0) = -0.09; T.F(2,  0, 1) = -0.08; T.F(2,  0, 2) = -0.07;
    T.F(2,  1, 0) =  0.01; T.F(2,  1, 1) =  0.00; T.F(2,  1, 2) =  0.00;
    T.F(2,  2, 0) = -0.08; T.F(2,  2, 1) = -0.06; T.F(2,  2, 2) =  0.00;
}

//------------------------------------------------------------------------------

static int usage(const char *exe)
{
    fprintf(stderr,
            "%s [-o output] [-b b] [-l l] [-m m] [-w w] [-h h]\n"
            "\t-o ... Output file name (out.tif)\n"
            "\t-b ... Output depth     (4)\n"
            "\t-l ... Harmonic degree  (0)\n"
            "\t-m ... Harmonic order   (0)\n"
            "\t-w ... Synthesis width  (64)\n"
            "\t-h ... Synthesis height (32)\n\n", exe);

    return -1;
}

int main(int argc, char **argv)
{
    // Set default options.

    const char *out = "out.tif";
    int         w   = 64;
    int         h   = 32;
    int         b   = 4;
    int         l   = 0;
    int         m   = 0;

    // Parse the options.

    int o;

    while ((o = getopt(argc, argv, "b:h:l:m:o:w:")) != -1)
        switch (o)
        {
            case 'b': b = strtol(optarg, 0, 0); break;
            case 'h': h = strtol(optarg, 0, 0); break;
            case 'l': l = strtol(optarg, 0, 0); break;
            case 'm': m = strtol(optarg, 0, 0); break;
            case 'w': w = strtol(optarg, 0, 0); break;
            case 'o': out = optarg;             break;

            default: return usage(argv[0]);
        }

    // Confirm a reasonable output request.

    if (w > 0 && h > 0 && b > 0 && m >= -l && m <= l)
    {
        float *src = 0;
        float *dst = 0;
        int n = l + 1;

        // Allocate source and destination buffers.

        if ((src = (float *) calloc(n * n * 3, sizeof (float))) &&
            (dst = (float *) calloc(w * h * 3, sizeof (float))))
        {
            // Construct the input.

            Flm<float> F(n, 3);

            F(l, m, 0) = -1;
            F(l, m, 1) =  1;

            F.get(src);

            // Instance the transformer and do the work.

            sht<long double> T(n, w, h, 3);

            T.F.set(src);
            T.syn();
            T.S.get(dst);

            image_write_float(out, w, h, 3, b, dst);
        }

        free(dst);
        free(src);
    }
    else usage(argv[0]);

    return 0;
}

//------------------------------------------------------------------------------

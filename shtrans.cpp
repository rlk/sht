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

static int usage(const char *exe)
{
    fprintf(stderr,
            "%s    [-o output] [-n n] [-b b] input\n"
            "\t-o ... Output file name (out.tif)\n"
            "\t-b ... Output depth     (4)\n"
            "\t-n ... Analysis degree  (3)\n\n"

            "%s -i [-o output] [-w w] [-h h] [-b b] input\n"
            "\t-o ... Output file name (out.tif)\n"
            "\t-b ... Output depth     (4)\n"
            "\t-w ... Synthesis width  (64)\n"
            "\t-h ... Synthesis height (32)\n\n",

            exe, exe);

    return -1;
}

template <typename real> void syn(const char *in, const char *out, int w, int h, int b)
{
    float *src = 0;
    float *dst = 0;
    int    W;
    int    H;
    int    B;
    int    C;

    // Load the image and cast it to floating point.

    if ((src = image_read_float(in, &W, &H, &C, &B)))
    {
        // The input must be square and the requested output size reasonable.

        if (W == H && w > 0 && h > 0 && b > 0)
        {
            // Allocate a destination buffer.

            if ((dst = (float *) calloc(w * h * C, sizeof (float))))
            {
                // Instance the transformer and do the work.

                sht<real> T(H, w, h, C);

                T.F.set(src);
                T.syn();
                T.S.get(dst);

                image_write_float(out, w, h, C, b, dst);
            }
        }
        else fprintf(stderr, "Bad parameters\n");
    }
    else fprintf(stderr, "Failed to load %s\n", in);

    free(dst);
    free(src);
}

template <typename real> void ana(const char *in, const char *out, int n, int b)
{
    float *src = 0;
    float *dst = 0;
    float *tmp = 0;
    int    W;
    int    H;
    int    B;
    int    C;

    // Load the image and cast it to floating point.

    if ((src = image_read_float(in, &W, &H, &C, &B)))
    {
        if (n > 0)
        {
            // Determine the smallest image meeting Nyquist's criterion.

            int w = std::max(n * 2, W);
            int h = std::max(n * 2, H);

            // If the input is too small, upsample it up before analysis.

            if (w != W || h != H)
            {
                tmp = src;
                src = image_scale_float(w, h, W, H, C, tmp);
            }

            // Allocate a destination buffer and do the work.

            if ((dst = (float *) calloc(n * n * C, sizeof (float))))
            {
                sht<real> T(n, w, h, C);

                T.S.set(src);
                T.ana();
                T.F.get(dst);

                image_write_float(out, n, n, C, b, dst);
            }
        }
        else fprintf(stderr, "Bad parameters\n");
    }
    else fprintf(stderr, "Failed to load %s\n", in);

    free(tmp);
    free(dst);
    free(src);
}

int main(int argc, char **argv)
{
    // Set default options.

    const char *out = "out.tif";
    bool        inv = false;
    int         w   = 64;
    int         h   = 32;
    int         b   = 4;
    int         n   = 3;
    int         t   = 'L';

    // Parse the options.

    int o;

    while ((o = getopt(argc, argv, "b:h:n:o:w:iFDL")) != -1)
        switch (o)
        {
            case 'b': b = strtol(optarg, 0, 0); break;
            case 'h': h = strtol(optarg, 0, 0); break;
            case 'n': n = strtol(optarg, 0, 0); break;
            case 'w': w = strtol(optarg, 0, 0); break;
            case 'i': inv = true;               break;
            case 'o': out = optarg;             break;
            case 'F': t = o;                    break;
            case 'D': t = o;                    break;
            case 'L': t = o;                    break;

            default: return usage(argv[0]);
        }

    // Run the job.

    if (argc > optind)
    {
        if (inv)
            switch (t)
            {
                case 'F':  syn<      float>(argv[optind], out, w, h, b); break;
                case 'D':  syn<     double>(argv[optind], out, w, h, b); break;
                case 'L':  syn<long double>(argv[optind], out, w, h, b); break;
            }
        else
            switch (t)
            {
                case 'F':  ana<      float>(argv[optind], out,    n, b); break;
                case 'D':  ana<     double>(argv[optind], out,    n, b); break;
                case 'L':  ana<long double>(argv[optind], out,    n, b); break;
            }
    }
    else return usage(argv[0]);

    return 0;
}

//------------------------------------------------------------------------------

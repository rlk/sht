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

static inline bool ispow2(int n)
{
    return n && !(n & (n - 1));
}

static inline int pow2(int n)
{
    n--;
    n |= n >> 1;
    n |= n >> 2;
    n |= n >> 4;
    n |= n >> 8;
    n |= n >> 16;
    n++;

    return n;
}

//------------------------------------------------------------------------------

static int usage(const char *exe)
{
    fprintf(stderr,
            "%s -i [-FDL] [-o output] [-n n] [-b b] input\n"
            "\t-F ..... Float precision\n"
            "\t-D ..... Double precision\n"
            "\t-L ..... Long double precision\n"
            "\t-o ..... Output file name    (out.tif)\n"
            "\t-b ..... Output depth        (4)\n\n"

            "%s    [-FDL] [-o output] [-n n] [-b b] [-hlg w] [-d] input\n"
            "\t-h width ... Hanning  filter (n)\n"
            "\t-l width ... Lanczos  filter (n)\n"
            "\t-g width ... Gaussian filter (n)\n"
            "\t-d ......... Diffuse convolution\n\n"

            , exe, exe);

    return -1;
}

template <typename real> void syn(const char *in, const char *out, int bb)
{
    float *src = 0;
    float *dst = 0;
    int    w;
    int    h;
    int    b;
    int    c;

    // Load the image and cast it to floating point.

    if ((src = image_read_float(in, &w, &h, &c, &b)))
    {
        // The input must be square.

        if (w == h)
        {
            int n = pow2(w);

            // Allocate a destination buffer.

            if ((dst = (float *) calloc(4 * n * n * c, sizeof (float))))
            {
                // Instance the transformer and do the work.

                sht<real> T(n, c);

                T.F.set(src, w);
                T.syn();
                T.S.get(dst, 2 * n);

                image_write_float(out, 2 * n, 2 * n, c, bb, dst);
            }
        }
        else fprintf(stderr, "Image must be square\n");
    }
    else fprintf(stderr, "Failed to load %s\n", in);

    free(dst);
    free(src);
}

template <typename real> void ana(const char *in,
                                  const char *out, int bb, int fo, int fn)
{
    float *src = 0;
    float *dst = 0;
    int    w;
    int    h;
    int    b;
    int    c;

    // Load the image and cast it to floating point.

    if ((src = image_read_float(in, &w, &h, &c, &b)))
    {
        if (ispow2(w) && ispow2(h))
        {
            int n = w / 2;

            // Allocate a destination buffer and do the work.

            if ((dst = (float *) calloc(n * n * c, sizeof (float))))
            {
                sht<real> T(n, c);

                T.S.set(src, h);
                T.ana();

                switch (fo)
                {
                    case 'h': T.F.hanning(fn ? fn : n); break;
                    case 'l': T.F.lanczos(fn ? fn : n); break;
                    case 'g': T.F.gauss  (fn ? fn : n); break;
                    case 'd': T.F.diffuse();            break;
                }
                T.F.get(dst, n);

                image_write_float(out, n, n, c, bb, dst);
            }
        }
        else fprintf(stderr, "Input must be power-of-two\n");
    }
    else fprintf(stderr, "Failed to load %s\n", in);

    free(dst);
    free(src);
}

int main(int argc, char **argv)
{
    // Set default options.

    const char *o = "out.tif";
    bool        i = false;
    int         b = 4;
    int         w = 0;
    int         t = 'L';
    int         f = ' ';

    // Parse the options.

    int c;

    while ((c = getopt(argc, argv, "b:o:h:l:g:diFDL")) != -1)
        switch (c)
        {
            case 'o': o = optarg;                      break;
            case 'i': i = true;                        break;
            case 'F':                           t = c; break;
            case 'D':                           t = c; break;
            case 'L':                           t = c; break;
            case 'b': b = strtol(optarg, 0, 0);        break;
            case 'h': w = strtol(optarg, 0, 0); f = c; break;
            case 'l': w = strtol(optarg, 0, 0); f = c; break;
            case 'g': w = strtol(optarg, 0, 0); f = c; break;
            case 'd':                           f = c; break;

            default: return usage(argv[0]);
        }

    // Run the job.

    if (argc > optind)
    {
        if (i)
            switch (t)
            {
                case 'F':  syn<      float>(argv[optind], o, b); break;
                case 'D':  syn<     double>(argv[optind], o, b); break;
                case 'L':  syn<long double>(argv[optind], o, b); break;
            }
        else
            switch (t)
            {
                case 'F':  ana<      float>(argv[optind], o, b, f, w); break;
                case 'D':  ana<     double>(argv[optind], o, b, f, w); break;
                case 'L':  ana<long double>(argv[optind], o, b, f, w); break;
            }
    }
    else return usage(argv[0]);

    return 0;
}

//------------------------------------------------------------------------------

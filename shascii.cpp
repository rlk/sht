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
            "%s [-o output] [-n n] input\n"
            "\t-o ... Output file name (out.tif)\n"
            "\t-n ... Degree           (256)\n\n", exe);

    return -1;
}

int main(int argc, char **argv)
{
    // Set default options.

    const char *out = "out.tif";
    int         n   = 256;

    // Parse the options.

    int o;

    while ((o = getopt(argc, argv, "n:o:")) != -1)
        switch (o)
        {
            case 'n': n = strtol(optarg, 0, 0); break;
            case 'o': out = optarg;             break;

            default: return usage(argv[0]);
        }

    // Confirm a reasonable output request.

    if (n > 0 && argc > optind)
    {
        // Allocate a destination buffer.

        if (float *dst = (float *) calloc(n * n * 1, sizeof (float)))
        {
            // Instance the transformer and parse the input.

            sht<long double> T(n, 1);

            if (FILE *fp = fopen(argv[optind], "r"))
            {
                long double c;
                long double s;
                int         l = 0;
                int         m = 0;
                char buf[256];

                while (fgets(buf, 256, fp) && l < n)
                {
                    if (fscanf(fp, "%d %d %Lf %Lf\n", &l, &m, &c, &s))
                    {
                        T.F(l,  m, 0) = c;
                        T.F(l, -m, 0) = s;
                    }
                }
                fclose(fp);
            }

            // Store the output.

            T.F.get(dst, n);

            image_write_float(out, n, n, 1, 4, dst);

            free(dst);
        }
    }
    else usage(argv[0]);

    return 0;
}

//------------------------------------------------------------------------------

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

#ifndef SH_HPP
#define SH_HPP

#include <cmath>

//------------------------------------------------------------------------------

// Plm is a real buffer indexed by degree l and order m where m is not negative.
// This is the structure of the associated Legendre functions.

template <typename real>
    class Plm
    {
    public:
        Plm(int n) : n(n)
        {
            P = new real[(n + 1) * n / 2];
        }
        real& operator()(int l, int m)
        {
//          assert(0 <= l && l <  n);
//          assert(0 <= m && m <= l);
            return P[(l + 1) * l / 2 + m];
        }
       ~Plm()
        {
            delete [] P;
        }
        int n;

    private:
        real *P;
    };

// Flm is a real buffer with c channels indexed by degree l and order m.
// This is the structure of a set of spherical harmonic coefficients.

template <typename real>
    class Flm
    {
    public:
        Flm(int n, int c) : n(n), c(c)
        {
            F = new real[n * n * c];
        }
        real& operator()(int l, int m, int k)
        {
            assert( 0 <= l && l <  n);
            assert(-l <= m && m <= l);
            assert( 0 <= k && k <  c);
            return F[((l + 1) * l + m) * c + k];
        }
        void set(const float *);
        void get(      float *);
       ~Flm()
        {
            delete [] F;
        }
        int n;
        int c;

    private:
        real *F;
    };

// Sij is a real buffer with c channels indexed by width w and height h.

template <typename real>
    class Sij
    {
    public:
        Sij(int w, int h, int c) : w(w), h(h), c(c)
        {
            S = new real[w * h * c];
        }
        real& operator()(int i, int j, int k)
        {
            assert(0 <= i && i < h);
            assert(0 <= j && j < w);
            assert(0 <= k && k < c);
            return S[(i * w + j) * c + k];
        }
        void set(const float *);
        void get(      float *);
       ~Sij()
        {
            delete [] S;
        }
        int w;
        int h;
        int c;

    private:
        real *S;
    };

// Class sht is the complete spherical harmonic transform implementation.

template <typename real>
    class sht
    {
    public:

        sht(int, int, int, int);

        void syn();
        void ana();

        Plm<real> P;
        Flm<real> F;
        Sij<real> S;

    private:

        real a(int)      const;
        real b(int, int) const;
        real c(int, int) const;

        void alf(real);

        real theta (int) const;
        real lambda(int) const;
        real delta (int) const;

        void anad(int, int, int, int, real);
        void synd(int, int, int, int);
    };

//------------------------------------------------------------------------------
// Construct a transformer. Allocate storage for ALFs as well as the frequency
// domain and spatial domain representations.

template <typename real> sht<real>::sht(int n, int w, int h, int c)
    : P(n), F(n, c), S(w, h, c)
{
    // Yep non-inlined empty constructor. Woo!
}

//------------------------------------------------------------------------------
// Compute orthonormalized ALF coefficients using the m-varying recurrence.

template <typename real> inline real sht<real>::a(int m) const
{
    return sqrt((2.0 * m + 1.0) /
                (2.0 * m      ));
}

template <typename real> inline real sht<real>::b(int l, int m) const
{
    return (2.0 * m + 2.0) / sqrt((l + m + 1.0) * (l - m));
}

template <typename real> inline real sht<real>::c(int l, int m) const
{
    return sqrt(((l + m + 2.0) * (l - m - 1.0)) /
                ((l + m + 1.0) * (l - m      )));
}

template <typename real> void sht<real>::alf(real x)
{
    const real y = sqrt(1.0 - x * x), t = x / y;
    int l;
    int m;

    // Begin with the recursive base case.

    P(0, 0) = 0.282094791773878143474040L;

    // Perform the diagonal recurrence.

    for (m = 1; m < P.n; m++)
        P(m, m) = y * a(m) * P(m - 1, m - 1);

    // Perform the triangular recurrence.

    for (l = 1; l < P.n; l++)
    {
        // The first element off the diagonal is simplified by a zero term.

        P(l, l - 1) = t * b(l, l - 1) * P(l, l);

        // Proceed from the diagonal toward column zero.

        for (m = l - 2; m >= 0; m--)
            P(l, m) = t * b(l, m) * P(l, m + 1)
                        - c(l, m) * P(l, m + 2);
    }
}

//------------------------------------------------------------------------------

// Compute spherical coordinates for a given image location.

template <typename real> inline real sht<real>::theta(int i) const
{
    return       M_PI * (real(i) + 0.5) / real(S.h);
}

template <typename real> inline real sht<real>::lambda(int j) const
{
    return 2.0 * M_PI * (real(j) + 0.5) / real(S.w);
}

// Calculate the Chebyshev quadrature weight for image row i.

template <typename real> inline real sht<real>::delta(int i) const
{
    const int n = F.n;
    real T = theta(i);
    real d = 0.0;

    for (int l = 0; l < n; ++l)
    {
        int x = l * 2 + 1;
        d += sin(T * x) / x;
    }

    return 2.0 * M_PI * d * sin(T) / (n * n);
}

//------------------------------------------------------------------------------

// Integrate one spatial sample Sij from one frequency sample Flm.

template <typename real> void sht<real>::synd(int i, int j, int l, int m)
{
    const real L = lambda(j);
    const real p = P(l, m);

    const real c = M_SQRT2 * cos(L * m);
    const real s = M_SQRT2 * sin(L * m);

    for (int k = 0; k < F.c; k++)
        if (m)
        {
            S(i, j, k) += F(l,  m, k) * p * c
                       +  F(l, -m, k) * p * s;
        }
        else
            S(i, j, k) += F(l,  0, k) * p;
}

// Integrate one frequency sample Flm from one spatial sample Sij.

template <typename real> void sht<real>::anad(int i, int j, int l, int m, real d)
{
    const real L = lambda(j);
    const real p = P(l, m) * d;

    const real c = M_SQRT2 * cos(L * m);
    const real s = M_SQRT2 * sin(L * m);

    for (int k = 0; k < F.c; k++)
        if (m)
        {
            F(l,  m, k) += S(i, j, k) * p * c;
            F(l, -m, k) += S(i, j, k) * p * s;
        }
        else
            F(l,  0, k) += S(i, j, k) * p;
}

// Perform a spherical harmonic synthesis.

template <typename real> void sht<real>::syn()
{
    const int n = F.n;
    const int w = S.w;
    const int h = S.h;

    int i, j, l, m;

    S.set(0);

    for (i = 0; i < h; i++)
    {
        alf(cos(theta(i)));

        #pragma omp parallel for private(l, m) schedule(static)
        for         (j = 0; j <  w; j++)
            for     (l = 0; l <  n; l++)
                for (m = 0; m <= l; m++)
                    synd(i, j, l, m);
    }
}

// Perform a spherical harmonic analysis.

template <typename real> void sht<real>::ana()
{
    const int n = F.n;
    const int w = S.w;
    const int h = S.h;

    int i, j, l, m;

    F.set(0);

    for (i = 0; i < h; i++)
    {
        const real d = delta(i);

        alf(cos(theta(i)));

        #pragma omp parallel for private(m, j) schedule(dynamic)
        for         (l = 0; l <  n; l++)
            for     (m = 0; m <= l; m++)
                for (j = 0; j <  w; j++)
                    anad(i, j, l, m, d);
    }
}

//------------------------------------------------------------------------------

// Map and cast the contents of the n * n * c floating point raster p into F.

template <typename real> void Flm<real>::set(const float *p)
{
    if (p == 0)
        memset(F, 0, n * n * c * sizeof (real));
    else
        for         (int l =  0; l <  n; l++)
            for     (int m = -l; m <= l; m++)
                for (int k =  0; k <  c; k++)
                {
                    int i = m < 0 ? l : l - m;
                    int j = m > 0 ? l : l + m;

                    F[((l + 1) * l + m) * c + k] = real(p[(i * n + j) * c + k]);
                }
}

// Map and cast the contents of F into the n * n * c floating point raster p.

template <typename real> void Flm<real>::get(float *p)
{
    for         (int l =  0; l <  n; l++)
        for     (int m = -l; m <= l; m++)
            for (int k =  0; k <  c; k++)
            {
                int i = m < 0 ? l : l - m;
                int j = m > 0 ? l : l + m;

                p[(i * n + j) * c + k] = float(F[((l + 1) * l + m) * c + k]);
            }
}

// Map and cast the contents of the w * h * c floating point raster p into S.

template <typename real> void Sij<real>::set(const float *p)
{
    if (p == 0)
         memset(S, 0, w * h * c * sizeof (real));
    else
        for         (int i = 0; i < h; i++)
            for     (int j = 0; j < w; j++)
                for (int k = 0; k < c; k++)

                    S[(i * w + j) * c + k] = real(p[(i * w + j) * c + k]);
}

// Map and cast the contents of S into the w * h * c floating point raster p.

template <typename real> void Sij<real>::get(float *p)
{
    for         (int i = 0; i < h; i++)
        for     (int j = 0; j < w; j++)
            for (int k = 0; k < c; k++)

                p[(i * w + j) * c + k] = float(S[(i * w + j) * c + k]);
}

//------------------------------------------------------------------------------

#endif

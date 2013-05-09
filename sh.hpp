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
#include <vector>
#include <complex>
#include <algorithm>

//------------------------------------------------------------------------------

// Cm is a complex vector indexed by width 2n. This is used as a working buffer
// for the Fourier transform and inverse.

template <typename real>
    class Cm
    {
    public:
        Cm(int n) : C(2 * n, 0) { }

        typedef std::complex<real> complex;

        complex& operator()(int j)
        {
            return C[j];
        }

        void clear();
        void fft(int);

    private:
        std::vector<complex> C;
    };

// Plm is a real vector indexed by degree l and order m where m is not negative.
// This is the structure of the associated Legendre functions.

template <typename real>
    class Plm
    {
    public:
        Plm(int n) : n(n), P((n + 1) * n / 2) { }

        real& operator()(int l, int m)
        {
            return p(l, m);
        }

        void set(real);

        int n;

    private:
        real a(int)      const;
        real b(int, int) const;
        real c(int, int) const;

        std::vector<real> P;

        real& p(int l, int m)
        {
            return P[(l + 1) * l / 2 + m];
        }
    };

// Flm is a real vector with c channels indexed by degree l and order m.
// This is the structure of a set of spherical harmonic coefficients.

template <typename real>
    class Flm
    {
    public:
        Flm(int n, int c) : n(n), c(c), F(n * n * c, 0) { }

        real& operator()(int l, int m, int k)
        {
            return f(l, m, k);
        }

        void set(const float *, int);
        void get(      float *, int);
        void acc(Flm<real>&);

        void hanning(int);
        void lanczos(int);
        void gauss  (int);
        void diffuse();

        int n;
        int c;

    private:
        std::vector<real> F;

        real& f(int l, int m, int k)
        {
            return F[((l + 1) * l + m) * c + k];
        }
        
        void filter(int, real);
    };

// Sij is a real vector with c channels indexed by width 2n and height 2n.
// This is the structure of an equirectangular spherical image.

template <typename real>
    class Sij
    {
    public:
        Sij(int n, int c) : n(n), c(c), S(4 * n * n * c, 0) { }

        real& operator()(int i, int j, int k)
        {
            return S[(2 * n * i + j) * c + k];
        }

        void set(const float *, int H=0);
        void get(      float *, int H=0);

        int n;
        int c;

    private:
        std::vector<real> S;
    };

// Class sht is the complete spherical harmonic transform implementation.

template <typename real>
    class sht
    {
    public:

        sht(int, int);

        void syn();
        void ana();

        Flm<real> F;
        Sij<real> S;

    private:

        std::vector<int> X;

        void rev(int, int, int, int *);

        real theta (int) const;
        real lambda(int) const;
        real delta (int) const;

        void syni(Plm<real>&,             Cm<real>&, int);
        void anai(Plm<real>&, Flm<real>&, Cm<real>&, int);

        int n;
    };

//------------------------------------------------------------------------------

// Construct a transformer. Allocate storage for ALFs as well as the frequency
// domain and spatial domain representations. Initialize the bit-reversal table.

template <typename real> sht<real>::sht(int n, int c)
    : F(n, c), S(n, c), X(2 * n), n(n)
{
    rev(0, 1, 2 * n, &X.front());
}

template <typename real> void sht<real>::rev(int a, int d, int n, int *x)
{
    if (n == 1)
        *x = a;
    else
    {
        rev(a,     d * 2, n / 2, x        );
        rev(a + d, d * 2, n / 2, x + n / 2);
    }
}

//------------------------------------------------------------------------------

// Compute spherical coordinates for a given image location.

template <typename real> inline real sht<real>::theta(int i) const
{
    return     M_PI * (real(i) + 0.5) / real(2 * n);
}

template <typename real> inline real sht<real>::lambda(int j) const
{
    return 2 * M_PI * (real(j)      ) / real(2 * n);
}

// Calculate the Chebyshev quadrature weight for image row i. Include pi / n
// to account for the post-analysis normalization.

template <typename real> inline real sht<real>::delta(int i) const
{
    real T = theta(i);
    real d = 0.0;

    for (int h = 0; h < n; h++)
    {
        int x = h * 2 + 1;
        d += sin(T * x) / x;
    }

//  return (M_PI / n) * (M_SQRT2 / n) * sin(T) * d;
    return (M_PI / n) * (    2.0 / n) * sin(T) * d;
//  return (M_PI / n) * (    1.0 / n) * sin(T) * d;
}

//------------------------------------------------------------------------------

// Perform the spherical harmonic synthesis of row i.

template <typename real> void sht<real>::syni(Plm<real>& P,
                                              Cm <real>& C, int i)
{
    P.set(cos(theta(i)));

    for (int k = 0; k < S.c; k++)
    {
        C.clear();

        // Acculumate each order in the working buffer, bit-reverse indexed.

        for     (int m = 0; m < n; m++)
            for (int l = m; l < n; l++)
                C(X[m]) += std::complex<real>(F(l,  m, k),
                                              F(l, -m, k)) * P(l, m);

        // Calculate the inverse Fourier transform.

        C.fft(-1);

        // Copy the result to the current row.

        for (int j = 0; j < 2 * n; j++)
            S(i, j, k) = std::real(C(j));
    }
}

// Perform the spherical harmonic synthesis.

template <typename real> void sht<real>::syn()
{
    #pragma omp parallel
    {
        // Each thread has local ALFs and FFT scratch space.

        Plm<real> P(n);
        Cm <real> C(n);

        // Synthesize all rows in parallel.

        int i;

        #pragma omp for
        for (i = 0; i < 2 * n; i++)
            syni(P, C, i);
    }
}

//------------------------------------------------------------------------------

// Perform the spherical harmonic analysis of row i.

template <typename real> void sht<real>::anai(Plm<real>& P,
                                              Flm<real>& F,
                                              Cm <real>& C, int i)
{
    const real d = delta(i);

    P.set(cos(theta(i)));

    for (int k = 0; k < S.c; k++)
    {
        // Copy the current row to the working buffer, bit-reverse indexed.

        for (int j = 0; j < 2 * n; j++)
            C(X[j]) = S(i, j, k);

        // Calculate the forward Fourier transform.

        C.fft(+1);

        // Accumulate the contribution to all amplitudes.

        for     (int m = 0; m < n; m++)
            for (int l = m; l < n; l++)
            {
                F(l,  m, k) += std::real(C(m)) * P(l, m) * d;
                F(l, -m, k) += std::imag(C(m)) * P(l, m) * d;
            }
    }
}

// Perform the spherical harmonic analysis.

template <typename real> void sht<real>::ana()
{
    F.set(0, 0);

    #pragma omp parallel
    {
        // Each thread has local ALFs, FFT scratch space, and results.

        Plm<real> P(n);
        Flm<real> G(n, S.c);
        Cm <real> C(n);

        // Analyze all rows in parallel.

        int i;

        #pragma omp for
        for (i = 0; i < 2 * n; i++)
            anai(P, G, C, i);

        // Accumulate all results.

        #pragma omp critical
        F.acc(G);
    }
}

//------------------------------------------------------------------------------

// Apply the fast Fourier transform or its inverse.

template <typename real> void Cm<real>::fft(int s)
{
    const int n = int(C.size());

    for (int M = 1; M < n; M *= 2)
    {
        const real a = sin(0.5 * s * M_PI / M);
        const real b = sin(      s * M_PI / M);

        complex k = complex(-2.0 * a * a, b);
        complex w = complex( 1.0);

        for (int m = 0; m < M; ++m)
        {
            for (int i = m; i < n; i += M * 2)
            {
                complex t = C[i + M] * w;
                C[i + M]  = C[i    ] - t;
                C[i    ]  = C[i    ] + t;
            }
            w = w + w * k;
        }
    }
}

// Zero the buffer.

template <typename real> void Cm<real>::clear()
{
    std::fill(C.begin(), C.end(), 0);
}

//------------------------------------------------------------------------------
// Compute orthonormalized ALF coefficients

template <typename real> inline real Plm<real>::a(int m) const
{
    return sqrt((2.0 * m + 1.0) /
                (2.0 * m      ));
}

template <typename real> inline real Plm<real>::b(int l, int m) const
{
    return (2.0 * m + 2.0) / sqrt((l + m + 1.0) * (l - m));
}

template <typename real> inline real Plm<real>::c(int l, int m) const
{
    return sqrt(((l + m + 2.0) * (l - m - 1.0)) /
                ((l + m + 1.0) * (l - m      )));
}

template <typename real> void Plm<real>::set(real x)
{
    const real y = sqrt(1.0 - x * x), t = x / y;
    int l;
    int m;

    // Begin with the base case.

    p(0, 0) = 0.282094791773878143474040L;

    // Perform the diagonal recurrence.

    for (m = 1; m < n; m++)
        p(m, m) = y * a(m) * p(m - 1, m - 1);

    // Perform the triangular recurrence.

    for (l = 1; l < n; l++)
    {
        // The first element off the diagonal is simplified by a zero term.

        p(l, l - 1) = t * b(l, l - 1) * p(l, l);

        // Proceed from the diagonal toward column zero.

        for (m = l - 2; m >= 0; m--)
            p(l, m) = t * b(l, m) * p(l, m + 1)
                        - c(l, m) * p(l, m + 2);
    }

    // Pre-multiply the real SH root 2 and Condon-Shortley phase.

    for     (l = 0; l <  n; l++)
        for (m = 1; m <= l; m++)
            p(l, m) *= M_SQRT2 * pow(-1.0, m);
}

//------------------------------------------------------------------------------

// Map and cast the contents of the N * N * c floating point raster p into F.

template <typename real> void Flm<real>::set(const float *p, int N)
{
    std::fill(F.begin(), F.end(), 0);

    if (p)
        for         (int l =  0; l <  N; l++)
            for     (int m = -l; m <= l; m++)
                for (int k =  0; k <  c; k++)
                {
                    int i = m < 0 ? l : l - m;
                    int j = m > 0 ? l : l + m;

                    F[((l + 1) * l + m) * c + k] = real(p[(i * N + j) * c + k]);
                }
}

// Map and cast the contents of F into the N * N * c floating point raster p.

template <typename real> void Flm<real>::get(float *p, int N)
{
    for         (int l =  0; l <  N; l++)
        for     (int m = -l; m <= l; m++)
            for (int k =  0; k <  c; k++)
            {
                int i = m < 0 ? l : l - m;
                int j = m > 0 ? l : l + m;

                p[(i * N + j) * c + k] = float(F[((l + 1) * l + m) * c + k]);
            }
}

// Accumulate the contents of G with F.

template <typename real> void Flm<real>::acc(Flm<real>& G)
{
    for (int i = 0; i < n * n * c; i++)
        F[i] += G.F[i];
}

// Apply a Hanning window of width w.

template <typename real> void Flm<real>::hanning(int w)
{
    for (int l = 0; l < n; l++)
        if (l > w)
            filter(l, 0);
        else
            filter(l, (cos(M_PI * real(l) / real(w)) + 1.0) * 0.5);
}

// Apply a Lanczos window of width w.

template <typename real> void Flm<real>::lanczos(int w)
{
    for (int l = 0; l < n; l++)
        if (l == 0)
            filter(l, 1);
        else
            filter(l, sin(M_PI * real(l) / real(w)) / 
                         (M_PI * real(l) / real(w)));
}

// Apply a Gaussian window of width w.

template <typename real> void Flm<real>::gauss(int w)
{
    for (int l = 0; l < n; l++)
        filter(l, exp(-pow(M_PI * real(l) / real(w), 2.0) / 2.0));
}

// Apply Ramamoorthi & Hanrahan's Lambertian convolution.

template <typename real> void Flm<real>::diffuse()
{
    real a = M_PI;

    if (0 < n) filter(0,     M_PI    );
    if (1 < n) filter(1, 2 * M_PI / 3);

    for (int l = 2; l < n; l += 2)
    {
        a *= 5.0 / (l + 2.0) - 1.0;
        filter(l + 0, a);
        filter(l + 1, 0);
    }
}

// Modulate all coefficients of degree l by scalar a.

template <typename real> void Flm<real>::filter(int l, real a)
{
    // printf("%d %.12lf\n", l, double(a));

    for     (int m = -l; m <= l; m++)
        for (int k =  0; k <  c; k++)
            f(l, m, k) *= a;
}

//------------------------------------------------------------------------------

// Map and cast the contents of the w * H * c floating point raster p into S.

template <typename real> void Sij<real>::set(const float *p, int H)
{
    const int w = 2 * n;
    const int h = 2 * n;

    std::fill(S.begin(), S.end(), 0);

    // Copy 2n rows into 2n rows.

    if (p && (H == h || H == 0))
        for         (int i = 0; i < h; i++)
            for     (int j = 0; j < w; j++)
                for (int k = 0; k < c; k++)

                    S[(i * w + j) * c + k] = real(p[(i * w + j) * c + k]);

    // Copy n rows into 2n rows, duplicating each row and scaling by half.

    if (p && (H == n))
        for         (int i = 0; i < n; i++)
            for     (int j = 0; j < w; j++)
                for (int k = 0; k < c; k++)
                {
                    float f = p[(i * w + j) * c + k];

                    S[((2 * i + 0) * w + j) * c + k] = real(f) / 2;
                    S[((2 * i + 1) * w + j) * c + k] = real(f) / 2;
                }
}

// Map and cast the contents of S into the w * h * c floating point raster p.

template <typename real> void Sij<real>::get(float *p, int H)
{
    const int w = 2 * n;
    const int h = 2 * n;

    // Copy 2n rows into 2n rows.

    if (H == h || H == 0)
        for         (int i = 0; i < h; i++)
            for     (int j = 0; j < w; j++)
                for (int k = 0; k < c; k++)

                    p[(i * w + j) * c + k] = float(S[(i * w + j) * c + k]);

    // Copy 2n rows int n rows, summing each pair of rows.

    if (H == n)
        for         (int i = 0; i < n; i++)
            for     (int j = 0; j < w; j++)
                for (int k = 0; k < c; k++)
                {
                    real r = S[((2 * i + 0) * w + j) * c + k];
                    real s = S[((2 * i + 1) * w + j) * c + k];

                    p[(i * w + j) * c + k] = float(r + s);
                }
}

//------------------------------------------------------------------------------

#endif

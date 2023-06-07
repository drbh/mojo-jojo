#include <iostream>
#include <vector>
#include <complex>
#include <cmath>

// function to reverse bits of a number
unsigned int bit_reverse(unsigned int x, int log2n)
{
    int n = 0;
    int mask = 0x1;
    for (int i = 0; i < log2n; i++)
    {
        n <<= 1;
        n |= (x & 1);
        x >>= 1;
    }
    return n;
}

void fft(std::vector<std::complex<float>> &a)
{
    // Bit-reversal permutation
    unsigned int n = a.size();
    unsigned int log2n = std::log2(n);
    // 👈
    std::vector<std::complex<float>> a2(n);
    for (unsigned int i = 0; i < n; ++i)
    {
        unsigned int rev = bit_reverse(i, log2n);
        a2[rev] = a[i];
    }

    // In-place Cooley-Tukey FFT
    std::complex<float> im(0, 1);
    for (unsigned int s = 1; s <= log2n; ++s)
    {
        unsigned int m = 1 << s;      // 2 power s
        unsigned int half_m = m >> 1; // m/2
        // Pre-compute twiddle factors
        std::vector<std::complex<float>> w(half_m);
        for (unsigned int j = 0; j < half_m; ++j)
        {
            float theta = -2.0 * M_PI * j / m;
            w[j] = std::exp(im * theta);
        }
        for (unsigned int k = 0; k < n; k += m)
        {
            for (unsigned int j = 0; j < half_m; ++j)
            {
                std::complex<float> t = w[j] * a2[k + j + half_m];
                std::complex<float> u = a2[k + j];
                a2[k + j] = u + t;
                a2[k + j + half_m] = u - t;
            }
        }
    }
    a = a2;
}

int main()
{
    std::vector<std::complex<float>> a = {1, 2, 3, 4, 5, 6, 7, 8};
    fft(a);
    for (auto i : a)
    {
        std::cout << i << std::endl;
    }
    return 0;
}

/// Comparable to numpy.fft.fft
// import numpy as np
// a = np.array([1,2,3,4,5,6,7,8])
// np.fft.fft(a)
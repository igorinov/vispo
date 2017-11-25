# pvft
Portable Vectorized Fourier Transform

This project is a cross-platform library for Fast Fourier Transform.

FFT algorithm used in this library allows full utilization of SIMD processor
extensions by optimized assembly code or automatic vectorization by compiler.

If the signal size is even but not power of 2, FFT synthesis is used to combine
multiple DFT results.

DFT functions allow arbitrary data size, using compensated summation (Kahan's
algorithm) to minimize error accumulation.
They can also be vectorized to improve compensated summation performance.


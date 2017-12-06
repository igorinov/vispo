/**
 *  Copyright (C) 2016 Ivan Gorinov
 *  License: BSD
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "pvft.h"

#if (__STDC_VERSION__ >= 199901L)
#define USE_FMA
#endif

#ifdef ASSEMBLY_FFT
int fft_combine_s(const complex_s *s, complex_s *data, int m, int n);
int fft_combine_d(const complex_d *s, complex_d *data, int m, int n);
#endif

/**
 *  dft_setup() - prepare complex sinusoid table to use with dft_... functions
 *  for fft_... functions, use fft_setup()
 *  @s - buffer for sinusoid table, @n elements of complex type
 *  @n - size of the transform
 *  @inverse - direction: 0 - forward, 1 - inverse
 */

int dft_setup_d(complex_d *s, int n, int inverse)
{
	double da;
	double t;
	int sign = inverse ? 1 : -1;
	int c, q, i;

	if (n == 0)
		return 0;

	da = 2 * M_PI / n;

	s[0].re = 1;
	s[0].im = 0;

	if ((n & 1) == 0) {
		c = n / 2;
		if ((c & 1) == 0) {
			q = c / 2;
			for (i = 1; i < q; i += 1) {
				t = i * da;
				s[i].re = cos(t);
				s[q - i].im = s[i].re * sign;
			}
			s[q].re = 0;
			s[q].im = sign;
			for (i = 1; i < q; i += 1) {
				s[q + i].re = - s[q - i].re;
				s[q + i].im = s[q - i].im;
			}
		} else {
			for (i = 1; i < c; i += 1) {
				t = i * da;
				s[i].re = cos(t);
				s[i].im = sin(t) * sign;
			}
		}
		for (i = 0; i < c; i += 1) {
			s[c + i].re = - s[i].re;
			s[c + i].im = - s[i].im;
		}
	} else {
		for (i = 1; i < n; i += 1) {
			t = i * da;
			s[i].re = cos(t);
			s[i].im = sin(t) * sign;
		}
	}

	return n;
}

int dft_setup_s(complex_s *s, int n, int inverse)
{
	double da;
	double t;
	int sign = inverse ? 1 : -1;
	int c, q, i;

	if (n == 0)
		return 0;

	da = 2 * M_PI / n;

	s[0].re = 1;
	s[0].im = 0;

	if ((n & 1) == 0) {
		c = n / 2;
		if ((c & 1) == 0) {
			q = c / 2;
			for (i = 1; i < q; i += 1) {
				t = i * da;
				s[i].re = cos(t);
				s[q - i].im = s[i].re * sign;
			}
			s[q].re = 0;
			s[q].im = sign;
			for (i = 1; i < q; i += 1) {
				s[q + i].re = - s[q - i].re;
				s[q + i].im = s[q - i].im;
			}
		} else {
			for (i = 1; i < c; i += 1) {
				t = i * da;
				s[i].re = cos(t);
				s[i].im = sin(t) * sign;
			}
		}
		for (i = 0; i < c; i += 1) {
			s[c + i].re = - s[i].re;
			s[c + i].im = - s[i].im;
		}
	} else {
		for (i = 1; i < n; i += 1) {
			t = i * da;
			s[i].re = cos(t);
			s[i].im = sin(t) * sign;
		}
	}

	return n;
}

/**
 *  fft_setup() - prepare complex sinusoid tables for fft_... functions
 *  for dft... functions, table (cs + n) can be used
 *  @cs - buffer for sinusoid table, @n * 2 elements of complex type
 *  @n - size of the transform
 *  @inverse - direction: 0 - forward, 1 - inverse
 */

int fft_setup_d(complex_d *cs, int n, int inverse)
{
	int i, j;

	if (n == 0)
		return 0;

	dft_setup_d(cs + n, n, inverse);

	while (!(n & 1)) {
		n >>= 1;
		j = n * 2;
		for (i = 0; i < n; i += 1) {
			cs[n + i] = cs[j];
			j += 2;
		}
	}

	for (i = 0; i < n; i += 1) {
		cs[i].re = 0;
		cs[i].im = 0;
	}

	return n;
}

int fft_setup_s(complex_s *cs, int n, int inverse)
{
	int i, j;

	if (n == 0)
		return 0;

	dft_setup_s(cs + n, n, inverse);

	while (!(n & 1)) {
		n >>= 1;
		j = n * 2;
		for (i = 0; i < n; i += 1) {
			cs[n + i] = cs[j];
			j += 2;
		}
	}

	for (i = 0; i < n; i += 1) {
		cs[i].re = 0;
		cs[i].im = 0;
	}

	return n;
}

/**
 *  dft_complex() - simple DFT, any size
 *  @s - sinusoid table prepared by dft_setup()
 *  @data - input data
 *  @out - output buffer
 *  @n - number of points in @data and @out
 */

int dft_complex_s(const complex_s *cs, complex_s *output, const complex_s *input, int n)
{
	return dft_complex_step_s(cs, output, input, n, 1);
}

int dft_complex_d(const complex_d *cs, complex_d *output, const complex_d *input, int n)
{
	return dft_complex_step_d(cs, output, input, n, 1);
}

#ifndef ASSEMBLY_DFT

/**
 *  dft_complex_step() - compute DFT of complex signal
 *  using compensated summation (Kahan's algorithm) to
 *  reduce numerical error
 *  @s - complex sinusoid table prepared by dft_setup()
 *  @input - input data
 *  @output - output buffer
 *  @n - number of rows
 *  @in_step - input data interleaving step
 */

int dft_complex_step_d(const complex_d *s, complex_d *output, const complex_d *input, int n, int in_step)
{
	const complex_d *px, *ps;
	double acc[4];
	double inc[4];
	double sum[4];
	complex_d x;
	int i, j, k, si;

	/* DC component */

	px = input;
	acc[0] = Re(*px);
	acc[1] = Im(*px);
	inc[0] = 0;
	inc[1] = 0;
	for (i = 1; i < n; i += 1) {
		px += in_step;
		inc[0] += Re(*px);
		inc[1] += Im(*px);
		for (j = 0; j < 2; j += 1) {
			sum[j] = acc[j] + inc[j];
			inc[j] -= (sum[j] - acc[j]);
			acc[j] = sum[j];
		}
	}
	output[0].re = acc[0];
	output[0].im = acc[1];

	for (k = 1; k < n; k += 1) {
		px = input;
		acc[0] = Re(*px);
		acc[1] = Im(*px);
		acc[2] = 0;
		acc[3] = 0;

		for (j = 0; j < 4; j += 1)
			inc[j] = 0;

        	si = k;
		for (i = 1; i < n; i += 1) {
			px += in_step;
			ps = s + si;
			x = *px;

#ifdef USE_FMA
			inc[0] = fma(Re(*ps), Re(x), inc[0]);
			inc[1] = fma(Re(*ps), Im(x), inc[1]);
			inc[2] = fma(Im(*ps), Re(x), inc[2]);
			inc[3] = fma(Im(*ps), Im(x), inc[3]);
#else
			inc[0] += x * ps->x;
			inc[1] += y * ps->x;
			inc[2] += x * ps->y;
			inc[3] += y * ps->y;
#endif
			/* this should be vectorized by SIMD-enabled compiler */

			for (j = 0; j < 4; j += 1)
				sum[j] = acc[j] + inc[j];

			for (j = 0; j < 4; j += 1)
				inc[j] -= (sum[j] - acc[j]);

			for (j = 0; j < 4; j += 1)
				acc[j] = sum[j];

			/*
			 * si = (i * k) % n
			 */
			si += k;
			if (si >= n)
				si -= n;
		}
		output[k].re = (acc[0] - acc[3]) + (inc[0] - inc[3]);
		output[k].im = (acc[1] + acc[2]) + (inc[1] + inc[2]);
	}

	return k;
}

int dft_complex_step_s(const complex_s *s, complex_s *output, const complex_s *input, int n, int in_step)
{
	const complex_s *px, *ps;
	float acc[4];
	float inc[4];
	float sum[4];
	complex_s x;
	int i, j, k, si;

	/* DC component */

	px = input;
	acc[0] = Re(*px);
	acc[1] = Im(*px);
	inc[0] = 0;
	inc[1] = 0;
	for (i = 1; i < n; i += 1) {
		px += in_step;
		inc[0] += Re(*px);
		inc[1] += Im(*px);
		for (j = 0; j < 2; j += 1) {
			sum[j] = acc[j] + inc[j];
			inc[j] -= (sum[j] - acc[j]);
			acc[j] = sum[j];
		}
	}
	output[0].re = acc[0];
	output[0].im = acc[1];

	for (k = 1; k < n; k += 1) {
		px = input;
		acc[0] = Re(*px);
		acc[1] = Im(*px);
		acc[2] = 0;
		acc[3] = 0;

		for (j = 0; j < 4; j += 1)
			inc[j] = 0;
        	si = k;
		for (i = 1; i < n; i += 1) {
			px += in_step;
			ps = s + si;
			x = *px;

#ifdef USE_FMA
			inc[0] = fma(Re(*ps), Re(x), inc[0]);
			inc[1] = fma(Re(*ps), Im(x), inc[1]);
			inc[2] = fma(Im(*ps), Re(x), inc[2]);
			inc[3] = fma(Im(*ps), Im(x), inc[3]);
#else
			inc[0] += Re(a) * Re(*ps);
			inc[1] += Im(a) * Re(*ps);
			inc[2] += Re(a) * Im(*ps);
			inc[3] += Im(a) * Im(*ps);
#endif
			/* this should be vectorized by SIMD-enabled compiler */

			for (j = 0; j < 4; j += 1)
				sum[j] = acc[j] + inc[j];

			for (j = 0; j < 4; j += 1)
				inc[j] -= (sum[j] - acc[j]);

			for (j = 0; j < 4; j += 1)
				acc[j] = sum[j];

			/*
			 * si = (i * k) % n
			 */
			si += k;
			if (si >= n)
				si -= n;
		}
		output[k].re = (acc[0] - acc[3]) + (inc[0] - inc[3]);
		output[k].im = (acc[1] + acc[2]) + (inc[1] + inc[2]);
	}

	return k;
}

#endif  /* ASSEMBLY_DFT */

/**
 *  dft_real_step() - compute DFT of real signal
 *  using Kahan's algorithm to reduce numerical error
 *  @s - complex sinusoid table prepared by dft_setup()
 *  @input - input data
 *  @output - output buffer
 *  @n - number of rows
 *  @in_step - interleaving step
 */

int dft_real_step_d(const complex_d *s, complex_d *output, const double *input, int n, int in_step)
{
	const complex_d *ps;
	const double *pt;
	double acc[2];
	double inc[2];
	double sum[2];
	double x, y;
	int i, j, k, si;

	/* DC component */

	pt = input;
	acc[0] = *pt;
	inc[0] = 0;
	for (i = 1; i < n; i += 1) {
		pt += in_step;
		inc[0] += *pt;
		sum[0] = acc[0] + inc[0];
		inc[0] -= (sum[0] - acc[0]);
		acc[0] = sum[0];
	}
	output[0].re = acc[0];
	output[0].im = 0;

	for (k = 1; k < n; k += 1) {
		pt = input;
		acc[0] = *pt;
		acc[1] = 0;
		inc[0] = 0;
		inc[1] = 0;

		for (j = 0; j < 2; j += 1)
			inc[j] = 0;

        	si = k;
		for (i = 1; i < n; i += 1) {
			pt += in_step;
			ps = s + si;
			x = *pt;

#ifdef USE_FMA
			inc[0] = fma(Re(*ps), x, inc[0]);
			inc[1] = fma(Im(*ps), x, inc[1]);
#else
			inc[0] += x * Re(*ps);
			inc[1] += x * Im(*ps);
#endif
			/* this should be vectorized by SIMD-enabled compiler */

			for (j = 0; j < 2; j += 1)
				sum[j] = acc[j] + inc[j];

			for (j = 0; j < 2; j += 1)
				inc[j] -= (sum[j] - acc[j]);

			for (j = 0; j < 2; j += 1)
				acc[j] = sum[j];

			/*
			 * si = (i * k) % n
			 */
			si += k;
			if (si >= n)
				si -= n;
		}
		output[k].re = acc[0];
		output[k].im = acc[1];
	}

	return k;
}

int dft_real_step_s(const complex_s *s, complex_s *output, const float *input, int n, int in_step)
{
	const complex_s *ps;
	const float *pt;
	float acc[2];
	float inc[2];
	float sum[2];
	float x, y;
	int i, j, k, si;

	/* DC component */

	pt = input;
	acc[0] = *pt;
	inc[0] = 0;
	for (i = 1; i < n; i += 1) {
		pt += in_step;
		inc[0] += *pt;
		sum[0] = acc[0] + inc[0];
		inc[0] -= (sum[0] - acc[0]);
		acc[0] = sum[0];
	}
	output[0].re = acc[0];
	output[0].im = 0;

	for (k = 1; k < n; k += 1) {
		pt = input;
		acc[0] = *pt;
		acc[1] = 0;
		inc[0] = 0;
		inc[1] = 0;

		for (j = 0; j < 2; j += 1)
			inc[j] = 0;

        	si = k;
		for (i = 1; i < n; i += 1) {
			pt += in_step;
			ps = s + si;
			x = *pt;

#ifdef USE_FMA
			inc[0] = fma(ps->re, x, inc[0]);
			inc[1] = fma(ps->im, x, inc[1]);
#else
			inc[0] += x * ps->re;
			inc[1] += x * ps->im;
#endif
			/* this should be vectorized by SIMD-enabled compiler */

			for (j = 0; j < 2; j += 1)
				sum[j] = acc[j] + inc[j];

			for (j = 0; j < 2; j += 1)
				inc[j] -= (sum[j] - acc[j]);

			for (j = 0; j < 2; j += 1)
				acc[j] = sum[j];

			/*
			 * si = (i * k) % n
			 */
			si += k;
			if (si >= n)
				si -= n;
		}
		output[k].re = acc[0];
		output[k].im = acc[1];
	}

	return k;
}

#ifndef ASSEMBLY_FFT

/**
 *  fft_combine() - combine M spectrum pairs, N points each
 *  into M spectra, N * 2 points each
 *
 *  @pt_w - complex sinusoid table for this stage
 *  @l - signal length, complex samples
 *  @m - number of signals
 */

int fft_combine_d(const complex_d *pt_w, complex_d *data, int l, int m)
{
	complex_d w, a, b, t;
	int i, j;

	for (i = 0; i < m; i += 1) {
		for (j = 0; j < l; j += 1) {
			w = pt_w[j];
			a = data[j];
			b = data[l + j];
			t.re = b.re * w.re - b.im * w.im;
			t.im = b.im * w.re + b.re * w.im;
			b.re = a.re - t.re;
			b.im = a.im - t.im;
			a.re += t.re;
			a.im += t.im;
			data[j] = a;
			data[l + j] = b;
		}
		data += l * 2;
	}
}

int fft_combine_s(const complex_s *pt_w, complex_s *data, int l, int m)
{
	complex_s w, a, b, t;
	int i, j;

	for (i = 0; i < m; i += 1) {
		for (j = 0; j < l; j += 1) {
			w = pt_w[j];
			a = data[j];
			b = data[l + j];
			t.re = b.re * w.re - b.im * w.im;
			t.im = b.im * w.re + b.re * w.im;
			b.re = a.re - t.re;
			b.im = a.im - t.im;
			a.re += t.re;
			a.im += t.im;
			data[j] = a;
			data[l + j] = b;
		}
		data += l * 2;
	}
}

#endif  /* ASSEMBLY_FFT */

/**
 *  fft_combine_interleaved() - one stage of interleaved transform
 *  combine M spectrum pairs, N points each
 *  into M spectra, N * 2 points each (interleaved data)
 *
 *  @pt_w - complex sinusoid table for this stage
 *  @l - signal length, complex samples
 *  @m - number of signals in one channel
 *  @n - number of interleaved channels
 */

int fft_combine_interleaved_d(const complex_d *pt_w, complex_d *data, int l, int m, int d)
{
	complex_d w, a, b, t;
	int ld = l * d;
	int i, j, k;

	for (i = 0; i < m; i += 1) {
		for (j = 0; j < l; j += 1) {
			w = pt_w[j];
			for (k = 0; k < d; k += 1) {
				a = data[k];
				b = data[ld + k];
				t.re = b.re * w.re - b.im * w.im;
				t.im = b.im * w.re + b.re * w.im;
				b.re = a.re - t.re;
				b.im = a.im - t.im;
				a.re += t.re;
				a.im += t.im;
				data[k] = a;
				data[ld + k] = b;
			}
			data += k;
		}
		data += ld;
	}
}

int fft_combine_interleaved_s(const complex_s *pt_w, complex_s *data, int l, int m, int d)
{
	complex_s w, a, b, t;
	int ld = l * d;
	int i, j, k;

	for (i = 0; i < m; i += 1) {
		for (j = 0; j < l; j += 1) {
			w = pt_w[j];
			for (k = 0; k < d; k += 1) {
				a = data[k];
				b = data[ld + k];
				t.re = b.re * w.re - b.im * w.im;
				t.im = b.im * w.re + b.re * w.im;
				b.re = a.re - t.re;
				b.im = a.im - t.im;
				a.re += t.re;
				a.im += t.im;
				data[k] = a;
				data[ld + k] = b;
			}
			data += k;
		}
		data += ld;
	}
}

/**
 *  fft_combine_loop() - combine M spectrum pairs, L points each
 *  into one spectrum, L * M * 2 points
 *
 *  @cs - complex sinusoid tables, prepared by fft_setup()
 */

int fft_combine_loop_s(const complex_s *ww, complex_s *data, int l, int m)
{
	while (!(m & 1)) {
		m >>= 1;
		fft_combine_s(ww + l * 2, data, l, m);
		l <<= 1;
	}

	return m;
}

int fft_combine_loop_d(const complex_d *ww, complex_d *data, int l, int m)
{
	while (!(m & 1)) {
		m >>= 1;
		fft_combine_d(ww + l * 2, data, l, m);
		l <<= 1;
	}

	return m;
}

/**
 *  fft_combine_loop_step() - combine M spectrum pairs, L points each
 *  into one spectrum, L * M * 2 points (D interleaved channels)
 *
 *  @cs - complex sinusoid tables, prepared by fft_setup()
 */

int fft_combine_loop_step_s(const complex_s *ww, complex_s *data, int l, int m, int d)
{
	while (!(m & 1)) {
		m >>= 1;
		fft_combine_interleaved_s(ww + l * 2, data, l, m, d);
		l <<= 1;
	}

	return m;
}

int fft_combine_loop_step_d(const complex_d *ww, complex_d *data, int l, int m, int d)
{
	while (!(m & 1)) {
		m >>= 1;
		fft_combine_interleaved_d(ww + l * 2, data, l, m, d);
		l <<= 1;
	}

	return m;
}

#ifndef BITREV

int bit_reverse(int i, int bits)
{
	unsigned int x = i;

	x = ((x & 0xffff0000) >> 16) | ((x & 0x0000ffff) << 16);
	x = ((x & 0xff00ff00) >>  8) | ((x & 0x00ff00ff) <<  8);
	x = ((x & 0xf0f0f0f0) >>  4) | ((x & 0x0f0f0f0f) <<  4);
	x = ((x & 0xcccccccc) >>  2) | ((x & 0x33333333) <<  2);
	x = ((x & 0xaaaaaaaa) >>  1) | ((x & 0x55555555) <<  1);

	return x >> (32 - bits);
}

#endif

/**
 *  fft_complex() - compute FFT from complex data
 */

int fft_complex_d(const complex_d *ww, complex_d *out, const complex_d *in, int size)
{
	int i, j, k;
	int bits = 0;
	int m = 1;
	int l = size;
	complex_d *dft_ptr;

	/* length = l * m,  m = 2 ^ bits */

	while (l > 1) {
		if (l & 1)
			break;
		l >>= 1;
		m <<= 1;
		bits += 1;
	}

	for (i = 0; i < m; i += 1) {
		j = bit_reverse(i, bits);

		if (l > 1) {
			dft_complex_step_d(ww + l, out + j * l, in + i, l, m);
		} else {
			out[j] = in[i];
		}
	}

	return fft_combine_loop_d(ww, out, l, m);
}

int fft_complex_s(const complex_s *ww, complex_s *out, const complex_s *in, int size)
{
	int i, j, k;
	int bits = 0;
	int m = 1;
	int l = size;
	complex_s *dft_ptr;

	/* length = l * m,  m = 2 ^ bits */

	while (l > 1) {
		if (l & 1)
			break;
		l >>= 1;
		m <<= 1;
		bits += 1;
	}

	for (i = 0; i < m; i += 1) {
		j = bit_reverse(i, bits);

		if (l > 1) {
			dft_complex_step_s(ww + l, out + j * l, in + i, l, m);
		} else {
			out[j] = in[i];
		}
	}

	return fft_combine_loop_s(ww, out, l, m);
}

/**
 *  fft_complex_step() - compute FFT from complex data, with interleaving
 */

int fft_complex_step_d(const complex_d *ww, complex_d *out, const complex_d *in, int size, int step)
{
	int i, j, k;
	int bits = 0;
	int m = 1;
	int l = size;
	complex_d *dft_ptr;

	/* length = m * n,  m = 2 ^ bits */

	while (l > 1) {
		if (l & 1)
			break;
		l >>= 1;
		m <<= 1;
		bits += 1;
	}

	for (i = 0; i < m; i += 1) {
		j = bit_reverse(i, bits) * step;

		if (l > 1) {
			dft_complex_step_d(ww + l, out + j * l, in + i * step, l, m * step);
		} else {
			out[j] = in[i * step];
		}
	}

	return fft_combine_loop_step_d(ww, out, l, m, step);
}

int fft_complex_step_s(const complex_s *ww, complex_s *out, const complex_s *in, int size, int step)
{
	int i, j, k;
	int bits = 0;
	int m = 1;
	int l = size;
	complex_s *dft_ptr;

	/* length = l * m,  m = 2 ^ bits */

	while (l > 1) {
		if (l & 1)
			break;
		l >>= 1;
		m <<= 1;
		bits += 1;
	}

	for (i = 0; i < m; i += 1) {
		j = bit_reverse(i, bits) * step;

		if (l > 1) {
			dft_complex_step_s(ww + l, out + j * l, in + i * step, l, m * step);
		} else {
			out[j] = in[i * step];
		}
	}

	return fft_combine_loop_step_s(ww, out, l, m, step);
}

/**
 *  fft_real_odd() - compute FFT from real data, odd size
 */

int fft_real_odd_d(const complex_d *cs, complex_d *out, const double *in, int size)
{
	int i, j, k;
	int bits = 0;
	int m = 1;
	int n = size;

	/* length = m * n,  m = 2 ^ bits */

	while (n > 1) {
		if (n & 1)
			break;
		n >>= 1;
		m <<= 1;
		bits += 1;
	}

	for (i = 0; i < m; i += 1) {
		j = bit_reverse(i, bits);

		if (n > 1) {
			dft_real_step_d(cs + n, out + j * n, in + i, n, m);
		} else {
			out[j].re = in[i];
			out[j].im = 0;
		}
	}

	return fft_combine_loop_d(cs, out, n, m);
}


int fft_real_odd_s(const complex_s *cs, complex_s *out, const float *in, int size)
{
	int i, j, k;
	int bits = 0;
	int m = 1;
	int n = size;

	/* length = m * n,  m = 2 ^ bits */

	while (n > 1) {
		if (n & 1)
			break;
		n >>= 1;
		m <<= 1;
		bits += 1;
	}

	for (i = 0; i < m; i += 1) {
		j = bit_reverse(i, bits);

		if (n > 1) {
			dft_real_step_s(cs + n, out + j * n, in + i, n, m);
		} else {
			out[j].re = in[i];
			out[j].im = 0;
		}
	}

	if (m == 1)
		return 1;

	return fft_combine_loop_s(cs, out, n, m);
}

/**
 *  fft_real() - compute FFT from real data
 *  @sig - real input signal
 *  @n - FFT size
 */

int fft_real_d(const complex_d *cs, complex_d *out, const double *input, int n)
{
	complex_d a, b;
	int c, quarter;
	int i, j;

	if (n & 3)
		return fft_real_odd_d(cs, out, input, n);

	c = n / 2;
	quarter = n / 4;

	/*
	 *  casting the pointer turns the real signal
	 *  into complex signal with size N/2
	 *  sig'[i].re = sig[i * 2]
	 *  sig'[i].im = sig[i * 2 + 1]
	 */

	fft_complex_d(cs, out, (const complex_d *) input, c);

	/* odd-even decomposition */

	i = 0;

	out[c + i].re = out[i].im;
	out[c + i].im = 0;
	out[i].im = 0;

	for (i = 1; i < quarter; i += 1) {
		j = c - i;
		a.re = out[i].re / 2;
		a.im = out[i].im / 2;
		b.re = out[j].re / 2;
		b.im = out[j].im / 2;
		out[c + i].re = b.im + a.im;
		out[c + i].im = b.re - a.re;
		out[c + j].re = a.im + b.im;
		out[c + j].im = a.re - b.re;
		out[j].re = b.re + a.re;
		out[j].im = b.im - a.im;
		out[i].re = a.re + b.re;
		out[i].im = a.im - b.im;
	}

	out[c + i].re = out[i].im;
	out[c + i].im = 0;
	out[i].im = 0;

	return fft_combine_loop_d(cs, out, c, 2);
}

int fft_real_s(const complex_s *cs, complex_s *out, const float *input, int n)
{
	complex_s a, b;
	int c, quarter;
	int i, j;

	if (n & 3)
		return fft_real_odd_s(cs, out, input, n);

	c = n / 2;
	quarter = n / 4;

	/*
	 *  casting the pointer turns the real signal
	 *  into complex signal with size N/2
	 *  sig'[i].re = sig[i * 2]
	 *  sig'[i].im = sig[i * 2 + 1]
	 */

	fft_complex_s(cs, out, (const complex_s *) input, c);

	/* odd-even decomposition */

	i = 0;

	out[c + i].re = out[i].im;
	out[c + i].im = 0;
	out[i].im = 0;

	for (i = 1; i < quarter; i += 1) {
		j = c - i;
		a.re = out[i].re / 2;
		a.im = out[i].im / 2;
		b.re = out[j].re / 2;
		b.im = out[j].im / 2;
		out[c + i].re = b.im + a.im;
		out[c + i].im = b.re - a.re;
		out[c + j].re = a.im + b.im;
		out[c + j].im = a.re - b.re;
		out[j].re = b.re + a.re;
		out[j].im = b.im - a.im;
		out[i].re = a.re + b.re;
		out[i].im = a.im - b.im;
	}

	out[c + i].re = out[i].im;
	out[c + i].im = 0;
	out[i].im = 0;

	return fft_combine_loop_s(cs, out, c, 2);
}


/**
 *  Copyright (C) 2016 Ivan Gorinov
 *
 *  SPDX-License-Identifier: BSD-2-Clause
 *
 *  This default implementation is built if $(ARCH)/combine.S is not present
 */

#include "pvft.h"

/**
 *  fft_combine() - combine M spectrum pairs, L points each
 *  into M spectra, L * 2 points each
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


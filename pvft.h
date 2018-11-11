/**
 *  Copyright (C) 2016 Ivan Gorinov
 *  License: BSD
 */

#ifndef __pvft_h
#define __pvft_h

#define Re(a) (a).re
#define Im(a) (a).im

typedef struct {
    double re;
    double im;
} complex_d;

typedef struct {
    float re;
    float im;
} complex_s;


#ifdef __cplusplus
extern "C" {
#endif

void pvft_init(void);

int fft_setup_s(complex_s *ww, int n, int inverse);
int fft_setup_d(complex_d *ww, int n, int inverse);

int fft_complex_s(const complex_s *ww, complex_s *out, const complex_s *in, int n);
int fft_complex_d(const complex_d *ww, complex_d *out, const complex_d *in, int n);

int fft_real_s(const complex_s *ww, complex_s *out, const float *in, int n);
int fft_real_d(const complex_d *ww, complex_d *out, const double *in, int n);

int dft_complex_step_s(const complex_s *w, complex_s *output, const complex_s *input, int n, int in_step);
int dft_complex_step_d(const complex_d *w, complex_d *output, const complex_d *input, int n, int in_step);

int dft_real_step_s(const complex_s *w, complex_s *output, const float *input, int n, int in_step, int f_step);
int dft_real_step_d(const complex_d *w, complex_d *output, const double *input, int n, int in_step, int f_step);

int dft_complex_s(const complex_s *w, complex_s *output, const complex_s *input, int n);
int dft_complex_d(const complex_d *w, complex_d *output, const complex_d *input, int n);

int bit_reverse(int i, int bits);
int fft_combine_s(const complex_s *s, complex_s *data, int m, int n);
int fft_combine_d(const complex_d *s, complex_d *data, int m, int n);

#ifdef __cplusplus
}
#endif

#endif  /* __pvft_h */


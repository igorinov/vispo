/**
 *  Copyright (C) Ivan Gorinov, 2016
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <fenv.h>
#include <time.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/mman.h>

#include "vispo.h"

typedef float single;

#define DFT_MAX 0x01000000

#define FFT_MAX_SINGLE 65536

int digits = 0;

void *alloc_large(size_t size)
{
    void *p;

    p = mmap(NULL, size,
        PROT_READ | PROT_WRITE,
        MAP_ANONYMOUS | MAP_PRIVATE,
        -1, 0);

    if (p == MAP_FAILED) {
        printf("Failed to allocate %lu bytes\n", (unsigned long) size);
        exit(1);
    }

    return p;
}

complex_d *cs_forward_d;
complex_d *cs_inverse_d;
complex_s *cs_forward_s;
complex_s *cs_inverse_s;

double *real_data_d;
single *real_data_s;
complex_d *data_d;
complex_s *data_s;
complex_d *data2_d;
complex_s *data2_s;
complex_d *spec1_d;
complex_s *spec1_s;
complex_d *spec2_d;
complex_s *spec2_s;

/**
 *  signal_power() - compute signal power
 *  Using Kahan's algorithm to reduce numerical error
 *  @signal - input complex signal
 */

double signal_power_d(complex_d *signal, int n)
{
    double inc[2], acc[2], sum[2];
    double x, y;
    int i;

    x = signal[0].re;
    y = signal[0].im;
    acc[0] = x * x;
    acc[1] = y * y;
    inc[0] = 0;
    inc[1] = 0;

    for (i = 1; i < n; i += 1) {
        x = signal[i].re;
        y = signal[i].im;
        inc[0] += x * x;
        inc[1] += y * y;
        sum[0] = acc[0] + inc[0];
        sum[1] = acc[1] + inc[1];
        inc[0] -= (sum[0] - acc[0]);
        inc[1] -= (sum[1] - acc[1]);
        acc[0] = sum[0];
        acc[1] = sum[1];
    }

    return (acc[0] + acc[1]) + (inc[0] + inc[1]);
}

double signal_power_s(complex_s *signal, int n)
{
    double inc[2], acc[2], sum[2];
    double x, y;
    int i;

    x = signal[0].re;
    y = signal[0].im;
    acc[0] = x * x;
    acc[1] = y * y;
    inc[0] = 0;
    inc[1] = 0;

    for (i = 1; i < n; i += 1) {
        x = signal[i].re;
        y = signal[i].im;
        inc[0] += x * x;
        inc[1] += y * y;
        sum[0] = acc[0] + inc[0];
        sum[0] = acc[1] + inc[1];
        inc[0] -= (sum[0] - acc[0]);
        inc[1] -= (sum[1] - acc[1]);
        acc[0] = sum[0];
        acc[1] = sum[1];
    }

    return (acc[0] + acc[1]) + (inc[0] + inc[1]);
}

int compare_arrays_d(complex_d *a, complex_d *b, int n, double eps)
{
    double dx, dy;
    double eps_sq = eps * eps;
    double err_sq;
    int i;

    for (i = 0; i < n; i += 1) {
        dx = a[i].re - b[i].re;
        dy = a[i].im - b[i].im;
        err_sq = dx * dx + dy * dy;
        if (err_sq > eps_sq)
            return i;
    }

    return n;
}

int compare_arrays_s(complex_s *a, complex_s *b, int n, float eps)
{
    float dx, dy;
    float eps_sq = eps * eps;
    float err_sq;
    int i;

    for (i = 0; i < n; i += 1) {
        dx = a[i].re - b[i].re;
        dy = a[i].im - b[i].im;
        err_sq = dx * dx + dy * dy;
        if (err_sq > eps_sq)
            return i;
    }

    return n;
}

double time_diff(struct timespec *t0, struct timespec *t1)
{
    return (t1->tv_sec - t0->tv_sec) + (t1->tv_nsec - t0->tv_nsec) * 1e-9;
}

int real_test_d(short *filedata, int m, int n, int show)
{
    struct timespec t0, t1;
    double r_m = 1.0 / m;
    double r_n = 1.0 / n;
    double eps_d = log(n) / (log(2) * 65536);
    double dt_complex_d = 0;
    double dt_real_d = 0;
    int pass = 1;
    double x, y;
    int off, i, j, k;
    long l;
    struct vispo_fft_d forward, inverse;
    int err;

    err = vispo_fft_alloc_d(&forward, n);
    if (err) {
        printf("Error allocating %d complex double\n", n);
        return err;
    }

    err = vispo_fft_alloc_d(&inverse, n);
    if (err) {
        printf("Error allocating %d complex double\n", n);
        return err;
    }

    vispo_fft_setup_d(&forward, 0);
    vispo_fft_setup_d(&inverse, 1);

    memset(spec1_d, 0, m * n * sizeof(complex_d));
    memset(spec2_d, 0, m * n * sizeof(complex_d));

    for (i = 0; i < m * n; i += 1) {
        x = filedata[i];
        data_d[i].re = x;
        data_d[i].im = 0;
        real_data_d[i] = x;
    }

    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t0);
    off = 0;
    for (i = 0; i < m; i += 1) {
        vispo_fft_complex_d(&forward, spec1_d + off, data_d + off);
        off += n;
    }
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t1);
    dt_complex_d = time_diff(&t0, &t1);

    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t0);
    off = 0;
    for (i = 0; i < m; i += 1) {
        vispo_fft_real_d(&forward, spec2_d + off, real_data_d + off);
        off += n;
    }
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t1);
    dt_real_d = time_diff(&t0, &t1);

    memset(data2_d, 0, m * n * sizeof(complex_d));

    off = 0;
    for (i = 0; i < m; i += 1) {
        vispo_fft_complex_d(&inverse, data2_d + off, spec2_d + off);
        off += n;
    }

    for (i = 0; i < m * n; i += 1) {
        l = round(data2_d[i].re * r_n);
        if (l != filedata[i]) {
            puts("inverse FFT (double) result does not match file data");
            pass = 0;
            break;
        }
    }

    if (compare_arrays_d(spec1_d, spec2_d, n, eps_d) != n) {
        puts("real and complex FFT (double) results do not match");
        pass = 0;
    }

    if (pass)
        puts("PASS");
    else
        puts("FAIL");

    dt_complex_d *= r_m;
    dt_real_d *= r_m;

    printf("transform time, complex / real: = %.*f / %.*f s\n",
            digits, dt_complex_d, digits, dt_real_d);

    vispo_fft_free_d(&forward);
    vispo_fft_free_d(&inverse);

    return pass;
}

int real_test_s(short *filedata, int m, int n, int show)
{
    struct timespec t0, t1;
    double r_m = 1.0 / m;
    double r_n = 1.0 / n;
    double eps_s = log(n) / log(2);
    double dt_complex_s = 0;
    double dt_real_s = 0;
    int pass = 1;
    float x, y;
    int off, i, j, k;
    long l;
    struct vispo_fft_s forward, inverse;
    int err;

    err = vispo_fft_alloc_s(&forward, n);
    if (err)
        return -1;

    err = vispo_fft_alloc_s(&inverse, n);
    if (err)
        return -1;

    vispo_fft_setup_s(&forward, 0);
    vispo_fft_setup_s(&inverse, 1);

    for (i = 0; i < n; i += 1) {
        x = filedata[i];
        data_s[i].re = x;
        data_s[i].im = 0;
        real_data_s[i] = x;
    }

    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t0);
    off = 0;
    for (i = 0; i < m; i += 1) {
        vispo_fft_complex_s(&forward, spec1_s + off, data_s + off);
        off += n;
    }
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t1);
    dt_complex_s = time_diff(&t0, &t1);

    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t0);
    off = 0;
    for (i = 0; i < m; i += 1) {
        vispo_fft_real_s(&forward, spec2_s, real_data_s);
    }
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t1);
    dt_real_s = time_diff(&t0, &t1);

    memset(data2_s, 0, n * sizeof(complex_s));

    off = 0;
    for (i = 0; i < m; i += 1) {
        vispo_fft_complex_s(&inverse, data2_s, spec2_s);
    }

    for (i = 0; i < n; i += 1) {
        l = round(data2_s[i].re * r_n);
        if (l != filedata[i]) {
            puts("inverse FFT (single) result does not match file data");
            pass = 0;
            break;
        }
    }

    if (compare_arrays_s(spec1_s, spec2_s, n, eps_s) != n) {
        puts("Real and complex FFT results do not match (single)");
        pass = 0;
    }

    if (pass)
        puts("PASS");
    else
        puts("FAIL");

    dt_complex_s *= r_m;
    dt_real_s *= r_m;

    printf("transform time, complex / real: = %.*f / %.*f s\n",
            digits, dt_complex_s, digits, dt_real_s);

    vispo_fft_free_s(&forward);
    vispo_fft_free_s(&inverse);

    return pass;
}

int main(int argc, char **argv)
{
    struct vispo_fft_s forward_s, inverse_s;
    struct vispo_fft_d forward_d, inverse_d;
    short *filedata;
    struct timespec t0, t1;
    double r_n;
    double a, da;
    double e_t, e_f;
    double eps = 1, eps_d, eps_s;
    double x, y;
    int pass;
    int m, n;
    int i, k, c;
    long l;
    int err;
/*
    vispo_init();
*/
    filedata = alloc_large(sizeof(short) * 2 * DFT_MAX);
    vispo_fft_alloc_d(&forward_d, DFT_MAX);
    vispo_fft_alloc_d(&inverse_d, DFT_MAX);
    vispo_fft_alloc_s(&forward_s, FFT_MAX_SINGLE);
    vispo_fft_alloc_s(&inverse_s, FFT_MAX_SINGLE);

    real_data_d = alloc_large(sizeof(double) * DFT_MAX);
    real_data_s = alloc_large(sizeof(single) * FFT_MAX_SINGLE);

    data_d = alloc_large(sizeof(complex_d) * DFT_MAX);
    data_s = alloc_large(sizeof(complex_s) * FFT_MAX_SINGLE);

    data2_d = alloc_large(sizeof(complex_d) * DFT_MAX);
    data2_s = alloc_large(sizeof(complex_s) * FFT_MAX_SINGLE);

    spec1_d = alloc_large(sizeof(complex_d) * DFT_MAX);
    spec1_s = alloc_large(sizeof(complex_s) * FFT_MAX_SINGLE);

    spec2_d = alloc_large(sizeof(complex_d) * DFT_MAX);
    spec2_s = alloc_large(sizeof(complex_s) * FFT_MAX_SINGLE);

    clock_getres(CLOCK_MONOTONIC, &t0);
    srand(t0.tv_sec * 1000000000ul + t0.tv_nsec);

    /* timer precision in decimal digits */

    clock_getres(CLOCK_PROCESS_CPUTIME_ID, &t1);
    digits = ceil(log10(1.0 / (t1.tv_sec + t1.tv_nsec * 1e-9)));

    fesetround(FE_TONEAREST);

    /*
     *  Random complex data test
     */

    pass = 1;
    eps = 1;

    for (n = 2; n <= 4096; n <<= 1) {
        r_n = 1.0 / n;
        eps_d = log(n) / (log(2) * 256);
        eps_s = log(n) / log(2);

        err = vispo_fft_alloc_d(&forward_d, n);
        if (err) {
            printf("Error allocating %d complex double\n", n);
            return err;
        }

        err = vispo_fft_alloc_d(&inverse_d, n);
        if (err) {
            printf("Error allocating %d complex double\n", n);
            return err;
        }

        vispo_fft_setup_d(&forward_d, 0);
        vispo_fft_setup_d(&inverse_d, 1);

        if (n <= FFT_MAX_SINGLE) {
            vispo_fft_alloc_s(&forward_s, n);
            vispo_fft_alloc_s(&inverse_s, n);
            vispo_fft_setup_s(&forward_s, 0);
            vispo_fft_setup_s(&inverse_s, 1);
        }

        for (c = 0; c * n < 4096; c += 1) {
            memset(spec1_d, 0, n * sizeof(complex_d));
            memset(spec2_d, 0, n * sizeof(complex_d));

            for (i = 0; i < n * 2; i += 1) {
                filedata[i] = rand();
            }

            e_t = 0;
            e_f = 0;
            for (i = 0; i < n; i += 1) {
                x = filedata[i * 2];
                y = filedata[i * 2 + 1];
                data_d[i].re = x;
                data_d[i].im = y;
                if (n <= FFT_MAX_SINGLE) {
                    data_s[i].re = x;
                    data_s[i].im = y;
                }
            }

            e_t = signal_power_d(data_d, n);

            vispo_dft_complex_d(forward_d.tables + n, spec1_d, data_d, n);
            vispo_fft_complex_d(&forward_d, spec2_d, data_d);
            vispo_fft_complex_d(&inverse_d, data2_d, spec2_d);

            if (compare_arrays_d(spec1_d, spec2_d, n, eps_d) != n) {
                puts("FFT and DFT results (double) do not match");
                pass = 0;
                break;
            }

            /* scale the forward-inverse conversion result */

            for (i = 0; i < n; i += 1) {
                data2_d[i].re *= r_n;
                data2_d[i].im *= r_n;
            }

            if (compare_arrays_d(data_d, data2_d, n, eps_d) != n) {
                puts("iFFT output (double) does not match input signal");
                pass = 0;
                break;
            }

            e_f = signal_power_d(spec2_d, n);

            /*
             *  Parseval's relation check
             */
            e_f *= r_n;
            if (fabs(e_t - e_f) > eps) {
                printf("power is different in time and frequency domain (n=%d): %lf, %lf\n", n, e_t, e_f);
                for (i = 0; i < n; i += 1) {
                    printf("  %d: (%lf, %lf) (%lf, %lf)\n", i,
                        data_d[i].re, data_d[i].im, spec2_d[i].re, spec2_d[i].im);
                }
                pass = 0;
            }

            if (n > FFT_MAX_SINGLE)
                continue;

            vispo_dft_complex_s(forward_s.tables + n, spec1_s, data_s, n);
            vispo_fft_complex_s(&forward_s, spec2_s, data_s);
            vispo_fft_complex_s(&inverse_s, data2_s, spec2_s);
            if (compare_arrays_s(spec1_s, spec2_s, n, eps_s) != n) {
                puts("FFT and DFT results (single) do not match");
                pass = 0;
                break;
            }

            /* scale the forward-inverse conversion result */

            for (i = 0; i < n; i += 1) {
                data2_s[i].re *= r_n;
                data2_s[i].im *= r_n;
            }

            if (compare_arrays_s(data_s, data2_s, n, eps_s) != n) {
                puts("iFFT output (single) does not match input signal");
                pass = 0;
                break;
            }
        }

        if (n <= FFT_MAX_SINGLE) {
            vispo_fft_free_s(&forward_s);
            vispo_fft_free_s(&inverse_s);
        }

        vispo_fft_free_d(&forward_d);
        vispo_fft_free_d(&inverse_d);

        printf("Random complex test, %d points: %s\n", n, pass ? "PASS" : "FAIL");

        if (!pass) {
            break;
        }
    }

    if (!pass) {
        return 1;
    }

    /*
     *  Random real data test
     */

    puts("\n*** random real data test (single, n = 1 .. 1024)");

    for (n = 1; n <= 1024; n += 1) {
        r_n = 1.0 / n;
        eps_s = log(n) / log(2);

        vispo_fft_alloc_s(&forward_s, n);
        vispo_fft_alloc_s(&inverse_s, n);
        vispo_fft_setup_s(&forward_s, 0);
        vispo_fft_setup_s(&inverse_s, 1);

        for (i = 0; i < n; i += 1) {
            filedata[i] = rand();
        }

        for (i = 0; i < n; i += 1) {
            x = filedata[i];
            data_s[i].re = x;
            data_s[i].im = 0;
            real_data_s[i] = x;
        }

        vispo_dft_complex_s(forward_s.tables + n, spec1_s, data_s, n);
        vispo_dft_complex_s(inverse_s.tables + n, data2_s, spec1_s, n);

        for (i = 0; i < n; i += 1) {
            data2_s[i].re *= r_n;
            data2_s[i].im *= r_n;
        }

        if (compare_arrays_s(data_s, data2_s, n, eps_s * 2) != n) {
            printf("single (%d) input data and inverse DFT results do not match", n);
            pass = 0;
            break;
        }

        memset(data2_s, 0, sizeof(complex_s) * n);

        vispo_fft_complex_s(&forward_s, spec2_s, data_s);
        vispo_fft_complex_s(&inverse_s, data2_s, spec2_s);

        if (compare_arrays_s(spec1_s, spec2_s, n, eps_s) != n) {
            printf("single (%d) DFT and complex FFT results do not match", n);
            pass = 0;
            break;
        }

        for (i = 0; i < n; i += 1) {
            data2_s[i].re *= r_n;
            data2_s[i].im *= r_n;
        }

        if (compare_arrays_s(data_s, data2_s, n, eps_s * 2) != n) {
            printf("single (%d) input data and inverse FFT results do not match", n);
            pass = 0;
            break;
        }

        memset(spec2_s, 0, sizeof(complex_s) * n);

        vispo_fft_real_s(&forward_s, spec2_s, real_data_s);

        vispo_fft_free_s(&forward_s);
        vispo_fft_free_s(&inverse_s);

        if (compare_arrays_s(spec1_s, spec2_s, n, eps_s) != n) {
            printf("single (%d) DFT and real FFT results do not match", n);
            pass = 0;
            break;
        }
    }

    puts("\n*** random real data test (double, n = 1 .. 1024)");

    for (n = 1; n <= 1024; n += 1) {
        r_n = 1.0 / n;
        eps_d = log(n) / (log(2) * 256);

        vispo_fft_alloc_d(&forward_d, n);
        vispo_fft_alloc_d(&inverse_d, n);
        vispo_fft_setup_d(&forward_d, 0);
        vispo_fft_setup_d(&inverse_d, 1);

        for (i = 0; i < n; i += 1) {
            filedata[i] = rand();
        }

        for (i = 0; i < n; i += 1) {
            x = filedata[i];
            data_d[i].re = x;
            data_d[i].im = 0;
            real_data_d[i] = x;
        }

        vispo_dft_complex_d(forward_d.tables + n, spec1_d, data_d, n);
        vispo_dft_complex_d(inverse_d.tables + n, data2_d, spec1_d, n);

        for (i = 0; i < n; i += 1) {
            data2_d[i].re *= r_n;
            data2_d[i].im *= r_n;
        }

        if (compare_arrays_d(data_d, data2_d, n, eps_d * 2) != n) {
            printf("double (%d) input data and inverse DFT results do not match", n);
            pass = 0;
            break;
        }

        memset(data2_d, 0, sizeof(complex_d) * n);

        vispo_fft_complex_d(&forward_d, spec2_d, data_d);
        vispo_fft_complex_d(&inverse_d, data2_d, spec2_d);

        if (compare_arrays_d(spec1_d, spec2_d, n, eps_d) != n) {
            printf("double (%d) DFT and complex FFT results do not match", n);
            pass = 0;
            break;
        }

        for (i = 0; i < n; i += 1) {
            data2_d[i].re *= r_n;
            data2_d[i].im *= r_n;
        }

        if (compare_arrays_d(data_d, data2_d, n, eps_d * 2) != n) {
            printf("double (%d) input data and inverse FFT results do not match", n);
            pass = 0;
            break;
        }

        memset(spec2_d, 0, sizeof(complex_d) * n);

        vispo_fft_real_d(&forward_d, spec2_d, real_data_d);

        vispo_fft_free_d(&forward_d);
        vispo_fft_free_d(&inverse_d);

        if (compare_arrays_d(spec1_d, spec2_d, n, eps_d) != n) {
            printf("double (%d) DFT and real FFT results do not match", n);
            pass = 0;
            break;
        }
    }

    puts("\n *** real signal, single precision ***\n\n");

    m = FFT_MAX_SINGLE / 256;
    n = 256;
    while (m != 0) {
        printf("Random real test (single), %d points\n", n);
        for (i = 0; i < m * n; i += 1) {
            filedata[i] = rand();
        }

        real_test_s(filedata, m, n, 1);

        m >>= 1;
        n <<= 1;
    }

    puts("\n *** real signal, double precision ***\n\n");

    m = DFT_MAX / 256;
    n = 256;
    while (m != 0) {
        printf("Random real test (double), %d points\n", n);

        for (i = 0; i < m * n; i += 1) {
            filedata[i] = rand();
        }

        real_test_d(filedata, m, n, 1);

        m >>= 1;
        n <<= 1;
    }

    return 0;
}


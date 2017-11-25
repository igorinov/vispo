package com.example.project;

/**
 * Created by igorinov on 11/22/16.
 */
public class FourierTransform {
    double ce[];
    int size = 0;

    public static native int fft(double[] ce, double[] out, double[] signal, int size);
    public static native int fftReal(double[] ce, double[] out, double[] signal, int size);

    public int fft(double[] out, double[] in, int size) {
        return fft(ce, out, in, size);
    }

    public int fft_real(double[] out, double[] in, int size) {
        return fftReal(ce, out, in, size);
    }

    FourierTransform(int n, boolean inverse) {
        double da = 2 * Math.PI / n;
        double a;
        int sign = inverse ? 1 : -1;
        int i;

        size = n;
        ce = new double[size * 2];

        for (i = 0; i < size; i += 1) {
            a = i * da;
            ce[i * 2] = StrictMath.cos(a);
            ce[i * 2 + 1] = StrictMath.sin(a) * sign;
        }
    }

    static {
        try {
            System.loadLibrary("pvft");
        } catch (UnsatisfiedLinkError e) {
            e.printStackTrace();
        }
    }
}

//
// Java interface for PVFT
//

#include <jni.h>

#include "pvft.h"

#define CLASS_NAME FourierTransform
#define CAT_NAMES(a, b, c, d) a ## _ ## b ## _ ## c ## _ ## d
#define JAVA_NAME(p, c, n) CAT_NAMES(Java, p, c, n)
#define FN_NAME(method) JAVA_NAME(PACKAGE_NAME, CLASS_NAME, method)

extern "C" JNIEXPORT jint JNICALL FN_NAME(fftSetupS)(JNIEnv *env, jclass jc,
    jfloatArray jce, jint size, jboolean inverse)
{
    float *ce  = env->GetFloatArrayElements(jce, 0);

    fft_setup_s((complex_s *) ce, size, inverse);

    env->ReleaseFloatArrayElements(jce, ce, 0);

    return 0;
}

extern "C" JNIEXPORT jint JNICALL FN_NAME(fftSetupD)(JNIEnv *env, jclass jc,
    jdoubleArray jce, jint size, jboolean inverse)
{
    double *ce  = env->GetDoubleArrayElements(jce, 0);

    fft_setup_d((complex_d *) ce, size, inverse);

    env->ReleaseDoubleArrayElements(jce, ce, 0);

    return 0;
}

extern "C" JNIEXPORT jint JNICALL FN_NAME(fftS)(JNIEnv *env, jclass jc,
    jfloatArray jcs, jfloatArray jout, jfloatArray jin, jint size)
{
    float *cs  = env->GetFloatArrayElements(jcs, 0);
    float *in  = env->GetFloatArrayElements(jin, 0);
    float *out  = env->GetFloatArrayElements(jout, 0);

    fft_complex_d((complex_d *) cs, (complex_d *) out, (complex_d *) in, size);

    env->ReleaseFloatArrayElements(jcs, cs, JNI_ABORT);  /* read only; not copied back */
    env->ReleaseFloatArrayElements(jin, in, JNI_ABORT);  /* read only; not copied back */
    env->ReleaseFloatArrayElements(jout, out, 0);

    return 0;
}

extern "C" JNIEXPORT jint JNICALL FN_NAME(fftD)(JNIEnv *env, jclass jc,
    jdoubleArray jce, jdoubleArray jout, jdoubleArray jin, jint size)
{
    double *ce  = env->GetDoubleArrayElements(jce, 0);
    double *in  = env->GetDoubleArrayElements(jin, 0);
    double *out  = env->GetDoubleArrayElements(jout, 0);

    fft_complex_d((complex_d *) ce, (complex_d *) out, (complex_d *) in, size);

    env->ReleaseDoubleArrayElements(jce, ce, JNI_ABORT);  /* read only; not copied back */
    env->ReleaseDoubleArrayElements(jin, in, JNI_ABORT);  /* read only; not copied back */
    env->ReleaseDoubleArrayElements(jout, out, 0);

    return 0;
}

extern "C" JNIEXPORT jint JNICALL FN_NAME(fftRealS)(JNIEnv *env, jclass jc,
    jfloatArray jce, jfloatArray jout, jfloatArray jin, jint size)
{
    float *ce  = env->GetFloatArrayElements(jce, 0);
    float *in  = env->GetFloatArrayElements(jin, 0);
    float *out  = env->GetFloatArrayElements(jout, 0);

    fft_real_s((complex_s *) ce, (complex_s *) out, in, size);

    env->ReleaseFloatArrayElements(jce, ce, JNI_ABORT);  /* read only; not copied back */
    env->ReleaseFloatArrayElements(jin, in, JNI_ABORT);  /* read only; not copied back */
    env->ReleaseFloatArrayElements(jout, out, 0);

    return 0;
}

extern "C" JNIEXPORT jint JNICALL FN_NAME(fftRealD)(JNIEnv *env, jclass jc,
    jdoubleArray jce, jdoubleArray jout, jdoubleArray jin, jint size)
{
    double *ce  = env->GetDoubleArrayElements(jce, 0);
    double *in  = env->GetDoubleArrayElements(jin, 0);
    double *out  = env->GetDoubleArrayElements(jout, 0);

    fft_real_d((complex_d *) ce, (complex_d *) out, in, size);

    env->ReleaseDoubleArrayElements(jce, ce, JNI_ABORT);  /* read only; not copied back */
    env->ReleaseDoubleArrayElements(jin, in, JNI_ABORT);  /* read only; not copied back */
    env->ReleaseDoubleArrayElements(jout, out, 0);

    return 0;
}

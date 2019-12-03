LOCAL_PATH := $(call my-dir)
include $(CLEAR_VARS)
LOCAL_MODULE := pvft
LOCAL_SRC_FILES := vispo-jni.cpp dft.c fft.c
LOCAL_CPP_EXTENSION := .cxx .cpp .cc

# replace with your project name
LOCAL_CFLAGS += -DPACKAGE_NAME=com_example_project

ifeq ($(TARGET_ARCH_ABI),armeabi-v7a)
LOCAL_CFLAGS += -DASSEMBLY_FFT
LOCAL_SRC_FILES += armeabi-v7a/combine.S
endif

ifeq ($(TARGET_ARCH_ABI),arm64-v8a)
LOCAL_CFLAGS += -DASSEMBLY_FFT
LOCAL_SRC_FILES += arm64-v8a/combine.S
endif

ifeq ($(TARGET_ARCH_ABI),x86)
LOCAL_CFLAGS += -DASSEMBLY_FFT
LOCAL_SRC_FILES += x86/combine.S
endif

ifeq ($(TARGET_ARCH_ABI),x86_64)
LOCAL_CFLAGS += -DASSEMBLY_FFT
LOCAL_SRC_FILES += x86_64/combine.S
endif

TARGET_ARCH := all
include $(BUILD_SHARED_LIBRARY)


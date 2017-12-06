VERSION_MAJOR = 1
VERSION = 1.0.0
PROJECT_NAME=pvft
PREFIX=/usr/local

LIB_NAME = lib$(PROJECT_NAME).so
SONAME = lib$(PROJECT_NAME).so.$(VERSION_MAJOR)

HEADERS = pvft.h

OBJECTS = pvft.o

TARGET_ARCH := $(shell uname -m)

CFLAGS = -fPIC -O2

ifeq ($(TARGET_ARCH),armeabi-v7a)
OBJECTS += $(TARGET_ARCH)/combine.o
CFLAGS += -DDASSEMBLY_FFT
CFLAGS += -DBITREV
CFLAGS += -mfpu=neon
endif

ifeq ($(TARGET_ARCH),arm64-v8a)
OBJECTS += $(TARGET_ARCH)/combine.o
CFLAGS += -DDASSEMBLY_FFT
CFLAGS += -DBITREV
endif

ifeq ($(TARGET_ARCH),x86)
OBJECTS += $(TARGET_ARCH)/combine.o
CFLAGS += -DASSEMBLY_FFT
CFLAGS += -m32
SFLAGS += -m32
endif

ifeq ($(TARGET_ARCH),x86_64)
OBJECTS += $(TARGET_ARCH)/combine.o
CFLAGS += -DASSEMBLY_FFT
CFLAGS += -m64
SFLAGS += -m64
endif

LIBS = -lm

# CC = aarch64-linux-gnu-gcc

TARGETS = $(SONAME) $(LIB_NAME)

$(SONAME): $(LIB_NAME)
	ln -sf $(LIB_NAME) $(SONAME)

all: $(TARGETS)

install: $(LIB_NAME)
	cp $(LIB_NAME) $(PREFIX)/lib/$(LIB_NAME).$(VERSION)
	ln -sf $(LIB_NAME).$(VERSION) $(PREFIX)/lib/$(LIB_NAME)
	ldconfig

$(LIB_NAME): $(OBJECTS)
	$(CC) -shared -Wl,-soname,$(SONAME) -o $@ $(OBJECTS) $(LIBS)

pvft.o: pvft.c
	$(CC) -c -o $@ $(CFLAGS) $(LIBS) pvft.c

$(TARGET_ARCH)/%.o: $(TARGET_ARCH)/%.S
	$(CC) -c -o $@ $(SFLAGS) $<

test: test.c $(OBJECTS)
	$(CC) -o $@ $(CFLAGS) test.c -l$(PROJECT_NAME) $(LIBS)

clean:
	rm -f $(TARGETS) $(OBJECTS)


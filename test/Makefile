all: test

CFLAGS += -I..

test: test.c
	$(CC) $(CFLAGS) -o $@ test.c ../fft.o ../dft.o ../combine.o ../alloc_pages.o -lm

clean:
	rm -f test


all: test

CFLAGS += -I..

test: test.c
	$(CC) $(CFLAGS) -o $@ test.c ../pvft.o ../combine.o -lm

clean:
	rm -f test


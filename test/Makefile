all: test

CFLAGS += -I..

test: test.c
	$(CC) $(CFLAGS) -o $@ test.c -lm -lpvft


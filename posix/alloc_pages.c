#include <stdlib.h>
#include <sys/mman.h>

void *vispo_alloc_pages(int length)
{
    int prot = PROT_READ | PROT_WRITE;
    int flags = MAP_ANONYMOUS | MAP_PRIVATE;
    void *p;

    p = mmap(NULL, length, prot, flags, -1, 0);
    if (p == MAP_FAILED)
        return NULL;

    return p;
}

int vispo_free_pages(void *address, int length)
{
    return munmap(address, length);
}


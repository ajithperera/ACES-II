
#include <stdio.h>
#include <stdlib.h>

void c_heap_core(void * icore)
{
    fflush(stdout);
    printf(" ==> C malloc vs. FORTRAN common // test\n\n");
    {
        void * heap = malloc(1024*1024*200); /* 200MB heap */
        size_t adr_icore = (size_t)icore;
        size_t adr_heap  = (size_t)heap;
        size_t offset    = adr_heap-adr_icore;
        long lloff = offset/sizeof(long long);
        long loff = offset/sizeof(long);
        int  ioff = offset/sizeof(int);
        printf("blank common is at %p\n"
               "heap         is at %p\n"
               "size_t offset   is %li bytes\n"
               "llong  offset   is %li long longs\n"
               "long   offset   is %li longs\n"
               "int    offset   is %i ints\n\n",
               icore,heap,offset,lloff,loff,ioff);
        free(heap);
    }
    fflush(stdout);
    return;
}

void c_heap (void * icore)
{ c_heap_core(icore); return; }

void c_heap_(void * icore)
{ c_heap_core(icore); return; }

void C_HEAP (void * icore)
{ c_heap_core(icore); return; }


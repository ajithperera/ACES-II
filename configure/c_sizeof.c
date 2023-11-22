
#include <stdio.h>

void c_sizeof_core()
{
    fflush(stdout);
    printf(" ==> C data type sizes test\n\n"
           " char   is %i bytes\n"
           " int    is %i bytes\n"
           " long   is %i bytes\n"
           " llong  is %i bytes\n"
           " float  is %i bytes\n"
           " double is %i bytes\n"
           " size_t is %i bytes\n"
           " void * is %i bytes\n\n",
           sizeof(char),
           sizeof(int),sizeof(long),sizeof(long long),
           sizeof(float),sizeof(double),
           sizeof(size_t),sizeof(void*)
          );
    if (sizeof(void*) > sizeof(int))
       printf(" > WARNING: ints cannot hold memory addresses\n");
    if (sizeof(void*) > sizeof(long))
       printf(" > WARNING: longs cannot hold memory addresses\n");
    if (sizeof(void*) > sizeof(long long))
       printf(" > WARNING: llongs cannot hold memory addresses\n");
    if (sizeof(void*) == 8)
       printf(
     " > WARNING: If the default FORTRAN integer is not 8 bytes (see above),\n"
     " >          then F_ADR must be explicitly defined as integer*8.\n\n");
    fflush(stdout);
    return;
}

void c_sizeof ()
{ c_sizeof_core(); return; }

void c_sizeof_()
{ c_sizeof_core(); return; }

void C_SIZEOF ()
{ c_sizeof_core(); return; }


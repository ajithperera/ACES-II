
#include <stdio.h>

void c_endian_core(unsigned char * c)
{
    fflush(stdout);
    printf(" ==> endian test\n\n");
    if (*c) printf(" Little endian\n\n");
    else    printf(" Big endian\n\n");
    fflush(stdout);
    return;
}

void c_endian ()
{ int i=1; c_endian_core((unsigned char *)&i); return; }

void c_endian_()
{ int i=1; c_endian_core((unsigned char *)&i); return; }

void C_ENDIAN ()
{ int i=1; c_endian_core((unsigned char *)&i); return; }


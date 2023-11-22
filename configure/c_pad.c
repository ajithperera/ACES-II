
#include <stdio.h>

struct { double d1, d2;
         long long L1, L2;
         long l1, l2;
         int i1, i2;
         short s1, s2;
         char c1, c2;
         void *p; } st;

void c_pad_core()
{
    void *d1 = &st.d1, *d2 = &st.d2;
    void *L1 = &st.L1, *L2 = &st.L2;
    void *l1 = &st.l1, *l2 = &st.l2;
    void *i1 = &st.i1, *i2 = &st.i2;
    void *s1 = &st.s1, *s2 = &st.s2;
    void *c1 = &st.c1, *c2 = &st.c2;
    void *p = &st.p;
    fflush(stdout);
    printf(" ==> C padding test\n\n");
    printf("%s are %i bytes apart\n","doubles",d2-d1);
    printf("%s are %i bytes apart\n","llongs ",L2-L1);
    printf("%s are %i bytes apart\n","longs  ",l2-l1);
    printf("%s are %i bytes apart\n","ints   ",i2-i1);
    printf("%s are %i bytes apart\n","shorts ",s2-s1);
    printf("%s are %i bytes apart\n","chars  ",c2-c1);
    printf("%s are padded by %i bytes\n","doubles",L1-d2-sizeof(double));
    printf("%s are padded by %i bytes\n","llongs ",l1-L2-sizeof(long long));
    printf("%s are padded by %i bytes\n","longs  ",i1-l2-sizeof(long));
    printf("%s are padded by %i bytes\n","ints   ",s1-i2-sizeof(int));
    printf("%s are padded by %i bytes\n","shorts ",c1-s2-sizeof(short));
    printf("%s are padded by %i bytes\n","chars  ",p-c2-sizeof(char));
    printf("\n");
    return;
}

void c_pad()
{ c_pad_core(); return; }

void c_pad_()
{ c_pad_core(); return; }

void C_PAD()
{ c_pad_core(); return; }


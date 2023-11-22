
#include <stdio.h>

struct { double d; } ftncom;
struct { double d; } ftncom_;
struct { double d; } FTNCOM;

void common_set()
{
  ftncom.d = 1.;
  ftncom_.d = 2.;
  FTNCOM.d = 3.;
}

void c_sub()
{ common_set();
  fflush(stdout);
  printf(" C function names are direct translation.\n"
         " > Remove '-DC_SUFFIX -DC_UPPER'\n"
        );
  fflush(stdout);
  return;
}

void c_sub_()
{ common_set();
  fflush(stdout);
  printf(" C function names are suffixed.\n"
         " > Add '-DC_SUFFIX', remove '-DC_UPPER'\n"
        );
  fflush(stdout);
  return;
}

void C_SUB()
{ common_set();
  fflush(stdout);
  printf(" C function names are uppercase.\n"
         " > Add '-DC_UPPER', remove '-DC_SUFFIX'\n"
        );
  fflush(stdout);
  return;
}


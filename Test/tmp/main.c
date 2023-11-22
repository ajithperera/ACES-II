
#include <ieeefp.h>
#include <stdio.h>
#include <signal.h>

volatile sig_atomic_t code;

void handler_siginfo(int sig, siginfo_t *info, void *nix)
{ code = info->si_code; }

int
main(void)
{
  double  x, y;
  int     i, j;
  //fp_except_t fpm = fpsetmask( 0xFFFFFF );

  struct sigaction act;
  act.sa_sigaction = handler_siginfo;
  act.sa_flags |= SA_SIGINFO;
  //sigemptyset(&act.sa_mask); sigaction(SIGFPE, &act, NULL);

  code = -1;

  //x = 1e19; i=(int)x;
  //i/=j;
  //x/=j;
  //x/=y;
  i=255; j=i*i*i*i*i*i; printf("j = %d\n",j);
  fprintf(stderr, "Code: %d\n", code);

  return 0;
}


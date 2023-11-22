
#include <signal.h>

#include "f77_name.h"

void sighandler(int sig)
{
    /*call fortran handler*/
#define f_ncycread F77_NAME(ncycread,NCYCREAD)
    f_ncycread();

    /*install handler again*/
    signal(sig,sighandler);
}

/* to install the handler first time for SIGUSR1 */

void
F77_NAME(installsig,INSTALLSIG)
()
{
    signal(SIGUSR1,sighandler);
}


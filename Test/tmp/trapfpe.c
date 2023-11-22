
/*
 * This routine will turn on IEEE floating-point exception trapping
 * if there are no compiler flags to do it automatically (e.g., g77).
 */

#include <ieeefp.h>
#include <signal.h>
#include <stdio.h>
extern void fpe_handler();
extern void stupid_backtrace();

extern void
#ifdef C_SUFFIX
  trapfpe_
#else
  trapfpe
#endif /* C_SUFFIX */
()
{
    /*
     * set the FPE mask
     * FP_X_INV  : invalid operation
     * FP_X_DNML : denormal
     * FP_X_DZ   : zero divide
     * FP_X_OFL  : overflow
     * FP_X_UFL  : underflow <- BEWARE, this will kill almost any ACES prog!
     * FP_X_IMP  : (im)precision <- BEWARE, this will kill almost ANY  prog!
     * FP_X_STK  : stack fault
     */
    //fp_except_t i = fpsetmask( FP_X_INV | FP_X_OFL | FP_X_DZ );
    fp_except_t i = fpsetmask( 0xFF01 );

    /* change the FPE handler */
    struct sigaction sa_fpe;
    sa_fpe.sa_flags     = SA_SIGINFO; /* use sigaction instead of handler */
    sa_fpe.sa_sigaction = fpe_handler;
    sigemptyset(&sa_fpe.sa_mask);
    if (sigaction(SIGFPE,&sa_fpe,NULL) < 0)
    {
        printf("@TRAPFPE: Floating-point exceptions cannot be caught.\n");
    }
    else
    {
        printf("@TRAPFPE: Lots of floating-point exceptions will be caught.\n");
    }
    return;
}

/******************************************************************************/

void
fpe_handler(int sig, siginfo_t * sip, void * uap)
{
    char *label;
    switch (sip->si_code)
    {
        case FPE_INTOVF: label = "integer overflow";                 break;
        case FPE_INTDIV: label = "integer divide by zero";           break;
        case FPE_FLTDIV: label = "floating point divide by zero";    break;
        case FPE_FLTOVF: label = "floating point overflow";          break;
        case FPE_FLTUND: label = "floating point underflow";         break;
        case FPE_FLTRES: label = "floating point inexact result";    break;
        case FPE_FLTINV: label = "invalid floating point operation"; break;
        case FPE_FLTSUB: label = "subscript out of range";           break;
        default:         label = "???";                              break;
    }
    if ( (sip->si_code == FPE_INTOVF) ||
         (sip->si_code == FPE_INTDIV) ||
         (sip->si_code == FPE_FLTDIV) ||
         (sip->si_code == FPE_FLTOVF) ||
         (sip->si_code == FPE_FLTRES) ||
         (sip->si_code == FPE_FLTINV) ||
         (sip->si_code == FPE_FLTSUB)    )
    {
        fprintf(stderr,
                "A floating-point exception has been raised.\n"
                "      type:    %s (0x%x)\n"
                "      address: %p\n"
                "      signal:  %d\n\n",
                label, sip->si_code, sip->si_addr, sig);
        stupid_backtrace();
        abort();
    }
    else
    {
        return;
    }
}

/******************************************************************************/

#include <sys/types.h>
#include <unistd.h>

void
stupid_backtrace()
{
    char bt[200];
    char command[100];
    char prgname[100];
    char elink[100];
    pid_t pid = getpid();
    int lsize;
    FILE * gdb;

    sprintf(elink,"/proc/%i/file",pid);
    lsize = readlink(elink,prgname,99);
    if (lsize < 0)
    {
        fprintf(stderr,"\n@stupid_backtrace: readlink failed\n");
        exit(1);
    }
    prgname[lsize] = 0;

    sprintf(command,
            "echo backtrace | gdb %s %i 2>/dev/null",prgname,pid);

    gdb = popen(command,"r");
    do
    {
        fgets(bt,200,gdb);
    } while ( ! (strncmp(bt,"#",1) == 0) );
    do
    {
        printf(bt);
        fgets(bt,200,gdb);
    } while ( (strncmp(bt,"#",1) == 0) || (strncmp(bt,"  ",2) == 0) );
    pclose(gdb);
}


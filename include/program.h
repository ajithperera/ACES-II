
#ifndef _PROGRM_H_
#define _PRGRAM_H_ /*  DEFINES which program to use

 * The following definitions define the blas getrec, putrec, getlst,
 * putlst lroutines to call depending on whether a2proc is run
 * after CFOUR or ACES II calculation 
 *
 */

#ifdef _USE_CFOUR

#  define B_GETREC    Getrec_c4
#  define B_PUTREC    Putrec_c4
#  define B_GETLST    Getlst_c4
#  define B_PUTLST    Putlst_c4

#else _USE_ACES2

#  define B_GETREC    Getrec
#  define B_PUTREC    Putrec
#  define B_GETLST    Getlst
#  define B_PUTLST    Putlst

#endif

#endif /* _PROGRAM_H_ */


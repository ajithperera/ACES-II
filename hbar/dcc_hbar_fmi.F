      Subroutine dcc_hbar_fmi(Work, Length, Iuhf)

      Implicit Double Precision (A-H, O-Z)

      Double Precision Mone, One

      Dimension Work(Length)

#include "machsp.com"
#include "syminf.com"
#include "sympop.com"
#include "parcc.com"
C
C Construct the DCC FMI intermediates.
C
      Irrepx = 1
#ifdef _DEBUG_LVL0
      Write(6,"(a)") " Checksums of FMI : At Entry"
      Write(6,*)
      IHHA_LENGTH = IRPDPD(IRREPX,21)
      IHHB_LENGTH = IRPDPD(IRREPX,22)
      CALL Getlst(Work,1,1,1,9,91)
      Call Checksum("FMI", Work, IHHA_LENGTH)
      IF (IUHF.gt.0) Call Getlst(Work,1,1,1,10,91)
      IF (IUHF.gt.0) Call Checksum("Fmi", Work, IHHB_LENGTH)
#endif

      Call pdcc_quad2(Work,Length,iuhf,Fmi_scale,8,1)

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") " Checksums of FMI : At Exit"
      Write(6,*)
      IHHA_LENGTH = IRPDPD(IRREPX,21)
      IHHB_LENGTH = IRPDPD(IRREPX,22)
      CALL Getlst(Work,1,1,1,9,91)
      Call Checksum("FMI", Work, IHHA_LENGTH)
      IF (IUHF.gt.0) Call Getlst(Work,1,1,1,10,91)
      IF (IUHF.gt.0) Call Checksum("Fmi", Work, IHHB_LENGTH)
      Write(6,*)
#endif

      Return
      End

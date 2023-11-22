      Subroutine dcc_hbar_fae(Work, Length, Iuhf)

      Implicit Double Precision (A-H, O-Z)

      Double Precision Mone, One

      Dimension Work(Length)

#include "machsp.com"
#include "syminf.com"
#include "sympop.com"
#include "parcc.com"
C
C Construct the DCC FAE intermediate
C
      Irrepx = 1 
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") " Checksums of FAE: At Entry"
      Write(6,*)
      IPPA_LENGTH = IRPDPD(IRREPX,19)
      IPPB_LENGTH = IRPDPD(IRREPX,20)
      CALL Getlst(Work,1,1,1,9,92)
      Call Checksum("FEA", Work, IPPA_LENGTH)
      IF (IUHF.gt.0) Call Getlst(Work,1,1,1,10,92)
      IF (IUHF.gt.0) Call Checksum("Fea", Work, IPPB_LENGTH)
#endif
      Call pdcc_quad3(Work,Length,iuhf,Fae_scale,8,1)

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") " Checksums of FAE: At Exit"
      Write(6,*)
      IPPA_LENGTH = IRPDPD(IRREPX,19)
      IPPB_LENGTH = IRPDPD(IRREPX,20)
      CALL Getlst(Work,1,1,1,9,92)
      Call Checksum("FEA", Work, IPPA_LENGTH)
      IF (IUHF.gt.0) Call Getlst(Work,1,1,1,10,92)
      IF (IUHF.gt.0) Call Checksum("Fea", Work, IPPB_LENGTH)
#endif
      Return 
      End


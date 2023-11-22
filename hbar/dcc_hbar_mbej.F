      Subroutine dcc_hbar_mbej(Work, Length, Iuhf)

      Implicit Double Precision (A-H, O-Z)

      Integer AAAA_LENGTH_MBEJ,BBBB_LENGTH_MBEJ,AABB_LENGTH_MBEJ
      Integer BBAA_LENGTH_MBEJ,ABAB_LENGTH_MBEJ,BABA_LENGTH_MBEJ
      Double Precision Mone, One

      Dimension Work(Length)

#include "machsp.com"
#include "syminf.com"
#include "sympop.com"
#include "parcc.com"
C                             _     =    
C In lambda, for DCC, we need W and W  along  with the DCC W(mb,ej). 
C                                              -     =
C Instead of standard CCSD/CCD construction of W and W using, here
C W(mb,ej) intermediate ar more direct route is employed here. 
C
      IRREPX = 1
      Imode  = 0
C
C Recreate the pointers 117,118,123-126 lists. These are used in lambda.
C                         =
C These list contains the W  lists that will be used in EOM code.

      Call Inipck(1,9,9,123,Imode,0,1)
      Call Inipck(1,9,10,118,Imode,0,1)
      Call Inipck(1,11,11,125,Imode,0,1)

      If (Iuhf .ne. 0) Then
         Call Inipck(1,10,10,124,Imode,0,1)
         Call Inipck(1,10,9,117,Imode,0,1)
         Call Inipck(1,12,12,126,Imode,0,1)
      Endif

      IF (IUhf .GT. 0) Then
         AAAA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,123),ISYTYP(2,123))
         BBBB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,124),ISYTYP(2,124))
         AABB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,118),ISYTYP(2,118))
         BBAA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,117),ISYTYP(2,117))
         ABAB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,125),ISYTYP(2,125))
         BABA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,126),ISYTYP(2,126))

         MAX_LENGTH = MAX(AAAA_LENGTH_MBEJ, BBBB_LENGTH_MBEJ)
         MAX_LENGTH = MAX(MAX_LENGTH, AABB_LENGTH_MBEJ)
         MAX_LENGTH = MAX(MAX_LENGTH, BBAA_LENGTH_MBEJ)
         MAX_LENGTH = MAX(MAX_LENGTH, ABAB_LENGTH_MBEJ)
         MAX_LENGTH = MAX(MAX_LENGTH, ABAB_LENGTH_MBEJ)

         Call Zero(Work, MAX_LENGTH)
         Call Putall(Work, AAAA_LENGTH_MBEJ, IRREPX, 123)
         Call Putall(Work, BBBB_LENGTH_MBEJ, IRREPX, 124)
         Call Putall(Work, AABB_LENGTH_MBEJ, IRREPX, 118)
         Call Putall(Work, BBAA_LENGTH_MBEJ, IRREPX, 117)
         Call Putall(Work, ABAB_LENGTH_MBEJ, IRREPX, 125)
         Call Putall(Work, BABA_LENGTH_MBEJ, IRREPX, 126)
      Else
         AAAA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,123),ISYTYP(2,123))
         AABB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,118),ISYTYP(2,118))
         ABAB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,125),ISYTYP(2,125))

         MAX_LENGTH = MAX(AAAA_LENGTH_MBEJ, AABB_LENGTH_MBEJ)
         MAX_LENGTH = MAX(MAX_LENGTH, ABAB_LENGTH_MBEJ)

         Call Zero(Work, MAX_LENGTH)
         Call Putall(Work, AAAA_LENGTH_MBEJ, IRREPX, 123)
         Call Putall(Work, AABB_LENGTH_MBEJ, IRREPX, 118)
         Call Putall(Work, ABAB_LENGTH_MBEJ, IRREPX, 125)
      Endif
C               =
C Construct the W(mb,ej) and store 123-126,117 and 118. Currently 54-59 
C             -
C has the DCC W(mb,ej) and the lists 254-259 has the I(mb,ej). 
C     = 
C The W(mb,ej) = <mb||ej> - 1/2 I(mb,ej)
C
C The purpose of the following call is to put <mb||ej> integrals to 
C 123-126 and 117 and 118 lists. 

      C1 =  1.0
      C2 =  0.0
      C3 =  0.0
      C4 =  0.0

      Call Pdcc_dwmbej(Work,Length,"ABAB",Iuhf,C1,C2,C3,C4,118,117,
     &                 .False.,.False.,Wmbej_scale)
      Call Pdcc_dwmbej(Work,Length,"ABBA",Iuhf,C1,C2,C3,C4,125,125,
     &                 .False.,.False.,Wmbej_scale)

      If (Iuhf .Gt. 0) Then
          Call Pdcc_dwmbej(Work,Length,"AAAA",Iuhf,C1,C2,C3,C4,122,122,
     &                 .False.,.False.,Wmbej_scale)
          Call Pdcc_dwmbej(Work,Length,"BBBB",Iuhf,C1,C2,C3,C4,122,122,
     &                 .False.,.False.,Wmbej_scale)
          Call Pdcc_dwmbej(Work,Length,"BABA",Iuhf,C1,C2,C3,C4,117,118,
     &                 .False.,.False.,Wmbej_scale)
          Call Pdcc_dwmbej(Work,Length,"BAAB",Iuhf,C1,C2,C3,C4,126,126,
     &                 .False.,.False.,Wmbej_scale)
      Endif

      If (Iuhf .eq. 0) Call Pdcc_reset(Work,Length,Iuhf,118,125,123)
C 
C                                                = 
C Add the -I(mb,ej) to complete the formation of W(mb,ej)
C
      Mone = -1.0D0
      One  =  1.0D0

      Call Saxlst(Work,Length,254,123,123,Mone,One)
      Call Saxlst(Work,Length,258,125,125,Mone,One)

      If (IUHF.NE.0)THEN
         CALL Saxlst(Work,Length,256,118,118,Mone,One)
         CALL Saxlst(Work,Length,255,124,124,Mone,One)
         CALL Saxlst(Work,Length,257,117,117,MOne,One)
         CALL Saxlst(Work,Length,259,126,126,Mone,One)
      Else
         Call modf_quikab(Work,Length,Iuhf,123,125,118)
      Endif 

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,*)     "                =         "
      Write(6,"(a)") "The checksums of W at exit:"
      Write(6,*)
      IF (Iuhf .Gt. 0) Then
         AAAA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,123),ISYTYP(2,123))
         BBBB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,124),ISYTYP(2,124))
         AABB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,118),ISYTYP(2,118))
         BBAA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,117),ISYTYP(2,117))
         ABAB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,125),ISYTYP(2,125))
         BABA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,126),ISYTYP(2,126))

         Call Getall(Work, AAAA_LENGTH_MBEJ, IRREPX, 123)
         Call Checksum("AAAA_MBEJ", Work, AAAA_LENGTH_MBEJ)
         Call Getall(Work, BBBB_LENGTH_MBEJ, IRREPX, 124)
         Call Checksum("BBBB_MBEJ", Work, BBBB_LENGTH_MBEJ)
         Call Getall(Work, AABB_LENGTH_MBEJ, IRREPX, 118)
         Call Checksum("AABB_MBEJ", Work, AABB_LENGTH_MBEJ)
         Call Getall(Work, BBAA_LENGTH_MBEJ, IRREPX, 117)
         Call Checksum("BBAA_MBEJ", Work, BBAA_LENGTH_MBEJ)
         Call Getall(Work, ABAB_LENGTH_MBEJ, IRREPX, 125)
         Call Checksum("ABAB_MBEJ", Work, AABB_LENGTH_MBEJ)
         Call Getall(Work, BABA_LENGTH_MBEJ, IRREPX, 126)
         Call Checksum("BABA_MBEJ", Work, BABA_LENGTH_MBEJ)
         Write(6,*) 
      Else
         AAAA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,123),ISYTYP(2,123))
         AABB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,118),ISYTYP(2,118))
         ABAB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,125),ISYTYP(2,125))

         Call Getall(Work, AAAA_LENGTH_MBEJ, IRREPX, 123)
         Call Checksum("AAAA_MBEJ", Work, AAAA_LENGTH_MBEJ)
         Call Getall(Work, AABB_LENGTH_MBEJ, IRREPX, 118)
         Call Checksum("AABB_MBEJ", Work, AABB_LENGTH_MBEJ)
         Call Getall(Work, ABAB_LENGTH_MBEJ, IRREPX, 125)
         Call Checksum("ABAB_MBEJ", Work, AABB_LENGTH_MBEJ)
         Write(6,*) 
      Endif
#endif

      Return
      End
 

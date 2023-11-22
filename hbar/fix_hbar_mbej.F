      Subroutine Fix_HBAR_mbej(Work, Maxcor, Iuhf)

      Implicit Double Precision (A-H, O-Z)

      Integer AAAA_LENGTH_MBEJ, BBBB_LENGTH_MBEJ, AABB_LENGTH_MBEJ
      Integer BBAA_LENGTH_MBEJ, ABAB_LENGTH_MBEJ, BABA_LENGTH_MBEJ
      Integer AAAA_LENGTH_AIBJ, BBBB_LENGTH_AIBJ
      Integer AA_FAE_LENGTH, BB_FAE_LENGTH

      Dimension Work(Maxcor)
      Logical UHF

#include "machsp.com"
#include "sympop.com"
#include "sym.com"

      UHF    = .FALSE.
      UHF    = (IUhf .EQ. 1)
      IRREPX = 1

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") "The checksums of Hbar(me,bj) at Entry:" 
      Write(6,*)
      IF (UHF) Then
         AAAA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,54),ISYTYP(2,54))
         BBBB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,55),ISYTYP(2,55))
         AABB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,56),ISYTYP(2,56))
         BBAA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,57),ISYTYP(2,57))
         ABAB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,58),ISYTYP(2,58))
         BABA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,59),ISYTYP(2,59))
         Call Getall(Work, AAAA_LENGTH_MBEJ, IRREPX, 54)
         Call Checksum("AAAA_MBEJ", Work, AAAA_LENGTH_MBEJ)
         Call Getall(Work, BBBB_LENGTH_MBEJ, IRREPX, 55)
         Call Checksum("BBBB_MBEJ", Work, BBBB_LENGTH_MBEJ)
         Call Getall(Work, AABB_LENGTH_MBEJ, IRREPX, 56)
         Call Checksum("AABB_MBEJ", Work, AABB_LENGTH_MBEJ)
         Call Getall(Work, BBAA_LENGTH_MBEJ, IRREPX, 57)
         Call Checksum("BBAA_MBEJ", Work, BBAA_LENGTH_MBEJ)
         Call Getall(Work, ABAB_LENGTH_MBEJ, IRREPX, 58)
         Call Checksum("ABAB_MBEJ", Work, AABB_LENGTH_MBEJ)
         Call Getall(Work, BABA_LENGTH_MBEJ, IRREPX, 59)
         Call Checksum("BABA_MBEJ", Work, BABA_LENGTH_MBEJ)
      Else
         AAAA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,54),ISYTYP(2,54))
         AABB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,56),ISYTYP(2,56))
         ABAB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,58),ISYTYP(2,58))
         Call Getall(Work, AAAA_LENGTH_MBEJ, IRREPX, 54)
         Call Checksum("AAAA_MBEJ", Work, AAAA_LENGTH_MBEJ)
         Call Getall(Work, AABB_LENGTH_MBEJ, IRREPX, 56)
         Call Checksum("AABB_MBEJ", Work, AABB_LENGTH_MBEJ)
         Call Getall(Work, ABAB_LENGTH_MBEJ, IRREPX, 58)
         Call Checksum("ABAB_MBEJ", Work, AABB_LENGTH_MBEJ)
      Endif
#endif 

C First initialize the contents in W(mb,ej) list to zero
     
      IF (UHF) Then
         AAAA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,54),ISYTYP(2,54))
         BBBB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,55),ISYTYP(2,55))
         AABB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,56),ISYTYP(2,56))
         BBAA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,57),ISYTYP(2,57))
         ABAB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,58),ISYTYP(2,58))
         BABA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,59),ISYTYP(2,59))

         MAX_LENGTH = MAX(AAAA_LENGTH_MBEJ, BBBB_LENGTH_MBEJ)
         MAX_LENGTH = MAX(MAX_LENGTH, AABB_LENGTH_MBEJ)
         MAX_LENGTH = MAX(MAX_LENGTH, BBAA_LENGTH_MBEJ)
         MAX_LENGTH = MAX(MAX_LENGTH, ABAB_LENGTH_MBEJ)
         MAX_LENGTH = MAX(MAX_LENGTH, ABAB_LENGTH_MBEJ)
          
         Call Zero(Work, MAX_LENGTH)
         Call Putall(Work, AAAA_LENGTH_MBEJ, IRREPX, 54)
         Call Putall(Work, BBBB_LENGTH_MBEJ, IRREPX, 55)
         Call Putall(Work, AABB_LENGTH_MBEJ, IRREPX, 56)
         Call Putall(Work, BBAA_LENGTH_MBEJ, IRREPX, 57)
         Call Putall(Work, ABAB_LENGTH_MBEJ, IRREPX, 58)
         Call Putall(Work, BABA_LENGTH_MBEJ, IRREPX, 59)
      Else
         AAAA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,54),ISYTYP(2,54))
         AABB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,56),ISYTYP(2,56))
         ABAB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,58),ISYTYP(2,58))

         MAX_LENGTH = MAX(AAAA_LENGTH_MBEJ, AABB_LENGTH_MBEJ)
         MAX_LENGTH = MAX(MAX_LENGTH, ABAB_LENGTH_MBEJ)

         Call Zero(Work, MAX_LENGTH)
         Call Putall(Work, AAAA_LENGTH_MBEJ, IRREPX, 54)
         Call Putall(Work, AABB_LENGTH_MBEJ, IRREPX, 56)
         Call Putall(Work, ABAB_LENGTH_MBEJ, IRREPX, 58)
      Endif

#ifdef _DEBUG_LVLM
      Write(6,"(a)") "Checking the integral list 19 and 20"
      AAAA_LENGTH_AIBJ = IDSYMSZ(IRREPX,ISYTYP(1,19),ISYTYP(2,19))
      Call Getall(Work, AAAA_LENGTH_AIBJ, IRREPX, 19)
      Call Checksum("AAAA_AIBJ", Work, AAAA_LENGTH_AIBJ)
      If (UHF) THEN
         BBBB_LENGTH_AIBJ = IDSYMSZ(IRREPX,ISYTYP(1,20),ISYTYP(2,20))
         Call Getall(Work, BBBB_LENGTH_AIBJ, IRREPX, 20)
         Call Checksum("BBBB_AIBJ", Work, BBBB_LENGTH_AIBJ)
      Endif 
#endif 

      C1 =  1.0
      C2 =  0.0
      C3 =  1.0
      C4 = -0.5

      AA_FAE_LENGTH = IRPDPD(1,19)
      BB_FAE_LENGTH = 0
      IF (Iuhf .gt. 0) BB_FAE_LENGTH = IRPDPD(1,20)

      I000 = 1
      I010 = I000 + AA_FAE_LENGTH
      I020 = I010 + BB_FAE_LENGTH

      Call Getlst(Work(I000), 1, 1, 1, 1, 92)
      
      If (Iuhf .Gt. 0) Call Getlst(Work(I010), 1, 1, 1, 2, 92)

#ifdef _DEBUG_LVL0
      Call Checksum("AA_FAE", Work(I000), AA_FAE_LENGTH,S)
      If (Iuhf .gt. 0) Call Checksum("BB_FAE", Work(I010),
     &                                BB_FAE_LENGTH)
      Write(6,*)
#endif 

      Call T1ring(Work(I020),(Maxcor-I020),Iuhf,.False.)

      Call Putlst(Work(I000), 1, 1, 1, 1, 92)
      If (Iuhf .Gt. 0) Call Putlst(Work(I010), 1, 1, 1, 2, 92)

      IF(UHF) THEN
            CALL Dwmbej(Work,MAXCOR,'AAAA',IUHF,C1,C2,C3,C4)
            CALL Dwmbej(Work,MAXCOR,'BBBB',IUHF,C1,C2,C3,C4)
            CALL Dwmbej(Work,MAXCOR,'BABA',IUHF,C1,C2,C3,C4)
            CALL Dwmbej(Work,MAXCOR,'BAAB',IUHF,C1,C2,C3,C4)
      ENDIF
C
      CALL Dwmbej(Work,MAXCOR,'ABAB',IUHF,C1,C2,C3,C4)
      CALL Dwmbej(WOrk,MAXCOR,'ABBA',IUHF,C1,C2,C3,C4)

      If (.NOT. UHF) Call Reset(Work,Maxcor,Iuhf)
C
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(*,"(a,a)") "The checksums of Hbar(me,bj) at Exit: Should", 
     &                 " not change for the standard" 
      Write(6,"(a,a)") "CCSD calculations; these",
     &              " will change (supposed to) for DCC calculations."
      Write(6,*)
      
      IF (UHF) Then
         AAAA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,54),ISYTYP(2,54))
         BBBB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,55),ISYTYP(2,55))
         AABB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,56),ISYTYP(2,56))
         BBAA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,57),ISYTYP(2,57))
         ABAB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,58),ISYTYP(2,58))
         BABA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,59),ISYTYP(2,59))
         Call Getall(Work, AAAA_LENGTH_MBEJ, IRREPX, 54)
         Call Checksum("AAAA_MBEJ", Work, AAAA_LENGTH_MBEJ)
         Call Getall(Work, BBBB_LENGTH_MBEJ, IRREPX, 55)
         Call Checksum("BBBB_MBEJ", Work, BBBB_LENGTH_MBEJ)
         Call Getall(Work, AABB_LENGTH_MBEJ, IRREPX, 56)
         Call Checksum("AABB_MBEJ", Work, AABB_LENGTH_MBEJ)
         Call Getall(Work, BBAA_LENGTH_MBEJ, IRREPX, 57)
         Call Checksum("BBAA_MBEJ", Work, BBAA_LENGTH_MBEJ)
         Call Getall(Work, ABAB_LENGTH_MBEJ, IRREPX, 58)
         Call Checksum("ABAB_MBEJ", Work, AABB_LENGTH_MBEJ)
         Call Getall(Work, BABA_LENGTH_MBEJ, IRREPX, 59)
         Call Checksum("BABA_MBEJ", Work, BABA_LENGTH_MBEJ)
      Else
         AAAA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,54),ISYTYP(2,54))
         AABB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,56),ISYTYP(2,56))
         ABAB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,58),ISYTYP(2,58))
         Call Getall(Work, AAAA_LENGTH_MBEJ, IRREPX, 54)
         Call Checksum("AAAA_MBEJ", Work, AAAA_LENGTH_MBEJ)
         Call Getall(Work, AABB_LENGTH_MBEJ, IRREPX, 56)
         Call Checksum("AABB_MBEJ", Work, AABB_LENGTH_MBEJ)
         Call Getall(Work, ABAB_LENGTH_MBEJ, IRREPX, 58)
         Call Checksum("ABAB_MBEJ", Work, AABB_LENGTH_MBEJ)
      Endif
#endif 

      Return
      End

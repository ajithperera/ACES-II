#include "flags.h"
      Subroutine post_dcc_mods(Work, Length, Iuhf)

      Implicit Double Precision (A-H, O-Z)
      Integer AAAA_LENGTH_MBEJ,BBBB_LENGTH_MBEJ,AABB_LENGTH_MBEJ
      Integer BBAA_LENGTH_MBEJ,ABAB_LENGTH_MBEJ,BABA_LENGTH_MBEJ
      Integer AAAA_LENGTH_MNIJ, BBBB_LENGTH_MNIJ,ABAB_LENGTH_MNIJ
      Integer AA_FAE_LENGTH, BB_FAE_LENGTH

      Dimension Work(Length)

#include "machsp.com"
#include "syminf.com"
#include "sympop.com"
#include "parcc.com"
#include "flags.com"

      Imode  = 0
      Irrepx = 1

      Write(6,"(a,a)") "  Post DCC reconstruction of W(mn,ij) and",
     &                 " I(mb,ej) = t2(mn,ef)*<nj||fb>." 
      Write(6,*)

      Call Inipck(1,9,9,254,Imode,0,1)
      Call Inipck(1,9,10,256,Imode,0,1)
      Call Inipck(1,11,11,258,Imode,0,1)

      If (Iuhf .ne. 0) Then
         Call Inipck(1,10,10,255,Imode,0,1)
         Call Inipck(1,10,9,257,Imode,0,1)
         Call Inipck(1,12,12,259,Imode,0,1)
      Endif 
C
C If Coulomb integrals are used in W(mb,ej), we need to have two versions
C of W(mb,ej) (ie. one with the Coulomb only and the other with standard
C antisymmetrized integrals). 
C

      Call Inipck(1,9,9,264,Imode,0,1)
      Call Inipck(1,9,10,266,Imode,0,1)
      Call Inipck(1,11,11,268,Imode,0,1)

      Call Inipck(1,9,9,164,Imode,0,1)
      Call Inipck(1,9,10,166,Imode,0,1)
      Call Inipck(1,11,11,168,Imode,0,1)

      If (Iuhf .ne. 0) Then
          Call Inipck(1,10,10,265,Imode,0,1)
          Call Inipck(1,10,9,267,Imode,0,1)
          Call Inipck(1,12,12,269,Imode,0,1)

          Call Inipck(1,10,10,165,Imode,0,1)
          Call Inipck(1,10,9,167,Imode,0,1)
          Call Inipck(1,12,12,169,Imode,0,1)
      Endif 
C
      IHHA_LENGTH = IRPDPD(IRREPX,21)
      IHHB_LENGTH = IRPDPD(IRREPX,22)

      IPPA_LENGTH = IRPDPD(IRREPX,19)
      IPPB_LENGTH = IRPDPD(IRREPX,20)

      Call Updmoi(1,IHHA_LENGTH,9,91, 0, 0)
      Call Updmoi(1,IPPA_LENGTH,9,92, 0, 0)
      If (Iuhf .ne.0) Then 
         Call Updmoi(1,IHHB_LENGTH,10,91, 0, 0)
         Call Updmoi(1,IPPB_LENGTH,10,92, 0, 0)
      Endif 

C Construct the piece of W(mn,ij)= W(mn,ij)'+1/4 T2(ij,ef)W(mn||ef). This term
C is removed from the W(min,ij) intermediate. Notice that this piece 
C is constructed in genint.F and it is intact. 
C
      If (Iflags(h_IFLAGS_calc) .gt. 9) Then

      If (IUHF .gt. 0) Then
         AAAA_LENGTH_MNIJ = IDSYMSZ(IRREPX,ISYTYP(1,51),ISYTYP(2,51))
         BBBB_LENGTH_MNIJ = IDSYMSZ(IRREPX,ISYTYP(1,52),ISYTYP(2,52))
         I000 = 1
         I010 = I000 + AAAA_LENGTH_MNIJ
         I020 = I010 + AAAA_LENGTH_MNIJ

         Call Getall(Work(I000), AAAA_LENGTH_MNIJ, IRREPX, 51)
         Call Getall(Work(I010), AAAA_LENGTH_MNIJ, IRREPX, 251)
         Call Daxpy(AAAA_LENGTH_MNIJ, 1.0D0, Work(I000), 1, Work(I010),
     &              1)
         Call Putall(Work(I010), AAAA_LENGTH_MNIJ, IRREPX, 51)
         Call Putall(Work(I000), AAAA_LENGTH_MNIJ, IRREPX, 251)

         I000 = 1
         I010 = I000 + BBBB_LENGTH_MNIJ
         I020 = I010 + BBBB_LENGTH_MNIJ

         Call Getall(Work(I000), BBBB_LENGTH_MNIJ, IRREPX, 52)
         Call Getall(Work(I010), BBBB_LENGTH_MNIJ, IRREPX, 252)
         Call Daxpy(BBBB_LENGTH_MNIJ, 1.0D0, Work(I000), 1, Work(I010),
     &              1)
         Call Putall(Work(I010), BBBB_LENGTH_MNIJ, IRREPX, 52)
         Call Putall(Work(I000), BBBB_LENGTH_MNIJ, IRREPX, 252)

         ABAB_LENGTH_MNIJ = IDSYMSZ(IRREPX,ISYTYP(1,53),ISYTYP(2,53))
         I000 = 1
         I010 = I000 + ABAB_LENGTH_MNIJ
         I020 = I010 + ABAB_LENGTH_MNIJ

         Call Getall(Work(I000), ABAB_LENGTH_MNIJ, IRREPX, 53)
         Call Getall(Work(I010), ABAB_LENGTH_MNIJ, IRREPX, 253)
         Call Daxpy(ABAB_LENGTH_MNIJ, 1.0D0, Work(I000), 1, Work(I010),
     &              1)
         Call Putall(Work(I010), ABAB_LENGTH_MNIJ, IRREPX, 53)
         Call Putall(Work(I000), ABAB_LENGTH_MNIJ, IRREPX, 253)

      Else
         ABAB_LENGTH_MNIJ = IDSYMSZ(IRREPX,ISYTYP(1,53),ISYTYP(2,53))
         I000 = 1
         I010 = I000 + ABAB_LENGTH_MNIJ
         I020 = I010 + ABAB_LENGTH_MNIJ

         Call Getall(Work(I000), ABAB_LENGTH_MNIJ, IRREPX, 53)
         Call Getall(Work(I010), ABAB_LENGTH_MNIJ, IRREPX, 253)
         Call Daxpy(ABAB_LENGTH_MNIJ, 1.0D0, Work(I000), 1, Work(I010),
     &              1)
         Call Putall(Work(I010), ABAB_LENGTH_MNIJ, IRREPX, 53)
         Call Putall(Work(I000), ABAB_LENGTH_MNIJ, IRREPX, 253)
      Endif

      Else 

      If (IUHF .gt. 0) Then
         
         AAAA_LENGTH_MNIJ = IDSYMSZ(IRREPX,ISYTYP(1,51),ISYTYP(2,51))
         BBBB_LENGTH_MNIJ = IDSYMSZ(IRREPX,ISYTYP(1,52),ISYTYP(2,52))

         I000 = 1
         I010 = I000 + AAAA_LENGTH_MNIJ
         Call Getall(Work(I000), AAAA_LENGTH_MNIJ, IRREPX, 51)
         Call Putall(Work(I000), AAAA_LENGTH_MNIJ, IRREPX, 251)

         I000 = 1
         I010 = I000 + BBBB_LENGTH_MNIJ
         Call Getall(Work(I000), BBBB_LENGTH_MNIJ, IRREPX, 52)
         Call Putall(Work(I000), BBBB_LENGTH_MNIJ, IRREPX, 252)

         ABAB_LENGTH_MNIJ = IDSYMSZ(IRREPX,ISYTYP(1,53),ISYTYP(2,53))
         I000 = 1
         I010 = I000 + ABAB_LENGTH_MNIJ
         Call Getall(Work(I000), ABAB_LENGTH_MNIJ, IRREPX, 53)
         Call Putall(Work(I000), ABAB_LENGTH_MNIJ, IRREPX, 253)

      Else

         ABAB_LENGTH_MNIJ = IDSYMSZ(IRREPX,ISYTYP(1,53),ISYTYP(2,53))
         I000 = 1
         I010 = I000 + ABAB_LENGTH_MNIJ
         Call Getall(Work(I000), ABAB_LENGTH_MNIJ, IRREPX, 53)
         Call Putall(Work(I000), ABAB_LENGTH_MNIJ, IRREPX, 253)
      Endif 

      Endif 
C
C The W(mb,ej) intermediate at this point is CC (pCC) or DCC depending 
C on wether we did a CC or DCC calculation. 
C
      AAAA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,54),ISYTYP(2,54))
      AABB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,56),ISYTYP(2,56))
      ABAB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,58),ISYTYP(2,58))

      If (Iuhf .NE. 0) Then
         AAAA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,54),ISYTYP(2,54))
         BBBB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,55),ISYTYP(2,55))
         AABB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,56),ISYTYP(2,56))
         BBAA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,57),ISYTYP(2,57))
         ABAB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,58),ISYTYP(2,58))
         BABA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,59),ISYTYP(2,59))
      End if 

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") "The checksums of Hbar(me,bj) at Entry:"
      Write(6,*)
      If (IUhf .Gt. 0) Then
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
         Call Getall(Work, AAAA_LENGTH_MBEJ, IRREPX, 54)
         Call Checksum("AAAA_MBEJ", Work, AAAA_LENGTH_MBEJ)
         Call Getall(Work, AABB_LENGTH_MBEJ, IRREPX, 56)
         Call Checksum("AABB_MBEJ", Work, AABB_LENGTH_MBEJ)
         Call Getall(Work, ABAB_LENGTH_MBEJ, IRREPX, 58)
         Call Checksum("ABAB_MBEJ", Work, ABAB_LENGTH_MBEJ)
      Endif
#endif
C
C Meaning of the C1,C2,C3,C of the Various forms of DWMBEJ.
C C1 > 0.0 Adds the bare integrals.
C C2 > 0.0 actually unused. 
C C3 > 0.0 The t1t1 is added to T2 before contracted with C the integrals.
C C3 > The spin-adapatation of T1ring for RHF is carried out. 
C c4 > Scalled ie C4*T2 is added. The pCC does 2.0*parbeta*c4*T2. 
C
C This block is reserved for when Coulomb is true (Ispar is always
C true since we will not be here otherwise). We need both the CC W(mb,ej)
C and p(D)CC W(mb,ej).
C                   -
C Construct the p(D)CC  W(mb,ej) here in 54-58 lists. Then below construct
C         -
C the CC W(mb,ej) in 264-269 lists. In lambda code they will be
C copied to 54-58 afte the Hbar(abci) and Hbar(iajk) is constructed.
C

C                  -
C Construct the CC W(mb,ej) here (overwrite the W(mb,ej)) 
C
      C1 =  1.0D0
      C2 =  0.0D0
      C3 =  0.0D0
      C4 = -1.0D0
      C5 = -1.0D0

      Do Ilist = 54, 59, 2-Iuhf
         Call Zersym(Work,Ilist)
      Enddo 

      IF (iflags(h_IFLAGS_calclevel) .gt. 9) Then
         AA_FAE_LENGTH = IRPDPD(1,19)
         BB_FAE_LENGTH = 0
         IF (Iuhf .gt. 0) BB_FAE_LENGTH = IRPDPD(1,20)
C
C Saving and restring FAE is not necessary here since pdcc_t1ring 
C does not change FAE (unlike t1ring), but leaving it alone since
C it does no harm.

         I000 = 1
         I010 = I000 + AA_FAE_LENGTH
         I020 = I010 + BB_FAE_LENGTH

         Call Getlst(Work(I000), 1, 1, 1, 1, 92)

         If (Iuhf .Gt. 0) Call Getlst(Work(I010), 1, 1, 1, 2, 92)

#ifdef _DEBUG_LVL0
      Call Checksum("AA_FAE   ", Work(I000), AA_FAE_LENGTH)
      If (Iuhf .gt. 0) Call Checksum("BB_FAE   ", Work(I010),
     &                                BB_FAE_LENGTH)
      Write(6,*)
#endif

         Call T1ring(Work(I020),(Length-I020),Iuhf, .False.)

         Call Putlst(Work(I000), 1, 1, 1, 1, 92)
         If (Iuhf .Gt. 0) Call Putlst(Work(I010), 1, 1, 1, 2, 92)

         If (Iuhf .eq. 0) call Spn_adapt_t1ring(Work,Length,Iuhf,
     &                                          56,58,56)

C This copies the T1ring lists to lists 264-269. The T2 contributions
C are added down below. 

         Call Getall(Work, AABB_LENGTH_MBEJ, IRREPX, 56)
         Call Putall(Work, AABB_LENGTH_MBEJ, IRREPX, 266)
         Call Getall(Work, ABAB_LENGTH_MBEJ, IRREPX, 58)
         Call Putall(Work, ABAB_LENGTH_MBEJ, IRREPX, 268)

         If (Iuhf .Gt. 0) Then
            Call Getall(Work, AAAA_LENGTH_MBEJ, IRREPX, 54)
            Call Putall(Work, AAAA_LENGTH_MBEJ, IRREPX, 264)
            Call Getall(Work, BBBB_LENGTH_MBEJ, IRREPX, 55)
            Call Putall(Work, BBBB_LENGTH_MBEJ, IRREPX, 265)
            Call Getall(Work, AABB_LENGTH_MBEJ, IRREPX, 56)
            Call putall(Work, AABB_LENGTH_MBEJ, IRREPX, 266)
            Call Getall(Work, BBAA_LENGTH_MBEJ, IRREPX, 57)
            Call Putall(Work, BBAA_LENGTH_MBEJ, IRREPX, 267)
            Call Getall(Work, ABAB_LENGTH_MBEJ, IRREPX, 58)
            Call Putall(Work, ABAB_LENGTH_MBEJ, IRREPX, 268)
            Call Getall(Work, BABA_LENGTH_MBEJ, IRREPX, 59)
            Call Putall(Work, BABA_LENGTH_MBEJ, IRREPX, 269)
         Endif

#ifdef _DEBUG_LVL0
      Write(6,"(a)") " In Post_dcc_mods C1, C2, C3 and C4"
      Write(6,"(4(1x,F5.2))") c1, c2, c3, c4
      Write(6,*)
#endif 

      Endif 
      
      Call Dwmbej(Work,Length,"ABAB",Iuhf,C1,C2,C3,C4)
      Call Dwmbej(Work,Length,"ABBA",Iuhf,C1,C2,C3,C4)

      If (Iuhf .Gt. 0) Then
          Call Dwmbej(Work,Length,"AAAA",Iuhf,C1,C2,C3,C4)
          Call Dwmbej(Work,Length,"BBBB",Iuhf,C1,C2,C3,C4)
          Call Dwmbej(Work,Length,"BABA",Iuhf,C1,C2,C3,C4)
          Call Dwmbej(Work,Length,"BAAB",Iuhf,C1,C2,C3,C4)
      Endif 

      Call T1t1_inwmbej(Work,Length,"ABAB",Iuhf,C5)
      Call T1t1_inwmbej(Work,Length,"ABBA",Iuhf,C5)

      If (Iuhf .Gt. 0) Then
          Call T1t1_inwmbej(Work,Length,"AAAA",Iuhf,C5)
          Call T1t1_inwmbej(Work,Length,"BBBB",Iuhf,C5)
          Call T1t1_inwmbej(Work,Length,"BABA",Iuhf,C5)
          Call T1t1_inwmbej(Work,Length,"BAAB",Iuhf,C5)
      Endif

      If (Iuhf .eq. 0) Call Reset(Work, Length, Iuhf)

#ifdef _DEBUG_LVL0
      Write(6,"(a)")   "                  -                "
      Write(6,"(a,a)") " The checksums of W reconstructed p(D)cc."
      Write(6,*)
      IF (Iuhf .Gt. 0) Then
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
         Call Getall(Work, AAAA_LENGTH_MBEJ, IRREPX, 54)
         Call Checksum("AAAA_MBEJ", Work, AAAA_LENGTH_MBEJ)
         Call Getall(Work, AABB_LENGTH_MBEJ, IRREPX, 56)
         Call Checksum("AABB_MBEJ", Work, AABB_LENGTH_MBEJ)
         Call Getall(Work, ABAB_LENGTH_MBEJ, IRREPX, 58)
         Call Checksum("ABAB_MBEJ", Work, ABAB_LENGTH_MBEJ)
      Endif
#endif
C                 -
C Save the p(D)CC W(mb,ej) lists in lists 54-59 to 164-169
C
      Call Getall(Work, AAAA_LENGTH_MBEJ, IRREPX, 54)
      Call Putall(Work, AAAA_LENGTH_MBEJ, IRREPX, 164)
      Call Getall(Work, AABB_LENGTH_MBEJ, IRREPX, 56)
      Call Putall(Work, AABB_LENGTH_MBEJ, IRREPX, 166)
      Call Getall(Work, ABAB_LENGTH_MBEJ, IRREPX, 58)
      Call Putall(Work, ABAB_LENGTH_MBEJ, IRREPX, 168)

      If (Iuhf .Gt. 0) Then
         Call Getall(Work, AAAA_LENGTH_MBEJ, IRREPX, 54)
         Call Putall(Work, AAAA_LENGTH_MBEJ, IRREPX, 164)
         Call Getall(Work, BBBB_LENGTH_MBEJ, IRREPX, 55)
         Call Putall(Work, BBBB_LENGTH_MBEJ, IRREPX, 165)
         Call Getall(Work, AABB_LENGTH_MBEJ, IRREPX, 56)
         Call putall(Work, AABB_LENGTH_MBEJ, IRREPX, 166)
         Call Getall(Work, BBAA_LENGTH_MBEJ, IRREPX, 57)
         Call Putall(Work, BBAA_LENGTH_MBEJ, IRREPX, 167)
         Call Getall(Work, ABAB_LENGTH_MBEJ, IRREPX, 58)
         Call Putall(Work, ABAB_LENGTH_MBEJ, IRREPX, 168)
         Call Getall(Work, BABA_LENGTH_MBEJ, IRREPX, 59)
         Call Putall(Work, BABA_LENGTH_MBEJ, IRREPX, 169)
      Endif

CSSS#ifdef _NO_SKIP
C                  -
C Construct the CC W(mb,ej) in 264-269 lists.
C
      C1 =  1.0D0
      C2 =  0.0D0
      C3 =  0.0D0
      C4 = -1.0D0
      C5 = -1.0D0

#ifdef _NO_NEED 
      Do Ilist = 54, 59, 2-Iuhf
         Call Zersym(Work,Ilist)
      Enddo 

      Do Ilist = 264, 269, 2-Iuhf
         Call Zersym(Work,Ilist)
      Enddo 

      IF (iflags(h_IFLAGS_calclevel) .gt. 9) Then
         AA_FAE_LENGTH = IRPDPD(1,19)
         BB_FAE_LENGTH = 0
         IF (Iuhf .gt. 0) BB_FAE_LENGTH = IRPDPD(1,20)
C
         I000 = 1
         I010 = I000 + AA_FAE_LENGTH
         I020 = I010 + BB_FAE_LENGTH

         Call Getlst(Work(I000), 1, 1, 1, 1, 92)

         If (Iuhf .Gt. 0) Call Getlst(Work(I010), 1, 1, 1, 2, 92)

#ifdef _DEBUG_LVL0
      Call Checksum("AA_FAE   ", Work(I000), AA_FAE_LENGTH)
      If (Iuhf .gt. 0) Call Checksum("BB_FAE   ", Work(I010),
     &                                BB_FAE_LENGTH)
      Write(6,*)
#endif
         Call T1ring(Work(I020),(Length-I020),Iuhf, .False.)

         Call Putlst(Work(I000), 1, 1, 1, 1, 92)
         If (Iuhf .Gt. 0) Call Putlst(Work(I010), 1, 1, 1, 2, 92)

          Call Getall(Work, AABB_LENGTH_MBEJ, IRREPX, 56)
          Call Putall(Work, AABB_LENGTH_MBEJ, IRREPX, 266)
          Call Getall(Work, ABAB_LENGTH_MBEJ, IRREPX, 58)
          Call Putall(Work, ABAB_LENGTH_MBEJ, IRREPX, 268)

          If (Iuhf .eq. 0) call Spn_adapt_t1ring(Work,Length,Iuhf,
     &                                          266,268,266)
          If (Iuhf .Gt. 0) Then
             Call Getall(Work, AAAA_LENGTH_MBEJ, IRREPX, 54)
             Call Putall(Work, AAAA_LENGTH_MBEJ, IRREPX, 264)
             Call Getall(Work, BBBB_LENGTH_MBEJ, IRREPX, 55)
             Call Putall(Work, BBBB_LENGTH_MBEJ, IRREPX, 265)
             Call Getall(Work, AABB_LENGTH_MBEJ, IRREPX, 56)
             Call putall(Work, AABB_LENGTH_MBEJ, IRREPX, 266)
             Call Getall(Work, BBAA_LENGTH_MBEJ, IRREPX, 57)
             Call Putall(Work, BBAA_LENGTH_MBEJ, IRREPX, 267)
             Call Getall(Work, ABAB_LENGTH_MBEJ, IRREPX, 58)
             Call Putall(Work, ABAB_LENGTH_MBEJ, IRREPX, 268)
             Call Getall(Work, BABA_LENGTH_MBEJ, IRREPX, 59)
             Call Putall(Work, BABA_LENGTH_MBEJ, IRREPX, 269)
          Endif
#ifdef _DEBUG_LVL0
      Write(6,"(a)") " In Post_dcc_mods C1, C2, C3 and C4"
      Write(6,"(4(1x,F5.2))") c1, c2, c3, c4
      Write(6,*)
#endif 
      Endif 
#endif
C Above endif is for the _NO_NEED block. 
 
      Call Pdcc_dwmbej(Work,Length,"ABAB",Iuhf,C1,C2,C3,C4,266,267,
     &                 .False.,.False.,1.0D0)
      Call Pdcc_dwmbej(Work,Length,"ABBA",Iuhf,C1,C2,C3,C4,268,268,
     &                .False.,.False.,1.0D0)

      If (Iuhf .Gt. 0) Then
          Call Pdcc_dwmbej(Work,Length,"AAAA",Iuhf,C1,C2,C3,C4,263,263,
     &                    .False.,.False.,1.0D0)
          Call Pdcc_dwmbej(Work,Length,"BBBB",Iuhf,C1,C2,C3,C4,263,263,
     &                    .False.,.False.,1.0D0)
          Call Pdcc_dwmbej(Work,Length,"BABA",Iuhf,C1,C2,C3,C4,267,266,
     &                    .False.,.False.,1.0D0)
          Call Pdcc_dwmbej(Work,Length,"BAAB",Iuhf,C1,C2,C3,C4,269,269,
     &                    .False.,.False.,1.0D0)
      Endif 

      Call Pdcc_t1t1_inwmbej(Work,Length,"ABAB",Iuhf,C5,266,267)
      Call Pdcc_t1t1_inwmbej(Work,Length,"ABBA",Iuhf,C5,268,268)

      If (Iuhf .Gt. 0) Then
          Call Pdcc_t1t1_inwmbej(Work,Length,"AAAA",Iuhf,C5,263,263)
          Call Pdcc_t1t1_inwmbej(Work,Length,"BBBB",Iuhf,C5,263,263)
          Call Pdcc_t1t1_inwmbej(Work,Length,"BABA",Iuhf,C5,267,266)
          Call Pdcc_t1t1_inwmbej(Work,Length,"BAAB",Iuhf,C5,269,269)
      Endif

      If (Iuhf .eq. 0) Call Pdcc_Reset(Work,Length,Iuhf,266,268,264)

#ifdef _DEBUG_LVL0
      Write(6,"(a)")   "                  -                "
      Write(6,"(a,a)") " The checksums of W reconstructed as (CC)."
      Write(6,*)
      IF (Iuhf .Gt. 0) Then
         Call Getall(Work, AAAA_LENGTH_MBEJ, IRREPX, 264)
         Call Checksum("AAAA_MBEJ", Work, AAAA_LENGTH_MBEJ)
         Call Getall(Work, BBBB_LENGTH_MBEJ, IRREPX, 265)
         Call Checksum("BBBB_MBEJ", Work, BBBB_LENGTH_MBEJ)
         Call Getall(Work, AABB_LENGTH_MBEJ, IRREPX, 266)
         Call Checksum("AABB_MBEJ", Work, AABB_LENGTH_MBEJ)
         Call Getall(Work, BBAA_LENGTH_MBEJ, IRREPX, 267)
         Call Checksum("BBAA_MBEJ", Work, BBAA_LENGTH_MBEJ)
         Call Getall(Work, ABAB_LENGTH_MBEJ, IRREPX, 268)
         Call Checksum("ABAB_MBEJ", Work, AABB_LENGTH_MBEJ)
         Call Getall(Work, BABA_LENGTH_MBEJ, IRREPX, 269)
         Call Checksum("BABA_MBEJ", Work, BABA_LENGTH_MBEJ)
      Else
         Call Getall(Work, AAAA_LENGTH_MBEJ, IRREPX, 264)
         Call Checksum("AAAA_MBEJ", Work, AAAA_LENGTH_MBEJ)
         Call Getall(Work, AABB_LENGTH_MBEJ, IRREPX, 266)
         Call Checksum("AABB_MBEJ", Work, AABB_LENGTH_MBEJ)
         Call Getall(Work, ABAB_LENGTH_MBEJ, IRREPX, 268)
         Call Checksum("ABAB_MBEJ", Work, ABAB_LENGTH_MBEJ)
      Endif
#endif
CSSS#endif 
C
C Compute a I(mb,ej) where I(mb,ej) is t(mn,ef)<nj||fb> and store
C them in 254-259 lists. In lambda code I would neeed -t(mn,ef)<nj||fb>. 
C
      C1 =  0.0D0
      C2 =  0.0D0
      C3 =  0.0D0
      C4 =  1.0D0

      Call Pdcc_dwmbej(Work,Length,"ABAB",Iuhf,C1,C2,C3,C4,256,257,
     &                 .False.,.False.,1.0D0)
      Call Pdcc_dwmbej(Work,Length,"ABBA",Iuhf,C1,C2,C3,C4,258,258,
     &                .False.,.False.,1.0D0)

      If (Iuhf .Gt. 0) Then
          Call Pdcc_dwmbej(Work,Length,"AAAA",Iuhf,C1,C2,C3,C4,253,253,
     &                    .False.,.False.,1.0D0)
          Call Pdcc_dwmbej(Work,Length,"BBBB",Iuhf,C1,C2,C3,C4,253,253,
     &                    .False.,.False.,1.0D0)
          Call Pdcc_dwmbej(Work,Length,"BABA",Iuhf,C1,C2,C3,C4,257,256,
     &                    .False.,.False.,1.0D0)
          Call Pdcc_dwmbej(Work,Length,"BAAB",Iuhf,C1,C2,C3,C4,259,259,
     &                    .False.,.False.,1.0D0)  
      Endif 

      If (Iuhf .eq. 0) Call Pdcc_Reset(Work,Length,Iuhf,256,258,254)
C
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a,a)") " The checksums of I(mb,ej)=t(mn,ef)<nj||fb> ",
     &                 "at exit:"
      Write(6,*)
      IF (Iuhf .Gt. 0) Then
         AAAA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,54),ISYTYP(2,254))
         BBBB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,55),ISYTYP(2,255))
         AABB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,56),ISYTYP(2,256))
         BBAA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,57),ISYTYP(2,257))
         ABAB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,58),ISYTYP(2,258))
         BABA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,59),ISYTYP(2,259))
         Call Getall(Work, AAAA_LENGTH_MBEJ, IRREPX, 254)
         Call Checksum("AAAA_MBEJ", Work, AAAA_LENGTH_MBEJ)
         Call Getall(Work, BBBB_LENGTH_MBEJ, IRREPX, 255)
         Call Checksum("BBBB_MBEJ", Work, BBBB_LENGTH_MBEJ)
         Call Getall(Work, AABB_LENGTH_MBEJ, IRREPX, 256)
         Call Checksum("AABB_MBEJ", Work, AABB_LENGTH_MBEJ)
         Call Getall(Work, BBAA_LENGTH_MBEJ, IRREPX, 257)
         Call Checksum("BBAA_MBEJ", Work, BBAA_LENGTH_MBEJ)
         Call Getall(Work, ABAB_LENGTH_MBEJ, IRREPX, 258)
         Call Checksum("ABAB_MBEJ", Work, AABB_LENGTH_MBEJ)
         Call Getall(Work, BABA_LENGTH_MBEJ, IRREPX, 259)
         Call Checksum("BABA_MBEJ", Work, BABA_LENGTH_MBEJ)
      Else
         AAAA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,254),ISYTYP(2,254))
         AABB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,256),ISYTYP(2,256))
         ABAB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,258),ISYTYP(2,258))
         Call Getall(Work, AAAA_LENGTH_MBEJ, IRREPX, 254)
         Call Checksum("AAAA_MBEJ", Work, AAAA_LENGTH_MBEJ)
         Call Getall(Work, AABB_LENGTH_MBEJ, IRREPX, 256)
         Call Checksum("AABB_MBEJ", Work, AABB_LENGTH_MBEJ)
         Call Getall(Work, ABAB_LENGTH_MBEJ, IRREPX, 258)
         Call Checksum("ABAB_MBEJ", Work, AABB_LENGTH_MBEJ)
      Endif
#endif
C
C Modify the FAE and FMI intermediates. 
C
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") " Checksums of FAE and FMI : At Entry"
      Write(6,*)
      IHHA_LENGTH = IRPDPD(IRREPX,21)
      IHHB_LENGTH = IRPDPD(IRREPX,22)
      CALL Getlst(Work,1,1,1,1,91)
      Call Checksum("FMI", Work, IHHA_LENGTH)
      IF (IUHF.gt.0) Call Getlst(Work,1,1,1,2,91)
      IF (IUHF.gt.0) Call Checksum("Fmi", Work, IHHB_LENGTH)
      IPPA_LENGTH = IRPDPD(IRREPX,19)
      IPPB_LENGTH = IRPDPD(IRREPX,20)
      CALL Getlst(Work,1,1,1,1,92)
      Call Checksum("FEA", Work, IPPA_LENGTH)
      IF (IUHF.gt.0) Call Getlst(Work,1,1,1,2,92)
      IF (IUHF.gt.0) Call Checksum("Fea", Work, IPPB_LENGTH)
#endif 

      Call zerlst(Work,IHHA_LENGTH,1,1,9,91)
      Call zerlst(Work,IPPA_LENGTH,1,1,9,92)
      If (iuhf .ne. 0) then
            call zerlst(Work,IHHB_LENGTH,1,1,10,91)
            call zerlst(Work,IPPB_LENGTH,1,1,10,92)
      Endif 

      Call pdcc_quad2(Work,Length,iuhf,Fmi_scale,8,0)
      Call pdcc_quad3(Work,Length,iuhf,Fae_scale,8,0)

#ifdef _DEBUG_LVLM
      Write(6,*)
      Write(6,"(a)") " The checksums of all the intms at VCC exit"
      Write(6,*)
      Call checkintms(Work,Length,Iuhf,2)
      Write(6,*)
#endif 
      Call Zersym(Work,119)
      Call Zersym(Work,120)

      Return
      End
 

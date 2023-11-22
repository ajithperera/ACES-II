      SUBROUTINE GENINT(ICORE,MAXCOR,IUHF,INTPCK,DO_HBAR_4LCCSD)
C
C DRIVER FOR INTERMEDIATE FORMATION.  
C
C   INTPCK = 0  , DUMP PACKED INTEGRALS ONLY FOR W INTERMEDIATES (NO F).
C   INTPCK = 1  , DUMP ONLY QUADRATIC PART OF W INTERMEDIATES AND F.
C   INTPCK = 2  , DUMP FULL W AND F INTERMEDIATES (For standard CCSD or beyond)
C   INTPCK = 3  , DUMP FIRST ORDER W AND F ONLY.  (For Linearized CC methods)
C
CEND
C
#include "flags.h" /* for h_IFLAGS_* and h_IFLAGS2_* */
C
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION TCPU,TSYS
      DOUBLE PRECISION C1,C2,C3,C4,FACT,C5
      LOGICAL DO_HBAR_4LCCSD
C
      DIMENSION ICORE(MAXCOR)
      LOGICAL MBPT3,MBPT4,CC,TRPEND,SNGEND,GRAD,MBPTT,SING1,
     &                QCISD,LAMBDA,NONHF,ROHF4,ITRFLG,UCC,
     &        NO_REDUNDANT,BRUECKNER
      EQUIVALENCE (METHOD,IFLAGS(2))
C
      COMMON/SWITCH/ MBPT3,MBPT4,CC,TRPEND,SNGEND,GRAD,MBPTT,SING1,
     &                QCISD,UCC
      COMMON/NHFREF/NONHF
      COMMON/ROHF/ROHF4,ITRFLG
      COMMON/FLAGS/IFLAGS(100)
      COMMON/FLAGS2/IFLAGS2(500)
      COMMON/SYM/POP(8,2),VRT(8,2),NT1AA,NT1BB,NF1AA,NF1BB,
     &           NF2AA,NF2BB
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)

#include "parcc.com"

      EQUIVALENCE(METHOD,IFLAGS(2))
C
C W INTERMEDIATES. SET COEFFICIENTS BASED ON CALCULATION TYPE.
C
      IF(IFLAGS(1).GE.10)THEN
         WRITE(6,100)
 100     FORMAT(T3,'@GENINT-I, Calculating intermediates.')
         CALL TIMER(1)
      ENDIF
C
      NO_REDUNDANT = iflags2(h_IFLAGS2_single_store).ne.0
      BRUECKNER    = iflags(h_IFLAGS_brueckner) .ne.0
C
      IF(INTPCK.EQ.2)THEN
         C1=1.0
         C2=0.0
         C3=0.0
         C4=-0.5
      ELSEIF(INTPCK.EQ.1)THEN
         C1=0.0
         C2=0.0
         C3=0.0
         C4=-0.5 
      ELSEIF(INTPCK.EQ.0.OR.INTPCK.EQ.3)THEN
         C1=1.0
         C2=0.0
         C3=0.0
         C4=0.0
      ENDIF
C
      FACT=1.D0
      IF(UCC) THEN
         FACT=0.5
         C4=C4*FACT
      ENDIF
C
      IF(METHOD.GT.9.AND.SING1.AND..NOT.QCISD.AND.INTPCK.NE.3) C3=1.0D0
C
C ---------------The first part of the W(mn,ij) begins Here -------------
C
C The <mn||ij> + 1/4 Sum Tau(ij,ef)<mn|ef> contribution to the W(mn,ij)
C intermediate. The other piece is calculated in T1INW2 (see below).
C
CSSS#ifdef _DCC_FLAG
C
C The DCC W(mn,ij) = <mn||ij> + 1/4 Sum T1*T1*<mn|ef>
C
      If (Ispar) write(6,"(a,a)") " ***The pCC/DCC W(mn,ij) is",
     &                            " constructed***"
CSSS#endif 
C
      CALL QUAD1(ICORE,MAXCOR,INTPCK,IUHF,FACT)

#ifdef _DEBUG_LVLM
      write(6,*) "After Q1"
      call checkintms(icore,Maxcor,iuhf,0)
#endif
CSSS#ifdef _DCC_FLAG
      If (Ispar) Then
C
C Note that the above W(mn,ij) leave out the term T1*T1*T1*T1*W(mn,ij)
C Construct it and add it to the T2 equation from here. 
C
         write(6,"(a)") " *** T1*T1*W(mn,ij) to T2 is constructed***"
         CALL CONSTRUCT_T2WMNIJ(ICORE,MaXCOR,IUHF)
       Endif 
CSS#endif 

C
C ---------------The fist part of the W(mn,ij) ends Here -----------------
C
C--------The W(mbej) and parts of F(ae) intermediate begins Here ---------
C
      DO 333 ILIST = 54,59,2-IUHF
         CALL ZERSYM(ICORE,ILIST) 
 333  CONTINUE
C
C The T1RING is responsible for calculating Sum T(jf)*<mb||ef> - 
C Sum T(nb)<mn||ej> contributions to W(mbej) intermediates. Also, calculates 
C the Sum T(mf)<ma||fe>  contributions to the F(ae) intermediate.
C All three contributinons must be included in the respective
C intermediates for LCCSD excitaion energy calculations. So, call genint
C with DO_HBAR_4LCCSD to true to acomplishe this since when the
C calc=lccsd, this does not goto this block (even when INTPCK=2 and
C SING=.TRUE. 05/2015, Ajith Perera.
C
      INCRF=0
      IF ((METHOD.GT.9.AND.SING1.AND.(.NOT.QCISD)).OR.
     &   (ROHF4.AND.INTPCK.NE.0) .OR. DO_HBAR_4LCCSD) THEN
C
C Non-zero T1 contributions to W intermediates for RHF/UHF ref., only if 
C CALCLEVEL is beyond CCSD, but not QCISD and SING1 which is true if there is
C non zero T1's available from a prior iteration or ROHF-MBPT(4) and 
C 
         CALL ZERLST(ICORE,NF2AA,1,1,1,92)
         IF (IUHF.NE.0) CALL ZERLST(ICORE,NF2BB,1,1,2,92)
         LAMBDA=.FALSE.
         CALL T1RING(ICORE,MAXCOR,IUHF,LAMBDA)
         INCRF = 1
      ENDIF
CSSS#ifdef _DCC_FLAG
      If (Ispar) Then
         Write(6,"(a)") " ***The pCC/DCC W(mb,ej) is constructed***" 
         Write(6,*)
         C1 =  1.0D0
         C2 =  0.0D0
         C3 =  0.0D0
         C4 = -0.50D0
         C5 = -1.0D0

#ifdef _DEBUG_LVL0
      Write(6,"(a)") " In GENINT C1, C2, C3 and C4"
      Write(6,"(4(1x,F5.2))") c1, c2, c3, c4
      Write(6,*)
#endif 
         If (Iflags(h_IFLAGS_calc) .gt. 9) Then

            If (Iuhf .eq. 0) call Spn_adapt_t1ring(Icore,Maxcor,Iuhf,
     &                                             56,58,56)  
         Endif

#ifdef _DEBUG_LVL0
        call scale_h4(icore,maxcor,iuhf,.false.,.true.)
#endif 
#ifdef _DEBUG_LVL0
        call scale_h4(icore,maxcor,iuhf,.true.,.false.)
#endif 
         Call Dwmbej(Icore,Maxcor,"ABAB",Iuhf,C1,C2,C3,C4)
#ifdef _DEBUG_LVL0
        call scale_h4(icore,maxcor,iuhf,.false.,.false.)
#endif 
         Call Dwmbej(Icore,Maxcor,"ABBA",Iuhf,C1,C2,C3,C4)

         If (Iuhf .Gt. 0) Then
             Call Dwmbej(Icore,Maxcor,"AAAA",Iuhf,C1,C2,C3,C4)
             Call Dwmbej(Icore,Maxcor,"BBBB",Iuhf,C1,C2,C3,C4)
             Call Dwmbej(Icore,Maxcor,"BABA",Iuhf,C1,C2,C3,C4)
             Call Dwmbej(Icore,Maxcor,"BAAB",Iuhf,C1,C2,C3,C4)
         Endif

         If (Iflags(h_IFLAGS_calc) .gt. 9) Then 
            Call T1t1_inwmbej(Icore,Maxcor,"ABAB",Iuhf,C5)
            Call T1t1_inwmbej(Icore,Maxcor,"ABBA",Iuhf,C5)

            If (Iuhf .Gt. 0) Then
               Call T1t1_inwmbej(Icore,Maxcor,"AAAA",Iuhf,C5)
               Call T1t1_inwmbej(Icore,Maxcor,"BBBB",Iuhf,C5)
               Call T1t1_inwmbej(Icore,Maxcor,"BABA",Iuhf,C5)
               Call T1t1_inwmbej(Icore,Maxcor,"BAAB",Iuhf,C5)
            Endif

         Endif 
CSSS#else
      Else
CSSS      Write(6,*) "Traditional pcc path"
C  
C The following calls calculate <mb||ej>-(1/2)*SUM [T(jn,fb)<mn||ef>] 
C contribution to W(mbej) intermediate. 
C
         IF (.NOT. NO_REDUNDANT) THEN

C#ifndef _T2T2L2_DEBUG
C#define _T2T2L2_DEBUG
C#endif 
C
#ifdef _T2T2L2_DEBUG
        Write(6,*)
        Write(6,"(a)") "!!!!!!DEBUGGING G(IJ,AB) FROM T2*T2*L2!!!!!!"
        If (Iuhf .eq. 0) call Spn_adapt_t1ring(Icore,Maxcor,Iuhf,
     &                                             56,58,56)  
         C1 =  1.0D0
         C2 =  0.0D0
         C3 =  0.0D0
         C4 = -0.50D0
         C5 = -1.0D0

#ifdef _DEBUG_LVL0
      Write(6,"(a)") " In GENINT C1, C2, C3 and C4"
      Write(6,"(4(1x,F5.2))") c1, c2, c3, c4
      call scale_h4(icore,maxcor,iuhf,.false.,.true.)
#endif
         Call Dwmbej(Icore,Maxcor,"ABAB",Iuhf,C1,C2,C3,C4)

#ifdef _DEBUG_LVL0
      call scale_h4(icore,maxcor,iuhf,.true.,.false.)
#endif
         Call Dwmbej(Icore,Maxcor,"ABBA",Iuhf,C1,C2,C3,C4)

#ifdef _DEBUG_LVL0
        call scale_h4(icore,maxcor,iuhf,.false.,.false.)
#endif

         If (Iuhf .Gt. 0) Then
             Call Dwmbej(Icore,Maxcor,"AAAA",Iuhf,C1,C2,C3,C4)
             Call Dwmbej(Icore,Maxcor,"BBBB",Iuhf,C1,C2,C3,C4)
             Call Dwmbej(Icore,Maxcor,"BABA",Iuhf,C1,C2,C3,C4)
             Call Dwmbej(Icore,Maxcor,"BAAB",Iuhf,C1,C2,C3,C4)
         Endif

         If (Iflags(h_IFLAGS_calc) .gt. 9) Then
            Call T1t1_inwmbej(Icore,Maxcor,"ABAB",Iuhf,C5)
            Call T1t1_inwmbej(Icore,Maxcor,"ABBA",Iuhf,C5)

            If (Iuhf .Gt. 0) Then
               Call T1t1_inwmbej(Icore,Maxcor,"AAAA",Iuhf,C5)
               Call T1t1_inwmbej(Icore,Maxcor,"BBBB",Iuhf,C5)
               Call T1t1_inwmbej(Icore,Maxcor,"BABA",Iuhf,C5)
               Call T1t1_inwmbej(Icore,Maxcor,"BAAB",Iuhf,C5)
            Endif

         Endif
#else 
            IF(IUHF.EQ.1)THEN
              CALL DWMBEJ(ICORE,MAXCOR,'AAAA',IUHF,C1,C2,C3,C4)
              CALL DWMBEJ(ICORE,MAXCOR,'BBBB',IUHF,C1,C2,C3,C4)
              CALL DWMBEJ(ICORE,MAXCOR,'BABA',IUHF,C1,C2,C3,C4)
              CALL DWMBEJ(ICORE,MAXCOR,'BAAB',IUHF,C1,C2,C3,C4)
            ENDIF
C
            CALL DWMBEJ(ICORE,MAXCOR,'ABAB',IUHF,C1,C2,C3,C4)
            CALL DWMBEJ(ICORE,MAXCOR,'ABBA',IUHF,C1,C2,C3,C4)
#endif 
C
CSSS         CALL CHECKLIST(ICORE, MAXCOR, 54)
CSSS         IF (IUHF .NE. 0) CALL CHECKLIST(ICORE, MAXCOR, 55)
CSSS         CALL CHECKLIST(ICORE, MAXCOR, 56)
CSSS         IF (IUHF. NE. 0) CALL CHECKLIST(ICORE, MAXCOR, 57)
CSSS         CALL CHECKLIST(ICORE, MAXCOR, 58)
CSSS         IF (IUHF. NE. 0) CALL CHECKLIST(ICORE, MAXCOR, 59)
C
         ELSE
            IF(IUHF.EQ.1)THEN
              CALL DWMBEJNR(ICORE,MAXCOR,'AAAA',IUHF,C1,C2,C3,C4)
              CALL DWMBEJNR(ICORE,MAXCOR,'BBBB',IUHF,C1,C2,C3,C4)
              CALL DWMBEJNR(ICORE,MAXCOR,'BABA',IUHF,C1,C2,C3,C4)
              CALL DWMBEJNR(ICORE,MAXCOR,'BAAB',IUHF,C1,C2,C3,C4)
         ENDIF

         CALL DWMBEJNR(ICORE,MAXCOR,'ABAB',IUHF,C1,C2,C3,C4)
         CALL DWMBEJNR(ICORE,MAXCOR,'ABBA',IUHF,C1,C2,C3,C4)
C
CSSS         CALL CHECKLIST(ICORE, MAXCOR, 54)
CSSS         IF (IUHF .NE. 0) CALL CHECKLIST(ICORE, MAXCOR, 55)
CSSS         CALL CHECKLIST(ICORE, MAXCOR, 56)
CSSS         IF (IUHF .NE. 0)  CALL CHECKLIST(ICORE, MAXCOR, 57)
CSSS         CALL CHECKLIST(ICORE, MAXCOR, 58)
CSSS         IF (IUHF. NE. 0)  CALL CHECKLIST(ICORE, MAXCOR, 59)

         ENDIF

      ENDIF

CSSS#endif 
#ifdef _DEBUG_LVLM
      Irrepx = 1
      Write(6,*)

      IF (iUHF .gt. 0) Then
         iAAAA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,54),ISYTYP(2,54))
         iBBBB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,55),ISYTYP(2,55))
         iAABB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,56),ISYTYP(2,56))
         iBBAA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,57),ISYTYP(2,57))
         iABAB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,58),ISYTYP(2,58))
         iBABA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,59),ISYTYP(2,59))
         Call Getall(Icore, iAAAA_LENGTH_MBEJ, IRREPX, 54)
         Call Checksum("AAAA_MBEJ", icore, iAAAA_LENGTH_MBEJ)
         Call Getall(icore, iBBBB_LENGTH_MBEJ, IRREPX, 55)
         Call Checksum("BBBB_MBEJ", icore, iBBBB_LENGTH_MBEJ)
         Call Getall(icore, iAABB_LENGTH_MBEJ, IRREPX, 56)
         Call Checksum("AABB_MBEJ", icore, iAABB_LENGTH_MBEJ)
         Call Getall(icore, iBBAA_LENGTH_MBEJ, IRREPX, 57)
         Call Checksum("BBAA_MBEJ", icore, iBBAA_LENGTH_MBEJ)
         Call Getall(icore, iABAB_LENGTH_MBEJ, IRREPX, 58)
         Call Checksum("ABAB_MBEJ", icore, iAABB_LENGTH_MBEJ)
         Call Getall(icore, iBABA_LENGTH_MBEJ, IRREPX, 59)
         Call Checksum("BABA_MBEJ", icore, iBABA_LENGTH_MBEJ)
      Else
         iAAAA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,54),ISYTYP(2,54))
         iAABB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,56),ISYTYP(2,56))
         iABAB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,58),ISYTYP(2,58))
         Call Getall(icore, iAAAA_LENGTH_MBEJ, IRREPX, 54)
         Call Checksum("AAAA_MBEJ", icore, iAAAA_LENGTH_MBEJ)
         Call Getall(icore, iAABB_LENGTH_MBEJ, IRREPX, 56)
         Call Checksum("AABB_MBEJ", icore, iAABB_LENGTH_MBEJ)
         Call Getall(icore, iABAB_LENGTH_MBEJ, IRREPX, 58)
         Call Checksum("ABAB_MBEJ", icore, iAABB_LENGTH_MBEJ)
      Endif
#endif

      IF(IFLAGS(1).GE.10)THEN
       CALL TIMER(1)
c       WRITE(6,101)TIMENEW
101    FORMAT(T3,'@GENINT-I, W intermediates required ',F9.3,
     &           ' seconds.')
      ENDIF
C
C ----------------------W(mb,ej) and parts of F(ae) ends here-----------
C
      IF (ROHF4.AND.INTPCK.EQ.0) RETURN
C
C  Linearized CC or ROHF-MBPT(4)
C
      IF(INTPCK.EQ.3.OR.ROHF4) THEN
C
         IF(NONHF)THEN
C
            CALL ZERLST(ICORE,NF1AA,1,1,1,91)
            CALL ZERLST(ICORE,NT1AA,1,1,1,93)
            IF(IUHF.NE.0)THEN
               CALL ZERLST(ICORE,NF1BB,1,1,2,91)
               CALL ZERLST(ICORE,NT1BB,1,1,2,93)
            ENDIF
C 
            IF(ROHF4) THEN
               CALL NHFFINT(ICORE,MAXCOR,IUHF,.TRUE.)
            ELSE
               CALL NHFFINT(ICORE,MAXCOR,IUHF,.FALSE.)
            ENDIF
C
            IF(.NOT.SING1)THEN
               CALL ZERLST(ICORE,NT1AA,1,1,1,90)
               CALL ZERLST(ICORE,NT1BB,1,1,2,90)
            ENDIF
            IF(.NOT.ROHF4)RETURN
C

         ELSEIF(SING1)THEN
C
            CALL ZERLST(ICORE,NF1AA,1,1,1,91)
            CALL ZERLST(ICORE,NF2AA,1,1,1,92)
            CALL ZERLST(ICORE,NT1AA,1,1,1,93)
C
            IF(IUHF.NE.0)THEN
               CALL ZERLST(ICORE,NF1BB,1,1,2,91)
               CALL ZERLST(ICORE,NF2BB,1,1,2,92)
               CALL ZERLST(ICORE,NT1BB,1,1,2,93)
            ENDIF
         ENDIF
C
         IF(.NOT.ROHF4)RETURN
C
      ENDIF
C
C ---- Pieces of F(AE), F(MI) AND F(ME) intermediates begin Here------------
C ---- Note that part of the F(AE) is already calculated in T1RING----------
C
      IF(INTPCK.EQ.0) RETURN  ! No need for any of the F intermediates.

      IF(IFLAGS(1).GE.10)CALL TIMER(1)
C
      IF(ROHF4) THEN
        NINCRF=INCRF
      ELSE
        NINCRF=0
      ENDIF
C
C Calculate 1/2 Sum Tau(in,ef)<mn|ef> contribution to the F(mi) intermediate
C
      CALL QUAD2(ICORE,MAXCOR,IUHF,NINCRF,FACT)
C
C Calculate 1/2 Sum Tau(mn,af)<mn|ef> contribution to the F(ae) intermediate
C The other piece (Sum T(mf) <ma||fe>) is calculated in T1RING (see above).

#ifdef _DEBUG_LVLM
      write(6,*) "After Q2"
      call checkintms(icore,Maxcor,iuhf,0)
#endif 
      CALL QUAD3(ICORE,MAXCOR,IUHF,INCRF,FACT)

#ifdef _DEBUG_LVLM
      call checkintms(icore,Maxcor,iuhf,0)
#endif 
C
C Need to call for methods other than CCSD when Brueckner reference 
C is used. Added for EOM-BMBPT(2). 07/2015, Ajith Perera 

      IF((METHOD.GE.9.AND.SING1).OR.ROHF4 .OR. DO_HBAR_4LCCSD .OR.
     &    BRUECKNER) THEN

C 
C Calculate F(me) + Sum T(nf)<mn||ef> contribution to the F(me). That 
C is all for the total F(me) intermediate. 
C
       CALL MAKFME(ICORE,MAXCOR,IUHF,1)
       IF(IUHF.NE.0) CALL MAKFME(ICORE,MAXCOR,IUHF,2)
C 
      ENDIF
C
      IF((NONHF.AND..NOT.ROHF4) .OR. (NONHF.AND.DO_HBAR_4LCCSD)) THEN
C
C Calculate (1-delta) f(ae) - 1/2 Sum f(me) t(am) and (1-delta)f(mi) - 1/2 Sum
C T(ie) f(me). These terms have non zero contributions only for NON-HF methods.
C ROHF-MBPT(4) does not need these contributions. 
C
        CALL NHFFINT(ICORE,MAXCOR,IUHF,SING1)
C
      ENDIF
C
      IF(IFLAGS(1).GE.10)THEN
       CALL TIMER(1)
c       WRITE(6,102)TIMENEW
102    FORMAT(T3,'@GENINT-I, F intermediates required ',F9.3,
     &           ' seconds.')
      ENDIF
C
C------The F(AE) AND F(ME) and parts of F(MI) intermediates end Here -------
C
C----- Rest of the W(mnij) and F(MI) and part of W(abef) contribution to
C-----------------------begins Here -----------------------------------
C
      IF((METHOD.GT.9.AND.SING1.AND.(.NOT.QCISD)).OR.ROHF4) THEN
C
C Calculate Sum T(ej)<mn||ie> contribution to W(mn,ij) intermediate. Also,
C add the Sum T(ne)<mn||ie> contriution to the F(MI) intermediate. At this 
C point F(MI) intermediate is complete. 
C
       CALL T1INW2(ICORE,MAXCOR,IUHF)
C
C Calculate the Sum T(bm) T(ef,ij) <am||ef> contribution to T2. Note that 
C this is not a contributio to the W(ab,ef) intermediate, instead it is a 
C contribution to the T2(ij,ab). This can avoid storing W(ab,ef) on the disk
C at this point. 
C
       CALL T1INW1(ICORE,MAXCOR,IUHF)
C
      ENDIF
C
C This should be invoked only in the case of LCCSD excitation energies
C are requested. This will add the T(ej)<mn||ie> to W(mn,ij) 
C 05/2015, Ajith Perera.
C 
      IF (DO_HBAR_4LCCSD) CALL T1INW2(ICORE,MAXCOR,IUHF)
C
C-------------------------W(mn,ij) and F(MI) is complete----------------
C
      IF(IFLAGS(1).GE.10)THEN
       CALL TIMER(1)
c       WRITE(6,103)TIMENEW
103    FORMAT(T3,'@GENINT-I, additional terms for CCSD required ',F9.3,
     &           ' seconds.')
      ENDIF
C
      RETURN
      END

C
      SUBROUTINE MKDBLAB1(Z, ZT, T ,TT, W, MAXSIZE, TA, TB, POP1, POP2,
     &                    VRT1, VRT2, DISSYZ, DISSYWA, DISSYWB, DISSYT, 
     &                    NUMSYZ, NUMSYWA, NUMSYWB, NUMSYT, NTASIZ, 
     &                    NTBSIZ, LISTT, LISTZ, LISTWA, LISTWB, 
     &                    IRREPWAL, IRREPWAR, IRREPWBL, IRREPWBR, 
     &                    IRREPTL, IRREPTR, IRREPQL, IRREPQR, IRREPX,
     &                    IUHF, TMP)
C
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      INTEGER DISSYZ,DISSYWA,DISSYWB,DIRPRD,POP1,POP2,VRT1,VRT2,DISSYT
      DIMENSION Z(DISSYZ,NUMSYZ),W(DISSYWA,1),TA(NTASIZ)
      DIMENSION T(DISSYT,NUMSYT),TB(NTBSIZ)
      DIMENSION ZT(NUMSYZ,DISSYZ),TT(NUMSYT,DISSYT)
      DIMENSION TMP(1),POP1(8),POP2(8),VRT1(8),VRT2(8),ITOFF1(8),
     &          ITOFF2(8),IZOFF1(8),IZOFF2(8)
C     
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/  LUOUT, MOINTS
      COMMON /FLAGS/ IFLAGS(100)
C   
C Common blocks used in the quadratic term
C 
      COMMON /QADINTI/ INTIMI, INTIAE, INTIME
      COMMON /QADINTG/ INGMNAA, INGMNBB, INGMNAB, INGABAA,
     &                 INGABBB, INGABAB, INGMCAA, INGMCBB,
     &                 INGMCABAB, INGMCBABA, INGMCBAAB,
     &                 INGMCABBA
      COMMON /APERTT1/ IAPRT1AA             
      COMMON /BPERTT1/ IBPRT1AA
      COMMON /APERTT2/ IAPRT2AA1, IAPRT2BB1, IAPRT2AB1, IAPRT2AB2, 
     &                 IAPRT2AB3, IAPRT2AB4, IAPRT2AA2, IAPRT2BB2,
     &                 IAPRT2AB5 
C
      DATA AZERO,ONE,ONEM /0.0D0,1.0D0,-1.0D0/
C
C Zero first output array (Neccesary because this is not necessarily
C done in the matrix multiplication
C
      CALL ZERO(Z, NUMSYZ*DISSYZ)
C     
C If there are no WA integrals there is nothing to do 
C     
      IF (MIN(NUMSYWA, DISSYWA) .NE. 0) THEN
C     
C Get T2 amplitudes (for quadratic term there is no neeed to form 
C the TAU amplitudes.
C     
         CALL GETLST(T, 1, NUMSYT, 1, IRREPTR, LISTT)
C     
C Process as many distributions as possible
C Get integrals Hbar(E<F,AM) from LISTW
C
         NINCOR = MAXSIZE/(DISSYWA*IINTFP)
         NLEFT  = NUMSYWA
         NFIRST = 1
         NPASS  = 0
 1       NREAD  = MIN(NLEFT, NINCOR)
C
         CALL GETLST(W, NFIRST, NREAD, 1, IRREPWAR, LISTWA)
C
C Perform multiplication 
C          +
C   T(Ef,Ij)*Hbar(Ef,Am) = Z(Ij,Am)
C     
         CALL XGEMM ('T', 'N', NUMSYT, NREAD, DISSYT, ONE, T,
     &                DISSYT, W, DISSYWA, AZERO, ZT(1,NFIRST), NUMSYZ)
C
         NFIRST = NFIRST + NREAD
         NLEFT  = NLEFT - NREAD
         NPASS  = NPASS + 1
         IF (NLEFT .NE. 0)GOTO 1
C     
C Do the second part of the multiplication
C     
         CALL ZERO(TT, NUMSYZ*DISSYZ)
C
         JOFFT = 1
         ITOFF1(1) = 1
         ITOFF2(1) = 1
C         
C Setup an offset for the T2 array
C
         DO 5000 IRREPM = 2, NIRREP
C
            IRREP = IRREPM - 1
            IRREPAO = DIRPRD(IRREP, IRREPWAR)
            IRREPBO = DIRPRD(IRREP, IRREPX)
C
            ITOFF1(IRREPM) = ITOFF1(IRREPM - 1) + VRT1(IRREPAO)*
     &                       POP2(IRREP)
            ITOFF2(IRREPM) = ITOFF2(IRREPM - 1) + POP2(IRREP)*
     &                       VRT2(IRREPBO)
C
 5000    CONTINUE
C
         DO 90 IRREPBO = 1, NIRREP
C
            IRREPM = DIRPRD(IRREPBO, IRREPX)
            IRREPAO = DIRPRD(IRREPM, IRREPWAR)
C     
            NOCCM = POP2(IRREPM)
            NVRTB = VRT2(IRREPBO)
            NVRTA = VRT1(IRREPAO)
C
            JOFFZ = ITOFF1(IRREPM)
            IOFF  = ITOFF2(IRREPM)
C     
            IF (NVRTB .NE. 0 .AND. NOCCM .NE. 0 .AND. NVRTA .NE. 0)
     &      THEN
C     
               CALL XGEMM('N', 'T', NUMSYZ*NVRTA, NVRTB, NOCCM, ONE,
     &                     ZT(1,JOFFZ), NUMSYZ*NVRTA, TB(IOFF), NVRTB, 
     &                     AZERO, TT(1,JOFFT), NUMSYZ*NVRTA)
C     
            ENDIF
C     
C Update pointers
C     
            JOFFT = JOFFT + NVRTA*NVRTB
C
 90      CONTINUE
C     
         CALL TRANSP(TT, Z, DISSYZ, NUMSYZ) 
C
C At the moment we do not have any Q(Ab,Ij) contributions stored 
C in the list. In our case we do not have to include the bare 
C integrals either. So let's put the first contribution of Q(Ab,Ij)
C to the disk. Fist negate the contribution to take care the overall
C sign.
C
         IF (IUHF .EQ. 1) THEN
            CALL VMINUS(Z, NUMSYZ*DISSYZ)
            CALL PUTLST(Z, 1, NUMSYZ, 1, IRREPQR, LISTZ)
         ENDIF
C
      ENDIF
C     
      IF (IUHF .EQ. 0) THEN
C     
C RHF case
C In RHF this is simply a transposition,
C     
         IF (MIN(NUMSYWA, DISSYWA) .EQ. 0) RETURN
C
C Copy Z vector into T before we transpose the indices.
C     
         CALL SCOPY(DISSYZ*NUMSYZ, Z, 1, T, 1)
C
         CALL SYMTR1(IRREPQR, POP1, POP2, DISSYZ, Z, TMP, 
     &               TMP(1 + DISSYZ), TMP(1 + 2*DISSYZ))
         CALL SYMTR3(IRREPQL, VRT1, VRT2, DISSYZ, NUMSYZ, Z, TMP,
     &               TMP(1 + NUMSYZ), TMP(1 + 2*NUMSYZ))
C
C Now add the transposition to the original and negate and then dump into 
C the disk.
C
         CALL VADD(T, T, Z, NUMSYZ*DISSYZ, ONE)
         CALL VMINUS(T, NUMSYZ*DISSYZ)
         CALL PUTLST(T, 1, NUMSYZ, 1, IRREPQR, LISTZ)
C
         IF (IFLAGS(1) .GE. 20) THEN
            NSIZE = NUMSYZ*DISSYZ
            CALL HEADER('Checksum @-GABEFIND per sym. block', 0, LUOUT)
            WRITE(LUOUT, *) ' ABAB =  ', SDOT(NSIZE, T, 1, T, 1)
         ENDIF
C
      ELSE
C     
C If there are no integrals skip multiplication
C     
         IF (MIN(NUMSYWB, DISSYWB) .EQ. 0) RETURN
C     
C Decide about the algorithm 
C     
         CALL ZERO(ZT, NUMSYT*NUMSYWB)
C     
C Get T2 amplitudes
C     
         CALL GETLST(T, 1, NUMSYT, 1, IRREPTR, LISTT)
C     
C Process as many Ef distributions at once as possible
C     
C Get integrals Hbar(Fe,Ma) from LISTWB
C     
         NINCOR = MAXSIZE/(DISSYWB*IINTFP)
         NLEFT  = NUMSYWB
         NFIRST = 1
         NPASS  = 0
C
 2       NREAD  = MIN(NLEFT, NINCOR)
C
         CALL GETLST(W, NFIRST, NREAD, 1, IRREPWBR, LISTWB)
C     
C Perform multiplication
C            +
C     T(Fe,Ij)*Hbar(Fe,Ma) = Z(Ij,Ma)
C     
         CALL XGEMM ('T', 'N', NUMSYT, NREAD, DISSYT, ONE, T,
     &                DISSYT, W, DISSYWB, AZERO, ZT(1,NFIRST),
     &                NUMSYZ)
C
         NFIRST = NFIRST + NREAD
         NLEFT  = NLEFT - NREAD
         NPASS  = NPASS + 1
         IF (NLEFT .NE. 0) GOTO 2
C     
C Do the second part of the multiplication
C     
         CALL ZERO(TT, NUMSYZ*DISSYZ)
         CALL SYMTR1(IRREPWBR, POP1, VRT2, NUMSYZ, ZT, TMP,
     &               TMP(1 + NUMSYZ), TMP(1 + 2*NUMSYZ))
C     
         JOFFT = 1
         IZOFF1(1) = 1
         IZOFF2(1) = 1 
C         
C Setup an offset for the T2 array
C
         DO 5100 IRREPM = 2, NIRREP
C
            IRREP = IRREPM - 1
            IRREPBO = DIRPRD(IRREP, IRREPWBR)
            IRREPAO = DIRPRD(IRREP, IRREPX)
C
            IZOFF1(IRREPM) = IZOFF1(IRREPM - 1) + VRT2(IRREPBO)*
     &                       POP1(IRREP)
            IZOFF2(IRREPM) = IZOFF2(IRREPM - 1) + POP1(IRREP)*
     &                       VRT1(IRREPAO)
 5100    CONTINUE
C
         DO 190 IRREPAO = 1, NIRREP
C
            IRREPM = DIRPRD(IRREPAO, IRREPX)
            IRREPBO = DIRPRD(IRREPM, IRREPWBR)
C
            NOCCM  = POP1(IRREPM)
            NVRTB  = VRT2(IRREPBO)
            NVRTA  = VRT1(IRREPAO)
C
            JOFFZ = IZOFF1(IRREPM)
            IOFF  = IZOFF2(IRREPM)
C
            IF (NVRTB .NE. 0 .AND. NOCCM .NE. 0 .AND. NVRTA .NE. 0) 
     &      THEN  
C     
               CALL XGEMM('N', 'T', NUMSYZ*NVRTB, NVRTA, NOCCM, ONE, 
     &                     ZT(1,JOFFZ), NUMSYZ*NVRTB, TA(IOFF), NVRTA,
     &                     ONE, TT(1,JOFFT), NUMSYZ*NVRTB)
            ENDIF
C
            JOFFT = JOFFT + NVRTB*NVRTA
C
 190     CONTINUE
C
         CALL SYMTRA(IRREPQL, VRT2, VRT1, NUMSYZ, TT, ZT)
         CALL TRANSP(ZT, T, DISSYZ, NUMSYZ)
C
C Now we have to update the List 63 with this contribution,
C since this term has the same structure as ones already
C in list 63 (Q(Ab,Ij)). 
         
         CALL GETLST(Z, 1, NUMSYZ, 2, IRREPQR, LISTZ)
         CALL VADD(Z, Z, T, NUMSYZ*DISSYZ, ONEM)
C
         IF (IFLAGS(1) .GE. 20) THEN
C
            NSIZE = NUMSYZ*DISSYZ
C
            CALL HEADER('Checksum @-GABEFIND per sym. block', 0, LUOUT)
C            
            WRITE(LUOUT, *) ' ABAB =  ', SDOT(NSIZE, Z, 1, Z, 1)
         ENDIF
C
         CALL PUTLST(Z, 1, NUMSYZ, 1, IRREPQR, LISTZ)
C
      ENDIF
C
      RETURN
      END

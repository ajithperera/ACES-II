C
      SUBROUTINE MKDBLAA1(Z, ZT, T, TT, W, MAXSIZE, TA, POP, VRT,
     &                    DISSYZ, DISSYW, DISSYT, NUMSYZ, NUMSYW,
     &                    NUMSYT, NVRTSQ, NTASIZ, LISTT, LISTZ, LISTW, 
     &                    IRREPTL, IRREPTR, IRREPWL, IRREPWR, IRREPQL, 
     &                    IRREPQR, IRREPX, IUHF, ISPIN, TMP)
C
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      INTEGER DISSYZ,DISSYW,DISSYT,DIRPRD,POP,VRT
      CHARACTER*8 SPCASE(2)
C
      DIMENSION Z(DISSYZ,NUMSYZ),W(DISSYW,1),TA(NTASIZ)
      DIMENSION T(DISSYT,NUMSYT)
      DIMENSION ZT(NUMSYZ,DISSYZ),TT(NUMSYT,DISSYT)
      DIMENSION TMP(1),POP(8),VRT(8),ITOFF1(8),ITOFF2(8)
C
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /FILES/ LUOUT, MOINTS
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
      DATA AZERO,ONE,ONEM /0.0D0,1.0D0 ,-1.0D0/
      DATA SPCASE /'AAAA =  ', 'BBBB =  '/
C
C Zero first output array (Neccesary because this is not necessarily
C done in the matrix multiplication
C
      CALL ZERO(Z, NUMSYZ*DISSYZ)
C
C If there are no W integrals there is nothing to do.
C     
      IF (MIN(NUMSYW, DISSYW) .NE. 0) THEN
C     
C Get T2 amplitudes (for quadratic term there is no need to form 
C the TAU amplitudes)
C
         CALL GETLST(T, 1, NUMSYT, 1, IRREPTR, LISTT)
C     
C Process as many distributions as possible
C Get integrals Hbar(E<F,AM) from LISTW
C     
         NINCOR = MAXSIZE/(DISSYW*IINTFP)
         NLEFT  = NUMSYW
         NFIRST = 1
         NPASS  = 0
C
 1       NREAD = MIN(NLEFT, NINCOR)
         CALL GETLST(W, NFIRST, NREAD, 1, IRREPWR, LISTW)
C     
C     Perform multiplication
C
C   T(E<F,I<J)*Hbar(E<F,AM) = Z(I<J,AM)
C     
         CALL XGEMM ('T', 'N', NUMSYT, NREAD, DISSYT, ONE, T,
     &                DISSYT, W, DISSYW, AZERO, ZT(1,NFIRST),
     &                NUMSYT)
C
         NFIRST = NFIRST + NREAD
         NLEFT  = NLEFT - NREAD
         NPASS  = NPASS + 1
         IF (NLEFT .NE. 0) GOTO 1
C     
C Do the second part of the multiplication 
C     
         CALL ZERO(TT, NUMSYZ*NVRTSQ)
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
            IRREPAO = DIRPRD(IRREP, IRREPWR)
            IRREPBO = DIRPRD(IRREP, IRREPX)
C
            ITOFF1(IRREPM) = ITOFF1(IRREPM - 1) + VRT(IRREPAO)*
     &                       POP(IRREP)
            ITOFF2(IRREPM) = ITOFF2(IRREPM -1) + POP(IRREP)*
     &                       VRT(IRREPBO)
C
 5000    CONTINUE
C
         DO 90 IRREPBO = 1, NIRREP
C
            IRREPM = DIRPRD(IRREPBO, IRREPX)
            IRREPAO = DIRPRD(IRREPM, IRREPWR)
C     
            NOCCM = POP(IRREPM)
            NVRTB = VRT(IRREPBO)
            NVRTA = VRT(IRREPAO)
C
            IOFF  = ITOFF2(IRREPM)
            JOFFZ = ITOFF1(IRREPM)
C     
            IF(NVRTB .NE.0 .AND. NOCCM .NE. 0 .AND. NVRTA .NE. 0) THEN
C     
               CALL XGEMM('N', 'T', NUMSYZ*NVRTA, NVRTB, NOCCM, ONE, 
     &                     ZT(1,JOFFZ), NUMSYZ*NVRTA, TA(IOFF), NVRTB,
     &                     AZERO,TT(1,JOFFT), NUMSYZ*NVRTA)
C     
            ENDIF
C     
C Update pointers
C     
            JOFFT = JOFFT + NVRTB*NVRTA
C
 90      CONTINUE
C     
         CALL ASSYM(IRREPQL, VRT, NUMSYZ, NUMSYZ, ZT, TT)
         CALL SCOPY(DISSYZ*NUMSYZ, ZT, 1, TT, 1)
         CALL TRANSP(TT, Z, DISSYZ, NUMSYZ) 
C
C At the moment we do not have any Q(A<B, I<J) contributions stored in the list.
C In our case we do not have to include the bare integrals either.
C So let's put the first contribution of Q(A<B, I<J) to the disk. First 
C we have to negate the contribution to take care of the sign.
C
         CALL VMINUS(Z, NUMSYZ*DISSYZ)
         CALL PUTLST(Z, 1, NUMSYZ, 1, IRREPQR, LISTZ)
C
         IF (IFLAGS(1) .GE. 20) THEN
C
            NSIZE = NUMSYZ*DISSYZ
            CALL HEADER('Checksum @-GABEFIND per sym. block', 0, LUOUT)
C            
            WRITE(LUOUT, *) SPCASE(ISPIN), SDOT(NSIZE, Z, 1, Z, 1)
C
         ENDIF
C
      ENDIF
C     
      RETURN
      END

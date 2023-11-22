C
      SUBROUTINE QTHRBAB1(YNI, T2, Q, NFOOY, DISSYT, NUMSYT, DISSYQ,
     &                    NUMSYQ, IRREPTR, IRREPQR, LIST2, LISTQ, 
     &                    POP1, POP2, VRT1, VRT2, TMP, IUHF, ISPIN)
C
C Take the product Q(Ab,Ij) = T2(Ab,Nj)*Y(N,I) in the calculation
C of three body contribution to the quadratic term.
C
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,ONEM,ZILCH,YNI,TMP, T2, Q, SDOT
      DIMENSION YNI(NFOOY), T2(DISSYT, NUMSYT), Q(DISSYQ, NUMSYQ)
      DIMENSION POP1(8), VRT1(8), TMP(1), POP2(8), VRT2(8)
C      
      COMMON /INFO/ NOCA,NOCB,NVRTA,NVRTB
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      COMMON /SYM/ POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF1BB,
     &             NF2AA,NF2BB
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
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
      DATA ONE /1.0D+00/ 
      DATA ONEM /-1.0D+00/
      DATA ZILCH /0.0D+00/
C
C Load the T2 vector in to the memory. T2 vector is ordered
C as T2(Ab,Nj) [ISPIN = 1] and T2(Ab,In) [ISPIN = 2]
C
      CALL GETLST(T2, 1, NUMSYT, 1, IRREPTR, LIST2)
C
C Now interchange the last two indices. [T2(Ab,Nj) ---> T2(Ab,jN)]
C in the case of [ISPIN = 1]. There is no need to do that for
C [ISPIN = 2].
C  
      IF (ISPIN .EQ. 1) THEN
         CALL SYMTR1(IRREPTR, POP1, POP2, DISSYT, T2, TMP, 
     &               TMP(1 + DISSYT), TMP(1 + 2*DISSYT))
      ENDIF
C
C Now carry out the multiplication, Q(Ab,Ij)= T2(Ab,jN)*Y(N,I) [ISPIN = 1]
C and Q(Ab,Ij) = T2(AB,In)*Y(n,i) [ISPIN = 2].
C  
      CALL ZERO(Q, DISSYQ*NUMSYQ)
C
      IOFFT1 = 1
      IOFFY1 = 1
      IOFFQ1 = 1
C
      DO 10 IRREPN = 1, NIRREP
C
         IRREPYN = IRREPN
         IRREPYI = IRREPYN
         IRREPTN = IRREPYN
         IRREPTJ = DIRPRD(IRREPTN, IRREPTR)
C     
         NOCCYN = POP1(IRREPYN)
         NOCCYI = POP1(IRREPYI)
C
         NOCCTN = POP1(IRREPTN)
         NOCCTJ = POP2(IRREPTJ)
C
         IF (NOCCYN .GT. 0 .AND. NOCCTJ .GT. 0) THEN
C
            CALL XGEMM('N', 'N', DISSYT*NOCCTJ, NOCCYI, NOCCTN,
     &                  ONE, T2(1, IOFFT1), DISSYT*NOCCTJ,
     &                  YNI(IOFFY1), NOCCYN, ZILCH, Q(1, IOFFQ1),
     &                  DISSYQ*NOCCTJ)
         ENDIF
C
C Update the offsets
C
         IOFFT1 = IOFFT1 + NOCCTJ*NOCCTN
         IOFFQ1 = IOFFQ1 + NOCCTJ*NOCCYI
         IOFFY1 = IOFFY1 + NOCCYI*NOCCYN
C     
 10   CONTINUE
C
C The resulting product is ordered as Q(Ab,jI) [ISPIN = 1] 
C and Q(Ab,Ij) [ISPIN = 2]. Now switch ordering for [ISPIN = 1].
C
      IF (ISPIN .EQ. 1) THEN
         CALL SYMTR1(IRREPQR, POP2, POP1, DISSYQ, Q, TMP,
     &              TMP(1 + DISSYQ), TMP(1 + 2*DISSYQ))
      ENDIF
C
C RHF cases transpose the I and J indices and add the resulting
C contribution.
C
      IF (IUHF .EQ. 0 .AND. ISPIN .EQ. 1) THEN
         CALL SYMRHF(IRREPQR, VRT1, POP1, DISSYQ, Q, TMP,
     &               TMP(1 + DISSYQ), TMP(1 + 2*DISSYQ))

      ENDIF
C
C Now we have what we want. Q(Ab,Ij) including the antisymmetrized
C piece which is the term caculate when [ISPIN = 2].
C Negate the whole contribution (This is the sign of the entire
C contribution)
C 
      CALL VMINUS(Q, DISSYQ*NUMSYQ)
C
C Update the Q(IJ,AB) contribution in the list.
C
      IF (IFLAGS(1) .GE. 20) THEN
C
         NSIZE = NUMSYQ*DISSYQ
C
         CALL HEADER('Checksum @-QTHRBDY-1 per sym. block', 0, LUOUT)
C            
         WRITE(LUOUT, *) ' ABAB =  ', SDOT(NSIZE, Q, 1, Q, 1)
C
      ENDIF
C      
      CALL SUMSYM2(Q, T2, DISSYQ*NUMSYQ, 1, IRREPQR, LISTQ)
C
      RETURN
      END

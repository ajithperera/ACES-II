C
      SUBROUTINE QTHRBAA1(YNI, T2, Q, NFOOY, DISSYT, NUMSYT, DISSYQ,
     &                    NUMSYQ, NOCC2SQ, IRREPTR, IRREPQR, LIST2,
     &                    LISTQ, POP1, VRT1, TMP, IUHF, ISPIN)
C
C Take the product Q(AB,IJ) = T2(AB,NJ)*Y(N,I) in the calculation
C of three body contribution to the quadratic term.
C
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,ONEM,ZILCH,YNI,TMP, T2, Q, SDOT
      DIMENSION YNI(NFOOY),T2(DISSYT, NOCC2SQ),Q(DISSYQ, NOCC2SQ)
      DIMENSION POP1(8), VRT1(8), TMP(1)
      CHARACTER*8 SPCASE(2)
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
      DATA SPCASE /'AAAA =  ', 'BBBB =  '/
C
C Load the T2 vector in to the memory. T2 vector is ordered
C as T2(A<B,N<J). 
C
      CALL GETLST(T2, 1, NUMSYT, 1, IRREPTR, LIST2)
C
C Transform T2(A<B,N<J) --> T2(A<B,NJ) (Uncompress the left hand side
C of the T2 array.
C
      CALL SYMEXP(IRREPTR, POP1, DISSYT, T2)
C
C Now interchange the last two indices. [T2(A<B,NJ) ---> T2(A<B,JN)]
C
      CALL SYMTR1(IRREPTR, POP1, POP1, DISSYT, T2, TMP, 
     &            TMP(1 + DISSYT), TMP(1 + 2*DISSYT))
C
C Now carry out the multiplication, Q(A<B,JI)= T2(A<B,JN)*Y(N,I)
C [ISPIN = 1] and Q(a<b,ji) = T2(a<b,jn)*Y(n,i) [ISPIN = 2].
C
      CALL ZERO(Q, DISSYQ*NOCC2SQ)
C
      IOFFT = 1
      IOFFY = 1
      IOFFQ = 1
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
         NOCCTJ = POP1(IRREPTJ)
         NOCCTN = POP1(IRREPTN)
C
         IF (NOCCTJ .GT. 0 .AND. NOCCYN .GT. 0) THEN
      
            CALL XGEMM('N', 'N', DISSYT*NOCCTJ, NOCCYI, NOCCTN, ONE, 
     &                  T2(1, IOFFT), DISSYT*NOCCTJ, YNI(IOFFY),
     &                  NOCCYN, ZILCH, Q(1, IOFFQ), DISSYQ*NOCCTJ)
         ENDIF
C
C Update the offsets
C
         IOFFT = IOFFT + NOCCTJ*NOCCTN
         IOFFQ = IOFFQ + NOCCTJ*NOCCYI
         IOFFY = IOFFY + NOCCYN*NOCCYI
C
 10   CONTINUE
C
C The resulting product is ordered as Q(A<B,JI) [ISPIN = 1] or
C Q(a<b,ji) [ISPIN = 2]. Switch the two left hand indices.
C
      CALL SYMTR1(IRREPQR, POP1, POP1, DISSYQ, Q, TMP,
     &            TMP(1 + DISSYQ), TMP(1 + 2*DISSYQ))
C
C      IF (IFLAGS(1) .GE. 20) THEN
C     
C         NSIZE = NUMSYQ*DISSYQ
C         CALL HEADER('Checksum @-QTHRBDY11 per sym. block', 0, LUOUT)
C            
C         WRITE(LUOUT, *) SPCASE(ISPIN), SDOT(NSIZE, Q, 1, Q, 1)
C      ENDIF
C
C Now we have something close to what we want. Q(A<B,IJ) or Q(a<b,ij) 
C depending the value of ISPIN. Now antisymmetrize the (P_(IJ)
C in algebraic expressions)
C
      CALL ASSYM(IRREPQR, POP1, DISSYQ, DISSYQ, T2, Q)
C
C Negate the whole contribution (This is the sign of the entire
C contribution)
C 
      CALL VMINUS(T2, DISSYQ*NUMSYQ)
C
C Update the Q(I<J,A<B) or Q(i<j,a<b) contribution in the list.
C
      IF (IFLAGS(1) .GE. 20) THEN
C     
         NSIZE = NUMSYQ*DISSYQ
         CALL HEADER('Checksum @-QTHRBDY1 per sym. block', 0, LUOUT)
C            
         WRITE(LUOUT, *) SPCASE(ISPIN), SDOT(NSIZE, T2, 1, T2, 1)
      ENDIF
C
      CALL SUMSYM2(T2, Q, DISSYQ*NUMSYQ, 1, IRREPQR, LISTQ)
C
C      CALL CHKSUM(T2, Q, DISSYQ*NUMSYQ, IRREPQR, LISTQ, SPCASE(ISPIN),
C     &            'QTHRBAA1')

C
      RETURN
      END

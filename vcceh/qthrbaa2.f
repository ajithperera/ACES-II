C
      SUBROUTINE QTHRBAA2(YFA, T2, Q, NFVVY, DISSYT, NUMSYT, DISSYQ,
     &                    NUMSYQ, NVRT2SQ, IRREPTL, IRREPTR, IRREPQL, 
     &                    IRREPQR, LIST2, LISTQ, POP1, VRT1, TMP,
     &                    IUHF, ISPIN)
C
C Take the product Q(AB,IJ) = T2(FB,IJ)*Y(F,A) in the calculation
C of three body contribution to the quadratic term.
C
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,ONEM,ZILCH,YFA,TMP, T2, Q, SDOT
      DIMENSION YFA(NFVVY),T2(NUMSYT,NVRT2SQ),Q(NUMSYQ, NVRT2SQ)
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
C as T2(F<B,I<J). 
C
      CALL GETLST(T2, 1, NUMSYT, 1, IRREPTR, LIST2)
C
C Take the transpose of the T2 vector. Now the new ordering of the T2
C vector is T2(I<J,F<B). Q has the same size as T2 and return the
C transpose matrix. Copy it back to the T2 vector.
C 
      CALL TRANSP(T2, Q, NUMSYT, DISSYT)
      CALL SCOPY(NUMSYT*DISSYT, Q, 1, T2, 1)
C
C Transform T2(I<J,F<B) --> T2(I<J,FB) (Uncompress the left hand side
C of the T2 array.
C
      CALL SYMEXP(IRREPTL, VRT1, NUMSYT, T2)
C
C Now interchange the last two indices. [T2(I<J,FB) ---> T2(I<J,BF)]
C
      CALL SYMTR1(IRREPTL, VRT1, VRT1, NUMSYT, T2, TMP, 
     &            TMP(1 + NUMSYT), TMP(1 + 2*NUMSYT))
C
C Now carry out the multiplication, Q(I<J,BA) = T2(I<J,BF)*Y(F,A) [ISPIN = 1]
C and Q(i<J,ba) = T2(i<j,bf)*Y(f,a) [ISPIN=2].
C
      IOFFT = 1
      IOFFY = 1
      IOFFQ = 1
C
      CALL ZERO(Q, NVRT2SQ*NUMSYQ)
C
      DO 10 IRREPF = 1, NIRREP
C
         IRREPYF = IRREPF
         IRREPYA = IRREPYF
         IRREPTF = IRREPYF
         IRREPTB = DIRPRD(IRREPTF, IRREPTL)
C     
         NVRTYF = VRT1(IRREPYF)
         NVRTYA = VRT1(IRREPYA)
C
         NVRTTB = VRT1(IRREPTB)
         NVRTTF = VRT1(IRREPTF)
C
         IF (NVRTTB .GT. 0 .AND. NVRTYF .GT. 0) THEN
C
            CALL XGEMM('N', 'N', NUMSYT*NVRTTB, NVRTYA, NVRTTF, ONE, 
     &                  T2(1, IOFFT), NUMSYT*NVRTTB, YFA(IOFFY),
     &                  NVRTYF, ZILCH, Q(1, IOFFQ), NUMSYT*NVRTTB)
         ENDIF
C     
C Update the offsets
C
         IOFFT = IOFFT + NVRTTB*NVRTTF
         IOFFQ = IOFFQ + NVRTTB*NVRTYA
         IOFFY = IOFFY + NVRTYF*NVRTYA
C
 10   CONTINUE
C
C The resulting product is ordered as Q(I<J,BA) [ISPIN = 1] or
C Q(i<j,ba) [ISPIN = 2]. Switch the two left hand indices.
C
      CALL SYMTR1(IRREPQL, VRT1, VRT1, NUMSYQ, Q, TMP,
     &            TMP(1 + NUMSYQ), TMP(1 + 2*NUMSYQ))
C
C Now we have something close to what we want. Q(I<J,AB) or Q(i<j,ab) 
C depending the value of ISPIN. Now antisymmetrize the (P_(IJ)
C in algebraic expressions)
C
      CALL ASSYM(IRREPQL, VRT1, NUMSYQ, NUMSYQ, T2, Q)
C
C Negate the whole contribution (This is the sign of the entire
C contribution) and take the transpose to have what we want, 
C Q(I<J,A<B) or Q(i<j,a<b) depending on ISPIN.
C 
      CALL VMINUS(T2, DISSYQ*NUMSYQ)
      CALL TRANSP(T2, Q, DISSYQ, NUMSYQ)
C
C Update the Q(I<J,A<B) or Q(i<j,a<b) contribution in the list.
C
      IF (IFLAGS(1) .GE. 20) THEN
C     
         NSIZE = NUMSYQ*DISSYQ
         CALL HEADER('Checksum @-QTHRBDY2 per sym. block', 0, LUOUT)
C            
         WRITE(LUOUT, *) SPCASE(ISPIN), SDOT(NSIZE, Q, 1, Q, 1)
      ENDIF
C
      CALL SUMSYM2(Q, T2, DISSYQ*NUMSYQ, 1, IRREPQR, LISTQ)
C     
      RETURN
      END

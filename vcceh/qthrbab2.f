C
      SUBROUTINE QTHRBAB2(YFA, T2, Q, NFVVY, DISSYT, NUMSYT, DISSYQ,
     &                    NUMSYQ, IRREPTL, IRREPTR, IRREPQL, IRREPQR,
     &                    LIST2, LISTQ, POP1, POP2, VRT1, VRT2, TMP, 
     &                    IUHF, ISPIN)
C
C Take the product Q(Ab,Ij) = T2(Fb,Ij)*Y(F,A) in the calculation
C of three body contribution to the quadratic term.

      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,ONEM,ZILCH,YFA,TMP, T2, Q, SDOT
      DIMENSION YFA(NFVVY),T2(NUMSYT, DISSYT),Q(NUMSYQ, DISSYQ)
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
C as T2(Fb,Ij) [ISPIN = 1] and T2(Af,Ij) [ISPIN=2]
C
      CALL GETLST(T2, 1, NUMSYT, 1, IRREPTR, LIST2)
C
C Take the transpose of the T2 vector. Now the new ordering of the T2
C vector is T2(Ij,Fb) [ISPIN = 1] and T2(Ij,Af) [ISPIN = 2]. Q has the same
C size as T2 and return the transpose matrix. Copy it back to the T2 vector.
C 
      CALL TRANSP(T2, Q, NUMSYT, DISSYT)
      CALL SCOPY(NUMSYT*DISSYT, Q, 1, T2, 1)
C
C Now interchange the last two indices. [T2(Ij,Fb) ---> T2(Ij,bF)]
C in the case of [ISPIN = 1]. There is no need to do that for
C [ISPIN = 2].
C  
      IF (ISPIN .EQ. 1) THEN
         CALL SYMTR1(IRREPTL, VRT1, VRT2, NUMSYT, T2, TMP, 
     &               TMP(1 + NUMSYT), TMP(1 + 2*NUMSYT))
      ENDIF
C
C Now carry out the multiplication, Q(Ij,bA)= T2(Ij,bF)*Y(F,A) [ISPIN = 1]
C and Q(Ij,Ab) = T2(Ij,Af)*Y(f,b) [ISPIN = 2].
C
      CALL ZERO(Q, DISSYQ*NUMSYQ)
C
      IOFFT = 1
      IOFFY = 1
      IOFFQ = 1 
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
         NVRTTF = VRT1(IRREPTF)
         NVRTTB = VRT2(IRREPTB)
C
         IF (NVRTTB .GT. 0 .AND. NVRTYF .GT. 0) THEN
C
            CALL XGEMM('N', 'N', NUMSYT*NVRTTB, NVRTYA, NVRTTF, ONE, 
     &                  T2(1, IOFFT), NUMSYT*NVRTTB, YFA(IOFFY),
     &                  NVRTYF,  ZILCH, Q(1, IOFFQ), NUMSYT*NVRTTB)
         ENDIF
C
C Update the offsets
C
         IOFFT = IOFFT + NVRTTB*NVRTTF
         IOFFQ = IOFFQ + NVRTYA*NVRTTB
         IOFFY = IOFFY + NVRTYF*NVRTYA
C
 10   CONTINUE
C
C The resulting product is ordered as Q(Ij,bA) [ISPIN = 1] 
C and Q(Ij,Ab) [ISPIN = 2]. Now switch ordering for [ISPIN = 2].
C
      IF (ISPIN .EQ. 1) THEN
         CALL SYMTR1(IRREPQL, VRT2, VRT1, NUMSYQ, Q, TMP, 
     &               TMP(1 + NUMSYQ), TMP(1 + 2*NUMSYQ))
      ENDIF
C
C RHF case add the transpose of I and J to the contribution
C
      IF (IUHF .EQ. 0 .AND. ISPIN .EQ. 1) THEN
         CALL SYMRHF(IRREPQL, POP1, VRT1, NUMSYQ, Q, TMP, 
     &               TMP(1 + NUMSYQ), TMP(1 + 2*NUMSYQ))
      ENDIF
C
C Now we have what we want. Q(Ab,Ij) including the antisymmetrized
C piece which is the term caculate when [ISPIN = 2].
C Negate the whole contribution (This is the sign of the entire
C contribution) and take the transpose.
C
      CALL VMINUS(Q, DISSYQ*NUMSYQ)
      CALL TRANSP(Q, T2, DISSYQ, NUMSYQ)
C
C Update the Q(IJ,AB) contribution in the list.
C
      IF (IFLAGS(1) .GE. 20) THEN

         NSIZE = NUMSYQ*DISSYQ

         CALL HEADER('Checksum @-QTHRBDY-2 per sym. block', 0, LUOUT)
            
         WRITE(LUOUT, *) ' ABAB =  ', SDOT(NSIZE, T2, 1, T2, 1)

      ENDIF
C
      CALL SUMSYM2(T2, Q, DISSYQ*NUMSYQ, 1, IRREPQR, LISTQ)
C
      RETURN
      END

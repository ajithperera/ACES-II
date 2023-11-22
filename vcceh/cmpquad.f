C
      SUBROUTINE CMPQUAD(ICORE, MAXCOR, QUAD, MAXPERT, IRREPX, IOFFSET, 
     &                   FACT, IUHF)
C
C Driver for the calculation of the full quadrtic contribution for given
C  A                 A,B                                       A   A,B
C Q(ALPHA,BETA) and Q(APHA,BETA) quadratic intermediates and  L , L    
C  I                 I,J                                       I   I,J
C lambda amplitudes.
C
C Algebraic expression for the quadratic contribution is as 
C follows
C                              A,B        A,B           A           A
C QUAD(ALPHA, BETA) = 1/4 SUM Q(ALP,BET)*L      +  SUM Q (ALP,BET)*L
C                         I,J  I,J        I,J      A,I  I           I
C                         A,B
C                                     A        B
C                   + 1/4 SUM [P_(AB)T(ALPHA)*T(BETA)]*Hbar(IJ,AB)
C                         I,J         I        J
C                         A,B
C Spin integrated formulas are as follows. Here alpha and beta is 
C implicit.
C
C UHF:
C
C  QUAD = SUM (A<B,I<J) Q(AB,IJ)*L(AB,IJ) + SUM (A,I) Q(A,I)*L(A,I)
C       + SUM (A<B,I<J) P_(AB)[T(A,I)*T(B,J)]*Hbar(IJ,AB)             [AAAA]
C
C  QUAD = SUM (a<b,i<j) Q(ab,ij)*L(ab,ij) + SUM (a,i) Q(a,i)*L(a,i)
C       + SUM (a<b,i<j) P_(ab)[T(a,i)*T(b,j)]*Hbar(ij,ab)             [BBBB]
C
C  QUAD = SUM (Ab,Ij) Q(Ab,Ij)*L(Ab,Ij) + T(A,I)*T(b,j)*Hbar(Ij,Ab)   [ABAB]
C
C RHF:
C
C  QUAD = 2[SUM (A<B,I<J) Q(AB,IJ)*L(AB,IJ) + SUM (A,I) Q(A,I)*L(A,I)
C       + SUM (A<B,I<J) P_(AB)[T(A,I)*T(B,J)]*Hbar(IJ,AB)]            [AAAA]
C
C  QUAD = SUM (Ab,Ij) Q(Ab,Ij)*L(Ab,Ij) + T(A,I)*T(b,j)*Hbar(Ij,Ab)   [ABAB]
C
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION FACTOR, QUAD, EDBL, ONE, ETOT, ESING, SDOT,
     &                 ONEM, FACT
C
      DIMENSION ICORE(MAXCOR), I0Q(2), I0L(2), IOFFT1A(8,2),
     &          IOFFT1B(8,2)
C
      CHARACTER*4 SPCASE(3)
      CHARACTER*6 SPCASE1(2)
      CHARACTER*8 SPCASE2(2)
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NF1(2),NF2(2)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
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
      DATA SPCASE /'AAAA', 'BBBB', 'ABAB'/
      DATA SPCASE1 /'AA =  ', 'BB =  '/
      DATA SPCASE2 /'AAAA =  ', 'BBBB =  '/
      DATA ONEM/-1.00D+00/
C
      MXCOR = MAXCOR
C
C Allocate memory for Q(A,I) ampltudes and load them in to the memory.
C Note that now Q(A,I) is totally symmetric.
C
      I0Q(1) = MXCOR + 1 - NT(1)*IINTFP
      MXCOR  = MXCOR - NT(1)*IINTFP
C
      CALL GETLST(ICORE(I0Q(1)), 1, 1, 1, 3, 90)
C
      IF (IUHF .EQ. 0) THEN
         I0Q(2) = I0Q(1)
      ELSE
         I0Q(2) = I0Q(1) - NT(2)*IINTFP
         MXCOR  = MXCOR - NT(2)*IINTFP
C
         CALL GETLST(ICORE(I0Q(2)), 1, 1, 1, 4, 90)
      ENDIF
C
C Allocate memory for L(A,I) and load them in to the memory.
C
      I0L(1) = I0Q(2) - NT(1)*IINTFP
      MXCOR  = MXCOR - NT(1)*IINTFP
C
      CALL GETLST(ICORE(I0L(1)), 1, 1, 1, 1, 190)
C
      IF (IUHF .EQ. 0) THEN  
         I0L(2) = I0L(1)
      ELSE
         I0L(2) = I0L(1) - NT(2)*IINTFP
         MXCOR  = MXCOR - NT(2)*IINTFP
C
         CALL GETLST(ICORE(I0L(2)), 1, 1, 1, 2, 190)
      ENDIF
C
      ETOT = 0.0D00
      ONE  = 1.0D00
C     
      FACTOR = 1.0D+00
      IF(IUHF .EQ. 0) FACTOR = 2.0D+00
C 
      DO 10 ISPIN = 1, IUHF + 1
C
         LISTQ = 60 +  ISPIN
         LISTL = 143 + ISPIN 
C
C Add here the Q(I,A) L(I,A) contribution. There are AA and BB
C contributions.
C
         ESING = SDOT(NT(ISPIN), ICORE(I0Q(ISPIN)), 1,
     &           ICORE(I0L(ISPIN)), 1)
C
         IF (IFLAGS(1) .GE. 20) THEN
            CALL HEADER('Singles contribution @-CMPQUAD', 0, LUOUT)
            WRITE(LUOUT, *) SPCASE1(ISPIN), ESING
         ENDIF
C
         ETOT = ETOT + FACTOR*ESING
C
C Add here Q(AB,IJ)*L(IJ,AB) contribution to the energy for AAAA and 
C BBBB spin cases
C
         IF (IUHF .NE. 0) THEN
C
            DO 100 IRREP = 1, NIRREP
C
               DISSYQ = IRPDPD(IRREP, ISYTYP(1, LISTQ))
               NUMSYQ = IRPDPD(IRREP, ISYTYP(2, LISTQ))
C     
               IF (MIN(NUMSYQ, DISSYQ) .NE. 0) THEN
C
                  I001 = 1
                  I002 = I001 + IINTFP*NUMSYQ*DISSYQ
                  I003 = I002 + IINTFP*NUMSYQ*DISSYQ
                  I004 = I003 + NUMSYQ
               
                  IF (I004 .LT. MXCOR) THEN
C
                     CALL GETLST(ICORE(I001), 1, NUMSYQ, 2, IRREP, 
     &                           LISTQ)
                     CALL GETLST(ICORE(I002), 1, NUMSYQ, 2, IRREP,
     &                           LISTL)
C
                     EDBL = SDOT(NUMSYQ*DISSYQ, ICORE(I001), 1,
     &                      ICORE(I002), 1)
C
                     IF (IFLAGS(1) .GE. 20) THEN
          CALL HEADER('Doubles contribution per sym. block @-CMPQUAD',
     &                 0, LUOUT)                
          WRITE(LUOUT, *) SPCASE2(ISPIN), EDBL
                     ENDIF
C     
                  ELSE
                     CALL INSMEM('CMPQUAD', I004, MXCOR)
                  ENDIF
C     
                  ETOT = ETOT + FACTOR*EDBL
C               
               ENDIF
C
 100        CONTINUE   
C
         ENDIF
C  
 10   CONTINUE
C
C Add here Q(AB,IJ)*L(IJ,AB) contribution to the energy for ABAB
C spin case.
C
      LISTQ = 63
      LISTL = 146
C
      DO 200 IRREP = 1, NIRREP
C
         DISSYQ = IRPDPD(IRREP, ISYTYP(1, LISTQ))
         NUMSYQ = IRPDPD(IRREP, ISYTYP(2, LISTQ))
C
         IF (MIN(NUMSYQ, DISSYQ) .NE. 0) THEN
C            
            I001 = 1
            I002 = I001 + IINTFP*NUMSYQ*DISSYQ
            I003 = I002 + IINTFP*NUMSYQ*DISSYQ
            I004 = I003 + IINTFP*MAX(DISSYQ, NUMSYQ)
            I005 = I004 + MAX(DISSYQ, NUMSYQ)
C
            IF (I005 .LT. MXCOR) THEN
C
               CALL GETLST(ICORE(I001), 1, NUMSYQ, 2, IRREP, LISTQ)
               CALL GETLST(ICORE(I002), 1, NUMSYQ, 2, IRREP, LISTL)
C
C Spin adapt for RHF caculations
C
               IF (IUHF .EQ. 0) THEN
                  CALL SPINAD1(IRREP, POP(1, 1), DISSYQ, ICORE(I001),
     &                         ICORE(I003), ICORE(I004))
               ENDIF
C
               EDBL = SDOT(NUMSYQ*DISSYQ, ICORE(I001), 1,
     &                     ICORE(I002), 1)
C
               IF (IFLAGS(1) .GE. 20) THEN
           CALL HEADER('Doubles contribution per sym. block @-CMPQUAD',
     &                  0, LUOUT)                
           WRITE(LUOUT, *) 'ABAB =  ',  EDBL
               ENDIF
C
            ELSE
               CALL INSMEM('CMPQUAD', I004, MXCOR)
            ENDIF
C
            ETOT = ETOT + EDBL
C     
         ENDIF
C
 200  CONTINUE
C
C Now let's do the SUM (A<B,I<J) [P_(AB)T(A,I)*T(B,J)]*Hbar(IJ,AB) term 
C to the quadratic contribution. This is moved out of the above doubles
C loop to avoid memory problems having to store L(A,I) and
C two T(A,I) vectors for both perturbations and quadratic double
C amplitudes and Hbar integrals at the same time.
C
C Allocate memory for T1(ALPHA) and T1(BETA) amplitudes and retrieve
C them from the disk.
C
      LISTA   = IAPRT1AA
      LISTB   = IBPRT1AA
C
      CALL GETPERT1(ICORE, MAXCOR, INTCOR, IUHF, IRREPX, LISTA,
     &              IOFFSET, IOFFT1A)

      CALL GETPERT1(ICORE, INTCOR, MXCOR, IUHF, IRREPX, LISTB,
     &              IOFFSET, IOFFT1B)
C 
C AAAA and BBBB spin cases
C
      IF (IUHF .NE. 0) THEN
         DO 20 ISPIN = 1, IUHF + 1
C         
            LISTH = 13 + ISPIN
        
            DO 300 IRREP = 1, NIRREP
C            
               DISSYH = IRPDPD(IRREP, ISYTYP(1, LISTH))
               NUMSYH = IRPDPD(IRREP, ISYTYP(2, LISTH))
               DISSYQ = DISSYH
               NUMSYQ = NUMSYH
C
               I001 = 1
               I002 = I001 + DISSYH*NUMSYH*IINTFP
               I003 = I002 + DISSYQ*NUMSYQ*IINTFP
C
               IF (MIN(NUMSYQ, DISSYQ) .NE. 0) THEN
C 
                  IF (I003 .LT. MAXCOR) THEN
C     
                     CALL GETLST(ICORE(I001), 1, NUMSYH, 2, IRREP,
     &                           LISTH)
C 
                     CALL MKT1TAU(ICORE(I002),
     &                            ICORE(IOFFT1A(IRREPX, ISPIN)),
     &                            ICORE(IOFFT1B(IRREPX, ISPIN)),
     &                            DISSYQ, NUMSYQ, POP(1, ISPIN),
     &                            POP(1, ISPIN), VRT(1, ISPIN),
     &                            VRT(1, ISPIN), IRREP, IRREPX, 
     &                            IRREPX, ISPIN, ONE)
C
                     EDBL = SDOT(NUMSYQ*DISSYQ, ICORE(I001), 1,
     &                           ICORE(I002), 1)
C     
                     IF (IFLAGS(1) .GE. 20) THEN
          CALL HEADER('Doubles contribution per sym. block @-CMPQUAD',
     &                0, LUOUT)                
          WRITE(LUOUT, *) SPCASE2(ISPIN), EDBL
                     ENDIF
C     
                  ELSE
                     CALL INSMEM('CMPQUAD', I004, MXCOR)
                  ENDIF
C
                  ETOT = ETOT + FACTOR*EDBL
C
               ENDIF
C
 300        CONTINUE
C
 20      CONTINUE
C
      ENDIF
C
C ABAB spin cases 
C
      LISTH = 16
      ISPIN = 3
         
      DO 30 IRREP = 1, NIRREP
C     
         DISSYH = IRPDPD(IRREP, ISYTYP(1, LISTH))
         NUMSYH = IRPDPD(IRREP, ISYTYP(2, LISTH))
         DISSYQ = DISSYH
         NUMSYQ = NUMSYH
C
         I001 = 1
         I002 = I001 + DISSYH*NUMSYH*IINTFP
         I003 = I002 + DISSYQ*NUMSYQ*IINTFP
         I004 = I003 + MAX(DISSYQ, NUMSYQ)*IINTFP
         I005 = I004 + MAX(DISSYQ, NUMSYQ)
C         
         IF(MIN(NUMSYQ, DISSYQ) .NE. 0) THEN
C
            IF (I005 .LT. MAXCOR) THEN
C     
               CALL GETLST(ICORE(I001), 1, NUMSYH, 2, IRREP, LISTH)
C     
               CALL MKT1TAU(ICORE(I002), ICORE(IOFFT1A(IRREPX, 1)), 
     &                      ICORE(IOFFT1B(IRREPX, 2)), DISSYQ, NUMSYQ,
     &                      POP(1, 1), POP(1, 2), VRT(1, 1), VRT(1, 2),
     &                      IRREP, IRREPX, IRREPX, ISPIN, ONE)
C
C Spin adapt for RHF caculations
C
               IF (IUHF .EQ. 0) THEN
                  CALL SPINAD1(IRREP, POP(1, 1), DISSYQ, ICORE(I002),
     &                         ICORE(I003), ICORE(I004))
               ENDIF
C
               EDBL = SDOT(NUMSYQ*DISSYQ, ICORE(I001), 1,
     &                     ICORE(I002), 1)
C
               IF (IFLAGS(1) .GE. 20) THEN
           CALL HEADER('Doubles contribution per sym. block @-CMPQUAD',
     &                  0, LUOUT)                
           WRITE(LUOUT, *) 'ABAB =  ',  EDBL
               ENDIF
C
            ELSE
               CALL INSMEM('CMPQUAD', I004, MXCOR)
            ENDIF
C
            ETOT = ETOT + EDBL
         ENDIF
C    
 30   CONTINUE
C  
C Notice the multiplication by minus one. The second-order properties
C are the second derivatives of the energy and the first-order perturb
C T amplitudes (T(alpha)) has to be negative by perturbation theory
C arguments. So the real quadratic contribution [[Hbar, T(alpha)],
C T(beta)] is always positive and it is the way it should always be. 
C John stanton, in his original implementation calculates the negative
C of the perturb T(alpha) amplitudes. Because of that the linear 
C terms (T(alpha)*Hbar()), the only terms included in John's original
C implementation comes out as positive values. The quadratic 
C contributions based on John's perturb T amplitudes calculates
C the negative of the actual quadratic contribution. So we need to
C negate the result to have the correct sign for the final quadratic
C contribution. Here, the factor (FACT) is comming from the frequency.
C
      QUAD = QUAD + ONEM*ETOT*FACT
C     
C Now finally, the quadratic contribution for this pair of 
C perturbations. Thank god. Is this the end?.
C
      IF (IFLAGS(1) .GE. 20) THEN
         WRITE(LUOUT, *)
         CALL HEADER('Final quadratic contribution', 1, LUOUT)
         WRITE (LUOUT, *) 'Q(ALPHA, BETA) =  ', QUAD
         WRITE(6,*)
      ENDIF
C     
      RETURN
      END

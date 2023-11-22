      SUBROUTINE GAAAA1(NIRREP,NBAS,IOFFL,
     &                  MULT,NRCSH,KHKTSH,INDGEN,ISTBSH,
     &                  ISHLD,NNIJ,IJADD,NSTSH,NENDSH,
     &                  BUF,BUF1,SORT)
C
C  GAAAA1 PERFORMS THE SORT OF THE GAMMA ELEMENTS OF SYMMETRY
C  TYPE AAAA. THIS IS THE INCORE VERSION.
C 
CEND
C
C  CODED DEC/90/JG
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER AND,OR,XOR
      PARAMETER (IBIT08 = 255)
C
      INTEGER P, Q, R, S, PQ, RS, PQRS
C
      DIMENSION NBAS(8),IOFFL(8)
      DIMENSION MULT(0:7),NRCSH(1),KHKTSH(1),INDGEN(1),ISTBSH(1)
      DIMENSION ISHLD(1),NNIJ(1),IJADD(1),NSTSH(8,1),NENDSH(8,1)
      DIMENSION BUF(1),BUF1(1),SORT(1)
C
C  COMMON /MACHSP/ CONTAINS MACHINE DEPENDENT STUFF
C
      COMMON /MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
C
C STATEMENT FUNCTIONS FOR PACKING, UNPACKING AND ADDRESSING 
C
      IBTAND(I,J) = AND(I,J)
      IBTSHL(I,J) = ISHFT(I,J)
      IBTSHR(I,J) = ISHFT(I,-J)
      IPART(I,J) = IBTAND(IBTSHR(I,8*(J-1)),IBIT08)
      ISHLL(I)   = IBTAND(I,IBIT08)
      ICONTR(I)  = IBTAND(IBTSHR(I,8),IBIT08)
      IANG(I)    = IBTAND(IBTSHR(I,16),IBIT08)
      JRREP(I)   = IBTAND(IBTSHR(I,24),IBIT08)
      ICR(I,J,K,L)  = NRCP*(NRCQ*(NRCR*(L-1)+K-1)+J-1)+I
      IAR(I,J,K,L)  = KHKTP*(KHKTQ*(KHKTR*(L-1)+K-1)+J-1)+I
      IRR(I,J,K,L)  = MULP*(MULQ*(MULR*(L)+K)+J)+I+1
      IADDRS(I,J,K) = NRCTOT*(MULTOT*(K-1)+J-1)+I
C
C  LOOP OVER ALL IRREPS :
C
      DO 100 IRREP=1,NIRREP
C
C  RESET COUNTER KCOUNT TO ZERO
C
       KCOUNT=0
C
C  DETERMINE NUMBER OF BASIS FUNCTIONS IN THIS IRREP
C
       NBAS1=NBAS(IRREP)
C
C  GET OFFSET FOR THIS IRREP IN THE LIST OF BASIS FUNCTIONS
C
       IOFFL1=IOFFL(IRREP)
C
C  LOOP OVER ALL P,Q,R,S IN CANONICAL ORDER
C
C   P GE Q P GE R R GE S AND IF P EQ R Q GE S
C
C   FOR EACH PQ A DISTRIBUTION WITH THE RS ELEMENT IS READ
C
C   LOOP OVER P
C
       DO 100 L1=IOFFL1+1,IOFFL1+NBAS1
C
C  GET PARAMETERS FOR THIS ORBITAL, NOTE THAT THE SHELLS
C  ARE ALSO IN CANONICAL ORDER WHEN THE ORBITALS ARE
C
        LAB1 = INDGEN(L1)
        P = ISHLL(LAB1)
        MULP = MULT(ISTBSH(P))
        NRCP = NRCSH(P)
        KHKTP = KHKTSH(P)
        IP = ICONTR(LAB1)
        KP = IANG(LAB1)
        NP = JRREP(LAB1)
C
C  SET MODIFIED LIMITS FOR L3
C
        NMAX3=NENDSH(IRREP,P)
C
C  LOOP OVER Q
C
       DO 100 L2=IOFFL1+1,L1
C
C  GET PARAMETER FOR THIS ORBITAL
C
        LAB2 = INDGEN(L2)
        Q = ISHLL(LAB2)
        MULQ = MULT(ISTBSH(Q))
        NRCQ = NRCSH(Q)
        KHKTQ = KHKTSH(Q)
        IQ = ICONTR(LAB2)
        KQ = IANG(LAB2)
        NQ = JRREP(LAB2)
C
C  SET MODIFIED POSSIBLE LIMIT FOR S
C
        NMAXQ=NENDSH(IRREP,Q)
C
C  ADDRESS OF SHELL COMBINATION P,Q
C
        PQ = (P*(P-1))/2 + Q
C
C GET ADDRESS WHICH CORRESPONDS TO THE SHELLS PQ ON THE LEFT SIDE
C
        ISHLAD1=0
        IF(PQ.NE.1) THEN
         MAXPQ=PQ-1
         DO 1000 M=1,MAXPQ
c          ISHLAD1=ISHLAD1+NNIJ(M)*IJADD(M)
          ISHLAD1=ISHLAD1+NNIJ(M)*IJADD(M+1)
1000     CONTINUE
        ENDIF
c         ISHLAD1=ISHLD(PQ)
C
C  PICK UP A LIST FROM THE AOGAM VERSION OF MOINTS
C  INCREMENT FIRST KCOUNT
C
        KCOUNT=KCOUNT+1   
C
        CALL GETLST(BUF,KCOUNT,1,1,1,IRREP)
C
C EXPAND THE TRIANGULAR MATRIX R>S FOR FIXED P,Q TO FULL SQUARE
C
        CALL EXPND(BUF,BUF1,NBAS1)
C
C LOOP OVER ALL ELEMENTS IN THIS BUFFER WHICH ARE NEEDED
C 
        DO 200 L3=IOFFL1+1,NMAX3
C
C  GET PARAMETER OF THIRD ORBITAL
C
         LAB3 = INDGEN(L3)
         R = ISHLL(LAB3)
         MULR = MULT(ISTBSH(R))
         NRCR = NRCSH(R)
         KHKTR = KHKTSH(R)
         IR = ICONTR(LAB3)
         KR = IANG(LAB3)
         NR = JRREP(LAB3)
C
C  LOOP OVER FOURTH ORBITALS, MAKE SURE THAT WE HAVE A CANONICAL ORDERING
C
         NMAX4=NENDSH(IRREP,R)
         IF(P.EQ.R) NMAX4=NMAXQ
C
C
         DO 200 L4=IOFFL1+1,NMAX4 
C
C GET ADRESS OF THE REQUIRED GAMMA ELEMENT 
C
          INDEX=(L4-IOFFL1)+(L3-IOFFL1-1)*NBAS1
C
C GET CHARACTERISTICS OF FOURTH ORBITAL
C
          LAB4 = INDGEN(L4)
          S = ISHLL(LAB4)
          MULS = MULT(ISTBSH(S))
          NRCS = NRCSH(S)
          KHKTS = KHKTSH(S)
          IS = ICONTR(LAB4)
          KS = IANG(LAB4)
c          NS = JRREP(LAB4)
C
C  ADDRESS FOR SHELl COMBINATION R,S
C
          RS = (R*(R-1))/2 + S
C
C  FINAL ADDRESS IN FULL LIST OF SHELLS P,Q,R,S
C
          PQRS = (PQ*(PQ-1))/2 + RS
C
C  STARTING ADDRESS OF SHELL COMBINATION P,Q,R,S IN THE LIST OF GAMMAS
C
c          ISHLAD=ISHLAD1
c          IF(RS.NE.1) THEN
c           ISHLAD=ISHLAD+NNIJ(PQ)*IJADD(RS-1)
c          ENDIF
          ISHLAD=ISHLAD1+NNIJ(PQ)*IJADD(RS)
C
C  EXTRACT SHELL COMPONENT INFORMATION FROM LABELS AND
C  COMPUTE SORT ADDRESS AND OFFSET
C
         NRCTOT = ICR(NRCP,NRCQ,NRCR,NRCS)
c         MULTOT = IRR(MULP-1,MULQ-1,MULR-1,MULS-1)
         MULTOT = IRR(MULP-1,MULQ-1,MULR-1,0)
C
C  ILOC IS THE ADDRESS OF THIS ELEMENT IN THE CHAINING AREA
C
         ILOC=ISHLAD+IADDRS(ICR(IP,IQ,IR,IS),
c     1                           IRR(NP,NQ,NR,NS),
     1                           IRR(NP,NQ,NR,0),
     2                           IAR(KP,KQ,KR,KS))
C
         SORT(ILOC)=BUF1(INDEX)
C
C  END OF RS LOOP
C
200      CONTINUE
C 
         IF(P.EQ.Q.AND.L1.NE.L2) THEN
C
C LOOP OVER ALL ELEMENTS IN THIS BUFFER WHICH ARE NEEDED
C 
        DO 300 L3=IOFFL1+1,NMAX3
C
C  GET PARAMETER OF THIRD ORBITAL
C
         LAB3 = INDGEN(L3)
         R = ISHLL(LAB3)
         MULR = MULT(ISTBSH(R))
         NRCR = NRCSH(R)
         KHKTR = KHKTSH(R)
         IR = ICONTR(LAB3)
         KR = IANG(LAB3)
         NR = JRREP(LAB3)
C
C  LOOP OVER FOURTH ORBITALS, MAKE SURE THAT WE HAVE A CANONICAL ORDERING
C
         NMAX4=NENDSH(IRREP,R)
         IF(P.EQ.R) NMAX4=NMAXQ
C
C
         DO 300 L4=IOFFL1+1,NMAX4 
C
C GET ADRESS OF THE REQUIRED GAMMA ELEMENT 
C
          INDEX=(L4-IOFFL1)+(L3-IOFFL1-1)*NBAS1
C
C GET CHARACTERISTICS OF FOURTH ORBITAL
C
          LAB4 = INDGEN(L4)
          S = ISHLL(LAB4)
          MULS = MULT(ISTBSH(S))
          NRCS = NRCSH(S)
          KHKTS = KHKTSH(S)
          IS = ICONTR(LAB4)
          KS = IANG(LAB4)
c          NS = JRREP(LAB4)
C
C  ADDRESS FOR SHELl COMBINATION R,S
C
          RS = (R*(R-1))/2 + S
C
C  FINAL ADDRESS IN FULL LIST OF SHELLS P,Q,R,S
C
          PQRS = (PQ*(PQ-1))/2 + RS
C
C  STARTING ADDRESS OF SHELL COMBINATION P,Q,R,S IN THE LIST OF GAMMAS
C
c          ISHLAD=ISHLAD1
c          IF(RS.NE.1) THEN
c           ISHLAD=ISHLAD+NNIJ(PQ)*IJADD(RS-1)
c          ENDIF
           ISHLAD=ISHLAD1+NNIJ(PQ)*IJADD(RS)
C
C  EXTRACT SHELL COMPONENT INFORMATION FROM LABELS AND
C  COMPUTE SORT ADDRESS AND OFFSET
C
         NRCTOT = ICR(NRCQ,NRCP,NRCR,NRCS)
c         MULTOT = IRR(MULQ-1,MULP-1,MULR-1,MULS-1)
         MULTOT = IRR(MULQ-1,MULP-1,MULR-1,0)
C
C  ILOC IS THE ADDRESS OF THIS ELEMENT IN THE CHAINING AREA
C
         ILOC=ISHLAD+IADDRS(ICR(IQ,IP,IR,IS),
c     1                      IRR(NQ,NP,NR,NS),
     1                      IRR(NQ,NP,NR,0),
     2                      IAR(KQ,KP,KR,KS))
C
          SORT(ILOC)=BUF1(INDEX)
C
300   CONTINUE
C
      ENDIF
C
C END OF ALL LOOPS
C
100   CONTINUE
C
      RETURN
      END

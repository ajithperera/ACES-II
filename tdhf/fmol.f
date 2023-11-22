      SUBROUTINE FMOL(D,F,BUF,IBUF,NINTMX,IA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*8 LABEL
      DIMENSION LABEL(52)
      DIMENSION D(2),F(2),IA(2)
      DIMENSION BUF(NINTMX),IBUF(NINTMX)
      COMMON/SWPPP/INP,IAMO,IFAMO,IINDO,IORTH                
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/INFOA/NBASIS,NUMSCF,NX,NMO2,NOC,NVT,NVO     
      DATA ZERO,ONE,TWO/0.D0,1.D0,2.D0/ 
C     WRITE(6,*) ' Now in FMOK '
      NTAPE=10
C  ...... F takes constant part from lower order .
      DO 1  M=2,NUMSCF
      MAX=M-1                                    
      DO 1  N=1,MAX                            
      NIJ=IA(M)+N                              
    1 F(NIJ)=F(NIJ)*TWO    
C     IF(IFAMO.EQ.0) THEN                   
      OPEN(UNIT=NTAPE,FORM='UNFORMATTED',FILE='HF2',ACCESS='SEQUENTIAL')
      READ(NTAPE) LABEL(1)
      READ(NTAPE) LABEL(1)
      READ(NTAPE) LABEL(1)
      READ(NTAPE) LABEL(1)
C     ELSE
C     OPEN(UNIT=NTAPE,FORM='UNFORMATTED',FILE='IIII',
C    X ACCESS='SEQUENTIAL')
C     CALL SEARCH('TWOELSUP',NTAPE)
C     END IF
C     WRITE(6,*) LABEL
C     WRITE(6,*) ' # , INTEGRAL I , J , K , L '
      NTOT=0
    5 READ(NTAPE,END=50) BUF,IBUF,NUT
      IF(NUT.LE.0) GO TO 30
      DO 10 IX=1,NUT
      NTOT = NTOT + 1
C     I = AND(ISHFT(IBUF(IX),-24),255)
C     J = AND(ISHFT(IBUF(IX),-16),255)
C     K = AND(ISHFT(IBUF(IX),- 8),255)
C     L = AND(IBUF(IX),255)
C     I = (SHIFTR(IBUF(IX),48).AND.65535)
C     J = (SHIFTR(IBUF(IX),32).AND.65535)
C     K = (SHIFTR(IBUF(IX),16).AND.65535)
C     L = IBUF(IX).AND.65535
      L = AND(IBUF(IX),IALONE)
      K = AND(ISHFT(IBUF(IX),-IBITWD),IALONE)
      J = AND(ISHFT(IBUF(IX),-2*IBITWD),IALONE)
      I = AND(ISHFT(IBUF(IX),-3*IBITWD),IALONE)
      IF(J.GT.I) THEN
      IJE=I
      I=J
      J=IJE
      END IF
      IF(L.GT.K) THEN
      KLE=K
      K=L
      L=KLE
      END IF
C  I>J, K>L and (I,J)>(K,L) 
      VAL=BUF(IX)
      IF(I.EQ.J) VAL= VAL/2.D0
      IF(K.EQ.L) VAL= VAL/2.D0
      NIJ=IA(I)+J
      NKL=IA(K)+L
      IF(NKL.GT.NIJ) THEN
      IKE=I
      I=K
      K=IKE
      JLE=J
      J=L
      L=JLE
      END IF
      IF(NIJ.EQ.NKL) VAL= VAL/2.D0
C     write(6,*) I,J,K,L,VAL
      VAL4 =(VAL+VAL)+(VAL+VAL)
      NIK=IA(I)+K             
      NIL=IA(I)+L           
      SIGNJK= ONE
      SIGNJL = ONE 
      IF(J.LT.K) GO TO 15  
      NJK=IA(J)+K          
      NJL=IA(J)+L         
      GO TO 25          
   15  NJK=IA(K)+J
      SIGNJK= -ONE       
      IF(J.LT.L) GO TO 20  
      NJL=IA(J)+L          
      GO TO 25           
   20 NJL=IA(L)+J 
      SIGNJL=-ONE       
   25 F(NIK)=F(NIK)-VAL*D(NJL)*SIGNJL   
      F(NIL)=F(NIL)-VAL*D(NJK)*SIGNJK  
      F(NJK)=F(NJK)-VAL*D(NIL)*SIGNJK 
      F(NJL)=F(NJL)-VAL*D(NIK)*SIGNJL      
C     IF(NIJ.LT.5.AND.NKL.LT.5)WRITE(6,*) BUF(IX),I,J,K,L
   10 CONTINUE
      GO TO  5
   30 CONTINUE
      IF(NUMSCF.EQ.1) GO TO 50
C .... when J=K or J=L SINGNJK and SIGNJL is no proper for F(NJK) and F(NJL)
C ... they are OK for F(NIK) and F(NIL) because D(NJL) and D(NJK) = 0
C .. therefor F(IJ)=0 must be forced for I=J for anti-hermitian
      DO 40  M=1,NUMSCF
      DO 40  N=1,M
      IJ=IA(M)+N 
      IF(M.EQ.N) THEN
      F(IJ)=ZERO
      ELSE
      F(IJ)=F(IJ)/TWO    
      END IF
   40 CONTINUE 
   50 CONTINUE       
      CLOSE(NTAPE) 
      RETURN 
      END

      SUBROUTINE ONEPRS(NDPROP,IVECT,X,BUFFER,INDX,NT,ISYMO,EVEC
     X ,HMO,DAP,SCR,LEN2,NSIZ1,NUMSCF)
C  .......... here NSIZ1=NBASIS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*8 NLAB,DABEL,LABEL(4),XYZ(12)
      COMMON/IPRNT/ IOPEV,IOPDA,IOPU,IOPFE,IOPPR,IWRPA
      COMMON/SWPPP/INP,IAMO,IFAMO,IINDO,IORTH
      COMMON/INFSYM/NSYMHF,NSO(8),NOCS(8),NVTS(8),IDPR(8,8),NVOS(8)
     X ,NIJS(8),NIJS2(8),NRDS(8),NIJSR(8)
      COMMON/INFPRS/IPRSYM(12),NPRSYM(8),JPRSYM(8,12)
      DIMENSION X(NSIZ1,NSIZ1,3),PTOT(3),HMO(IVECT,3),
     X EVEC(NSIZ1,NSIZ1),DAP(1),SCR(1),ISYMO(NUMSCF)
      DIMENSION BUFFER(LEN2),INDX(LEN2)
      DATA NLAB/'********'/
      DATA DABEL/'   DEN  '/
      DATA XYZ/'     x  ','     y  ','     z  ','    xx  ','    yy  ',
     X '    zz  ','    xy  ','    yz  ','    zx  ','    yx  ','    zy  '
     X ,'    xz  '/
      DATA ZERO/0.D00/,ONE/1.D0/
      DATA TWO/2.D00/
      DATA THRPR/1.D-7/
      OPEN(UNIT=NDPROP,FILE='VPOUT',STATUS='UNKNOWN',FORM='UNFORMATTED')
      NT=1
   10 READ (NDPROP,END=99,ERR=99) LABEL,PTOT(NT)
      IF (LABEL(1).NE.NLAB) GO TO 10
      WRITE(6,*) LABEL(4),NT,PTOT(NT)
      DO 20 I=1,IVECT
   20 HMO(I,NT)=ZERO
   30 READ (NDPROP,END=99,ERR=99) BUFFER,INDX,LIMIT
C     write(6,*) ' limit = ',LIMIT
C     write(6,1000) (BUFFER(II),INDX(II),II=1,LIMIT)
C1000 format(1h ,6(f10.4,I3))
      IF (LIMIT.LE.0) GO TO 40
      DO 50 I=1,LIMIT
      K=INDX(I)
   50 HMO(K,NT)=BUFFER(I)
      GO TO 30
   40 CONTINUE
      CALL RTPUT(X(1,1,NT),HMO(1,NT),NSIZ1)
      IF(IOPPR.NE.0) CALL OUTMXD(X(1,1,NT),NSIZ1,NSIZ1,NSIZ1)
      NT=NT+1
      IF(NT.GT.3) GO TO 99
      GO TO 10
   99 NT = NT - 1
  110 READ (NDPROP,END=199,ERR=199) LABEL
      IF (LABEL(1).NE.NLAB) GO TO 110
C  ...... in the case of OAO basis ............................
      IF(IORTH.NE.0) THEN
      WRITE(6,*) 'label =',LABEL(4)
      DO 120 I=1,IVECT
  120 HMO(I,3)=ZERO
  130 READ (NDPROP,END=199,ERR=199) BUFFER,INDX,LIMIT
      IF (LIMIT.LE.0) GO TO 140
      DO 150 I=1,LIMIT
      K=INDX(I)
  150 HMO(K,3)=BUFFER(I)
      GO TO 130
  140 CONTINUE
      CALL TRPUT(X(1,1,1),HMO(1,1),NSIZ1)
      CALL TRPUT(X(1,1,2),HMO(1,2),NSIZ1)
      CALL RTPUT(X(1,1,1),HMO(1,3),NSIZ1)
      CALL TRPUT(X(1,1,3),HMO(1,3),NSIZ1)
      IF(IOPPR.NE.0) THEN
      WRITE(6,*) ' Overlap Matrix '
      CALL OUTMXD(X,NSIZ1,NSIZ1,NSIZ1)
      END IF
C     CALL MHOUSE(NSIZ1,NSIZ1,X,X(1,1,2))
      CALL MHOUSE(NSIZ1,NSIZ1,NSIZ1+1,X,X(1,1,2),SCR)
CSS      WRITE(6,*) ' Eigen values of S-matrix '
CSS      WRITE(6,*) (X(J,J,1),J=1,NSIZ1)
      DO 159 I=1,NSIZ1
      IF(X(I,I,1).LE.0.D0) WRITE(6,*) ' Minus Eigen value !!!!!! '
  159 X(I,I,3)=SQRT(X(I,I,1))
C     write(6,*)  ' sqrt(eigen value)'
C     WRITE(6,*) (X(J,J,3),J=1,NSIZ1)
      DO 160 J=1,NSIZ1
      DO 160 I=1,NSIZ1
C     write(6,*) I,J,X(I,J,2),X(J,J,3)
  160 X(I,J,1)=X(I,J,2)/X(J,J,3)
C     write(6,*) ' before TRANSQ '
      CALL TRANSQ(X(1,1,2),NSIZ1)
C     write(6,*) ' X(1,1,1)'
C     CALL OUTMXD(X(1,1,1),NSIZ1,NSIZ1,NSIZ1)
C     write(6,*) ' X(1,1,2)'
C     CALL OUTMXD(X(1,1,2),NSIZ1,NSIZ1,NSIZ1)
C     write(6,*) ' X(1,1,3)'
C     CALL OUTMXD(X(1,1,3),NSIZ1,NSIZ1,NSIZ1)
C     write(6,*) ' before MATMUL NSIZ1',NSIZ1
      CALL MATMUL(X,X(1,1,2),X(1,1,3),NSIZ1,NSIZ1,NSIZ1,1,0)
C     write(6,*) 'after MATMUL '
      IF(IOPPR.NE.0) THEN
      WRITE(6,*) ' S**(-1/2) '
      CALL OUTMXD(X(1,1,3),NSIZ1,NSIZ1,NSIZ1)
      END IF
C ............................................................
      DO 161 M=1,3
      CALL RTPUT(X,HMO(1,M),NSIZ1)
      CALL MATMUL(X,X(1,1,3),X(1,1,2),NSIZ1,NSIZ1,NSIZ1,1,0)
      CALL MATMUL(X(1,1,3),X(1,1,2),X,NSIZ1,NSIZ1,NSIZ1,1,0)
      CALL TRPUT(X,HMO(1,M),NSIZ1)
  161 CONTINUE
      END IF
  199 CONTINUE
      IF(IAMO.EQ.0) THEN
C  External perturbation over MO basis
      DO 19 M=1,3
      WRITE(6,*) ' Properties over MO ,Component = ',M
      CALL RTPUT(X,HMO(1,M),NSIZ1)
C     CALL MATMUL(X,EVEC,X(1,1,2),NSIZ1,NUMSCF,NSIZ1,1,0)
      CALL XGEMM('N','N',NSIZ1,NUMSCF,NSIZ1,ONE,X,NSIZ1,EVEC,NSIZ1
     X ,ZERO,X(1,1,2),NSIZ1)
      CALL TRANSQ(EVEC,NSIZ1)
C     CALL MATMUL(EVEC,X(1,1,2),X,NUMSCF,NUMSCF,NSIZ1,1,0)
      CALL XGEMM('N','N',NUMSCF,NUMSCF,NSIZ1,ONE,EVEC,NSIZ1,X(1,1,2)
     X ,NSIZ1,ZERO,X,NSIZ1)
      CALL TRANSQ(EVEC,NSIZ1)
      CALL TRPUT(X,HMO(1,M),NSIZ1)
C  ......  Properties are packed in HMO trianglly up to NUMSCF
C  ......  Space over NUMSCF in HMO is garbage ....................
      IF(IOPPR.NE.0) CALL OUTMXD(X,NSIZ1,NSIZ1,NSIZ1)
   19 CONTINUE
      END IF
C  ..............  Symmetry Informations of properties .......
C  ..............   x,y,z ....................................
      NISO=NSO(1)
      DO 46 M=1,3
      CALL RTPUT(X,HMO(1,M),NSIZ1)
C     CALL OUTMXD(X,NSIZ1,NSIZ1,NSIZ1)
      IPRSYM(M)=0
C     NJSUM=0
C     DO 41 ISS=1,NSYMHF
C     NJSO=NSO(ISS)
C     DO 45 II=1,NISO
C     DO 45 JJ=1,NJSO
C     IF(ABS(X(II,NJSUM+JJ,M)).GT.THRPR) IPRSYM(M)=ISS
C  45 CONTINUE
C     NJSUM=NJSUM+NJSO
C     write(6,*) ' niso,njso,iss,njsum ',NISO,NJSO,ISS,NJSUM
C  41 CONTINUE
      DO 41 II=1,NUMSCF
      DO 41 JJ=1,NUMSCF
C     IF(ABS(X(II,JJ,1)).GT.THRPR) WRITE(6,*) M,II,JJ,ISYMO(II),
C    X ISYMO(JJ),IDPR(ISYMO(II),ISYMO(JJ))
      IF(ABS(X(II,JJ,1)).GT.THRPR) IPRSYM(M)=IDPR(ISYMO(II),ISYMO(JJ))
   41 CONTINUE
      write(6,*) ' m,iprsym = ',M,IPRSYM(M)
   46 CONTINUE
C  .... Symmetry information on xx,yy,zz,xy,yz,zx,yx,zy,xz
      IPRSYM(4)=IDPR(IPRSYM(1),IPRSYM(1))
      IPRSYM(5)=IDPR(IPRSYM(2),IPRSYM(2))
      IPRSYM(6)=IDPR(IPRSYM(3),IPRSYM(3))
      IPRSYM(7)=IDPR(IPRSYM(1),IPRSYM(2))
      IPRSYM(8)=IDPR(IPRSYM(2),IPRSYM(3))
      IPRSYM(9)=IDPR(IPRSYM(3),IPRSYM(1))
      IPRSYM(10)=IDPR(IPRSYM(2),IPRSYM(1))
      IPRSYM(11)=IDPR(IPRSYM(3),IPRSYM(2))
      IPRSYM(12)=IDPR(IPRSYM(1),IPRSYM(3))
      WRITE(6,*) ' Symmetry Information '
      DO 43 II=1,8
   43 NPRSYM(II)=0
      DO 42 III=1,12
      IIIS=IPRSYM(III)
      NPRSYM(IIIS)=NPRSYM(IIIS)+1
      JPRSYM(IIIS,NPRSYM(IIIS))=III
   42 WRITE(6,*) XYZ(III),' = ',IPRSYM(III)
      DO 44 II=1,8
      NNPRS=NPRSYM(II)
   44 WRITE(6,*) 'Symmetry ',II,' includes ',
     X  (XYZ(JPRSYM(II,JJ)),JJ=1,NNPRS)
      WRITE(6,*) ' FIRST-ORDER ENERGY '
      DO 152 MT=1,3
      EHF1=TRACEP(DAP,HMO(1,MT),NUMSCF)
      EHFT1=EHF1 + PTOT(MT)
      WRITE(6,*) MT,'-TH electronic component ',EHF1,
     X ' Total ',EHFT1
  152 CONTINUE
      REWIND NDPROP
      RETURN
      END

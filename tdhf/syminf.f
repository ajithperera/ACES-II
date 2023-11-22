C                                                         
      SUBROUTINE SYM_INF(NIREP,NSIZVO,NSIZ1,NOC,ISYMO,IA,IAS,IVRT,
     X IOCC,IVRTS,IOCCS,IVO,IVOS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     CHARACTER*8 LABEL
      DIMENSION TITLE(24)
      DIMENSION IA(NSIZ1),IAS(NSIZ1,NIREP),IVRT(NSIZVO),IOCC(NSIZVO),
     X IVRTS(NSIZVO,NIREP),IOCCS(NSIZVO,NIREP),IVO(NSIZ1,NOC),
     X IVOS(NSIZ1,NOC,NIREP),ISYMO(NSIZ1)
      COMMON /INFSYM/NSYMHF,NSO(8),NOCS(8),NVTS(8),IDPR(8,8),NVOS(8)
     X ,NIJS(8),NIJS2(8),NRDS(8),NIJSR(8)
C     NTAP=8
C     OPEN (UNIT=NTAP,FILE='IIII',STATUS='UNKNOWN',FORM='UNFORMATTED')
C     READ(NTAP) TITLE,NSYMHF,(NSO(I),I= 1,NSYMHF),POTNUC
      WRITE(6,11)NSYMHF,(NSO(I),I = 1,NSYMHF)
   11 FORMAT(/,10X,' The number of symmetry types is ',I3,
     *       /,10X,' The number of each sym. type are',
     *       /,10X,8I6)
      WRITE(6,*) 'Group multiplication table '
      DO 45 I=1,8
   45 WRITE(6,*) (IDPR(I,J),J=1,8)
      NSTO=0
      IVECT=0
      DO 12 I=1,NSYMHF
      IVECT=IVECT+((NSO(I)+1)*NSO(I))/2
      NNSO=NSO(I)
   12 NSTO=NSTO+NSO(I)
      WRITE(6,13)NSTO,IVECT
   13 FORMAT(/,10X,' The number of MO is ',I3,
     *       /,10X,' The number of onel elements is ',I5)
      NVT=NSTO-NOC
      WRITE(6,*) ' NVT,NOC ',NVT,NOC
      DO 46 I=1,NSTO
   46 ISYMO(I)=0
      IJ=0
      DO 47 I=1,NSYMHF
      NNOC=NOCS(I)
      IF(NNOC.GT.0) THEN
      DO 48 J=1,NNOC
      IJ=IJ+1
   48 ISYMO(IJ)=I
      END IF
   47 CONTINUE
      DO 49 I=1,NSYMHF
      NNVT=NSO(I)-NOCS(I)
      IF(NNVT.GT.0) THEN
      DO 51 J=1,NNVT
      IJ=IJ+1
   51 ISYMO(IJ)=I
      END IF
   49 CONTINUE
      WRITE(6,*) ' ISYMO ',(ISYMO(II),II=1,NSTO)
      NIJ=(NSTO*(NSTO+1))/2
      NSTO2=NSTO*NSTO
      DO 21 I=1,NSTO
      IA(I)=(I*(I-1))/2
   21 CONTINUE
C     WRITE(6,*) ' IA = ',(IA(II),II=1,NSTO)
      DO 50 I=1,NSTO
      DO 50 J=1,NOC
   50 IVO(I,J)=0
      DO 25 I=1,NSYMHF
      DO 25 J=1,NSTO
      DO 25 K=1,NOC
      IVOS(J,K,I)=0
   25 NVOS(I)=0
      IJ=0
      DO 22 I=1,NOC
      DO 22 J=1,NVT
      IS= ISYMO(I)
      J1=J+NOC
      JS = ISYMO(J1)
      IJSYM= IDPR(IS,JS)
      NVOS(IJSYM)=NVOS(IJSYM) + 1
      IJ=IJ+1
      IVRT(IJ)=J1
      IOCC(IJ)=I
      IVRTS(NVOS(IJSYM),IJSYM)=J1
      IOCCS(NVOS(IJSYM),IJSYM)=I
      JI=(I-1)*NSTO+J1
      IVOS(J1,I,IJSYM)=NVOS(IJSYM)
   22 IVO(J1,I)=IJ
      NIJVO=IJ
C     WRITE(6,*) ' IJ,JI ',IJ,JI
      WRITE(6,*) ' V-O ',NIJVO,(NVOS(II),II=1,NSYMHF)
C     WRITE(6,*) ' IVRT,IOCC ',(IVRT(II),IOCC(II),II=1,IJ)
C     WRITE(6,*) ' IVO '
C     DO 52 I=1,NSTO
C  52 WRITE(6,*) (IVO(I,J),J=1,NOC)
      DO 27 I=1,NSYMHF
C     WRITE(6,*) I,' IVOS '
C     DO 53 II=1,NSTO
C  53 WRITE(6,*) (IVOS(II,J,I),J=1,NOC) 
      IF(I.EQ.1) THEN 
      NIJS(I) = 0 
      NIJS2(I) = 0 
      NIJSR(I) = 0 
      ELSE 
      NIJS(I) = NIJS(I-1) + NVOS(I-1)
      NIJS2(I) = NIJS2(I-1) + NVOS(I-1)*NVOS(I-1)
      NIJSR(I) = NIJSR(I-1) + NVOS(I-1)*NRDS(I-1)
      END IF
   27 CONTINUE
      WRITE(6,*) ' NVOS,NIJS ',(NVOS(II),NIJS(II),II=1,NSYMHF)
      WRITE(6,*) ' NIJS2,NIJSR ',(NIJS2(II),NIJSR(II),II=1,NSYMHF)
      WRITE(6,*) ' IVRTS,IOCCS '
      DO 26 I=1,NSYMHF
      NIJ=NVOS(I)
C     WRITE(6,*) I,' sym ', (IVRTS(II,I),IOCCS(II,I),II=1,NIJ)
   26 CONTINUE
      DO 30 ISYM=1,NSYMHF
      NSTOS=NSO(ISYM)
      IF(NSTOS.GT.0) THEN
      DO 31 I=1,NSTOS
      IAS(I,ISYM)=(I*(I-1))/2
   31 CONTINUE
C     WRITE(6,*) ' IAS = ',NSTOS,(IAS(II,ISYM),II=1,NSTOS)
      END IF
   30 CONTINUE
C     CLOSE(NTAP)
      RETURN
      END

      SUBROUTINE WRDIR_C4(IUNIT,IREC,IVEC,LENGTH,IMOD)
      IMPLICIT INTEGER (A-Z)
      DIMENSION IVEC(LENGTH)
      COMMON/FSTRIPE/NUMSPLIT(5)
      IMOD=0
      nsplit=max(1,numsplit(iunit-49))
      IF(IREC.LE.0)RETURN
      ioff=mod(irec-1,nsplit)
      irealunit=iunit+ioff*5
      irealrec =1+(irec-1)/nsplit
      WRITE(IrealUNIT,REC=IrealREC,ERR=555,IOSTAT=IER)IVEC
      RETURN
555   CALL IOERR('WRDIR',IUNIT,IER)
      END

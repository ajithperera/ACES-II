      SUBROUTINE RDDIR_c4(IUNIT,IREC,IVEC,LENGTH,IMOD)
      IMPLICIT INTEGER (A-Z)
      DIMENSION IVEC(LENGTH)
      COMMON/FSTRIPE/NUMSPLIT(5)
      nsplit=max(1,numsplit(iunit-49))
      ioff=mod(irec-1,nsplit)
      irealunit=iunit+ioff*5
      irealrec =1+(irec-1)/nsplit
      READ(IrealUNIT,REC=IrealREC,ERR=555,IOSTAT=IER)IVEC
      IMOD=0
      RETURN
555   CALL IOERR('RDDIR',IUNIT,IER)
      END

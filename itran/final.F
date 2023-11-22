      SUBROUTINE FINAL(LUHF2,BUF,IBUF,ILNBUF,INUT)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION BUF(ILNBUF),IBUF(ILNBUF)
      
      WRITE(LUHF2) BUF,IBUF,INUT
C
#ifdef _DEBUG_LVLM1
      Write(6,*) 
      Write(6,*) "Writing in final", inut
      Write(*,'(4(1X,F12.7)))'), (BUF(I), I=1, ilnbuf)
#endif 
      RETURN
      END

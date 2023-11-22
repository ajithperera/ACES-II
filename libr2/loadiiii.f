      SUBROUTINE LOADIIII(W,BUF,IBUF,NAO,ILNBUF,LUINT,NSTART,
     &                    NEND)
C
C THIS ROUTINE LOADS THE AO INTEGRALS FROM THE IIII FILE. 
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER AND,OR
      DIMENSION BUF(ILNBUF),IBUF(ILNBUF),W(NAO*NAO,*)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      IUPKI(INT)=AND(INT,IALONE)
      IUPKJ(INT)=AND(ISHFT(INT,-IBITWD),IALONE)
      IUPKK(INT)=AND(ISHFT(INT,-2*IBITWD),IALONE)
      IUPKL(INT)=AND(ISHFT(INT,-3*IBITWD),IALONE)
      IPACK(I,J,K,L)=OR(OR(OR(I,ISHFT(J,IBITWD)),ISHFT(K,2*IBITWD)),
     &                  ISHFT(L,3*IBITWD))
      INDX(I,J)=J+(I*(I-1))/2
      INDX3(I,J)=I+(J*(J-1))/2
      INDX2(I,J,N)=I+(J-1)*N
C
      NAOBUF=0
      NUMINT=0
      NUMDIS=(NEND-NSTART+1)
      CALL LOCATE(LUINT,'TWOELSUP')
      CALL ZERO(W,NAO*NAO*NUMDIS)
1     READ(LUINT)BUF,IBUF,NUT
      NAOBUF=NAOBUF+1      
      DO 10 INT=1,NUT
       X=BUF(INT)
       I=IUPKI(IBUF(INT))
       J=IUPKJ(IBUF(INT))
       K=IUPKK(IBUF(INT))
       L=IUPKL(IBUF(INT))
       IK=INDX2(I,K,NAO)
       IL=INDX2(I,L,NAO)
       JK=INDX2(J,K,NAO)
       JL=INDX2(J,L,NAO)
       KI=INDX2(K,I,NAO)
       KJ=INDX2(K,J,NAO)
       LI=INDX2(L,I,NAO)
       LJ=INDX2(L,J,NAO)
       if(jl-nstart.ge.0.and.jl-nend.le.0)W(IK,JL-nstart+1)=X
       if(il-nstart.ge.0.and.il-nend.le.0)W(JK,IL-nstart+1)=X
       if(jk-nstart.ge.0.and.jk-nend.le.0)W(IL,JK-nstart+1)=X
       if(ik-nstart.ge.0.and.ik-nend.le.0)W(JL,IK-nstart+1)=X

       if(lj-nstart.ge.0.and.lj-nend.le.0)W(KI,LJ-nstart+1)=X
       if(li-nstart.ge.0.and.li-nend.le.0)W(KJ,LI-nstart+1)=X
       if(kj-nstart.ge.0.and.kj-nend.le.0)W(LI,KJ-nstart+1)=X
       if(ki-nstart.ge.0.and.ki-nend.le.0)W(LJ,KI-nstart+1)=X
10    CONTINUE
      NUMINT=NUMINT+NUT
C
      IF(NUT.NE.-1)GOTO 1
C
c      WRITE(6,*)' processed ',numint,' ao basis integrals ',
c     &          'from ',naobuf,' buffers.'
C
      RETURN
      END

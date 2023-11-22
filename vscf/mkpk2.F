      SUBROUTINE MKPK2(PK,BUF,IBUF,ILNBUF,IPKSIZ,NBAS,IUHF)
C
C  This routine reads in the AO integrals and places them in the integral
C  vector W based on what kind of symmetry integral it is, i.e., (aa|aa),
C  (ab|ab), (aa|bb), or (ab|cd).  
C
C      PK - PK integral vector; size is IPKSIZ
C
C     BUF - AO integral buffer; dimension is ILNBUF
C
C    IBUF - AO integral index buffer; dimension is ILNBUF
C
C  Other variables used, such as offsets and the like, are contained in
C  the common block SYMM2.  Most of these are determined in the
C  subroutine SYMSIZ, and their definitions can be found in the 
C  description of SYMSIZ.
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER(LUINT=10, LUNITSE=25)
C
      INTEGER AND,OR
      CHARACTER*80 FNAME
      DIMENSION PK(IPKSIZ),BUF(ILNBUF),IBUF(ILNBUF)
      DIMENSION IPKOFF(73)
C
#include "molcas.com"
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /FLAGS/ IFLAGS(100)
#include "symm2.com"
      COMMON /PKOFF/ IPKOFF
      COMMON /TIMEINFO/ TIMEIN, TIMENOW, TIMETOT, TIMENEW
C
C  Statement functions for packing and unpacking indices.
C
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
      CALL TIMER(1)
      CALL ZERO(PK,IPKSIZ)
      NAOBUF=0
      NAOINT=0
C
      IF(IUHF.EQ.1) THEN
        SCALE=0.5
      ELSE
        SCALE=0.25
      ENDIF
C
C  Read in a buffer full of integrals.  Make sure they have i>=j, k>=l and
C  ij >= kl.  Then place them in the appropriate location in the integral
C  vector W.
C
      CALL GFNAME('IIII    ',FNAME,ILENGTH)
      OPEN(LUINT,FILE=FNAME(1:ILENGTH),STATUS='OLD',FORM='UNFORMATTED',
     &     ACCESS='SEQUENTIAL')
      CALL LOCATE(LUINT,'TWOELSUP')
    1 READ(LUINT)BUF,IBUF,NUT
      NAOBUF=NAOBUF+1
      DO 10 INT=1,NUT
        IX=IUPKI(IBUF(INT))
        JX=IUPKJ(IBUF(INT))
        KX=IUPKK(IBUF(INT))
        LX=IUPKL(IBUF(INT))
        X=BUF(INT)
        I=MAX(IX,JX)
        J=MIN(IX,JX)
        K=MAX(KX,LX)
        L=MIN(KX,LX)
        IND1=INDX(I,J)
        IND2=INDX(K,L)
        IF(IND1.LT.IND2) THEN
          IX=I
          JX=J
          I=K
          J=L
          K=IX
          L=JX
        ENDIF
C
        scale1=1.
        scale2=1.
        IJUNK=INDX(IDXVEC(I),IDXVEC(J))
        IND1=INDX(INEWVC(I),INEWVC(J))
        IND2=INDX(INEWVC(K),INEWVC(L))
C
C  Add the J piece into PK.
C
        PK(IPKOFF(IJUNK)+INDX2(IND1,IND2,IRPDS1(IJUNK)))=
     &    PK(IPKOFF(IJUNK)+INDX2(IND1,IND2,IRPDS1(IJUNK)))+X
C
C  If this is a UHF calculation, put the J integral in the proper
C  place.
C
        IF(IUHF.NE.0) THEN
          PK(IPKOFF(IJUNK+36)+INDX2(IND1,IND2,IRPDS1(IJUNK)))=
     &      PK(IPKOFF(IJUNK+36)+INDX2(IND1,IND2,IRPDS1(IJUNK)))+X
        ENDIF
C
C  Now add the K piece into PK.
C
        IND3=INDX(INEWVC(I),INEWVC(K))
        IND4=INDX(MAX(INEWVC(J),INEWVC(L)),MIN(INEWVC(J),
     &            INEWVC(L)))
        if(ind4.gt.ind3) then
          icrap=ind3
          ind3=ind4
          ind4=icrap
        endif
        IND5=INDX(INEWVC(I),INEWVC(L))
        IND6=INDX(MAX(INEWVC(J),INEWVC(K)),MIN(INEWVC(J),
     &            INEWVC(K)))
        if(ind6.gt.ind5) then
          icrap=ind5
          ind5=ind6
          ind6=icrap
        endif
        if((i.eq.k.or.j.eq.l).and.(i.ne.j)) scale1=2.
        if((i.eq.j.or.k.eq.l).and.(i.ne.k)) scale2=0.
        if((j.eq.k).and.(i.ne.j.and.j.ne.l)) scale2=2.
        PK(IPKOFF(IJUNK)+INDX2(IND3,IND4,IRPDS1(IJUNK)))=
     &    PK(IPKOFF(IJUNK)+INDX2(IND3,IND4,IRPDS1(IJUNK)))-
     &    scale1*SCALE*X
        PK(IPKOFF(IJUNK)+INDX2(IND5,IND6,IRPDS1(IJUNK)))=
     &    PK(IPKOFF(IJUNK)+INDX2(IND5,IND6,IRPDS1(IJUNK)))-
     &    scale2*SCALE*X
C
   10 CONTINUE
      IF(NUT.EQ.ILNBUF) GOTO 1
C
      NAOINT=(NAOBUF-1)*ILNBUF+NUT
      CLOSE(LUINT,STATUS='KEEP')
C
      IF(NIRREP.EQ.1) GOTO 200
C
      CALL GFNAME('IIJJ    ',FNAME,ILENGTH)
      OPEN(LUINT,FILE=FNAME(1:ILENGTH),STATUS='OLD',FORM='UNFORMATTED',
     &     ACCESS='SEQUENTIAL')
      CALL LOCATE(LUINT,'TWOELSUP')
      NAOBUF=0
    2 READ(LUINT)BUF,IBUF,NUT
      NAOBUF=NAOBUF+1
      DO 20 INT=1,NUT
        IX=IUPKI(IBUF(INT))
        JX=IUPKJ(IBUF(INT))
        KX=IUPKK(IBUF(INT))
        LX=IUPKL(IBUF(INT))
        X=BUF(INT)
        I=MAX(IX,JX)
        J=MIN(IX,JX)
        K=MAX(KX,LX)
        L=MIN(KX,LX)
        IND1=INDX(I,J)
        IND2=INDX(K,L)
        IF(IND1.LT.IND2) THEN
          IX=I
          JX=J
          I=K
          J=L
          K=IX
          L=JX
        ENDIF
C
C  The following conditional handles the (aa|bb) symmetry integrals.
C
        JJUNK=INDX(IDXVEC(I)-1,IDXVEC(K))
        IJUNK=INDX(IDXVEC(I),IDXVEC(K))
        IND1=INDX(INEWVC(I),INEWVC(J))
        IND2=INDX(INEWVC(K),INEWVC(L))
        PK(IPKOFF(IJUNK)+INDX2(IND1,IND2,IRPDS2(2*JJUNK-1)))=
     &   PK(IPKOFF(IJUNK)+INDX2(IND1,IND2,IRPDS2(2*JJUNK-1)))+X
C
C  If this is a UHF calculation, then put the J integral in the proper
C  place.
C
        IF(IUHF.NE.0) THEN
          PK(IPKOFF(IJUNK+36)+INDX2(IND1,IND2,IRPDS2(2*JJUNK-1)))=
     &       PK(IPKOFF(IJUNK+36)+INDX2(IND1,IND2,
     &                                 IRPDS2(2*JJUNK-1)))+X
        ENDIF
   20 CONTINUE
C
      IF(NUT.EQ.ILNBUF) GOTO 2
C
      NAOINT=NAOINT+(NAOBUF-1)*ILNBUF+NUT
      CLOSE(LUINT,STATUS='KEEP')
C
C  The following conditional handles the (ab|ab) symmetry integrals.  They
C  are stored in square canonical order within W, but the indexing within
C  a distribution is rectangular (IND1 & IND2).  The starting point in
C  W is given by IWOFF1.  The value of IJUNK is as for the (aa|aa) types.
C
      CALL GFNAME('IJIJ    ',FNAME,ILENGTH)
      OPEN(LUINT,FILE=FNAME(1:ILENGTH),STATUS='OLD',FORM='UNFORMATTED',
     &     ACCESS='SEQUENTIAL')
      CALL LOCATE(LUINT,'TWOELSUP')
      NAOBUF=0
    3 READ(LUINT)BUF,IBUF,NUT
      NAOBUF=NAOBUF+1
      DO 30 INT=1,NUT
        IX=IUPKI(IBUF(INT))
        JX=IUPKJ(IBUF(INT))
        KX=IUPKK(IBUF(INT))
        LX=IUPKL(IBUF(INT))
        X=BUF(INT)
        I=MAX(IX,JX)
        J=MIN(IX,JX)
        K=MAX(KX,LX)
        L=MIN(KX,LX)
        IND1=INDX(I,J)
        IND2=INDX(K,L)
        IF(IND1.LT.IND2) THEN
          IX=I
          JX=J
          I=K
          J=L
          K=IX
          L=JX
        ENDIF
C
        scale1=1.0
        JJUNK=INDX(IDXVEC(I)-1,IDXVEC(J))
        IJUNK=INDX(IDXVEC(I),IDXVEC(J))
        ind1=indx(inewvc(i),inewvc(k))
        ind2=indx(max(inewvc(j),inewvc(l)),min(inewvc(j),
     &                inewvc(l)))
        if(i.eq.k.or.j.eq.l) then
          scale1=2.0
        endif
        PK(IPKOFF(IJUNK)+INDX2(IND1,IND2,IRPDS2(2*JJUNK-1)))=
     &    PK(IPKOFF(IJUNK)+INDX2(IND1,IND2,IRPDS2(2*JJUNK-1)))-
     &    scale1*SCALE*X
C
   30 CONTINUE
C
      IF(NUT.EQ.ILNBUF) GOTO 3
C
      NAOINT=NAOINT+(NAOBUF-1)*ILNBUF+NUT
      CLOSE(LUINT,STATUS='KEEP')
C
C  Now we need to change the (aa|aa) PK symmetry super-matrices from
C  canonical to square canonical.
C
 200  DO 99 ISPIN=1,(IUHF+1)
       DO 100 IRREP=1,NIRREP
        IJUNK=INDX(IRREP,IRREP)
        DO 101 I=1,NBFIRR(IRREP)
          DO 102 J=1,I
            DO 103 K=1,I
              IF(I.EQ.K) THEN
                ITOP=J
              ELSE
                ITOP=K
              ENDIF
              DO 104 L=1,ITOP
                IND1=INDX(I,J)
                IND2=INDX(K,L)
                PK(IPKOFF(IJUNK+(ISPIN-1)*36)+INDX2(IND2,IND1,
     &                                        IRPDS1(IJUNK)))=
     &           PK(IPKOFF(IJUNK+(ISPIN-1)*36)+INDX2(IND1,IND2,
     &                                         IRPDS1(IJUNK)))
  104         CONTINUE
  103       CONTINUE
  102     CONTINUE
  101   CONTINUE
  100  CONTINUE
   99 CONTINUE
C
      CALL TIMER(1)
cYAU      WRITE(LUOUT,8999)NAOINT,TIMENEW
 8999 FORMAT(T3,'@MKPK2-I, There are ',I10,' unique AO integrals.',/,
     &       T14,'AO integral reading and sorting required ',F10.3,
     &           ' seconds.',/)

C  All done getting stuff. Let's get on with the feature presentation.
C  Need to reopen the IIII/ONE_INT file because it contains stuff needed later.
      IF (SEWARD) THEN
         CALL GFNAME('ONE_INT ',FNAME,ILENGTH)
         OPEN(LUNITSE,FILE=FNAME(1:ILENGTH),
     &        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      ELSE
         CALL GFNAME('IIII    ',FNAME,ILENGTH)
         OPEN(LUINT,FILE=FNAME(1:ILENGTH),STATUS='OLD',
     &        FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      END IF

      RETURN
      END


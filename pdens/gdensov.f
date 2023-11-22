      SUBROUTINE GDENSOV(IRREPXL,IRREPXR,IRREPX,DOV,R0,SCR,
     &                  MAXCOR,IUHF,
     &                  LSTGRL,LSTGTL,LSTGRLOF,LSTGTLOF,LSTTMP,LSTTMPOF,
     &                  LSTR1,LSTL1,LSTR1OFF,LSTL1OFF,LISTR2,LISTL2,
     &                  LSTR2RS,LSTL2RS)
C
C THIS ROUTINE CALCULATES THE OCCUPIED-VIRTUAL COMPONENTS OF THE
C ONE-PARTICLE EXCITED STATE DENSITY MATRIX
C
C                         +
C     DE= <0| L exp(-T) p q exp(T) R |0>
C
C DROV(ia) =    R(em) * L(ae,im) + R0 * L(ia)
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION SCR(MAXCOR),DOV(*)
      DIMENSION I0R1(2),I0L1(2)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
C
      DATA ONE /1.0D0/
C
      I0R1(1)=1
      I0L1(1)=I0R1(1)+IRPDPD(IRREPXR,9)
      IF(LSTR1.NE.-1)CALL GETLST(SCR(I0R1(1)),1,1,1,1+LSTR1OFF,LSTR1)
      CALL GETLST(SCR(I0L1(1)),1,1,1,1+LSTL1OFF,LSTL1)
      IF(IUHF.NE.0)THEN
       I0R1(2)=I0L1(1)+IRPDPD(IRREPXL,9)
       I0L1(2)=I0R1(2)+IRPDPD(IRREPXR,10)
       IF(LSTR1.NE.-1)CALL GETLST(SCR(I0R1(2)),1,1,1,2+LSTR1OFF,LSTR1)
       CALL GETLST(SCR(I0L1(2)),1,1,1,2+LSTL1OFF,LSTL1)
      ELSE
       I0R1(2)=I0R1(1)
       I0L1(2)=I0L1(1)
      ENDIF
C
      ITOP=I0L1(2)+IRPDPD(IRREPXL,10)
      ITOP2=ITOP+IRPDPD(IRREPX,9)+IUHF*IRPDPD(IRREPX,10)
C
      MXCOR=MAXCOR-ITOP+1
C
      IF(LSTR1.NE.-1)THEN
       CALL GT2XF(SCR(ITOP),IRREPXL,IRREPXR,LSTL2RS,
     &            LSTR1,LSTR1OFF,
     &            SCR(ITOP2),MXCOR,IUHF,0)
       CALL SAXPY(IRPDPD(IRREPX,9)+IUHF*IRPDPD(IRREPX,10),ONE,
     &            SCR(ITOP),1,DOV,1)
      ENDIF
      RETURN
      END

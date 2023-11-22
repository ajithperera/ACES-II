      SUBROUTINE DRVSRT(ICORE,MEMSIZ,ILNBUF,NBAS,IUHF,YESNO2)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER(MAXBUK=5000)
C
      LOGICAL YESNO2
C
      DIMENSION ICORE(MEMSIZ)
      DIMENSION IBKNMOF(37),IBKOF(MAXBUK),IBKSTRT(MAXBUK)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
 
#include "icdacc.com"
 
      COMMON /FILES/ LUOUT,MOINTS
#include "symm2.com"
      COMMON /TIMEINFO/ TIMEIN, TIMENOW, TIMETOT, TIMENEW
C
      DATA I37 /37/
C
      INDX(I,J)=J+(I*(I-1))/2
C
      LENINT=1+IINTFP
      NBKINT=(32768-(2*IINTLN))/(IINTLN*(1+IINTFP))
c YAU - old
cC      IRECL=IINTLN*(NBKINT*(1+IINTFP)+2)
cC Nevin
c      if(iintfp.eq.1)then
c        IRECL=IINTLN*(NBKINT*(1+IINTFP)+2)
c      else
c        IRECL=idaccm*(NBKINT*(1+IINTFP)+2)
c      endif
c YAU - new
      IRECL=idaccm*(NBKINT*(1+IINTFP)+2)
c YAU - end
C
C  Calculate the number of buckets required and the bucket lookup
C  vector for placing integrals into buckets.
C
      CALL IZERO(IBKNMOF,I37)
      IBKDIS=0
      MAXDIS=0
      DO 100 I=1,NIRREP
        DO 101 J=1,I
          IF(I.EQ.J) THEN
            K=INDX(I,I)
            IBKNMOF(K+1)=IBKNMOF(K)+IRPDS1(K)
            IBKDIS=IBKDIS+IRPDS1(K)
            MAXDIS=MAX(MAXDIS,IRPDS1(K))
          ELSE
            K=INDX(I,J)
            L=INDX(I-1,J)
            IBKNMOF(K+1)=IBKNMOF(K)+IRPDS2(2*L)
            IBKDIS=IBKDIS+IRPDS2(2*L)
            MAXDIS=MAX(MAXDIS,IRPDS2(2*L))
            MAXDIS=MAX(MAXDIS,IRPDS2(2*L-1))
          ENDIF
  101   CONTINUE
  100 CONTINUE
C
C  Determine how much memory we really have for holding the
C  integrals; i.e., the size of the bucket.
C
      MEMSZ=MEMSIZ
      MEMSZ=MEMSZ-2*MAX(ILNBUF,NBKINT)*LENINT-4*IBKDIS-MAXDIS*IINTFP
C
      I000=1
      I010=I000+IBKDIS
      I020=I010+2*IBKDIS
      I030=I020+IBKDIS
      IF(I030-I000.GT.MEMSIZ) THEN
        CALL NOMEM('Calculate number of PK sort buckets',
     &             '{SRTALC} <-- DRVSRT <-- OOFCOR <-- VSCF',
     &             I030-I000,MEMSIZ)
      ENDIF
C
      CALL SRTALC(ICORE(I020),ICORE(I000),ICORE(I010),IBKOF,IBKSTRT,
     &            IBKDIS,MAXBUK,NBUCK,MEMSZ,IBKSIZ)
C
      WRITE(LUOUT,5000)MEMSZ,IBKDIS,NBUCK,IBKSIZ
 5000 FORMAT(T3,'@DRVSRT-I, PK sort information:',/,
     &       T8,'Total memory available (integer words): ',I10,/,
     &       T8,'         Total number of distributions: ',I10,/,
     &       T8,'                     Number of buckets: ',I10,/,
     &       T8,'            Number of integrals/bucket: ',I10,/)
C
C  I000 = Address for bucket <=> distribution lookup vector
C  I010 = Address for distribution size lookup vector
C  I020 = Address for distribution number and symmetry combination lookup 
C         vector
C  I030 = Address for chain vector
C  I040 = Address for integral bucket buffers
C  I050 = Address for integral indices bucket buffers
C  I060 = Address for integrals-in-bucket-buffer counter
C  I070 = Address for VMOL integral buffer
C  I080 = Address for VMOL integral indices buffer
C
      I040=I030+3*NBUCK
      i040=i040+mod(i040+1,2)
      I050=I040+3*NBUCK*NBKINT*IINTFP
      I060=I050+3*NBUCK*NBKINT
      I070=I060+3*NBUCK+mod(nbuck,2)
      I080=I070+ILNBUF*IINTFP
      I090=I080+ILNBUF
      IF(I090-I000.GT.MEMSIZ) THEN
        CALL NOMEM('Sort AO integrals to PK lists','{PKSORT21} <-- '//
     &             'DRVSRT <-- OOFCOR <-- VSCF',I090-I000,MEMSIZ)
      ENDIF
C
      IF(YESNO2) THEN
        CALL PKSORT21(ICORE(I020),ICORE(I030),ICORE(I040),ICORE(I050),
     &                ICORE(I060),ICORE(I070),ICORE(I080),IBKNMOF,I37,
     &                IRECL,NBKINT,ILNBUF,IBKDIS,NBUCK,IUHF)
C
C  I030 = Address for chain vector
C  I040 = Address for bucket
C  I050 = Address for sortfile bucket integral buffer
C  I060 = Address for sortfile bucket indices buffer
C
        I050=I040+IBKSIZ*IINTFP
        I060=I050+NBKINT*IINTFP
        I070=I060+NBKINT
        IF(I070-I000.GT.MEMSIZ) THEN
          CALL NOMEM('Chain AO integrals to PK lists','{PKCHAN1} <-- '//
     &               'DRVSRT <-- OOFCOR <-- VSCF',I070-I000,MEMSIZ)
      ENDIF
C
        CALL PKCHAN1(ICORE(I030),ICORE(I040),ICORE(I000),ICORE(I010),
     &               ICORE(I050),ICORE(I060),IBKSTRT,IBKOF,MAXBUK,
     &               IBKDIS,IBKSIZ,NBUCK,NBKINT,IUHF)
C
        I040=I030+3*NBUCK
        i040=i040+mod(i040+1,2)
        I050=I040+3*NBUCK*NBKINT*IINTFP
        I060=I050+3*NBUCK*NBKINT
        I070=I060+3*NBUCK+mod(nbuck,2)
        I080=I070+ILNBUF*IINTFP
        I090=I080+ILNBUF
C
        CALL PKSORT22(ICORE(I020),ICORE(I030),ICORE(I040),ICORE(I050),
     &                ICORE(I060),ICORE(I070),ICORE(I080),IBKNMOF,I37,
     &                IRECL,NBKINT,ILNBUF,IBKDIS,NBUCK,IUHF)
C
C  I030 = Address for chain vector
C  I040 = Address for bucket
C  I050 = Address for sortfile bucket integral buffer
C  I060 = Address for sortfile bucket indices buffer
C  I070 = Address for scratch vector for summation of PK lists
C
        I050=I040+IBKSIZ*IINTFP
        I060=I050+NBKINT*IINTFP
        I070=I060+NBKINT
        I080=I070+MAXDIS*IINTFP
        IF(I080-I000.GT.MEMSIZ) THEN
          CALL NOMEM('Chain AO integrals to PK lists','{PKCHAN2} <-- '//
     &               'DRVSRT <-- OOFCOR <-- VSCF',I080-I000,MEMSIZ)
        ENDIF
C
        CALL PKCHAN2(ICORE(I030),ICORE(I040),ICORE(I000),ICORE(I010),
     &               ICORE(I050),ICORE(I060),ICORE(I070),IBKSTRT,IBKOF,
     &               MAXBUK,IBKDIS,IBKSIZ,NBUCK,NBKINT,IUHF,MAXDIS)
C
        I040=I030+3*NBUCK
        i040=i040+mod(i040+1,2)
        I050=I040+3*NBUCK*NBKINT*IINTFP
        I060=I050+3*NBUCK*NBKINT
        I070=I060+3*NBUCK+mod(nbuck,2)
        I080=I070+ILNBUF*IINTFP
        I090=I080+ILNBUF
C
        CALL PKSORT23(ICORE(I020),ICORE(I030),ICORE(I040),ICORE(I050),
     &                ICORE(I060),ICORE(I070),ICORE(I080),IBKNMOF,I37,
     &                IRECL,NBKINT,ILNBUF,IBKDIS,NBUCK,IUHF)
C
C  I030 = Address for chain vector
C  I040 = Address for bucket
C  I050 = Address for sortfile bucket integral buffer
C  I060 = Address for sortfile bucket indices buffer
C  I070 = Address for scratch vector for summation of PK lists
C
        I050=I040+IBKSIZ*IINTFP
        I060=I050+NBKINT*IINTFP
        I070=I060+NBKINT
        I080=I070+MAXDIS*IINTFP
C
        CALL PKCHAN3(ICORE(I030),ICORE(I040),ICORE(I000),ICORE(I010),
     &               ICORE(I050),ICORE(I060),ICORE(I070),IBKSTRT,IBKOF,
     &               MAXBUK,IBKDIS,IBKSIZ,NBUCK,NBKINT,IUHF,MAXDIS)
C
C  Now symmetrize the PK matrices which are square.
C
        CALL TIMER(1)
        DO 300 IRREP=1,NIRREP
          ILIST=(IRREP*(IRREP+1))/2
          CALL SYMMETW3(1,ILIST,IRPDS1(ILIST),ICORE,MEMSIZ/IINTFP)
          IF(IUHF.NE.0) THEN
            CALL SYMMETW3(2,ILIST,IRPDS1(ILIST),ICORE,MEMSIZ/IINTFP)
          ENDIF
  300   CONTINUE
C
        CALL TIMER(1)
        WRITE(LUOUT,4000)TIMENEW
 4000   FORMAT(T3,'@DRVSRT-I, Symmetrization of PK matrix required ',
     &            F10.3,' seconds.',/)
C
      ELSE
        CALL PKSORT(ICORE(I020),ICORE(I030),ICORE(I040),ICORE(I050),
     &              ICORE(I060),ICORE(I070),ICORE(I080),IBKNMOF,I37,
     &              IRECL,NBKINT,ILNBUF,IBKDIS,NBUCK,IUHF)
C
C  I030 = Address for chain vector
C  I040 = Address for bucket
C  I050 = Address for sortfile bucket integral buffer
C  I060 = Address for sortfile bucket indices buffer
C
      I050=I040+IBKSIZ*IINTFP
      I060=I050+NBKINT*IINTFP
      I070=I060+NBKINT
      IF(I070-I000.GT.MEMSIZ) THEN
        CALL NOMEM('Chain AO integrals to PK lists','{PKCHAN} <-- '//
     &             'DRVSRT <-- OOFCOR <-- VSCF',I070-I000,MEMSIZ)
      ENDIF
C
      CALL PKCHAN(ICORE(I030),ICORE(I040),ICORE(I000),ICORE(I010),
     &            ICORE(I050),ICORE(I060),IBKSTRT,IBKOF,MAXBUK,IBKDIS,
     &            IBKSIZ,NBUCK,NBKINT,IUHF)
C
      ENDIF
C
      RETURN
      END

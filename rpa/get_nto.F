      SUBROUTINE GET_NTO(DOO,DVV,EVECOO,EVECVV,CMO,NTO,IRREPC,
     &                  ISPIN,IUHF)
C
C GET NATURAL TRANSITION ORBITALS (NTO) FROM DENSITY MATRICES
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "maxbasfn.par"
      INTEGER DIRPRD,POP,VRT
      INTEGER INDOO(8,2),INDVV(8,2)
      CHARACTER*1 BLANK
      DIMENSION DOO(*),DVV(*),EVECOO(*),EVECVV(*),CMO(*),NTO(*)
      DIMENSION IMAP_A(MAXBASFN),IMAP_B(MAXBASFN)
      COMMON/SYMINF/NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/SYMLOC/ISYMOFF(8,8,25)
      DATA ONE,ONEM,ZILCH /1.0D0,-1.0D0,0.0D0/
      DATA BLANK /" "/
C
      LOGICAL YESNO
#include "machsp.com"
#include "info.com"

#ifdef _DEBUG_NTO
      WRITE(6,*) "GET_NTO START"
#endif

C
C Read SCF MO
C
      CALL GETREC(20,'JOBARC','NBASTOT ',1,NBAS)
C     WRITE(6,*) "NBAS=",NBAS
      IF(ISPIN.EQ.1)THEN
       CALL GETREC(20,'JOBARC','SCFEVCA0',NBAS*NBAS*IINTFP,CMO)
      ELSE
       CALL GETREC(20,'JOBARC','SCFEVCB0',NBAS*NBAS*IINTFP,CMO)
      ENDIF
C
      CALL GETREC(20,"JOBARC","ORBMAP_A",NBAS,IMAP_A)
      CALL GETREC(20,"JOBARC","ORBMAP_A",NBAS,IMAP_B)
#ifdef  _DEBUG_NTO
      WRITE(6,*) "ORBMAP_A" 
      DO I = 1, NBAS
         WRITE(6,*) I, IMAP_A(I)
      ENDDO
#endif 
C
#ifdef _DEBUG_NTO
      WRITE(6,*) 'SCF MO'
      CALL OUTPUT(CMO,1,NBAS,1,NBAS,NBAS,NBAS,1)
#endif
C 
C OPEN FILE 
C
      LUNITN = 30
      INQUIRE(FILE='NTO_EXCIT',EXIST=YESNO)
      IF (YESNO) THEN
       OPEN(LUNITN,FILE='NTO_EXCIT',STATUS='OLD',FORM='FORMATTED',
     &     POSITION='APPEND')
      ELSE
       OPEN(LUNITN,FILE='NTO_EXCIT',STATUS='NEW',FORM='FORMATTED',
     &     POSITION='REWIND')
      ENDIF
      WRITE(LUNITN,'(A,1X,I1)') "Symmetry of the excited State: ",
     &                           IRREPC
      WRITE(LUNITN,'(A)') BLANK
C
      IRREPX=1
      DO I=1,ISPIN
       INDOO(1,I)=0
       INDVV(1,I)=POP(1,I)
       DO IRREPJ=1,NIRREP-1
        INDOO(IRREPJ+1,I)=INDVV(IRREPJ,I)+VRT(IRREPJ,I)
        INDVV(IRREPJ+1,I)=INDOO(IRREPJ+1,I)+POP(IRREPJ+1,I)
       ENDDO
      ENDDO
C
C NTO FROM OCCUPIED-OCCUPIED CONTRIBUTION
C
      ICOLL=1
      IROWL=1
      DO 10 IRREPJ=1,NIRREP
C      Since IRREPX=1, IRREPJ=IRREPI, NUMJ=NUMI
       IRREPI=DIRPRD(IRREPJ,IRREPX)
       NUMJ=POP(IRREPJ,ISPIN)
       NUMI=POP(IRREPI,ISPIN)
       IOFFD =ISYMOFF(IRREPJ,IRREPX,20+ISPIN)
C
C      Get eigenvalues(DOO) and eigenvectors(EVECOO) from DOO matrix
C      DOO(IOFFD) has NUMJ*NUMJ elements
C
       IF (NUMJ.GT.0) THEN
#ifdef _DEBUG_NTO
        WRITE(6,'(A,I1)') "DOO before digaonalization, Block:",
     &   IRREPJ
        CALL OUTPUT(DOO(IOFFD),1,NUMJ,1,NUMJ,NUMJ,NUMJ,1)
#endif
        CALL EIG(DOO(IOFFD),EVECOO(IOFFD),1,NUMJ,1)
#ifdef _DEBUG_NTO
        
        WRITE(6,'(A,I1)') "Eigenvalues after digaonalization, Block:",
     &   IRREPJ
        CALL OUTPUT(DOO(IOFFD),1,NUMJ,1,NUMJ,NUMJ,NUMJ,1)
        DIAGVAL = 0.0D0
        DO I = 1, NUMJ
           IND = NUMJ*(I-1) + I-1
           DIAGVAL = DIAGVAL + DOO(IOFFD+IND)
        ENDDO
        WRITE(6,'(A,I1,A,F15.7)') "Sum of diagonal elements for block(",
     &       IRREPJ,"):",DIAGVAL
        WRITE(6,'(A,I1)') "Eigenvectors after digaonalization, Block:",
     &   IRREPJ
        CALL OUTPUT(EVECOO(IOFFD),1,NUMJ,1,NUMJ,NUMJ,NUMJ,1)
#endif
C       Check corresponding SCF MO
CSS        WRITE(6,*) "NTO OCCUPIED"
        IND = 0
        DO I = 1, NUMJ
         TEMP = 0.0D0 
         IND2 = 0
         DO J = 1, NUMJ
          IF (ABS(EVECOO(IOFFD+IND)).GT.TEMP) THEN
           TEMP = ABS(EVECOO(IOFFD+IND))
           IND2 = J
          ENDIF
          IND = IND + 1
         ENDDO
CSSS         WRITE(6,'(2(A,I4))') "NTO:",I," corresponds to SCFMO:",IND2 
C        IF (I.EQ.1) THEN
C           WRITE(LUNITN,'(A,I7,1X,F15.10)') "OCC ",IND2, 
C    &       EVECOO((I-1)*NUMJ+IND2)
C        ENDIF
        ENDDO
C
C       Rotate HF orbitals to get NTO orbitals
C
        ICOLH = ICOLL+POP(IRREPJ,ISPIN)-1
        IROWH = IROWL+POP(IRREPJ,ISPIN)+VRT(IRREPJ,ISPIN)-1
#ifdef _DEBUG_NTO
        WRITE(6,'(A,I1)') "SCF MO for block:",IRREPJ
        CALL OUTPUT(CMO,IROWL,IROWH,ICOLL,ICOLH,NBAS,NBAS,1)
#endif
C       NTO(NBAS,ICOLL~ICOLH) = CMO(NBAS,ICOLL~ICOLH)*EVECOO(NUMJ,NUMJ)
        IOFFN = (ICOLL-1)*NBAS+1
        IOFFN2= INDOO(IRREPJ,ISPIN)*NBAS+1
        CALL XGEMM('N','N',NBAS,NUMJ,NUMJ,ONE,CMO(IOFFN),NBAS,
     &              EVECOO(IOFFD),NUMJ,ZILCH,NTO(IOFFN2),NBAS)
#ifdef _DEBUG_NTO
C       WRITE(6,'(A,3I4)') "ICOLL/ICOLH/IOFFN",ICOLL,ICOLH,IOFFN
        ICOLL2=INDOO(IRREPJ,ISPIN)+1
        ICOLH2=INDOO(IRREPJ,ISPIN)+POP(IRREPJ,ISPIN)
C       WRITE(6,'(A,3I4)') "ICOLL2/ICOLH2/IOFFN2",ICOLL2,ICOLH2,IOFFN2
C        WRITE(6,'(A,I1)') "NTO for block:",IRREPJ
C       CALL OUTPUT(NTO,IROWL,IROWH,ICOLL,ICOLH,NBAS,NBAS,1)
        CALL OUTPUT(NTO,1,NBAS,ICOLL2,ICOLH2,NBAS,NBAS,1)
#endif
C
        WRITE(LUNITN,'(A,1X,I1)') "Symmetry of occ. orbitals:",IRREPI
        WRITE(LUNITN,*) "IOFFN2/IOFFN/IOFFD",IOFFN2,IOFFN,IOFFD
        WRITE(LUNITN,'(A)') BLANK
        CALL TRACK_OCC_NTOS(NTO(IOFFN2),EVECOO(IOFFD),CMO(IOFFN),
     +                      DOO(IOFFD),NUMI,NBAS,ICOLL,IMAP_A,IMAP_B,
     +                      LUNITN) 
        WRITE(LUNITN,'(A)') BLANK
C
        ICOLL = ICOLH+1
        IROWL = IROWH+1
C
       ENDIF
10    CONTINUE
C
C NTO FROM VIRTUAL-VIRTUAL CONTRIBUTION 
C
      NOCS = NOCCO(1)
      ICOLL=ICOLH+1
      IROWL=1
      DO 20 IRREPB=1,NIRREP
       IRREPA=DIRPRD(IRREPB,IRREPX)
       NUMB=VRT(IRREPB,ISPIN)
       NUMA=VRT(IRREPA,ISPIN)
       IOFFD =ISYMOFF(IRREPB,IRREPX,18+ISPIN)
C
C      Get eigenvalues(DVV) and eigenvectors(EVECVV) from DVV matrix
C      DVV(IOFFD) has NUMB*NUMB elements
C
       IF (NUMB.GT.0) THEN
#ifdef _DEBUG_NTO
        WRITE(6,'(A,I1)') "DVV before digaonalization, Block:",
     &   IRREPB
        CALL OUTPUT(DVV(IOFFD),1,NUMB,1,NUMB,NUMB,NUMB,1)
#endif
C
C       
C
        CALL EIG(DVV(IOFFD),EVECVV(IOFFD),1,NUMB,-1)
#ifdef _DEBUG_NTO
        WRITE(6,'(A,I1)') "Eigenvalues after digaonalization, Block:",
     &   IRREPB
        CALL OUTPUT(DVV(IOFFD),1,NUMB,1,NUMB,NUMB,NUMB,1)
        DIAGVAL = 0.0D0
        DO I = 1, NUMB
           IND = NUMB*(I-1) + I-1
           DIAGVAL = DIAGVAL + DVV(IOFFD+IND)
        ENDDO
        WRITE(6,'(A,I1,A,F15.7)') "Sum of diagonal elements for block(",
     &       IRREPB,"):",DIAGVAL
        WRITE(6,'(A,I1)') "Eigenvectors after digaonalization, Block:",
     &   IRREPB
        CALL OUTPUT(EVECVV(IOFFD),1,NUMB,1,NUMB,NUMB,NUMB,1)
#endif
C       Check corresponding SCF MO
CSSS        WRITE(6,*) "NTO VIRTUAL"
        IND = 0
        DO I = 1, NUMB
         TEMP = 0.0D0
         IND2 = 0
         DO J = 1, NUMB
          IF (ABS(EVECVV(IOFFD+IND)).GT.TEMP) THEN
           TEMP = ABS(EVECVV(IOFFD+IND))
           IND2 = J
          ENDIF
          IND = IND + 1
         ENDDO
CSSS         WRITE(6,'(2(A,I4))') "NTO:",I," corresponds to SCFMO:",IND2
C        IF (I.EQ.1) THEN
C           WRITE(LUNITN,'(A,I7,1X,F15.10)') "VRT ",IND2+NOCS,
C    &                     EVECVV((I-1)*NUMB+IND2)
C        ENDIF
        ENDDO
C
C       Rotate HF orbitals to get NTO orbitals
C
        ICOLH = ICOLL+VRT(IRREPB,ISPIN)-1
        IROWH = IROWL+VRT(IRREPB,ISPIN)+POP(IRREPB,ISPIN)-1
#ifdef _DEBUG_NTO
        WRITE(6,'(A,I1)') "SCF MO for block:",IRREPB
        CALL OUTPUT(CMO,IROWL,IROWH,ICOLL,ICOLH,NBAS,NBAS,1)
#endif
C       NTO(NBAS,ICOLL~ICOLH) = CMO(NBAS,ICOLL~ICOLH)*EVECVV(NUMB,NUMB)
        IOFFN = (ICOLL-1)*NBAS+1
        IOFFN2= INDVV(IRREPB,ISPIN)*NBAS+1
        CALL XGEMM('N','N',NBAS,NUMB,NUMB,ONE,CMO(IOFFN),NBAS,
     &              EVECVV(IOFFD),NUMB,ZILCH,NTO(IOFFN2),NBAS)
#ifdef _DEBUG_NTO
        WRITE(6,'(A,3I4)') "ICOLL/ICOLH/IOFFN",ICOLL,ICOLH,IOFFN
        ICOLL2=INDVV(IRREPB,ISPIN)+1
        ICOLH2=INDVV(IRREPB,ISPIN)+VRT(IRREPB,ISPIN)
        WRITE(6,'(A,3I4)') "ICOLL2/ICOLH2/IOFFN2",ICOLL2,ICOLH2,IOFFN2
        WRITE(6,'(A,I1)') "NTO for block:",IRREPB
C       CALL OUTPUT(NTO,IROWL,IROWH,ICOLL,ICOLH,NBAS,NBAS,1)
        CALL OUTPUT(NTO,1,NBAS,ICOLL2,ICOLH2,NBAS,NBAS,1)
#endif
C
C        WRITE(6,'(A,3I4)') "ICOLL/ICOLH/IOFFN",ICOLL,ICOLH,IOFFN
        WRITE(LUNITN,'(A,1X,I1)') "Symmetry of vrt. orbitals:",IRREPB
        WRITE(LUNITN,*) "IOFFN2/IOFFN/IOFFD",IOFFN2,IOFFN,IOFFD
        WRITE(LUNITN,'(A)') BLANK
        CALL TRACK_VRT_NTOS(NTO(IOFFN2),EVECVV(IOFFD),CMO(IOFFN),
     +                      DVV(IOFFD),NUMB,NBAS,ICOLL,IMAP_A,IMAP_B,
     +                      LUNITN) 
        WRITE(LUNITN,'(A)') BLANK
C
        ICOLL = ICOLH+1
        IROWL = IROWH+1
C
       ENDIF
20    CONTINUE
C
C     Overwrite SCFMO with NTO
C
       IF(ISPIN.EQ.1)THEN
        CALL PUTREC(20,'JOBARC','SCFEVCA0',NBAS*NBAS*IINTFP,NTO)
        CALL PUTREC(20,'JOBARC','SCFEVCB0',NBAS*NBAS*IINTFP,NTO)
       ENDIF
C     ENDIF

      CLOSE(LUNITN)
#ifdef _DEBUG_NTO
      WRITE(6,*) "Full NTO"
      CALL OUTPUT(NTO,1,NBAS,1,NBAS,NBAS,NBAS,1)

      WRITE(6,*) "GET_NTO END"
#endif
      RETURN
      END

      SUBROUTINE MEMUSE
C
C COMPUTES MEMORY REQUIRED TO PERFORM CERTAIN OPERATIONS
C
      IMPLICIT INTEGER (A-Z)
      CHARACTER*30 XXX(15)
      CHARACTER*8 PTGRP
      DIMENSION MEMREQ(15,2)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPY(255,2),DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),D(18) 
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /FILES/ LUOUT,MOINTS
      DATA XXX /'Particle-particle ladder (AB) ',
     &          'Hole-hole ladder (AB)         ',
     &          'Fourth-order singles (AB)     ',
     &          'Rings (AB)                    ',
     &          'T2 vector (AB)                ',
     &          'T2 vector (AA)                ',
     &          'Particle-particle ladder (AA) ',
     &          'Hole-hole ladder (AA)         ',
     &          'Rings (AA)                    ',
     &          'Fourth-order singles (AA)     ',
     &          'Particle-particle ladder (BB) ',
     &          'Hole-hole ladder (BB)         ',
     &          'Rings (BB)                    ',
     &          'Fourth-order singles (BB)     ',
     &          'T2 vector (BB)                '/
      NNM1O2(I)=(I*(I-1))/2
      CALL IZERO(MEMREQ,30)
      NOCA=NOCCO(1)
      NOCB=NOCCO(2)
      NVRTA=NVRTO(1)
      NVRTB=NVRTO(2)
      DO 10 IRREP=1,NIRREP
C
C AB PPL.
C
       MEMT2=IRPDPD(IRREP,13)*IRPDPD(IRREP,14)
       MEMW =IRPDPD(IRREP,13)*IRPDPD(IRREP,13)
       MEMPPL=(MEMW+2*MEMT2)
       MEMREQ(1,1)=MAX(MEMREQ(1,1),MEMPPL)
C
C AB HHL.
C
       MEMW =IRPDPD(IRREP,14)*IRPDPD(IRREP,14)
       MEMHHL=(MEMW+2*MEMT2)
       MEMREQ(2,1)=MAX(MEMREQ(2,1),MEMHHL)

C
C AB RING.
C
       MEMTXX=MAX(IRPDPD(IRREP,9),IRPDPD(IRREP,11),IRPDPD(IRREP,10),
     &            IRPDPD(IRREP,12))
       MEMT2A=MEMTXX*MEMTXX
       MEMWA =MEMT2A
       MEMRNG=(MEMWA+2*MEMT2A)
       MEMREQ(4,1)=MAX(MEMRNG,MEMREQ(4,1))
C
C FULL T2 VECTOR.
C
       MEMREQ(5,1)=ISYMSZ(13,14)
       MEMREQ(6,1)=ISYMSZ(1,3)
10    CONTINUE
      MEMREQ(1,2)=NVRTA*NVRTA*NVRTB*NVRTB+2*(NOCA*NOCB*NVRTA*NVRTB)
      MEMREQ(2,2)=NOCA*NOCA*NOCB*NOCB+2*(NOCA*NOCB*NVRTA*NVRTB)
      MEMREQ(4,2)=3*NOCA*NOCB*NVRTA*NVRTB*IINTFP     
      MEMREQ(5,2)=NVRTA*NVRTB*NOCA*NOCB
      MEMREQ(6,2)=NNM1O2(NVRTA)*NNM1O2(NOCA)
      WRITE(LUOUT,999)
999   FORMAT(72('-'))
      WRITE(LUOUT,1000)
1000  FORMAT(T11,'Core memory requirement for full matrix ',
     &       'algorithms',//,T5,'Contraction or quantity',
     &       T41,'Memory needed (words)')
      CALL GETCREC(20,'JOBARC','COMPPTGP',8,PTGRP)
      WRITE(LUOUT,1001)PTGRP(1:3)
      WRITE(LUOUT,999)
1001  FORMAT(T39,'C1',T62,A3)
      DO 20 ITYPE=1,6
       WRITE(LUOUT,1002)XXX(ITYPE),IINTFP*MEMREQ(ITYPE,2),
     &                  IINTFP*MEMREQ(ITYPE,1)
1002   FORMAT(T3,A30,T34,I10,T58,I10)
20    CONTINUE
      WRITE(LUOUT,999)
      RETURN
      END

      SUBROUTINE SETUPB(EVEC,SCR,NBFATM,IMEMCMP,IPOPCMP,
     &                  IMAP,NBAS,NSIZE,NATOMS,LABELIN,
     &                  LABELOUT)
C
C THIS ROUTINE REORDERS ELEMENTS OF MOLECULAR ORBITAL EIGENVECTORS
C  SO THAT THE ATOMIC ORDERING MATCHES THAT SPECIFIED IN THE ZMAT 
C  FILE.  THE EIGENVECTORS ARE ASSUMED TO BE IN THE FULL AO BASIS
C  (NOT THE SO!) AND ORDERED ACCORDING TO THE VMOL ATOM ORDERING ON 
C  INPUT.  THEY ARE WRITTEN TO JOBARC AT THE END OF THE ROUTINE.
C
CEND 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "maxbasfn.par"
      CHARACTER*4 CMPPGP
      CHARACTER*8 LABELIN,LABELOUT
      DIMENSION EVEC(NBAS,NBAS),SCR(NBAS,NBAS),NBFATM(NATOMS)
      DIMENSION IMEMCMP(NSIZE),IPOPCMP(NSIZE)
      DIMENSION IMAP(NSIZE),IBASOFF(MAXBASFN,2),IDUMMY(MAXBASFN)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      common /centre/ numcen0
C
      ISTART=0
      IONE=1
      ITHR=3*IINTFP
C
C PICK UP SOME INFORMATION FROM THE JOBARC FILE
C
      CALL GETREC(20,'JOBARC',LABELIN,IINTFP*NBAS*NBAS,EVEC)
c      CALL GETCREC(-1,'JOBARC','COMPPTGP',4,CMPPGP)
      CALL GETREC(-1,'JOBARC','COMPORDR',IONE,IORDGP)
      CALL GETREC(-1,'JOBARC','COMPNORB',IONE,NORBIT)
      CALL GETREC(-1,'JOBARC','COMPPOPV',NORBIT,IPOPCMP)
      CALL GETREC(-1,'JOBARC','COMPMEMB',NATOMS,IMEMCMP)
      CALL GETREC(-1,'JOBARC','NBASATOM',NORBIT,NBFATM)
      CALL GETREC(20,'JOBARC','MAP2ZMAT',NATOMS,IMAP)
c
c loop over z-matrix centers and determine which orbit they belong to
c  and compute basis function offsets for zmat ordering.  the ilocate
c  stuff is only to deal with troublesome dummy atoms
c
      ibasoff(1,1)=1
      ibasoff(1,2)=1
      do 1000 icenter=1,natoms-1
       ioff=0
       ilocate=0
       do 1001 iorbit=1,norbit
        numcen=ipopcmp(iorbit)
        do 1002 ninorb=1,ipopcmp(iorbit)
         if(imemcmp(ioff+ninorb).eq.icenter)then
          ibasoff(icenter+1,1)=ibasoff(icenter,1)+
     &    nbfatm(iorbit)/numcen
          ilocate=1
         endif
1002    continue
        ioff=ioff+ipopcmp(iorbit)
1001   continue
       if(ilocate.eq.0)ibasoff(icenter+1,1)=ibasoff(icenter,1)
1000  continue
c
c now determine basis function offsets for vmol ordering
c
      icenter=1
      do 2000 iorbit=1,norbit
       npop=ipopcmp(iorbit)
       nbasx=nbfatm(iorbit)/npop
       do 2001 ninorb=1,npop
         ibasoff(icenter+1,2)=ibasoff(icenter,2)+nbasx
         icenter=icenter+1
2001   continue
2000  continue
C
C
C REORDER THE EIGENVECTORS BY EIGENVALUE
C
c      DO 410 I=1,NBAS
c
c       IDUMMY(I)=I
c410   CONTINUE
c      CALL PIKSR2(NBAS,EVAL,IDUMMY)
c      DO 411 I=1,NBAS
c       IPOS1=IDUMMY(I)
c       IPOS2=I
c       CALL SCOPY(NBAS,EVEC(1,IPOS1),1,SCR(1,IPOS2),1)
c411   CONTINUE
c      CALL SCOPY(NBAS*NBAS,SCR,1,EVEC,1)
C
c      CALL PUTREC(20,'JOBARC','EVALORDR',IINTFP*NBAS,EVAL)
C
C DO REORDERING PROCEDURE NUMBER 2.
C
      CALL GETREC(20,'JOBARC','ANGMOMBF',NBAS,IDUMMY)
      IATOM=0
      DO 500 IORBIT=1,NORBIT
       NPOP=IPOPCMP(IORBIT)
       NBASX=NBFATM(IORBIT)/NPOP
       DO 501 NINORB=1,NPOP
        IATOM=IATOM+1
        IF(NINORB.EQ.1)THEN
         ISTART=IBASOFF(IATOM,2)
        ELSE
         ISTART=ISTART+1
        ENDIF
        IOFF2=IBASOFF(IATOM,2)-1
        DO 502 IBAS=1,NBASX
         IOFF1=ISTART+(IBAS-1)*NPOP
         IOFF2=IOFF2+1
         CALL SCOPY(NBAS,EVEC(IOFF1,1),NBAS,SCR(IOFF2,1),NBAS) 
502     CONTINUE
501    CONTINUE
500   CONTINUE
      CALL SCOPY(NBAS*NBAS,SCR,1,EVEC,1)
C
C NOW REORDER THE EIGENVECTORS
C
C FIRST LOOP OVER THE NUMBER OF COMPUTATIONAL ORBITS
C  
      CALL ZERO(SCR,NBAS*NBAS)
      IVMOL=0
      numcen0=0
      DO 100 IORBIT=1,NORBIT
       NUMCEN=IPOPCMP(IORBIT)
       numcen0=numcen0+numcen
       NBASFN=NBFATM(IORBIT)/NUMCEN
       DO 120 ICENT=1,NUMCEN
        IVMOL=IVMOL+1
        IZMAT=IMAP(IVMOL)
        IOFF=0
C
C NOW LOOP OVER ALL BASIS FUNCTIONS ON THIS CENTER AND COPY
C  EIGENVECTORS TO A SCRATCH ARRAY AND BUILD ANGULAR MOMENTUM
C  AND ATOMIC CENTER ARRAYS ON THE FLY
C
        DO 130 IBAS=1,NBASFN
         ILOC1=IBASOFF(IVMOL,2)+IOFF
         ILOC2=IBASOFF(IZMAT,1)+IOFF
         CALL SCOPY(NBAS,EVEC(ILOC1,1),NBAS,SCR(ILOC2,1),NBAS)
         IOFF=IOFF+1
130     CONTINUE
120    CONTINUE
100   CONTINUE
      CALL PUTREC(20,'JOBARC',LABELOUT,IINTFP*NBAS*NBAS,SCR)
C
      RETURN 
      END

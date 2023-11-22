C      SUBROUTINE OSPHASE(EVEC,SOVRLP,SCR1,SCR2,LDIM1,LDIM2,LDIM3,NBAS,
C     $   ICORE,MAXMEM,IUHF,DENS,PK,ONEH,FOCK,IPKSIZ,INCORE)

      SUBROUTINE OSPHASE(PK,ONEH,DENS,FOCK,EVAL,EVEC,ICORE,MAXMEM,LDIM1,
     &                  LDIM2,NBAS,IPKSIZ,REPULS,ETOT,IUHF,LUINT,IOS,
     &                  SOVRLP,SCR1,SCR2,INCORE,
     &                  naobasfn,scfks,scfksexact,scfkslastiter,
     &                  V,z1,ksa,ksb,
     &                  screxc,scr3,scr4,
     &                  valao,valgradao,totwt,
     &                  max_angpts,natoms,
     &                  intnumradpts,ncount,kshf)

C
C     THIS ROUTINE ENSURES THAT ALPHA AND BETA ORBITALS HAVE THE SAME PHASE
C     AS REQUIRED BY OPEN-SHELL SINGLET CALCULATIONS
C
CEND
C
CPROGRAMED BY PS JAN 93
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
#include <maxbasfn.par>
#include "symm2.com"
C
      LOGICAL CHNGFLG,INCORE,AOFIL
      DIMENSION ONEH(LDIM1),DENS((IUHF+1)*LDIM1),FOCK((IUHF+1)*LDIM1)
      DIMENSION ICORE(MAXMEM),PK(IPKSIZ)
      DIMENSION EVEC(2*LDIM1),SOVRLP(NBAS*NBAS)
      DIMENSION SCR1(LDIM2),SCR2(LDIM2), DOCC(MAXBASFN*2)
      DIMENSION NOCC(16), NSUM(16)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /POPUL/ NOCC
      COMMON /FOCK/AOFIL
C
      DATA ONE /1.0D0/
      DATA ONEM /-1.0D0/
      DATA TWO /2.0D0/
      DATA ZILCH /0.0D0/
C
      INDX2(I,J,N)=I+(J-1)*N
      CHNGFLG=.FALSE.
C
C  Form overlap of alpha and beta density matrices.
C
      CALL GETREC(20,'JOBARC','AOOVRLAP',NBAS*NBAS*IINTFP,SOVRLP)
      DO 100 I=1,NIRREP
        IF(NOCC(I).EQ.0) GOTO 100
        CALL GETBLK(SOVRLP,SCR1,NBFIRR(I),NBAS,IREPS(I))
        CALL XGEMM('N','N',NBFIRR(I),NBFIRR(I),NBFIRR(I),ONE,
     &             SCR1,NBFIRR(I),EVEC(ISQROF(I)),NBFIRR(I),
     &             ZILCH,SCR2,NBFIRR(I))
        CALL XGEMM('T','N',NBFIRR(I),NBFIRR(I),NBFIRR(I),ONE,
     &             EVEC(LDIM1+ISQROF(I)),NBFIRR(I),SCR2,NBFIRR(I),
     &             ZILCH,SCR1,NBFIRR(I))
C
        IOFF=LDIM1+ISQROF(I)
        DO 110 J=1,NBFIRR(I)
           IF(SCR1((J-1)*NBFIRR(I)+J).LT.0) THEN
              CALL SSCAL(NBFIRR(I),ONEM,EVEC(IOFF),1)
              CHNGFLG=.TRUE.
           ENDIF
           IOFF=IOFF+NBFIRR(I)
  110   CONTINUE
  100 CONTINUE
C
      CHNGFLG = .FALSE. 
      IF(CHNGFLG) THEN
C
        write(6,"(a)") ' Warning! Some orbitals change sign in osphase.'
        write(6,"(a)") ' The Fock matrices are being rebulit. '
C
C  Now construct a new density matrix.
C
         I000=1
         I010=I000+MXIRR2*IINTFP
         I015=I010+MXIRR2*IINTFP
         IF(I015-I000.GT.MAXMEM) THEN
          CALL NOMEM('Construct density','{MKDEN} <-- OSPHASE <-- VSCF',
     &          I010-I000,MAXMEM)
         ENDIF
C
#ifdef _DPOCC
        CALL MKDP_OCCNOS(DOCC, NOCC, NSUM, NBAS, IUHF, .TRUE.)
#endif
         DO 1000 ISPIN=1,(IUHF+1)
            IROFF = (ISPIN - 1)*NBAS
            DO 1001 I=1,NIRREP
               IF(NBFIRR(I).EQ.0) GOTO 1001
#ifdef _DPOCC
              IOFF = IROFF + 1
              CALL MKDDEN_4IRREP(EVEC((ISPIN-1)*LDIM1+ISQROF(I)),
     &                      ICORE(I000),ICORE(I010),
     &                      NSUM(I+8*(ISPIN-1)),
     &                      DOCC(IOFF),NBFIRR(I),ISPIN,IUHF)
              IROFF = IROFF + NBFIRR(I)
#else
               CALL MKDEN(EVEC((ISPIN-1)*LDIM1+ISQROF(I)),ICORE(I000),
     &            NOCC(I+8*(ISPIN-1)),NBFIRR(I),IUHF)
#endif
               CALL SQUEZ2(ICORE(I000),DENS((ISPIN-1)*LDIM3+ITRIOF(I)),
     &            NBFIRR(I))
 1001       CONTINUE
 1000    CONTINUE
C
C  Now construct the new Fock matrix.
C
      
      IF (AOFIL) THEN
         ilnbuf=600
         i000 = 1
         i010=i000+itriln(nirrep+1)*iintfp
         i020=i010+iintfp*ilnbuf
         i030=i020+ilnbuf
         i040=i030+nbas*nbas
         if(i040-i000.gt.maxmem) then
            call nomem('make fock matrix',
     &           '{mk.hff} <-- uno_ref <-- vscf',
     &           i040-i000,maxmem)
         endif
         if (iuhf .eq. 0) then
            call mkrhff(fock, dens, icore(i000),oneh,
     &           icore(i010),icore(i020),itriln(nirrep+1),
     &           nbas,nbfirr,icore(i030),ilnbuf,luint,.true.,
     &           naobasfn,0,scfks,scfksexact,scfkslastiter,
     &           V,z1,ksa,ksb,screxc,scr,scr2,valao,valgradao,
     &           totwt,max_angpts,natoms,intnumradpts,ncount,kshf)
         else
            call mkuhff(fock,fock(1+itriln(nirrep+1)),dens,
     &           dens(1+itriln(nirrep+1)),icore(i000),oneh,
     &           icore(i010),icore(i020),itriln(nirrep+1),
     &           nbas,nbfirr,icore(i030),ilnbuf,luint,.true.,
     &           naobasfn,iuhf,scfks,scfksexact,scfkslastiter,
     &           V,z1,ksa,ksb,screxc,scr,scr2,valao,valgradao,
     &           totwt,max_angpts,natoms,intnumradpts,ncount,kshf)
         endif

      ELSE

         I010=I000+(IUHF+1)*ITRILN(NIRREP+1)*IINTFP
         I020=I010+MXIRR2*IINTFP
         I030=I020+NBAS*NBAS*IINTFP
         IMEM=MAXMEM-(I030-I000)
         IF(I030-I000.GT.MAXMEM) THEN
          CALL NOMEM('Make Fock matrix','{MKFOCK} <-- OSPHASE <-- VSCF',
     &         I030-I000,MAXMEM)
         ENDIF
C
         IF(INCORE) THEN 
            CALL MKFOCK(PK,ONEH,FOCK,DENS,ICORE(I000),ICORE(I010),
     &         ICORE(I020),IPKSIZ,ITRILN(NIRREP+1),
     &         MXIRR2,NBAS,IUHF)
         ELSE
            CALL MKFOCK2(ONEH,FOCK,DENS,ICORE(I000),ICORE(I010),
     &         ICORE(I020),ICORE(I030),ITRILN(NIRREP+1),
     &         MXIRR2,IMEM,NBAS,IUHF)
         ENDIF

      ENDIF
      ENDIF 
C
      RETURN
      END

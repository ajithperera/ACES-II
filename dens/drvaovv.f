C
C DRIVER FOR AO-BASED I(ab) = G(aecd)*<be||cd> contraction
C modified slightly from DRAOLAD
c
c Started by PR 12/2003, but left incomplete. Substantial 
C amount of work was need to salvage this to a correctly functioning
C driver for the AO-based formation of I(ab).  Ajith Perera, 04/2005. 
C 
      SUBROUTINE DRVAOVV(AIVV,ICORE,MAXCOR,IUHF,LISTFLAG,IRREPX,
     &                   LSTMO,LSTMOINC,LSTAO,LSTAOINC)
C
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
C
      LOGICAL TAU,MBPT3,MBPT4,TRPEND,SNGEND,GRAD,MBPTT,SING1,
     &        QCISD,ROHF4,ITRFLG,LAMBDA,UCC,CC,M4DQ,M4SDQ,
     &        M4SDTQ,CCSD,CCD
      INTEGER POP ,VRT
      LOGICAL SING
      DIMENSION AIVV(1), IT1OFF(2)
      DIMENSION ICORE(MAXCOR)
C
      COMMON /SYMINF/NSTART,NIRREP,IRREPY(255,2),DIRPRD(8,8)
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/ LUOUT,MOINTS
      COMMON/METH/MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC
      COMMON /SYM/POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF1BB,NF2AA,NF2BB
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),NJUNK(18)
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /ROHF/ ROHF4,ITRFLG
C
      CC   = CCD .OR. CCSD .OR. QCISD
      MBPT4= M4DQ .OR. M4SDQ .OR. M4SDTQ  
C   
      IF (MBPT3 .OR. UCC) THEN
         NTIMES = 1
      ELSE 
         NTIMES = 2
      ENDIF
C
      TAU = .FALSE. 
      IF (CC) THEN
         TAU = .TRUE. 
         IT1OFF(1) = 1
         IT1OFF(2) = IT1OFF(1) + IINTFP*NTAA*IUHF
         IT1END    = IT1OFF(2) + IINTFP*NTBB 
      ELSE
         IT1END    = 1 
      ENDIF 
C
C In this call to AOLADLST, LISTFLAG=0,IRREPX=1. 
C
      CALL AOLADLST(IUHF,LISTFLAG,IRREPX) 
C
C Here in this call to T2TOAO, LSTAO is 243. 
C 
      DO ITIME = 1, NTIMES
C
C ITIME=1 transform T, and ITIME=2 transform Lambda
C
          IF (ITIME .EQ. 1) CALL T2TOAO(ICORE,MAXCOR,IUHF,TAU,43,
     &                                  LSTAO,IRREPX)
          IF (ITIME .EQ. 2) CALL T2TOAO(ICORE,MAXCOR,IUHF,.FALSE.,
     &                                  143, LSTAO,IRREPX)
c
C Do the contraction of Tau(ij,munu)<munu||lamsig> = Z(ij,lamsig)
C when ITIME=1 and when ITIME=2, do the L(ij,munu)<munu||lamsig> =
C Z(ij,lamsig). LSTAOINC=260. AOLAD2 is the single pass algorithm
C and the AOLAD3 is the multipass algorithm. 
C
          IF (IFLAGS(95).EQ.1)THEN
             CALL AOLAD2(ICORE,MAXCOR,IUHF,.FALSE.,IRREPX,LSTAO,
     &                   LSTAOINC)
          ELSE
             CALL AOLAD3(ICORE,MAXCOR,IUHF,.FALSE.,IRREPX,LSTAO,
     &                   LSTAOINC)
          ENDIF
C
C Now we should have a intermediate quantity of the kind Z(mn,lamsig)
C Convert the Z(mn,lamsig) ---> Z(mn,ab), The LSTMOINC (63-66) contains 
C the Z(mn,ab) spin combinations. 
C 
          CALL Z2TOMO(ICORE,MAXCOR,IUHF,.FALSE.,IRREPX,LSTMOINC,
     &                LSTMOINC,LSTAOINC,.TRUE.)
C
C Do the Z(mn,af)L(or T)(mn,bf) to get the I(ab) contribution. 
C 
          CALL Z2INIVV(AIVV,ICORE,MAXCOR,IUHF,TAU,IT1OFF,IT1END, 
     &                 ITIME)
      ENDDO

      RETURN
      END

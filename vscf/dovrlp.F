
#include "flags.h"

      SUBROUTINE DOVRLP(OVRLP,KINETIC,SCR1,SCR2,SCR3,SCRTMP,SCRSQR,
     &                  ZDEP,ICOL,LDIM1,LDIM2,NBAS,ITOTAL)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INTEGER SCFPRNT
      DOUBLE PRECISION KINETIC
      DIMENSION OVRLP(LDIM1),SCR1(LDIM2),SCR2(LDIM2),SCR3(LDIM2)
      DIMENSION SCRTMP(LDIM2),SCRSQR(NBAS*NBAS),ZDEP(NBAS),ICOL(NBAS)
      DIMENSION KINETIC(LDIM1)
      DIMENSION NLINDEP_PER_IRREP(8)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/ LUOUT,MOINTS
#include "flags.com"
#include "flags2.com"
      COMMON /KILLS/ KILLFLG(3)
#include "symm2.com"
      COMMON /SCFPRINT/ SCFPRNT
C
      DATA ZILCH /0.0/
      DATA ONE /1.0/
C
      INDX2(I,J,N)=I+(J-1)*N
C
C  Dump out the overlap matrix before we overwrite it.
C
      CALL ZERO(SCRSQR,NBAS*NBAS)
      CALL MKFULL(SCRSQR,OVRLP,SCR1,LDIM1,LDIM2,NBAS,1)
      CALL PUTREC(20,'JOBARC','AOOVRLAP',NBAS*NBAS*IINTFP,SCRSQR)
#ifdef _DEBUG_LVLM
      Write(6,"(a)") "The overlap matrix"
      call output(SCRSQR,1,nbas, 1,nbas, nbas, nbas, 1)
      call checksum("@-VSCF:Ovlrp",SCRSQR,nbas*nbas)
#endif
      CALL ZERO(SCRSQR,NBAS*NBAS)
      CALL MKFULL(SCRSQR,KINETIC,SCR1,LDIM1,LDIM2,NBAS,1)
      CALL PUTREC(20,'JOBARC','KINETINT',NBAS*NBAS*IINTFP,SCRSQR)
C
      IEXP=-ABS(IFLAGS(8))
      TOL=10.**IEXP
      ITOTAL=0
      CALL ZERO(SCRSQR,NBAS*NBAS)
C
      IF(SCFPRNT.GE.10) THEN
        WRITE(LUOUT,5300)
 5300   FORMAT(T3,'@DOVRLP-I, Eigenvalues of the diagonalized ',
     &            'overlap matrix:',/)
      ENDIF
C
      DO 100 I=1,NIRREP
        NLINDEP=0
        IF(NBFIRR(I).EQ.0)GOTO 100
        CALL EXPND2(OVRLP(ITRIOF(I)),SCR1,NBFIRR(I))
        CALL EIG(SCR1,SCR2,NBFIRR(I),NBFIRR(I),1)
        DO 101 J=1,NBFIRR(I)
          DIAG=SCR1(INDX2(J,J,NBFIRR(I)))
          IF(DIAG.GT.TOL) THEN
            SCALE=ONE/SQRT(DIAG)
          ELSE
            SCALE=ZILCH
            NLINDEP=NLINDEP+1
            ZDEP(NLINDEP)=DIAG
            ICOL(NLINDEP)=J
          ENDIF
          SCR1(INDX2(J,J,NBFIRR(I)))=SCALE
          SCRSQR(J)=DIAG
  101   CONTINUE
        NLINDEP_PER_IRREP(I) = NLINDEP 
C
        IF(SCFPRNT.GE.10) THEN
          WRITE(LUOUT,5050)I
 5050     FORMAT(T7,'Eigenvalues of irrep ',I1,'.',/)
          WRITE(LUOUT,5051)(SCRSQR(K),K=1,NBFIRR(I))
 5051     FORMAT((T9,4(F14.10,3X)))
          WRITE(LUOUT,5101)
        ENDIF
C
C  Now print out the info on the dependent basis functions.
C
        IF(NLINDEP.NE.0) THEN
          IF(IFLAGS( 3).GT.0) KILLFLG(2)=1
          IF(IFLAGS(18).EQ.1 .OR.
     &       IFLAGS(18).EQ.2 .OR.
     &       IFLAGS(18).EQ.3) KILLFLG(2)=1
          IF(IFLAGS2(h_IFLAGS2_geom_opt).NE.0) KILLFLG(2)=1
          WRITE(LUOUT,5000)NLINDEP,I
 5000     FORMAT(T3,'@DOVRLP-W, ',I3,' linear dependencies found for ',
     &              'irrep ',I1,'.',/)
          DO 200 K=1,NLINDEP
            IOFF=INDX2(1,ICOL(K),NBAS)
            WRITE(LUOUT,5100)ZDEP(K),(SCR2(IOFF+L-1),L=1,NBAS)
 5100       FORMAT(T7,'Eigenvalue: ',F20.10,/,/,(T7,5(F10.6,4X)))
            WRITE(LUOUT,5101)
 5101       FORMAT(/)
  200     CONTINUE
        ENDIF
        ITOTAL=ITOTAL+NLINDEP
C
C  Now form the orthogonalization matrix.
C
        CALL TRANSP(SCR2,SCR3,NBFIRR(I),NBFIRR(I))
        CALL MXM(SCR2,NBFIRR(I),SCR1,NBFIRR(I),SCRTMP,NBFIRR(I))
        CALL MXM(SCRTMP,NBFIRR(I),SCR3,NBFIRR(I),SCR1,NBFIRR(I))
        CALL SQUEZ2(SCR1,OVRLP(ITRIOF(I)),NBFIRR(I))
  100 CONTINUE
        CALL PUTREC(20,"JOBARC","NLINDEP ",8,NLINDEP_PER_IRREP) 
C
C  Write out some more info to the user so he's not completely lost
C  when he sees eigenvectors and eigenvalues which are zeros.
C
      IF(ITOTAL.NE.0) WRITE(LUOUT,5200)ITOTAL
 5200 FORMAT(T3,'@DOVRLP-I, A total of ',I3,' linear dependencies ',
     &          'were found in the basis.',/,
     &       T14,'The dependencies were eliminated from the ',
     &           'calculation.',/,
     &       T14,'There will be eigenvalues and eigenvectors of ',
     &           'zero.',/)
C
      RETURN
      END

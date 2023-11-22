      SUBROUTINE MKFOCK(PK,ONEH,FOCK,DENS,SCR1,SCR2,scr3,IPKSIZ,LDIM1,
     &                  LDIM2,nbas,IUHF)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      character*8 cfock(2)
      DIMENSION ONEH(LDIM1),FOCK((IUHF+1)*LDIM1),DENS((IUHF+1)*LDIM1)
      DIMENSION PK(IPKSIZ),scr3(nbas*nbas)
      DIMENSION SCR1(LDIM1),SCR2(LDIM2)
      DIMENSION IPKOFF(73)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/ LUOUT,MOINTS
      common /flags/ iflags(100)
#include "symm2.com"
      COMMON /PKOFF/ IPKOFF
      COMMON /TIMEINFO/ TIMEIN, TIMENOW, TIMETOT, TIMENEW
C
      DATA ONE /1.0/
      DATA TWO /2.0/
      DATA HALF /0.5/
      DATA ZILCH /0.0/
      data cfock /'FOCKA   ','FOCKB   '/
C
      INDX(I,J)=J+I*(I-1)/2
C
C  Here we must appropriately scale the density matrix for use
C  in its symmetry packed form.
C
      CALL TIMER(1)
      CALL ZERO(SCR1,LDIM1)
C
      DO 100 ISPIN=1,(IUHF+1)
c YAU : old
c       CALL ICOPY(ITRILN(NIRREP+1)*IINTFP,ONEH,1,FOCK(((ISPIN-1)*LDIM1)+1),1)
c       CALL ICOPY(ITRILN(NIRREP+1)*IINTFP,DENS(((ISPIN-1)*LDIM1)+1),1,SCR1,1)
c YAU : new
        CALL DCOPY(ITRILN(NIRREP+1),ONEH,1,FOCK(((ISPIN-1)*LDIM1)+1),1)
        CALL DCOPY(ITRILN(NIRREP+1),DENS(((ISPIN-1)*LDIM1)+1),1,SCR1,1)
c YAU : end
        DO 101 I=1,NIRREP
          IF(NBFIRR(I).EQ.0) GOTO 101
          CALL SSCAL(ITRILN(I),TWO,SCR1(ITRIOF(I)),1)
          CALL SCLDEN(SCR1(ITRIOF(I)),ITRILN(I),SCR2,NBFIRR(I),HALF)
          DO 102 J=1,I
            IF(NBFIRR(J).EQ.0) GOTO 102
            IF(I.EQ.J) THEN
              IJUNK=INDX(I,J)
              CALL XGEMM('N','N',ITRILN(I),1,IRPDS1(IJUNK),ONE,
     &                   PK(IPKOFF(IJUNK)+1),IRPDS1(IJUNK),
     &                   SCR1(ITRIOF(I)),ITRILN(I),
     &                   ONE,FOCK(((ISPIN-1)*LDIM1)+ITRIOF(I)),
     &                   ITRILN(I))
            ELSE
              IJUNK=INDX(I,J)
              JJUNK=INDX(I-1,J)
              CALL XGEMM('T','N',ITRILN(J),1,IRPDS2(2*JJUNK-1),ONE,
     &                   PK(IPKOFF(IJUNK)+1),IRPDS2(2*JJUNK-1),
     &                   SCR1(ITRIOF(I)),ITRILN(I),
     &                   ONE,FOCK(((ISPIN-1)*LDIM1)+ITRIOF(J)),
     &                   ITRILN(J))
              CALL XGEMM('N','N',ITRILN(I),1,IRPDS2(2*JJUNK),ONE,
     &                   PK(IPKOFF(IJUNK)+1),IRPDS2(2*JJUNK-1),
     &                   SCR1(ITRIOF(J)),ITRILN(J),
     &                   ONE,FOCK(((ISPIN-1)*LDIM1)+ITRIOF(I)),
     &                   ITRILN(I))
            ENDIF
  102     CONTINUE
  101   CONTINUE
  100 CONTINUE
      IF(IUHF.EQ.1) THEN
        DO 200 ISPIN=1,(IUHF+1)
c YAU : old
c         CALL ICOPY(ITRILN(NIRREP+1)*IINTFP,DENS(((2-ISPIN)*LDIM1)+1),1,SCR1,1)
c YAU : new
          CALL DCOPY(ITRILN(NIRREP+1),
     &               DENS(((2-ISPIN)*LDIM1)+1),1,SCR1,1)
c YAU : end
          DO 201 I=1,NIRREP
            IF(NBFIRR(I).EQ.0) GOTO 201
            CALL SSCAL(ITRILN(I),TWO,SCR1(ITRIOF(I)),1)
            CALL SCLDEN(SCR1(ITRIOF(I)),ITRILN(I),SCR2,NBFIRR(I),HALF)
            DO 202 J=1,I
              IF(NBFIRR(J).EQ.0) GOTO 202
              IF(I.EQ.J) THEN
                IJUNK=INDX(I,J)
                CALL XGEMM('N','N',ITRILN(I),1,IRPDS1(IJUNK),ONE,
     &                     PK(IPKOFF(36+IJUNK)+1),IRPDS1(IJUNK),
     &                     SCR1(ITRIOF(I)),ITRILN(I),
     &                     ONE,FOCK(((ISPIN-1)*LDIM1)+ITRIOF(I)),
     &                     ITRILN(I))
              ELSE
                IJUNK=INDX(I,J)
                JJUNK=INDX(I-1,J)
                CALL XGEMM('T','N',ITRILN(J),1,IRPDS2(2*JJUNK-1),ONE,
     &                     PK(IPKOFF(36+IJUNK)+1),IRPDS2(2*JJUNK-1),
     &                     SCR1(ITRIOF(I)),ITRILN(I),ONE,
     &                     FOCK(((ISPIN-1)*LDIM1)+ITRIOF(J)),
     &                     ITRILN(J))
                CALL XGEMM('N','N',ITRILN(I),1,IRPDS2(2*JJUNK),ONE,
     &                     PK(IPKOFF(36+IJUNK)+1),IRPDS2(2*JJUNK-1),
     &                     SCR1(ITRIOF(J)),ITRILN(J),ONE,
     &                     FOCK(((ISPIN-1)*LDIM1)+ITRIOF(I)),
     &                     ITRILN(I))
              ENDIF
  202       CONTINUE
  201     CONTINUE
  200   CONTINUE
      ENDIF
C
      CALL TIMER(1)
C
      if(iflags(16).eq.0) then
        do 500 ispin=1,(iuhf+1)
          call zero(scr3,nbas*nbas)
          call mkfull(scr3,fock((ispin-1)*ldim1+1),scr2,ldim1,ldim3,
     &                nbas,1)
          call putrec(20,'JOBARC',cfock(ispin),nbas*nbas*iintfp,scr3)
  500   continue
        write(luout,5010)
 5010   format(t3,'@MKFOCK-I, Construction of new Fock matrix using ',
     &            'correlated density complete.',/)
      endif
      RETURN
      END

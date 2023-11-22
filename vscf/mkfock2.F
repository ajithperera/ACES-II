      SUBROUTINE MKFOCK2(ONEH,FOCK,DENS,SCR1,SCR2,scr3,PK,LDIM1,LDIM2,
     &                   MEMSIZ,nbas,IUHF)
C  
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      character*8 cfock(2)
      DIMENSION ONEH(LDIM1),FOCK((IUHF+1)*LDIM1),DENS((IUHF+1)*LDIM1)
      DIMENSION PK(MEMSIZ),scr3(nbas*nbas)
      DIMENSION SCR1(LDIM1),SCR2(LDIM2)

      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/ LUOUT,MOINTS
      common /flags/ iflags(100)
#include "symm2.com"
      COMMON /TIMEINFO/ TIMEIN, TIMENOW, TIMETOT, TIMENEW
C
      DATA ONE /1.0/
      DATA TWO /2.0/
      DATA HALF /0.5/
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
              IJUNK=INDX(I,I)
C
C  Now we determine how many distributions of PK lists we can hold
C  in-core at once, and loop until we've done them all.
C
              ISTART=1
              INCR=MEMSIZ/(IRPDS1(IJUNK)*IINTFP)
              IF(INCR.LT.IRPDS1(IJUNK)) THEN
                IGET=INCR
              ELSE
                IGET=IRPDS1(IJUNK)
              ENDIF
C
C             CALL GETLST(PK,1,IRPDS1(IJUNK),2,1,IJUNK)
C
 3000         CALL GETLST(PK,ISTART,IGET,2,1,IJUNK)
              Write(6,*) IGET, IJUNK
              CALL XGEMM('N','N',ITRILN(I),1,IGET,ONE,PK,IRPDS1(IJUNK),
     &                   SCR1(ITRIOF(I)+ISTART-1),ITRILN(I),
     &                   ONE,FOCK(((ISPIN-1)*LDIM1)+ITRIOF(I)),
     &                   ITRILN(I))
C
              ISTART=ISTART+IGET
              IF(ISTART.LE.IRPDS1(IJUNK)) THEN
                NREM=IRPDS1(IJUNK)-ISTART+1
                IF(NREM.GT.INCR) THEN
                  IGET=INCR
                ELSE
                  IGET=NREM
                ENDIF
                GOTO 3000
              ENDIF
C
            ELSE
              JJUNK=INDX(I-1,J)
C
C
              ISTART=1
              INCR=MEMSIZ/(IRPDS2(2*JJUNK-1)*IINTFP)
              IF(INCR.LT.IRPDS2(2*JJUNK)) THEN
                IGET=INCR
              ELSE
                IGET=IRPDS2(2*JJUNK)
              ENDIF
C
C             CALL GETLST(PK,1,IRPDS2(2*JJUNK),2,1,36+JJUNK)
C
 3001         CALL GETLST(PK,ISTART,IGET,2,1,36+JJUNK)
              CALL XGEMM('T','N',IGET,1,IRPDS2(2*JJUNK-1),ONE,
     &                   PK,IRPDS2(2*JJUNK-1),SCR1(ITRIOF(I)),ITRILN(I),
     &                   ONE,FOCK(((ISPIN-1)*LDIM1)+ITRIOF(J)+ISTART-1),
     &                   ITRILN(J))
              CALL XGEMM('N','N',ITRILN(I),1,IGET,ONE,
     &                   PK,IRPDS2(2*JJUNK-1),
     &                   SCR1(ITRIOF(J)+ISTART-1),ITRILN(J),
     &                   ONE,FOCK(((ISPIN-1)*LDIM1)+ITRIOF(I)),
     &                   ITRILN(I))
C
C
              ISTART=ISTART+IGET
              IF(ISTART.LE.IRPDS2(2*JJUNK)) THEN
                NREM=IRPDS2(2*JJUNK)-ISTART+1
                IF(NREM.GT.INCR) THEN
                  IGET=INCR
                ELSE
                  IGET=NREM
                ENDIF
                GOTO 3001
              ENDIF
C
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
                IJUNK=INDX(I,I)
C
C
                ISTART=1
                INCR=MEMSIZ/(IRPDS1(IJUNK)*IINTFP)
                IF(INCR.LT.IRPDS1(IJUNK)) THEN
                  IGET=INCR
                ELSE
                  IGET=IRPDS1(IJUNK)
                ENDIF
C
C               CALL GETLST(PK,1,IRPDS1(IJUNK),2,2,IJUNK)
C
 3100           CALL GETLST(PK,ISTART,IGET,2,2,IJUNK)
                CALL XGEMM('N','N',ITRILN(I),1,IGET,ONE,PK,
     &                     IRPDS1(IJUNK),SCR1(ITRIOF(I)+ISTART-1),
     &                     ITRILN(I),ONE,
     &                     FOCK(((ISPIN-1)*LDIM1)+ITRIOF(I)),
     &                     ITRILN(I))
C
C
                ISTART=ISTART+IGET
                IF(ISTART.LE.IRPDS1(IJUNK)) THEN
                  NREM=IRPDS1(IJUNK)-ISTART+1
                  IF(NREM.GT.INCR) THEN
                    IGET=INCR
                  ELSE
                    IGET=NREM
                  ENDIF
                  GOTO 3100
                ENDIF
C
C
              ELSE
                JJUNK=INDX(I-1,J)
C
C
                ISTART=1
                INCR=MEMSIZ/(IRPDS2(2*JJUNK-1)*IINTFP)
                IF(INCR.LT.IRPDS2(2*JJUNK)) THEN
                  IGET=INCR
                ELSE
                  IGET=IRPDS2(2*JJUNK)
                ENDIF
C
C               CALL GETLST(PK,1,IRPDS2(2*JJUNK),2,2,36+JJUNK)
C
 3101           CALL GETLST(PK,ISTART,IGET,2,2,36+JJUNK)
                CALL XGEMM('T','N',IGET,1,IRPDS2(2*JJUNK-1),ONE,
     &                     PK,IRPDS2(2*JJUNK-1),SCR1(ITRIOF(I)),
     &                     ITRILN(I),ONE,
     &                     FOCK(((ISPIN-1)*LDIM1)+ITRIOF(J)+ISTART-1),
     &                     ITRILN(J))
                CALL XGEMM('N','N',ITRILN(I),1,IGET,ONE,PK,
     &                     IRPDS2(2*JJUNK-1),SCR1(ITRIOF(J)+ISTART-1),
     &                     ITRILN(J),ONE,
     &                     FOCK(((ISPIN-1)*LDIM1)+ITRIOF(I)),
     &                     ITRILN(I))
C
C
                ISTART=ISTART+IGET
                IF(ISTART.LE.IRPDS2(2*JJUNK)) THEN
                  NREM=IRPDS2(2*JJUNK)-ISTART+1
                  IF(NREM.GT.INCR) THEN
                    IGET=INCR
                  ELSE
                    IGET=NREM
                  ENDIF
                  GOTO 3101
                ENDIF
C
C
              ENDIF
  202       CONTINUE
  201     CONTINUE
  200   CONTINUE
      ENDIF
C
      CALL TIMER(1)
C     WRITE(LUOUT,5000) TIMENEW
C5000 FORMAT(T3,'@MKFOCK2-I, Construction of Fock matrix required ',
C    &          F10.3,' seconds.',/)
C
      if(iflags(16).eq.0) then
        do 500 ispin=1,(iuhf+1)
          call zero(scr3,nbas*nbas)
          call mkfull(scr3,fock((ispin-1)*ldim1+1),scr2,ldim1,ldim3,
     &                nbas,1)
          call putrec(20,'JOBARC',cfock(ispin),nbas*nbas*iintfp,scr3)
  500   continue
        write(luout,5010)
 5010   format(t3,'@MKFOCK2-I, Construction of new Fock matrix using ',
     &            'correlated density complete.',/)
      endif
      RETURN
      END

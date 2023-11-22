      SUBROUTINE MODF_AMPS(ICORE,MAXCOR,IUHF,NLIST,SING1,TYPE,ACT,
     &                     MODE)
C
C This routine drives number of routine which will be used
C to change the nature of the CC/EOM vectors during iterations.
C 
c  NLIST -- list from which to read amplitudes:
C            0; T amplitudes, 100 -- lambda, 433, R and L
C
c  SING1 -- .TRUE. if singles are included in calculation
C  TYPE  -- Character*1 variable (T, L, R..)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
#include <maxbasfn.par>
C
      LOGICAL SING1, ACTIVE_AMPS_EXSIST, ACTISFULL
      INTEGER DIRPRD, POP, VRT, T1LIST, T2LIST
      CHARACTER*1 TYPE
      CHARACTER*80 FNAME
      CHARACTER*4 ACT, MODE
      CHARACTER*10 GET_FILE
C
      DIMENSION ICORE(MAXCOR), NDRPOP(8), NDRVRT(8)
      COMMON /SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMINF/ NSTART, NIRREP, IRREPS(255,2), DIRPRD(8,8)
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /IRREPS_OF_ORB/ IRREP_ORB(MAXBASFN, 2)
C
      IUNIT = 99
      ACTISFULL = .FALSE.
      IF (TYPE .EQ. "T") GET_FILE = "ACT_AMPS_T"
      IF (TYPE .EQ. "L") GET_FILE = "ACT_AMPS_L"
      IF (TYPE .EQ. "R") GET_FILE = "ACT_AMPS_R"
      IF (ACT .EQ. "SAVE") THEN
         CALL GFNAME(GET_FILE, FNAME, ILENGTH)
         INQUIRE (FILE=FNAME(1:ILENGTH), EXIST=ACTIVE_AMPS_EXSIST)
         IF (ACTIVE_AMPS_EXSIST) THEN
             IF (TYPE .EQ. "T") istat = ishell("rm -rf ACT_AMPS_T")
             IF (TYPE .EQ. "L") istat = ishell("rm -rf ACT_AMPS_L")
             IF (TYPE .EQ. "R") istat = ishell("rm -rf ACT_AMPS_R")
             OPEN (UNIT=IUNIT, FILE=GET_FILE, FORM="FORMATTED",
     &             STATUS="NEW")
         ELSE
            OPEN (UNIT=IUNIT, FILE=GET_FILE, FORM="FORMATTED",
     &             STATUS="NEW")
         ENDIF
 
      ELSE IF (ACT .EQ. "READ") THEN

         CALL GFNAME(GET_FILE, FNAME, ILENGTH)
         INQUIRE (FILE=FNAME(1:ILENGTH), EXIST=ACTIVE_AMPS_EXSIST)
         IF (ACTIVE_AMPS_EXSIST) THEN
            CLOSE(IUNIT)
            OPEN (UNIT=IUNIT, FILE=GET_FILE, FORM="FORMATTED",
     &            STATUS="OLD")
         ELSE
           WRITE(6, "(2x,a,10a,a)") "Active space amplitude file ", 
     &                         GET_FILE ," is not found: standard",
     &                         " CC calculation "
           Write(6,"(2x,a)") "is performed" 
           RETURN
         ENDIF
      ELSE IF (ACT .EQ. "ZERO") THEN
         Write(6, "(2X, a)") "Zeroing selected doubles Amplitudes."
      ENDIF
C
      CALL GETREC(-1,'JOBARC','NDROPPOP',NIRREP,NDRPOP)
      CALL GETREC(-1,'JOBARC','NDROPVRT',NIRREP,NDRVRT)
      CALL GETREC(-1,'JOBARC','NDROTPOP',1,NDRPOPT)
      CALL GETREC(-1,'JOBARC','NDROTVRT',1,NDRVRTT)
C

      DO ISPIN = 1, 2
         IOFF = 1
         DO IRREP = 1, NIRREP
            NOCC = POP(IRREP, ISPIN)
            DO IOCC=1, NOCC
               IRREP_ORB(IOFF, ISPIN) = IRREP
               IOFF = IOFF + 1
            ENDDO
         ENDDO
C
         DO IRREP = 1, NIRREP
            NVRT = VRT(IRREP, ISPIN)
            DO IVRT =1, NVRT
               IRREP_ORB(IOFF, ISPIN) = IRREP
               IOFF = IOFF + 1
            ENDDO
        ENDDO
C
      ENDDO
C
      IF (MODE .EQ. "AMPS") THEN
          IF (TYPE .EQ. "T") THEN
             T1LIST = 90+NLIST
             T2LIST = 43+NLIST
          ELSEIF (TYPE .EQ. "L") THEN
             T1LIST = 190+NLIST
             T2LIST = 143+NLIST
          ELSEIF (TYPE .EQ. "R") THEN
             T1LIST = 490+NLIST
             T2LIST = 443+NLIST
          ENDIF 
      ELSE IF (MODE .EQ. "RESD") THEN
          T1LIST = 90+NLIST
          T2LIST = 60+NLIST
      ENDIF
C          
      IF (IUHF.NE.0) THEN
C
         DO I=1,2
            IF(SING1 .AND. ACT .EQ. "READ") CALL READ_ACT1AMPS(IUNIT,
     &                                      T1LIST, I, MODE, ACTISFULL)

            IF(SING1 .AND. .NOT. ACTISFULL) CALL SETUP_MODF_SINGLES
     &                                           (ICORE,MAXCOR,I,
     &                                            T1LIST,TYPE,
     &                                            ACT,IUNIT,MODE)
C
            IF (ACT .EQ. "READ") CALL READ_ACT2AMPS(IUNIT, T2LIST+I,
     &                                I, ACT, ACTISFULL)
C
            IF (.NOT. ACTISFULL) CALL SETUP_MODF_DOUBLES(ICORE,MAXCOR,
     &                                I,T2LIST+I,TYPE,ACT,IUNIT)
         ENDDO 
C
         I=3
         IF (ACT .EQ. "READ") CALL READ_ACT2AMPS(IUNIT, T2LIST+I, I,
     &                             ACT, ACTISFULL)
         IF  (.NOT. ACTISFULL) CALL SETUP_MODF_DOUBLES(ICORE,MAXCOR,I,
     &                              T2LIST+I,TYPE, ACT,IUNIT)
      ELSE
C
         I = 1
         IF(SING1 .AND. ACT .EQ. "READ") CALL READ_ACT1AMPS(IUNIT,
     &                                   T1LIST, I, MODE, ACTISFULL)
C
         IF(SING1 .AND. .NOT. ACTISFULL) CALL SETUP_MODF_SINGLES
     &                                        (ICORE,MAXCOR,1,
     &                                        T1LIST,TYPE,ACT,IUNIT,
     &                                        MODE)
C
         IF (ACT .EQ. "READ") CALL READ_ACT2AMPS(IUNIT, T2LIST+I,
     &                             I, ACT, ACTISFULL)
C
         IF (.NOT. ACTISFULL) CALL SETUP_MODF_DOUBLES(ICORE,MAXCOR,1,
     &                             T2LIST+1,TYPE,ACT,IUNIT)
C
         I = 3
         IF (ACT .EQ. "READ") CALL READ_ACT2AMPS(IUNIT, T2LIST+I,
     &                             I, ACT, ACTISFULL)
          IF (.NOT. ACTISFULL) CALL SETUP_MODF_DOUBLES(ICORE,MAXCOR,3,
     &                              T2LIST+3,TYPE, ACT,IUNIT)
      ENDIF
C
      CLOSE (IUNIT,STATUS="KEEP")
      RETURN
      END

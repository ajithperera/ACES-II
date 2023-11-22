
      SUBROUTINE TDEE_HBAR_MULT(ICORE,MAXCOR,IUHF,IRREPX,ISIDE)
C
C THIS ROUTINE DRIVES THE CALCULATION OF THE MATRIX-VECTOR
C PRODUCT BETWEEN A SINGLES AND DOUBLES HAMILTONIAN MATRIX 
C AND A VECTOR.
C
CEND
      IMPLICIT INTEGER (A-Z)
      LOGICAL DOUBLE, SS, SD, DS, DD
      LOGICAL Source, Target
      DOUBLE PRECISION TIN, TOUT, TIMDUM, TDAVID, TMULT
      DIMENSION ICORE(MAXCOR)
      DIMENSION LISTIJKL(2),LISTIJKA(2),LISTAIBJ(2),LISTABCI(2)
      DIMENSION LISTABCD(2)
      DIMENSION LISTR1(2),LISTR2(2),LISTZ0(2),LISTT2RS(2)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/EXTRAP/MAXEXP,NREDUCE,NTOL,NSIZEC
      COMMON/TIMSUB/TDAVID, TMULT
      COMMON/DRVHBAR/SS, SD, DS, DD
      COMMON /TIMEINFO/ TIMEIN, TIMENOW, TIMETOT, TIMENEW
C
#include "parcc.com"
C
      DATA LISTABCD /231, 231/
      DATA LISTABCI /127, 27/
      DATA LISTIJKA /107, 7/
      DATA LISTIJKL /51,  51/
      DATA LISTAIBJ /54,  54/

      DATA LISTT2RS /350,360/
      DATA LISTR1   /390,392/
      DATA LISTR2   /314,324/
      DATA LISTZ0   /334,344/
      DATA LSTT2RNG /304/

#ifdef _DEBUG_LVL0
      write(6,"(a)") "Entered tdee_hbar_mult"
      Write(6,"(a,3i10)") "The memory left: ",Maxcor
      Call check_source_lists(icore,maxcor,Irrepx,Iuhf,Iside)
#endif
      Source = .False.
      Target = .True.
      Call Tdee_init_lists(icore,maxcor,Iuhf,Irrepx,Source,
     +                     Target,Iside)
      CALL TIMER(1)
C
C CALCULATE H x C1 -> C1
C
      IF (SS) THEN

        If (Ispar) Then 
           Call RESTORE_CC_WMBEJ(ICORE,MAXCOR,IUHF)
           if (iuhf .eq. 0) CALL TDEE_MAKESS(ICORE,MAXCOR,IUHF)
        Endif 

        CALL TDEE_DR1INR1(ICORE,MAXCOR,IUHF,IRREPX, 
     &                  IUHF.EQ.0,LISTAIBJ(1),
     &                   LISTR1(ISIDE),ISIDE)
        If (Ispar) Then 
           Call RESTORE_PDCC_WMBEJ(ICORE,MAXCOR,IUHF)
           if (iuhf .eq. 0) CALL TDEE_MAKESS(ICORE,MAXCOR,IUHF)
        Endif 

        CALL TDEE_DFR1INR1(ICORE,MAXCOR,IUHF,IRREPX,LISTR1(ISIDE),
     &                     ISIDE)

      ENDIF
C
      DOUBLE = .TRUE.
  
      IF (DOUBLE) THEN

C CALCULATE H x C2 -> C1
C
        IF (DS) THEN

         CALL TDEE_DR2INR1A(ICORE,MAXCOR,IUHF,IRREPX,ISIDE,
     &                      LISTABCI(3-ISIDE),LISTR2(ISIDE))

         CALL TDEE_DR2INR1B(ICORE,MAXCOR,IUHF,IRREPX,ISIDE,
     &                      LISTIJKA(3-ISIDE),LISTR2(ISIDE))

         CALL TDEE_RESORT(ICORE,MAXCOR,IUHF,IRREPX,LISTR2(ISIDE),
     &                    LSTT2RNG)

         IF (ISIDE.EQ.1) THEN
            CALL TDEE_DFR2INR1(ICORE,MAXCOR,IUHF,1,IRREPX)
         ENDIF
C
        ENDIF

#ifdef _DEBUG_LVL0
C      Call check_source_lists(icore,maxcor,Irrepx,Iuhf,Iside)
      Call check_target_lists(icore,maxcor,Irrepx,Iuhf,Iside)
#endif
C
C CALCULATE H x C1 -> C2
C
        IF (SD) THEN

          CALL TDEE_DR1INR2A(ICORE,MAXCOR,IUHF,IRREPX,
     &                       LISTIJKA(ISIDE),LISTZ0(ISIDE),
     &                       LISTR1(ISIDE))

          CALL TDEE_DR1INR2B(ICORE,MAXCOR,IUHF,IRREPX,
     &                       LISTABCI(ISIDE),LISTZ0(ISIDE),
     &                       LISTR1(ISIDE))
        ENDIF
#ifdef _DEBUG_LVL0
C      Call check_source_lists(icore,maxcor,Irrepx,Iuhf,Iside)
      Call check_target_lists(icore,maxcor,Irrepx,Iuhf,Iside)

#endif
C
C CALCULATE H x C2 -> C2
C
        IF (DD) THEN
C
          CALL TDEE_DR2INR2(ICORE,MAXCOR,IUHF,IRREPX,ISIDE,LISTIJKL,
     &                      LISTABCD,LISTAIBJ,LISTR2(ISIDE),
     &                      LISTZ0(ISIDE),LISTT2RS,LSTT2RNG)
C
        ELSE
C
          CALL TDEE_RSDSNSD(ICORE,MAXCOR,IUHF,IRREPX,ISIDE,
     &                      LISTZ0(ISIDE),LISTT2RS)
          CALL HDIAGXDD(ICORE,MAXCOR/IINTFP,IUHF,IRREPX,NSIZEC)
C 
        ENDIF
C
      ENDIF 
C
      CALL TIMER(1)      
      TMULT = TMULT + TIMENEW

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a)") "At Exit"
      Call check_target_lists(icore,maxcor,Irrepx,Iuhf,Iside)
#endif
C
      RETURN
      END

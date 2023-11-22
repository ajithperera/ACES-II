C
      SUBROUTINE PRTLIST(SCR, NSIZE, ISPIN, LIST, IRREPX, IUNIT, 
     &                   T1T2FLG, DIRPRD, IRPDPD, ISYTYP, NIRREP)
C
C This utility routine prints the content of a list into an output.  
C  
C Arguments
C
C  SCR     : Scratch area to store the T1 or T2 amplitudes.      
C  SIZE    : Size of the scratch array.
C  ISPIN   : Spin case of the list. 
C            ISPIN = 1 [AAAA]
C            ISPIN = 2 [BBBB]
C            ISPIN = 3 [ABAB] Extension to other spin cases is stright
C            forward.
C  LIST    : List number to be printed.
C  IRREPX  : The overall symmetry of the list.
C  IUNIT   : Unit number of the out-put file
C  T1T2FLG : Flag to identify whether we are printing a T1 or T2 list.
C            .TRUE.  T1 lsit
C            .FALSE. T2 list
C  DIRPRD  : Dirrect product tabel. To be passed from 
C            the calling routine  
C  IRPDPD  : Sizes of the left and right sides of the list.
C          : To be passed from the calling routine.
C   ISYTYP : Spin type of the left and right hand size of the 
C            list. To be passed from the calling program.
C   NIRREP : The total number of irreps in the particular
C            point group.
C
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION SCR
      LOGICAL T1T2FLG
C    
      DIMENSION SCR(NSIZE), DIRPRD(8, 8), IRPDPD(8, 22), ISYTYP(2, 500)
C
      IF (T1T2FLG) THEN
C         
C We are trying to print a T1 list
C
         IF (ISPIN .EQ. 1) THEN
C
            CALL HEADER('T1(AA) Amplitudes', 0, IUNIT)
            LENT1 = IRPDPD(IRREPX, 9)
            CALL GETLST(SCR, 1, 1, 1, ISPIN, LIST)
            CALL TAB(IUNIT, SCR, LENT1, 1, LENT1, 1)
C
         ELSE IF (ISPIN .EQ. 2) THEN
C
            CALL HEADER('T1(BB) Amplitudes', 0, IUNIT)
            LENT1 = IRPDPD(IRREPX, 10)
            CALL GETLST(SCR, 1, 1, 1, ISPIN, LIST)
            CALL TAB(IUNIT, SCR, LENT1, 1, LENT1, 1)
C
         ENDIF
C
      ELSE
C
C We are trying to print a T2 list.
C
         IF (ISPIN .EQ. 1) THEN
            CALL HEADER('T2(AAAA) Amplitudes for each irrep', 0, IUNIT)
         ELSE IF (ISPIN .EQ. 2) THEN
            CALL HEADER('T2(BBBB) Amplitudes for each irrep', 0, IUNIT)
         ELSE IF (ISPIN .EQ. 3) THEN
            CALL HEADER('T2(ABAB) Amplitudes for each irrep', 0, IUNIT)
         ENDIF
C
         DO 10 IRREP = 1, NIRREP
C     
            IRREPR = IRREP
            IRREPL = DIRPRD(IRREPR, IRREPX)
C 
            DISSYT2 = IRPDPD(IRREPL, ISYTYP(1, LIST))
            NUMSYT2 = IRPDPD(IRREPR, ISYTYP(2, LIST))
C
            CALL GETLST(SCR, 1, NUMSYT2, 1, IRREPR, LIST)
C
            CALL TAB(IUNIT, SCR, DISSYT2, NUMSYT2, DISSYT2,
     &               NUMSYT2)
 10      CONTINUE
C
      ENDIF
C
C Everything is done. Return to the calling program.
C
      RETURN
      END

C
      SUBROUTINE PRTEXTLST(SCR, IRREPX, LEN, IUHF, IFCOMPR, IUNIT,
     &                     DIRPRD, IRPDPD, NIRREP)
C
C This utility routine is used to print out a list which has 
C the folowing structure. 
C
C  T1(AA) - T1(BB) - T2(ABAB) - T2(BBBB) - T2(ABAB) [UHF]
C  T1(AA) - T2(AB)                                  [RHF]
C
C The T1 and T2 vectors are printed seperatley.
C
C Arguments:
C
C  SCR     : Input array which keeps the full vector.
C  IRREPX  : Overall symmetry of the list.
C  IUHF    : UHF/RHF flag
C  IFCOMPR : Flag to indicate whether the T2 lists are compressed.
C             .TRUE.  - If symmetry compressed
C             .FALSE. - If not symmetry compressed
C  IUNIT   : Unit number of the output file
C  DIRPRD  : Dirrect product tabel. To be passed from 
C            the calling routine  
C  IRPDPD  : Sizes of the left and right sides of the list.
C          : To be passed from the calling routine.
C  NIRREP  : The number of irreps in the particular point group.
C
      IMPLICIT INTEGER (A-Z)
C
      LOGICAL IFCOMPR
      CHARACTER*40 STRNG
      DOUBLE PRECISION SCR
      DIMENSION DIRPRD(8,8),IRPDPD(8,22),SCR(LEN)
C
      I000 = 1
C
C First let's determine the sizes of individual lengths of T1 and
C T2 pieces of full vector
C
      IF (IFCOMPR) THEN
         SYTYPAAL = 1
         SYTYPAAR = 3
         SYTYPBBL = 2
         SYTYPBBR = 4
         SYTYPABL = 13
         SYTYPABR = 14
      ELSE
         SYTYPAAL = 19
         SYTYPAAR = 21
         SYTYPBBL = 20
         SYTYPBBR = 22
         SYTYPABL = 13
         SYTYPABR = 14
      ENDIF
C
      IF (IUHF .NE. 0) THEN
C
         LENT1AA = IRPDPD(IRREPX, 9)
         LENT1BB = IRPDPD(IRREPX, 10)
         LENT2AA = IDSYMSZ(IRREPX, SYTYPAAL, SYTYPAAR)
         LENT2BB = IDSYMSZ(IRREPX, SYTYPBBL, SYTYPBBR)
         LENT2AB = IDSYMSZ(IRREPX, SYTYPABL, SYTYPABR)
C
      ELSE
C
         LENT1AA = IRPDPD(IRREPX, 9)
         LENT2AB = IDSYMSZ(IRREPX, SYTYPABL, SYTYPABR)
C
      ENDIF         
C
C First write the T1 vectors into the output file.
C
      DO 10 ISPIN = 1, IUHF + 1
C
         IF (ISPIN .EQ. 1) THEN
            LENT1 = LENT1AA
            STRNG = 'Perturb T1(AA) amplitudes'
         ELSE
            LENT1 = LENT1BB
            STRNG = 'Perturb T1(BB) amplitudes'
         ENDIF
C
         IOFF = I000 + (ISPIN - 1)*LENT1AA
C
         CALL HEADER(STRNG, 0, IUNIT)
         CALL TAB(IUNIT, SCR(IOFF), LENT1, 1, LENT1, 1)
C
 10   CONTINUE
C     
C First do the ABAB spin case (Both UHF and RHF cases).
C
      IF (IUHF .NE. 0) THEN
         LENT1 = (LENT1AA + LENT1BB) 
      ELSE 
         LENT1 =  LENT1AA 
      ENDIF
C
      IOFFT2 = LENT1 + 1
C
      STRNG = 'Perturb T2(ABAB) amplitudes per irrep'
C
      CALL HEADER(STRNG, 0, IUNIT)
C     
      DO 60 IRREP = 1, NIRREP
C
         IRREPR = IRREP
         IRREPL = DIRPRD(IRREPR, IRREPX)
C     
         DISSYT2 = IRPDPD(IRREPL, 13)
         NUMSYT2 = IRPDPD(IRREPR, 14)
C
         CALL TAB(IUNIT, SCR(IOFFT2), DISSYT2, NUMSYT2, DISSYT2, 
     &              NUMSYT2)
C
         IOFFT2 = IOFFT2 + DISSYT2*NUMSYT2
C     
 60   CONTINUE
C
C Now, write the T2 vectors into the output file for AAAA and BBBB
C spin cases.
C
      LENT2 = LENT1 + LENT2AB
C
      IF (IUHF .NE. 0) THEN
C
         DO 40 ISPIN = 2, IUHF, -1
C         
            IOFFT2 =  LENT2 - (ISPIN - 2)*LENT2BB + 1
C
            IF (ISPIN .EQ. 1) THEN
               STRNG = 'Perturb T2(AAAA) amplitudes per irrep'
            ELSE
               STRNG = 'Perturb T2(BBBB) amplitudes per irrep'
            ENDIF
C     
            CALL HEADER(STRNG, 0, IUNIT)
C
            DO 50 IRREP = 1, NIRREP
C     
               IRREPR = IRREP
               IRREPL = DIRPRD(IRREPR, IRREPX)
C     
               DISSYT2 = IRPDPD(IRREPL, ISPIN)
               NUMSYT2 = IRPDPD(IRREPR, 2 + ISPIN)
C
               CALL TAB(IUNIT, SCR(IOFFT2), DISSYT2, NUMSYT2, DISSYT2, 
     &                  NUMSYT2)
C
               IOFFT2 = IOFFT2 + DISSYT2*NUMSYT2
C     
 50         CONTINUE
C   
 40      CONTINUE
C
      ENDIF
C      
      RETURN
      END

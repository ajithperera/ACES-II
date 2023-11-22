C
      SUBROUTINE QRESORT(SCR, WOUT, MAXCOR, POP1, POP2, POP3, POP4,
     &                   IRREPX, NSIZEOUT, SPCASE, LISTFROM, LHSEXP, 
     &                   LFTTYP1, SYTYPL, SYTYPR)
C
C This routine drives the resort of a list on a disk. Adapted from 
C the ACESII library routines and modified to be used in the calculation
C of quadratic contributions to EOM-CCSD properties. The original
C list is labeled as PQRS or 1234.
C
C  SCR     - Working scratch area used in this routine
C  WOUT    - Resorted list; Output argument. 
C  MAXCOR  - Total core memory available in double words
C  POP1    - Population vector for P
C  POP2    - Population vector for Q
C  POP3    - Population vector for R
C  POP4    - Population vector for S
C  IRREPX  - The overall symetry of the list (Usually 1)
C  LISTIN  - The list holding the PQRS quantity
C  LHSEXP  - Logical flag set to .TRUE. if left-hand side needs 
C            to be expanded before resorting.
C  LFTTYP1 - Left hand symmetry type of the PQRS list 
C            [Used only if LHSEXP = .TRUE.]
C  SYTYPL - Left hand symmetry type of the resorted list.
C  SYTYPR - Right hand symmetry type of the resorted list
C
      IMPLICIT INTEGER (A-Z)
      CHARACTER*4 SPCASE
      LOGICAL LHSEXP
      DOUBLE PRECISION SCR, WOUT
C
      DIMENSION SCR(MAXCOR), WOUT(NSIZEOUT)
      DIMENSION POP1(8),POP2(8),POP3(8),POP4(8)
C
      COMMON /SYMINF/ NSTART,NIRREP,IRRVEC(255,2),DIRPRD(8,8)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /FLAGS/ IFLAGS(100)
      COMMON/ FILES/ LUOUT, MOINTS
C
      PRTLEVL = IFLAGS(1)
C
C Calculate size of input lists
C
      
      IF (.NOT. LHSEXP) THEN

         NSIZEIN = IDSYMSZ(IRREPX, ISYTYP(1, LISTFROM), 
     &             ISYTYP(2, LISTFROM))
C
         MAXDSZ  = -1
C
         DO 10 I = 1, NIRREP
C
            DSZIN = IRPDPD(I, ISYTYP(1, LISTFROM))
            DSZOUT  = IRPDPD(I, SYTYPL)
C
            INTDSZ  = MAX(DSZIN, DSZOUT)  
            MAXDSZ  = MAX(MAXDSZ, INTDSZ)
C     
 10      CONTINUE
C
C See if there is sufficient core to do this the easy way,
C If there is sufficient core memory allocate the memory and load
C the full vector in to the memory.
C
         ITOP = 0
         ITOP = NSIZEIN + MAXDSZ
C
         IF (ITOP .LE. MAXCOR) THEN
C     
            IF (IFLAGS(1) .GE. 20) THEN
               WRITE(LUOUT, *)
               WRITE(LUOUT, 1000) LISTFROM
            ENDIF
C     
            I000 = 1
            I010 = I000 + NSIZEIN
            I020 = I010 + MAXDSZ
C
C Now resort the list using an in-core algorithm. 
C SPCASE determines the structure of the resorted list and it
C is an input argument to the subroutine.
C
            CALL GETALL(SCR(I000), NSIZEIN, IRREPX, LISTFROM)
            CALL SSTGEN(SCR(I000), WOUT, NSIZEIN, POP1, POP2, POP3,
     &                  POP4, SCR(I010), IRREPX, SPCASE)
         ELSE
C
            I000 = 1
            I010 = I000 + MAXDSZ
            I020 = I010 + MAXDSZ
C
            IF (I020 .GT. MAXCOR) THEN
               CALL INSMEM('QRESORT', I020, MAXCOR)
            ELSE
C     
               IF (IFLAGS(1) .GE. 20) THEN
                  WRITE(LUOUT, *)
                  WRITE(LUOUT, 2000) LISTFROM
               ENDIF
C
C Actual reordering is carried out in QSSTDSK. The QSSTDSK is 
C corrected version of ACESII library routine SSTDSK. The 
C ACES2 library routine is a disaster.
C
               CALL QSSTDSK(SCR(I000), WOUT, NSIZEIN, POP1, POP2, POP3,
     &                      POP4, SCR(I010), IRREPX, SPCASE, LISTFROM,
     &                      LISTOUT, SYTYPL, SYTYPR, LHSEXP, LFTTYP1)
            ENDIF
C
         ENDIF
C
      ELSE 
C
C Calculate the size of input lists.
C
         NSIZEIN = IDSYMSZ(IRREPX, LFTTYP1, ISYTYP(2, LISTFROM))
C
         MAXNEED = -1
         MAXDSZ  = -1
C
         DO 20 I = 1, NIRREP
C     
            DSZIN    = IRPDPD(I, LFTTYP1)
            NUMDISIN = IRPDPD(I, ISYTYP(2, LISTFROM))
C         
            MAXDSZ  = MAX(MAXDSZ, DSZIN)
            MAXNEED = MAX(MAXNEED, DSZIN*NUMDISIN)
C     
 20      CONTINUE
C
C See if there is sufficient core to do this the easy way,
C If there is sufficient core memory allocate the memory and load
C the full vector in to the memory.
C
         ITOP = 0
         ITOP = NSIZEIN + MAXNEED
C
         IF (ITOP .LT. MAXCOR) THEN
C     
            IF (IFLAGS(1) .GE. 20) THEN
               WRITE(LUOUT, *)
               WRITE(LUOUT, 1000) LISTFROM
            ENDIF
C     
            IOFF = 1
            I000 = 1
            I010 = I000 + NSIZEIN
C
C Now resort the list using an in-core algorithm. 
C SPCASE determines the structure of the resorted list and it
C is an input argument to the subroutine.
C
            DO 30 IRREPWR = 1, NIRREP
C
               IRREPWL = DIRPRD(IRREPWR, IRREPX)
C     
               DISSYW = IRPDPD(IRREPWL, ISYTYP(1, LISTFROM))
               NUMSYW = IRPDPD(IRREPWR, ISYTYP(2, LISTFROM))
               DISEXP = IRPDPD(IRREPWL, LFTTYP1)
C
               I020 = I010 + DISEXP*NUMSYW
C
               CALL GETLST(SCR(I010), 1, NUMSYW, 2, IRREPWR, LISTFROM)
C     
C Expand the left side and accumulate on the  SCR(I000) array which
C has to reordered using SSTGEN.
C
               CALL SYMEXP2(IRREPWL, POP1, DISEXP, DISSYW, NUMSYW,
     &                      SCR(I010), SCR(I010))
               CALL SCOPY(DISEXP*NUMSYW, SCR(I010), 1, SCR(IOFF), 1)
C
               IOFF = IOFF + DISEXP*NUMSYW
C
 30         CONTINUE
C
C Now we have the full expanded list in SCR(I000). Resort the
C by calling SSTGEN.
C
            I020 = I010 + MAXDSZ
C     
            CALL SSTGEN(SCR(I000), WOUT, NSIZEIN, POP1, POP2, POP3,
     &                  POP4, SCR(I010), IRREPX, SPCASE)
C     
         ELSE
C     
            I000 = 1
            I010 = I000 + MAXDSZ
            I020 = I010 + MAXDSZ
C
            IF (I020 .GT. MAXCOR) THEN
               CALL INSMEM('QRESORT', I020, MAXCOR)
            ELSE
C     
               IF (IFLAGS(1) .GE. 20) THEN
                  WRITE(LUOUT, *)
                  WRITE(LUOUT, 2000) LISTFROM
               ENDIF
C     
C Actual reordering is carried out in QSSTDSK. The QSSTDSK is 
C corrected version of ACESII library routine SSTDSK. The 
C ACES2 library routine is a disaster.
C     
               CALL QSSTDSK(SCR(I000), WOUT, NSIZEIN, POP1, POP2, POP3,
     &                      POP4, SCR(I010), IRREPX, SPCASE, LISTFROM,
     &                      LISTOUT, SYTYPL, SYTYPR, LHSEXP, LFTTYP1)
            ENDIF
C     
         ENDIF
C
      ENDIF
C
1000  FORMAT(T3,'@QRESORT-I, Reordering list',I3,' using in-core ',
     &          'algorithm.')
2000  FORMAT(T3,'@QRESORT-I, Reordering list',I3,' using out-of-core '
     &          ,'algorithm.')
C
      RETURN
      END

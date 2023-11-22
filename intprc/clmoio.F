C
C THIS ROUTINE CALCULATES THE ADDRESSES OF THE VARIOUS LISTS
C  IN THE MOINTS FILE.  
C
C THIS ROUTINE IS CALLED IN INTPRC TO INITIALIZE THE MOINTS FILE
C AND TO CREATE ALL THE INTEGRAL AND AMPLITUDE LISTS ON MOINTS
C (AT LEAST THE STANDARD LISTS)
C
#include "flags.h" /* for h_IFLAGS_* and h_IFLAGS2_* */
C 
      SUBROUTINE CLMOIO(ICORE,MAXCOR,NOCCO,NVRTO,IUHF)
      IMPLICIT INTEGER(A-H,O-Z)
      LOGICAL DOALL,ABIJ,NOABCD,NO_REDUNDANT,SYMMETRIC,RCCD
      LOGICAL DRCCD
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /FLAGS2/ IFLAGS2(500)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /SYMINF/ NSTART,NIRREP,IRREPY(255,2),DIRPRD(8,8)
      COMMON /CONTROL/ IPRNT,IXXX,IXXX2
      COMMON /DOINTS/ DOALL,ABIJ,NOABCD
      COMMON /SHIFT/  ISHIFT,NDRGEO

C Parameters for parametrized CC
#include "parcc.com"

      DIMENSION NOCCO(2),NVRTO(2),ICORE(MAXCOR)
      IX1(I)=(I*(I-1))/2
      IX2(I)=(I*(I+1))/2
      NOCCA=NOCCO(1)
      NOCCB=NOCCO(2)
      NVRTA=NVRTO(1)
      NVRTB=NVRTO(2)
      NINT=0
      NAMP=0
      NO_REDUNDANT = iflags2(h_IFLAGS2_single_store).ne.0 
      DRCCD        = (IFLAGS(h_IFLAGS_calclevel) .EQ. 49)
c      IMODE=1
c      IMODE3=1
      call aces_io_reset
      IMODE=0
      IMODE3=0
      NFULL=NOCCA*NOCCB*NVRTA*NVRTB
C***********************************************************************
C PHHH INTEGRALS.
C
      IF(.NOT.ABIJ)THEN
c      CALL INIPCK(1,3,16,7,IMODE,0,1)
       CALL INIPCK(1,3,16,7+ISHIFT,IMODE,0,1)
       IF(IUHF.NE.0)THEN
c       CALL INIPCK(1,4,17,8,IMODE,0,1)
        CALL INIPCK(1,4,17,8+ISHIFT,IMODE,0,1)
c       CALL INIPCK(1,14,11,9,IMODE,0,1)
        CALL INIPCK(1,14,11,9+ISHIFT,IMODE,0,1)
       ENDIF
c      CALL INIPCK(1,14,18,10,IMODE,0,1)
       CALL INIPCK(1,14,18,10+ISHIFT,IMODE,0,1)
C***********************************************************************
C PHPH INTEGRALS.
C
c      CALL INIPCK(1,9,9,23,IMODE,0,1)
       CALL INIPCK(1,9,9,23+ISHIFT,IMODE,0,1)
C
C This 123 and 124 lists store <IA||JB> (BI,AJ ordered) store the
C Coulomb only contributions for DRPA (perhaps parametrized CC) etc.
C Notice that this list is resorted in VDENS. 

       IF ((ISPAR .AND. COULOMB) .OR. DRCCD) CALL INIPCK(1,9,9,
     &                                            123,IMODE,0,1)

c      IF(IUHF.NE.0)CALL INIPCK(1,10,10,24,IMODE,0,1)
       IF(IUHF.NE.0)CALL INIPCK(1,10,10,24+ISHIFT,IMODE,0,1)
      
       IF (IUHF. NE.0 .AND. ((ISPAR .AND. COULOMB) .OR. DRCCD)) 
     &     CALL INIPCK(1,10,10,124,IMODE,0,1)

c      CALL INIPCK(1,11,11,25,IMODE,0,1)
       CALL INIPCK(1,11,11,25+ISHIFT,IMODE,0,1)
c      IF(IUHF.NE.0)CALL INIPCK(1,12,12,26,IMODE,0,1)
       IF(IUHF.NE.0)CALL INIPCK(1,12,12,26+ISHIFT,IMODE,0,1)
C***********************************************************************
C PPPH INTEGRALS.
C
       IF(IUHF.NE.0)THEN
c       CALL INIPCK(1,1,9,27,IMODE,0,1)
        CALL INIPCK(1,1,9,27+ISHIFT,IMODE,0,1)
c       CALL INIPCK(1,2,10,28,IMODE,0,1)
        CALL INIPCK(1,2,10,28+ISHIFT,IMODE,0,1)
c       CALL INIPCK(1,13,18,29,IMODE,0,1)
        CALL INIPCK(1,13,18,29+ISHIFT,IMODE,0,1)
       ENDIF
c      CALL INIPCK(1,13,11,30,IMODE,0,1)
       CALL INIPCK(1,13,11,30+ISHIFT,IMODE,0,1)
      ENDIF
C**********************************************************************
C PPPP INTEGRALS.  AA TYPE ARE STORED SYMMETRY PACKED SO WE NEED
C  TO GET THE APPROPRIATE SYMMETRY INFORMATION FIRST.  THESE AND
C  THE HHHH INTEGRALS ARE DONE ONLY IF THE CALCLEVEL IS HIGHER THAN
C  MBPT(2).
C
      IF(.NOT.NOABCD)THEN
c       IF(IUHF.NE.0)THEN
c        CALL INIPCK(1,1,1,231,IMODE3,0,1)
c        CALL INIPCK(1,2,2,232,IMODE3,0,1)
c        CALL INIPCK(1,13,13,233,IMODE3,0,1)
c       ELSE
C
C  CURRENTLY WE ARE USING FOR GRADIENT RUNS A DIFFERENT THAN IN
C  ENERGY CALCULATIONS
C
c        IF(IFLAGS(3).NE.0.OR.IFLAGS(18).NE.0.OR.
c     &     IFLAGS(19).NE.0) THEN
c         CALL INIPCK(1,13,13,233,IMODE3,0,1)
c        ELSE 
c         CALL INIPCK(1,5,13,233,IMODE3,0,1)
c        ENDIF
c       ENDIF
      ENDIF
C**********************************************************************
C HHHH INTEGRALS.  AA TYPE ARE STORED SYMMETRY PACKED SO WE NEED
C  TO GET THE APPROPRIATE SYMMETRY INFORMATION FIRST.
C
      IF(.NOT.ABIJ)THEN
       IF(IUHF.NE.0)THEN
c       CALL INIPCK(1,3,3,11,IMODE,0,1)
        CALL INIPCK(1,3,3,11+ISHIFT,IMODE,0,1)
c       CALL INIPCK(1,4,4,12,IMODE,0,1)
        CALL INIPCK(1,4,4,12+ISHIFT,IMODE,0,1)
       ENDIF
c      CALL INIPCK(1,14,14,13,IMODE,0,1)
       CALL INIPCK(1,14,14,13+ISHIFT,IMODE,0,1)
      ENDIF
C
C***********************************************************************
C PPHH INTEGRALS.  THESE SIT NEXT TO T AMPLITUDES ON DISK.
C
      CALL INIPCK(1,1,3,14+ISHIFT,IMODE,0,1)

C This 114 list store <AB|IJ> (AB,IJ ordered) store the
C Coulomb only contributions for parameterized CCSD, CCD etc.

            IF ((ISPAR .AND. COULOMB) .OR. DRCCD) 
CSSS     &          CALL INIPCK(1,19,21,114,IMODE,0,1)
     &          CALL INIPCK(1,19,21,214,IMODE,0,1)

      IF (IUHF.NE.0)THEN
            CALL INIPCK(1,2,4,15+ISHIFT,IMODE,0,1)

C This 115 list store <ab|ij> (ab,ij ordered) store the
C Coulomb only contributions for parameterized CCSD, CCD etc.

            IF ((ISPAR .AND. COULOMB) .OR. DRCCD)
CSSS     &          CALL INIPCK(1,20,22,115,IMODE,0,1)
     &          CALL INIPCK(1,20,22,215,IMODE,0,1)
 
            CALL INIPCK(1,10,10,6+ISHIFT,IMODE,0,1)

         IF (.NOT. NO_REDUNDANT) THEN
            CALL INIPCK(1,10,9,17+ISHIFT,IMODE,0,1)
            CALL INIPCK(1,10,10,20+ISHIFT,IMODE,0,1)

C This 120 list store <ab|ij> (ai,bj ordered) store the
C Coulomb only contributions for parameterized CCSD, CCD etc.

            IF ((ISPAR .AND. COULOMB) .OR. DRCCD) THEN
                CALL INIPCK(1,10,10,120,IMODE,0,1)
            ENDIF 
 
            CALL INIPCK(1,12,11,22+ISHIFT,IMODE,0,1)
         ELSE
            CALL INIT_NOREDNT_LIST(1,10,9,17+ISHIFT)
            CALL INIT_NOREDNT_LIST(1,10,10,20+ISHIFT)
            CALL INIT_NOREDNT_LIST(1,12,11,22+ISHIFT)
         ENDIF 
      ENDIF
C
         CALL INIPCK(1,13,14,16+ISHIFT,IMODE,0,1)
         CALL INIPCK(1,9,9,5+ISHIFT,IMODE,0,1)
      IF (.NOT. NO_REDUNDANT) THEN
         CALL INIPCK(1,9,10,18+ISHIFT,IMODE,0,1)
         CALL INIPCK(1,9,9,19+ISHIFT,IMODE,0,1)

C This 119 list store <AB|IJ> (AI,BJ ordered) store the
C Coulomb only contributions for parameterized CCSD, CCD etc.

         IF ((ISPAR .AND. COULOMB) .OR. DRCCD) THEN
            CALL INIPCK(1,9,9,119,IMODE,0,1)
         ENDIF 

         CALL INIPCK(1,11,12,21+ISHIFT,IMODE,0,1)
      ELSE
         CALL INIT_NOREDNT_LIST(1,9,10,18+ISHIFT)
         CALL INIT_NOREDNT_LIST(1,9,9,19+ISHIFT)
         CALL INIT_NOREDNT_LIST(1,11,12,21+ISHIFT)
      ENDIF
C
c   o do not recreate unshifted lists (in a dropped MO basis)
      IF (ISHIFT.NE.0) RETURN
C
C***********************************************************************
C DENOMINATOR ARRAYS
C
      SYMMETRIC = (IFLAGS(h_IFLAGS_calclevel) .EQ. 48 .OR.
     +             IFLAGS(h_IFLAGS_calclevel) .EQ. 49)
      RCCD      = (IFLAGS(h_IFLAGS_calclevel) .EQ. 48) 
      CALL INIPCK(1,13,14,50,IMODE,0,1)
      IF (SYMMETRIC) THEN
CSSS         CALL INIPCK(1,1,3,48,IMODE,0,1)
CSSS         IF(IUHF.NE.0) CALL INIPCK(1,2,4,49,IMODE,0,1)
CSSS         CALL INIPCK(1,19,21,48,IMODE,0,1)

         IF (RCCD .AND. IUHF .EQ.0) THEN
            CALL INIPCK(1,13,14,48,IMODE,0,1)
         ELSE 
            CALL INIPCK(1,19,21,48,IMODE,0,1)
         ENDIF 
         IF(IUHF.NE.0) CALL INIPCK(1,20,22,49,IMODE,0,1)
      ELSE
         CALL INIPCK(1,1,3,48,IMODE,0,1)
         IF(IUHF.NE.0) CALL INIPCK(1,2,4,49,IMODE,0,1)
      ENDIF

      IF(IFLAGS(21).NE.3)THEN
       CALL INIPCK(1,13,14,66,IMODE,0,1)
       IF (SYMMETRIC) THEN
CSSS          CALL INIPCK(1,1,3,64,IMODE,0,1)
CSSS          IF(IUHF.NE.0) CALL INIPCK(1,2,4,65,IMODE,0,1)
CSSS          CALL INIPCK(1,19,21,64,IMODE,0,1)
          IF (RCCD .AND. IUHF .EQ.0) THEN
             CALL INIPCK(1,13,14,64,IMODE,0,1)
          ELSE
             CALL INIPCK(1,19,21,64,IMODE,0,1)
          ENDIF 
          IF(IUHF.NE.0) CALL INIPCK(1,20,22,65,IMODE,0,1)
       ELSE
          CALL INIPCK(1,1,3,64,IMODE,0,1)
          IF(IUHF.NE.0)CALL INIPCK(1,2,4,65,IMODE,0,1)
       ENDIF 

      ENDIF
C
C***********************************************************************
C T2 AMPLITUDES (START ON PHYSICAL RECORD BOUNDARY)
C
      CALL INIPCK(1,13,14,46,IMODE,-1,1)
      IF (SYMMETRIC) THEN
CSSS          CALL INIPCK(1,1,3,44,IMODE,0,1)
CSSS          IF(IUHF.NE.0) CALL INIPCK(1,2,4,45,IMODE,0,1)
CSSS          CALL INIPCK(1,19,21,44,IMODE,0,1)

         IF (RCCD .AND. IUHF .EQ.0) THEN
            CALL INIPCK(1,13,14,44,IMODE,0,1)
         ELSE 
          CALL INIPCK(1,19,21,44,IMODE,0,1)
         ENDIF 
         IF(IUHF.NE.0) CALL INIPCK(1,20,22,45,IMODE,0,1)
      ELSE
        CALL INIPCK(1,1,3,44,IMODE,0,1)
        IF(IUHF.NE.0)CALL INIPCK(1,2,4,45,IMODE,0,1)
      ENDIF 
C
C ONLY NEED THE REST FOR POST-MBPT(2) STUFF.
C
C SG 4/24/96 Remodified by Ajith 03/20/97
C      IF(IFLAGS(2).GT.1.OR.IFLAGS(22).NE.0)THEN
C
      IF (IFLAGS(2).GT.1.OR.IFLAGS(22).NE.0.OR.IFLAGS2(117).GT.0) THEN
         IF (.NOT. NO_REDUNDANT) THEN 
            CALL INIPCK(1,9,9,34,IMODE,0,1)
            CALL INIPCK(1,9,10,37,IMODE,0,1)
            CALL INIPCK(1,11,12,39,IMODE,0,1)
            IF (IUHF.NE.0)THEN
               CALL INIPCK(1,10,9,36,IMODE,0,1)
               CALL INIPCK(1,12,11,38,IMODE,0,1)
               CALL INIPCK(1,10,10,35,IMODE,0,1)
            ENDIF
         ELSE
            CALL INIT_NOREDNT_LIST(1,9,9,34)
            CALL INIT_NOREDNT_LIST(1,9,10,37)
            CALL INIT_NOREDNT_LIST(1,11,12,39)
cSSS            CALL INIT_NOREDNT_LIST(1,9,9,134)
cSSS            CALL INIT_NOREDNT_LIST(1,9,10,137)
cSSS            CALL INIT_NOREDNT_LIST(1,11,12,139)
            IF (IUHF.NE.0)THEN
               CALL INIT_NOREDNT_LIST(1,10,9,36)
               CALL INIT_NOREDNT_LIST(1,12,11,38)
               CALL INIT_NOREDNT_LIST(1,10,10,35)
cSSS               CALL INIT_NOREDNT_LIST(1,10,9,136)
cSSS               CALL INIT_NOREDNT_LIST(1,12,11,138)
cSSS               CALL INIT_NOREDNT_LIST(1,10,10,135)
            ENDIF
         ENDIF
      ENDIF
C
      IF(IFLAGS(1).GT.100)CALL aces_io_summary
C
      RETURN
      END

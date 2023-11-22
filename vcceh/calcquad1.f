


      SUBROUTINE CALCQUAD1(ICORE, QUAD, MAXCOR, ISYMPRT, IOFFPRT, 
     &                     NUMPERT, IUHF, FACT, MAXPRT, IFREQ, ITOP,
     $     image_quad)
C
C This routine initiate the calculation of quadratic contribution
C to EOM-CCSD second-order properties. This assign list numbers
C required in quadratic term.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z) 
C     
      DIMENSION ICORE(MAXCOR), QUAD(MAXPRT, MAXPRT), ISYMPRT(MAXPRT),
     &          IOFFPRT(MAXPRT)
      INTEGER DIRPRD
      logical image_quad
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPY(255,2),DIRPRD(8,8)
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /FILES/ LUOUT, MOINTS
      COMMON/TIMSUB/TDAVID, TMULT, TQUAD
      common /TIMEINFO/ TIMEIN, TIMENOW, TIMETOT, TIMENEW
C
      CALL TIMER(1)
C
      I000 = 1 
      IOFFSET = 1
C
C Print the description of the program, What it does, authors etc.
C to the user.    
C
      IF (IFREQ .EQ. 1) CALL BANNER(LUOUT)
C
C Initialize the common blocks used in the quadratic term
C
      CALL INITQUAD(ICORE(I000), MAXCOR, IUHF)      
C
C Create list Hbar(MB,EJ), Hbar(Mb,Ej) and Hbar(Ab,Cd) for RHF
C calculations.
C
      IF (IUHF .EQ. 0 .AND. IFREQ .EQ. 1) THEN
         CALL GENRHFLST(ICORE, MAXCOR, IUHF)
      ENDIF
C     
C First-order perturb T amplitudes are stored in list 373 according
C to the following order,
C T1(AA)-T1(BB)-T2(AB)-T2(AA)-T2(BB) [UHF]
C T1[AA]-T2[AB]                      [RHF]
C Unwind them and put them in separate lists to make life simpler.
C
      DO 10 IALPHA = 1, NUMPERT
C
         IOFFST = IOFFPRT(IALPHA)
         IRREPX = ISYMPRT(IALPHA)
C
         IF (IRREPX .NE. 0) THEN
C     
C Kill the lists 400-500 and reinitialize them to be used in
C quadratic term. Also create the lists used in the quadratic
C term.
C
c YAU : old
c           CALL KILL400 (ICORE(I000), MAXCOR, IUHF)
c YAU : new
            CALL ACES_IO_REMOVE(54,'DERGAM')
c YAU : end
            CALL MKQUDLST(ICORE(I000), MAXCOR, IUHF, IRREPX)
C
            LISTIN = (IFREQ - 1)*2 + 373
c                  write(6,*) ' ialpha, ioffst, listin ',
c     $                 ialpha, ioffst, listin
            CALL UNWIND (ICORE(I000), MAXCOR, IRREPX, IUHF, LISTIN,
     &                   NUMPERT, IOFFST, .FALSE.)
C            
C Generate I(ALPA) and G(ALPA) quadratic intermediates. Equations 
C for the I(ALP) and G(BETA) are given in each routines called
C by GENQINT.
C
            CALL GENQINT(ICORE(I000), MAXCOR, IRREPX, IOFFSET, IUHF)
C
            DO 20 IBETA = 1, NUMPERT
C
               IRREPZ = ISYMPRT(IBETA)
               JOFFST = IOFFPRT(IBETA)
C
               IF (IRREPZ .EQ. IRREPX) THEN
C
                  LISTIN = (ITOP - IFREQ)*2 + 373
                  if (image_quad) then
                     LISTIN = (IFREQ - 1)*2 + 373
                  endif
c
c                  write(6,*) ' ibeta, joffst, listin ',
c     $                 ibeta, joffst, listin
c                  write(6,*) ' wrong list in calcquad1'
                  CALL UNWIND(ICORE(I000), MAXCOR, IRREPX, IUHF, LISTIN,
     &                        NUMPERT, JOFFST, .TRUE.)
C 
C Calculate the contribution of I(ALP) and G(ALP) quadratic
C intermediates to the singles and doubles equation (Q(ALPA,BETA)).
C The singles and doubles eqautins are given in GENQSD.
C
                  CALL GENQSD(ICORE(I000), MAXCOR, IRREPX, IOFFSET,
     &                        IUHF)
C
C Calculate quadratic contribution for every pair of perturbations  
C and store them in the array QUAD.
C
                  CALL CMPQUAD(ICORE(I000), MAXCOR, QUAD(IALPHA,
     &                         IBETA), MAXPRT, IRREPX, IOFFSET,
     &                         FACT, IUHF)
C
               ENDIF
C
 20         CONTINUE
C
         ENDIF
C
 10   CONTINUE
C
      CALL TIMER(1)      
      TQUAD = TQUAD + TIMENEW
C
      RETURN
      END

      SUBROUTINE RCC_DODIIS0(SCR,MAXCOR,IUHF,IRREPX,ICYCLE,
     &                       ICONVG,ICONTL,SING,LISTT20,LISTT21,
     &                       LISTT10,IOFF0,LISTT11,IOFF1,LSTERR,
     &                       TYPE,DAMP_PARAMETER)

C
C DRIVER FOR DIIS CONVERGENCE ACCELERATION FOR DERIVATIVE T AMPLITUDE
C EQUATIONS
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL SING
C
      DIMENSION SCR(MAXCOR)
C
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/FLAGS/IFLAGS(100)
cYAU - This is the same in DIISLST and should probably be an include file.
cjp
cjp separate data for individual references
cjp probably would be better to treat all amplitudes togerher
cjp as one large vector, but it seems to work like this nicely
#include "bwcc.com"
      integer ld_diis
      parameter (ld_diis=25)
      double precision r(ld_diis,ld_diis,maxref)
      common/diisdat1/r,idimdiis(maxref)
      save  /diisdat1/
cjp
C
      DATA ONE,ONEM/1.D0,-1.D0/
C
      INDXF(I,J,N)=I+(J-1)*N
C
cjp   initialize
      lenvec=idimdiis(iref)
      ndimr=min(icycle,iflags(12))
      idimr=1+mod(icycle-1,iflags(12))
cjp
      TOL=10.D0**(-ICONTL)

C
C DO BOOKKEEPING AND COMPUTE RESIDUAL VECTOR
C
      MAXDM1=ld_diis+2
      I000=1
      I010=I000+MAX(LENVEC,MAXDM1*MAXDM1)
      I020=I010+MAX(LENVEC,MAXDM1*MAXDM1)
      IEND=I020+MAXDM1
      IF(IEND.GE.MAXCOR) CALL INSMEM('DODIIS0',IEND,MAXCOR)
C
C Symmetrize the T2 amplitudes for RHF cases.
C
      IF (IUHF.EQ.0) THEN
         CALL SYMETRIZE_T2(SCR, MAXCOR, 1, LISTT20+2)
         CALL SYMETRIZE_T2(SCR, MAXCOR, 1, LISTT21+2)
      END IF
C
      IOFFT=I000
      IOFFHT=I010
      IF(SING)THEN
       DO 10 ISPIN=1,1+IUHF
        CALL GETLST(SCR(IOFFT) ,1,1,1,IOFF0+ISPIN,LISTT10)
        CALL GETLST(SCR(IOFFHT),1,1,1,IOFF1+ISPIN,LISTT11)
#ifdef _DEBUG_LVLM
       call checksum("T1-DIIS1",scr(iofft), IRPDPD(IRREPX,8+ISPIN))
       call checksum("T1-DIIS1",scr(ioffht),IRPDPD(IRREPX,8+ISPIN))
#endif
        IOFFT=IOFFT+IRPDPD(IRREPX,8+ISPIN)
        IOFFHT=IOFFHT+IRPDPD(IRREPX,8+ISPIN)
10     CONTINUE
      ENDIF
       DO 11 ISPIN=3,1,IUHF-2
       LISTT=LISTT20-1+ISPIN
       NSIZE=IDSYMSZ(IRREPX,ISYTYP(1,LISTT),ISYTYP(2,LISTT))
       CALL GETALL(SCR(IOFFT) ,NSIZE,IRREPX,LISTT20-1+ISPIN)
       CALL GETALL(SCR(IOFFHT),NSIZE,IRREPX,LISTT21-1+ISPIN)

#ifdef _DEBUG_LVLM
       call checksum("T2-DIIS1",scr(iofft), nsize)
       call checksum("T2-DIIS1",scr(ioffht), nsize)
#endif
       IOFFT=IOFFT+NSIZE
       IOFFHT=IOFFHT+NSIZE
11    CONTINUE
C
C SAVE OLD T AS NEW APPROXIMATE SOLUTION VECTOR ON 1,LSTERR
C
      CALL PUTLST(SCR(I000),IDIMR,1,1,1,LSTERR)
      CALL SAXPY (LENVEC,ONEM,SCR(I000),1,SCR(I010),1)

      IF (DAMP_PARAMETER .NE. 1.0D0) THEN

         CALL DSCAL(LENVEC, 1.0D0/DAMP_PARAMETER, SCR(I010), 1)
         RESMITER=SCR(I010-1+ISAMAX(LENVEC,SCR(I010),1))
         CALL DSCAL(LENVEC, DAMP_PARAMETER, SCR(I010), 1)
         
      ElSE

         RESMITER=SCR(I010-1+ISAMAX(LENVEC,SCR(I010),1))

      ENDIF 
C
C CALCULATE ERROR MATRIX UP TO THIS POINT
C
      DO I=1,IDIMR-1
         CALL GETLST(SCR(I000),I,1,1,2,LSTERR)
         X=SDOT(LENVEC,SCR(I000),1,SCR(I010),1)
         R(I,IDIMR,iref)=X
         R(IDIMR,I,iref)=X
      END DO
C
C SAVE T(NEW)-T(OLD) AS CORRESPONDING ERROR VECTOR ON 2,LSTERR
C
      CALL PUTLST(SCR(I010),IDIMR,1,1,2,LSTERR)
      X=SDOT(LENVEC,SCR(I010),1,SCR(I010),1)
      
      R(IDIMR,IDIMR,iref)=X
C
C CALCULATE ERROR MATRIX UP TO THIS POINT (cont.)
C
      DO I=IDIMR+1,NDIMR
         CALL GETLST(SCR(I000),I,1,1,2,LSTERR)
         X=SDOT(LENVEC,SCR(I000),1,SCR(I010),1)
         R(I,IDIMR,iref)=X
         R(IDIMR,I,iref)=X
      END DO
C
C CALCULAtE DIIS EXPANSION COEFFICIENTS
C
cjp
      I=ld_diis
      CALL DODIIS(R(1,1,iref),SCR(I000),SCR(I020),I,NDIMR)
C
      IF(isbwcc.or.IFLAGS(1).GE.10)THEN
       WRITE(6,1000)
       CALL PRVECR(SCR(I020),NDIMR)
      ENDIF
C
C EVALUATE DIIS RESIDUAL
C
      fact = scr(i020-1+idimr)
      do i = 0, lenvec-1
         scr(i000+i) = fact*scr(i010+i)
      end do
cYAU - fix this to load as many vectors as possible and use dgemm
      IOFF=I020-1+NDIMR
      DO I=NDIMR,IDIMR+1,-1
         CALL GETLST(SCR(I010),I,1,1,2,LSTERR)
         CALL SAXPY (LENVEC,SCR(IOFF),SCR(I010),1,SCR(I000),1)
         IOFF=IOFF-1
      END DO
      IOFF=IOFF-1
      DO I=IDIMR-1,1,-1
         CALL GETLST(SCR(I010),I,1,1,2,LSTERR)
         CALL SAXPY (LENVEC,SCR(IOFF),SCR(I010),1,SCR(I000),1)
         IOFF=IOFF-1
      END DO
      RESMDIIS=SCR(ISAMAX(LENVEC,SCR(I000),1))
C
C WRITE OUT RESIDUALS AND CHECK CONVERGENCE
C
      WRITE(6,1001)ICYCLE
      WRITE(6,1002)RESMITER
      WRITE(6,1003)RESMDIIS
      IF(ABS(RESMITER).LT.TOL)THEN
       ICONVG=0
      WRITE(6,1004)ICYCLE
C
C This is commented to avoid problems where DIIS converge
C faster than the amplitudes. Be/STO-3G, Ajith /0308/2000 
C 
C      ELSEIF(ABS(RESMDIIS).LT.TOL)THEN
C      ICONVG=0
C       WRITE(6,1004)ICYCLE
      ELSE
       ICONVG=1
      ENDIF
C SG 7/23/98
      IF (IFLAGS(22) .EQ. 1) THEN
        CALL PUTREC(20, 'JOBARC', 'T2CNVCRT', IINTFP,
     &     MIN(ABS(RESMITER), ABS(RESMDIIS)))
      ENDIF
C
C EXTRAPOLATE, SKIP ONLY IN THE FIRST ITERATION
C
      IF(NDIMR.GE.2) THEN
COLD
C
C PUT DIIS RESIDUAL ON DISK - THIS CORRESPONDS TO DIIS EXTRAPOLANT
C
c       CALL PUTLST(SCR(I000),1,1,1,2,LSTERR)
COLD
C
C USE (DIIS EXTRAPOLANT + JACOBI UPDATE) AS NEXT VECTOR
C
C GENERATE UPDATED GUESS, THE JACOBI UPDATE IS ALREADY ON SCR(I000)
C ADD THE SOLUTION VECTORS WITH APPROPRIATE COEFFICIENTS
C
c       CALL ZERO(SCR(I000),LENVEC)
       IOFF=I020
       DO 501 I=1,NDIMR
        FACT=SCR(IOFF)
        CALL GETLST(SCR(I010),I,1,1,1,LSTERR)
        CALL SAXPY (LENVEC,FACT,SCR(I010),1,SCR(I000),1)
        IOFF=IOFF+1
501    CONTINUE
COLD
c       CALL PUTLST(SCR(I000),1,1,1,1,LSTERR)
C
C GENERATE UPDATED JACOBI GUESS
C
c       CALL GETLST(SCR(I010),1,1,1,2,LSTERR)
c       CALL SAXPY (LENVEC,ONE,SCR(I010),1,SCR(I000),1)
COLD
C
C WRITE IT OVER CURRENT INCREMENT LISTS
C
       IOFFT=I000
       IF(SING)THEN
        DO 110 ISPIN=1,1+IUHF
         CALL PUTLST(SCR(IOFFT),1,1,1,IOFF1+ISPIN,LISTT11)
#ifdef _DEBUG_LVLM
         call checksum("T1-DISS2",scr(iofft), IRPDPD(IRREPX,8+ISPIN))
#endif 
         IOFFT=IOFFT+IRPDPD(IRREPX,8+ISPIN)
110     CONTINUE
       ENDIF
       DO 111 ISPIN=3,3-2*IUHF,-1
        LISTT=LISTT20-1+ISPIN
        NSIZE=IDSYMSZ(IRREPX,ISYTYP(1,LISTT),ISYTYP(2,LISTT))
        CALL PUTALL(SCR(IOFFT),NSIZE,IRREPX,LISTT21-1+ISPIN)
#ifdef _DEBUG_LVLM
        call checksum("T2-DISS2",scr(iofft), nsize)
#endif 
        IOFFT=IOFFT+NSIZE
111    CONTINUE
C
C Create the AA list from AB list for RHF calcualtions. For DRCCD
C Calculations, we do not need to do this for every iteration. After
C conergence, we can do this once to get the correct energies. 
C
CSSS       IF (IUHF.EQ.0) CALL RCC_ABTOAA(SCR, MAXCOR*IINTFP, 
CSSS     &                                IUHF, LISTT21)
C
      ENDIF
C
      RETURN
1000  FORMAT(T3,'Current DIIS expansion coefficients : ')
1001  FORMAT(T3,'Convergence information after ',I5,' iterations: ')
1002  FORMAT(T3,'Largest element of residual vector : ',E15.8,'.')
1003  FORMAT(T3,'Largest element of DIIS residual   : ',E15.8,'.')
1004  FORMAT(T3,'Amplitude equations converged in ',I5, 
     &          ' iterations.')
      END

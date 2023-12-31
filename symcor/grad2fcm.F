
C THIS ROUTINE CONSTRUCTS THE FORCE CONSTANT MATRIX FOR SYMMETRY
C BLOCK NDIM FROM NUMERICAL DIFFERENTIATION OF INTERNAL COORDINATE
C (SYMMETRY ADAPTED AND MASS WEIGHTED) GRADIENTS

C 07/98, Extensions for Raman Intensities, Ajith and John

c INPUT
c integer NATOM     :
c integer NDIM      :
c double  GRDPNT(*) :
c double  DIPPNT(*) :
c double  POLPNT(*) :
c integer INVOP(*)  :
c double  SYOP(*)   :
c double  STPSIZ    :
c integer NDSCR     :

c OUTPUT
c double FCM(*)      :
c double DIPDER(*)   :
c double POLDER(*)   :
c double DSCR(NDSCR) :

      SUBROUTINE GRAD2FCM(NATOM,NDIM,
     &                    GRDPNT,FCM,
     &                    DIPPNT,DIPDER,
     &                    POLPNT,POLDER,
     &                    INVOP,SYOP,STPSIZ,DSCR,NDSCR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION GRDPNT(*),FCM(NDIM,NDIM)
      DIMENSION DIPPNT(*),DIPDER(3,NDIM)
      DIMENSION POLPNT(*),POLDER(9,NDIM)
      DIMENSION INVOP(NDIM),SYOP(*)
      double precision dscr(ndscr)

      DIMENSION DA3(3),DA3x3(3,3),WORK3x3(3,3)
      LOGICAL PRINTQ

      COMMON /FLAGS/ IFLAGS(100)
#include "control.com" /* for RAMAN */

      DATA ZILCH /0.0/
      DATA FACT /5.14048D03/

      INDXF(I,J,N)=I+(J-1)*N

      IF (NDSCR.LT.NDIM*(1+2*NDIM)) THEN
         print *, '@GRAD2FCM: Insufficient memory.'
         print *, '           have ',NDSCR,' doubles'
         print *, '           need ',NDIM*(1+2*NDIM),' doubles'
         call aces_exit(1)
      END IF

      PRINTQ=(IFLAGS(1).GE.1)

      NSIZE=3*NATOM

C IF THE REPRESENTATION IS TOTALLY SYMMETRIC, THEN WE HAVE DONE ONLY
C ONE DISPLACEMENT FOR EACH DEGREE OF FREEDOM.
C NOW BUILD THE FCM ROW BY ROW.  EASY STUFF.

      IF (INVOP(1).LE.0) THEN

c      o cases with positive and negative displacements
         IOFF1=0
         IOFF2=NSIZE
         IOFF1D=0
         IOFF2D=3
         IOFF1P=0
         IOFF2P=9
         DO IDIM=1,NDIM
            DO I = 1, NDIM
               FCM(I,IDIM) = GRDPNT(IOFF1+I)-GRDPNT(IOFF2+I)
            END DO
            IOFF1=IOFF1+2*NSIZE
            IOFF2=IOFF2+2*NSIZE
            DIPDER(1,IDIM) = DIPPNT(IOFF1D+1)-DIPPNT(IOFF2D+1)
            DIPDER(2,IDIM) = DIPPNT(IOFF1D+2)-DIPPNT(IOFF2D+2)
            DIPDER(3,IDIM) = DIPPNT(IOFF1D+3)-DIPPNT(IOFF2D+3)
            IOFF1D=IOFF1D+6
            IOFF2D=IOFF2D+6
            DO I = 1, 9
               POLDER(I,IDIM) = POLPNT(IOFF1P+I)-POLPNT(IOFF2P+I)
            END DO
            IOFF1P=IOFF1P+18
            IOFF2P=IOFF2P+18
         END DO
         FFACT=0.5d0/STPSIZ
         CALL DSCAL(NDIM*NDIM,FFACT,FCM,   1)
         CALL DSCAL(   3*NDIM,FFACT,DIPDER,1)
         CALL DSCAL(   9*NDIM,FFACT,POLDER,1)

      ELSE

c      o cases with only one displacement
         IOFF1=1
         IOFF1D=1
         IOFF1P=1
         DO IDIM=1,NDIM
            IOP=INVOP(IDIM)
            IPOS=1+(IOP-1)*9
            CALL DCOPY(NDIM,GRDPNT(IOFF1),1,FCM(1,IDIM),1)
            IOFF1=IOFF1+NSIZE
c         o generate image dipole
            CALL XGEMM('N','N',3,1,3,
     &                 1.d0,SYOP(IPOS),    3,
     &                      DIPPNT(IOFF1D),3,
     &                 0.d0,DIPDER(1,IDIM),3)
            CALL DAXPY(3,-1.d0,DIPPNT(IOFF1D),1,DIPDER(1,IDIM),1)
            IOFF1D=IOFF1D+3
C 07/98 Ajith and John, generate polarizabilty image.
            CALL XGEMM('N','N',3,3,3,
     &                 1.d0,SYOP(IPOS),    3,
     &                      POLPNT(IOFF1P),3,
     &                 0.d0,DA3x3,         3)
            CALL XGEMM('N','T',3,3,3,
     &                 1.d0,DA3x3,         3,
     &                      SYOP(IPOS),    3,
     &                 0.d0,POLDER(1,IDIM),3)
            CALL DAXPY(9,-1.d0,POLPNT(IOFF1P),1,POLDER(1,IDIM),1)
            IOFF1P=IOFF1P+9
         END DO
         FFACT=0.5d0/STPSIZ
         CALL DSCAL(NDIM*NDIM,2.d0*FFACT,FCM,   1)
         CALL DSCAL(   3*NDIM,     FFACT,DIPDER,1)
         CALL DSCAL(   9*NDIM,     FFACT,POLDER,1)

      END IF

CSSS      IF (PRINTQ) THEN
         WRITE(6,1010)
1010     FORMAT(T3,'Force constant matrix : ')
         WRITE(6,1011)((I,J,FCM(I,J),J=1,NDIM),I=1,NDIM)
1011     FORMAT(3('[',I3,',',I3,']',1X,F9.6,1X),'[',I3,',',I3,']',1X,
     &          F9.6)
CSSS      END IF

c   o initialize DSCR pointers for FCM e.vectors and dipder in norm coords
      NFCEV =     1 + NDIM*NDIM
      NDDNC = NFCEV + NDIM*NDIM
C
C Anthony has changed the logic here to call DSYEV instead of EIG
C (EIG itself call DSYEV). I do not find anything really wrong in
C what Anthony did, but some compilers for some frequency calculations
C DSYEV fails. So, I have changed this to the orignal version with
C calls to EIG. This required changes to offsets. The commented line
C are Anthony's offsets. 

c   o diagonalize force constant matrix
      CALL SYMMET2(FCM,NDIM)

CSSS      CALL DCOPY(NDIM*NDIM,FCM,1,DSCR(NFCEV),1)

      CALL DCOPY(NDIM*NDIM,FCM,1,DSCR(1),1)

CSSS      L = NDSCR - NDIM*(1+NDIM)
CSSS      ITMP = NFCEV + NDIM*NDIM
C
      CALL EIG(DSCR(1), DSCR(NFCEV), NDIM, NDIM, 1)

CSSS      IF (NDIM .LE. 4) THEN
CSSS      CALL DSYEV('V','U',NDIM,DSCR(NFCEV),NDIM,DSCR(1),DSCR(NDDNC),9,I)
CSSS      ELSE
CSSS      L = NDIM*NDIM -NDIM
CSSS      CALL DSYEV('V','U',NDIM,DSCR(NFCEV),NDIM,DSCR(1),DSCR(NDDNC),L,I)
CSSS      ENDIF
      
CSSS      if (i.ne.0) then
CSSS         print *, '@GRAD2FCM: Eigensolver failed.'
CSSS         print *, '           dsyev returned error ',i
CSSS         call aces_exit(i)
cSSS      end if

c   o transform dipole derivatives to normal coordinate representation
      CALL XGEMM('N','N',3,NDIM,NDIM,
     &           1.d0,DIPDER,     3,
     &                DSCR(NFCEV),NDIM,
     &           0.d0,DSCR(NDDNC),3)

      IF (RAMAN) THEN
         Write(6, *) "Raman Intensities and depolarization ratios"
         NDX = NFCEV 
         DO IMODE=1,NDIM
c         o transform polarizability derivatives to normal coordinates
            CALL XGEMM('N','N',9,1,NDIM,
     &                 1.d0,POLDER,   9,
     &                      DSCR(NDX),NDIM,
     &                 0.d0,DA3x3,    9)
            NDX = NDX+NDIM
c         o convert the tensor into its principal axis orientation
            CALL DSYEV('N','U',3,DA3x3,3,DA3,WORK3x3,9,I)
            ALPHA = ( DA3(1)
     &               +DA3(2)
     &               +DA3(3))/3.0D0
            BETAS = ( (DA3(1)-DA3(2))**2
     &               +(DA3(2)-DA3(3))**2
     &               +(DA3(3)-DA3(1))**2)*0.5D0
            IF (ALPHA .NE. 0.0D0) THEN
               DTMP = 45.0D0*ALPHA*ALPHA + 4.0D0*BETAS
               RMNIN = 0.0784159D0*(DTMP + 3.0D0*BETAS)
               DEPOL = 3.0D0*BETAS/DTMP
            ELSE
               RMNIN = 0.0784159D0*7.0D0*BETAS
               DEPOL = 0.75D0
            END IF
            Write(6,*) IMODE,RMNIN,DEPOL
         END DO
      END IF

c   o write out frequencies and intensities
      WRITE(6,1000)
1000  FORMAT(T3,'Vibrational frequencies  (cm-1) :')
      NDX = NDDNC
      DO IMODE=1,NDIM
         W=DNRM2(3,DSCR(NDX),1)
         ZINT=W*W*974.868
         NDX=NDX+3
CSSS         X=DSCR(INDXF(IMODE,IMODE,NDIM)+NDIM)
         X=DSCR(INDXF(IMODE,IMODE,NDIM))
         IF (X.LT.ZILCH) THEN
            Z=SQRT(-X)*FACT
            WRITE(6,1001)IMODE,Z,ZINT
1001        FORMAT(T3,I5,1X,F10.5,'i ',F10.5)
         ELSE
            Z=SQRT(X)*FACT
            WRITE(6,1002)IMODE,Z,ZINT
1002        FORMAT(T3,I5,1X,F10.5,2X,F10.5)
         END IF
      END DO

      RETURN
      END


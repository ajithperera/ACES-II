      SUBROUTINE PROCESTEP_XYZ(SCRATCH, AMATRX, TOTREDNCO)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL CONVERGED, CHANGED, NO_ITER
      INTEGER TOTREDNCO, TOTNFBND
      PARAMETER (MAXITER = 50)
C
#include "mxatms.par"
#include "machsp.com"
C
      PARAMETER (EPSILON = 1.0D-5)
      COMMON /USINT/ NX, NXM6, IARCH, NCYCLE, NUNIQUE, NOPT
      COMMON /COORD/ Q(3*MXATMS), R(3*MXATMS), NCON(3*MXATMS),
     &               NR(MXATMS), ISQUASH(3*MXATMS), IATNUM(MXATMS),
     &               ATMASS(MXATMS), IUNIQUE(3*MXATMS), NEQ(3*MXATMS),
     &               IEQUIV(3*MXATMS,3*MXATMS),NOPTI(3*MXATMS),
     &               NATOMS

      DIMENSION SCRATCH(NX*NX), IQFIX(3*MXATMS, 3*MXATMS),
     &          AMATRX(3*NATOMS*TOTREDNCO), STATS(6), 
     &          REDUNCO(MAXREDUNCO), SCRATCH_LCL(3*MAXREDUNCO+6*
     &          MXATMS+15*MXATMS) 
C
      DO I = 1, NUNIQUE
         DO J = 1, NEQ(IUNIQUE(I))
            IQFIX(IUNIQUE(I), J) = IEQUIV(I,J)
         ENDDO
      ENDDO
C
      DO J = 1, NOPT
C
         SCRATCH_LCL(NOPTI(J)) = SCRATCH(J + NOPT)
         SCRATCH_LCL(NOPT + NOPTI(J)) = SCRATCH(J + NOPT)
C
         DO K = 1, NEQ(NOPTI(J))
            SCRATCH_LCL(IQFIX(NOPTI(J), K)) = SCRATCH(J + NOPT)
            SCRATCH_LCL(NOPT + IQFIX(NOPTI(J), K)) = SCRATCH(J + NOPT)  
         ENDDO
      ENDDO
C
#ifdef _DEBUG_LVL0
      WRITE(6,*) "The optimization parameters:nx, nopt, totrednco",
     &            nx, totrednco, nopt
      WRITE(6,*) "The delta q0"
      WRITE(6,*) (SCRATCH_LCL(I), I = 1, TOTREDNCO)
#endif
C
C#ifdef _DNO_SKIP
      IOFF4q0 = 3*NATOMS + 4*NOPT
      IOFF4X0 = IOFF4q0  + TOTREDNCO
C
      ISCRATCH = 3*MAXREDUNCO + 6*MXATMS           
      ITMPQ    = ISCRATCH + 9*NATOMS
      INEWQ    = ITMPQ    + 3*NATOMS
      ILAST    = INEWQ    + 3*NATOMS
C
      CALL DAXPY(TOTREDNCO, -1.0D0, SCRATCH_LCL, 1, R, 1)
      CALL DCOPY(TOTREDNCO, R, 1, SCRATCH_LCL(2*NOPT + 1), 1)
      CALL DCOPY(TOTREDNCO, SCRATCH_LCL, 1, SCRATCH_LCL(IOFF4q0 + 1),
     &           1)
      CALL DCOPY(3*NATOMS, Q, 1, SCRATCH_LCL(IOFF4X0 +1), 1)
C
#ifdef _DEBUG_LVL0
      WRITE(6,*) "The starting Cartesian coords. x0"
      WRITE(6,*) (Q(I), I = 1, 3*NATOMS)
      Write(6,*) "The starting RICs  (q0)"
      WRITE(6,*) (SCRATCH_LCL(2*NOPT + I), I=1, TOTREDNCO)
#endif
C     
      CONVERGED = .FALSE. 
        CHANGED = .FALSE.
        NO_ITER = .FALSE.
          NITER = 0 
C
      DO WHILE (.NOT. CONVERGED .AND. NITER .LE. MAXITER .AND.
     &          .NOT. CHANGED .AND. .NOT. NO_ITER) 
C
         CALL DCOPY(3*NATOMS, Q, 1, SCRATCH_LCL(3*NOPT + 1), 1)
#ifdef _DEBUG_LVL0
         Write(6,*)
         Write(6,*) "The A matrix"
         CALL OUTPUT(AMATRX, 1, 3*NATOMS, 1, TOTREDNCO, 3*NATOMS,
     &               TOTREDNCO, 1) 
         WRITE(6,*) "The starting Cartesian during iteration; x_K"
         WRITE(6,*) (SCRATCH_LCL(3*NOPT+I), I = 1, 3*NATOMS)

#endif
         CALL XGEMM('N', 'N', 3*NATOMS, 1, TOTREDNCO, 1.0D0,
     &               AMATRX, 3*NATOMS, SCRATCH_LCL, TOTREDNCO,
     &               1.0D0, Q, 3*NATOMS)
C
#ifdef _DEBUG_LVL0
      WRITE(6,*) "The new Cartesian Coord., x+delta x0"
      WRITE(6,*) (Q(I), I = 1, 3*NATOMS)
C            ISCRATCH = 3*MAXREDUNCO + 6*MXATMS 
C            ITMPQ    = ISCRATCH + 9*NATOMS
C            INEWQ    = ITMPQ    + 3*NATOMS
C            ILAST    = INEWQ    + 3*NATOMS
C            CALL SYMMETRY(Z(ISCRATCH), Z(ITMPQ), Z(INEWQ), .TRUE.)

#endif
         CALL GEN_NEW_RIC(Q, REDUNCO, IATNUM, NATOMS, TOTNFBND,
     &                    NW_TOTREDNCO)
C
C In rare cases, it is possible that the # of RIC coordinates
C can change. If that happens leave the iterative procedure and
C and use the linear update.
C
         IF (NW_TOREDNCO .GT. TOTREDNCO) CHANGED = .TRUE.
C
#ifdef _DEBUG_LVL0  
      WRITE(6,*) "The new q: q_k"
      WRITE(6,*) (REDUNCO(I), I = 1, TOTREDNCO)
      WRITE(6,*) "The strating q: q0"
      WRITE(6,*) (SCRATCH_LCL(2*NOPT+I), I = 1, TOTREDNCO)
#endif
C
         CALL DCOPY(TOTREDNCO, REDUNCO, 1, R, 1)
         CALL REMOVE_360(R, TOTNFBND, TOTREDNCO)
         CALL DAXPY(TOTREDNCO, -1.0D0, SCRATCH_LCL(2*NOPT + 1), 1, 
     &              R, 1)
#ifdef _DEBUG_LVL0 
      WRITE(6,*) "The (q_k - q0)"
      WRITE(6,*) (R(I), I = 1, TOTREDNCO)
      Write(6,*) "The delta q0" 
      write(6,*) (SCRATCH_LCL(IOFF4q0 + I), I=1, TOTREDNCO)
#endif
C
C Remove +- 360 dgree changes from angel coordinates.
C
CSSS         CALL REMOVE_360(R, TOTNFBND, TOTREDNCO)
C
#ifdef _DEBUG_LVL0
      WRITE(6,*) "after removing -+360 angles"
      WRITE(6,*) (R(I), I = 1, TOTREDNCO)
#endif

         CALL DCOPY(TOTREDNCO, SCRATCH_LCL(IOFF4q0 + 1), 1, 
     &              SCRATCH_LCL, 1)
         CALL DAXPY(TOTREDNCO, -1.0D0, R, 1, SCRATCH_LCL, 1)
#ifdef _DEBUG_LVL0
      WRITE(6,*) "The delta delta q_k" 
      WRITE(6,*) (SCRATCH_LCL(I), I=1, TOTREDNCO)
#endif
C
         IF (NITER .GE. 0) CALL VSTAT(SCRATCH_LCL, STATS, TOTREDNCO)
         IF (NITER .EQ. 0) THEN
             RMSOF_DDQ1 = STATS(5)
         ELSE
             IF (STATS(5) .GE. RMSOF_DDQ1) NO_ITER = .TRUE.
         ENDIF
#ifdef _DEBUG_LVL0
         Print*, "RMS RIC change:", STATS(5)
#endif
C
         CALL DAXPY(3*NATOMS, -1.0D0, Q, 1, SCRATCH_LCL(3*NOPT+1), 1)
         CALL VSTAT(SCRATCH_LCL(3*NOPT+1), STATS, 3*NATOMS)
C
#ifdef _DEBUG_LVL0
         Write(6,*) "Change in Cartesians x_k+1 - x_k"
         WRITE(6,*) (SCRATCH_LCL(3*NOPT+I), I=1, 3*NATOMS)
         Print*, "The root mean sq change of x:", STATS(5)
#endif
         IF (STATS(5)  .LE. EPSILON) CONVERGED = .TRUE.
         NITER = NITER + 1 

      ENDDO 
      IF (CONVERGED)  Print*, "Iterative update conveged in", NITER
C
C If the iterative procedure fail to converge in maximum iterations
C alowed, then do a simple linear update.
C
C#endif
C
      IF (.NOT. CONVERGED) THEN
         CALL DCOPY(3*NATOMS, SCRATCH_LCL(IOFF4X0 +1), 1, Q, 1)
         CALL XGEMM('N', 'N', 3*NATOMS, 1, TOTREDNCO, 1.0D0,
     &              AMATRX, 3*NATOMS, SCRATCH_LCL(IOFF4q0 + 1), 
     &              TOTREDNCO, 1.0D0, Q, 3*NATOMS)
      ENDIF
C
C Store the Cartesian coordinates with a different record label
C from COORD since finite difference calculations use COORD
C for the grid points, not the current point.
C
      CALL PUTREC(1,'JOBARC','COORD_OP',IINTFP*3*NATOMS,Q)
C
#ifdef _DEBUG_LVL0
      WRITE(6,*) "The new Cartesian Coord., x+delta x0"
      WRITE(6,*) (Q(I), I = 1, 3*NATOMS)
#endif
      RETURN
      END

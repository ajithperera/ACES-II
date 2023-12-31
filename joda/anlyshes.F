#include "flags.h"
      SUBROUTINE ANLYSHES(HESMOD, DIAGHES, QST_TANGENT, LST_TANGENT,
     &                    HES, GRDMOD, SCRATCH, PMAT, VEC, NOPT, NX,
     &                    NCYCLE, INR, IVEC, IMODE, IDIE, TS, 
     &                    NRORMANR, RFA, EVFTS, IGTS, QSD, LUOUT,
     &                    QSTLST_CLIMB, NATOMS,
     &                    IPRNT)
C
C Do some printing, and based on the character of the Hessian 
C plus user input decide the what we are going to do.
C
C HES    : The symmetrized Hessian
C HESMOD : The eigenvalues of the symmetrized Hessian
C DIAGHES: The eigenvectors of the symmetrized Hessian
C GRDMOD : The symmetrized gradients
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      LOGICAL TS, NRORMANR, RFA, EVFTS, IGTS, QSD
      LOGICAL LST, QST, QSTORLST, QSTLST_CLIMB, FAILED, RX_PR_EXIST
C
#include "jodaflags.com"
C
      DIMENSION HESMOD(NOPT, NOPT), DIAGHES(NOPT, NOPT), 
     &          SCRATCH(NX*NX), VEC(NOPT), PMAT(NX*NX)
      DIMENSION QST_TANGENT(NOPT), LST_TANGENT(NOPT), GRDMOD(NOPT),
     &          HES(NOPT, NOPT)
C
      INEG = 0
      ZOVLP = 0.0D00
C
      TS       = (INR .EQ. 4 .OR. INR .EQ. 5)
      NRORMANR = (INR .EQ. 1 .OR. INR .EQ. 3)
      RFA      = (INR .EQ. 2)
      EVFTS    = (INR .EQ. 4 .OR. INR .EQ. 5)
      IGTS     = (INR .EQ. 6)
      QSD      = (INR .EQ. 7)
      QSTLST_CLIMB = .FALSE. 
      RX_PR_EXIST  = .FALSE. 
      QST          = .FALSE.
      LST          = .FALSE.
C     
C Let's Print out the eigenvalues and eigenvectors of the Hessian
C
#ifdef _DEBUG_LVL0
          Write(6,*)
          WRITE(LUOUT,2110)
 2110     FORMAT (T3,' The eigenvectors of the Hessian matrix: ')
          CALL HESSOUT(DIAGHES, NOPT, NOPT, 1)
#endif
#ifdef _DEBUG_LVL0
          Write(6,*)
          WRITE(LUOUT,2103)
2103      FORMAT(T3,' The eigenvalues of the Hessian matrix: ')
          WRITE(LUOUT,2101)((HESMOD(I,J),J=I,I),I=1,NOPT)
2101      FORMAT((T3, 6(F10.5,1X)))
#endif
C
      Do Imode = 1, NOPT
         If (Dabs(HesMoD(Imode, Imode)) .LT.  1.0D-09)
     &           HesMOD(Imode, Imode) = 0.0D0
         If (Dabs(GrdMod(Imode)) .LT. 1.0D-09) GrdMod(Imode) = 0.0D0
      Enddo 
       
      DO 1102 J = 1, NOPT
         SCRATCH(J) = HESMOD(J,J)
         IF (HESMOD(J,J).LT.0.0D0) INEG = INEG + 1
 1102 CONTINUE
C
      Write(6,*)
      WRITE(LUOUT, 2102) INEG
 2102 FORMAT(T3,' There are ',i2,' Negative Eigenvalues.')
C
      CALL GETREC(0,'JOBARC','RXSTRUCT',ILENGTH1,TMP)
      CALL GETREC(0,'JOBARC','PRSTRUCT',ILENGTH2,TMP)
      IF (ILENGTH1 .GT. 0 .AND. ILENGTH2 .GT. 0) RX_PR_EXIST = .TRUE. 

      IF (IFLAGS2(h_IFLAGS2_opt_control) .EQ. 1 
     &                                   .AND. RX_PR_EXIST) THEN
         LST = .TRUE. 
      ELSE IF (IFLAGS2(h_IFLAGS2_opt_control) .EQ. 2 
     &                                   .AND. RX_PR_EXIST) THEN
         QST = .TRUE.
      ENDIF 
 
      QSTORLST = (NCYCLE.LE.4.AND.(QST.OR.LST).AND.TS)

#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a,3l)") "Logicals;LST,QST,QSTORLST :", LST, QST, 
     &                                                 QSTORLST
      Write(6,*)
#endif

C I doubt that if we are following an eigenvector that does not
C correspond to the lowest eiegenvalue, the QST or LST is
C useful. That is the reason for test IVEC > 1.
C
      IF (TS.AND.IVEC.GE.1.AND.QSTORLST) THEN

C Get the LST and QST direction vectors. The formulas directly taken from the
C Schlegel, Isreal Journal of Chemistry 33, 449, 1993.

         If (iFlags2(h_IFLAGS2_geom_opt) .EQ. 1) Then
            CALL EVAL_QSTLST_PATH_UDI(QST_TANGENT, LST_TANGENT, 
     &                                SCRATCH, LST, QST)
         Else if (iFlags2(h_IFLAGS2_geom_opt) .EQ. 2) then
            CALL EVAL_QSTLST_PATH_XYZ(QST_TANGENT, LST_TANGENT, 
     &                                SCRATCH, LST, QST)
         Else if (iFlags2(h_IFLAGS2_geom_opt) .Gt. 2) Then
            CALL EVAL_QSTLST_PATH_RIC(QST_TANGENT, LST_TANGENT, 
     &                                SCRATCH, PMAT, LST, QST)   
         Endif
         
C We need to decide whether we choose QST or LST climbing or eigenvector
C following.
C The recommendation in Schlegel et al. is "estimated displacement
C along the tangent vector is greater than 0.05 au. to use climbing.
C Otherwise it is a standard eigenvector following. With regards to
C eigenvector following Schlegel recommend that when the overlap
C between LST or QST tangent with the Hessian eigenvectors is greater
C than 0.08 au, then switch the eigenvectors to follow the one with
C the maximum overlap with the tangent; otherwise, the eigenvector
C with the smallest eigenvalue is followed.
C
C Also, Schlegel recommend to cap the QST or LST climbing steps to 4.
C Since I do not know an empirical formula to estimate the displacement,
C I am going to use simple Newton-Raphson to estimate the stepsize to
C decide whether we are going to climb. That is what's decided in
C QSTLST_OR_EVEC. There is one aspect that is not clear from the papers:
C the direction of the climb. Since the tangent is not an eigenvector
C of the Hessian, I presume when he says that "the tangent to the
C path is used to choose the best eigenvectors for ascent..." what he
C meant is that the eigenvector that has maximum overlap with the tangent.

         IF (LST) THEN
            CALL QSTLST_OR_EVEC(LST_TANGENT, GRDMOD, HESMOD, DIAGHES,
     &                          HES, SCRATCH, IMODE, QSTLST_CLIMB)
            IF (QSTLST_CLIMB) THEN
               CALL FOLOWTANGENT(HESMOD, DIAGHES, LST_TANGENT,
     &                           SCRATCH, NX, NOPT, IMODE)
            ELSE
               CALL WHAT2FOLLOW(HESMOD, DIAGHES, LST_TANGENT,
     &                          SCRATCH, NX, NOPT, IMODE)
   
               CALL FOLOWDFLT(HESMOD, DIAGHES, SCRATCH, VEC, TS,
     &                        NRORMANR, RFA, IVEC, IMODE, NCYCLE,
     &                        NX, NOPT)
            END IF
         ELSE IF (QST) THEN
            CALL QSTLST_OR_EVEC(QST_TANGENT, GRDMOD, HESMOD, DIAGHES,
     &                          HES, SCRATCH, IMODE, QSTLST_CLIMB)
            IF (QSTLST_CLIMB) THEN
               CALL FOLOWTANGENT(HESMOD, DIAGHES, QST_TANGENT,
     &                           SCRATCH, NX, NOPT, IMODE)
            ELSE
               CALL WHAT2FOLLOW(HESMOD, DIAGHES, QST_TANGENT,
     &                          SCRATCH, NX, NOPT, IMODE)

               CALL FOLOWDFLT(HESMOD, DIAGHES, SCRATCH, VEC, TS,
     &                        NRORMANR, RFA, IVEC, IMODE, NCYCLE,
     &                        NX, NOPT)
            END IF
         END IF
         RETURN

C     ELSE IF (.NOT.TS.OR.IVEC.LT.1.OR..NOT.QSTORLST) THEN
      ELSE

C The following logic (bit cleaned up) is what we have been using. It
C is retained just in case the user did not want to do any QST or LST
C eigenvector following. The logic pertains to identifying the eigenvector
C that corresponds to the lowest (IVEC=1, the default), next to the lowest
C (IVEC=2), and so on (user's choice) during the first cycle and then for
C the subsequent steps use the eigenvector that has the largest overlap with
C the previous vector. I think we should add the capability to follow the
C lowest, next to lowest, and so on regardless of what happens in the previous
C cycle. It might just even work better!!!
C Ajith Perera 07/04.
C
  10       CONTINUE
C
           QSTLST_CLIMB = .FALSE. 
C
           CALL FOLOWDFLT(HESMOD, DIAGHES, SCRATCH, VEC, TS, NRORMANR,
     &                    RFA, IVEC, IMODE, NCYCLE, NX, NOPT)

      END IF
C
      RETURN
      END
      

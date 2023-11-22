
#include "flags.h"

      SUBROUTINE DO_RHF_RPA(WORK,MAXCOR,EIGVAL,NSIZ,IUHF)
C
      IMPLICIT INTEGER (A-Z)

      CHARACTER*34 TYPE
      CHARACTER*9 DIFSAME
      DOUBLE PRECISION ONE,ONEM,TWO,FACT,SNRM2,X,ZILCH
      DOUBLE PRECISION EIGVAL,EMIN,TOLER2,TOLER,WORK,EVAL, 
     &                 DET, ECOR, ZSCALE, YSCALE, DSCALE, 
     &                 DDOT
      INTEGER SCR_SIZE, INFO
      logical bUHF, bStop, bWarn
C
      PARAMETER (TOLER=1.D-5, TOLER2=1.D-8)
      INTEGER IONE
      PARAMETER (IONE=1)
      LOGICAL FULLDIAG
      PARAMETER (FULLDIAG=.TRUE.)
C
      DIMENSION WORK(MAXCOR),EIGVAL(2*NSIZ)
#include "machsp.com"
#include "jodaflags.com"

      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NF(2,2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPY(255,2),DIRPRD(8,8)
C
      DATA TWO /2.0/
      DATA ONE /1.0/
      DATA ZILCH/0.0/
      DATA ONEM/-1.0/
      DATA FACT/1.0/

c ----------------------------------------------------------------------

C GET SINGLE EXCITATION DENOMINATORS FOR DIAGONAL ELEMENTS

      ECOR = 0.0D0
      CALL UPDMOI(1,NT(1),3,90,0,0)
C
      CALL GETREC(20, 'JOBARC', 'NBASTOT ', 1, NBAS)
C
      DO IRREP = 1, NIRREP
C
         NUMDIS=IRPDPD(IRREP,ISYTYP(1,25))
         DISSIZ=IRPDPD(IRREP,ISYTYP(2,25))
         SIZE=NUMDIS*DISSIZ

         I000M = 1
         I010M = I000M + 4*SIZE
         I020M = I010M + MAX(DISSIZ,NUMDIS)*2
         I030M = I020M + MAX(DISSIZ,NUMDIS)
         I040M = I030M + NBAS

         CALL MKAMBRHF(WORK(I000M),WORK(I010M),IRREP)
         CALL SINDEN(WORK(I010M),WORK(I030M),IRREP,1)
C
         I000P = I040M 
         I010P = I000P + 4*SIZE
         I020P = I010P + MAX(DISSIZ,NUMDIS)*2
C
         CALL MKAPBRHF(WORK(I000P),WORK(I010P),IRREP)

         DAISIZ = DISSIZ
         ROWDIM = DISSIZ

         IOFFAAM = I000M
         IOFFBBM = (SIZE*2+DISSIZ)+IOFFAAM
#ifdef _DEBUG_LVLM1
        Write(6,*) 
        Write(6,*) "The A-B matrix"
        CALL OUTPUT(WORK(I000M), 1,ROWDIM,1,ROWDIM,ROWDIM,
     &                  ROWDIM, 1)
        Write(6,*) "The A+B matrix"
        CALL OUTPUT(WORK(I000P), 1,ROWDIM,1,ROWDIM,ROWDIM,
     &                  ROWDIM, 1)
#endif
C
         CALL DAXPY(DAISIZ,ONE,WORK(I010M),1,WORK(IOFFAAM),
     &              ROWDIM+1)
         CALL DAXPY(DAISIZ,ONE,WORK(I010M),1,WORK(IOFFBBM),
     &              ROWDIM+1)
C
         IOFFAAP = I000P
         IOFFBBP = (SIZE*2+DISSIZ) + IOFFAAP

         CALL DAXPY(DAISIZ,ONE,WORK(I010M),1,WORK(IOFFAAP),
     &              ROWDIM+1)
         CALL DAXPY(DAISIZ,ONE,WORK(I010M),1,WORK(IOFFBBP),
     &              ROWDIM+1)
c
#ifdef _DEBUG_LVLM1
        Write(6,*)
        Write(6,*) "The A-B matrix"
        CALL OUTPUT(WORK(I000M), 1,ROWDIM,1,ROWDIM,ROWDIM,
     &                  ROWDIM, 1)
        Write(6,*) "The A+B matrix"
        CALL OUTPUT(WORK(I000P), 1,ROWDIM,1,ROWDIM,ROWDIM,
     &                  ROWDIM, 1)
#endif
C
        I030P = I020P + ROWDIM*ROWDIM
        I040P = I030P + ROWDIM*ROWDIM
        IEND  = I040P + ROWDIM*ROWDIM
        CALL EIG(WORK(I000P), WORK(I020P), 1, ROWDIM, 1)
C
#ifdef _DEBUG_LVLM1
        Write(6,*)
        Write(6,*) "The eigenvalues of A+B"
        CALL OUTPUT(WORK(I000P), 1,ROWDIM,1,ROWDIM,ROWDIM,
     &                  ROWDIM, 1)
        Write(6,*) "The eigenvectors of A+B"
        CALL OUTPUT(WORK(I020P), 1,ROWDIM,1,ROWDIM,ROWDIM,
     &                  ROWDIM, 1)
#endif 
C
        DO ILOOP = 1, ROWDIM
           IOFF =  (ILOOP - 1)*ROWDIM + ILOOP - 1
           IF (WORK(I000P + IOFF) .LT. 0.0D0) Then
              WRITE(6,"(a)")"The A-B matrix is not positive definite"
              CALL ERREX
           ENDIF
           WORK(I000P + IOFF) = DSQRT(WORK(I000P + IOFF))
        ENDDO
#ifdef _DEBUG_LVLM1
        Write(6,*) 
        Write(6,*) "The sqrt-eigenvalues of A+B"
        CALL OUTPUT(WORK(I000P), 1,ROWDIM,1,ROWDIM,ROWDIM,
     &                  ROWDIM, 1)
#endif
        CALL DZERO(WORK(I030P), ROWDIM*ROWDIM)
        CALL DZERO(WORK(I040P), ROWDIM*ROWDIM)
        DO ILOOP = 1, ROWDIM
           IDIAG = (ILOOP - 1)*ROWDIM + ILOOP - 1 + I000P
           IELEM = (ILOOP - 1)*ROWDIM + I020P
           IDQRT = (ILOOP - 1)*ROWDIM + I030P
           ISQRT = (ILOOP - 1)*ROWDIM + I040P
           CALL DAXPY(ROWDIM, 1.0D0/WORK(IDIAG), WORK(IELEM), 1,
     &                WORK(IDQRT), 1)
           CALL DAXPY(ROWDIM, WORK(IDIAG), WORK(IELEM), 1, 
     &                WORK(ISQRT), 1)
        ENDDO
#ifdef _DEBUG_LVLM1
        Write(6,*) 
        Write(6,*) "(A+B)*X1/2"
        CALL OUTPUT(WORK(I040P), 1,ROWDIM,1,ROWDIM,ROWDIM,
     &                  ROWDIM, 1)
#endif
        CALL XGEMM("T", "N", ROWDIM, ROWDIM, ROWDIM, 1.0D0, 
     &              WORK(I040P), ROWDIM, WORK(I000M), ROWDIM,
     &              0.0D0, WORK(I020P), ROWDIM)

        CALL XGEMM("N", "N", ROWDIM, ROWDIM, ROWDIM, 1.0D0, 
     &              WORK(I020P), ROWDIM, WORK(I040P), ROWDIM,
     &              0.0D0, WORK(I000M), ROWDIM)
#ifdef _DEBUG_LVLM1
        Write(6,*) 
        Write(6,*) "The scaled of A-B"
        CALL OUTPUT(WORK(I000M), 1,ROWDIM,1,ROWDIM,ROWDIM,
     &                  ROWDIM, 1)
#endif
         IF (FULLDIAG) THEN
CSSS            CALL DSYEV('V','L',ROWDIM,WORK(I000M),ROWDIM,
CSSS     &                 EIGVAL,WORK(I030P),SCR_SIZE,INFO)
            IEVEC = IEND
            IZVEC = IEVEC + ROWDIM*ROWDIM
            IYVEC = IZVEC + ROWDIM*ROWDIM
            ITVEC = IYVEC + ROWDIM*ROWDIM
            ISVEC = ITVEC + ROWDIM*ROWDIM
            IT2AM = ISVEC + ROWDIM*ROWDIM
            IDUMY = IT2AM + ROWDIM*ROWDIM
            IENDD = IDUMY + ROWDIM

            CALL EIG(WORK(I000M), WORK(IEVEC), 1, ROWDIM, 1)
            DO ILOOP = 1, ROWDIM
               IDIAG = (ILOOP - 1)*ROWDIM + ILOOP - 1 + I000M
               EIGVAL(ILOOP) = WORK(IDIAG)
            ENDDO

C#ifdef _DEBUG_LVLM1
            Write(6,*) "The eigenvalues of RPA matrix"
            Write(6, "(4F10.5)") (Dsqrt(EIGVAL(I))*27.2113961, 
     &               I=1, ROWDIM)
CSSS            Write(6, "(4F10.5)") (DSqrt(EIGVAL(I)), I=1, ROWDIM)
            Write(6,*) "The eigen vectors"
        CALL OUTPUT(WORK(IEVEC), 1,ROWDIM,1,ROWDIM,ROWDIM,
     &                  ROWDIM, 1)
C#endif
            CALL DZERO(WORK(IZVEC), ROWDIM*ROWDIM)
            CALL DZERO(WORK(IYVEC), ROWDIM*ROWDIM)
C
            DO ILOOP = 1, ROWDIM
               IELEM = (ILOOP - 1)*ROWDIM + IEVEC
               IOFFZ = (ILOOP - 1)*ROWDIM + IZVEC
               IOFFY = (ILOOP - 1)*ROWDIM + IYVEC
               EVAL = DSQRT(DSQRT(EIGVAL(ILOOP)))
               CALL DAXPY(ROWDIM, EVAL, WORK(IELEM), 1, 
     &                    WORK(IOFFZ), 1)
               CALL DAXPY(ROWDIM, 1.0D0/EVAL, WORK(IELEM), 1,
     &                    WORK(IOFFY), 1)
            ENDDO
#ifdef _DEBUG_LVLM1
        Write(6,*)
        Write(6,*) "The ZVEC"
        CALL OUTPUT(WORK(IZVEC), 1,ROWDIM,1,ROWDIM,ROWDIM,
     &                  ROWDIM, 1)
        Write(6,*) "The YVEC"
        CALL OUTPUT(WORK(IYVEC), 1,ROWDIM,1,ROWDIM,ROWDIM,
     &                  ROWDIM, 1)
       
#endif 
               CALL XGEMM("N", "N", ROWDIM, ROWDIM, ROWDIM, 1.0D0, 
     &                     WORK(I040P), ROWDIM, WORK(IYVEC), ROWDIM,
     &                     0.0D0, WORK(ITVEC), ROWDIM)
               CALL XGEMM("N", "N", ROWDIM, ROWDIM, ROWDIM, 1.0D0,
     &                     WORK(I030P), ROWDIM, WORK(IZVEC), ROWDIM,
     &                     0.0D0, WORK(ISVEC), ROWDIM)
#ifdef _DEBUG_LVLM1
        Write(6,*)
        Write(6,*) "The TVEC"
        CALL OUTPUT(WORK(ITVEC), 1,ROWDIM,1,ROWDIM,ROWDIM,
     &                  ROWDIM, 1)
        Write(6,*) "The SVEC"
        CALL OUTPUT(WORK(ISVEC), 1,ROWDIM,1,ROWDIM,ROWDIM,
     &                  ROWDIM, 1)

#endif
               DO J = 1, ROWDIM
                  DO I = 1, ROWDIM
                     IYOFF = I + (J-1)*ROWDIM + IYVEC - 1
                     IzOFF = I + (J-1)*ROWDIM + IZVEC - 1
                     ITOFF = I + (J-1)*ROWDIM + ITVEC - 1
                     ISOFF = I + (J-1)*ROWDIM + ISVEC - 1

                     WORK(IZOFF) = (WORK(ISOFF) + WORK(ITOFF))/2.0D0
                     WORK(IYOFF) = (WORK(ISOFF) - WORK(ITOFF))/2.0D0
                  ENDDO
               ENDDO
#ifdef _DEBUG_LVLM1
        Write(6,*)
        Write(6,*) "The ZVEC"
        CALL OUTPUT(WORK(IZVEC), 1,ROWDIM,1,ROWDIM,ROWDIM,
     &                  ROWDIM, 1)
        Write(6,*) "The YVEC"
        CALL OUTPUT(WORK(IYVEC), 1,ROWDIM,1,ROWDIM,ROWDIM,
     &                  ROWDIM, 1)
#endif
C
              DO IVEC = 1, ROWIDIM
                 IZOFF = (IVEC - 1)*ROWDIM + 1 + IZVEC
                 IYOFF = (IVEC - 1)*ROWDIM + 1 + IYVEC
                 ZSCALE= DDOT(ROWDIM, WORK(IZOFF), 1, WORK(IZOFF),1)
                 YSCALE= DDOT(ROWDIM, WORK(IYOFF), 1, WORK(IYOFF),1)
                 DSCALE= 1.0D0/DSQRT(ZSCALE+YSCALE)
                 CALL DSCAL(ROWDIM, DSCALE, WORK(IZOFF), 1)
                 CALL DSCAL(ROWDIM, DSCALE, WORK(IYOFF), 1)
              ENDDO
C
              CALL MINV(WORK(IZVEC), ROWDIM, ROWDIM, WORK(IDUMY), 
     &                  DET, 1.0D-8, 0, 1)
              CALL XGEMM("N", "N", ROWDIM, ROWDIM, ROWDIM, 1.0D0,
     &                    WORK(IYVEC), ROWDIM, WORK(IZVEC), ROWDIM,
     &                    0.0D0, WORK(IT2AM), ROWDIM)
#ifdef _DEBUG_LVLM1
        Write(6,*)
        Write(6,*) "The TVEC"
        CALL OUTPUT(WORK(IT2AM), 1,ROWDIM,1,ROWDIM,ROWDIM,
     &                  ROWDIM, 1)
#endif
               CALL MKB_RHF(WORK(I000M), WORK(I010M), IRREP)
#ifdef _DEBUG_LVL0
        Write(6,*)
        Write(6,*) "The B matrix"
        CALL OUTPUT(WORK(I000M), 1,ROWDIM,1,ROWDIM,ROWDIM,
     &                  ROWDIM, 1)
#endif

               ECOR = EcOR + DDOT(ROWDIM*ROWDIM,WORK(I000M), 1,
     &                            WORK(IT2AM), 1)

               Write(6,*)
               Write(6, "(a,F12.10)")"The r-CCD Cooreation energy = ",
     &                                  1.0D0/4.0D0*ECOR
            INFO = 0
            IF (INFO.NE.0) THEN
               WRITE(*,*)
     & '@RHFSTAB: There was a problem diagonalizing the matrix.'
               CALL ERREX
            END IF
            ILOCEVC=I000P
            ILOCSCR=I010P
         ELSE
            NROOT=MAX(48/NIRREP,IFLAGS(75)+5)
            NROOT=MIN(ROWDIM,NROOT,15)
            MAXIT=IFLAGS2(105)
            I015=I010+MAXIT*MAXIT
            I020=I015+ROWDIM*NROOT
            I030=I020+ROWDIM*(MAXIT+NROOT)
            I040=I030+MAX(ROWDIM,NROOT)
            I050=I040+ROWDIM
            I060=I050+MAXIT*MAX(ROWDIM,MAXIT)
            I070=I060+ROWDIM
            CALL DAVID1(WORK(I000P), WORK(I010P),
     &                  EIGVAL,      WORK(I015),
     &                  WORK(I020), WORK(I030),
     &                  WORK(I040), WORK(I050),
     &                  WORK(I060), MAXIT,
     &                  NROOT,       IROOT,
     &                  ROWDIM,      DISSIZ,
     &                  TOLER,       IERR)
            ILOCEVC=I015
            ILOCSCR=I030
            WRITE(*,*) 'Davidson procedure found ',IROOT,' roots.'
            CALL SRTEVAL(IROOT,         ROWDIM,     EIGVAL,
     &                   WORK(ILOCEVC),WORK(I020),WORK(I030))
            NROOT=IROOT
         END IF
      ENDDO

      RETURN
      END





      SUBROUTINE DO_UHF_RPA(WORK,MAXCOR,EIGVAL,NSIZ,IUHF)
      IMPLICIT INTEGER (A-Z)

      CHARACTER*34 TYPE
      CHARACTER*9 DIFSAME
      DOUBLE PRECISION ONE,ONEM,TWO,FACT,SNRM2,X,ZILCH
      DOUBLE PRECISION EIGVAL,EMIN,TOLER2,TOLER,WORK,EVAL,
     &       DET, ECOR, ZSCALE, YSCALE, DSCALE,DDOT
      INTEGER SCR_SIZE, INFO
      logical bWarn

      PARAMETER (TOLER=1.D-5, TOLER2=1.0D-8 )
      INTEGER IONE
      PARAMETER (IONE=1)
      LOGICAL FULLDIAG
      PARAMETER (FULLDIAG=.TRUE.)

      DIMENSION WORK(MAXCOR),EIGVAL(2*NSIZ)
C
#include "machsp.com"
#include "jodaflags.com"

      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NF(2,2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPY(255,2),DIRPRD(8,8)

      DATA TWO /2.0/
      DATA ONE /1.0/
      DATA ZILCH/0.0/
      DATA ONEM/-1.0/
      DATA FACT/1.0/
C
C LOOP OVER IRREPS OF SPIN ORBITAL A MATRIX
C
C GET SINGLE EXCITATION DENOMINATORS FOR DIAGONAL ELEMENTS
C
      ISIZEA=NT(1)
      ISIZEB=NT(2)
      CALL UPDMOI(1,NT(1),3,90,0,0)
      CALL UPDMOI(1,NT(2),4,90,0,0)
C
      CALL INIPCK(1, 9,  9 , 200, 0, 0, 1)
      CALL INIPCK(1, 10, 10, 201, 0, 0, 1)
      CALL INIPCK(1, 10, 9,  202, 0, 0, 1)
CSSS      CALL INIPCK(1, 9, 10,  202, 0, 0, 1)
C
      CALL GETREC(20, 'JOBARC', 'NBASTOT ', 1, NBAS)

      DO IRREP=1,NIRREP
         SIZEAA=IRPDPD(IRREP,ISYTYP(1,19))
         SIZEBB=IRPDPD(IRREP,ISYTYP(1,20))
         ROWDIM=SIZEAA+SIZEBB
       IF (ROWDIM .NE. ZILCH) THEN
C
         IOFFAAAAM = 1
         IOFFBBBBM = IOFFAAAAM + (SIZEAA*ROWDIM+SIZEAA)
         IOFFBBAAM = IOFFAAAAM + ISIZEA
         IOFFAABBM = IOFFAAAAM + (SIZEAA*ROWDIM)
         SIZTOT    = ROWDIM*ROWDIM

         I000M = 1
         I010M = I000M + SIZTOT
         I020M = I010M + 2*ROWDIM
 
         CALL MKAMBUHF(WORK(I000M),WORK(I010M),IRREP)

         I000P = I020M
         I010P = I000P + SIZTOT
         I020P = I010P + 2*ROWDIM

         CALL MKAPBUHF(WORK(I000P),WORK(I010P),IRREP)

#ifdef _DEBUG_LVLM1
        Write(6,*) 
        Write(6,*) "The A-B matrix"
        CALL OUTPUT(WORK(I000M), 1,ROWDIM,1,ROWDIM,ROWDIM,
     &                  ROWDIM, 1)
        Write(6,*) "The A+B matrix"
        CALL OUTPUT(WORK(I000P), 1,ROWDIM,1,ROWDIM,ROWDIM,
     &                  ROWDIM, 1)
#endif
         IOFFA = I020P
         IOFFB = IOFFA + SIZEAA
         IDENA = IOFFB + SIZEBB
         IDENB = IDENA + NBAS
         IEND  = IDENB + NBAS
C
         CALL SINDEN(WORK(IOFFA),WORK(IDENA),IRREP,1)
         CALL SINDEN(WORK(IOFFB),WORK(IDENB),IRREP,2)
C
         CALL DAXPY(SIZEAA,ONE,WORK(IOFFA),1,WORK(IOFFAAAAM),ROWDIM+1)
         CALL DAXPY(SIZEBB,ONE,WORK(IOFFB),1,WORK(IOFFBBBBM),ROWDIM+1)
C
         IOFFAAAAP = I000P
         IOFFBBBBP = IOFFAAAAP + (SIZEAA*ROWDIM+SIZEAA)
         IOFFBBAAP = IOFFBBBBP + ISIZEA

         CALL DAXPY(SIZEAA,ONE,WORK(IOFFA),1,WORK(IOFFAAAAP),ROWDIM+1)
         CALL DAXPY(SIZEAA,ONE,WORK(IOFFB),1,WORK(IOFFBBBBP),ROWDIM+1)
C
#ifdef _DEBUG_LVLM1
        Write(6,*) "The A-B matrix"
        CALL OUTPUT(WORK(I000M), 1,ROWDIM,1,ROWDIM,ROWDIM,
     &                  ROWDIM, 1)
        Write(6,*) "The A+B matrix"
        CALL OUTPUT(WORK(I000P), 1,ROWDIM,1,ROWDIM,ROWDIM,
     &                  ROWDIM, 1)
#endif
         I020P  = IEND
         I030P  = I020P + ROWDIM*ROWDIM
         I040P  = I030P + ROWDIM*ROWDIM
         IEND   = I040P + ROWDIM*ROWDIM
C
CSS         SCR_SIZE = MAXCOR - IEND
CSS         IF(IEND.GT.MAXCOR)CALL INSMEM('DO_UHF_RPA', IEND, MAXCOR)

         CALL EIG(WORK(I000P), WORK(I020P), 1, ROWDIM, 1)

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
               WRITE(6,"(a)")"The A+B matrix is not positive definite"
CSSS               CALL ERREX
               GO TO 100
            ENDIF
            WORK(I000P + IOFF) = DSQRT(WORK(I000P + IOFF))
         ENDDO
C
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
C
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

        IF (FULLDIAG)THEn
C
            IEVEC = IEND
            IZVEC = IEVEC + ROWDIM*ROWDIM
            IYVEC = IZVEC + ROWDIM*ROWDIM
            ITVEC = IYVEC + ROWDIM*ROWDIM
            ISVEC = ITVEC + ROWDIM*ROWDIM
            IT2AM = ISVEC + ROWDIM*ROWDIM
            IDUMY = IT2AM + ROWDIM*ROWDIM
            IENDD = IDUMY + ROWDIM
            IF(IENDD.GT.MAXCOR)CALL INSMEM('DO_UHF_RPA', IEND, MAXCOR)

            CALL EIG(WORK(I000M), WORK(IEVEC), 1, ROWDIM, 1)
            DO ILOOP = 1, ROWDIM
               IDIAG = (ILOOP - 1)*ROWDIM + ILOOP - 1 + I000M
               EIGVAL(ILOOP) = WORK(IDIAG)
               IF (EIGVAL(ILOOP) .LT. 0.0D0) THEN
                  WRITE(6,"(a)")"There are negative excitations"
CSSS               CALL ERREX
                   GO TO 100
               ENDIF
            ENDDO

#ifdef _DEBUG_LVL0
            Write(6,*) "The eigenvalues of RPA matrix"
            Write(6, "(4F10.5)") (DSQRT(EIGVAL(I))*27.2113961, 
     &               I=1, ROWDIM)
CSSS        CALL OUTPUT(WORK(IEVEC), 1,ROWDIM,1,ROWDIM,ROWDIM,
CSSS     &                  ROWDIM, 1)
#endif
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
     &                   WORK(I040P), ROWDIM, WORK(IYVEC), ROWDIM,
     &                   0.0D0, WORK(ITVEC), ROWDIM)
             CALL XGEMM("N", "N", ROWDIM, ROWDIM, ROWDIM, 1.0D0,
     &                   WORK(I030P), ROWDIM, WORK(IZVEC), ROWDIM,
     &                   0.0D0, WORK(ISVEC), ROWDIM)
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
                   IZOFF = I + (J-1)*ROWDIM + IZVEC - 1
                   ITOFF = I + (J-1)*ROWDIM + ITVEC - 1
                   ISOFF = I + (J-1)*ROWDIM + ISVEC - 1

                   WORK(IZOFF) = (WORK(ISOFF) + WORK(ITOFF))/2.0D0
                   WORK(IYOFF) = (WORK(ISOFF) - WORK(ITOFF))/2.0D0
                ENDDO
             ENDDO
#ifdef _DEBUG_LVM1
        Write(6,*)
        Write(6,*) "The ZVEC"
        CALL OUTPUT(WORK(IZVEC), 1,ROWDIM,1,ROWDIM,ROWDIM,
     &                  ROWDIM, 1)
        Write(6,*) "The YVEC"
        CALL OUTPUT(WORK(IYVEC), 1,ROWDIM,1,ROWDIM,ROWDIM,
     &                  ROWDIM, 1)

#endif
             ZSCALE=0.0D0
             DO IVEC = 1, ROWDIM
                IZOFF  = (IVEC - 1)*ROWDIM  + IZVEC
                IYOFF  = (IVEC - 1)*ROWDIM  + IYVEC
                ZSCALE = DDOT(ROWDIM, WORK(IZOFF), 1, WORK(IZOFF),1)
                YSCALE = DDOT(ROWDIM, WORK(IYOFF), 1, WORK(IYOFF),1)
                DSCALE = 1.0D0/DSQRT(ZSCALE+YSCALE)
                CALL DSCAL(ROWDIM, DSCALE, WORK(IZOFF), 1)
               CALL DSCAL(ROWDIM, DSCALE, WORK(IYOFF), 1)
             ENDDO
#ifdef _DEBUG_LVLM1
        Write(6,*)
        Write(6,*) "Normalized ZVEC"
        CALL OUTPUT(WORK(IZVEC), 1,ROWDIM,1,ROWDIM,ROWDIM,
     &                  ROWDIM, 1)
        Write(6,*) "Normalized YVEC"
        CALL OUTPUT(WORK(IYVEC), 1,ROWDIM,1,ROWDIM,ROWDIM,
     &                  ROWDIM, 1)
      
#endif
             CALL MINV(WORK(IZVEC), ROWDIM, ROWDIM, WORK(IDUMY),
     &                 DET, 1.0D-8, 0, 1)
#ifdef _DEBUG_LVLM1
        Write(6,*)
        Write(6,*) "The T1 VEC"
        CALL OUTPUT(WORK(IZVEC), 1,ROWDIM,1,ROWDIM,ROWDIM,
     &                  ROWDIM, 1)
#endif
             CALL XGEMM("N", "N", ROWDIM, ROWDIM, ROWDIM, 1.0D0,
     &                   WORK(IYVEC), ROWDIM, WORK(IZVEC), ROWDIM,
     &                   0.0D0, WORK(IT2AM), ROWDIM)
C
#ifdef _DEBUG_LVLM1
        Write(6,*)
        Write(6,*) "The T2 VEC"
        CALL OUTPUT(WORK(IT2AM), 1,ROWDIM,1,ROWDIM,ROWDIM,
     &                  ROWDIM, 1)
#endif
             CALL DZERO(WORK(I000M), ROWDIM*ROWDIM)
             CALL MKB_UHF(WORK(I000M), WORK(I010M), IRREP)
#ifdef _DEBUG_LVLM1
        Write(6,*)
        Write(6,*) "The B Matrix"
        CALL OUTPUT(WORK(I000M), 1,ROWDIM,1,ROWDIM,ROWDIM,
     &                  ROWDIM, 1)
#endif

             ECOR = ECOR + DDOT(ROWDIM*ROWDIM,WORK(I000M), 1,
     &                            WORK(IT2AM), 1)
  
             Write(6,*)
#ifdef _RCC
             Write(6, "(a,F12.10)")"The r-CCD Correlation energy = ",
     &                                  1.0D0/4.0D0*ECOR
#else
             Write(6, "(a,F12.10)")"The dr-CCD Corrleation energy = ",
     &                                  1.0D0/2.0D0*ECOR
#endif
             LOAVA = IRPDPD(IRREP,  9)
             LOBVB = IRPDPD(IRREP, 10)
CSSS             LOAVB = IRPDPD(IRREP, 18)
             LOAVB = IRPDPD(IRREP, 9)
             LOBVA = IRPDPD(IRREP, 12)
               
             LT2AAAA = LOAVA*LOAVA
             LT2BBBB = LOBVB*LOBVB
             LT2ABAB = LOAVB*LOAVB
             LT2BABA = LOBVA*LOBVA
C
             IT2AA = IENDD 
             IT2BB = IT2AA + LT2AAAA
             IT2AB = IT2BB + LT2BBBB
             IT2BA = IT2AB + LT2ABAB
             IENDD = IT2BA + LT2BABA
             IF(IENDD.GT.MAXCOR)CALL INSMEM('DO_UHF_RPA', IEND, 
     &                                       MAXCOR)

             DO ILOOP = 1, LOAVA
                IORGNAA = (ILOOP - 1)*ROWDIM + IT2AM
                ITRGTAA = (ILOOP - 1)*LOAVA  + IT2AA
C
                IORGNBA =  LOAVA + (ILOOP - 1)*ROWDIM + IT2AM
                ITRGTBA = (ILOOP - 1)*LOBVA  + IT2BA
C
                CALL DCOPY(LOAVA, WORK(IORGNAA), 1, WORK(ITRGTAA),
     &                     1)             
                CALL DCOPY(LOBVA, WORK(IORGNBA), 1, WORK(ITRGTBA),
     &                     1)             
             ENDDO
             DO ILOOP = 1, LOBVB
                IORGNBB = ROWDIM*LOAVA + LOBVB + (ILOOP - 1)*ROWDIM
     &                      + IT2AM
                ITRGTBB = (ILOOP - 1)*LOBVB + IT2BB
C
                IORGNAB = ROWDIM*LOAVA +  (ILOOP - 1)*ROWDIM 
     &                    + IT2AM 
                ITRGTAB = (ILOOP - 1)*LOAVB  + IT2AB
                CALL DCOPY(LOBVB, WORK(IORGNBB), 1, WORK(ITRGTBB),
     &                     1)             
                CALL DCOPY(LOAVB, WORK(IORGNAB), 1, WORK(ITRGTAB),
     &                     1)             
             ENDDO

#ifdef _DEBUG_LVLM1
        Write(6,*)
        Write(6,*) "The T2AA"
        CALL OUTPUT(WORK(IT2AA), 1, LOAVA, 1, LOAVA, LOAVA, LOAVA, 1)
        Write(6,*)
        Write(6,*) "The T2BB"
        CALL OUTPUT(WORK(IT2BB), 1, LOBVB, 1, LOBVB, LOBVB, LOBVB, 1)
        Write(6,*)
        Write(6,*) "The T2AB"
        CALL OUTPUT(WORK(IT2AB), 1, LOAVB, 1, LOAVB, LOAVB, LOAVB, 1)
        Write(6,*)
        Write(6,*) "The T2BA"
        CALL OUTPUT(WORK(IT2BA), 1, LOBVA, 1, LOBVA, LOBVA, LOBVA, 1)
#endif
            CALL PUTLST(WORK(IT2AA), 1, LOAVA, 1, IRREP, 200)
            CALL PUTLST(WORK(IT2BB), 1, LOBVB, 1, IRREP, 201)
            CALL PUTLST(WORK(IT2AB), 1, LOAVB, 1, IRREP, 202)
            
CSSS            CALL MKT1AT2(WORK(IT2AA), WORK(IT2BB), WORK(IT2AB),
CSSS     &                   LT2AAAA, LT2BBBB, LT2ABAB, WORK(IENDD),
CSSS     &                   (MAXCOR - IENDD), IRREP)
        ELSE
C
            NROOT = MAX(48/NIRREP,IFLAGS(75)+5)
            NROOT = MIN(ROWDIM,NROOT,15)
            MAXIT = IFLAGS2(105)
C
            I015 = I010+MAXIT*MAXIT
            I020 = I015+ROWDIM*NROOT
            I030 = I020+ROWDIM*(MAXIT+NROOT)
            I040 = I030+MAX(NROOT,ROWDIM)
            I050 = I040+ROWDIM
            I060 = I050+MAXIT*MAX(ROWDIM,MAXIT)
            I070 = I060+ROWDIM
            CALL DAVID1(WORK(I000),WORK(I010),EIGVAL,WORK(I015),
     &                  WORK(I020),WORK(I030),WORK(I040),WORK(I050),
     &                  WORK(I060),MAXIT,NROOT,IROOT,ROWDIM,ROWDIM,
     &                  TOLER,IERR)
            ILOCEVC=I015
            ILOCSCR=I030
            WRITE(6,1010)IROOT
1010        FORMAT(T3,'Davidson procedure found ',I5,' roots.')
            CALL SRTEVAL(IROOT,ROWDIM,EIGVAL,WORK(ILOCEVC),
     &                   WORK(I020),WORK(I030))
            NROOT=IROOT
        ENDIF
       ENDIF
  100  CONTINUE
C
      ENDDO
      CALL MKT1AT2(WORK, MAXCOR)
C
      RETURN
      END

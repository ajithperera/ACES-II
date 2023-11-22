      SUBROUTINE BUILDX(X, FAB, FIJ, IRREP0, NOCC, NVRT, ISPIN)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER POP,VRT,DIRPRD
      DIMENSION X(1),IOFFO(8),IOFFV(8),SCR(1), FAB(NVRT,NVRT), 
     &          FIJ(NOCC,NOCC)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYM/ POP(8,2),VRT(8,2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPY(255,2),DIRPRD(8,8)
C
      IOFFO(1) = 1
      IOFFV(1) = 1
      IOFF     = 1
C
      DO IRREP=2,NIRREP
         IOFFO(IRREP)=IOFFO(IRREP-1)+POP(IRREP-1,ISPIN)
         IOFFV(IRREP)=IOFFV(IRREP-1)+VRT(IRREP-1,ISPIN)
      ENDDO
C
      DO IRREPR = 1, NIRREP
         IRREPL = DIRPRD(IRREPR, IRREP0)
         IOFFO0 = IOFFO(IRREPR)
         IOFFV0 = IOFFV(IRREPL)
         JOFFO0 = IOFFO(IRREPR)
         JOFFV0 = IOFFV(IRREPL)
       
         IOFFV00 = IOFFV0 
         JOFFV00 = JOFFV0
         IOFFO00 = IOFFO0
         JOFFO00 = IOFFO0

         Write(6,"(2(I2))") IRREPL,IRREPR
         Write(6,"(2(1x,I2))") POP(IRREPR, ISPIN), VRT(IRREPL, ISPIN)
         Write(6,*)
         DO JOCCR = 1, POP(IRREPR, ISPIN)
            DO JVRTL = 1, VRT(IRREPL, ISPIN)

               DO IOCCR = 1, POP(IRREPR, ISPIN)
                  DO IVRTL = 1, VRT(IRREPL, ISPIN)
                  Write(6,"(4(1x,I2))") IOFFV0, IOFFO0, JOFFV0, JOFFO0
          
                     IF (IOFFV0 .EQ. JOFFV0 .AND. IOFFO0 .EQ. JOFFO0)
     &                THEN
                      X(IOFF) = FAB(IOFFV0, JOFFV0) - 
     &                         FIJ(IOFFO0, JOFFO0)
                     Write(6, "(F10.5)") FAB(IOFFV0, JOFFV0) -
     &                         FIJ(IOFFO0, JOFFO0)

                     ELSE IF (IOFFV0 .EQ. JOFFV0) THEN
                      X(IOFF) = -FIJ(IOFFO0, JOFFO0)
                      Write(6, "(F10.5)") -FIJ(IOFFO0, JOFFO0)
                      ELSE IF (IOFFO0 .EQ. JOFFO0) THEN
                      X(IOFF) = FAB(IOFFV0, JOFFV0)
                      Write(6, "(F10.5)")  FAB(IOFFV0, JOFFV0)
                     ELSE  
                      X(IOFF) = 0.0D0
                     ENDIF 
                     IOFF = IOFF   + 1
                   IOFFV0 = IOFFV0 + 1
                  ENDDO
                  IOFFO0 = IOFFO0 + 1
                  IOFFV0 = IOFFV00
               ENDDO
               JOFFV0 = JOFFV0 + 1 
               IOFFO0 = IOFFO00
            ENDDO
            JOFFO0 = JOFFO0 + 1
            JOFFV0 = JOFFV00
          ENDDO
      ENDDO
C
      Write(6,*) ioff
      Write(6, "(5(1x,F10.5))") (X(I), I=1, IOff)
      RETURN
      END

      SUBROUTINE Y10C(CORE,MAXCOR,IUHF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER A,C,B,D,E,F,BD,EF,DF,BF,EB,ED,ELF,BLD,BDEF,BLDELF,
     1        EBDF,EDBF,
     1        DISSIZ,DISSIZA,DISSIZB
      INTEGER POP,VRT,DIRPRD
      DIMENSION CORE(1),DISSIZ(8),DISSIZA(8),DISSIZB(8),NDIS(8)
      DIMENSION IADY5A(8),LENY5A(8),IOFFW(8)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYM/    POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /INFO/   NOCA,NOCB,NVRTA,NVRTB
      COMMON /T3OFF/  IOFFVV(8,8,10),IOFFOO(8,8,10),IOFFVO(8,8,4)
C
cYAU      INDEX(I) = I * (I-1)/2
      INDEX(I) = ishft((I*(I-1)),-1)

#ifdef _DEBUG_VCC5Q
      WRITE(*,*) '@VCC5Q: Entered ', __FILE__
c1000 FORMAT(' @Y10C-I, Welcome to one of the black holes ! ')
#endif /* _DEBUG */

C     AAAA/BBBB = AAAA * AAAA - AAAA * AAAA / BBBB * BBBB - BBBB * BBBB
C
      DO  200 ISPIN=1,IUHF+1
C
C     Read and sort y3 ring intermediates.
C
C     We read KBEJ as (E,K,B,J). Sort to EBKJ.
C
      LEN = 0
      DO IRREP = 1,NIRREP
         LEN = LEN + IRPDPD(IRREP,8+ISPIN) * IRPDPD(IRREP,8+ISPIN)
      END DO
C
      I000 = 1
      I010 = I000 + LEN
      I020 = I010 + LEN
      I030 = I020 +   2 * (NOCA * NOCB + NVRTA * NVRTB)

cYAU - this format is used many places
 1020 FORMAT('@Y10C: Insufficient memory. Need ',I15,'. Have ',I15,'.')

      NEED = IINTFP * I030
      IF (NEED.GT.MAXCOR) THEN
         WRITE(6,1020) NEED,MAXCOR
         STOP 'Y10C'
      END IF
C
C     Get the whole list at I010.
C
      CALL GETALL(CORE(I010),LEN,1,53+ISPIN)
C
C     Sort to different order : EBKJ is at I000.
C
      CALL SSTGEN(CORE(I010),CORE(I000),LEN,
     1            VRT(1,ISPIN),POP(1,ISPIN),VRT(1,ISPIN),POP(1,ISPIN),
     1            CORE(I020),1,'1324')
C
C     Note that y3 has same symmetry as integrals at present. this is
C     exploited here, but it should be remembered that this may not be
C     in the case of Q (rather than Q*) corrections.
C
C     BDEF = KBEJ * JDFK - KDEJ * JBFK
C                   KFDJ          KFBJ
C
C     KBEJ was             (E,K,B,J)   now is (E,B,K,J)
C     JDFK was (F,J,D,K) = (D,K,F,J)   now is (D,F,K,J)
C     KDEJ was             (E,K,D,J)   now is (E,D,K,J)
C     JBFK was (F,J,B,K) = (B,K,F,J)   now is (B,F,K,J)
C
C     Therefore we are trying to form
C
C     W(B,D,E,F) = y3(E,B,KJ) * y3(D,F,KJ) - y3(E,D,KJ) * y3(B,F,KJ)
C
C     Somewhat crude, in core implementation as a trial (assumes that
C     all symmetry allowed (B,D,E,F) can be held in core).
C
      LNABCD = 0
      DO IRREP = 1,NIRREP
         DISSIZ(IRREP) = 0
         NDIS(IRREP)   = 0
         DO JRREP=1,NIRREP
            KRREP = DIRPRD(IRREP,JRREP)
            DISSIZ(IRREP) = DISSIZ(IRREP)
     1                    + VRT(JRREP,ISPIN) * VRT(KRREP,ISPIN)
            NDIS(IRREP)   = NDIS(IRREP)
     1                    + POP(JRREP,ISPIN) * POP(KRREP,ISPIN)
         END DO
         LNABCD = LNABCD + DISSIZ(IRREP) * DISSIZ(IRREP)
      END DO
C
      I020 = I010 + LNABCD
C
      NEED = IINTFP * I020
      IF (NEED.GT.MAXCOR) THEN
         WRITE(6,1020) NEED,MAXCOR
         STOP 'Y10C'
      END IF
C
      IOFFY3 = I000
      DO IRPJK=1,NIRREP
         IF (IRPJK.EQ.1) THEN
            IOFFW(IRPJK) = I010
         ELSE
            IOFFW(IRPJK) = IOFFW(IRPJK-1)
     &                   + DISSIZ(IRPJK-1)*DISSIZ(IRPJK-1)
         END IF
         CALL XGEMM('N','T',DISSIZ(IRPJK),DISSIZ(IRPJK),NDIS(IRPJK),
     &              -4.0D+00,CORE(IOFFY3),DISSIZ(IRPJK),
     &                       CORE(IOFFY3),DISSIZ(IRPJK),
     &              0.0D+00, CORE(IOFFW(IRPJK)),DISSIZ(IRPJK))
         IOFFY3 = IOFFY3 + DISSIZ(IRPJK) * NDIS(IRPJK)
      END DO
C
      DO  190 IRPEF=1,NIRREP
      IRPBD = IRPEF
C
C     Get a symmetry block of Y10(BD,EF) at I020.
C
      I030 = I020 + IRPDPD(IRPEF,ISPIN) * IRPDPD(IRPEF,ISPIN)
      NEED = IINTFP * I030
      IF (NEED.GT.MAXCOR) THEN
         WRITE(6,1020) NEED,MAXCOR
         STOP 'Y10C'
      END IF

cYAU - This huge loop was rearranged to allow the GNU optimizer unfettered
c      access to the nirvana we call "optimum code". The main change was
c      setting MAX?=VRT(IRP?,ISPIN). After that, I noticed the loops were
c      really over (?-1) instead of (?) so I replaced the index and shifted
c      the loop bounds. Ah, the joys of hand-optimized code. BTW, I only
c      found this when a GNU-compiled version crashed in this routine.
cYAU - anchor - search for 'anchor' to hit all of my changes
      CALL GETLIST(CORE(I020),1,IRPDPD(IRPEF,ISPIN),1,IRPEF,3+ISPIN)
      DO IRPF = 1,NIRREP
         IRPE = DIRPRD(IRPEF,IRPF)
         IF (IRPE.LE.IRPF) THEN
            MAXE = VRT(IRPE,ISPIN)
            MAXF = VRT(IRPF,ISPIN)
            DO IRPD = 1,NIRREP
               IRPB = DIRPRD(IRPBD,IRPD)
               IF (IRPB.LE.IRPD) THEN
                  MAXB = VRT(IRPB,ISPIN)
                  MAXD = VRT(IRPD,ISPIN)
                  IRPDF = DIRPRD(IRPD,IRPF)
                  IRPBF = DIRPRD(IRPB,IRPF)
                  IRPEB = DIRPRD(IRPE,IRPB)
                  IRPED = DIRPRD(IRPE,IRPD)
                  IF (IRPEF.EQ.1) THEN
                     DO F = 1,MAXF-1
                        DO E = 0,F-1
                           DO D = 1,MAXD-1
                              DO B = 0,D-1
      DF = IOFFVV(IRPF,IRPDF,2+ISPIN) + F*MAXD + D
      BF = IOFFVV(IRPF,IRPBF,2+ISPIN) + F*MAXB + B
      ED = IOFFVV(IRPD,IRPED,2+ISPIN) + D*MAXE + E
      EB = IOFFVV(IRPB,IRPEB,2+ISPIN) + B*MAXE + E
      EBDF = IOFFW(IRPDF) + DF*DISSIZ(IRPDF) + EB
      EDBF = IOFFW(IRPBF) + BF*DISSIZ(IRPBF) + ED
      ELF = IOFFVV(IRPF,IRPEF,ISPIN) + ishft((F*(F-1)),-1) + E
      BLD = IOFFVV(IRPD,IRPBD,ISPIN) + ishft((D*(D-1)),-1) + B
      BLDELF = ELF*IRPDPD(IRPBD,ISPIN) + BLD
      CORE(I020 + BLDELF) = CORE(I020 + BLDELF)
     &                    + CORE(EBDF) - CORE(EDBF)
c                             END DO B = 0,D-1
                              END DO
c                          END DO D = 1,MAXD-1
                           END DO
c                       END DO E = 0,F-1
                        END DO
c                    END DO F = 1,MAXF-1
                     END DO
                  ELSE
                     DO F = 0,MAXF-1
                        DO E = 0,MAXE-1
                           DO D = 0,MAXD-1
                              DO B = 0,MAXB-1
      DF = IOFFVV(IRPF,IRPDF,2+ISPIN) + F*MAXD + D
      BF = IOFFVV(IRPF,IRPBF,2+ISPIN) + F*MAXB + B
      ED = IOFFVV(IRPD,IRPED,2+ISPIN) + D*MAXE + E
      EB = IOFFVV(IRPB,IRPEB,2+ISPIN) + B*MAXE + E
      EBDF = IOFFW(IRPDF) + DF*DISSIZ(IRPDF) + EB
      EDBF = IOFFW(IRPBF) + BF*DISSIZ(IRPBF) + ED
      ELF = IOFFVV(IRPF,IRPEF,ISPIN) + F*MAXE + E
      BLD = IOFFVV(IRPD,IRPBD,ISPIN) + D*MAXB + B
      BLDELF = ELF*IRPDPD(IRPBD,ISPIN) + BLD
      CORE(I020 + BLDELF) = CORE(I020 + BLDELF)
     &                    + CORE(EBDF) - CORE(EDBF)
c                             END DO B = 0,MAXB-1
                              END DO
c                          END DO D = 0,MAXD-1
                           END DO
c                       END DO E = 0,MAXE-1
                        END DO
c                    END DO F = 0,MAXF-1
                     END DO
c                 END IF (IRPEF.EQ.1)
                  END IF
c              END IF (IRPB.LE.IRPD)
               END IF
c           END DO IRPD = 1,NIRREP
            END DO
c        END IF (IRPE.LE.IRPF)
         END IF
c     END DO IRPF = 1,NIRREP
      END DO
      CALL PUTLIST(CORE(I020),1,IRPDPD(IRPEF,ISPIN),1,IRPEF,3+ISPIN)

c     END DO 190 IRPEF=1,NIRREP
  190 CONTINUE
c     END DO 200 ISPIN=1,IUHF+1
  200 CONTINUE

c ----------------------------------------------------------------------

C                 kBEj   jDFk   kDEj   jBFk   KbeJ   JdfK   KdeJ   JbfK
C
C     AAAA/BBBB = BAAB * BAAB - BAAB * BAAB / ABBA * ABBA - ABBA * ABBA
C
      DO  400 ISPIN=1,IUHF+1
C
      IF (ISPIN.EQ.1) THEN
         ISPIN1 = 1
         ISPIN2 = 2
      ELSE
         ISPIN1 = 2
         ISPIN2 = 1
      END IF
C
C     Read and sort y3 ring intermediates.
C
C     We read KBEJ as (E,K,B,J). Sort to EBKJ.
C
      LEN = 0
      DO IRREP=1,NIRREP
         LEN = LEN + IRPDPD(IRREP,10+ISPIN) * IRPDPD(IRREP,10+ISPIN)
      END DO
C
      I000 = 1
      I010 = I000 + LEN
      I020 = I010 + LEN
      I030 = I020 +   2 * (NOCA * NOCB + NVRTA * NVRTB)
C
      NEED = IINTFP * I030
      IF (NEED.GT.MAXCOR) THEN
         WRITE(6,1020) NEED,MAXCOR
         STOP 'Y10C'
      END IF
C
C     Get the whole list at I010.
C
      CALL GETALL(CORE(I010),LEN,1,57+ISPIN)
C
C     Sort to different order : EBKJ is at I000.
C
      CALL SSTGEN(CORE(I010),CORE(I000),LEN,
     1        VRT(1,ISPIN1),POP(1,ISPIN2),VRT(1,ISPIN1),POP(1,ISPIN2),
     1        CORE(I020),1,'1324')
C
C     Note that y3 has same symmetry as integrals at present. this is
C     exploited here, but it should be remembered that this may not be
C     in the case of Q (rather than Q*) corrections.
C
C     BDEF = KBEJ * JDFK - KDEJ * JBFK
C                   KFDJ          KFBJ
C
C     KBEJ was             (E,K,B,J)   now is (E,B,K,J)
C     JDFK was (F,J,D,K) = (D,K,F,J)   now is (D,F,K,J)
C     KDEJ was             (E,K,D,J)   now is (E,D,K,J)
C     JBFK was (F,J,B,K) = (B,K,F,J)   now is (B,F,K,J)
C
C     Therefore we are trying to form
C
C     W(B,D,E,F) = y3(E,B,KJ) * y3(D,F,KJ) - y3(E,D,KJ) * y3(B,F,KJ)
C
C     Somewhat crude, in core implementation as a trial (assumes that
C     all symmetry allowed (B,D,E,F) can be held in core).
C
      LNABCD = 0
      DO IRREP=1,NIRREP
         DISSIZ(IRREP) = 0
         NDIS(IRREP)   = 0
         DO JRREP=1,NIRREP
            KRREP = DIRPRD(IRREP,JRREP)
            DISSIZ(IRREP) = DISSIZ(IRREP)
     &                    + VRT(JRREP,ISPIN) * VRT(KRREP,ISPIN)
            NDIS(IRREP)   = NDIS(IRREP)
     &                    + POP(JRREP,ISPIN2) * POP(KRREP,ISPIN2)
         END DO
         LNABCD = LNABCD + DISSIZ(IRREP) * DISSIZ(IRREP)
      END DO
C
      I020 = I010 + LNABCD
C
      NEED = IINTFP * I020
      IF (NEED.GT.MAXCOR) THEN
         WRITE(6,1020) NEED,MAXCOR
         STOP 'Y10C'
      END IF
C
      IOFFY3 = I000
      DO IRPJK=1,NIRREP
         IF (IRPJK.EQ.1) THEN
            IOFFW(IRPJK) = I010
         ELSE
            IOFFW(IRPJK) = IOFFW(IRPJK-1)
     &                   + DISSIZ(IRPJK-1)*DISSIZ(IRPJK-1)
         END IF
         CALL XGEMM('N','T',DISSIZ(IRPJK),DISSIZ(IRPJK),NDIS(IRPJK),
     &              -4.0D+00,
     &              CORE(IOFFY3),DISSIZ(IRPJK),
     &              CORE(IOFFY3),DISSIZ(IRPJK),0.0D+00,
     &              CORE(IOFFW(IRPJK)),DISSIZ(IRPJK))
         IOFFY3 = IOFFY3 + DISSIZ(IRPJK) *   NDIS(IRPJK)
      END DO
C
      DO  390 IRPEF=1,NIRREP
      IRPBD = IRPEF
C
C     Get a symmetry block of Y10(BD,EF) at I020.
C
      I030 = I020 + IRPDPD(IRPEF,ISPIN) * IRPDPD(IRPEF,ISPIN)
      NEED = IINTFP * I030
      IF (NEED.GT.MAXCOR) THEN
         WRITE(6,1020) NEED,MAXCOR
         STOP 'Y10C'
      END IF
C
cYAU - anchor
      CALL GETLIST(CORE(I020),1,IRPDPD(IRPEF,ISPIN),1,IRPEF,3+ISPIN)
      DO IRPF = 1,NIRREP
         IRPE = DIRPRD(IRPEF,IRPF)
         IF (IRPE.LE.IRPF) THEN
            MAXE = VRT(IRPE,ISPIN)
            MAXF = VRT(IRPF,ISPIN)
            DO IRPD = 1,NIRREP
               IRPB = DIRPRD(IRPBD,IRPD)
               IF (IRPB.LE.IRPD) THEN
                  MAXB = VRT(IRPB,ISPIN)
                  MAXD = VRT(IRPD,ISPIN)
                  IRPDF = DIRPRD(IRPD,IRPF)
                  IRPBF = DIRPRD(IRPB,IRPF)
                  IRPEB = DIRPRD(IRPE,IRPB)
                  IRPED = DIRPRD(IRPE,IRPD)
                  IF (IRPEF.EQ.1) THEN
                     DO F = 1,MAXF-1
                        DO E = 0,F-1
                           DO D = 1,MAXD-1
                              DO B = 0,D-1
      DF = IOFFVV(IRPF,IRPDF,2+ISPIN) + F*MAXD + D
      BF = IOFFVV(IRPF,IRPBF,2+ISPIN) + F*MAXB + B
      ED = IOFFVV(IRPD,IRPED,2+ISPIN) + D*MAXE + E
      EB = IOFFVV(IRPB,IRPEB,2+ISPIN) + B*MAXE + E
      EBDF = IOFFW(IRPDF) + DF*DISSIZ(IRPDF) + EB
      EDBF = IOFFW(IRPBF) + BF*DISSIZ(IRPBF) + ED
      ELF = IOFFVV(IRPF,IRPEF,ISPIN) + ishft((F*(F-1)),-1) + E
      BLD = IOFFVV(IRPD,IRPBD,ISPIN) + ishft((D*(D-1)),-1) + B
      BLDELF = ELF*IRPDPD(IRPBD,ISPIN) + BLD
      CORE(I020 + BLDELF) = CORE(I020 + BLDELF)
     &                    + CORE(EBDF) - CORE(EDBF)
c                             END DO B = 0,D-1
                              END DO
c                          END DO D = 1,MAXD-1
                           END DO
c                       END DO E = 0,F-1
                        END DO
c                    END DO F = 1,MAXF-1
                     END DO
                  ELSE
                     DO F = 0,MAXF-1
                        DO E = 0,MAXE-1
                           DO D = 0,MAXD-1
                              DO B = 0,MAXB-1
      DF = IOFFVV(IRPF,IRPDF,2+ISPIN) + F*MAXD + D
      BF = IOFFVV(IRPF,IRPBF,2+ISPIN) + F*MAXB + B
      ED = IOFFVV(IRPD,IRPED,2+ISPIN) + D*MAXE + E
      EB = IOFFVV(IRPB,IRPEB,2+ISPIN) + B*MAXE + E
      EBDF = IOFFW(IRPDF) + DF*DISSIZ(IRPDF) + EB
      EDBF = IOFFW(IRPBF) + BF*DISSIZ(IRPBF) + ED
      ELF = IOFFVV(IRPF,IRPEF,ISPIN) + F*MAXE + E
      BLD = IOFFVV(IRPD,IRPBD,ISPIN) + D*MAXB + B
      BLDELF = ELF*IRPDPD(IRPBD,ISPIN) + BLD
      CORE(I020 + BLDELF) = CORE(I020 + BLDELF)
     &                    + CORE(EBDF) - CORE(EDBF)
c                             END DO B = 0,MAXB-1
                              END DO
c                          END DO D = 0,MAXD-1
                           END DO
c                       END DO E = 0,MAXE-1
                        END DO
c                    END DO F = 0,MAXF-1
                     END DO
c                 END IF (IRPEF.EQ.1)
                  END IF
c              END IF (IRPB.LE.IRPD)
               END IF
c           END DO IRPD = 1,NIRREP
            END DO
c        END IF (IRPE.LE.IRPF)
         END IF
c     END DO IRPF = 1,NIRREP
      END DO
      CALL PUTLIST(CORE(I020),1,IRPDPD(IRPEF,ISPIN),1,IRPEF,3+ISPIN)

c     END DO 390 IRPEF=1,NIRREP
  390 CONTINUE
c     END DO 400 ISPIN=1,IUHF+1
  400 CONTINUE

c ----------------------------------------------------------------------

C            kBEj   kfdj
C
C     ABAB = BAAB * BBBB
C
C     Read and sort y3 ring intermediates.
C
C     We read kBEj as (E,k,B,j). Sort to EBkj.
C
      LEN = 0
      DO  410 IRREP=1,NIRREP
      LEN = LEN + IRPDPD(IRREP,11) * IRPDPD(IRREP,11)
  410 CONTINUE
C
      I000 = 1
      I010 = I000 + LEN
      I020 = I010 + LEN
      I030 = I020 +   2 * (NOCA * NOCB + NVRTA * NVRTB)
C
      NEED = IINTFP * I030
      IF(NEED.GT.MAXCOR)THEN
      WRITE(6,1020) NEED,MAXCOR
      STOP 'Y10C'
      ENDIF
C
C     Get the whole list at I010.
C
      CALL GETALL(CORE(I010),LEN,1,58)
C
C     Sort to different order : EBkj is at I000.
C
      CALL SSTGEN(CORE(I010),CORE(I000),LEN,
     1        VRT(1,1),POP(1,2),VRT(1,1),POP(1,2),
     1        CORE(I020),1,'1324')
C
C     We read kfdj as (f,k,d,j). Sort to fdkj.
C
      LEN = 0
      DO  420 IRREP=1,NIRREP
      LEN = LEN + IRPDPD(IRREP,10) * IRPDPD(IRREP,10)
  420 CONTINUE
C
      I020 = I010 + LEN
      I030 = I020 + LEN
      I040 = I030 +   2 * (NOCA * NOCB + NVRTA * NVRTB)
C
      NEED = IINTFP * I040
      IF(NEED.GT.MAXCOR)THEN
      WRITE(6,1020) NEED,MAXCOR
      STOP 'Y10C'
      ENDIF
C
C     Get the whole list at I020.
C
      CALL GETALL(CORE(I020),LEN,1,54+IUHF)
C
C     Sort to different order : fdkj is at I010.
C
      CALL SSTGEN(CORE(I020),CORE(I010),LEN,
     1        VRT(1,2),POP(1,2),VRT(1,2),POP(1,2),
     1        CORE(I030),1,'1324')
C
C     Note that y3 has same symmetry as integrals at present. this is
C     exploited here, but it should be remembered that this may not be
C     in the case of Q (rather than Q*) corrections.
C
C     BDEF = KBEJ * JDFK - KDEJ * JBFK
C                   KFDJ          KFBJ
C
C     KBEJ was             (E,K,B,J)   now is (E,B,K,J)
C     JDFK was (F,J,D,K) = (D,K,F,J)   now is (D,F,K,J)
C     KDEJ was             (E,K,D,J)   now is (E,D,K,J)
C     JBFK was (F,J,B,K) = (B,K,F,J)   now is (B,F,K,J)
C
C     Therefore we are trying to form
C
C     W(B,D,E,F) = y3(E,B,KJ) * y3(D,F,KJ) - y3(E,D,KJ) * y3(B,F,KJ)
C
C     Somewhat crude, in core implementation as a trial (assumes that
C     all symmetry allowed (B,D,E,F) can be held in core).
C
      LNABCD = 0
      DO  440 IRREP=1,NIRREP
      DISSIZA(IRREP) = 0
      DISSIZB(IRREP) = 0
      NDIS(IRREP)   = 0
      DO  430 JRREP=1,NIRREP
      KRREP = DIRPRD(IRREP,JRREP)
      DISSIZA(IRREP) = DISSIZA(IRREP)
     1                 + VRT(JRREP,1) * VRT(KRREP,1)
      DISSIZB(IRREP) = DISSIZB(IRREP)
     1                 + VRT(JRREP,2) * VRT(KRREP,2)
      NDIS(IRREP)   = NDIS(IRREP)
     1                 + POP(JRREP,2) * POP(KRREP,2)
  430 CONTINUE
      LNABCD = LNABCD + DISSIZA(IRREP) * DISSIZB(IRREP)
  440 CONTINUE
C
      I030 = I020 + LNABCD
C
      NEED = IINTFP * I030
      IF(NEED.GT.MAXCOR)THEN
      WRITE(6,1020) NEED,MAXCOR
      STOP 'Y10C'
      ENDIF
C
      IOFFY3A = I000
      IOFFY3B = I010
      DO  450 IRPJK=1,NIRREP
C
      IF(IRPJK.EQ.1)THEN
      IOFFW(IRPJK) = I020
      ELSE
      IOFFW(IRPJK) = IOFFW(IRPJK-1) + DISSIZA(IRPJK-1)*DISSIZB(IRPJK-1)
      ENDIF
C
      CALL XGEMM('N','T',DISSIZA(IRPJK),DISSIZB(IRPJK),NDIS(IRPJK),
     1            -4.0D+00,
     1            CORE(IOFFY3A),DISSIZA(IRPJK),
     1            CORE(IOFFY3B),DISSIZB(IRPJK),0.0D+00,
     1            CORE(IOFFW(IRPJK)),DISSIZA(IRPJK))
C
      IOFFY3A = IOFFY3A + DISSIZA(IRPJK) *   NDIS(IRPJK)
      IOFFY3B = IOFFY3B + DISSIZB(IRPJK) *   NDIS(IRPJK)
C
  450 CONTINUE
C
      DO  590 IRPEF=1,NIRREP
      IRPBD = IRPEF
C
C     Get a symmetry block of Y10(BD,EF) at I030.
C
      I040 = I030 + IRPDPD(IRPEF,13) * IRPDPD(IRPEF,13)
      NEED = IINTFP * I040
      IF (NEED.GT.MAXCOR) THEN
         WRITE(6,1020) NEED,MAXCOR
         STOP 'Y10C'
      END IF
C
cYAU - anchor
      CALL GETLIST(CORE(I030),1,IRPDPD(IRPEF,13),1,IRPEF,6)
      DO IRPF = 1,NIRREP
         IRPE = DIRPRD(IRPEF,IRPF)
         MAXE = VRT(IRPE,1)
         MAXF = VRT(IRPF,2)
         DO IRPD = 1,NIRREP
            IRPB = DIRPRD(IRPBD,IRPD)
            MAXB = VRT(IRPB,1)
            MAXD = VRT(IRPD,2)
            IRPDF = DIRPRD(IRPD,IRPF)
            IRPEB = DIRPRD(IRPE,IRPB)
            DO F = 0,MAXF-1
               DO E = 0,MAXE-1
                  DO D = 0,MAXD-1
                     DO B = 0,MAXB-1
      DF = IOFFVV(IRPF,IRPDF,4) + F*MAXD + D
      EB = IOFFVV(IRPB,IRPEB,3) + B*MAXE + E
      EBDF = IOFFW(IRPDF) + DF*DISSIZA(IRPEB) + EB
      EF = IOFFVV(IRPF,IRPEF,5) + F*MAXE + E
      BD = IOFFVV(IRPD,IRPBD,5) + D*MAXB + B
      BDEF = EF*IRPDPD(IRPBD,13) + BD
      CORE(I030 + BDEF) = CORE(I030 + BDEF) + CORE(EBDF)
c                    END DO B = 0,MAXB-1
                     END DO
c                 END DO D = 0,MAXD-1
                  END DO
c              END DO E = 0,MAXE-1
               END DO
c           END DO F = 0,MAXF-1
            END DO
c        END DO IRPD = 1,NIRREP
         END DO
c     END DO IRPF = 1,NIRREP
      END DO
      CALL PUTLIST(CORE(I030),1,IRPDPD(IRPEF,13),1,IRPEF,6)
C
  590 CONTINUE
C
C
C            KBEJ   KfdJ
C
C     ABAB = AAAA * ABBA
C
C     Read and sort y3 ring intermediates.
C
C     We read KBEJ as (E,K,B,J). Sort to EBKJ.
C
      LEN = 0
      DO  610 IRREP=1,NIRREP
      LEN = LEN + IRPDPD(IRREP,9) * IRPDPD(IRREP,9)
  610 CONTINUE
C
      I000 = 1
      I010 = I000 + LEN
      I020 = I010 + LEN
      I030 = I020 +   2 * (NOCA * NOCB + NVRTA * NVRTB)
C
      NEED = IINTFP * I030
      IF(NEED.GT.MAXCOR)THEN
      WRITE(6,1020) NEED,MAXCOR
      STOP 'Y10C'
      ENDIF
C
C     Get the whole list at I010.
C
      CALL GETALL(CORE(I010),LEN,1,54)
C
C     Sort to different order : EBKJ is at I000.
C
      CALL SSTGEN(CORE(I010),CORE(I000),LEN,
     1        VRT(1,1),POP(1,1),VRT(1,1),POP(1,1),
     1        CORE(I020),1,'1324')
C
C     We read KfdJ as (f,K,d,J). Sort to fdKJ.
C
      LEN = 0
      DO  620 IRREP=1,NIRREP
      LEN = LEN + IRPDPD(IRREP,12) * IRPDPD(IRREP,12)
  620 CONTINUE
C
      I020 = I010 + LEN
      I030 = I020 + LEN
      I040 = I030 +   2 * (NOCA * NOCB + NVRTA * NVRTB)
C
      NEED = IINTFP * I040
      IF(NEED.GT.MAXCOR)THEN
      WRITE(6,1020) NEED,MAXCOR
      STOP 'Y10C'
      ENDIF
C
C     Get the whole list at I020.
C
      CALL GETALL(CORE(I020),LEN,1,58+IUHF)
C
C     Sort to different order : fdKJ is at I010.
C
      CALL SSTGEN(CORE(I020),CORE(I010),LEN,
     1        VRT(1,2),POP(1,1),VRT(1,2),POP(1,1),
     1        CORE(I030),1,'1324')
C
C     Note that y3 has same symmetry as integrals at present. this is
C     exploited here, but it should be remembered that this may not be
C     in the case of Q (rather than Q*) corrections.
C
C     BDEF = KBEJ * JDFK - KDEJ * JBFK
C                   KFDJ          KFBJ
C
C     KBEJ was             (E,K,B,J)   now is (E,B,K,J)
C     JDFK was (F,J,D,K) = (D,K,F,J)   now is (D,F,K,J)
C     KDEJ was             (E,K,D,J)   now is (E,D,K,J)
C     JBFK was (F,J,B,K) = (B,K,F,J)   now is (B,F,K,J)
C
C     Therefore we are trying to form
C
C     W(B,D,E,F) = y3(E,B,KJ) * y3(D,F,KJ) - y3(E,D,KJ) * y3(B,F,KJ)
C
C     Somewhat crude, in core implementation as a trial (assumes that
C     all symmetry allowed (B,D,E,F) can be held in core).
C
      LNABCD = 0
      DO  640 IRREP=1,NIRREP
      DISSIZA(IRREP) = 0
      DISSIZB(IRREP) = 0
      NDIS(IRREP)   = 0
      DO  630 JRREP=1,NIRREP
      KRREP = DIRPRD(IRREP,JRREP)
      DISSIZA(IRREP) = DISSIZA(IRREP)
     1                 + VRT(JRREP,1) * VRT(KRREP,1)
      DISSIZB(IRREP) = DISSIZB(IRREP)
     1                 + VRT(JRREP,2) * VRT(KRREP,2)
      NDIS(IRREP)   = NDIS(IRREP)
     1                 + POP(JRREP,1) * POP(KRREP,1)
  630 CONTINUE
      LNABCD = LNABCD + DISSIZA(IRREP) * DISSIZB(IRREP)
  640 CONTINUE
C
      I030 = I020 + LNABCD
C
      NEED = IINTFP * I030
      IF(NEED.GT.MAXCOR)THEN
      WRITE(6,1020) NEED,MAXCOR
      STOP 'Y10C'
      ENDIF
C
      IOFFY3A = I000
      IOFFY3B = I010
      DO  650 IRPJK=1,NIRREP
C
      IF(IRPJK.EQ.1)THEN
      IOFFW(IRPJK) = I020
      ELSE
      IOFFW(IRPJK) = IOFFW(IRPJK-1) + DISSIZA(IRPJK-1)*DISSIZB(IRPJK-1)
      ENDIF
C
      CALL XGEMM('N','T',DISSIZA(IRPJK),DISSIZB(IRPJK),NDIS(IRPJK),
     1            -4.0D+00,
     1            CORE(IOFFY3A),DISSIZA(IRPJK),
     1            CORE(IOFFY3B),DISSIZB(IRPJK),0.0D+00,
     1            CORE(IOFFW(IRPJK)),DISSIZA(IRPJK))
C
      IOFFY3A = IOFFY3A + DISSIZA(IRPJK) *   NDIS(IRPJK)
      IOFFY3B = IOFFY3B + DISSIZB(IRPJK) *   NDIS(IRPJK)
C
  650 CONTINUE
C
      DO  790 IRPEF=1,NIRREP
      IRPBD = IRPEF
C
C     Get a symmetry block of Y10(BD,EF) at I030.
C
      I040 = I030 + IRPDPD(IRPEF,13) * IRPDPD(IRPEF,13)
      NEED = IINTFP * I040
      IF(NEED.GT.MAXCOR)THEN
      WRITE(6,1020) NEED,MAXCOR
      STOP 'Y10C'
      ENDIF
C
cYAU - anchor
      CALL GETLIST(CORE(I030),1,IRPDPD(IRPEF,13),1,IRPEF,6)
      DO IRPF = 1,NIRREP
         IRPE = DIRPRD(IRPEF,IRPF)
         MAXE = VRT(IRPE,1)
         MAXF = VRT(IRPF,2)
         DO IRPD = 1,NIRREP
            IRPB = DIRPRD(IRPBD,IRPD)
            MAXB = VRT(IRPB,1)
            MAXD = VRT(IRPD,2)
            IRPDF = DIRPRD(IRPD,IRPF)
            IRPEB = DIRPRD(IRPE,IRPB)
            DO F = 0,MAXF-1
               DO E = 0,MAXE-1
                  DO D = 0,MAXD-1
                     DO B = 0,MAXB-1
      DF = IOFFVV(IRPF,IRPDF,4) + F*MAXD + D
      EB = IOFFVV(IRPB,IRPEB,3) + B*MAXE + E
      EBDF = IOFFW(IRPDF) + DF*DISSIZA(IRPEB) + EB
      EF = IOFFVV(IRPF,IRPEF,5) + F*MAXE + E
      BD = IOFFVV(IRPD,IRPBD,5) + D*MAXB + B
      BDEF = (EF)*IRPDPD(IRPBD,13) + BD
      CORE(I030 + BDEF) = CORE(I030 + BDEF) + CORE(EBDF)
c                    END DO B = 0,MAXB-1
                     END DO
c                 END DO D = 0,MAXD-1
                  END DO
c              END DO E = 0,MAXE-1
               END DO
c           END DO F = 0,MAXF-1
            END DO
c        END DO IRPD = 1,NIRREP
         END DO
c     END DO IRPF = 1,NIRREP
      END DO
      CALL PUTLIST(CORE(I030),1,IRPDPD(IRPEF,13),1,IRPEF,6)
C
  790 CONTINUE
C
C
C
C              KdEj   KfBj
C
C     ABAB = - ABAB * ABAB
C
C     Read and sort y3 ring intermediates.
C
C     We read KdEj as (E,K,d,j). Sort to EdKj.
C
      LEN = 0
      DO  810 IRREP=1,NIRREP
      LEN = LEN + IRPDPD(IRREP,9) * IRPDPD(IRREP,10)
  810 CONTINUE
C
      I000 = 1
      I010 = I000 + LEN
      I020 = I010 + LEN
      I030 = I020 +   2 * (NOCA * NOCB + NVRTA * NVRTB)
C
      NEED = IINTFP * I030
      IF(NEED.GT.MAXCOR)THEN
      WRITE(6,1020) NEED,MAXCOR
      STOP 'Y10C'
      ENDIF
C
C     Get the whole list at I010.
C
      CALL GETALL(CORE(I010),LEN,1,56)
C
C     Sort to different order : EdKj/BfKj is at I000.
C
      CALL SSTGEN(CORE(I010),CORE(I000),LEN,
     1        VRT(1,1),POP(1,1),VRT(1,2),POP(1,2),
     1        CORE(I020),1,'1324')
C
C     We read KfdJ as (f,K,d,J). Sort to fdKJ.
C
C      LEN = 0
C      DO  620 IRREP=1,NIRREP
C      LEN = LEN + IRPDPD(IRREP,12) * IRPDPD(IRREP,12)
C  620 CONTINUE
CC
C      I020 = I010 + LEN
C      I030 = I020 + LEN
C      I040 = I030 +   2 * (NOCA * NOCB + NVRTA * NVRTB)
CC
C      NEED = IINTFP * I040
C      IF(NEED.GT.MAXCOR)THEN
C      WRITE(6,1020) NEED,MAXCOR
C      STOP 'Y10C'
C      ENDIF
C
C     Get the whole list at I020.
C
C      CALL GETALL(CORE(I020),LEN,1,58+IUHF)
CC
CC     Sort to different order : fdKJ is at I010.
CC
C      CALL SSTGEN(CORE(I020),CORE(I010),LEN,
C     1        VRT(1,2),POP(1,1),VRT(1,2),POP(1,1),
C     1        CORE(I030),1,'1324')
C     1
C
C     Note that y3 has same symmetry as integrals at present. this is
C     exploited here, but it should be remembered that this may not be
C     in the case of Q (rather than Q*) corrections.
C
C     BDEF = KBEJ * JDFK - KDEJ * JBFK
C                   KFDJ          KFBJ
C
C     KBEJ was             (E,K,B,J)   now is (E,B,K,J)
C     JDFK was (F,J,D,K) = (D,K,F,J)   now is (D,F,K,J)
C     KDEJ was             (E,K,D,J)   now is (E,D,K,J)
C     JBFK was (F,J,B,K) = (B,K,F,J)   now is (B,F,K,J)
C
C     Therefore we are trying to form
C
C     W(B,D,E,F) = y3(E,B,KJ) * y3(D,F,KJ) - y3(E,D,KJ) * y3(B,F,KJ)
C
C     Somewhat crude, in core implementation as a trial (assumes that
C     all symmetry allowed (B,D,E,F) can be held in core).
C
      LNABCD = 0
      DO  840 IRREP=1,NIRREP
      NDIS(IRREP)   = 0
      DISSIZ(IRREP) = IRPDPD(IRREP,13)
      DO  830 JRREP=1,NIRREP
      KRREP = DIRPRD(IRREP,JRREP)
      NDIS(IRREP)   = NDIS(IRREP)
     1                 + POP(JRREP,1) * POP(KRREP,2)
  830 CONTINUE
      LNABCD = LNABCD + IRPDPD(IRREP,13) * IRPDPD(IRREP,13)
  840 CONTINUE
C
      I020 = I010 + LNABCD
C
      NEED = IINTFP * I020
      IF(NEED.GT.MAXCOR)THEN
      WRITE(6,1020) NEED,MAXCOR
      STOP 'Y10C'
      ENDIF
C
      IOFFY3 = I000
C     IOFFY3B = I010
      DO  850 IRPJK=1,NIRREP
C
      IF(IRPJK.EQ.1)THEN
      IOFFW(IRPJK) = I010
      ELSE
      IOFFW(IRPJK) = IOFFW(IRPJK-1) + DISSIZ(IRPJK-1)*DISSIZ(IRPJK-1)
      ENDIF
C
      CALL XGEMM('N','T',DISSIZ(IRPJK),DISSIZ(IRPJK),NDIS(IRPJK),
     1            -4.0D+00,
     1            CORE(IOFFY3),DISSIZ(IRPJK),
     1            CORE(IOFFY3),DISSIZ(IRPJK),0.0D+00,
     1            CORE(IOFFW(IRPJK)),DISSIZ(IRPJK))
C
      IOFFY3 = IOFFY3 + DISSIZ(IRPJK) *   NDIS(IRPJK)
C     IOFFY3B = IOFFY3B + DISSIZB(IRPJK) *   NDIS(IRPJK)
C
  850 CONTINUE
C
      DO  990 IRPEF=1,NIRREP
      IRPBD = IRPEF
C
C     Get a symmetry block of Y10(BD,EF) at I020.
C
      I030 = I020 + IRPDPD(IRPEF,13) * IRPDPD(IRPEF,13)
      NEED = IINTFP * I030
      IF(NEED.GT.MAXCOR)THEN
      WRITE(6,1020) NEED,MAXCOR
      STOP 'Y10C'
      ENDIF
C
cYAU - anchor
      CALL GETLIST(CORE(I020),1,IRPDPD(IRPEF,13),1,IRPEF,6)
      DO IRPF = 1,NIRREP
         IRPE = DIRPRD(IRPEF,IRPF)
         MAXE = VRT(IRPE,1)
         MAXF = VRT(IRPF,2)
         DO IRPD = 1,NIRREP
            IRPB = DIRPRD(IRPBD,IRPD)
            MAXB = VRT(IRPB,1)
            MAXD = VRT(IRPD,2)
            IRPED = DIRPRD(IRPE,IRPD)
            IRPBF = DIRPRD(IRPB,IRPF)
            DO F = 0,MAXF-1
               DO E = 0,MAXE-1
                  DO D = 0,MAXD-1
                     DO B = 0,MAXB-1
      ED = IOFFVV(IRPD,IRPED,5) + D*MAXE + E
      BF = IOFFVV(IRPF,IRPBF,5) + F*MAXB + B
      EDBF = IOFFW(IRPED) + BF*DISSIZ(IRPBF) + ED
      EF = IOFFVV(IRPF,IRPEF,5) + F*MAXE + E
      BD = IOFFVV(IRPD,IRPBD,5) + D*MAXB + B
      BDEF = EF*IRPDPD(IRPBD,13) + BD
      CORE(I020 + BDEF) = CORE(I020 + BDEF) - CORE(EDBF)
c                    END DO B = 0,MAXB-1
                     END DO
c                 END DO D = 0,MAXD-1
                  END DO
c              END DO E = 0,MAXE-1
               END DO
c           END DO F = 0,MAXF-1
            END DO
c        END DO IRPD = 1,NIRREP
         END DO
c     END DO IRPF = 1,NIRREP
      END DO
      CALL PUTLIST(CORE(I020),1,IRPDPD(IRPEF,13),1,IRPEF,6)
C
  990 CONTINUE
      RETURN
      END


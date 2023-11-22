       SUBROUTINE T2IN_DENSVV(DVV,T2AMPS,ICORE,MAXCOR,T2LN,T2LN_AA,
     +                        T2LN_BB,T2LN_AB,IRREPX,IUHF)
C
C THE FORMULAS ARE
C D(A,B) = 1/2 SUM M,N,E T[1](MN,AE) T[1](MN,BE)              (MBPT(2))
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DIRPRD,DISSYT,POP,VRT
      INTEGER T2LN,T2LN_AA,T2LN_BB,T2LN_AB,T2LEN
      INTEGER T2OFF_AA,T2OFF_AB
C
      COMMON/SYM/POP(8,2),VRT(8,2),NTAA,NTBB,NF1(2),NF2(2)
      COMMON/INFO/NOCCO(2),NVRTO(2)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/FLAGS/IFLAGS(100)
      COMMON/SYMINF/NSTART,NIRREP,IRREPA(255),IRREPB(255),DIRPRD(8,8)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),NTOT(18)

      DIMENSION ICORE(MAXCOR)
      DIMENSION DVV(NF2(1)+IUHF*NF2(2))
      DIMENSION T2AMPS(T2LN)
C
      DATA ONE,TWO,ONEM,HALF /1.0D0,2.0D0,-1.0D0,0.5D0/
C
C LOOP OVER SPIN CASES (ISPIN=1: ALPHA DENSITY MATRIX, =2: BETA DENSITY MATRIX)
C
      T2OFF_AA = 1
      DO ISPIN=1,IUHF+1
C
         IF(ISPIN.EQ.1) THEN
          IOFF=1
          JOFF=0
         ELSE
           IOFF=NF2(1)+1
           JOFF=0
         ENDIF

         IF (IUHF.EQ.1) THEN
C
            FACT     = ONE
            T2OFF_AA = T2OFF_AA + (ISPIN-1)*T2LN_AA

            DO IRREPR=1,NIRREP

               IRREPL = DIRPRD(IRREPR,IRREPX)

               NVRT2SQ=0
               DO IRREPJ=1,NIRREP
                  NVRT2SQ=NVRT2SQ+VRT(IRREPJ,ISPIN)*
     &                    VRT(DIRPRD(IRREPJ,IRREPR),ISPIN)
               ENDDO
C
               DISSYT=IRPDPD(IRREPR,ISPIN)
               NUMSYT=IRPDPD(IRREPL,2+ISPIN) 
               T2LEN = DISSYT*NUMSYT

               I000=1
               I010=I000+IINTFP*NUMSYT*NVRT2SQ
               I020=I010+IINTFP*NUMSYT*NVRT2SQ
               IEND=I020+IINTFP*NUMSYT*DISSYT
               IF (IEND .GE. MAXCOR) CALL INSMEM("densvv_mp2",
     &                                           IEND,MAXCOR)
               IF (MIN(NUMSYT,DISSYT).NE.0) THEN

                  CALL DVVAA_MP2(T2AMPS(T2OFF_AA+JOFF),ICORE(I000),
     &                           ICORE(I010),ICORE(I020),DVV(IOFF),
     &                           FACT,ISPIN,POP(1,ISPIN),VRT(1,ISPIN),
     &                           DISSYT,NUMSYT,NVRT2SQ,T2LEN,
     &                           IRREPR,IUHF)
               ENDIF
               JOFF = JOFF + DISSYT*NUMSYT
            ENDDO
         ENDIF
#ifdef _DEBUG_LVL0
        Call checksum("T(im,ef)*L(jm,ef)->DVV:",DVV(1+(ISPIN-1)*
     &                 NF2(ISPIN)),NF2(ISPIN))
#endif

      FACT     = ONE
      T2OFF_AB = T2LN_AA + T2LN_BB + 1
      JOFF     = 0

      DO IRREPR=1,NIRREP
C
         IRREPL = DIRPRD(IRREPR,IRREPX)

         DISSYT=IRPDPD(IRREPL,13)
         NUMSYT=IRPDPD(IRREPR,14)
         T2LEN = DISSYT*NUMSYT

         I000=1
         I010=I000+IINTFP*NUMSYT*DISSYT
         I020=I010+IINTFP*NUMSYT*DISSYT
C
         IF (MIN(NUMSYT,DISSYT).NE.0) THEN

           IEND=I020+3*IINTFP*MAX(DISSYT,NUMSYT)

           IF (IEND.GE.MAXCOR) CALL INSMEM("densvv_mp2",
     &                                      IEND,MAXCOR)

             CALL DVVAB_MP2(T2AMPS(T2OFF_AB+JOFF),ICORE(I000),
     &                      ICORE(I010),DVV(IOFF),FACT,ISPIN,
     &                      POP(1,ISPIN),POP(1,3-ISPIN),VRT(1,ISPIN),
     &                      VRT(1,3-ISPIN),DISSYT,NUMSYT,
     &                      T2LEN,IRREPL,ICORE(I020),IUHF)
         ENDIF
         JOFF = JOFF + DISSYT*NUMSYT
      ENDDO
#ifdef _DEBUG_LVL0
        Call checksum("T(im,ef)*L(jm,ef)->DVV:",DVV(1+(ISPIN-1)*
     &                 NF2(ISPIN)),NF2(ISPIN))
#endif
      ENDDO

      RETURN
      END

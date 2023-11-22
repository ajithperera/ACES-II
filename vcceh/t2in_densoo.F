       SUBROUTINE T2IN_DENSOO(DOO,T2AMPS,ICORE,MAXCOR,T2LN,T2LN_AA,
     +                        T2LN_BB,T2LN_AB,IRREPX,IUHF)     
C
C THE GENERAL FORMULA IS
C
C D(I,J) = - 1/2 SUM M,E,F T[1](IM,EF) T[1](JM,EF)      (MBPT(2))
C
C THE RHF PART USES EXPLICITELY SPIN ADAPTED CODE 
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DIRPRD,DISSYT,POP,VRT
      INTEGER T2LN,T2LN_AA,T2LN_BB,T2LN_AB
      INTEGER T2OFF_AA,T2OFF_AB,T2LEN
C
      COMMON/SYM/POP(8,2),VRT(8,2),NTAA,NTBB,NF1(2),NF2(2)
      COMMON/INFO/NOCCO(2),NVRTO(2)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/FLAGS/IFLAGS(100)
      COMMON/SYMINF/NSTART,NIRREP,IRREPA(255),IRREPB(255),DIRPRD(8,8)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),NTOT(18)

      DIMENSION ICORE(MAXCOR),DOO(NF1(1)+IUHF*NF1(2))
      DIMENSION T2AMPS(T2LN)
C
      DATA ONE,ONEM,TWOM,HALF /1.0D0,-1.0D0,-2.0D0,0.5D0/
C
      T2OFF_AA = 1
      DO ISPIN=1,IUHF+1 
C
         IF(ISPIN.EQ.1) THEN
          IOFF=1
          JOFF=0
         ELSE
          IOFF=NF1(1)+1
          JOFF=0
         ENDIF

         T2OFF_AA = T2OFF_AA + (ISPIN-1)*T2LN_AA

         IF(IUHF.EQ.1) THEN

            FACT=ONEM

            DO IRREPR=1,NIRREP
C
               IRREPL = DIRPRD(IRREPR,IRREPX)

               NOCC2SQ=0 
               DO IRREPJ=1,NIRREP
                  NOCC2SQ=NOCC2SQ+POP(IRREPJ,ISPIN)
     &                   *POP(DIRPRD(IRREPJ,IRREPR),ISPIN)
               ENDDO 
 
               DISSYT=IRPDPD(IRREPL,ISPIN)
               NUMSYT=IRPDPD(IRREPR,2+ISPIN)  
               T2LEN = DISSYT*NUMSYT

               I000=1
               I010=I000+IINTFP*NOCC2SQ*DISSYT
               IEND=I010+IINTFP*NOCC2SQ*DISSYT

               IF(MIN(NUMSYT,DISSYT).NE.0) THEN

               IF(IEND .GT. MAXCOR) CALL INSMEM('densoo_mp2',IEND,
     &                                           MAXCOR) 
                 CALL DOOAA_MP2(T2AMPS(T2OFF_AA+JOFF),ICORE(I000),
     &                          ICORE(I010),DOO(IOFF),FACT,
     &                          ISPIN,POP(1,ISPIN),VRT(1,ISPIN),
     &                          DISSYT,NUMSYT,NOCC2SQ,T2LEN,IRREPR,
     &                          IUHF)
               ENDIF 
               JOFF = JOFF + DISSYT*NUMSYT
            ENDDO 
         ENDIF

#ifdef _DEBUG_LVL0
        Call checksum("T(mn,ae)*L(mn,be)->DOO:",DOO(1+(ISPIN-1)*
     &               NF1(ISPIN)),NF1(ISPIN))
#endif
         FACT     = ONEM
         T2OFF_AB = T2LN_AA + T2LN_BB + 1
         T2LEN    = T2LN_AB
         JOFF     = 0

         DO IRREPR=1,NIRREP
C
            IRREPL = DIRPRD(IRREPR,IRREPX)

            DISSYT=IRPDPD(IRREPL,13)
            NUMSYT=IRPDPD(IRREPR,14)
            T2LEN = DISSYT*NUMSYT

            I000 = 1
            I010 = I000 + NUMSYT*DISSYT*IINTFP
            I020 = I010 + NUMSYT*DISSYT*IINTFP

            IF (MIN(NUMSYT,DISSYT).NE.0) THEN
               IEND=I020 + 3*IINTFP*MAX(NUMSYT,DISSYT)
               IF(IEND .GT. MAXCOR) CALL INSMEM('densoo_mp2',IEND,
     &                                           MAXCOR)

               CALL DOOAB_MP2(T2AMPS(T2OFF_AB+JOFF),ICORE(I000),
     &                        ICORE(I010),DOO(IOFF),FACT,ISPIN,
     &                        POP(1,ISPIN),POP(1,3-ISPIN),VRT(1,ISPIN),
     &                        VRT(1,3-ISPIN),DISSYT,NUMSYT,T2LEN,
     &                        IRREPR,ICORE(I020),IUHF)
            ENDIF
            JOFF = JOFF + DISSYT*NUMSYT
         ENDDO 

#ifdef _DEBUG_LVL0
        Call checksum("T(mn,ae)*L(mn,be)->DOO:",DOO(1+(ISPIN-1)*
     &               NF1(ISPIN)),NF1(ISPIN))
#endif
      ENDDO
C
      RETURN
      END

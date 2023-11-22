       SUBROUTINE T1IN_DENSOO(DOO,T1AMPS,WORK,MAXCOR,T1LN,T1LN_AA,
     +                        T1LN_BB,IRREPX,IUHF)     
C
C THE GENERAL FORMULA IS
C
C D(I,J) = - 1/2 SUM M,E,F T[1](A,I) T[1](A,J)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DIRPRD,DISSYT,POP,VRT
      INTEGER T1LN,T1LN_AA,T1LN_BB
      INTEGER T1OFF_AA
C
      COMMON/SYM/POP(8,2),VRT(8,2),NTAA,NTBB,NF1(2),NF2(2)
      COMMON/INFO/NOCCO(2),NVRTO(2)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/FLAGS/IFLAGS(100)
      COMMON/SYMINF/NSTART,NIRREP,IRREPA(255),IRREPB(255),DIRPRD(8,8)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),NTOT(18)

      DIMENSION WORK(MAXCOR),DOO(NF1(1)+IUHF*NF1(2))
      DIMENSION T1AMPS(T1LN)
C
      DATA ONE,ONEM,TWOM,HALF /1.0D0,-1.0D0,-2.0D0,0.5D0/
C
      T1OFF_AA = 1
      DO ISPIN=1, IUHF+1 
C
         T1OFF_AA = T1OFF_AA + (ISPIN-1)*T1LN_AA

        IF (ISPIN .EQ. 1) THEN
            IOFF = 0
            JOFF = 1
         ELSE
            IOFF = 0
            JOFF = NF1(1) + 1
         ENDIF 

         FACT=ONEM

         DO IRREPR=1,NIRREP

             IRREPL = DIRPRD(IRREPR,IRREPX) 

             NOCC = POP(IRREPR,ISPIN)
             NVRT = VRT(IRREPL,ISPIN)

             IF (MIN(NVRT,NOCC).GT.0) THEN

               CALL XGEMM('T','N',NOCC,NOCC,NVRT,FACT,
     &                     T1AMPS(T1OFF_AA+IOFF),NVRT,
     &                     T1AMPS(T1OFF_AA+IOFF),NVRT,
     &                     ONE,DOO(JOFF),NOCC)
             ENDIF 
             IOFF = IOFF + NOCC*NVRT
             JOFF = JOFF + NOCC*NOCC
         ENDDO 

#ifdef _DEBUG_LVL0
        Call checksum("T(a,i)*L(a,j)->DOO    :",DOO(1+(ISPIN-1)*
     &               NF1(ISPIN)),NF1(ISPIN))
#endif
      ENDDO

      RETURN
      END

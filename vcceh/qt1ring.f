C
      SUBROUTINE QT1RING(ICORE, MAXCOR, IRREPX, IOFFSET, IUHF)
C
C Driver for G(MB,EJ) <-- T1 contributions. Selects between incore
C and out of core algorithms and call the appropriate routines.
C Spin orbital expression for the term evaluated by dependencies
C of this routine can be written as,
C
C   G(MB,EJ) = SUM T(J,F)*Hbar(MB,EF) - SUM T(N,B)*Hbar(MN,EJ) 
C 
C Individual spin integrated formulas (6 of them) can be found in
C each subroutine which evaluate each contribution.
   
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE, ONEM, ZILCH
      DIMENSION ICORE(MAXCOR)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYM/ POP(8,2),VRT(8,2),NT(2),NFEA(2),NFMI(2)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /FLAGS/ IFLAGS(100)
C
C Common blocks used in the quadratic term
C
      COMMON /QADINTI/ INTIMI, INTIAE, INTIME
      COMMON /QADINTG/ INGMNAA, INGMNBB, INGMNAB, INGABAA,
     &                 INGABBB, INGABAB, INGMCAA, INGMCBB,
     &                 INGMCABAB, INGMCBABA, INGMCBAAB,
     &                 INGMCABBA
      COMMON /APERTT1/ IAPRT1AA
      COMMON /BPERTT1/ IBPRT1AA
      COMMON /APERTT2/ IAPRT2AA1, IAPRT2BB1, IAPRT2AB1, IAPRT2AB2, 
     &                 IAPRT2AB3, IAPRT2AB4, IAPRT2AA2, IAPRT2BB2,
     &                 IAPRT2AB5 
C
      DATA ONE  /1.0D+00/   
      DATA ONEM /-1.0D+00/
      DATA ZILCH /0.0D+00/
C
C Spin cases AAAA and BBBB
C
      INEED = -1 
C
      DO 10 ISPIN = 1, 1 + IUHF
         LSTGTL   = 8 + ISPIN
         LSTGTR   = 8 + ISPIN
         LSTWINT1 = 6 + ISPIN
         LSTWINT2 = 26 + ISPIN
         FULTYP1  = 18 + ISPIN
         FULTYP2  = 20 + ISPIN
         T1SIZE   = IRPDPD(IRREPX, 9) + IUHF*IRPDPD(IRREPX, 10)
         ABFULL   = IRPDPD(IRREPX, FULTYP1)
         MNFULL   = IRPDPD(IRREPX, FULTYP2)
C
         DO 20 IRREPBM = 1, NIRREP
            IRREPWBM = IRREPBM
            IRREPWEF = IRREPWBM            
            IRREPGBM = IRREPWBM
            IRREPGEJ = DIRPRD(IRREPGBM, IRREPX)
            IRREPWEJ = IRREPGEJ
            IRREPWMN = IRREPWEJ
C
            TARDSZ = IRPDPD(IRREPGBM, LSTGTL)
            TARDIS = IRPDPD(IRREPGEJ, LSTGTR)
            TARSIZ = TARDSZ*TARDIS
C
            FULDSZ1 = IRPDPD(IRREPWEF, FULTYP1)
            FULDSZ2 = IRPDPD(IRREPWMN, FULTYP2)
C
            INTDSZ2 = IRPDPD(IRREPWMN, ISYTYP(1, LSTWINT1))
            INTDIS2 = IRPDPD(IRREPWEJ, ISYTYP(2, LSTWINT1))
            INTDSZ1 = IRPDPD(IRREPWEF, ISYTYP(1, LSTWINT2))
            INTDIS1 = IRPDPD(IRREPWBM, ISYTYP(2, LSTWINT2))
C
 1          I1 = TARDIS*TARSIZ + T1SIZE
            I2 = 3*MAX(TARDSZ, TARDIS)
C     
            I3A = FULDSZ1*INTDIS1 + INTDSZ1
            I3B = FULDSZ2*INTDIS2 + INTDSZ2
C     
            I4A = 2*ABFULL
            I4B = 2*MNFULL
C
            INEEDA = I1 + MAX(I2, I3A, I4A)
            INEEDB = I1 + MAX(I2, I3B, I4B)
            INEED = MAX(INEED, INEEDA, INEEDB)
 20      CONTINUE
 10   CONTINUE
C
      TOTSIZ = INEED*IINTFP
C 
c      IF (TOTSIZ .LT. MAXCOR) THEN
      IF (INEED .LT. MAXCOR/IINTFP) THEN
         IF (IFLAGS(1) .GE. 20) THEN
            WRITE(LUOUT, *)
            WRITE(LUOUT, 1000) 
         ENDIF
         CALL MKT1GMBEJA1(ICORE, MAXCOR, IUHF, IRREPX, IOFFSET)
      ELSE
         IF (IFLAGS(1) .GE. 20) THEN
            WRITE(LUOUT, *)
            WRITE(LUOUT, 2000) 
         ENDIF
         CALL MKT1GMBEJA2(ICORE, MAXCOR, IUHF, IRREPX, IOFFSET)
      ENDIF
C
C Spin cases BAAB and ABBA 
C     
      INEED = -1
      DO 110 ISPIN = 1, 1 + IUHF
         LSTWINT1 = 31 - ISPIN
         LSTWINT2 = 8 + ISPIN + (1 - IUHF)
         FULTYP1  = 18 + ISPIN
         FULTYP2  = 20 + ISPIN
         T1SIZE  = IRPDPD(IRREPX, 9) + IUHF*IRPDPD(IRREPX, 10)
         ABFULL  = IRPDPD(IRREPX, FULTYP1)
C         
         DO 120 IRREPBM = 1, NIRREP
            IRREPWBM = IRREPBM
            IRREPWEF = IRREPWBM
            IRREPGBM = IRREPWBM
            IRREPGEJ = DIRPRD(IRREPGBM, IRREPX) 
            IRREPWEJ = IRREPGEJ
            IRREPWMN = IRREPWEJ
C
            IF (ISPIN .EQ. 1) THEN
               TARDSZ = IRPDPD(IRREPGBM, 11)
               TARDIS = IRPDPD(IRREPGEJ, 11)
            ELSE
               TARDSZ = IRPDPD(IRREPGBM, 12)
               TARDIS = IRPDPD(IRREPGEJ, 12)
            ENDIF
            TARSIZ = TARDSZ*TARDIS
C
            INTDSZ1 = IRPDPD(IRREPWEF, ISYTYP(1, LSTWINT1))
            INTDIS1 = IRPDPD(IRREPWBM, ISYTYP(2, LSTWINT1))
            INTDSZ2 = IRPDPD(IRREPWMN, ISYTYP(1, LSTWINT2))
            INTDIS2 = IRPDPD(IRREPWEJ, ISYTYP(2, LSTWINT2))
C            	
            I1 = TARDIS*TARSIZ + T1SIZE
            I2 = 3*MAX(TARDSZ, TARDIS)
C     
            I3A = INTDSZ1*INTDIS1 + 3*MAX(INTDSZ1, INTDIS1)
            I3B = INTDSZ2*INTDIS2 + 3*MAX(INTDSZ2, INTDIS2)
C     
            I4 = 2*ABFULL
C     
            I5 = TARDIS*TARDSZ + 3*TARDSZ
C
            INEEDA = I1 + MAX(I2, I3A, I4, I5)
            INEEDB = I1 + MAX(I2, I3B, I4, I5)
            INEED  = MAX(INEED, INEEDA, INEEDB)
 120     CONTINUE
 110  CONTINUE
C
      TOTSIZ = INEED*IINTFP
C
c      IF (TOTSIZ .LT. MAXCOR) THEN
      IF (INEED .LT. MAXCOR/IINTFP) THEN
         IF (IFLAGS(1) .GE. 20) THEN
            WRITE(LUOUT, *)
            WRITE(LUOUT, 1000)
         ENDIF
         CALL MKT1GMBEJB1(ICORE, MAXCOR, IUHF, IRREPX, IOFFSET)
      ELSE
         IF (IFLAGS(1) .GE. 20) THEN
            WRITE(LUOUT, *)
            WRITE(LUOUT, 2000)
         ENDIF
         CALL MKT1GMBEJB2(ICORE, MAXCOR, IUHF, IRREPX, IOFFSET)
      ENDIF
C     
C SPIN CASES ABAB AND BABA (FOR UHF ONLY)
C     
      INEED = -1
      DO 210 ISPIN = 1, 2
         LSTTAR   = (INGMCABAB - 1) + ISPIN
         LSTTMP1  = 37 + ISPIN
         LSTTMP2  = 40 - ISPIN
         LSTWINT1 = 28 + ISPIN
         LSTWINT2 = 8 +  ISPIN
         T1SIZE   =  IRPDPD(IRREPX, 9) + IUHF*IRPDPD(IRREPX, 10)
C
         DO 220 IRREPBM = 1, NIRREP
            IRREPWBM = IRREPBM
            IRREPWEF = IRREPWBM
            IRREPGBM = IRREPWBM
            IRREPGEJ = DIRPRD(IRREPGBM, IRREPX)
            IRREPWEJ = IRREPGEJ
            IRREPWMN = IRREPWEJ
C
            IRREPGBJ = IRREPBM
            IRREPGEM = DIRPRD(IRREPGBJ, IRREPX)
C
            TARDSZ = IRPDPD(IRREPGEM, ISYTYP(1, LSTTAR))
            TARDIS = IRPDPD(IRREPGBJ, ISYTYP(2, LSTTAR))
            TARSIZ = TARDSZ*TARDIS
C
            TMPDSZ1 = IRPDPD(IRREPGBM, ISYTYP(1, LSTTMP1))
            TMPDIS1 = IRPDPD(IRREPGEJ, ISYTYP(2, LSTTMP1))
            TMPDSZ2 = IRPDPD(IRREPGEJ, ISYTYP(1, LSTTMP2))
            TMPDIS2 = IRPDPD(IRREPGBM, ISYTYP(2, LSTTMP2))
C     
            INTDSZ1 = IRPDPD(IRREPWEF, ISYTYP(1, LSTWINT1))
            INTDIS1 = IRPDPD(IRREPWBM, ISYTYP(2, LSTWINT1))
            INTDSZ2 = IRPDPD(IRREPWMN, ISYTYP(1, LSTWINT2))
            INTDIS2 = IRPDPD(IRREPWEJ, ISYTYP(2, LSTWINT2))
C  
            I1 = TARDIS*TARSIZ + T1SIZE
            I2 = 3*MAX(TARDSZ, TARDIS)
C     
            I3A = INTDSZ1*INTDIS1 + 3*MAX(INTDSZ1, INTDIS1)
            I3B = INTDSZ2*INTDIS2 + 3*MAX(INTDSZ2, INTDIS2)
C     
            I4 = TARDIS*TARDSZ + 3*TARDSZ
C     
            I5A = 3*MAX(TMPDSZ1, TMPDIS1)
            I5B = 3*MAX(TMPDSZ2, TMPDIS2)

            INEEDA = I1 + MAX(I2, I3A, I4, I5A)
            INEEDB = I1 + MAX(I2, I3B, I4, I5B)
            INEED = MAX(INEED, INEEDA, INEEDB)
 220     CONTINUE
 210  CONTINUE
C      
      TOTSIZ = INEED*IINTFP      
C
c      IF(TOTSIZ .LT. MAXCOR) THEN
      IF (INEED .LT. MAXCOR/IINTFP) THEN
         IF (IFLAGS(1) .GE. 20) THEN
            WRITE(LUOUT, *)
            WRITE(LUOUT, 1000)
         ENDIF
         CALL MKT1GMBEJC1(ICORE, MAXCOR, IUHF, IRREPX, IOFFSET)
      ELSE
         IF (IFLAGS(1) .GE. 20) THEN
            WRITE(LUOUT, *)
            WRITE(LUOUT, 2000)
         ENDIF
         CALL MKT1GMBEJC2(ICORE, MAXCOR, IUHF, IRREPX, IOFFSET)
      ENDIF
C
1000  FORMAT(T3, '@-QT1RING-I, T1-Ring contributions by',
     &           ' in-core algorithm.')
2000  FORMAT(T3, '@-QT1RING-I, T1-Ring contributions by',
     &           ' out-of-core algorithm.')
C
      RETURN
      END

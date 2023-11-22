
C     Eventually supposed to become a decent program for noniterative T4
C     corrections to CCSD and other methods. At the moment it is a piece
C     of junk (and that is being nice).

      PROGRAM T4
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER POP1,POP2,VRT1,VRT2,DIRPRD,DISTSZ
      INTEGER RLECYC
      LOGICAL MBPT3,MBPT4,CC,TRPEND,SNGEND,GRAD,MBPTT,SING1,QCISD,LINCC
      LOGICAL CICALC,NONHF,TRIPIT,TRIPNI,TRIPNI1,UCC,RESTART,BRUECK
      LOGICAL ROHF4,ITRFLG
      DIMENSION ECORR(3),E(10)
      DIMENSION E45T(20)
      COMMON / / ICORE(1)
      COMMON /ISTART/ I0,ICRSIZ
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /FILES/ LUOUT,MOINTS
      COMMON /INFO/ NOCCO(2),NVRTO(2)
      COMMON /CONTROL/ IPRNT,IXXX,IXXX2
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /SWITCH/ MBPT3,MBPT4,CC,TRPEND,SNGEND,GRAD,MBPTT,SING1,
     &                QCISD,UCC
      COMMON/TRIPLES/TRIPNI,TRIPNI1,TRIPIT
      COMMON /LINEAR / LINCC,CICALC
C
C HP is giving problems with this common blok. However,
C this is not being used in the program. Ajith 05/01/97
C
C      COMMON /NHFREF/ NONHF
      COMMON /ROHF/ ROHF4,ITRFLG
      COMMON /CORENG/ ELAST
      EQUIVALENCE (IFLAGS(1),IPRINT)
      EQUIVALENCE (IFLAGS(2),METHOD)
C      
      COMMON /AUXIO/ DISTSZ(8,100),NDISTS(8,100),INIWRD(8,100),LNPHYR,
     1               LUAUX
C
C
C         POP1   ...... NUMBER OF OCCUPIED ALPHA ORBITALS WITHIN EACH IRREP
C         POP2   ...... NUMBER OF OCCUPIED BETA ORBITALS WITIN EACH IRREP
C         VRT1   ...... NUMBER OF VIRTUAL ALPHA ORBITALS WITHIN EACH IRREP
C         VRT2   ...... NUMBER OF VIRTUAL BETA ORBITALS WITHIN EACH IRREP
C         NTAA   ...... LENGTH OF T1(A,I) ALPHA
C         NTBB   ...... LENGTH OF T1(A,I) BETA
C         NF1AA  ...... LENGTH OF F(M,I) ALPHA ( NOTE ALL Fs ARE NOT SYMMETRIC)
C         NF1BB  ...... LENGTH OF F(M,I) BETA
C         NF2AA  ...... LENGTH OF F(A,E) ALPHA
C         NF2BB  ...... LENGTH OF F(A,E) BETA
C
      COMMON/SYM/POP1(8),POP2(8),VRT1(8),VRT2(8),
     &           NTAA,NTBB,NF1AA,NF1BB,NF2AA,NF2BB
      COMMON/SYMINF/NSTART,NIRREP,IRREPA(255),IRREPB(255),DIRPRD(8,8)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),NJUNK(18)
C
      COMMON /T3OFF/  IOFFVV(8,8,10),IOFFOO(8,8,10),IOFFVO(8,8,4)
C
      CALL CRAPSI(ICORE,IUHF,0)
C
      MAXCOR=ICRSIZ
C
C     Set this so that particle-particle ladder contractions are handled
C     correctly for RHF.
C
      GRAD =.TRUE.
C
      WRITE(6,1000)
 1000 FORMAT(/,72('-'),/,T10,
     1       ' Noniterative 5th Order Quadruples Calculation (E5QQ) ',
     1       /,72('-'),/)
C
C     Initialize auxiliary i/o stuff.
C
      CALL AUXIOI
C
C     Initialize offset arrays.
C
      CALL MKOFOO
      CALL MKOFVV
      CALL MKOFVO
C
C     Initialize lists for y1 and y2 quantities. These reside in the space
C     formerly occupied by the F intermediates.
C
        CALL ZERLST(ICORE(I0),NF1AA,1,1,1,91)
        CALL ZERLST(ICORE(I0),NF2AA,1,1,1,92)
C       CALL ZERLST(ICORE(I0),NT1AA,1,1,1,93)
        IF(IUHF.NE.0)THEN
         CALL ZERLST(ICORE(I0),NF1BB,1,1,2,91)
         CALL ZERLST(ICORE(I0),NF2BB,1,1,2,92)
C        CALL ZERLST(ICORE(I0),NT1BB,1,1,2,93)
        ENDIF
        call zerlst(ICORE(I0),ntaa,1,1,1,90)
        call zerlst(ICORE(I0),ntaa,1,1,3,90)
        IF(IUHF.NE.0)THEN
        call zerlst(ICORE(I0),ntbb,1,1,2,90)
        call zerlst(ICORE(I0),ntbb,1,1,4,90)
        ENDIF
C
C     Form y1 and y2 and write to lists.
C
      INCRF = 0
      FACT = 1.0D+00
C
C     y1
      WRITE(6,1010)
 1010 FORMAT(' @T4-I, Computing Y1 intermediate. ')
      CALL QUAD3(ICORE(I0),MAXCOR,IUHF,INCRF,FACT)
C
C     y2
      WRITE(6,1020)
 1020 FORMAT(' @T4-I, Computing Y2 intermediate. ')
      CALL QUAD2(ICORE(I0),MAXCOR,IUHF,INCRF,FACT)
C
C     y3
C
C     Initialize ring intermediate lists.
C
      CALL ZERSYM(ICORE(I0),54)
      CALL ZERSYM(ICORE(I0),56)
      CALL ZERSYM(ICORE(I0),58)
      IF(IUHF.EQ.1)THEN
      CALL ZERSYM(ICORE(I0),55)
      CALL ZERSYM(ICORE(I0),57)
      CALL ZERSYM(ICORE(I0),59)
      ENDIF
       C1= 0.0D+00
       C2= 0.0D+00
       C3= 0.0D+00
       C4=-0.5D+00
      WRITE(6,1030)
 1030 FORMAT(' @T4-I, Computing Y3 intermediate. ')
      CALL DWMBEJ(ICORE(I0),MAXCOR,'ABAB',IUHF,C1,C2,C3,C4)
      CALL DWMBEJ(ICORE(I0),MAXCOR,'ABBA',IUHF,C1,C2,C3,C4)
      IF(IUHF.EQ.1)THEN
       CALL DWMBEJ(ICORE(I0),MAXCOR,'AAAA',IUHF,C1,C2,C3,C4)
       CALL DWMBEJ(ICORE(I0),MAXCOR,'BBBB',IUHF,C1,C2,C3,C4)
       CALL DWMBEJ(ICORE(I0),MAXCOR,'BABA',IUHF,C1,C2,C3,C4)
       CALL DWMBEJ(ICORE(I0),MAXCOR,'BAAB',IUHF,C1,C2,C3,C4)
      ENDIF
C
C     y6 (particle-particle ladder)
C
C      Initialize lists.
C
      CALL ZERSYM(ICORE(I0),48)
      IF(IUHF.EQ.1)THEN
      CALL ZERSYM(ICORE(I0),49)
      ENDIF
      CALL ZERSYM(ICORE(I0),50)
      ITYPE = 6
      WRITE(6,1040)
 1040 FORMAT(' @T4-I, Computing Y6 intermediate. ')
      CALL DRLAD(ICORE(I0),MAXCOR,IUHF,ITYPE)
C
C     y9 (f*t2 and ring contractions).
C
      CALL ZERSYM(ICORE(I0),61)
      IF(IUHF.EQ.1)THEN
      CALL ZERSYM(ICORE(I0),62)
      ENDIF
      CALL ZERSYM(ICORE(I0),63)
C
      WRITE(6,1050)
 1050 FORMAT(' @T4-I, Computing Y9 intermediate. ')
      CALL FEACONT(ICORE(I0),MAXCOR,IUHF)
      IF(IUHF.EQ.0)THEN
      IBOT=3
      ELSE
      IBOT=1
      ENDIF
      if(iuhf.gt.0)then
      call zersym(icore(i0),40)
      call zersym(icore(i0),41)
      endif
      call zersym(icore(i0),42)
      call zersym(icore(i0),43)
      DO 10 ISPIN=IBOT,3
      CALL DRRNG(ICORE(I0),MAXCOR,ISPIN,IUHF)
   10 CONTINUE
C
C     compute energy contribution from y6 and y9.
C
      WRITE(6,1060)
 1060 FORMAT(' @T4-I, Computing Y6-Y9 energy contribution. ')
C
      CALL ENABIJ(ICORE(I0),MAXCOR,IUHF,1,E,3,ESPAD)
      IF(IUHF.EQ.0)THEN
      E(2) = E(1)
      ENDIF
      IF(IUHF.EQ.0)THEN
      E69 = ESPAD
      ELSE
      E69 = E(1) + E(2) + E(3)
      ENDIF
      WRITE(6,2010) E(1),E(2),E(3),E69
 2010 FORMAT('        AAAA ',F20.12,/,
     1       '        BBBB ',F20.12,/,
     1       '        ABAB ',F20.12,/,
     1       '        Tot. ',F20.12)
C
C     y7 (hole-hole ladder; overwrites y6)
C
      CALL ZERSYM(ICORE(I0),48)
      IF(IUHF.EQ.1)THEN
      CALL ZERSYM(ICORE(I0),49)
      ENDIF
      CALL ZERSYM(ICORE(I0),50)
      ITYPE=1
      WRITE(6,1070)
 1070 FORMAT(' @T4-I, Computing Y7 intermediate. ')
      CALL DRLAD(ICORE(I0),MAXCOR,IUHF,ITYPE)
C
C     y11 (f*t2 and ring contractions; overwrites y9).
C
      CALL ZERSYM(ICORE(I0),61)
      IF(IUHF.EQ.1)THEN
      CALL ZERSYM(ICORE(I0),62)
      ENDIF
      CALL ZERSYM(ICORE(I0),63)
C
      WRITE(6,1080)
 1080 FORMAT(' @T4-I, Computing Y11 intermediate. ')
      CALL FMICONT(ICORE(I0),MAXCOR,IUHF)
      if(iuhf.gt.0)then
      call zersym(icore(i0),40)
      call zersym(icore(i0),41)
      endif
      call zersym(icore(i0),42)
      call zersym(icore(i0),43)
      DO 20 ISPIN=IBOT,3
      CALL DRRNG(ICORE(I0),MAXCOR,ISPIN,IUHF)
   20 CONTINUE
C
C     compute energy contribution from y7 and y11.
C
      WRITE(6,1090)
 1090 FORMAT(' @T4-I, Computing Y7-Y11 energy contribution. ')
      CALL ENABIJ(ICORE(I0),MAXCOR,IUHF,1,E,3,ESPAD)
      IF(IUHF.EQ.0)THEN
      E(2) = E(1)
      ENDIF
      IF(IUHF.EQ.0)THEN
      E711 = ESPAD
      ELSE
      E711 = E(1) + E(2) + E(3)
      ENDIF
      WRITE(6,2010) E(1),E(2),E(3),E711
C
C     Massage mbej lists.
C
      if(iuhf.eq.0) call reset(ICORE(I0),MAXCOR,IUHF)
      CALL ALTMBEJ(ICORE(I0),MAXCOR,IUHF)
C
C     Attempt to form Y4 and Y5.
C
      CALL DRY4Y5(ICORE(I0),MAXCOR,IUHF,1)
      CALL DRY4Y5(ICORE(I0),MAXCOR,IUHF,0)
C
C     attempt to compute Y10 and Y12.
C
      CALL Y10A(ICORE(I0),MAXCOR,IUHF)
      CALL Y10B(ICORE(I0),MAXCOR,IUHF)
      CALL Y10C(ICORE(I0),MAXCOR,IUHF)
      CALL Y12A(ICORE(I0),MAXCOR,IUHF)
      CALL Y12B(ICORE(I0),MAXCOR,IUHF)
      CALL Y12C(ICORE(I0),MAXCOR,IUHF)
C
      WRITE(6,2020)
 2020 FORMAT(' @T4-I, Computing <mn//jl>-Y12(mnjl) energy. ')
      E(1) = 0.0D+00
      E(2) = 0.0D+00
      E(3) = 0.0D+00
      ESPAD = 0.0D+00
      CALL ENABIJ(ICORE(I0),MAXCOR,IUHF,2,E,3,ESPAD)
      IF(IUHF.EQ.0)THEN
      E(2) = E(1)
      ENDIF
      IF(IUHF.EQ.0)THEN
      EIJKL = ESPAD
      ELSE
      EIJKL = E(1) + E(2) + E(3)
      ENDIF
      WRITE(6,2010) E(1),E(2),E(3),EIJKL
C
      WRITE(6,2030)
 2030 FORMAT(' @T4-I, Computing <bd//ef>-Y10(bdef) energy. ')
      E(1) = 0.0D+00
      E(2) = 0.0D+00
      E(3) = 0.0D+00
      ESPAD = 0.0D+00
      CALL ENABIJ(ICORE(I0),MAXCOR,IUHF,3,E,3,ESPAD)
      IF(IUHF.EQ.0)THEN
      E(2) = E(1)
      ENDIF
      IF(IUHF.EQ.0)THEN
      EABCD = ESPAD
      ELSE
      EABCD = E(1) + E(2) + E(3)
      ENDIF
      WRITE(6,2010) E(1),E(2),E(3),EABCD
C
      CALL   Y8(ICORE(I0),MAXCOR,IUHF)
      CALL Y13A(ICORE(I0),MAXCOR,IUHF)
      CALL Y13B(ICORE(I0),MAXCOR,IUHF)
      CALL Y13C(ICORE(I0),MAXCOR,IUHF)
      CALL Y13D(ICORE(I0),MAXCOR,IUHF)
C
      WRITE(6,2040)
 2040 FORMAT(' @T4-I, Computing Y8(cekm)-Y13(cekm) energy. ')
      E(1) = 0.0D+00
      E(2) = 0.0D+00
      E(3) = 0.0D+00
      E(4) = 0.0D+00
      E(5) = 0.0D+00
      E(6) = 0.0D+00
      CALL ENABIJ(ICORE(I0),MAXCOR,IUHF,4,E,6,ESPAD)
      IF(IUHF.EQ.0)THEN
      E(2) = E(1)
      E(4) = E(3)
      E(6) = E(5)
      ENDIF
      E813 = E(1) + E(2) + E(3) + E(4) + E(5) + E(6)
      WRITE(6,2050) E(1),E(2),E(3),E(4),E(5),E(6),E813
 2050 FORMAT('        AAAA ',F20.12,/,
     1       '        BBBB ',F20.12,/,
     1       '        ABAB ',F20.12,/,
     1       '        BABA ',F20.12,/,
     1       '        ABBA ',F20.12,/,
     1       '        BAAB ',F20.12,/,
     1       '        Tot. ',F20.12)
C
      CALL Y14A(ICORE(I0),MAXCOR,IUHF)
      CALL Y14B(ICORE(I0),MAXCOR,IUHF)
      CALL Y14C(ICORE(I0),MAXCOR,IUHF)
      CALL Y14D(ICORE(I0),MAXCOR,IUHF)
C
      WRITE(6,2060)
 2060 FORMAT(' @T4-I, Computing <mb//le>-Y14(b,e,l,m) energy. ')
      E(1) = 0.0D+00
      E(2) = 0.0D+00
      E(3) = 0.0D+00
      E(4) = 0.0D+00
      E(5) = 0.0D+00
      E(6) = 0.0D+00
      CALL ENABIJ(ICORE(I0),MAXCOR,IUHF,5,E,6,ESPAD)
      IF(IUHF.EQ.0)THEN
      E(2) = E(1)
      E(4) = E(3)
      E(6) = E(5)
      ENDIF
      E14 = E(1) + E(2) + E(3) + E(4) + E(5) + E(6)
      WRITE(6,2070) E(1),E(2),E(3),E(4),E(5),E(6),E14
 2070 FORMAT('        AAAA ',F20.12,/,
     1       '        BBBB ',F20.12,/,
     1       '        ABAB ',F20.12,/,
     1       '        BABA ',F20.12,/,
     1       '        AABB ',F20.12,/,
     1       '        BBAA ',F20.12,/,
     1       '        Tot. ',F20.12)
C
      E5QQ = E69 + E711 + EIJKL + EABCD + E813 + E14
      WRITE(6,3010) E5QQ
 3010 FORMAT(/,72('-'),/,T20,' @T4-I, E5QQ ',F20.12,/,72('-'),/)
C
C     --- Summarize the results of the whole monster calculation ---
C
      WRITE(LUOUT,3090)
 3090 FORMAT(/,'                    ******************** ',/,
     1         '                          Summary        ',/,
     1         '                    ******************** ',/)
      CALL GETREC(20,'JOBARC','NONITT3',20*IINTFP,E45T)
      E4T     = E45T(1)
      E5ST    = E45T(2)
      E4DT    = E45T(3)
      E5TD    = E45T(4)
      E5TT    = E45T(5)
      E5QT    = E45T(6)
      E5QTC   = E5TD
      E5QTD   = E5QT - E5TD
      E6TT    = E45T(7)
      E45T(8) = E5QQ
      CALL PUTREC(20,'JOBARC','NONITT3',20*IINTFP,E45T)
      WRITE(LUOUT,4000) E4T,E5ST,E4DT,E5TD,E5TT,E5QTC,E5QTD,E5QT,E6TT,
     1                  E5QQ
 4000 FORMAT(' Noniterative triple & quadruple excitation energies :',/,
     1       '           E4T                         ',F20.12,/,
     1       '           E5ST                        ',F20.12,/,
     1       '           E4DT                        ',F20.12,/,
     1       '           E5TD                        ',F20.12,/,
     1       '           E5TT                        ',F20.12,/,
     1       '           E5QT(c)                     ',F20.12,/,
     1       '           E5QT(d)                     ',F20.12,/,
     1       '           E5QT                        ',F20.12,/,
     1       '           E6TT                        ',F20.12,/,
     1       '           E5QQ                        ',F20.12)
C
C SCF
      CALL GETREC(20,'JOBARC','SCFENEG',IINTFP,EREF)
      WRITE(LUOUT,4001) EREF
 4001 FORMAT(/,'        Reference        energy ',F20.12,' a.u. ')
C
      CALL GETREC(20,'JOBARC','TOTENERG',IINTFP,ETOT)
C
      IF(METHOD.EQ.12.OR.METHOD.EQ.27)THEN
C
C     ETOT at this point is CC5SD[T] energy.
C CCSD
      ETEMP = ETOT - E4T - E5ST - E5TD - E5TT
      WRITE(LUOUT,4002) ETEMP
 4002 FORMAT('        CCSD             energy ',F20.12,' a.u. ')
C CCSD+T(CCSD)
      ETEMP = ETEMP + E4T
      WRITE(LUOUT,4003) ETEMP
 4003 FORMAT('        CCSD+T(CCSD)     energy ',F20.12,' a.u. ')
C CCSD(T)
      ETEMP = ETEMP + E5ST
      WRITE(LUOUT,4004) ETEMP
 4004 FORMAT('        CCSD(T)          energy ',F20.12,' a.u. ')
C CCSD+T*(CCSD)
      ETEMP = ETEMP - E5ST + E5TD
      WRITE(LUOUT,4005) ETEMP
 4005 FORMAT('        CCSD+T*(CCSD)    energy ',F20.12,' a.u. ')
      WRITE(LUOUT,4010) ETOT
 4010 FORMAT('        CC5SD[T]         energy ',F20.12,' a.u. ')
      ETOT = ETOT + E5QT + E5QQ
      WRITE(LUOUT,4020) ETOT
 4020 FORMAT('        CCSD+TQ*(CCSD)   energy ',F20.12,' a.u. ')
      ETOT = ETOT + E6TT
      WRITE(LUOUT,4030) ETOT
 4030 FORMAT('        CCSD(TQ)         energy ',F20.12,' a.u. ')
      CALL PUTREC(20,'JOBARC','TOTENERG',IINTFP,ETOT)
      ENDIF
C
C     QCISD(TQ)
C
      IF(METHOD.EQ.26)THEN
C
C     ETOT at this point is the QCISD(T) energy.
C
C QCISD
      ETEMP = ETOT - E4T - E5ST
      WRITE(LUOUT, 4031) ETEMP
 4031 FORMAT('        QCISD            energy ',F20.12,' a.u. ')
C QCISD+T(QCISD)
      ETEMP = ETOT - E5ST
      WRITE(LUOUT, 4032) ETEMP
 4032 FORMAT('        QCISD+T(QCISD)   energy ',F20.12,' a.u. ')
      WRITE(LUOUT,4040) ETOT
 4040 FORMAT('        QCISD(T)         energy ',F20.12,' a.u. ')
      ETOT = ETOT + E5TD + E5TT + E5QT + E5QQ
      WRITE(LUOUT,4050) ETOT
 4050 FORMAT('        QCISD+TQ*(QCISD) energy ',F20.12,' a.u. ')
      ETOT = ETOT + E5QTD + 2.0D+00*E6TT
      WRITE(LUOUT,4060) ETOT
 4060 FORMAT('        QCISD(TQ)        energy ',F20.12,' a.u. ')
      CALL PUTREC(20,'JOBARC','TOTENERG',IINTFP,ETOT)
C
      IF(IPRINT.GE.5)THEN
C
C     Debug/test aids : compute the combined quantities which are used
C     by Raghavachari et al., J. Phys. Chem. 94, 5579 (1990).
C
C     TQ(I) ~ E5QTD; TQ(II) ~ E5QTC; QQ(II) ~ E5QQ.
C     Note that E5ST already has the factor of 2.
C
      ET  = E4T + E5ST
      ETT = E5TT + 2.0D+00 * E6TT
      ETQ = 2.0D+00 * E5QT
      EQQ = E5QQ
      WRITE(LUOUT,4070) ET,ETT,ETQ,EQQ
 4070 FORMAT('        E(T)   ',F20.12,' a.u. ',/,
     1       '        E(TT)  ',F20.12,' a.u. ',/,
     1       '       2E(TQ)  ',F20.12,' a.u. ',/,
     1       '        E(QQ)  ',F20.12,' a.u. ')
      ENDIF
C
      ENDIF
C
      IF(METHOD.EQ.29)THEN
C
      ETEMP = ETOT - E5QT
      ETOT  = ETOT + E5QQ
C
      WRITE(LUOUT,4080) ETEMP,ETOT
 4080 FORMAT('        CCSDT            energy ',F20.12,' a.u. ',/,
     &       '        CCSDT+Q*(CCSDT)  energy ',F20.12,' a.u. ')
      CALL PUTREC(20,'JOBARC','TOTENERG',IINTFP,ETOT)
C
      ENDIF
C
C     Close auxiliary i/o file.
C
      CALL AUXIOO
C
      call aces_fin
      STOP
      END

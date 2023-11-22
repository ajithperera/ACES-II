      SUBROUTINE EOM_CCSDPT(ICORE,MAXCOR,IUHF,IRREPX,ROOT,ISIDEPT,
     &                      EL1R3,EL2R3,EL3R2,EL3R3,TOTL3R3,
     &                      EL2HR1T)
C
C     This subroutine drives the noniterative triples corrections to
C     EOM-CCSD excitation energies.
C
      IMPLICIT NONE
C-----------------------------------------------------------------------
C     Arguments.
C-----------------------------------------------------------------------
      INTEGER ICORE,MAXCOR,IUHF,IRREPX,ISIDEPT
      DOUBLE PRECISION ROOT,EL1R3,EL2R3,EL3R2,EL3R3,TOTL3R3,EL2HR1T,
     &                 DDOT
C-----------------------------------------------------------------------
C     Common blocks and parameters.
C-----------------------------------------------------------------------
      INTEGER IINTLN,IFLTLN,IINTFP,IALONE,IBITWD,
     &        NSTART,NIRREP,IRREPA,IRREPB,DIRPRD,
     &        IRPDPD,ISYTYP,ID,
     &        POP,VRT,NTAA,NTBB,NF1AA,NF2AA,NF1BB,NF2BB,
     &        NOCCO,NVRTO,LUOUT,MOINTS,IFLAGS,IFLAGS2,
     &        ABSVRT,ABSOCC,MAXNBF,
     &        IOFFVV,IOFFOO,IOFFVO,LIST2I1,LIST2I2,LIST2I3,LISTT2,
     &        LWIC11,LWIC12,LWIC13,LWIC14,LWIC15,LWIC16,LWIC17,LWIC18,
     &        LWIC21,LWIC22,LWIC23,LWIC24,LWIC25,LWIC26,LWIC27,LWIC28,
     &        LWIC31,LWIC32,LWIC33,LWIC34,LWIC35,LWIC36,
     &        LWIC37,LWIC38,LWIC39,LWIC40,LWIC41,LWIC42,
     &        DISTSZ,NDISTS,INIWRD,LNPHYR,NRECS,LUAUX,IDSYMSZ
C-----------------------------------------------------------------------
C     Local variables.
C-----------------------------------------------------------------------
      DOUBLE PRECISION AAAL3R3,BBBL3R3,AABL3R3,BBAL3R3,EL2HR1TEV,sum
      INTEGER ICLLVL,IREFNC,IQRHFP,IQRHFM,INONHF,NOCA,NOCB,NVRTA,NVRTB,
     &        MAXPOP,MAXVRT,MAX1,MAX2,MAX3,MAX4,MAX9,MAX10,MAX11,MAX12,
     &        MAX13,MAX14,MAXOO,MAXVV,MAXOV,SCRDIM,
     &        LNOOOV1,LNOOOV2,LNOOOV3,LNOOOV4,LNOOOV,IRREP,IRPIJ,
     &        NO,NV,ISPIN,LABCI0,LMCJK0,NFREE,NALLOC,
     &        I000,I010,I020,I030,I040,I050,I060,I070,I080,I090,I100,
     &        I110,I120,ITOP,length,kk
      LOGICAL INT1,INT2,NONHF,HF,ISKP14,ISKP2,ISKP3,INCR2,threebod,
     &        hackup
C-----------------------------------------------------------------------
      PARAMETER(MAXNBF=1000)
C-----------------------------------------------------------------------
      DIMENSION ICORE(MAXCOR)
C-----------------------------------------------------------------------
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYM/    POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF2AA,
     1                NF1BB,NF2BB
      COMMON /INFO/   NOCCO(2),NVRTO(2)
      COMMON /FILES/  LUOUT,MOINTS
      COMMON /FLAGS/  IFLAGS(100)
C
      EQUIVALENCE(IREFNC,IFLAGS(11))
      EQUIVALENCE(IQRHFP,IFLAGS(32))
      EQUIVALENCE(IQRHFM,IFLAGS(33))
      EQUIVALENCE(INONHF,IFLAGS(38))
C
      COMMON /FLAGS2/ IFLAGS2(500)
      COMMON /ACTORB/ ABSVRT(MAXNBF,8,2),ABSOCC(MAXNBF,8,2)
      COMMON /T3OFF/  IOFFVV(8,8,10),IOFFOO(8,8,10),IOFFVO(8,8,4)
      COMMON /T2ILIS/ LIST2I1,LIST2I2,LIST2I3
      COMMON /T2LIST/ LISTT2(3)
      COMMON /LISWI/  LWIC11,LWIC12,LWIC13,LWIC14,
     &                LWIC15,LWIC16,LWIC17,LWIC18,
     &                LWIC21,LWIC22,LWIC23,LWIC24,
     &                LWIC25,LWIC26,LWIC27,LWIC28,
     &                LWIC31,LWIC32,LWIC33,
     &                LWIC34,LWIC35,LWIC36,
     &                LWIC37,LWIC38,LWIC39,LWIC40,LWIC41,LWIC42
      COMMON /AUXIO / DISTSZ(8,100),NDISTS(8,100),INIWRD(8,100),LNPHYR,
     &                NRECS,LUAUX
C-----------------------------------------------------------------------
C
      Write(6,*)
      WRITE(LUOUT,1010)
 1010 FORMAT(' @EOM_CCSDPT-I, Evaluating noniterative triple terms. ')
      WRITE(LUOUT,1015) ROOT
 1015 FORMAT(' @EOM_CCSDPT-I, CCSD Root used is ',F20.12)
C
C-----------------------------------------------------------------------
C     Some traps.
C-----------------------------------------------------------------------
C Eleminate the setting of ICLLVL based on CALC flag. Lets leave old 
C numbering system intact to avoid editing many files, but set the
C value of ICLLVL based of EXCITE flag.
C
      IF (IFLAGS(87) .EQ. 11) ICLLVL = 22

      IF(ICLLVL.NE.22)THEN
        WRITE(6,1020) ICLLVL
        CALL ERREX
      ENDIF
 1020 FORMAT(' @EOM_CCSDPT, ICLLVL is ',I3,'. It must be 22. ')
C
      IF(ISIDEPT.NE.1.AND.ISIDEPT.NE.2)THEN
        WRITE(6,1030) ISIDEPT
        CALL ERREX
      ENDIF
 1030 FORMAT(' @EOM_CCSDPT, ISIDEPT is ',I3,'. It must be 1 or 2. ')
C-----------------------------------------------------------------------
C
      ISKP14 = .false.
      ISKP2  = .false.
      ISKP3  = .false.
C
      threebod = .false.
C
      NOCA  = NOCCO(1)
      NOCB  = NOCCO(2)
      NVRTA = NVRTO(1)
      NVRTB = NVRTO(2)
C
C-----------------------------------------------------------------------
C     Determine dimensions of ijk/abc triples.
C-----------------------------------------------------------------------
      CALL DIMT3(DISTSZ(1,1),NDISTS(1,1))
      CALL DIMT3(DISTSZ(1,5),NDISTS(1,5))
C
C-----------------------------------------------------------------------
C     Determine kind of reference function.
C-----------------------------------------------------------------------
      NONHF = INONHF.NE.0.OR.
     &       (IREFNC.EQ.1.AND.(IQRHFP.NE.0.OR.IQRHFM.NE.0)).OR.
     &        IREFNC.EQ.2
      HF = .NOT.NONHF
C
      WRITE(6,1061)
      IF(NONHF)THEN
      WRITE(6,1062)
      ELSE
      WRITE(6,1063)
      ENDIF
C
C-----------------------------------------------------------------------
C     Set inclusion parameters.
C-----------------------------------------------------------------------
      IF(ISIDEPT.EQ.1)THEN
       INT1 = .TRUE.
       INT2 = .TRUE.
      ELSE
       INT1 = .FALSE.
       INT2 = .FALSE.
      ENDIF
C
      IF(     INT1) WRITE(6,1072)
      IF(.NOT.INT1) WRITE(6,1073)
      IF(     INT2) WRITE(6,1074)
      IF(.NOT.INT2) WRITE(6,1075)
C
C-----------------------------------------------------------------------
C     Calculate maximum values for various two index quantities. Get
C     overall maximum for scratch array dimension.
C-----------------------------------------------------------------------
      MAXPOP = 0
      MAXVRT = 0
      MAX1   = 0
      MAX2   = 0
      MAX3   = 0
      MAX4   = 0
      MAX9   = 0
      MAX10  = 0
      MAX11  = 0
      MAX12  = 0
      MAX13  = 0
      MAX14  = 0
      DO   20 IRREP=1,NIRREP
      IF(POP(IRREP,1).GT.MAXPOP) MAXPOP = POP(IRREP,1)
      IF(POP(IRREP,2).GT.MAXPOP) MAXPOP = POP(IRREP,2)
      IF(VRT(IRREP,1).GT.MAXVRT) MAXVRT = VRT(IRREP,1)
      IF(VRT(IRREP,2).GT.MAXVRT) MAXVRT = VRT(IRREP,2)
C
      IF(IRPDPD(IRREP, 1).GT.MAX1 ) MAX1   = IRPDPD(IRREP, 1)
      IF(IRPDPD(IRREP, 2).GT.MAX2 ) MAX2   = IRPDPD(IRREP, 2)
      IF(IRPDPD(IRREP, 3).GT.MAX3 ) MAX3   = IRPDPD(IRREP, 3)
      IF(IRPDPD(IRREP, 4).GT.MAX4 ) MAX4   = IRPDPD(IRREP, 4)
      IF(IRPDPD(IRREP, 9).GT.MAX9 ) MAX9   = IRPDPD(IRREP, 9)
      IF(IRPDPD(IRREP,10).GT.MAX10) MAX10  = IRPDPD(IRREP,10)
      IF(IRPDPD(IRREP,11).GT.MAX11) MAX11  = IRPDPD(IRREP,11)
      IF(IRPDPD(IRREP,12).GT.MAX12) MAX12  = IRPDPD(IRREP,12)
      IF(IRPDPD(IRREP,13).GT.MAX13) MAX13  = IRPDPD(IRREP,13)
      IF(IRPDPD(IRREP,14).GT.MAX14) MAX14  = IRPDPD(IRREP,14)
   20 CONTINUE
      MAXOO  = MAX(MAX3,MAX4 ,MAX14)
      MAXVV  = MAX(MAX1,MAX2 ,MAX13)
      MAXOV  = MAX(MAX9,MAX10,MAX11,MAX12)
      SCRDIM = MAX(MAXOV,MAXVV,MAXOO)
C
C-----------------------------------------------------------------------
C     Evaluate symmetry packed lengths of <ij//ka> integrals/intermediates.
C-----------------------------------------------------------------------
C
      LNOOOV1 = 0
      LNOOOV2 = 0
      LNOOOV3 = 0
      LNOOOV4 = 0
      DO   30 IRPIJ=1,NIRREP
      LNOOOV1 = LNOOOV1 + IRPDPD(IRPIJ, 3) * IRPDPD(IRPIJ, 9)
      LNOOOV2 = LNOOOV2 + IRPDPD(IRPIJ, 4) * IRPDPD(IRPIJ,10)
      LNOOOV3 = LNOOOV3 + IRPDPD(IRPIJ,14) * IRPDPD(IRPIJ,11)
      LNOOOV4 = LNOOOV4 + IRPDPD(IRPIJ,14) * IRPDPD(IRPIJ,12)
   30 CONTINUE
C
      LNOOOV = MAX0(LNOOOV1,LNOOOV2,LNOOOV3,LNOOOV4)
      WRITE(6,1081) LNOOOV1,LNOOOV2,LNOOOV3,LNOOOV4,LNOOOV
C
C-----------------------------------------------------------------------
C     Set up orbital index tables
C-----------------------------------------------------------------------
C
      CALL ABSORB(ABSVRT,ABSOCC,MAXNBF)
C
C-----------------------------------------------------------------------
C     Calculate VV, OO, VO offsets
C-----------------------------------------------------------------------
C
      CALL MKOFVV
      CALL MKOFOO
      CALL MKOFVO
C
C-----------------------------------------------------------------------
C     Set W intermediate list numbers
C-----------------------------------------------------------------------
C
      CALL E_INITRP
C
C-----------------------------------------------------------------------
C     Initialize R1/L1/R2/L2 lists and increment lists.
C
C     R1 : Parts 1,2 of 410; L1 : Parts 3,4 of 410; R1 increment : 5,6.
C     R2 : 401-403; L2 : 404-406; R2/L2 Increments : 407-409.
C     L increments are currently not needed in noniterative code ---
C     INT1 and INT2 are .FALSE. when ISIDEPT=2.
C-----------------------------------------------------------------------
C
      IF(IUHF.NE.0)THEN
       CALL ZEROLIST(ICORE,MAXCOR,407)
       CALL ZEROLIST(ICORE,MAXCOR,408)
       CALL ZERO(ICORE,IRPDPD(IRREPX,10))
       CALL PUTLST(ICORE,1,1,1,6,410)
      ENDIF
      CALL ZEROLIST(ICORE,MAXCOR,409)
      CALL ZERO(ICORE,IRPDPD(IRREPX,9))
      CALL PUTLST(ICORE,1,1,1,5,410)
C
C     Increments.
C
      LIST2I1 = 407
      LIST2I2 = 408
      LIST2I3 = 409
C
C     R2/L2 vectors.
C
      IF(ISIDEPT.EQ.1)THEN
       LISTT2(1) = 401
       LISTT2(2) = 402
       LISTT2(3) = 403
      ELSE
       LISTT2(1) = 404
       LISTT2(2) = 405
       LISTT2(3) = 406
      ENDIF

C Check the bi-orthognality
C
#ifdef _DEBUG_LVL0
      If (IRREPX .NE. 1) Then
          If (Iuhf .EQ. 0) Then
              Length = IDSYMSZ(IRREPX,ISYTYP(1,46),ISYTYP(2,46)) +
     &                 IRPDPD(IRREPX,9)
              i000 = 1
              i010 = i000 + length
              i020 = i010 + length

              CALL GETLST(ICORE(i000),2,1,1,1,472)
              CALL GETLST(ICORE(i010),3,1,1,1,472)

              sum = ddot(length, icore(i000), 1, icore(i000),1)
              Write(6, "(a,F15.8)") "<R|R> =", sum 
              sum = ddot(length, icore(i010), 1, icore(i010),1)
              Write(6, "(a,F15.8)") "<L|L> =", sum
              sum = ddot(length, icore(i000), 1, icore(i010),1)
              Write(6, "(a,F15.8)") "<L|R> =", sum
          Else 
              Length = IDSYMSZ(IRREPX,ISYTYP(1,44),ISYTYP(2,44)) +
     &                 IDSYMSZ(IRREPX,ISYTYP(1,45),ISYTYP(2,45)) +
     &                 IDSYMSZ(IRREPX,ISYTYP(1,46),ISYTYP(2,46)) +
     &                 IRPDPD(IRREPX,9) + IRPDPD(IRREPX,10)
              i000 = 1
              i010 = i000 + length
              i020 = i010 + length

              CALL GETLST(ICORE(i000),2,1,1,1,472)
              CALL GETLST(ICORE(i010),3,1,1,1,472)

              sum = ddot(length, icore(i000), 1, icore(i000),1)
              Write(6, "(a,F15.8)") "<R|R> =", sum
              sum = ddot(length, icore(i010), 1, icore(i010),1)
              Write(6, "(a,F15.8)") "<L|L> =", sum
              sum = ddot(length, icore(i000), 1, icore(i010),1)
              Write(6, "(a,F15.8)") "<L|R> =", sum
          Endif  
        
      Endif 

C Prepare lists
#endif 

      CALL GETLST(ICORE,2,1,1,1,472)
      CALL MKSD(IRREPX,ICORE,401,410,0,IUHF)
      CALL GETLST(ICORE,3,1,1,1,472)
      CALL MKSD(IRREPX,ICORE,404,410,2,IUHF)

C-----------------------------------------------------------------------
C     Compute "derivatives" of W(abci), W(mcjk).
C-----------------------------------------------------------------------
C
c       if(  iflags(2).eq.33 .or.
c     &     (iflags(2).eq.22 .and. iflags2(124).eq.5) )then
c         hackup = .true.
c       else
c         hackup = .false.
c       endif
C
       LMCJK0 = 376
       LABCI0 = 380
      IF(ICLLVL.EQ.22.AND.IFLAGS2(124).GE.5)THEN
      WRITE(6,*) ' @EOM_CCSDPT-I, Doing additional CC3 mjck,abci terms'
C
       write(6,*) '@EOM_CCSDPT-I, Debug write A '
       IF(IFLAGS2(124).GE.11)THEN
         INCR2 = .TRUE.
       ELSE
         INCR2 = .FALSE.
       ENDIF
       LMCJK0 = 376
       LABCI0 = 380
       write(6,*) '@EOM_CCSDPT-I, Debug write B '
         IF(IUHF.NE.0)THEN
           CALL ZEROLIST(ICORE,MAXCOR,377)
           CALL ZEROLIST(ICORE,MAXCOR,378)
           CALL ZEROLIST(ICORE,MAXCOR,379)
           CALL ZEROLIST(ICORE,MAXCOR,381)
           CALL ZEROLIST(ICORE,MAXCOR,382)
           CALL ZEROLIST(ICORE,MAXCOR,383)
         ENDIF
         CALL ZEROLIST(ICORE,MAXCOR,380)
       write(6,*) '@EOM_CCSDPT-I, Debug write C '
         CALL ZEROLIST(ICORE,MAXCOR,384)
       write(6,*) '@EOM_CCSDPT-I, Debug write D '
C
C abci stuff.
C
       IF(IFLAGS2(124).GE.11)THEN
        CALL  GW5T1ABCD(ICORE,MAXCOR,IUHF,IRREPX,  490,0,230,LABCI0)
cfudge  CALL  GW5T1ABCD(ICORE,MAXCOR,IUHF,IRREPX,  490,0,330,LABCI0)
       write(6,*) '@EOM_CCSDPT-I, Debug write E '
       ELSE
        CALL  GW5T1ABCD(ICORE,MAXCOR,IUHF,IRREPX,  490,0,330,LABCI0)
       write(6,*) '@EOM_CCSDPT-I, Debug write F '
       ENDIF
c      CALL GCC3W5RING(ICORE,MAXCOR,IUHF,IRREPX,LABCI0,363,490,0)
       CALL GCC3W5RING(ICORE,MAXCOR,IUHF,IRREPX,LABCI0,353,490,0)
       write(6,*) '@EOM_CCSDPT-I, Debug write G '
       IF(IFLAGS2(124).GE.7)THEN
       CALL  GCC3W5TAU(ICORE,MAXCOR,IUHF,IRREPX,LABCI0,
     &                  6,443,90,0,490,0,INCR2)
       write(6,*) '@EOM_CCSDPT-I, Debug write H '
       ENDIF
C
C abci piece for nccsdt-3. [cut out for testing]
C
c       if(threebod)then
       IF(IFLAGS2(124).GE.11)THEN
coldvee
        IF(IUHF.EQ.0)THEN
       write(6,*) '@EOM_CCSDPT-I, Debug write I '
          CALL GW5INRHF(ICORE,MAXCOR,IUHF,IRREPX,LABCI0,26,400)
       write(6,*) '@EOM_CCSDPT-I, Debug write J '
        ELSE
          CALL GW5INUHF1(ICORE,MAXCOR,IUHF,IRREPX,LABCI0,26,400)
          CALL GW5INUHF2(ICORE,MAXCOR,IUHF,IRREPX,LABCI0,26,400)
          CALL GW5INUHF3(ICORE,MAXCOR,IUHF,IRREPX,LABCI0,26,400)
        ENDIF
c the +20 is for new vee
c        IF(IUHF.EQ.0)THEN
c          CALL GW5INRHF(ICORE,MAXCOR,IUHF,IRREPX,LABCI0,26,400+20)
c        ELSE
c          CALL GW5INUHF1(ICORE,MAXCOR,IUHF,IRREPX,LABCI0,26,400+20)
c          CALL GW5INUHF2(ICORE,MAXCOR,IUHF,IRREPX,LABCI0,26,400+20)
c          CALL GW5INUHF3(ICORE,MAXCOR,IUHF,IRREPX,LABCI0,26,400+20)
cc          WRITE(6,*) ' @EOM_CCSDPT-F, UHF not coded yet. '
cc          CALL ERREX
c        ENDIF
       ENDIF
c       endif
C
C mcjk stuff.
C
C ijka piece really only for nccsdt-3 [note that dhbiajk2 must be called
C before anything else. gdhbiajk2 does same job as dhbiajk2, but it can
C be called anywhere].
C
c       if(threebod)then
      IF(IFLAGS2(124).GE.11)THEN
       IF(IUHF.EQ.0)THEN
coldvee
       write(6,*) '@EOM_CCSDPT-I, Debug write K '
        CALL    DHBIAJK2(ICORE,MAXCOR,IUHF,IRREPX,1,437,439,10,LMCJK0+4)
cnewvee
c        CALL    DHBIAJK2(ICORE,MAXCOR,IUHF,IRREPX,1,457,459,10,LMCJK0+4)
calternative for old vee
c       CALL   GDHBIAJK2(ICORE,MAXCOR,IUHF,IRREPX,1,437,439,10,LMCJK0+4)
       ELSE
coldvee
       write(6,*) '@EOM_CCSDPT-I, Debug write L '
        CALL DHBIAJKA(ICORE,MAXCOR,IUHF,IRREPX,1,433,6,LMCJK0)
       write(6,*) '@EOM_CCSDPT-I, Debug write M '
        CALL DHBIAJKB(ICORE,MAXCOR,IUHF,IRREPX,1,433,6,LMCJK0)
       write(6,*) '@EOM_CCSDPT-I, Debug write N '
        CALL DHBIAJKC(ICORE,MAXCOR,IUHF,IRREPX,1,433,6,LMCJK0)
       write(6,*) '@EOM_CCSDPT-I, Debug write O '
cnewvee
c        CALL DHBIAJKA(ICORE,MAXCOR,IUHF,IRREPX,1,453,6,LMCJK0)
c        CALL DHBIAJKB(ICORE,MAXCOR,IUHF,IRREPX,1,453,6,LMCJK0)
c        CALL DHBIAJKC(ICORE,MAXCOR,IUHF,IRREPX,1,453,6,LMCJK0)
       ENDIF
      ENDIF
c       endif
C
       CALL   DHBIAJK4(ICORE,MAXCOR,IUHF,IRREPX,1,490,  350,LMCJK0) 
       CALL GCC3W4RING(ICORE,MAXCOR,IUHF,IRREPX,LMCJK0,363,490,0)
       IF(IFLAGS2(124).GE.7)THEN
coldvee ?
        CALL  GCC3W4TAU(ICORE,MAXCOR,IUHF,IRREPX,LMCJK0,
     &                  26,443,90,0,490,0,INCR2)
cnewvee
c       CALL  GCC3W4TAU(ICORE,MAXCOR,IUHF,IRREPX,LMCJK0,
c     &                 26,460,90,0,490,0,INCR2)
       ENDIF
C
C another NCCSDT-3 piece for both abci and mcjk.
C
c       if(threebod)then
       IF(IFLAGS2(124).GE.11)THEN
        CALL GWR1T2(ICORE,MAXCOR,IUHF,IRREPX,490,LABCI0,LMCJK0,
     &              .TRUE.,.TRUE.,.FALSE.,0.5D+00,0.5D+00)
c    &              .TRUE.,.TRUE.,.FALSE.,0.5D+00,0.5D+00)
       ENDIF
c       endif
C
      ENDIF
C-----------------------------------------------------------------------
C     End of computation of "derivatives" of W(abci), W(mcjk).
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
C     Transpose records of LWIC17 and LWIC27 (the latter only if it is
C     not the same list). Worry about position --- move to after W(abci)
C     and W(mcjk) derivative stuff ?
C-----------------------------------------------------------------------
C
      IF(IUHF.EQ.1)THEN
CSSS        CALL E_TSPABCI(ICORE,MAXCOR,LWIC17,1,1)
        CALL E_TSPABCI(ICORE,MAXCOR,1)
        IF(LWIC27.NE.LWIC17)THEN
CSSS          CALL E_TSPABCI(ICORE,MAXCOR,LWIC27,1,1)
          CALL E_TSPABCI(ICORE,MAXCOR,1)
        ENDIF
        IF(IFLAGS(2).EQ.22 .AND. IFLAGS2(124).GE.5)THEN
CSSS          CALL E_TSPABCI(ICORE,MAXCOR,LABCI0+3,1,IRREPX)
          CALL E_TSPABCI(ICORE,MAXCOR,1)
        ENDIF
      ENDIF
C
C-----------------------------------------------------------------------
C     Transpose ooov lists to ovoo.
C-----------------------------------------------------------------------
      Write(*,*) LWIC11, LWIC12, LWIC13, LWIC14, LWIC21, LWIC22,
     &           LWIC23, LWIC24, LMCJK0+1, LMCJK0+2, LMCJK0+3,
     &           LMCJK0+4

      IF(IUHF.NE.0)THEN
        CALL TSPLST(ICORE,MAXCOR,LWIC11,1)
        CALL TSPLST(ICORE,MAXCOR,LWIC12,1)
        CALL TSPLST(ICORE,MAXCOR,LWIC13,1)
      ENDIF
        CALL TSPLST(ICORE,MAXCOR,LWIC14,1)
C
      IF(LWIC14.NE.LWIC24)THEN
        IF(IUHF.NE.0)THEN
          CALL TSPLST(ICORE,MAXCOR,LWIC21,1)
          CALL TSPLST(ICORE,MAXCOR,LWIC22,1)
          CALL TSPLST(ICORE,MAXCOR,LWIC23,1)
        ENDIF
          CALL TSPLST(ICORE,MAXCOR,LWIC24,1)
      ENDIF
C
      IF(ISIDEPT.EQ.1.AND.IFLAGS2(124).GE.5)THEN
        IF(IUHF.NE.0)THEN
          CALL TSPLST(ICORE,MAXCOR,LMCJK0+1,IRREPX)
          CALL TSPLST(ICORE,MAXCOR,LMCJK0+2,IRREPX)
          CALL TSPLST(ICORE,MAXCOR,LMCJK0+3,IRREPX)
        ENDIF
          CALL TSPLST(ICORE,MAXCOR,LMCJK0+4,IRREPX)
      ENDIF
C-----------------------------------------------------------------------
C     End of transposition.
C-----------------------------------------------------------------------
C
      EL3R3   = 0.0D+00
      AAAL3R3 = 0.0D+00
      BBBL3R3 = 0.0D+00
      AABL3R3 = 0.0D+00
      BBAL3R3 = 0.0D+00
      TOTL3R3 = 0.0D+00
C
      IF(ISKP14.OR.IUHF.EQ.0) GOTO 1000
C
C-----------------------------------------------------------------------
C     AAA/BBB spin case.
C-----------------------------------------------------------------------
C
      DO  110 ISPIN=1,2
      NO = NOCCO(ISPIN)
      NV = NVRTO(ISPIN)
C
      IF(NO.LT.3.OR.NV.LT.3) GOTO 110
C
C     I000   OOOV INTEGRALS
C     I010   OOOV intermediates
C     I020   BUF  BUFFER FOR OOOV INTEGRAL LIST READ
C     I030   IADBLK
C     I040   LENBLK
C     I050   IADT2
C     I060   LENT2
C     I070   IADV
C     I080   LENV
C     I090   OOOVAD
C     I100   LENINT
C     I110   EVAL
C     I120   ICORE (THE REST)
C
      I000 = 1
      I010 = I000 + IINTFP * LNOOOV
      I020 = I010 + IINTFP * LNOOOV
      I030 = I020 + IINTFP * MAXOO
      I040 = I030 + IINTFP * 8
      I050 = I040 + IINTFP * 8
      I060 = I050 + IINTFP * 8
      I070 = I060 + IINTFP * 8
      I080 = I070 + IINTFP * 8
      I090 = I080 + IINTFP * 8
      I100 = I090 + IINTFP * 8
      I110 = I100 + IINTFP * 8
      I120 = I110 + IINTFP * (NO + NV)
      ITOP = I120
C     COMPUTE NUMBER OF FREE INTEGER WORDS
      NFREE = MAXCOR + 1 - ITOP
      NALLOC = ITOP - 1
      WRITE(6,2001) MAXCOR,NALLOC,NFREE
      IF(NFREE.LE.0)THEN
      WRITE(6,2002)
      CALL INSMEM('TRPS',ITOP,MAXCOR)
      ENDIF
CSSS      Write(6,*) "entering NT3AAA"
      CALL NT3AAA(ICORE(I000),ICORE(I010),ICORE(I020),
     1            ICORE(I030),ICORE(I040),ICORE(I050),ICORE(I060),
     1            ICORE(I070),ICORE(I080),ICORE(I090),ICORE(I100),
     1            ICORE(I110),
     1            LNOOOV,NO,NV,
     1            ICORE(ITOP),
     1            NFREE,ISPIN,INT1,INT2,NONHF,IRREPX,ISIDEPT,ROOT,
     1            EL3R3,AAAL3R3,BBBL3R3)
  110 CONTINUE
 1000 CONTINUE
C
      IF(NOCA .LT.2.OR.NOCB .LT.1.OR.
     1   NVRTA.LT.2.OR.NVRTB.LT.1.OR.ISKP2) GOTO 2000
C
C-----------------------------------------------------------------------
C     AAB spin case.
C-----------------------------------------------------------------------
C
C
C     I000   D1T1 alpha
C     I010   D1T1 beta
C     I020   S1   alpha
C     I030   S2   beta
C     I040   SCR1
C     I050   SCR2
C     I060   SCR3
C     I070   EVAL
C     I080   ICORE (THE REST)
C
      I000 =    1
      I010 = I000 + IINTFP * MAXOV
      I020 = I010 + IINTFP * MAXOV
      I030 = I020 + IINTFP * MAXOV
      I040 = I030 + IINTFP * MAXOV
      I050 = I040 + IINTFP * SCRDIM
      I060 = I050 + IINTFP * SCRDIM
      I070 = I060 + IINTFP * SCRDIM
      I080 = I070 + IINTFP * 2 * (NOCA + NVRTA)
      ITOP = I080
C
C     Compute number of free integer words.
C
      NFREE  = MAXCOR + 1 - ITOP
      NALLOC = ITOP   - 1
      WRITE(6,2011) MAXCOR,NALLOC,NFREE
      IF(NFREE.LE.0)THEN
      WRITE(6,2012)
      CALL INSMEM('TRPS',ITOP,MAXCOR)
      ENDIF
CSSS      Write(6,*) "Entering Nt3AAB"
      CALL NT3AAB(ICORE(I000),ICORE(I010),ICORE(I020),ICORE(I030),
     &            ICORE(I040),ICORE(I050),ICORE(I060),
     &            ICORE(I070),
     &            ICORE(ITOP),
     &            NFREE,INT1,INT2,NONHF,IUHF,IRREPX,ISIDEPT,ROOT,
     &            NOCA,NOCB,NVRTA,NVRTB,EL3R3,AABL3R3)
 2000 CONTINUE
C
      IF(IUHF.EQ.0.OR.
     1   NOCA .LT.1.OR.NOCB .LT.2.OR.
     1   NVRTA.LT.1.OR.NVRTB.LT.2.OR.ISKP3) GOTO 3000
C
C-----------------------------------------------------------------------
C     BBA spin case.
C-----------------------------------------------------------------------
C
C     Reorganize abci/ijka intermediate lists, t2 ab lists and,
C     if appropriate, d2t2 ab increment lists.
C
      IF(IUHF.EQ.1)THEN
CSSS        CALL E_TSPABCI2(ICORE,MAXCOR,LWIC17,LWIC18,1,1)
        CALL E_TSPABCI2(ICORE,MAXCOR,1)
        IF(LWIC17.NE.LWIC27)THEN
CSSS          CALL E_TSPABCI2(ICORE,MAXCOR,LWIC27,LWIC28,1,1)
          CALL E_TSPABCI2(ICORE,MAXCOR,1)
        ENDIF
        IF(IFLAGS(2).EQ.22 .AND. IFLAGS2(124).GE.5)THEN
          CALL E_TSPABCI2(ICORE,MAXCOR,LABCI0+3,LABCI0+4,1,IRREPX)
        ENDIF
        IF(ISIDEPT.EQ.2) CALL AT2IAB(ICORE,      403,IRREPX,0)
                       CALL AT2IAB(ICORE,LISTT2(3),IRREPX,0)
                       CALL AT2IAB(ICORE,       46,     1,0)
        IF(INT2)       CALL AT2IAB(ICORE,LIST2I3  ,IRREPX,0)
        IF(INT1.OR.ISIDEPT.EQ.2) CALL AT2IAB(ICORE,       16,     1,0)
      ENDIF
C
      I000 =    1
      I010 = I000 + IINTFP * MAXOV
      I020 = I010 + IINTFP * MAXOV
      I030 = I020 + IINTFP * MAXOV
      I040 = I030 + IINTFP * MAXOV
      I050 = I040 + IINTFP * SCRDIM
      I060 = I050 + IINTFP * SCRDIM
      I070 = I060 + IINTFP * SCRDIM
      I080 = I070 + IINTFP * 2 * (NOCA + NVRTA)
      ITOP = I080
C
C     Compute number of free integer words.
C
      NFREE = MAXCOR + 1 - ITOP
      NALLOC = ITOP - 1
      WRITE(6,2021) MAXCOR,NALLOC,NFREE
      IF(NFREE.LE.0)THEN
      WRITE(6,2022)
      CALL INSMEM('EOM_CCSDPT',ITOP,MAXCOR)
      ENDIF
      Write(6,*) "Entering NT3BBA"
      CALL NT3BBA(ICORE(I000),ICORE(I010),ICORE(I020),ICORE(I030),
     1            ICORE(I040),ICORE(I050),ICORE(I060),
     1            ICORE(I070),ICORE(ITOP),NFREE,
     1            INT1,INT2,NONHF,IUHF,IRREPX,ISIDEPT,ROOT,
     1            NOCA,NOCB,NVRTA,NVRTB,EL3R3,BBAL3R3)
C
      IF(IUHF.EQ.1)THEN
c      CALL TSPABCI2(ICORE,MAXCOR,2)
CSSS        CALL E_TSPABCI2(ICORE,MAXCOR,LWIC17,LWIC18,2,1)
        CALL E_TSPABCI2(ICORE,MAXCOR,2)
        IF(LWIC17.NE.LWIC27)THEN
CSSS          CALL E_TSPABCI2(ICORE,MAXCOR,LWIC27,LWIC28,2,1)
          CALL E_TSPABCI2(ICORE,MAXCOR,2)
        ENDIF
        IF(IFLAGS(2).EQ.22 .AND. IFLAGS2(124).GE.5)THEN
          CALL E_TSPABCI2(ICORE,MAXCOR,LABCI0+3,LABCI0+4,2,IRREPX)
        ENDIF
        IF(ISIDEPT.EQ.2) CALL AT2IAB(ICORE,      403,IRREPX,1)
                       CALL AT2IAB(ICORE,LISTT2(3),IRREPX,1)
                       CALL AT2IAB(ICORE,       46,     1,1)
        IF(INT2)       CALL AT2IAB(ICORE,LIST2I3  ,IRREPX,1)
        IF(INT1.OR.ISIDEPT.EQ.2) CALL AT2IAB(ICORE,       16,     1,1)
      ENDIF
C
 3000 CONTINUE
C
C-----------------------------------------------------------------------
C     Transpose ovoo lists back to ooov.
C-----------------------------------------------------------------------
C
      IF(IUHF.NE.0)THEN
        CALL TSPLST(ICORE,MAXCOR,LWIC11,1)
        CALL TSPLST(ICORE,MAXCOR,LWIC12,1)
        CALL TSPLST(ICORE,MAXCOR,LWIC13,1)
      ENDIF
        CALL TSPLST(ICORE,MAXCOR,LWIC14,1)
C
      IF(LWIC14.NE.LWIC24)THEN
        IF(IUHF.NE.0)THEN
          CALL TSPLST(ICORE,MAXCOR,LWIC21,1)
          CALL TSPLST(ICORE,MAXCOR,LWIC22,1)
          CALL TSPLST(ICORE,MAXCOR,LWIC23,1)
        ENDIF
          CALL TSPLST(ICORE,MAXCOR,LWIC24,1)
      ENDIF
C
      IF(ISIDEPT.EQ.1.AND.IFLAGS2(124).GE.5)THEN
        IF(IUHF.NE.0)THEN
          CALL TSPLST(ICORE,MAXCOR,LMCJK0+1,IRREPX)
          CALL TSPLST(ICORE,MAXCOR,LMCJK0+2,IRREPX)
          CALL TSPLST(ICORE,MAXCOR,LMCJK0+3,IRREPX)
        ENDIF
          CALL TSPLST(ICORE,MAXCOR,LMCJK0+4,IRREPX)
      ENDIF
C
C-----------------------------------------------------------------------
C     Put 29 type of list back to AbIc order.
C-----------------------------------------------------------------------
C
c      IF(IUHF.EQ.1)THEN
c        CALL TSPABCI(ICORE,MAXCOR,2)
c      ENDIF
      IF(IUHF.EQ.1)THEN
CSSS        CALL E_TSPABCI(ICORE,MAXCOR,LWIC17,2,1)
        CALL E_TSPABCI(ICORE,MAXCOR,2)
        IF(LWIC27.NE.LWIC17)THEN
CSSS          CALL E_TSPABCI(ICORE,MAXCOR,LWIC27,2,1)
          CALL E_TSPABCI(ICORE,MAXCOR,2)

        ENDIF
        IF(IFLAGS(2).EQ.22 .AND. IFLAGS2(124).GE.5)THEN
          CALL E_TSPABCI(ICORE,MAXCOR,LABCI0+3,2,IRREPX)
        ENDIF
      ENDIF
C
C-----------------------------------------------------------------------
C     Evaluate EL1R3, RL2R3.
C-----------------------------------------------------------------------
C
      IF(ISIDEPT.EQ.1)THEN
        CALL CLCNT3(ICORE,MAXCOR,IUHF,ISIDEPT,IRREPX,EL1R3,EL2R3,EL3R2)
      ENDIF
C-----------------------------------------------------------------------
C     Set total overlap term.
C-----------------------------------------------------------------------
C
      IF(IUHF.EQ.0)THEN
        TOTL3R3 = AABL3R3
      ELSE
        TOTL3R3 = AAAL3R3 + BBBL3R3 + AABL3R3 + BBAL3R3
      ENDIF
C
C-----------------------------------------------------------------------
C     Now compute the L2*A_DS(T)*R1 energy. First zero doubles increment
C     lists.
C-----------------------------------------------------------------------
C
      IF(IFLAGS2(124).EQ.4.OR.IFLAGS2(124).EQ.6.OR.IFLAGS2(124).EQ.8.OR.
     &   IFLAGS2(124).EQ.10.OR.IFLAGS2(124).EQ.12)THEN
       IF(ISIDEPT.EQ.1)THEN
        IF(IUHF.NE.0)THEN
         CALL ZEROLIST(ICORE,MAXCOR,407)
         CALL ZEROLIST(ICORE,MAXCOR,408)
        ENDIF
        CALL ZEROLIST(ICORE,MAXCOR,409)
        CALL EL2R1(ICORE,MAXCOR,IUHF,1,IRREPX,
     &             388,384,410,0,403,406,EL2HR1T)
        write(6,*) ' eomnt3-i, el2hr1t is ',el2hr1t
        el2hr1tev = 27.2116d+00 * el2hr1t
        write(6,*) ' eomnt3-i, el2hr1t is ',el2hr1tev
       ENDIF
      ENDIF
C-----------------------------------------------------------------------
C
      RETURN
 1061 FORMAT(' @TRPS-I, Reference function information : ')
 1062 FORMAT('          Reference function is general non-HF ')
 1063 FORMAT('          Reference function is Hartree-Fock ')
 1064 FORMAT('          Reference function is RHF Hartree-Fock ')
 1065 FORMAT('          Reference function is RHF non-HF ')
 1066 FORMAT('          Reference function is UHF Hartree-Fock ')
 1067 FORMAT(' @TRPS-I, Reference function is UHF non-HF ')
 1068 FORMAT(' @TRPS-I, Reference function is ROHF ')
 1069 FORMAT(' @TRPS-I, Reference function is QRHF ')
 1072 FORMAT(' @TRPS-I, INT1 is true ')
 1073 FORMAT(' @TRPS-I, INT1 is false ')
 1074 FORMAT(' @TRPS-I, INT2 is true ')
 1075 FORMAT(' @TRPS-I, INT2 is false ')
 1081 FORMAT(' @TRPS-I, Lengths of ijka integrals. #1 ',I10,/,
     1       '                                     #2 ',I10,/,
     1       '                                     #3 ',I10,/,
     1       '                                     #4 ',I10,/,
     1       '                                    Max ',I10)
 2001 FORMAT(' @TRPS-I, Available memory in integer words ',I12,/,
     1 '          Number of pre-allocated words for TRPS14 ',I12,/,
     1 '          Number of free words for TRPS14 ',I12)
 2002 FORMAT(' TRPS-I, Insufficient memory for TRPS14 pre-allocation. ')
 2011 FORMAT(' @TRPS-I, Available memory in integer words ',I12,/,
     1 '          Number of pre-allocated words for TRPS2 ',I12,/,
     1 '          Number of free words for TRPS2 ',I12)
 2012 FORMAT(' TRPS-I, Insufficient memory for TRPS2  pre-allocation. ')
 2021 FORMAT(' @TRPS-I, Available memory in integer words ',I12,/,
     1 '          Number of pre-allocated words for TRPS3 ',I12,/,
     1 '          Number of free words for TRPS3 ',I12)
 2022 FORMAT(' TRPS-I, Insufficient memory for TRPS3  pre-allocation. ')
      END

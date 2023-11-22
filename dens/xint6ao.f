      SUBROUTINE XINT6AO(DIVO,DIOO,ICORE,MAXCOR,IRREPX,IUHF)
C
C THIS ROUTINE CALCULATES THE SIXTH TERM OF THE OCCUPIED-VIRTUAL BLOCK
C OF THE INTERMEDIATE XIA
C
C THE GENERAL FORMULA IS
C
C 1. TERM : + 1/2 SUM E,F,G G(EF,GI) <EF||GA>
C
C          SPIN TYPES : AA   AAAA AAAA    (UHF)
C                            ABAB ABAB    (UHF + RHF)
C                       BB   BBBB BBBB    (UHF)
C                            BABA BABA    (UHF)
C
C  FOR RHF A SPIN ADAPTED CODE IS USED SO ONLY ONE TERM
C  HAS TO BE CALCULATED.
C
C  THIS TERM IS ONLY REQUIRED FOR METHODS WHICH INCLUDE
C  SINGLE EXCITATION, E.G. MBPT(4), QCISD, CCSD
C
C The ACES II Gainsville version was put together by A. Perera 
C liberaly using routines written by Jurgen Gauss. I must say
C that some of the routines used in here are so complicated that
C only Jurgen can write them!! 
C 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DIRPRD,POP,VRT,IXIA(2)
      CHARACTER*80 FNAME
C
      DIMENSION ICORE(MAXCOR),DIVO(*),DIOO(*)
      LOGICAL DONE
C
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMINF/NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON/AOSYM/IAOPOP(8),IOFFAO(8),IOFFV(8,2),IOFFO(8,2),
     &             IRPDPDAO(8),IRPDPDAOS(8),ISTART(8,8),ISTARTMO(8,3)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NF1(2),NF2(2)
      COMMON/INFO/NOCCO(2),NVRTO(2)
      COMMON/G6OOC/IPOPS(8,8),IPOPE(8,8),NDONE(8,8),NSTARTM,
     &             ISTARTI,NSIZG,IOFFSET,DONE
C
      DATA ZILCH,ONE,ONEM/0.D0,1.D0,-1.D0/
C
      IONE=1
C
C Unit for integral files
C
      LUINT=10
C
C Get the number of MOS
C
      NMO = NOCCO(1)+NVRTO(1)
C
C Initialize /AOSYM/
C
      CALL GETAOINF(IUHF, IRREPX)
C
C Get the number of AOs in the calculation. 
C
      NAO = 0
      DO IRREP = 1, NIRREP
         NAO = NAO + IAOPOP(IRREP)
      ENDDO
C
C List number for the G(ab,ci) 
C
      LISTMO = 126
C
      IXIA(1) = 1
      IXIA(2) = IXIA(1) + NOCCO(1)*NAO*IUHF*IINTFP
      ISTARTC = IXIA(2) + NOCCO(2)*NAO*IINTFP
      MXCOR   = MAXCOR - IINTFP*(NOCCO(1)+IUHF*NOCCO(2))*NAO
      CALL ZERO(ICORE(IXIA(1)),(NOCCO(1)+IUHF*NOCCO(2))*NAO)

      DO ISPIN = 4, 4-3*IUHF,-1
C
C Information for multiple passes 
C
       DONE = .FALSE.
       CALL IZERO(NDONE,64)
C
C NSTARTM keeps the information about the starting MOs.
C
       ISTARTI = 1 
       NSTARTM = 1
       IOFFSET = 0
C
C Spin of the occupied orbitals.
C
       IF (ISPIN .LE. 2) THEN
          ISPINI = ISPIN
       ELSE
          ISPINI= ISPIN - 2
       ENDIF
C
C Reentry point for multiple passes of G(ab,ci)
C
 1     CONTINUE
C
       DO  IRREP1 = 1, NIRREP
           DO  IRREP2 = 1, NIRREP
               IPOPS(IRREP1,IRREP2) = 1
           ENDDO
       ENDDO
C
       CALL IZERO(IPOPE, 64)
C
C Transform G(ab,ci) to AO-basis, IRREPX is symmetry of G(ab,ci)    
C
       CALL DG6TOAO(ICORE(ISTARTC), MXCOR, IUHF, ISPIN, LISTMO,
     &              NAO, IRREPX)
C
       ILNBUF = 600
       I140   = ISTARTC
C
C Memory for G(MU NU, SIGMA I)
C
       I000 = I140 + IINTFP*NSIZG
C
C Memory for BUF
C
       I010 = I000 + ILNBUF*IINTFP
C
C Memory for IBUF
C
       I120 = I010 + ILNBUF*2
C
C Memory for IAOSYM
C
       I130 = I120 + NAO
C
C Memory for symmetry vector IMAP
C
       IEND = I130 + NAO*NAO
C
       IOFF = I120
C
       DO IRREP = 1, NIRREP
          DO I = 1, IAOPOP(IRREP)
               ICORE(IOFF) = IRREP
               IOFF        = IOFF + 1
          ENDDO
       ENDDO
C
       CALL AOSYMVEC(ICORE(I130), ICORE(I120), NAO)
C
       CALL GFNAME('IIII    ',FNAME,ILENGTH)
       OPEN(UNIT=LUINT,FILE=FNAME(1:ILENGTH),
     &      FORM='UNFORMATTED',STATUS='OLD')
       CALL RDAOIJKL3(ICORE(I140),ICORE(IXIA(ISPINI)),
     &                ICORE(I000),ICORE(I010),
     &                ICORE(I120),ICORE(I130),POP(1,ISPINI),
     &                ILNBUF,LUINT,IUHF,NAO,IRREPX)
       CLOSE(UNIT=LUINT,STATUS='KEEP')
C
       IF (NIRREP .GT. 1) THEN
          CALL GFNAME('IJIJ    ',FNAME,ILENGTH)
          OPEN(UNIT=LUINT,FILE=FNAME(1:ILENGTH),
     &         FORM='UNFORMATTED',STATUS='OLD')
          CALL RDAOIJKL3(ICORE(I140),ICORE(IXIA(ISPINI)),
     &                   ICORE(I000),ICORE(I010),
     &                   ICORE(I120),ICORE(I130),POP(1,ISPINI),
     &                   ILNBUF,LUINT,IUHF,NAO,IRREPX)
          CLOSE(UNIT=LUINT,STATUS='KEEP')
C
          CALL GFNAME('IIJJ    ',FNAME,ILENGTH)
          OPEN(UNIT=LUINT,FILE=FNAME(1:ILENGTH),
     &         FORM='UNFORMATTED',STATUS='OLD')
          CALL RDAOIJKL3(ICORE(I140),ICORE(IXIA(ISPINI)),
     &                   ICORE(I000),ICORE(I010),
     &                   ICORE(I120),ICORE(I130),POP(1,ISPINI),
     &                   ILNBUF,LUINT,IUHF,NAO,IRREPX)
          CLOSE(UNIT=LUINT,STATUS='KEEP')
C
         IF (NIRREP .GT. 2) THEN
            CALL GFNAME('IJKL    ',FNAME,ILENGTH)
            OPEN(UNIT=LUINT,FILE=FNAME(1:ILENGTH),
     &           FORM='UNFORMATTED',STATUS='OLD')
            CALL RDAOIJKL3(ICORE(I140),ICORE(IXIA(ISPINI)),
     &                     ICORE(I000),ICORE(I010),
     &                     ICORE(I120),ICORE(I130),POP(1,ISPINI),
     &                     ILNBUF,LUINT,IUHF,NAO,IRREPX)
            CLOSE(UNIT=LUINT,STATUS='KEEP')

         ENDIF
      ENDIF
C
      IOFFSET = IOFFSET + NSIZG
      IF(.NOT. DONE) GO TO 1
C
C
C This endo is for the loop over spins types. 
C
      ENDDO
C
C Transform DIVO to MO representation 
C
      I010 = ISTARTC
      IEND = I010 + IINTFP*NAO*NMO
C
       FACT=ONEM
       IRREPXAO = DIRPRD(IRREPX, IRREPX)
       CALL DTRANDIVO(ICORE(IXIA(1)),ICORE(I010),ICORE(IEND),
     &                DIVO,DIOO,FACT,NAO,NMO,IUHF,IRREPXAO)
C
      RETURN
      END

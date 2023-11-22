      SUBROUTINE PCCD_GFOCK_DEBUG(ICORE,MAXCOR,IUHF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      INTEGER POP,VRT,DIRPRD,POP2,VRT2
      LOGICAL QRHF,NONHF,ROHF,RELAXED
      LOGICAL CCD,CANON,NONHF_TERMS_EXIST,TRULY_NONHF
      LOGICAL MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,QCISD,CCSD,UCC
C
      DIMENSION ICORE(MAXCOR)

      COMMON/ADD/SUM
      COMMON/INFO/NOCCO(2), NVRTO(2)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMLOC/ISYMOFF(8,8,25)
#include "flags.h" /* for h_IFLAGS_* and h_IFLAGS2_* */
#include "flags.com"
#include "flags2.com"

      COMMON/QRHFINF/POPRHF(8),VRTRHF(8),NOSH1(8),NOSH2(8),
c&two line mod
     &               POPDOC(8),VRTDOC(8),NAI,N1I,NA2,
     &               NUMISCF,NUMASCF,ISPINP,ISPINM,IQRHF
CJDW KKB stuff
      COMMON /SHIFT/   ISHIFT
      COMMON /SYMDROP/ NDRPOP(8),NDRVRT(8)
      COMMON /DROPGEO/ NDRGEO
      COMMON /SYM/ POP(8,2),VRT(8,2),MTAA,MTBB,MD1AA,MD1BB,
     &             MD2AA,MD2BB
      COMMON /SYM2/ POP2(8,2),VRT2(8,2),NTAA,NTBB,ND1AA,ND1BB,
     &              ND2AA,ND2BB
      COMMON /SYMPOP/  IRP_DM(8,22),ISYTYP(2,500),ID(18)
      COMMON /SYMPOP2/ IRPDPD(8,22)
      COMMON/METH/MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC
C
      DATA HALF,ONE,IONE /0.5D0,1.0D0,1/
C
      SUM = 0.0D0
      I0  = IONE
      NDRGEO = 0
      ISHIFT = 0

      NONHF = IFLAGS(h_IFLAGS_nonhf).EQ.1
      QRHF  = IFLAGS(h_IFLAGS_qrhf_gen).NE.0
      ROHF  = IFLAGS(h_IFLAGS_ref).EQ.2
      CCD   = IFLAGS(h_IFLAGS_calc) .EQ. 5 .OR.
     &        IFLAGS(h_IFLAGS_calc) .EQ. 8 .OR.
     &        IFLAGS(h_IFLAGS_calc) .EQ. 52 .OR.
     &        IFLAGS(h_IFLAGS_calc) .EQ. 53 .OR.
     &        IFLAGS(h_IFLAGS_calc) .EQ. 54
      CANON = IFLAGS(h_IFLAGS_pert_orb).EQ.1
      TRULY_NONHF = (IFLAGS(h_IFLAGS_ref).GE.2 .OR.
     &               IFLAGS(h_IFLAGS_qrhf_gen) .NE. 0)
      RELAXED    = .TRUE.
      MBPT2      = .FALSE.
      MBPT3      = .FALSE.
      M4DQ       = .FALSE.
      M4DSQ      = .FALSE.
      M4SDTQ     = .FALSE.
      CCD        = .TRUE.
      QCID       = .FALSE.
      CCSD       = .FALSE.
      UCC        = .FALSE.
      Bredundant = .TRUE.
      ISHIFT     = 0
C
      CALL DCOPY(8*22,IRP_DM,1,IRPDPD,1)
      CALL DCOPY(16,POP,1,POP2,1)
      CALL DCOPY(16,VRT,1,VRT2,1)
      DO I = 1, 8
         NDRPOP(I)= 0
         NDRVRT(I)= 0
      ENDDO 
      NTAA  = MTAA
      NTBB  = MTBB
      ND1AA = MD1AA
      ND1BB = MD1BB
      ND2AA = MD2AA
      ND2BB = MD2BB

      MXCOR=MAXCOR
C
C MEMORY FOR THE OCCUPIED-OCCUPIED BLOCK
C
      IDOO=I0+MXCOR-(ND1AA+IUHF*ND1BB)*IINTFP
      MXCOR=MXCOR-(ND1AA+IUHF*ND1BB)*IINTFP
C
C MEMORY FOR THE VIRTUAL-VIRTUAL BLOCK
C
      IDVV=IDOO-(ND2AA+IUHF*ND2BB)*IINTFP
      MXCOR=MXCOR-(ND2AA+IUHF*ND2BB)*IINTFP
C
C MEMORY FOR THE VIRTUAL-OCCUPIED BLOCK
C
      IDOV=IDVV-(NTAA+IUHF*NTBB)*IINTFP
      MXCOR=MXCOR-(NTAA+IUHF*NTBB)*IINTFP
C
C MEMORY FOR THE XOV INTERMEDIATE
C
C FOR NON HF WE HAVE TO ALLOCATE ADDITIONAL MEMORY FOR THE XOV INTERMEDIATE
C FOR HF CASES EQUIVALENCE MEMORY FOR XIA WITH THAT FOR DOV
C
      IF(QRHF.OR.NONHF.OR.ROHF) THEN
       IXOV=IDOV-(MAX(NTAA,NTBB)+IUHF*MAX(NTAA,NTBB))*IINTFP
       MXCOR=MXCOR-(MAX(NTAA,NTBB)+IUHF*MAX(NTAA,NTBB))*IINTFP
       LXOV=(1+IUHF)*MAX(NTAA,NTBB)
      ELSE
       IXOV=IDOV
       LXOV=0
      ENDIF
C
C MEMORY FOR THE IOV INTERMEDIATE
C
      IIOV=IXOV-(NTAA+IUHF*NTBB)*IINTFP
      MXCOR=MXCOR-(NTAA+IUHF*NTBB)*IINTFP
C
C MEMORY FOR THE IVO INTERMEDIATE
C
      IIVO=IIOV-(NTAA+IUHF*NTBB)*IINTFP
      MXCOR=MXCOR-(NTAA+IUHF*NTBB)*IINTFP
C
C MEMORY FOR THE IOO INTERMEDIATE
C
      IIOO=IIVO-(ND1AA+IUHF*ND1BB)*IINTFP
      MXCOR=MXCOR-(ND1AA+IUHF*ND1BB)*IINTFP
C
C MEMORY FOR THE IVV INTERMEDIATE
C
      IIVV=IIOO-(ND2AA+IUHF*ND2BB)*IINTFP
      MXCOR=MXCOR-(ND2AA+IUHF*ND2BB)*IINTFP
C
C SET ITOP SO THAT NON-HF WILL WORK
C
      ITOP=IIVV
C
C ZERO THE ALLOCATED SPACE FOR THE INTERMEDIATES IOO AND IVV
C
      LENGTH=LXOV+3*NTAA+2*ND1AA+2*ND2AA+
     &                IUHF*(3*NTBB+2*ND1BB+2*ND2BB)

      IF (.NOT.TRULY_NONHF .AND. NONHF) THEN
         CALL ANALYZE_FOCK(ICORE(IDOO), ICORE(IDVV), ICORE(IDOV),
     &                     ND1AA, ND2AA, NTAA, NONHF_TERMS_EXIST)
      ELSE
         NONHF_TERMS_EXIST = TRULY_NONHF
      END IF

      CALL DZERO(ICORE(IIOO),(ND1AA+IUHF*ND1BB))
      CALL DZERO(ICORE(IIVV),(ND2AA+IUHF*ND2BB))
      CALL DZERO(ICORE(IIVO),(NTAA+IUHF*NTBB))
      CALL DZERO(ICORE(IIOV),(NTAA+IUHF*NTBB))
      CALL DZERO(ICORE(IDOO),(ND1AA+IUHF*ND1BB))
      CALL DZERO(ICORE(IDVV),(ND2AA+IUHF*ND2BB))
C
      CALL PCCD_DENSOO(ICORE(IDOO),ICORE(i0),MXCOR,IUHF)
      CALL PCCD_DENSVV(ICORE(IDVV),ICORE(i0),MXCOR,IUHF)
      CALL PCCD_DENSVO(ICORE(IDOV),ICORE(i0),MXCOR,IUHF,.TRUE.)
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,*) "The checksums of the Density blocks"
      Call checksum("DENSOO ", ICORE(IDOO),ND1AA+IUHF*ND1BB)
      Call checksum("DENSVV ", ICORE(IDVV),ND2AA+IUHF*ND2BB)
      Call checksum("DENSOV ", ICORE(IDOV),NTAA+IUHF*NTBB)
      Write(6,*)
      Write(6,"(a)") "Density matrices"
       CALL PCCD_DUMPIT(ICORE(IDOO),ICORE(IDVV),ICORE(IDOV),
     &                  ICORE(IIOV),IUHF,"AL")
      Write(6,*)
#endif
      IF(IFLAGS(1).GE.15)THEN
       write(6,*) ' I intermediates after  INTO0 '
       CALL PCCD_DUMPIT(ICORE(IIOO),ICORE(IIVV),ICORE(IIVO),
     &                  ICORE(IIOV),IUHF,"OO")
      ENDIF
C
C NOTE THE INTERMEDIATES ARE CALCULATED HERE IN ORDER TO
C BE ABLE TO RESORT SOME OF THE GAMMA INTERMEDIATES TO
C A NEW ORDERING AS REQUIRED FOR THE CALCULATION OF IOV
C
C JCP-95-2623-91, Eqn. 43, 3-6 terms when (p,q = i,j; INTO1;(r,s,t = m,e,f
C INTO2; (r,s,t = m,n,o), INTO3; (r,s,t = m,e,f), INTO4 (r,s,t = e,f,g)
C INTO5; (r,s,t = m,n,e), INTO6; (r,s,t = m,n,e). Note that when r,s, t
C is m,e,f or m,n,e the i and j can permute to give two distinct terms.
C That is the reason behind INTO1,INTO3 and INTO5, INTO6 involve same
C r,s,t indices. 09/2004, Ajith Perera
C
      CALL INTO1(ICORE(IIOO),ICORE(i0),MXCOR,IUHF)

      IF(IFLAGS(1).GE.15)THEN
       write(6,*) ' I intermediates after  INTO1 '
       CALL PCCD_DUMPIT(ICORE(IIOO),ICORE(IIVV),ICORE(IIVO),
     &                  ICORE(IIOV),IUHF,"OO")
      ENDIF

      CALL INTO2(ICORE(IIOO),ICORE(i0),MXCOR,IUHF)

      IF(IFLAGS(1).GE.15)THEN
       write(6,*) ' I intermediates after  INTO2 '
       CALL PCCD_DUMPIT(ICORE(IIOO),ICORE(IIVV),ICORE(IIVO),
     &                  ICORE(IIOV),IUHF,"OO")
      ENDIf

      CALL INTO3(ICORE(IIOO),ICORE(i0),MXCOR,IUHF,.TRUE.)

      IF(IFLAGS(1).GE.15)THEN
       write(6,*) ' I intermediates after  INTO3 '
       CALL PCCD_DUMPIT(ICORE(IIOO),ICORE(IIVV),ICORE(IIVO),
     &                  ICORE(IIOV),IUHF,"OO")
      ENDIF

      IF(.NOT.CCD) THEN
C
C IF, THE <AB||CD> CONTRIBUTION TO XIA IS CALCULATED
C WITH AO INTEGRALS, SKIP INTO4. THIS CONTRIBUTION
C IS THEN A BY-PRODUCT IN XINT6AO
C
      CALL INTO4(ICORE(IIOO),ICORE(i0),MXCOR,IUHF)

      IF(IFLAGS(1).GE.15)THEN
       write(6,*) ' I intermediates after  INTO4-MO '
       CALL PCCD_DUMPIT(ICORE(IIOO),ICORE(IIVV),ICORE(IIVO),
     &                  ICORE(IIOV),IUHF,"OO")
      ENDIF

      CALL INTO5(ICORE(IIOO),ICORE(i0),MXCOR,IUHF)

      IF(IFLAGS(1).GE.15)THEN
       write(6,*) ' I intermediates after  INTO5 '
       CALL PCCD_DUMPIT(ICORE(IIOO),ICORE(IIVV),ICORE(IIVO),
     &                  ICORE(IIOV),IUHF,"OO")
      ENDIF

      CALL INTO6(ICORE(IIOO),ICORE(i0),MXCOR,IUHF)

      IF(IFLAGS(1).GE.15)THEN
       write(6,*) ' I intermediates after  INTO6 '
       CALL PCCD_DUMPIT(ICORE(IIOO),ICORE(IIVV),ICORE(IIVO),
     &                  ICORE(IIOV),IUHF,"OO")
      ENDIF

      ENDIF

      IF(.NOT.CANON) THEN
CSSS       Call dzero(icore(iioo),ND1AA)
       CALL DXAINI(ICORE(IDOO),ICORE(IDOV),ICORE(IDVV),
     &             ICORE(IIOO),ICORE(i0),MXCOR,IUHF,NONHF_TERMS_EXIST)
      ENDIF
      IF(IFLAGS(1).GE.15)THEN
       write(6,*) ' I intermediates after  DXAINI '
       CALL PCCD_DUMPIT(ICORE(IIOO),ICORE(IIVV),ICORE(IIVO),
     &                  ICORE(IIOV),IUHF,"OO")
      ENDIF
C
C JCP-95-2623-91, Eqn. 43, 3-6 terms when (p,q = a,b; INTV1;(r,s,t = m,n,e
C INTV2; (r,s,t = e,f,g), INTV3; (r,s,t = m,n,e,), INTV4 (r,s,t = m,e,f)
C INTV5; (r,s,t = m,e,f), INTV6; (r,s,t = m,n,o). Note that when r,s,t
C is m,e,f or m,n,e the a and b can permute to give two distinct terms.
C That is the reason behind INTV1,INTO3 and INTV4, INTV5 involve same
C r,s,t indices.

      CALL INTV1(ICORE(IIVV),ICORE(i0),MXCOR,IUHF)

      IF(IFLAGS(1).GE.15)THEN
       write(6,*)
       write(6,*) ' I intermediates after  INTV1 '
       CALL PCCD_DUMPIT(ICORE(IIOO),ICORE(IIVV),ICORE(IIVO),
     &                  ICORE(IIOV),IUHF,"VV")
      ENDIF

      CALL INTV2(ICORE(IIVV),ICORE(i0),MXCOR,IUHF)

      IF (IFLAGS(1).GE.15) THEN
          write(6,*) ' I intermediates after INTV2-MO '
          CALL PCCD_DUMPIT(ICORE(IIOO),ICORE(IIVV),ICORE(IIVO),
     &                     ICORE(IIOV),IUHF,"VV")
      ENDIF

      CALL INTV3(ICORE(IIVV),ICORE(i0),MXCOR,IUHF,.TRUE.)

      IF(IFLAGS(1).GE.15)THEN
          write(6,*) ' I intermediates after  INTV3 '
          CALL PCCD_DUMPIT(ICORE(IIOO),ICORE(IIVV),ICORE(IIVO),
     &                     ICORE(IIOV),IUHF,"VV")
      ENDIF

      IF (.NOT.CCD) THEN
        CALL INTV4(ICORE(IIVV),ICORE(i0),MXCOR,IUHF)

        IF(IFLAGS(1).GE.15)THEN
          write(6,*) ' I intermediates after  INTV4 '
          CALL PCCD_DUMPIT(ICORE(IIOO),ICORE(IIVV),ICORE(IIVO),
     &                     ICORE(IIOV),IUHF,"VV")
         ENDIF

        CALL INTV5(ICORE(IIVV),ICORE(i0),MXCOR,IUHF)

        IF(IFLAGS(1).GE.15)THEN
          write(6,*) ' I intermediates after  INTV5 '
          CALL PCCD_DUMPIT(ICORE(IIOO),ICORE(IIVV),ICORE(IIVO),
     &                     ICORE(IIOV),IUHF,"VV")
         ENDIF

        CALL INTV6(ICORE(IIVV),ICORE(i0),MXCOR,IUHF)

        IF(IFLAGS(1).GE.15)THEN
          write(6,*) ' I intermediates after  INTV6 '
          CALL PCCD_DUMPIT(ICORE(IIOO),ICORE(IIVV),ICORE(IIVO),
     &                     ICORE(IIOV),IUHF,"VV")
        ENDIF
      ENDIF 
C
      CALL INTOV1(ICORE(IIOV),ICORE(i0),MXCOR,IUHF)

      IF (IFLAGS(1).GE.15) THEN
         Write(6,*)
         Write(6,*)' I intermediates after INTOV1 '
          CALL PCCD_DUMPIT(ICORE(IIOO),ICORE(IIVV),ICORE(IIVO),
     &                     ICORE(IIOV),IUHF,"OV")
      ENDIF
C
C RESORT FIRST THE G(IA,JB) AMPLITUDES
C
      CALL SORTGAM(ICORE(i0),MXCOR,IUHF)

      CALL INTOV2(ICORE(IIOV),ICORE(i0),MXCOR,IUHF)

      IF (IFLAGS(1).GE.15) THEN
         Write(6,*)' I intermediates after INTOV2 '
         CALL PCCD_DUMPIT(ICORE(IIOO),ICORE(IIVV),ICORE(IIVO),
     &                    ICORE(IIOV),IUHF,"OV")
      ENDIF

      CALL INTOV3(ICORE(IIOV),ICORE(i0),MXCOR,IUHF)
      IF (IFLAGS(1).GE.15) THEN
         Write(6,*)' I intermediates after INTOV3 '
         CALL PCCD_DUMPIT(ICORE(IIOO),ICORE(IIVV),ICORE(IIVO),
     &                    ICORE(IIOV),IUHF,"OV")
      ENDIF

      IF (.NOT.CCD) THEN
        CALL INTOV4(ICORE(IIOV),ICORE(i0),MXCOR,IUHF)
        IF (IFLAGS(1).GE.15) THEN
          Write(6,*)' I intermediates after INTOV4 '
         CALL PCCD_DUMPIT(ICORE(IIOO),ICORE(IIVV),ICORE(IIVO),
     &                    ICORE(IIOV),IUHF,"OV")
        ENDIF
        CALL INTOV5(ICORE(IIOV),ICORE(i0),MXCOR,IUHF,bRedundant)
        IF (IFLAGS(1).GE.15) THEN
         Write(6,*)' I intermediates after INTOV5 '
         CALL PCCD_DUMPIT(ICORE(IIOO),ICORE(IIVV),ICORE(IIVO),
     &                    ICORE(IIOV),IUHF,"OV")
        ENDIF
        CALL INTOV6(ICORE(IIOV),ICORE(i0),MXCOR,IUHF)
        IF (IFLAGS(1).GE.15) THEN
         Write(6,*)' I intermediates after INTOV6 '
         CALL PCCD_DUMPIT(ICORE(IIOO),ICORE(IIVV),ICORE(IIVO),
     &                    ICORE(IIOV),IUHF,"OV")
        ENDIF
      ENDIF 
C
C JCP-95-2623-91, The 43 the 7th term. When the perturb orbitals
C are set to canonical (ROHF non iterative triples), the Eqn, 70
C of JCP-98-8718-93 applies(skip at this point to do latter).
C A very important to note that what we are calculating is I(a,i)
C contirbutions. They will eventually turned into the X(a,i) intermediate.
C The XINT1; (r,s = m,n), XINT2; (r,s = e,f). 09/2004, Ajith Perera.
C
      IF(.NOT.CANON) THEN

       CALL XINT1(ICORE(IIVO),ICORE(IDOO),ICORE(i0),MXCOR,IUHF)

       IF (IFLAGS(1).GE.15) THEN
       Write(6,*)
       Write(6,*)  ' XI intermediates after  XINT1'
       CALL PCCD_DUMPIT(ICORE(IIOO),ICORE(IIVV),ICORE(IIVO),
     &                  ICORE(IIOV),IUHF,"VO")
       ENDIF

       CALL XINT2(ICORE(IIVO),ICORE(IDVV),ICORE(i0),MXCOR,IUHF)

       IF (IFLAGS(1).GE.15) THEN
       Write(6,*)  ' XI intermediates after  XINT2'
       CALL PCCD_DUMPIT(ICORE(IIOO),ICORE(IIVV),ICORE(IIVO),
     &                  ICORE(IIOV),IUHF,"VO")
       ENDIF
      ENDIF

        CALL XINT3(ICORE(IIVO),ICORE(i0),MXCOR,IUHF)

        IF (IFLAGS(1).GE.15) THEN
            Write(6,*)  ' XI intermediates after  XINT3'
            CALL PCCD_DUMPIT(ICORE(IIOO),ICORE(IIVV),ICORE(IIVO),
     &                       ICORE(IIOV),IUHF,"VO")
        ENDIF

C JCP-95-2623-91, The 43 the 7th term. Only when D(O,V) is involoved.
C XINT9 (r,s = e,m). The direct contribution of D(O,V) is applicable
C only when F(a,i) is non zero.
C
      IF(NONHF_TERMS_EXIST) THEN
      CALL XINT9(ICORE(IIVO),ICORE(IDOV),ICORE(i0),MXCOR,IUHF,
     &            .TRUE.)
       IF (IFLAGS(1).GE.15) THEN
       Write(6,*)  ' XI intermediates after  XINT9'
       CALL PCCD_DUMPIT(ICORE(IIOO),ICORE(IIVV),ICORE(IIVO),
     &                  ICORE(IIOV),IUHF,"VO")
       ENDIF
      ENDIF
C
C The formation of X(a,i) intermediates. See Eqns. 49 in  JCP-95-2623-91
C followed by Eqn. 43. The X(a,i) is needed for Z-vector method
C for orbital relaxations.
C
C JCP-95-2623-91, Eqn. 43, 3-6 terms when (p,q = a,i; XINT4;(r,s,t = m,n,o
C XINT5; (r,s,t = n,e,f), XINT6; (r,s,t = e,f,g), XINT7 (r,s,t = m,n,e)
C XINT8; (r,s,t = m,n,e), XINT3; (r,s,t = m,e,f). Notice the symmetry
C in XINT5, XINT3 and XINT4 and XINT7. 09/2004, Ajith Perera.
C
       CALL XINT4(ICORE(IIVO),ICORE(i0),MXCOR,IUHF)

       IF (IFLAGS(1).GE.15) THEN
       Write(6,*)  ' XI intermediates after  XINT4'
       CALL PCCD_DUMPIT(ICORE(IIOO),ICORE(IIVV),ICORE(IIVO),
     &                  ICORE(IIOV),IUHF,"VO")
       ENDIF

       CALL XINT5(ICORE(IIVO),ICORE(i0),MXCOR,IUHF)

       IF (IFLAGS(1).GE.15) THEN
       Write(6,*)  ' XI intermediates after  XINT5'
       CALL PCCD_DUMPIT(ICORE(IIOO),ICORE(IIVV),ICORE(IIVO),
     &                  ICORE(IIOV),IUHF,"VO")
      ENDIF

      IF(.NOT.CCD) THEN

         CALL XINT6(ICORE(IIVO),ICORE(i0),MXCOR,IUHF)

         IF (IFLAGS(1).GE.15) THEN
            Write(6,*)  ' XI intermediates after  XINT6MO'
            CALL PCCD_DUMPIT(ICORE(IIOO),ICORE(IIVV),ICORE(IIVO),
     &                       ICORE(IIOV),IUHF,"VO")
         ENDIF

        CALL XINT7(ICORE(IIVO),ICORE(i0),MXCOR,IUHF,.TRUE.)

        IF (IFLAGS(1).GE.15) THEN
            Write(6,*)  ' XI intermediates after  XINT7'
            CALL PCCD_DUMPIT(ICORE(IIOO),ICORE(IIVV),ICORE(IIVO),
     &                       ICORE(IIOV),IUHF,"VO")
        ENDIF

        CALL XINT8(ICORE(IIVO),ICORE(i0),MXCOR,IUHF)

        IF (IFLAGS(1).GE.15) THEN
            Write(6,*)  ' XI intermediates after  XINT8'
            CALL PCCD_DUMPIT(ICORE(IIOO),ICORE(IIVV),ICORE(IIVO),
     &                       ICORE(IIOV),IUHF,"VO")
        ENDIF

      ENDIF
C
C        CALL XINT3(ICORE(IIVO),ICORE(i0),MXCOR,IUHF)
C
C        IF (IFLAGS(1).GE.15) THEN
C            Write(6,*)  ' XI intermediates after  XINT3'
C            CALL PCCD_DUMPIT(ICORE(IIOO),ICORE(IIVV),ICORE(IIVO),
C     &                       ICORE(IIOV),IUHF,"VO")
C        ENDIF
C
C ADD FINALLY THE F(P,R)*D(Q,R) CONTRIBUTION TO THE  INTERMEDIATES
C
C This term is JCP-95-2623-91, Eqn. 43s first two terms for (p,q = a,b; i,j; a,i) 
C and what is calculated here  is contributions to the I(a,b), I(i,j) and I(i,a). 
C 09/2004, Ajith Perera.
C
      IF(IFLAGS(1).GE.10)THEN
       write(6,*)
       write(6,*) ' D matrices '
       CALL PCCD_DUMPIT(ICORE(IDOO),ICORE(IDVV),ICORE(IDOV),ICORE(IDOV),
     &                  IUHF,"AL")
      ENDIF
C
      CALL PCCD_DFINI(ICORE(IDOO),ICORE(IDVV),ICORE(IDOV),
     &                ICORE(IIOO),ICORE(IIVV),ICORE(IIOV),
     &                ICORE(IIVO),ICORE(i0),MXCOR,IUHF,
     &               .TRUE.)
CSS     &                NONHF_TERMS_EXIST)
C
      IF(IFLAGS(1).GE.15)THEN
       Write(6,*)
       write(6,*) ' I intermediates after F(P,R)*D(Q,R) '
       CALL PCCD_DUMPIT(ICORE(IIOO),ICORE(IIVV),ICORE(IIVO),ICORE(IIOV),
     &                  IUHF,"AL")
      ENDIF

C IN THE CASE OF CANONICAL PERTURBED ORBITALS, FORM HERE
C X(IJ) AND X(AB) AS WELL AS THE FINAL I(IJ) AND I(AB).
C SOLVE FOR D(IJ) AND D(AB) AND AUGMENT X(AI)
C
C Once again the perturb canonical orbitals are choosen for
C ROHF non-iterative triples along with the semicanonical option.
C What that does is to slightly change the formulas for I and as
C a result X intermediates.  09/2004, Ajith Perera.
C
C#ifdef _NOSKIP
      IF(CANON) THEN
C
       IEVAL=i0
       IXOO=IEVAL+IINTFP*(NOCCO(1)+NVRTO(1))
C
C
      IF(IFLAGS(1).GE.15)THEN
       write(6,*)
       Write(6,*) "I, D, X before formxij"
       CALL PCCD_DUMPIT(ICORE(IIOO),ICORE(IDVV),ICORE(IIVO),
     &                  ICORE(IIOV),IUHF,"AL")
       CALL PCCD_DUMPIT(ICORE(IDOO),ICORE(IDVV),ICORE(IIVO),
     &                  ICORE(IIOV),IUHF,"AL")
       CALL PCCD_DUMPIT(ICORE(IXOO),ICORE(IEVAL),ICORE(IIVO),
     &                  ICORE(IIOV),IUHF,"AL")
      ENDIF
C
C The Eqn. 66 of JCP-98-8718-93, first term. Form the D(i,j) contribution.
C
       CALL FORMXIJ(ICORE(IIOO),ICORE(IDOO),
     &              ICORE(IXOO),ICORE(IEVAL),IUHF,QRHF)

      IF(IFLAGS(1).GE.15)THEN
       write(6,*)
       Write(6,*) "I, D, X after formxij"
       CALL PCCD_DUMPIT(ICORE(IIOO),ICORE(IDVV),ICORE(IIVO),
     &                  ICORE(IIOV),IUHF,"AL")
       CALL PCCD_DUMPIT(ICORE(IDOO),ICORE(IDVV),ICORE(IIVO),
     &                  ICORE(IIOV),IUHF,"AL")
       CALL PCCD_DUMPIT(ICORE(IXOO),ICORE(IEVAL),ICORE(IIVO),
     &                  ICORE(IIOV),IUHF,"AL")
      ENDIF
C
C The Eqn. 70 of JCP-98-8718-93, part of the second contribution.
C A contribution to the I(a,i) intermediate. Once again notice that
C only so called non-HF methods contribute.
C
       IF(NONHF_TERMS_EXIST) THEN
        IFAI=IXOO+IINTFP*(ND1AA+ND1BB)
        CALL ADDIA1(ICORE(IXOO),ICORE(IIVO),ICORE(IFAI))
       ENDIF
C
       IXVV=IXOO
C
C The Eqn. 66 of JCP-98-8718-93, second  term. Form the D(a,b) contribution
C
       CALL FORMXAB(ICORE(IIVV),ICORE(IDVV),
     &              ICORE(IXVV),ICORE(IEVAL),IUHF,QRHF)
C
C The Eqn. 70 of JCP-98-8718-93, part of the second contribution.
C A contribution to the I(a,i) intermediate. Once again notice that
C only so called non-HF methods contribute.
C
       IF(NONHF_TERMS_EXIST) THEN
        IFAI=IXVV+IINTFP*(ND2AA+ND2BB)
        CALL ADDIA2(ICORE(IXVV),ICORE(IIOV),ICORE(IFAI))
       ENDIF
C
C These are the contributions to the I(a,i) intermediate that we
C did not included previously (see above). The reason is that the
C form of D(i,j) and D(a,b) is different (see Eqn. 66 in JCP-98-8718-93).
C Since now we have the correct form of D(i,j) and D(a,b) (after
C calls to FORMXIJ and FORMXAB, we can do the now. What is calculated
C is contribution to X(a,i) intermediate via as in Eqn. 70, the
C second term in the first Equation. The XINT1; (r,s = m,n),
C XINT2; (r,s = e,f).
C
       IF(IFLAGS(1).GE.15)THEN
       write(6,*)
       Write(6,*) "I, and D before XINT1 and XINT2"
       CALL PCCD_DUMPIT(ICORE(IIOO),ICORE(IDVV),ICORE(IIVO),
     &                  ICORE(IIOV),IUHF,"AL")
       CALL PCCD_DUMPIT(ICORE(IDOO),ICORE(IDVV),ICORE(IIVO),
     &                  ICORE(IIOV),IUHF,"AL")
       ENDIF

       CALL XINT1(ICORE(IIVO),ICORE(IDOO),ICORE(i0),MXCOR,IUHF)
       CALL XINT2(ICORE(IIVO),ICORE(IDVV),ICORE(i0),MXCOR,IUHF)

       IF(IFLAGS(1).GE.15)THEN
       Write(6,*) "I, and D after XINT1 and XINT2"
       CALL PCCD_DUMPIT(ICORE(IIOO),ICORE(IDVV),ICORE(IIVO),
     &                  ICORE(IIOV),IUHF,"AL")
       ENDIF
C
C  ADD THE CONTRIBUTION DUE TO DOO AND DVV TO IOO
C
C This is the 7th term in the Eqn. 43 of JCP-95-2623-91.
C As discussed above, now we have the correct form of the
C D(i,j) and D(a,b) for perturb orbitals canonical (FORMXIJ,
C FORMXAB), so we can do them here now.
C
       CALL DXAINI(ICORE(IDOO),ICORE(IDOV),ICORE(IDVV), 
     &             ICORE(IIOO),ICORE(i0),MXCOR,IUHF,
     &             NONHF_TERMS_EXIST)
C
      ENDIF
C
      IF(IFLAGS(1).GE.10)THEN
       write(6,*)
       Write(6,*) "Intermediates before solving Z vec eqns"
       CALL PCCD_DUMPIT(ICORE(IIOO),ICORE(IIVV),ICORE(IIVO),ICORE(IIOV),
     &                  IUHF,"AL")
      ENDIF
C
C FORM X(A,I) REQUIRED FOR SOLVING THE Z-VECTOR EQUATIONS
C
C The relevent paper is JCP-95-2639-91 and the equations are 15 and
C 16 for the QRHFX. The equations 24 applies for QRHFZ. Beware that
C QRHFX and QRHFZ are doing specilized tasks. Unless you know
C exactly what you doing do not turn on the QRHF flag. All the
C other methods, the equation 48 of JCP-95-2623-91 applies. Note
C the memory location IIVO has the I(a,i) contributions and the
C and the  IIOV has the I(i,a) contributions. Now, form the
C I(i,a) - I(a,i) and copy it to the IXOV memory location. A comment
C about the test with NONHF_TERMS_EXIST. If true what does this tells
C is that the user has only specified NONHF flag, but analysis of the
C Fock matrix revels that it indeed a NONHF type fock matrix. Then
C follow the usual path. If that is not the case, user is doing
C something strange but not invalid. In that case, memory for the D(a,i)
C (direct contributions is not required). Since memory is allocated
C based on NONHF=.TRUE. before the analysis of the Fock matrix, proper
C location for IXOV is IDOV (see the memory allocation section for D(a,i) at the
C begining). Obviously, if this happen memory allocated for D(a,i) is simply a
C waste, but it is such small amount and can be safely ignored.
C 09/2004, Ajith Perera.
C
      IF(QRHF) THEN
C
C FOR NON-HF WAVE FUNCTION THIS IS A LITTLE BIT COMPLICATED
C NONHFX EXTRACTS THE RELEVANT TERMS FOR QRHF-CCSD CALCULATIONS
C
C ALLOCATE FIRST ADDITIONAL MEMORY FOR THE X1A AND X1B ARRAYS
C
C
       IX1=IIVV-2*IINTFP*N1I
       MXCOR=MXCOR-IINTFP*2*N1I
       IX2=IX1-2*IINTFP*NA2
       MXCOR=MXCOR-IINTFP*2*NA2
       ITOP=IX2
       CALL IZERO(ICORE(IX1),2*IINTFP*N1I)
       CALL IZERO(ICORE(IX2),2*IINTFP*NA2)
C
       CALL QRHFX(ICORE(IIOO),ICORE(IIVV),ICORE(IIOV),ICORE(IIVO),
     &            ICORE(IIOO+IINTFP*ND1AA),ICORE(IIVV+IINTFP*ND2AA),
     &            ICORE(IIOV+IINTFP*NTAA),ICORE(IIVO+IINTFP*NTAA),
     &            ICORE(IXOV),ICORE(IXOV+IINTFP*NAI),ICORE(IX1),
     &            ICORE(IX1+IINTFP*N1I),ICORE(IX2),
     &            ICORE(IX2+IINTFP*NA2),ICORE(i0))
C
C  FORM HERE NOW THE Z1I (OR ZA1) CONTRIBUTION TO THE DENSITY
C  AND CALCULATE XAI(TWIDDLE)
C
       IEV=ITOP-IINTFP*(NOCCO(1)+NVRTO(1))
       CALL QRHFZ(ICORE(IXOV),ICORE(IXOV+IINTFP*NAI),ICORE(IX1),
     &            ICORE(IX1+IINTFP*N1I),ICORE(IX2),
     &            ICORE(IX2+IINTFP*NA2),ICORE(IEV),ICORE(i0),
     &            MXCOR,IUHF)
C
      ELSE
       IF (.NOT.NONHF_TERMS_EXIST) IXOV=IDOV
       CALL DCOPY((NTAA+IUHF*NTBB),ICORE(IIVO),1,ICORE(IXOV),1)
       CALL SAXPY((NTAA+IUHF*NTBB),ONE,ICORE(IIOV),1,ICORE(IXOV),1)
       ITOP=IIVV
      ENDIF
#ifdef _DEBUG_LVL0
      nocc=nocco(1)
      nvrt=nvrto(1)
      Write(6,"(a)") "The antisymmetrized I(pq):The OV block"
      call output(icore(ixov),1,nocc,1,nvrt,nocc,nvrt,1)
#endif

      IF (NONHF_TERMS_EXIST .AND. .NOT. TRULY_NONHF) THEN 
 
          Write(6,"(2x,a,a)")"Warning!! Enter into exclusively reserved"
     &                     ," block that is only for non-digonal"
          Write(6,"(12x,a)") "elements in FOCK but not ROHF or QRHF." 

          CALL MAKDIA1_4STD_NONCAN(ICORE(IXOV),ICORE(i0),MXCOR,IUHF,
     &                             bRedundant)
      ELSE
          CALL MAKDIA1(ICORE(IXOV),ICORE(i0),MXCOR,IUHF,bRedundant)
      ENDIF
C
C  FOR HF CASES MAKDIA1 RETURNS THE OCCUPIED-VIRTUAL BLOCK OF THE
C  DENSITY MATRIX AT ADDRESS  IDOV (=IXOV), FOR NON-HF CASES 
C  IDOV NE IXOV
C
C  FORM NOW FINAL INTERMEDIATES AND SAVE QUANTITIES
C
C  NOTE, THAT DAINI TAKES CARE OF ALL THE ROHF STUFF
C
       CALL DAIINI(ICORE(IXOV),ICORE(IIOO),ICORE(IIOV),ICORE(i0),
     &             MXCOR,IUHF,NONHF_TERMS_EXIST)
C
C FOR ROHF, ADD DIRECT AND IN-DIRECT CONTRIBUTIONS
C
      IF(ROHF .OR. NONHF_TERMS_EXIST) THEN
        IF (IUHF .NE. 0) THEN
            CALL SAXPY(NTAA+NTBB,ONE,ICORE(IXOV),1,ICORE(IDOV),1)
        ELSE
            CALL SAXPY(NTAA,ONE,ICORE(IXOV),1,ICORE(IDOV),1)
        ENDIF 
      ENDIF
C
C The following two FRMFUL calls built and write the complete
C density matrix and the I intermediates, depending on the last
C argument. The argument 'DEN' will write out RELDENSA and RELDENSB records.
C The argument INT will write out the I intermediate. The I intermediate
C is required only for the gradient calculations. 09/2004, Ajith Perera.
C
C#endif 

      CALL PUTREC(20,"JOBARC","ROTGRDOO",ND1AA+IUHF*ND1BB,
     &            ICORE(IIOO))
      CALL PUTREC(20,"JOBARC","ROTGRDVV",ND2AA+IUHF*ND2BB,
     &            ICORE(IIVV))
      CALL PUTREC(20,"JOBARC","ROTGRDOV",NTAA+IUHF*NTBB,
     &            ICORE(IIOV))
      CALL PUTREC(20,"JOBARC","ROTGRDVO",NTAA+IUHF*NTBB,
     &            ICORE(IIVO))

      NBAS=NOCCO(1)+NVRTO(1)
      IDFULL=ITOP-(1+IUHF)*NBAS*NBAS*IINTFP
      MXCOR1=MXCOR-(1+IUHF)*NBAS*NBAS*IINTFP
      CALL FRMFUL(ICORE(IDOO),ICORE(IDVV),ICORE(IDOV),
     &            ICORE(IDFULL),ICORE(i0),MXCOR,NBAS,IUHF,'DEN')

      IIFULL=ITOP-(1+IUHF)*NBAS*NBAS*IINTFP
      MXCOR=MXCOR1-(1+IUHF)*NBAS*NBAS*IINTFP
      CALL FRMFUL(ICORE(IIOO),ICORE(IIVV),ICORE(IIOV),
     &            ICORE(IIFULL),ICORE(i0),MXCOR,NBAS,IUHF,'INT')

      RETURN
      END 

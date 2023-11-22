      SUBROUTINE ONEINT(PROPTY,CORE,MAXCOR,
     &                  ENERKE,GRADKE,HESSKE,ENERNA,GRADNA,
     &                  HESSNA,GRADFS,HESFS2,DIPME,DIPMN,
     &                  DDIPE,DDIPN,QUADME,QUADMN,
     &                  DQUADE,DQUADN,NCOORD,DOSCF)
C
C  DRIVER FOR EVALUATION OF ONE-ELECTRON INTEGRALS
C  AND INTEGRAL DERIVATIVES
C  MODIFIED TO ALLOW IN ADDITION CALCULATION OF ONE
C  ELECTRON INTEGRALS NEEDED FOR NMR SHIFT CALCULATIONS
C  AND INDIRECT SPIN-SPIN COUPLING CONSTANTS
C
CEND
C
C  ADAPTED FROM ABACUS, SEPTEMBER/90 JG
C  GIAO INTEGRALS, AUGUST/91 JG
C  FC,SD, AND DSO INTEGRALS,  APRIL/93 JG
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL PROPTY,DOSCF
      LOGICAL MOLGRD,MOLHES,DIPDER,POLAR,INPTES,VIB,RESTAR,
     &        DOWALK,GDALL,FOCK1,H2MO
      LOGICAL         SKIP, DTEST, DIFINT, NODC, NODV, DIFDIP
      LOGICAL ELECT,GEOM,MAGNET,GRAD,INTWRIT,FOCK,MAGNET2,GIAO
      LOGICAL JSO,JFC,JSD
C
      DIMENSION CORE(MAXCOR)
      DIMENSION GRADKE(NCOORD),HESSKE(NCOORD,NCOORD),
     &          GRADNA(NCOORD),HESSNA(NCOORD,NCOORD),
     &          GRADFS(NCOORD),HESFS2(NCOORD,NCOORD),
     &          DIPME(3),DDIPE(3,NCOORD),QUADME(6),
     &          DQUADE(6,NCOORD),DIPMN(3),DDIPN(3,NCOORD),
     &          QUADMN(6),DQUADN(6,NCOORD)
    
C
      COMMON/ABAINF/IPRDEF,MOLGRD,MOLHES,DIPDER,POLAR,
     &              INPTES,VIB,RESTAR,DOWALK,GDALL,FOCK1,H2MO
      COMMON /CBIONE/ IPRINT, SKIP, MAXDIF, DTEST, IDATOM,
     &                IDCOOR, DIFINT, NODC, NODV, DIFDIP
      COMMON/PROP/ELECT,MAGNET,GEOM,MAGNET2,GIAO
      COMMON/OPTION/GRAD,INTWRIT,FOCK
      COMMON/JNMR/JSO,JFC,JSD
C
      LMAX   = 0
C
      IF(ELECT.OR.GEOM.OR.JFC.OR.JSD) THEN
C
C  CALCULATE INTEGRALS NEEDED FOR GEOMETRIC AND ELECTRIC PERTURBATIONS
C
       CALL ONEDRV(PROPTY,MAXDIF,FOCK,IDATOM,IDCOOR,
     &             DIFINT,NODC,NODV,DIFDIP,LMAX,CORE,
     &             MAXCOR, ENERKE,GRADKE,HESSKE,ENERNA,
     &             GRADNA,HESSNA,GRADFS,HESFS2,DIPME,DIPMN,
     &             DDIPE,DDIPN,QUADME,DQUADE,QUADMN,
     &             DQUADN,NCOORD,DOSCF)
C
      ELSE IF(MAGNET) THEN
C
C  CALCULATE INTEGRALS NEEDED FOR MAGNETIC PERTURBATIONS
C
       CALL NMR1DR(PROPTY,MAXDIF,FOCK,IDATOM,IDCOOR,
     &             DIFINT,NODC,NODV,DIFDIP,LMAX,CORE,
     &             MAXCOR,HESSKE,DIPME,DIPMN,NCOORD)
C
       IF(JSO) THEN
C
C  CALCULATE INTEGRALS NEEDED FOR DIAMAGNETIC SO CONTRIBUTION
C
c       CALL DSODRV(PROPTY,MAXDIF,FOCK,IDATOM,IDCOOR,
c     &             DIFINT,NODC,NODV,DIFDIP,LMAX,CORE,
c     &             MAXCOR,HESSNA,NCOORD)
C
       ENDIF
      ELSE
C
       write(*,*) ' type of perturbation has not been specified !'
C
      ENDIF       
C
C ALL DONE, RETURN
C
      RETURN
      END

      SUBROUTINE Scrnc_built_respn_den(ICORE, MAXCOR, IUHF) 
C
C Calculate the response density matrix for second-order response
C properties.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      INTEGER POP, VRT
      DIMENSION ICORE(MAXCOR)
      LOGICAL MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD_DUPL,UCC
      LOGICAL DENS,GRAD,QRHF,NONHF,ROHF,SEMI,CANON,TRIP1,
     &             TRIP2
      LOGICAL CCSD,MBPT,PARTEOM,NODAVID,AOLADDER,RCCD,DRCCD
      LOGICAL CC
      COMMON /FILES/ LUOUT, MOINTS
      COMMON /INFO/ NOCCO(2), NVRTO(2)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYM/ POP(8,2),VRT(8,2),NTAA,NTBB,ND1AA,ND1BB,
     &            ND2AA,ND2BB
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /EOMINFO/CCSD,MBPT,PARTEOM,NODAVID,AOLADDER,RCCD,DRCCD
      COMMON/METH/MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD_DUPL,
     &            UCC
      COMMON/DERIV/DENS,GRAD,QRHF,NONHF,ROHF,SEMI,CANON,TRIP1,
     &             TRIP2
      COMMON/LISTDENS/LDENS
C
      DATA ONE /1.00D+00/
C
      Write(6,"(a)") " --------Warning---------"
      Write(6,"(a,a,a)") " Calculation level is set to CC in",
     &                 " scrnc_built_respn_den.F for"
      Write(6,"(a)")   " CALC={CCSD,RCCD,DRCCD}"
      Write(6,"(a,a)") " Ground state denisty is valid only for HF",
     &                 " (F(pq)=delta(p,q)) methods"
      Write(6,"(a)") " --------Warning---------"
      Write(6,*)

      CALL SCRNC_SET_DENS_VARS(IUHF)

      CC    =  (CCSD .OR. RCCD .OR. DRCCD)

      MXCOR = MAXCOR - MOD(MAXCOR, 2)
      LDENS = 160
C     
C Memory for the occupied-occupied block.
C     
      IO = 1
      IDOO  = MXCOR - (ND1AA+IUHF*ND1BB)*IINTFP + 1
      MXCOR = MXCOR - (ND1AA + IUHF*ND1BB)*IINTFP
C     
C Memory for the virtual-virtual block.
C     
      IDVV  = IDOO - (ND2AA + IUHF*ND2BB)*IINTFP 
      MXCOR = MXCOR - (ND2AA + IUHF*ND2BB)*IINTFP
C     
C Memory for the virtual-occupied block.
C     
      IDOV  = IDVV - (NTAA + IUHF*NTBB)*IINTFP
      MXCOR = MXCOR - (NTAA + IUHF*NTBB)*IINTFP
C     
C Memory for the XOV intermediate
C For non-hf we have to allocate additional memory for the XOV intermediate.
C
      IF (QRHF .OR. NONHF .OR. ROHF) THEN
         IXOV  = IDOV - (MAX(NTAA,NTBB) + IUHF*MAX(NTAA,NTBB))*IINTFP
         MXCOR = MXCOR - (MAX(NTAA,NTBB) + IUHF*MAX(NTAA,NTBB))*IINTFP
         LXOV  = (1 + IUHF)*MAX(NTAA,NTBB)
      ELSE
         IXOV = IDOV
         LXOV = 0
      ENDIF
C     
C Memory for the IOV intermediate 
C     
      IIOV  = IXOV - (NTAA+IUHF*NTBB)*IINTFP
      MXCOR = MXCOR - (NTAA+IUHF*NTBB)*IINTFP
C     
C Memory for the IVO intermediate
C     
      IIVO  = IIOV -  (NTAA + IUHF*NTBB)*IINTFP
      MXCOR = MXCOR - (NTAA + IUHF*NTBB)*IINTFP
C     
C Memory for the IOO intermediate
C     
      IIOO  = IIVO - (ND1AA + IUHF*ND1BB)*IINTFP
      MXCOR = MXCOR - (ND1AA + IUHF*ND1BB)*IINTFP
C     
C     MEMORY FOR THE IVV INTERMEDIATE
C     
      IIVV = IIOO - (ND2AA + IUHF*ND2BB)*IINTFP
      MXCOR = MXCOR - (ND2AA+IUHF*ND2BB)*IINTFP
C     
C Set ITOP so that non-HF will work
C     
      ITOP = IIVV
C     
C Zero the allocated space for the intermediates IOO and IVV
C     
      LENGTH = LXOV+3*NTAA + 2*ND1AA + 2*ND2AA +
     &         IUHF*(3*NTBB+2*ND1BB+2*ND2BB)
C     
      CALL IZERO(ICORE(IIVV), LENGTH*IINTFP)
C     
C Compute the occupied-occupied block of the density matrix.
C
      CALL SCRNC_DENS_OO(ICORE(IDOO), ICORE(IO), MXCOR, IUHF)

C Compute the virtual-virtual block of the density matrix.
C      
       CALL SCRNC_DENS_VV(ICORE(IDVV),ICORE(IO),MXCOR,IUHF)

C
C Compute the direct contribution to the virtual-occupied block of the 
C density matrix
C
       CALL SCRNC_DENS_VO(ICORE(IDOV),ICORE(IO),MXCOR,IUHF)

#ifdef _DEBUG_LVL0
CSSS       CALL DUMPIT(ICORE(IDOO),ICORE(IDVV),ICORE(IDOV),ICORE(IIOV),
CSSS     &              IUHF)
       call checksum("Doo(AA)",ICORE(IDOO),ND1AA)
       if (iuhf.ne.0) call checksum("Doo(BB)",ICORE(IDOO+ND1AA),ND1BB)
       call checksum("Dvv(AA)",ICORE(IDVV),ND2AA)
       if (iuhf.ne.0) call checksum("Dvv(BB)",ICORE(IDVV+ND2AA),ND2BB)
       call checksum("Dvo(AA)",ICORE(IDOV),NTAA)
       if (iuhf.ne.0) call checksum("Dvo(BB)",ICORE(IDOV+NTAA),NTBB)
#endif 
C 
C Write out the reduced density matrix to the disk.
C
       CALL SCRNC_SAVD(ICORE(IDOO),ICORE(IDOO+IINTFP*ND1AA),
     &                 ICORE(IDVV),ICORE(IDVV+IINTFP*ND2AA),
     &                 ICORE(IDOV),ICORE(IDOV+IINTFP*NTAA),
     &                 IUHF)
C

      RETURN
      END


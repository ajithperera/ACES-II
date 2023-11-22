      SUBROUTINE GENRSPNSDEN(ICORE, MAXCOR, IUHF) 
C
C Calculate the response density matrix for second-order response
C properties.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      INTEGER POP, VRT
      DIMENSION ICORE(MAXCOR)
C
      LOGICAL MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC
      LOGICAL DENS,GRAD,QRHF,NONHF,ROHF,SEMI,CANON,TRIP1,TRIP2

C      
      COMMON /FILES/ LUOUT, MOINTS
      COMMON /INFO/ NOCCO(2), NVRTO(2)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYM/ POP(8,2),VRT(8,2),NTAA,NTBB,ND1AA,ND1BB,
     &            ND2AA,ND2BB
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),DIRPRD(8,8)
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON /METH/ MBPT2,MBPT3,M4DQ,M4SDQ,M4SDTQ,CCD,QCISD,CCSD,UCC
      COMMON /DERIV/DENS,GRAD,QRHF,NONHF,ROHF,SEMI,CANON,TRIP1,
     &              TRIP2
C
      DATA ONE /1.00D+00/
C
C Fill the common block /METH/
C
      CALL GETDENMETHDS(IUHF)
C
      MXCOR = MAXCOR - MOD(MAXCOR, 2)
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
      CALL DENSOO(ICORE(IDOO), ICORE(IO), MXCOR, IUHF)
C      CALL DUMPIT(ICORE(IDOO),ICORE(IDVV),ICORE(IDOV),ICORE(IIOV),
C     &              IUHF)
#ifdef _DEBUG_LVL0
C      call output(icore(idoo),1,nd1aa,1,nd1aa,nd1aa,nd1aa,1)
#endif 
C      
C Compute the virtual-virtual block of the density matrix.
C      
       CALL DENSVV(ICORE(IDVV),ICORE(IO),MXCOR,IUHF)
C       CALL DUMPIT(ICORE(IDOO),ICORE(IDVV),ICORE(IDOV),ICORE(IIOV),
C     &              IUHF)
#ifdef _DEBUG_LVL0
C      call output(icore(idvv),1,nd2aa, 1,nd2aa,nd2aa,nd2aa,1)
#endif
C
C Compute the direct contribution to the virtual-occupied block of the 
C density matrix
C
       CALL DENSVO(ICORE(IDOV),ICORE(IO),MXCOR,IUHF)
C       CALL DUMPIT(ICORE(IDOO),ICORE(IDVV),ICORE(IDOV),ICORE(IIOV),
C     &              IUHF)
C
#ifdef _DEBUG_LVL0
C      call output(icore(idov),1,ntaa, 1,ntaa,ntaa,ntaa,1)
#endif
C 
C Write out the reduced density matrix to the disk.
C
       CALL SAVD(ICORE(IDOO),ICORE(IDOO+IINTFP*ND1AA),
     &           ICORE(IDVV),ICORE(IDVV+IINTFP*ND2AA),
     &           ICORE(IDOV),ICORE(IDOV+IINTFP*NTAA),
     &           IUHF)
C

      RETURN
      END


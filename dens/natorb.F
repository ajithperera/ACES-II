      SUBROUTINE NATORB(DENS,BUF,NMO,MAXCOR,IUHF,NDROP0)
C
C DIAGONALIZES RELAXED DENSITY MATRIX TO OBTAIN OCCUPATION NUMBERS
C  AND NATURAL ORBITALS
C
C
C  DENS - NMOxNMO MATRIX TO HOLD RELAXED DENSITY
C   BUF - A SCRATCH VECTOR OF LENGTH 3*NBAS*NBAS, WHERE NBAS IS THE
C            NUMBER OF BASIS FUNCTIONS
C   NMO - THE NUMBER OF ACTIVE MOLECULAR ORBITALS
C  IUHF - RHF/UHF FLAG 
C
CEND
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C  
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
C
      DIMENSION DENS(NMO*NMO),BUF(1), NBF4IRREP(8),NOCC(8,2)
C
      DATA ONE   /1.0/
      DATA ZILCH /0.0/ 
C
      IONE=1
      WRITE(6,1000)
      WRITE(6,1001)
      WRITE(6,1000)
C SG 2/96
      IF (IUHF.EQ.1) THEN
       WRITE(6,1002)
       WRITE(6,1000)
      ENDIF
C
C This record is simply to indicate that natural orbitals are formed
C and available to read from any subsequent program(s). See vscf.F 
C for one such scenario, 11/09, Ajith Perera. 

      CALL PUTREC(0, "JOBARC", "EXTORBS", IONE, IONE)
C
      CALL GETREC(20,'JOBARC','RELDENSA',NMO*NMO*IINTFP,DENS)
      CALL GETREC(20,'JOBARC','NBASTOT ',IONE,NBAS)
      CALL GETREC(20,'JOBARC','NAOBASFN',IONE,NAOBFNS)
#ifdef _DEBUG_LVL0
      Write(6,"(a,(1x,I4,1x,I4))") "NMO, NBAS :", NMO, NBAS
      Write(6, "(a)") "The Alpha density"
      call output(dens, 1, nbas, 1, nbas, nbas, nbas, 1)
#endif 
      IOFFB=NBAS*NBAS+1
      IOFFB2=NBAS*NBAS+IOFFB
      IOFFB3=NBAS*NBAS+IOFFB2
      IOFFB4=NBAS*NBAS+IOFFB3

      CALL EIG(DENS,BUF,NMO,NMO,-1)
      WRITE(6,*)
      WRITE(6,2000)(DENS(I),I=1,NMO*NMO,NMO+1)
      TRACE=SSUM(NMO,DENS,NMO+1)
      WRITE(6,2001)TRACE
      CALL PUTNATORBOCCS(DENS,NBAS,"Alpha") 
      CALL GETREC(20,'JOBARC','SCFEVECA',NMO*NBAS*IINTFP,BUF(IOFFB))
      CALL XGEMM('N','N',NBAS,NMO,NMO,ONE,BUF(IOFFB),NBAS,
     &           BUF,NMO,ZILCH,BUF(IOFFB2),NBAS)
#ifdef _DEBUG_LVLM
      Write(6, "(a)") "The Alpha natural orbitals"
      call output(buf(Ioffb2), 1, nbas, 1, nbas, nbas, nbas, 1)
      call checksum("natorba:",buf(Ioffb2),nbas*nbas)
#endif 
      CALL PUTREC(20,'JOBARC','AONTORBA',NBAS*NMO*IINTFP,BUF(IOFFB2))
      CALL SCOPY (NMO,DENS,NMO+1,BUF(IOFFB),1)
      CALL PUTREC(20,'JOBARC','OCCNUM_A',NMO*IINTFP,BUF(IOFFB)) 
c antiquated.
c      IF(IFLAGS(42).NE.0)THEN
c        CALL NEWEVC(BUF(IOFFB2),BUF(1),NBAS,1)
c        CALL PUTREC(20,'JOBARC','SCFEVCA0',NBAS*NMO*IINTFP,BUF(IOFFB2))
c        IFLAGS(16)=0
c        IFLAGS(38)=1
c      ENDIF
C SG 2/96
C changed so that Beta spin natural orbitals are not printed out for
C   RHF calculation
      IF (IUHF.EQ.1) THEN
       WRITE(6,1000)
       WRITE(6,1003)
       WRITE(6,1000)
       CALL GETREC(20,'JOBARC','RELDENSB',NMO*NMO*IINTFP,DENS)
#ifdef _DEBUG_LVLM
      Write(6,"(a,(1x,I4,1x,I4))") "NMO, NBAS :", NMO, NBAS
      Write(6, "(a)") "The Beta density"
      call output(dens, 1, nbas, 1, nbas, nbas, nbas, 1)
#endif
       CALL EIG(DENS,BUF,NMO,NMO,-1)
       CALL GETREC(20,'JOBARC','SCFEVECB',NMO*NBAS*IINTFP,BUF(IOFFB))
       CALL XGEMM('N','N',NBAS,NMO,NMO,ONE,BUF(IOFFB),NBAS,
     &            BUF,NMO,ZILCH,BUF(IOFFB3),NBAS)
#ifdef _DEBUG_LVLM
       write(6,*)
       write(6, "(a)") "The Beta natural orbitals"
       call output(buf(ioffb3), 1, nbas, 1, nbas, nbas, nbas, 1)
      call checksum("natorbb:",buf(Ioffb3),nbas*nbas)
#endif 
       CALL PUTREC(20,'JOBARC','AONTORBB',NBAS*NMO*IINTFP,BUF(IOFFB3))
       CALL SCOPY (NMO,DENS,NMO+1,BUF(IOFFB),1)
       CALL PUTREC(20,'JOBARC','OCCNUM_B',NMO*IINTFP,BUF(IOFFB)) 
c antiquated.
c       IF(IFLAGS(42).NE.0)THEN
c        CALL NEWEVC(BUF(IOFFB2),BUF(1),NBAS,2)
c        CALL PUTREC(20,'JOBARC','SCFEVCB0',NBAS*NMO*IINTFP,BUF(IOFFB2))
c        IFLAGS(16)=0
c        IFLAGS(38)=1
c       ENDIF
       WRITE(6,*)
       WRITE(6,2000)(DENS(I),I=1,NMO*NMO,NMO+1)
       TRACE=SSUM(NMO,DENS,NMO+1)
       WRITE(6,2001)TRACE
      ENDIF

      CALL GETREC(20, "JOBARC", 'NIRREP  ',  1,      NIRREP)
      CALL GETREC(20, "JOBARC", 'NUMBASIR',  NIRREP, NBF4IRREP)
      Call GETREC(20, "JOBARC", 'OCCUPYA0',  NIRREP, NOCC(1,1))
      IF (IUHF .NE. 0) THEN
         CALL GETREC(20, "JOBARC", 'OCCUPYB0',  NIRREP, NOCC(1,2))
      ELSE
         CALL ICOPY(8, NOCC(1,1), 1, NOCC(1,2), 1)
      ENDIF 

      WRITE(6,1000)

      If (Ndrop0 .Eq. 0) Then
      Do Ispin = 1, IUhf+1

C Reorder the alpha natural orbital to match with the order that NEWMOS
C in SCF is written

         Iscr1 = Ioffb4
         IScr2 = Iscr1 + Nmo*Nmo
         Iscr3 = Iscr2 + Nmo*Nmo
         Iscr4 = Iscr3 + Nmo+1
         Ikeep = Iscr4 + Nmo+1
         Iend  = Ikeep + 2*Nmo*Nmo
         If (Iend .Gt. Maxcor) Call Insmem("-Natorb",Iend,Maxcor)

         If (Ispin .EQ. 1) Ioff1 = Ioffb2 
         If (Ispin .EQ. 2) Ioff1 = Ioffb3
         Call Process_natorbs(Buf(Ioff1),Buf(Iscr1),Buf(Iscr2),
     +                        Buf(Iscr3),Buf(Iscr4),Nmo,Length,
     +                        Ispin)

         Call Dcopy(Length,Buf(Ioff1),1,Buf(Ikeep),1)
         If (Ispin .Eq. 2) Then
            Ioff = Length + Ikeep 
            Call Dcopy(Nmo*Nmo,Buf(Ioff1),1,Buf(ioff),1)
         Endif 
      Enddo 

      CALL PUTNATORBS(BUF(Ikeep), NBAS*NBAS, IUHF, NIRREP, NBF4IRREP,
     &                NOCC)

      CALL PUTNATORBOCCS(DENS,NBAS,"Beta ") 
      Else
      Write(6,"(2a)") "          ------------------Alert--------------",
     &                "----"
      Write(6,"(2x,2a)") "Natural orbitals and occupations are not",
     &                   " written for frozen orbital calculations."   
      Write(6,"(a,a)") "          ------------------------------------",
     &                "-------"
      Write(6,*)
      Endif 
C
      ISIZE=100
      CALL PUTREC(20,'JOBARC','IFLAGS  ',ISIZE,IFLAGS)
      RETURN
1000  FORMAT(T3,70('-'))
1001  FORMAT(T18,'Natural orbital occupation numbers')
1002  FORMAT(T29,'Alpha spin')
1003  FORMAT(T29,'Beta spin')
2000  FORMAT((8(2X,F7.5)))
2001  FORMAT(T3,'Trace of density matrix : ',F14.10,'.')
      END

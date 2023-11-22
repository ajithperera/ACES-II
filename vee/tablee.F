#include "flags.h"
      SUBROUTINE TABLEE(EIGVAL,EIGVAL_T,OSCSTR,NATURE,BGN,BGN_IRP,
     &                  END,END_IRP,METHOD)
C
C  THE CALCULATION OF EXCITATION ENERGIES IS SUMMARIZED
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      INTEGER DIRPRD,BGN,END,BGN_IRP,END_IRP,S,T
      Logical Triplet,G_TDENS,CORE_SEARCH,DOUBLE_EES,EXT_PRINTING
      Logical pCCD
      CHARACTER*1 NATURE(100,8)
      DIMENSION OSC(1024)
      DIMENSION OSC_QM(1024)
      DIMENSION OSC_OM(1024)
      DIMENSION OSC_AM(1024)
      DIMENSION OSC_MQ(1024)
      DIMENSION OSC_DM(1024)

      DIMENSION EIGVAL(100,8),EIGVAL_T(100,8),OSCSTR(100,8),BGN(100,8),
     &          BGN_IRP(100,8),END(100,8),END_IRP(100,8)

      LOGICAL MBPT2, CC,CCD,RCCD,DRCCD,LCCD,LCCSD,EOM_TRPS
C
      COMMON /SYMINF/ NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON/CALCINFO/NROOT(8)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/FLAGS/IFLAGS(100)
      COMMON/FLAGS2/IFLAGS2(500)
      COMMON/REFTYPE/MBPT2,CC,CCD,RCCD,DRCCD,LCCD,LCCSD
C
      DATA FACTEV, FACTCM /27.2113957D0,  2.19474625D5/
      DATA ZILCH /0.0D0/
C
C
#ifdef _DEBUG_LVLM
      do irrep = 1, nirrep
      do iroot = 1, nroot(irrep)
      Write(6,"(a,1x,a)") "Singlet or triplet?:", Nature(Iroot,Irrep)
      Write(6,"(a,i3)") "The origin      irrep:", Bgn_irp(Iroot,Irrep)
      Write(6,"(a,i3)") "The destination irrep:", end_irp(Iroot,Irrep)
      Write(6,"(a,i3)") "The origin        orb:", Bgn(Iroot,Irrep)
      Write(6,"(a,i3)") "The destination   orb:", end (Iroot,Irrep)
      enddo 
      enddo
#endif
      Triplet = .FALSE.
      EOM_TRPS = (iflags(87) .EQ. 11)

      CALL GETREC(0,"JOBARC","PCCD_RUN",ILENGTH,JUNK)
      pCCD = .FALSE.
      IF (ILENGTH .GT. 0) pCCD = .TRUE.

      WRITE(6,*)
      IF (METHOD .EQ. 3) THEN
        IF (CC .AND. (.NOT. CCD)) THEN
          IF (PCCD) THEN
             WRITE(6,*) ' Summary of pCCD-EOM-CCSD excitation energies'
          ELSE
             WRITE(6,*) ' Summary of EOM-CCSD excitation energies'
          ENDIF 
        ELSE IF (MBPT2) THEN
          WRITE(6,*) ' Summary of EOM-MBPT(2) excitation energies'
        ELSE IF (CCD) THEN
          IF (PCCD) THEN
             WRITE(6,*) ' Summary of pCCD-EOM-CCSD excitation energies'
          ELSE
             WRITE(6,*) ' Summary of CCD excitation energies'
          ENDIF 
        ELSE IF (RCCD) THEN
          WRITE(6,*) ' Summary of rCCD excitation energies'
        ELSE IF (DRCCD) THEN
          WRITE(6,*) ' Summary of drCCD excitation energies'
        ELSE IF (LCCD) THEN
          WRITE(6,*) ' Summary of LCCD excitation energies'
        ELSE IF (LCCSD) THEN
          WRITE(6,*) ' Summary of LCCSD excitation energies'
        ENDIF
      ELSE IF (METHOD .EQ. 7) THEN
        IF (CC) THEN
       WRITE(6,*) ' Summary of P-EOM-CCSD excitation energies'
        ELSE
       WRITE(6,*) ' Summary of P-EOM-MBPT(2) excitation energies'
        ENDIF
      ELSE IF (METHOD .EQ. 8) THEN
        IF (CC) THEN
          WRITE(6,*) ' Summary of BWPT(2)-EOM-CCSD excitation energies'
        ELSE
          WRITE(6,*) ' Summary of BWPT(2)-EOM-MBPT excitation energies'
        ENDIF
      ENDIF

#ifdef _NOSKIP
        WRITE(6,1300)
         WRITE(6, 2015)
 2015    FORMAT(2X,'SYM. ','Origin',1x,'Destination',1x,'Nature', 
     &      3X, 'EE(eV)',
     $      3x, 'EE(cm-1)',4X,'OSC. STR.',3X, 'TOTAL ENERGY')
         WRITE(6,*)'  ------------------------------------',
     $      '---------------------------------------'
         WRITE(6,*)

         DO IRREP = 1, NIRREP
           DO  IROOT = 1, NROOT(IRREP)
             EE = EIGVAL(IROOT, IRREP) * FACTEV
             EECM = EIGVAL(IROOT, IRREP) * FACTCM
             ET = EIGVAL(IROOT, IRREP) + ECC
             F = OSCSTR(IROOT,IRREP)
             
             IF (IFLAGS(91) .EQ. 0) THEN
               WRITE(6,2031)  IRREP,BGN(IROOT,IRREP),END(IROOT,IRREP),
     &                         NATURE(IROOT,IRREP), EE, EECM, ET
             ELSE
               IF (F .LT. 1.0D-15) F = ZILCH
               WRITE(6,2030)  IRREP, BGN(IROOT,IRREP),
     &                         END(IROOT,IRREP),NATURE(IROOT,IRREP),
     &                         EE,EECM, F, ET
             ENDIF
           ENDDO
         ENDDO
         WRITE(6,*)
         WRITE(6,*)'  ------------------------------------',
     $      '--------------------------------------'
         WRITE(6,*)
C
 2030    FORMAT(1X,I3,2X,I4,4x,I4,9x,A,1X, 
     &          F12.4,F12.2,3X,E9.4,2X,F15.8)
 2031    FORMAT(1X,I3,2X,I4,4x,I4,9x,A,1X,F12.4,F12.2,8X,F15.8)
C
#endif 
      G_TDENS      = (IFLAGS2(h_IFLAGS2_genrlizd_tdens) .Gt. 0) 
      CORE_SEARCH  = (Iflags2(h_IFLAGS2_ee_search) .EQ. 1)
      DOUBLE_EES   = (Iflags2(h_IFLAGS2_ee_search) .EQ. 4)
      EXT_PRINTING = (G_TDENS .OR. CORE_SEARCH .OR. DOUBLE_EES)

      CALL GETREC(20,'JOBARC','TOTENERG',IINTFP,ECC)

      IF (G_TDENS) THEN
         Call GETREC(0,"JOBARC","GEN_OSCS",Length,Dummy) 
         Call GETREC(20,"JOBARC","GEN_OSCS",Length*Iintfp,Osc)
      ENDIF 

      IF (CORE_SEARCH .OR. DOUBLE_EES .OR. G_TDENS) THEN
          Call GETREC(0,"JOBARC","GEN_OSCS",Length,Dummy) 
          Call GETREC(20,"JOBARC","QUADMOMC",Length*Iintfp,Osc_Qm)
          Call GETREC(20,"JOBARC","OCTPMOMC",Length*Iintfp,Osc_Om)
          Call GETREC(20,"JOBARC","ANGLMOMC",Length*Iintfp,Osc_Am)
          Call GETREC(20,"JOBARC","MAGQUADC",Length*Iintfp,Osc_Mq)
          Call GETREC(20,"JOBARC","MAGQUADC",Length*Iintfp,Osc_Dm)
      ENDIF 

      CALL GETREC(20,'JOBARC','TOTENERG',IINTFP,ECC)
      If (EOM_Trps) Then

      Write(6,*)
      Write(6,100)
      Write(6,101)
 100  Format(2x,'Sym',2x,'Origin',4x,'Destination',1x,'EOM-',7x,
     &       'EOM-',8x,'Osc. Str.',4x,'Total Energy')
 101  Format(28x,' CCSD(eV)',3x,'CCSD(T)(eV)')
      Write(6,*)' ----------------------------------------------------',
     &          '-------------------------'
      Write(6,*)
      Write(6,"(2x,a)") "The Singlet Excited States:"
      Write(6,*)
         Do Irrep = 1, Nirrep
            Do Iroot = 1, Nroot(Irrep)
               If (Nature(Iroot,Irrep) .EQ. "S") Then
                  Eed   = Eigval(Iroot,Irrep) * Factev
                  Eet   = Eigval_t(Iroot,Irrep) * Factev
                  Et    = Eigval(Iroot,Irrep) + Ecc
                  F     = Oscstr(Iroot,Irrep)

                  Write(6,99) Irrep,Bgn(iroot,irrep),
     &                       Bgn_irp(iroot,irrep),
     &                       End(iroot,irrep),End_irp(iroot,irrep),
     &                       Eed,Eet,F,Et
               Else
                  Triplet = .True.
               Endif
            Enddo
         Enddo

         If (Triplet) Then
         Write(6,*)
         Write(6,"(2x,a)") "The Triplet Excited States:"
         Write(6,*)
         Do Irrep = 1, Nirrep
            Do Iroot = 1, Nroot(Irrep)
               If (Nature(Iroot,Irrep) .EQ. "T") Then
                  Eed   = Eigval(Iroot,Irrep) * Factev
                  Eet   = Eigval_t(Iroot,Irrep) * Factev
                  Et    = Eigval(Iroot,Irrep) + Ecc
                  F     = Oscstr(Iroot,Irrep)

                  Write(6,99) Irrep,Bgn(iroot,irrep),
     &                       Bgn_irp(iroot,irrep),
     &                       End(iroot,irrep),End_irp(iroot,irrep),
     &                       Eed,Eet,F,Et
               Endif
            Enddo
          Enddo
         Endif

  99   Format(1x,i3,1x,i3,"[",i1,"]",4x,i3,"[",i1,"]",4x,F12.4,1x,F12.4,
     &        3x,E9.4,1x,F15.8)
      Else

      Ndone = 0
      Write(6,*) 
      Write(6,10)
 10   Format(2x,'Sym',2x,'Origin',4x,'Destination',3x,'EE(eV)',5x,
     &       'EE(cm-1)',4x,'Osc. Str.',4x,'Total Energy')
      Write(6,*)' ----------------------------------------------------',
     &          '-------------------------'
      Write(6,*)
      Write(6,"(2x,a)") "The Singlet Excited States:" 
      Write(6,*)
      S = 0
      T = 0
         Do Irrep = 1, Nirrep 
            Do Iroot = 1, Nroot(Irrep)
               If (Nature(Iroot,Irrep) .EQ. "S") Then
                  Ee    = Eigval(Iroot,Irrep) * Factev
                  Eecm  = Eigval(Iroot,Irrep) * Factcm
                  Et    = Eigval(Iroot,Irrep) + Ecc
                  F     = Oscstr(Iroot,Irrep) 
                  S     = S +  1

                  Write(6,9) Irrep,Bgn(iroot,irrep),
     &                       Bgn_irp(iroot,irrep),
     &                       End(iroot,irrep),End_irp(iroot,irrep),
     &                       Ee,Eecm,F,Et 
               Else 
                  Triplet = .True.
               Endif 
            Enddo
         Enddo

          IF (EXT_PRINTING .And. S .Gt. 0) THEN
              Write(6,*)
              Write(6,"(3a)")"  -------------------Oscillator",
     &                       " strengths between states-----",
     &                       "------------------"
              Write(6,"(25x,3a)")"  -(Beypnd Dipole",
     &                       " approximation)-"
              Write(6,*)
              Write(6,19)
              Write(6,20)
              Do Irrep = 1, Nirrep
                 Do Iroot = 1, Nroot(Irrep)
                  Es    = Eigval(Iroot,Irrep) + Ecc
                   Do Jrrep = 1, Nirrep
                      Do Jroot = 1,  Nroot(Jrrep)
                         If (Nature(Iroot,Irrep) .EQ. "S" .and.
     &                       Nature(Jroot,Jrrep) .EQ. "S") Then
                             Ndone = Ndone + 1
                             Ef    = Eigval(Jroot,Jrrep) + Ecc
      Write(6,*)
      Write(6,200) "Electric dipole contribution                   :",
     +              Osc_dm(Ndone)
      Write(6,200) "Magnetic dipole contribution                   :",
     +              OSc_am(Ndone)
      Write(6,200) "Electric quadrupole contribution               :",
     +              Osc_qm(Ndone)
      Write(6,200) "Electric dipole/electric octupole contribution :",
     +              Osc_om(Ndone)
      Write(6,200) "Electric dipole/Magnetic octupole contribution :",
     +              Osc_mq(Ndone)
      Os =  Osc(Ndone)+OSc_am(Ndone)+Osc_qm(Ndone)+ Osc_om(Ndone) + 
     +      Osc_mq(Ndone)
      Write(6,*)
      Write(6,29) Irrep,Iroot,Jrrep,jroot,Es,Ef,Os
                         Endif
                      Enddo
                   Enddo
                 Enddo
              Enddo
          Endif

         If (Triplet) Then
         Write(6,*)
         Write(6,"(2x,a)") "The Triplet Excited States:"
         Write(6,*)
         Do Irrep = 1, Nirrep
            Do Iroot = 1, Nroot(Irrep)
               If (Nature(Iroot,Irrep) .EQ. "T") Then
                  Ee    = Eigval(Iroot,Irrep) * Factev
                  Eecm  = Eigval(Iroot,Irrep) * Factcm
                  Et    = Eigval(Iroot,Irrep) + Ecc
                  F     = Oscstr(Iroot,Irrep)
                  T     = T + 1

                  Write(6,9) Irrep,Bgn(iroot,irrep),
     &                       Bgn_irp(iroot,irrep),
     &                       End(iroot,irrep),End_irp(iroot,irrep),
     &                       Ee,Eecm,F,Et
               Endif
            Enddo
          Enddo
          Endif 

          IF (EXT_PRINTING .And. T .Gt. 0) THEN
              Write(6,*)
              Write(6,"(3a)")"  -------------------Oscillator",
     &                       " strengths between states-----",
     &                       "------------------"
              Write(6,"(25x,3a)")"  -(Beypnd Dipole",
     &                       " approximation)-"
              Write(6,*)
              Write(6,19)
              Do Irrep = 1, Nirrep
                 Do Iroot = 1, Nroot(Irrep)
                  Es    = Eigval(Iroot,Irrep) + Ecc
                   Do Jrrep = 1, Nirrep
                      Do Jroot = 1,  Nroot(Jrrep)
                         If (Nature(Iroot,Irrep) .EQ. "T" .and.
     &                       Nature(Jroot,Jrrep) .EQ. "T") Then
                             Ndone = Ndone + 1
                             Ef    = Eigval(Jroot,Jrrep) + Ecc
      Write(6,*)
      Write(6,200) "Electric dipole contribution                   :",
     +              Osc_dm(Ndone)
      Write(6,200) "Magnetic dipole contribution                   :",
     +              OSc_am(Ndone)
      Write(6,200) "Electric quadrupole contribution               :",
     +              Osc_qm(Ndone)
      Write(6,200) "Electric dipole/electric octupole contribution :",
     +              Osc_om(Ndone)
      Write(6,200) "Electric dipole/Magnetic octupole contribution :",
     +              Osc_mq(Ndone)
      Os =  Osc(Ndone)+OSc_am(Ndone)+Osc_qm(Ndone)+ Osc_om(Ndone) + 
     +      Osc_mq(Ndone)
      Write(6,*)
      Write(6,29) Irrep,Iroot,Jrrep,jroot,Es,Ef,Os

                         Endif 
                      Enddo
                   Enddo 
                 Enddo
              Enddo
          Endif 
      Endif 

      Write(6,*)
      Write(6,*)' ----------------------------------------------------',
     &          '-------------------------'
  9   Format(1x,i3,1x,i3,"[",i1,"]",4x,i3,"[",i1,"]",4x,F12.4,1x,F12.2,
     &       3x,E9.4,1x,F15.8)


 19   Format(2x,'State-1',12x,'State-2',15x,'E1',12x,'E2',12x,
     &      'Osc. Str.')
 20   Format(2x,'Irrep1[root No.]',3x,'Irrep2[root No.]')

 29   Format(1x,i3,"[",i1,"]",14x,i3,"[",i1,"]",7x,F15.8,1x,F15.8,
     &       5x,E9.4)
 1300    FORMAT(/)
 200  Format(1x,a,1x,es12.4e2)
      Write(6,*)
      Write(6,"(a,a)") " Comments! If there are empty rows all that ",
     &                 "means assigining singlet or triplet"
      Write(6,"(a,a)") " characted to a state has failed. Most likely",
     &                 " scenario is that they are not" 
      Write(6,"(a,a)") " pure singlets or triplets. Also, the zero",
     &                 " excitation energy indicates that the"
      Write(6,"(a)")   " root is not converged."

      RETURN
      END

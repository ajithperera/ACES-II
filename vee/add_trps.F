#include <flags.h>
      Subroutine Add_trps(Scr, Maxcor, Irrepx, Root, Iuhf, Iside,
     &                    Iroot)

      Implicit Double Precision (A-H, O-Z)
      Logical T_TILDE
      Character*1 Nature(100,8)
      Integer Bgn,Bgn_irp,End,End_irp

#include <machsp.com>
#include <jodaflags.com>

      COMMON/ROOTS/EIGVAL(100,8),EIGVAL_T(100,8),OSCSTR(100,8),
     &             BGN(100,8),BGN_IRP(100,8),END(100,8),
     &             END_IRP(100,8),NATURE
      
      Dimension Scr(Maxcor)
      DATA AU2EV /27.2113957D0/

C Non iterative models (there are no triples in H_bar*R). The Iside_PT
C controls whether R3 or L3 being built. 

      Iside_PT = 1
      T_TILDE  = .FALSE.
      Rootin   = 0.0D0

      If (Iside .EQ. 2) Call Eom_ccsdpt(Scr, Maxcor, Iuhf, Irrepx, 
     &                                  Rootin, ISide_PT, EL1R3, 
     &                                  EL2R3,
     &                                  EL3R2, EL3R3, OL3R3, EL2TR1)
C
      ECCSD_PT_0             = ROOT + EL1R3 + EL2R3
      EIGVAL_T(Iroot,Irrepx) = ECCSD_PT_0 

      Write(6, "(a)") " Omega = 0.0 is used for The triples correction."

      Write(6, "(a, F20.10,a)") " The triples correction : ", 
     &                            EL1R3 + EL2R3," a.u."

      Write(6, "(a, F20.10,a)") " The EOM-CCSD energy    : ", 
     &                            ROOT*AU2EV, " eV"
      Write(6, "(a, F20.10,a)") " The EOM-CCSD(T) energy : ", 
     &                            ECCSD_PT_0*AU2EV, " eV"
C
      If (T_TILDE) THEN
      Rootin   =  Root
      If (Iside .EQ. 2) Call Eom_ccsdpt(Scr, Maxcor, Iuhf, Irrepx, 
     &                                  Rootin, ISide_PT, EL1R3, 
     &                                  EL2R3,
     &                                  EL3R2, EL3R3, OL3R3, EL2TR1)
C
      ECCSD_PT_0             = ROOT + EL1R3 + EL2R3
      EIGVAL_T(Iroot,Irrepx) = ECCSD_PT_0

      Write(6, "(a,a)") " Omega = EOM-CCSD is used for The triples",
     &                  " correction."

      Write(6, "(a, F20.10,a)") " The triples correction :", 
     &                            EL1R3 + EL2R3," a.u."

      Write(6, "(a, F20.10,a)") " The EOM-CCSD energy    : ", 
     &                            ROOT*AU2EV, " eV"
      Write(6, "(a, F20.10,a)") " The EOM-CCSD(T) energy :", 
     &                            ECCSD_PT_0*AU2EV, " eV"
      Endif
C

C
      Return
      End

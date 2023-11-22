#include "flags.h"

      Subroutine Scrnc_set_vars(Iuhf)

      Implicit Double Precision (A-H,O-Z)
      Logical Ccsd,Mbpt,Parteom,Nodavid,Rccd,Drccd
      Logical Freq_exist,Realfreq
      Logical Ss,Sd,Ds,Dd
      Logical Respn_d
      Logical Aoladder 
      Double Precision XFreq(200)

      Common /Extrap/Maxexp,Nreduce,Ntol,Nsizec
      Common /Eominfo/Ccsd,Mbpt,Parteom,Nodavid,Aoladder,Rccd,
     +                Drccd 
      Common /Freqinfo/Freq,Realfreq,Xfreq,No_freqs 
      Common /Drvhbar/Ss,Sd,Ds,dd
      Common /eomprop/Ieomprop 

#include "flags.com"
#include "flags2.com"
#include "flags2.com"
#include "files.com"

      Ntol    = Iflags(h_IFLAGS_estate_tol) 
      Maxexp  = Iflags(h_IFLAGS_estate_maxc)
      Nreduce = 0

      Ccsd    = .False.
      Parteom = .False.
      Aoladder= .False.
      Rccd    = .False. 
      Drccd   = .False. 

      If (Iflags2(h_IOPPAR_parteom) .Gt.0)  Parteom = .True.
      If (Iflags(h_IOPPAR_calc)    .Eq.10)  Ccsd    = .True.
      If (Iflags(h_IOPPAR_calc)    .Eq.48)  Rccd    = .True.
      If (Iflags(h_IOPPAR_calc)    .Eq.49)  Drccd   = .True.
      If (Iflags(h_IOPPAR_density) .Eq.1)   Respn_d = .True.
      Mbpt  = .Not. (Ccsd .OR. Rccd .OR. Drccd)
      Luout = 6

      Aoladder = (Iflags(h_IFLAGS_abcdtype).eq. 2)
      Ieomprop =  Iflags2(h_IFLAGS2_eomprop)

      If (Iuhf .EQ. 0) Then
         Write(6,*)
         Write(6,"(a,a)") " The approach implemented here is not",
     +                    " supported for RHF reference."
         Write(6,*)
         Call Errex
      Endif 
      If (.NOT. Respn_d) Then
         Write(6,"(a,a)")" The Density=Response is not set in the",
     +                   " ZMAT and calculation can not proceed."
         Write(6,*)
         Call Errex
      Endif 
      Inquire(File="frequency",Exist=Freq_exist)
      If (Freq_Exist) Then
         OPen(Unit=10,File="Frequency",Form="Formatted",Status="OLD")
         Junk = 0
         Read(10,*,End=800) No_Freqs
         Read(10,*,End=800) JunK
         Realfreq = (Junk.eq.0)
         Do Ifreq = 1, No_Freqs
            Read(10,*,End=800) Freq
            Xfreq(Ifreq) = Freq
         Enddo 
 800     Write(6,900)
      Else
         No_freqs = 1
         Freq     = 0.0D0
         RealFreq = .True.  
         Xfreq(1) = Freq
      Endif 

      If (Rccd .or. Drccd) Then
         Nodavid  = .True. 
      Else
         Nodavid  = Parteom .and. Realfreq 
      Endif 

      Ss = .True.
      Sd = .True.
      Ds = .True.
      If (Iflags2(h_IOPPAR_parteom) .Eq.1) Then
         Dd = .False.
      Else 
         Dd = .True.
      Endif 

      Write(6,*)
      If (Nodavid) Then
          If (Rccd .Or. Drccd) Then

          write(6,"(a)") " @-scrnc_set_vars, Linear equations solved",
     +                   " by matrix inversion (AC=Y => C=A^(-1)Y)"
          ELse
          write(6,"(a)") " @-scrnc_set_vars, Linear equations solved",
     +                   " using partitioning."
          Endif 
      Else 
          write(6,*) " @-scrnc_set_vars, Linear equations solved",
     +               " using Davidson."
      Endif 

     
  900 Format(T3,"@-Scrnc_set_vars; ", "file containing listing of",
     +       " frequencies (i.e. frequency)",//, "is incomplete.")
      Return
      End


#include "flags.h"
      Subroutine Check_1body_dens(Doo,Dvv,Dvo,Ioo,Ivv,Ivo,Work,
     &                            Maxcor,Nbas,IUhf) 

      Implicit Integer (A-Z)

#include "machsp.com"
#include "sympop.com"
#include "sym.com"
#include "parcc.com"
#include "maxbasfn.par"
#include "info.com"
#include "flags.com"

      Double Precision Work(Maxcor),Ioo(1),Ivv(1),Ivo(1)
      Double Precision Doo(1),Dvv(1),Dvo(1)
      Double Precision Eoo_aa,Eoo_bb,Evv_aa,Evv_bb,Evo_aa,Evo_bb
      Double Precision scr1(Maxbasfn),Ddot 
      Integer Iscr2(Maxbasfn)
      Logical UHF
C
C This is a routine written for debugging purposes. What this
C is to add delta to one electron terms (fock matrices etc.).

      UHF    = .False.
      UHF    = (Iuhf .EQ. 1)
      Irrepx = 1

#ifdef _DFULL
      I000 = 1
      I010 = I000 + Nbas*Nbas
      I020 = I010 + Nbas*Nbas
      I030 = I020 + Nbas*Nbas
      I040 = I030 + Nbas*Nbas
      If (I040  .Gt. Maxcor) 
     &   Call Insmem("@-Check_1body_dens",I040,Maxcor)

      Call Getrec(20,"JOBARC","ONEINTAO",Nbas*Nbas,Work(I000))
      Call Getrec(20,"JOABRC","RELDENSA",Nbas*Nbas,Work(I010))
      Call mo2ao2(Work(I010),Work(I020),Work(I030),Work(I040),
     &            Nbas, 1)
      Eoo_aa = Ddot(Nbas*Nbas, Work(I000), 1, Work(I020),1)
      Write(6,*)
      Write(6,"(a,F15.9)") "The AA Dpq*hpq energy = ", Eoo_aa
      If (UHF) Then
         I010 = I000 + Nbas*Nbas
         I020 = I010 + Nbas*Nbas
         I030 = I020 + Nbas*Nbas
         I040 = I030 + Nbas*Nbas
         If (I040  .Gt. Maxcor)
     &       Call Insmem("@-Check_1body_dens",I040,Maxcor)

         Call Getrec(20,"JOABRC","RELDENSA",Nbas*Nbas,Work(I010))
         Call mo2ao2(Work(I010),Work(I020),Work(I030),Work(I040),
     &                   Nbas, 1)
         Eoo_bb = Ddot(Nbas*Nbas, Work(I000), 1, Work(I020),1)
         Write(6,"(a,F15.9)") "The BB Dpq*hpq energy = ", Eoo_aa
      Endif 
#endif 

C The DOO check 
#ifdef _OOD_BLOCK
      I000 = 1
      Call Getrec(20,"JOBARC","ONEINTSA",Nbas*Nbas,Work(I000))
      Lenooa      = Nocco(1)*Nocco(1)
      Lenoob      = Nocco(2)*Nocco(2)
      I010 = I000 + Nbas*Nbas
      I020 = I010 + Nbas*Nbas
      I030 = I020 + Nbas*Nbas
      I040 = I030 + Nbas*Nbas
      If (I040  .Gt. Maxcor) 
     &   Call Insmem("@-Check_1body_dens",I040,Maxcor)
      
      Call Dzero(Work(I010),Nbas*Nbas)
      Call Putblk(Work(I010),Doo(1),"OO",scr1,Iscr2,Nbas,
     &            Irpdpd(1,21),Lenooa,1)
cSSS      Write(6,*) "The occ-occ block in MO basis"
CsSS      call output(Work(I010),1,nocco(1),1,nocco(1),nbas,
CSSS     &            nbas,1)
      Call mo2ao2(Work(I010),Work(I020),Work(I030),Work(I040),
     &            Nbas, 1)
CSSS      Write(6,*) "The occ-occ block in AO basis"
CSsS      call output(Work(I020),1,nbas,1,nbas,nbas,nbas,1)
C
      Eoo_aa = Ddot(Nbas*Nbas, Work(I000), 1, Work(I020),1)
      Write(6,*)
      Write(6,"(a,F15.9)") "The AA DOO*hOO energy = ", Eoo_aa

      IF (UHF) Then
         Call Getrec(20,"JOBARC","ONEINTSB",Nbas*Nbas,Work(I000))
         I010 = I000 + Nbas*Nbas
         I020 = I010 + Nbas*Nbas 
         I030 = I020 + Nbas*Nbas
         I040 = I030 + Nbas*Nbas
         If (I040  .Gt. Maxcor) 
     &      Call Insmem("@-Check_1body_dens",I040,Maxcor)

         Call Dzero(Work(I010),Nbas*Nbas)
         Call Putblk(Work(I010),Doo(nfmi(1)+1),"OO",scr1,Iscr2,Nbas,
     &               Irpdpd(1,22),Lenoob,2)
CSSS      Write(6,*) "The occ-occ block in MO basis"
CSSS      call output(Work(I010),1,nocco(2),1,nocco(2),nbas,
CssS     &            nbas,1)

         Call mo2ao2(Work(I010),Work(I020),Work(I030),Work(I040),
     &               Nbas, 2)
CSSS      Write(6,*) "The occ-occ block in AO basis"
CSSS      call output(Work(I020),1,nbas,1,nbas,nbas,nbas,1)

         Eoo_bb = Ddot(Nbas*Nbas, Work(I000), 1, Work(I020),1)
         Write(6,"(a,F15.9)") "The BB DOO*hOO energy = ", Eoo_bb
      Endif 
#endif 

C The DVV check 
#ifdef _VVD_BLOCK

      I000 = 1
      Call Getrec(20,"JOBARC","ONEINTSA",Nbas*Nbas,Work(I000))
      Lenvva      = Nvrto(1)*Nvrto(1)
      Lenvvb      = Nvrto(2)*Nvrto(2)
      I010 = I000 + Nbas*Nbas 
      I020 = I010 + Nbas*Nbas 
      I030 = I020 + Nbas*Nbas
      I040 = I030 + Nbas*Nbas
      If (I040  .Gt. Maxcor) 
     &   Call Insmem("@-Check_1body_dens",I040,Maxcor)

      Call Dzero(Work(I010),Nbas*Nbas)
      Call Putblk(Work(I010),Dvv(1),"VV",scr1,Iscr2,Nbas,
     &            Irpdpd(1,19),Lenvva,1)
CSSS      Write(6,*) "The vrt-vrt block in MO basis"
CSSS      call output(Work(I010),1,nvrto(1),1,nvrto(1),nbas,nbas,1)

      Call mo2ao2(Work(I010),Work(I020),Work(I030),Work(I040),
     &            Nbas, 1)
CSSS      Write(6,*) "The vrt-vrt block in AO basis"
CSSS      call output(Work(I020),1,nbas,1,nbas,nbas,nbas,1)
      Evv_bb = Ddot(Nbas*Nbas, Work(I000), 1, Work(I020),1)
      Write(6,"(a,F15.9)") "The AA DVV*hVV energy = ", Evv_bb

      IF (UHF) Then
         Call Getrec(20,"JOBARC","ONEINTSB",Nbas*Nbas,Work(I000))
         I010 = I000 + Nbas*Nbas
         I020 = I010 + Nbas*Nbas 
	 I030 = I020 + Nbas*Nbas
	 I040 = I030 + Nbas*Nbas
	 If (I040  .Gt. Maxcor) 
     &      Call Insmem("@-Check_1body_dens",I040,Maxcor)

         Call Dzero(Work(I010),Nbas*Nbas)
         Call Putblk(Work(I010),Dvv(nfea(1)+1),"VV",scr1,Iscr2,Nbas,
     &               Irpdpd(1,20),Lenvvb,2)
         Call mo2ao2(Work(I010),Work(I020),Work(I030),Work(I040),
     &               Nbas, 2)
         Evv_bb = Ddot(Nbas*Nbas, Work(I000), 1, Work(I020),1)
         Write(6,"(a,F15.9)") "The BB Dvv*hVV energy = ", Evv_bb
      Endif 
#endif 

C The DVO check 
#ifdef _OVD_BLOCK

      I000 = 1
      Call Getrec(20,"JOBARC","ONEINTSA",Nbas*Nbas,Work(I000))
      If (iflags(h_IFLAGS_nonhf) .GT. 0) Then
      Lenvoa      = Nocco(1)*Nvrto(1)
      LenvOb      = Nocco(2)*Nvrto(2)
      I010 = I000 + Nbas*Nbas
      I020 = I010 + Nbas*Nbas
      I030 = I020 + Nbas*Nbas
      I040 = I030 + Nbas*Nbas
      If (I040  .Gt. Maxcor) 
     &   Call Insmem("@-Check_1body_dens",I040,Maxcor)

      Call Dzero(Work(I010),Nbas*Nbas)
      Call Putblk(Work(I010),Dvo(1),"VO",scr1,Iscr2,Nbas,
     &            Irpdpd(1,9),Lenvoa,1)
      Call mo2ao2(Work(I010),Work(I020),Work(I030),Work(I040),
     &            Nbas, 1)
      Evo_aa = Ddot(Nbas*Nbas, Work(I000), 1, Work(I020),1)
      Write(6,"(a,F15.9)") "The AA DVO*hVO energy = ", Evo_aa

      If (UHF) Then
         Call Getrec(20,"JOBARC","ONEINTSB",Nbas*Nbas,Work(I000))
         I010 = I000 + Nbas*Nbas
         I020 = I010 + Nbas*Nbas
         I030 = I020 + Nbas*Nbas
         I040 = I030 + Nbas*Nbas
         If (I040  .Gt. Maxcor) 
     &   Call Insmem("@-Check_1body_dens",I040,Maxcor)

         Call Dzero(Work(I010),Nbas*Nbas)
         Call Putblk(Work(I010),Dvo(nt(1)+1),"VO",scr1,Iscr2,Nbas,
     &            Irpdpd(1,10),Lenvoa,1)
         Call mo2ao2(Work(I010),Work(I020),Work(I030),Work(I040),
     &               Nbas, 2)
         Evo_bb = Ddot(Nbas*Nbas, Work(I000), 1, Work(I020),1)
         Write(6,"(a,F15.9)") "The AA DVo*hVO energy = ", Evo_bb

      Endif 
      Endif 
#endif 

C The IOO check 
#ifdef _OOI_BLOCK
      I000 = 1
      Call Getrec(20,"JOBARC","ONEINTSA",Nbas*Nbas,Work(I000))
      I010 = I000 + Nbas*Nbas 
      I020 = I010 + Nbas*Nbas 
      I030 = I020 + Nbas*Nbas
      I040 = I030 + Nbas*Nbas
      If (I040  .Gt. Maxcor) 
     &   Call Insmem("@-Check_1body_dens",I040,Maxcor)

      Call Dzero(Work(I010),Nbas*Nbas)
      Call Putblk(Work(I010),Ioo(1),"OO",scr1,Iscr2,Nbas,
     &            Irpdpd(1,21),Lenooa,1)

CSSS      Write(6,*) "The occ-occ block in MO basis"
CSSS      call output(Work(I010),1,nocco(2),1,nocco(2),nbas,
CSSS     &            nbas,1)

      Call mo2ao2(Work(I010),Work(I020),Work(I030),Work(I040),
     &            Nbas, 1)
      Eoo_aa = Ddot(Nbas*Nbas, Work(I000), 1, Work(I020),1)
      Write(6,"(a,F15.9)") "The AA IOO*hOO energy = ", Eoo_aa

      IF (UHF) Then
         Call Getrec(20,"JOBARC","ONEINTSB",Nbas*Nbas,Work(I000))
         I010 = I000 + Nbas*Nbas 
         I020 = I010 + Nbas*Nbas
         I030 = I020 + Nbas*Nbas
         I040 = I030 + Nbas*Nbas
         If (I040  .Gt. Maxcor) 
     &   Call Insmem("@-Check_1body_dens",I040,Maxcor)

         Call Dzero(Work(I010),Nbas*Nbas)
         Call Putblk(Work(I010),Ioo(nfmi(1)+1),"OO",scr1,Iscr2,Nbas,
     &               Irpdpd(1,22),Lenoob,2)

CSSS      Write(6,*) "The occ-occ block in MO basis"
CSSS      call output(Work(I010),1,nocco(2),1,nocco(2),nbas,
CSSS     &            nbas,1)

         Call mo2ao2(Work(I010),Work(I020),Work(I030),Work(I040),
     &               Nbas, 2)
         Eoo_bb = Ddot(Nbas*Nbas, Work(I000), 1, Work(I020),1)
         Write(6,"(a,F15.9)") "The BB IOO*hOO energy = ", Eoo_bb
      Endif 
#endif 

C The IVV check 
#ifdef _VVI_BLOCK
      I000 = 1
      Call Getrec(20,"JOBARC","ONEINTSA",Nbas*Nbas,Work(I000))
      I010 = I000 + Nbas*Nbas 
      I020 = I010 + Nbas*Nbas
      I030 = I020 + Nbas*Nbas
      I040 = I030 + Nbas*Nbas
      If (I040  .Gt. Maxcor) 
     &   Call Insmem("@-Check_1body_dens",I040,Maxcor)

      Call Dzero(Work(I010),Nbas*Nbas)
      Call Putblk(Work(I010),Ivv(1),"VV",scr1,Iscr2,Nbas,
     &            Irpdpd(1,19),Lenvva,1)
      Call mo2ao2(Work(I010),Work(I020),Work(I030),Work(I040),
     &            Nbas, 1)
      Evv_bb = Ddot(Nbas*Nbas, Work(I000), 1, Work(I020),1)
      Write(6,"(a,F15.8)") "The AA IVV*hVV energy = ", Evv_bb

      IF (UHF) Then
         Call Getrec(20,"JOBARC","ONEINTSA",Nbas*Nbas,Work(I000))
         I010 = I000 + Nbas*Nbas 
         I020 = I010 + Nbas*Nbas
         I030 = I020 + Nbas*Nbas
         I040 = I030 + Nbas*Nbas
         If (I040  .Gt. Maxcor) 
     &   Call Insmem("@-Check_1body_dens",I040,Maxcor)

         Call Dzero(Work(I010),Nbas*Nbas)
         Call Putblk(Work(I010),Ivv(nfea(1)+1),"VV",scr1,Iscr2,Nbas,
     &               Irpdpd(1,20),Lenvvb,2)
         Call mo2ao2(Work(I010),Work(I020),Work(I030),Work(I040),
     &               Nbas, 22)
         Evv_bb = Ddot(Nbas*Nbas, Work(I000), 1, Work(I020),1)
         Write(6,"(a,F15.8)") "The BB Ivv*hVV energy = ", Evv_bb
      Endif 
#endif 

C The IVO check 
#ifdef _OVI_BLOCK
      If (iflags(h_IFLAGS_nonhf) .GT. 0) Then
      I000 = 1
      Call Getrec(20,"JOBARC","ONEINTSA",Nbas*Nbas,Work(I000))
      I010 = I000 + Nbas*Nbas
      I020 = I010 + Nbas*Nbas
      I030 = I020 + Nbas*Nbas
      I040 = I030 + Nbas*Nbas
      If (I040  .Gt. Maxcor) 
     &   Call Insmem("@-Check_1body_dens",I040,Maxcor)

      Call Dzero(Work(I010),Nbas*Nbas)
      Call Putblk(Work(I010),Ivo(1),"VO",scr1,Iscr2,Nbas,
     &            Irpdpd(1,9),Lenvoa,1)
      Call mo2ao2(Work(I010),Work(I020),Work(I030),Work(I040),
     &            Nbas, 1)
      Evo_aa = Ddot(Nbas*Nbas, Work(I000), 1, Work(I020),1)
      Write(6,"(a,F15.8)") "The AA IVO*hVO energy = ", Evo_aa

      If (UHF) Then
         Call Getrec(20,"JOBARC","ONEINTSB",Nbas*Nbas,Work(I000))
         I010 = I000 + Nbas*Nbas
         I020 = I010 + Nbas*Nbas
         I030 = I020 + Nbas*Nbas
         I040 = I030 + Nbas*Nbas
         If (I040  .Gt. Maxcor) 
     &   Call Insmem("@-Check_1body_dens",I040,Maxcor)

         Call Dzero(Work(I010),Nbas*Nbas)
         Call Putblk(Work(I010),Ivo(nt(1)+1),"VO",scr1,Iscr2,Nbas,
     &            Irpdpd(1,10),Lenvoa,1)
         Call mo2ao2(Work(I010),Work(I020),Work(I030),Work(I040),
     &               Nbas, 1)
         Evo_bb = Ddot(Nbas*Nbas, Work(I000), 1, Work(I020),1)
         Write(6,"(a,F15.8)") "The BB IVO*hVO energy = ", Evo_bb

      Endif 
      Endif 
#endif
      Return
      End 


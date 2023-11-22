#include "flags.h"
      Subroutine Hbar_mult_l(Work, Length, Iuhf)

      Implicit Double Precision (A-H, O-Z)

      Dimension Work(Length)

#include "machsp.com"
#include "syminf.com"
#include "sympop.com"
#include "flags.com"
#include "parcc.com"
C
      Irrepx = 1
C
C All the comments refers to Hbar elements (unless otherwise
C stated) not to CC intermediates.
C
C Do L2xF(ae), L2xF(mi), L1xF(ae) and L1xF(mi) terms to L1 and L2.
C Noticed that for pCC we need scalled F(ae) and F(mi) intermediates
c for the doubles term. 
C The lists 91,92 (1,2) contains the unscalled F elements and
C the 91,92(10,11) contains the scalled elements. The F in L2 
c requires the scalled and F in L1 requires the unscalled. Notice
C that this routine is called after Hbar formation. So, the 
C F(m,e) T1 contribution is allready added to the Habr (unlike
c the lambda code). In vee, these correspond to DFT1INT1 (to singles)
C DF2T2INT2 (doubles) from Hbarxc.
C
      Call F1inl2_leom(Work,Length,Iuhf)
      Call F2inl2_leom(Work,Length,Iuhf)
C
#ifdef _DEBUG_LVL0
      Write(6,*) "The lambda residuals after F"
      Call check_leom(Work,length,Iuhf)
#endif
C
C This is -1/2 L(af,mn)T(mn,ef) or 1/2 T(mn,ef)L(ef,in). For parametrized
C CC this need to be scalled by 1/2 when we do these contributions to the
C doubles equation. In vee, GINC2L. 
C
CSSS#ifdef _DCC_FLAG
      If (Ispar) Then
         Call formg_leom(Work,Length,Iuhf,100,Gae_scale,Gmi_scale)
      Else
CSSS#else
         Call formg_leom(Work,Length,Iuhf,100,1.0D0,1.0D0)
      Endif 
CSSS#endif 
C 
C Do the G(be)<ij||ae> and G(mi)<im||ab> terms to L2.
C
      Call G1inl2_leom(Work,Length,Iuhf)
      Call G2inl2_leom(Work,Length,Iuhf)
C
#ifdef _DEBUG_LVL0
      Write(6,*) "The lambda residuals after G"
      Call check_leom(Work,length,Iuhf)
#endif
C
      If (iflags(h_IFLAGS_calclevel) .gt. 9) Then
C
C L1xW(je,bm) and L1xW(ij,mb) terms to l2
C
          Call L1inl2a_leom(Work,Length,Iuhf,.False.)
          Call L1inl2b_leom(Work,Length,Iuhf,.False.)
C
#ifdef _DEBUG_LVL0
      Write(6,*) "The lambda residuals after l1inl2"
      Call check_leom(Work,length,Iuhf)
#endif 
         If (Ispar) Call Restore_cc_wmbej(Work,Length,Iuhf)
         Call L1inl1_leom(Work,Length,Iuhf,1)
         If (Iuhf .NE. 0) Call L1inl1_leom(Work,Length,Iuhf,2)
         If (Ispar) Call Restore_pdcc_wmbej(Work,Length,Iuhf) 

#ifdef _DEBUG_LVL0
      Write(6,*) "The lambda residuals after l1inl1"
      Call check_leom(Work,length,Iuhf)
#endif
C
C G(ef)xW(ei,fa) and G(mn)xW(mi,na). This is the G intermediate
C contributions to the L1 equation. For this we need unscalled 
C G intermediate and lets redo them unscaled. In vee this is 
C done in GINC1L. 
C
CSSS#ifdef _DCC_FLAG
         If (Ispar) Call formg_leom(Work,Length,Iuhf,100,1.0D0,1.0D0)
CSSS#endif
         Call Ginc1l_leom(Irrepx,Work,Length,Iuhf)

#ifdef _DEBUG_LVL0
      Write(6,*) "The lambda residuals after Ginl1"
      Call check_leom(Work,length,Iuhf)
      Call checkintms(Work,length,Iuhf)
#endif
      Endif
C
C The following there calls does the ladder and ring terms. The first
C one does the HH ladder. Here we need the W(mn,ij) PCc which is formed
C by vcc (post_vcc_mods) and stored in 251-253 list. The W(mn,ij) is
C stored in list 51-53. In vee this DLADDER with the 251-253 lists. 
C The PP ladder we do in pieces. The L2LAD does the <ab||cd> integral
C contribution. The W1LAD terms does the 1/2 L2(ef,j)Tau(mn,ef)<ef|ij> 
C In lambda code this is done ising the V(mn,ij) intermediate. This
C requires T2=0 in the Tau formation. 
C
C Do L2 x W(ij,mn) + L2 x W(efab) terms
C
      If (Ispar) Then
         Call Pdcc_formv1(Work,Length,Iuhf)
      Else
         Call Formv1(Work,Length,Iuhf)
      Endif 
 
      If (Ispar) Then
         Call L2lad_leom(Work,Length,Iuhf,1,250)
         If (Iflags(93) .Eq. 2) Then
            Call Draolad(Work,Length,Iuhf,.True.,1,0,143,60,243,260)
         Else
            Call L2lad_leom(Work,Length,Iuhf,6,250)
         Endif 
      Else 
         Call L2lad_leom(Work,Length,Iuhf,1,50)
         If (Iflags(93) .Eq. 2) Then
            Call Draolad(Work,Length,Iuhf,.True.,1,0,143,60,243,260)
         Else
            Call L2lad_leom(Work,Length,Iuhf,6,50)
         Endif 
      Endif 
#ifdef _DEBUG_LVL0
      Write(6,*) "The lambda residuals after l2lads"
      Call check_leom(Work,length,Iuhf)
#endif

CSSS#ifdef _DCC_FLAG
      If (Ispar) Then
         Call W1lad_leom(Work,Length,Iuhf)
      Endif 
CSSS#endif 
      If (Iuhf .Eq. 0) Then
          Ibgn = 3
      ELse
          Ibgn = 1
      Endif 
 
      Do Ispin = Ibgn, 3
         Call l2rng_leom(Work,Length,Ispin,Iuhf)
      Enddo 

#ifdef _DEBUG_LVL0
      Write(6,*) "The lambda residuals after l2lads/rings"
      Call check_leom(Work,length,Iuhf)
#endif

      If (iflags(h_IFLAGS_calclevel) .gt. 9) Then
C
C This is the remaining piece of the HBAR (Hbar(abef) -t(e,m)<mf||ab>)
C In vee this in done in L1W1.
C
          Call modaibc(Work,length,Iuhf,-1.0D0)
          Call l1w1inl2_leom(Work,Length,Iuhf)
          Call modaibc(Work,length,Iuhf, 1.0D0)
C
C L(ef,im)xW(ie,am) and L(ae,mn)xW(ie,mn). These are Habr(abci)
C and Hbar(ijka) contrubtions to L1. No changes to these Hbar elements
C for pCC and in vee these are done in DT2INT1A and B. 
C
          Call L2inl1b_leom(Work,Length,Iuhf)
          Call L2inl1a_leom(Work,Length,Iuhf)

       Endif 

#ifdef _DEBUG_LVL0
      Write(6,*) "The lambda residuals after l2inl1"
      Call check_leom(Work,length,Iuhf)
#endif

      Return
      End

      

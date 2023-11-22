      Subroutine Scale_h4(Work,Maxcor,IUhf,Scale_yesno,Form_dlists)
   
      Implicit Integer (A-Z)

      Double Precision Work(Maxcor), Delta, Scale_factor,Scale 
      Double Precision Buf(600)
      Dimension Ibuf(600)
      Logical UHF,Scale_yesno,Form_dlists
      Character*80 FNAME

#include "syminf.com"
#include "machsp.com"
#include "sympop.com"
#include "parcc.com"
#include "sym.com"
#include "info.com"
#include "parcc.com"
C
C This is a routine written for debugging purposes. What this 
C is to add delat to type of MO two electon integrals.

      UHF    = .False.
      UHF    = (Iuhf .EQ. 1)
      Irrepx = 1 
      Imode  = 0 
      Write(6,*) 
      Return 

      If (Form_dlists) Then
C Lets scale all the <pq||rs> integrals 
C
C Remove the dependency of t1int2, dwmbej, and t1int1 on integral
C list 17 and 18. These terms also dependent on 23,24,25 and 26
C but those do not interfer with debugging G(IJ,AB) in the denisty
C code. Dependecy on 17 and 18 cause G(IA,JB) to intefere with
C G(IJ,AB). Make copies of 17 and 18 and store them in 121 and
C 122. Also, make duplicates copies of 14,15 and 16 in 147,148
C and 149
     
      Write(6,*) 
      Write(6,"(a)") "-----Warning! Only do this for debugging-----" 
      Write(6,"(a,a)") " Duplicate copies of W list 17 and 18 are",
     +               " formed on 121 and 122"
      Write(6,"(a,a)") " Duplicate copies of W list 14-16 are",
     +               " formed on 147-149"
      Write(6,*)

      If (Iuhf .Ne. 0) Then
         Call Inipck(1,10,9,121,Imode,0,1)
         Length_121 = IDSYMSZ(IRREPX,ISYTYP(1,17),ISYTYP(2,17))
         Call Getall(Work, Length_121, Irrepx, 17)
         Call Putall(Work, Length_121, Irrepx, 121)

         Call Inipck(1,12,12, 148,Imode,0,1)
         Length_148 = IDSYMSZ(IRREPX,ISYTYP(1,26),ISYTYP(2,26))
         Call Getall(Work, Length_147, Irrepx, 26)
         Call Putall(Work, Length_147, Irrepx, 148)
      Endif 

      Call Inipck(1,10,10,147,Imode,0,1)
      Length_147 = IDSYMSZ(IRREPX,ISYTYP(1,25),ISYTYP(2,25))
      Call Getall(Work, Length_147, Irrepx, 25)
      Call Putall(Work, Length_147, Irrepx, 147)

      Call Inipck(1,9,10,122,Imode,0,1)
      Length_122 = IDSYMSZ(IRREPX,ISYTYP(1,18),ISYTYP(2,18))
      Call Getall(Work, Length_122, Irrepx, 18)
      Call Putall(Work, Length_122, Irrepx, 122)

      Call Inipck(1,11,12,123,Imode,0,1)
      Length_123 = IDSYMSZ(IRREPX,ISYTYP(1,21),ISYTYP(2,21))
      Call Getall(Work, Length_123, Irrepx, 21)
      Call Putall(Work, Length_123, Irrepx, 123)
      Return 

      Endif 

      IF (Scale_yesno) Then

      Delta        = 0.001D0
      Scale        = 0.0D0
      Scale_factor = 1.0D0 + Delta*Scale

C#ifdef _SCALE_IJAB 

      Write(6,"(a,F7.5)") "<pp||hh> integrals are incremented by ", 
     &                    delta*scale 

      If (Uhf) Then 

C Caution: The scalling of <Ab|Ij> integrals also scales the
C integrals <IA|jb>  and so on. Therfore the numerical gradients 
C correpond to sum of the G(Ij|Ab)+G(Ib,Ja).
C In UHF to test <Ij|Ab> Gamma scalles the lists 16,17,18,21,22
C To test <IJ|AB> scale 19,20 (14 and 15 may not be necessary)
C Also 17 and 18 should not be used in G(IA,JB) contributions.

         Length_14 = IDSYMSZ(IRREPX,ISYTYP(1,14),ISYTYP(2,14))
         Length_15 = IDSYMSZ(IRREPX,ISYTYP(1,15),ISYTYP(2,15))
         Length_16 = IDSYMSZ(IRREPX,ISYTYP(1,16),ISYTYP(2,16))
         Length_17 = IDSYMSZ(IRREPX,ISYTYP(1,17),ISYTYP(2,17))
         Length_18 = IDSYMSZ(IRREPX,ISYTYP(1,18),ISYTYP(2,18))
         Length_19 = IDSYMSZ(IRREPX,ISYTYP(1,19),ISYTYP(2,19))
         Length_20 = IDSYMSZ(IRREPX,ISYTYP(1,20),ISYTYP(2,20))
         Length_21 = IDSYMSZ(IRREPX,ISYTYP(1,21),ISYTYP(2,21))
         Length_22 = IDSYMSZ(IRREPX,ISYTYP(1,22),ISYTYP(2,22))
         Length_23 = IDSYMSZ(IRREPX,ISYTYP(1,23),ISYTYP(2,23))
         Length_24 = IDSYMSZ(IRREPX,ISYTYP(1,24),ISYTYP(2,24))
         Length_25 = IDSYMSZ(IRREPX,ISYTYP(1,25),ISYTYP(2,25))
         Length_26 = IDSYMSZ(IRREPX,ISYTYP(1,26),ISYTYP(2,26))

#ifdef _NOSKIP
         Call Getall(Work, Length_14, Irrepx, 14)
         Call Dscal(Length_14, Scale_factor, Work, 1)
CSSS         Call Dscal(Length_14, Zlich, Work, 1)
         Call Putall(Work, Length_14, Irrepx, 14)
#endif 

#ifdef _NOSKIP
         Call Getall(Work, Length_15, Irrepx, 15)
         Call Dscal(Length_15, Scale_factor, Work, 1)
CSS         Call Dscal(Length_15, Zlich, Work, 1)
        Call Putall(Work, Length_15, Irrepx, 15)
#endif 

C#ifdef _NOSKIP
         Call Getall(Work, Length_16, Irrepx, 16)
         Call Dscal(Length_16, Scale_factor, Work, 1)
         Call Putall(Work, Length_16, Irrepx, 16)
C#endif 

C#ifdef _NOSKIP 
         Call Getall(Work, Length_17, Irrepx, 17)
         Call Dscal(Length_17, Scale_factor, Work, 1)
         Call Putall(Work, Length_17, Irrepx, 17)
C#endif
C#ifdef _NOSKIP
         Call Getall(Work, Length_18, Irrepx, 18)
         Call Dscal(Length_18, Scale_factor, Work, 1)
         Call Putall(Work, Length_18, Irrepx, 18)
C#endif 

#ifdef _NOSKIP
         Call Getall(Work, Length_19, Irrepx, 19)
         Call Dscal(Length_19, Scale_factor, Work, 1)
CSSS         Call Dscal(Length_19, Zlich, Work, 1)
         Call Putall(Work, Length_19, Irrepx, 19)
#endif 

#ifdef _NOSKIP
         Call Getall(Work, Length_20, Irrepx, 20)
         Call Dscal(Length_20, Scale_factor, Work, 1)
CSSS         Call Dscal(Length_20, Zlich, Work, 1)
         Call Putall(Work, Length_20, Irrepx, 20)
#endif 

C#ifdef _NOSKIP 
         Call Getall(Work, Length_21, Irrepx, 21)
         Call Dscal(Length_21, Scale_factor, Work, 1)
CSSS         Call Dscal(Length_21, Zlich, Work, 1)
         Call Putall(Work, Length_21, Irrepx, 21)
C#endif 

C#ifdef _NOSKIP 
         Call Getall(Work, Length_22, Irrepx, 22)
         Call Dscal(Length_22, Scale_factor, Work, 1)
CSSS         Call Dscal(Length_21, Zlich, Work, 1)
         Call Putall(Work, Length_22, Irrepx, 22)
C#endif
         If (Ispar .and. Coulomb) Then

            Write(6,"(a,a,F7.5)")"<pp|hh> integrals for DCC are also",
     &                      " incremented by ", Delta*scale
            Length_119 = IDSYMSZ(IRREPX,ISYTYP(1,119),ISYTYP(2,119))
            Length_120 = IDSYMSZ(IRREPX,ISYTYP(1,120),ISYTYP(2,120))

#ifdef _NOSKIP
            Call Getall(Work, Length_119, Irrepx, 119)
            Call Dscal(Length_119, Scale_factor, Work, 1)
            Call Putall(Work, Length_119, Irrepx, 119)
#endif

#ifdef _NOSKIP
           Call Getall(Work, Length_120, Irrepx, 120)
           Call Dscal(Length_120, Scale_factor, Work, 1)
           Call Putall(Work, Length_120, Irrepx, 120)
#endif
         Endif

      Else

C For RHF Gamma(Ij|Ab)

         Length_14 = IDSYMSZ(IRREPX,ISYTYP(1,14),ISYTYP(2,14))
         Length_16 = IDSYMSZ(IRREPX,ISYTYP(1,16),ISYTYP(2,16))
         Length_18 = IDSYMSZ(IRREPX,ISYTYP(1,18),ISYTYP(2,18))
         Length_19 = IDSYMSZ(IRREPX,ISYTYP(1,19),ISYTYP(2,19))
         Length_21 = IDSYMSZ(IRREPX,ISYTYP(1,21),ISYTYP(2,21))

C For RHF, scale the lists 14,16,18,19 and 21. Make sure that
C 18 is not used (IA,JB) contributions 

#ifdef _NOSKIP
         Call Getall(Work, Length_14, Irrepx, 14)
CSSS         Write(6, "(6(1x,F10.6))") (Work(i),i=1,Length_14)
         Call Dscal(Length_14, Scale_factor, Work, 1)
         Call Putall(Work, Length_14, Irrepx, 14)
#endif 

C#ifdef _NOSKIP
         Call Getall(Work, Length_16, Irrepx, 16)
         Call Dscal(Length_16, Scale_factor, Work, 1)
         Call Putall(Work, Length_16, Irrepx, 16)
C#endif 

#ifdef _NOSKIP
         Call Getall(Work, Length_18, Irrepx, 18)
         Call Dscal(Length_18, Scale_factor, Work, 1)
         Call Putall(Work, Length_18, Irrepx, 18)
#endif 

#ifdef _NOSKIP
         Call Getall(Work, Length_19, Irrepx, 19)
         Call Dscal(Length_19, Scale_factor, Work, 1)
         Call Putall(Work, Length_19, Irrepx, 19)
#endif 

C#ifdef _NOSKIP
         Call Getall(Work, Length_21, Irrepx, 21)
         Call Dscal(Length_21, Scale_factor, Work, 1)
         Call Putall(Work, Length_21, Irrepx, 21)
C#endif 
         If (Ispar .and. Coulomb) Then
            Write(6,"(a,a,F7.5)")"<pp|hh> integrals for DCC are also",
     &                      " incremented by ", Delta*scale
             Length_119 = IDSYMSZ(IRREPX,ISYTYP(1,119),ISYTYP(2,119))

C#ifdef _NOSKIP
             Call Getall(Work, Length_119, Irrepx, 119)
             Call Dscal(Length_119, Scale_factor, Work, 1)
             Call Putall(Work, Length_119, Irrepx, 119)
C#endif
         Endif
      Endif

      Else 

      Delta        = 0.001D0
      Scale        = 0.0D0
      Scale_factor = 1.0D0 + Delta*Scale
      Scale_factor = 1.0D0/Scale_factor

      UHF    = .False.
      UHF    = (Iuhf .EQ. 1)
      Irrepx = 1
      Imode  = 0
      Write(6,*)

C#ifdef _SCALE_IJAB 

      Write(6,"(a,a)") "Restoring scalled 14-22 lists to orignal",
     &                 " values"
      Write(6,"(a,F7.5)") "<pp||hh> integrals are scalled by ", 
     &                     Scale_factor

      If (Uhf) Then 

C Caution: The scalling of <Ab|Ij> integrals also scales the
C integrals <IA|jb>  and so on. Therfore the numerical gradients 
C correpond to sum of the G(Ij|Ab)+G(Ib,Ja).
C In UHF to test <Ij|Ab> Gamma scalles the lists 16,17,18,21,22
C To test <IJ|AB> scale 19,20 (14 and 15 may not be necessary)
C Also 17 and 18 should not be used in G(IA,JB) contributions.

         Length_14 = IDSYMSZ(IRREPX,ISYTYP(1,14),ISYTYP(2,14))
         Length_15 = IDSYMSZ(IRREPX,ISYTYP(1,15),ISYTYP(2,15))
         Length_16 = IDSYMSZ(IRREPX,ISYTYP(1,16),ISYTYP(2,16))
         Length_17 = IDSYMSZ(IRREPX,ISYTYP(1,17),ISYTYP(2,17))
         Length_18 = IDSYMSZ(IRREPX,ISYTYP(1,18),ISYTYP(2,18))
         Length_19 = IDSYMSZ(IRREPX,ISYTYP(1,19),ISYTYP(2,19))
         Length_20 = IDSYMSZ(IRREPX,ISYTYP(1,20),ISYTYP(2,20))
         Length_21 = IDSYMSZ(IRREPX,ISYTYP(1,21),ISYTYP(2,21))
         Length_22 = IDSYMSZ(IRREPX,ISYTYP(1,22),ISYTYP(2,22))
         
#ifdef _NOSKIP
         Call Getall(Work, Length_14, Irrepx, 14)
         Call Dscal(Length_14, Scale_factor, Work, 1)
CSSS         Call Dscal(Length_14, Zlich, Work, 1)
         Call Putall(Work, Length_14, Irrepx, 14)
#endif 

#ifdef _NOSKIP
         Call Getall(Work, Length_15, Irrepx, 15)
         Call Dscal(Length_15, Scale_factor, Work, 1)
CSS         Call Dscal(Length_15, Zlich, Work, 1)
        Call Putall(Work, Length_15, Irrepx, 15)
#endif 

C#ifdef _NOSKIP
         Call Getall(Work, Length_16, Irrepx, 16)
         Call Dscal(Length_16, Scale_factor, Work, 1)
         Call Putall(Work, Length_16, Irrepx, 16)
C#endif 

C#ifdef _NOSKIP 
         Call Getall(Work, Length_17, Irrepx, 17)
         Call Dscal(Length_17, Scale_factor, Work, 1)
         Call Putall(Work, Length_17, Irrepx, 17)
C#endif
C#ifdef _NOSKIP
         Call Getall(Work, Length_18, Irrepx, 18)
         Call Dscal(Length_18, Scale_factor, Work, 1)
         Call Putall(Work, Length_18, Irrepx, 18)
C#endif 

#ifdef _NOSKIP
         Call Getall(Work, Length_19, Irrepx, 19)
         Call Dscal(Length_19, Scale_factor, Work, 1)
CSSS         Call Dscal(Length_19, Zlich, Work, 1)
         Call Putall(Work, Length_19, Irrepx, 19)
#endif 

#ifdef _NOSKIP
         Call Getall(Work, Length_20, Irrepx, 20)
         Call Dscal(Length_20, Scale_factor, Work, 1)
CSSS         Call Dscal(Length_20, Zlich, Work, 1)
         Call Putall(Work, Length_20, Irrepx, 20)
#endif 

C#ifdef _NOSKIP 
         Call Getall(Work, Length_21, Irrepx, 21)
         Call Dscal(Length_21, Scale_factor, Work, 1)
CSSS         Call Dscal(Length_21, Zlich, Work, 1)
         Call Putall(Work, Length_21, Irrepx, 21)
C#endif 

C#ifdef _NOSKIP 
         Call Getall(Work, Length_22, Irrepx, 22)
         Call Dscal(Length_22, Scale_factor, Work, 1)
CSSS         Call Dscal(Length_21, Zlich, Work, 1)
         Call Putall(Work, Length_22, Irrepx, 22)
C#endif

         If (Ispar .and. Coulomb) Then
   
            Write(6,"(a,a,F7.5)") "<pp|hh> integrals for DCC are also",
     &                      " incremented by ", scale_factor 
            Length_119 = IDSYMSZ(IRREPX,ISYTYP(1,119),ISYTYP(2,119))
            Length_120 = IDSYMSZ(IRREPX,ISYTYP(1,120),ISYTYP(2,120))

C#ifdef _NOSKIP
            Call Getall(Work, Length_119, Irrepx, 119)
            Call Dscal(Length_119, Scale_factor, Work, 1)
            Call Putall(Work, Length_119, Irrepx, 119)
C#endif

C#ifdef _NOSKIP
           Call Getall(Work, Length_120, Irrepx, 120)
           Call Dscal(Length_120, Scale_factor, Work, 1)
           Call Putall(Work, Length_120, Irrepx, 120)
C#endif
         Endif

      Else

C For RHF Gamma(Ij|Ab)

         Length_14 = IDSYMSZ(IRREPX,ISYTYP(1,14),ISYTYP(2,14))
         Length_16 = IDSYMSZ(IRREPX,ISYTYP(1,16),ISYTYP(2,16))
         Length_18 = IDSYMSZ(IRREPX,ISYTYP(1,18),ISYTYP(2,18))
         Length_19 = IDSYMSZ(IRREPX,ISYTYP(1,19),ISYTYP(2,19))
         Length_21 = IDSYMSZ(IRREPX,ISYTYP(1,21),ISYTYP(2,21))

C For RHF, scale the lists 14,16,18,19 and 21. Make sure that
C 18 is not used (IA,JB) contributions 

#ifdef _NOSKIP
         Call Getall(Work, Length_14, Irrepx, 14)
         Call Dscal(Length_14, Scale_factor, Work, 1)
         Call Putall(Work, Length_14, Irrepx, 14)
#endif 

C#ifdef _NOSKIP
         Call Getall(Work, Length_16, Irrepx, 16)
         Call Dscal(Length_16, Scale_factor, Work, 1)
         Call Putall(Work, Length_16, Irrepx, 16)
C#endif 

#ifdef _NOSKIP
         Call Getall(Work, Length_18, Irrepx, 18)
         Call Dscal(Length_18, Scale_factor, Work, 1)
         Call Putall(Work, Length_18, Irrepx, 18)
#endif 

#ifdef _NOSKIP
         Call Getall(Work, Length_19, Irrepx, 19)
         Call Dscal(Length_19, Scale_factor, Work, 1)
         Call Putall(Work, Length_19, Irrepx, 19)
#endif 

C#ifdef _NOSKIP
         Call Getall(Work, Length_21, Irrepx, 21)
         Call Dscal(Length_21, Scale_factor, Work, 1)
         Call Putall(Work, Length_21, Irrepx, 21)
C#endif 

         If (Ispar .and. Coulomb) Then
            Write(6,"(a,a,F7.5)")"<pp|hh> integrals for DCC are also",
     &                      " incremented by ", scale_factor 
             Length_119 = IDSYMSZ(IRREPX,ISYTYP(1,119),ISYTYP(2,119))

C#ifdef _NOSKIP
             Call Getall(Work, Length_119, Irrepx, 119)
             Call Dscal(Length_119, Scale_factor, Work, 1)
             Call Putall(Work, Length_119, Irrepx, 119)
C#endif
         Endif

      Endif 
      Endif 

      Return 
      End


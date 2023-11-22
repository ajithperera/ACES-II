      Subroutine Scale_ints_driver(Work,Maxcor,IUhf)
   
      Implicit Integer (A-Z)

      Double Precision Work(Maxcor), Delta, Scale_factor,Scale 
      Double Precision Buf(600)
      Dimension Ibuf(600)
      Logical UHF
      Character*80 FNAME

#include "syminf.com"
#include "machsp.com"
#include "sympop.com"
#include "parcc.com"
#include "sym.com"
#include "info.com"
C
C This is a routine written for debugging purposes. What this 
C is to add delat to type of MO two electon integrals.

CSSS      Return
C
      Delta        = 0.001D0
      Scale        = 0.0D0
      Scale_factor = 1.0D0 + Delta*Scale 
      Zlich        = 0.0D0

      UHF    = .False.
      UHF    = (Iuhf .EQ. 1)
      Irrepx = 1 
      Imode  = 0 
      Write(6,*) 

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

C         Call Inipck(1,2,4, 148,Imode,0,1)
C         Call Inipck(1,13,14,149,Imode,0,1)
C         Length_147 = IDSYMSZ(IRREPX,ISYTYP(1,147),ISYTYP(2,147))
C         Length_148 = IDSYMSZ(IRREPX,ISYTYP(1,148),ISYTYP(2,148))
C         Length_149 = IDSYMSZ(IRREPX,ISYTYP(1,149),ISYTYP(2,149))
C         Call Getall(Work, Length_147, Irrepx, 14)
C         Call Putall(Work, Length_147, Irrepx, 147)
C         Call Getall(Work, Length_148, Irrepx, 15)
C         Call Putall(Work, Length_148, Irrepx, 148)
C         Call Getall(Work, Length_149, Irrepx, 16)
C         Call Putall(Work, Length_149, Irrepx, 149)
      Endif 

      Call Inipck(1,10,10,147,Imode,0,1)
      Length_147 = IDSYMSZ(IRREPX,ISYTYP(1,25),ISYTYP(2,25))
      Call Getall(Work, Length_147, Irrepx, 25)
      Call Putall(Work, Length_147, Irrepx, 147)
      Call Inipck(1,9,10,122,Imode,0,1)
      Length_122 = IDSYMSZ(IRREPX,ISYTYP(1,18),ISYTYP(2,18))
      Call Getall(Work, Length_122, Irrepx, 18)
      Call Putall(Work, Length_122, Irrepx, 122)
      Call Inipck(1,9,10,148,Imode,0,1)

C      Call Inipck(1,13,14,149,Imode,0,1)
C      Length_147 = IDSYMSZ(IRREPX,ISYTYP(1,147),ISYTYP(2,147))
C      Length_149 = IDSYMSZ(IRREPX,ISYTYP(1,149),ISYTYP(2,149))
C      Call Getall(Work, Length_147, Irrepx, 14)
C      Call Putall(Work, Length_147, Irrepx, 147)
C      Call Getall(Work, Length_149, Irrepx, 16)
C      Call Putall(Work, Length_149, Irrepx, 149)

#ifdef _SCALE_IJAB 

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

         Write(6,"(a)") " Lists 23,24,25 and 26 before scalling 14-22"
         
C         Call Dzero(Work,Length_23)
C         Call Putall(Work, Length_23, Irrepx, 23)
C         Call checksum("List-23",Work,Length_23)
C         Call Dzero(Work,Length_24)
C         Call Putall(Work, Length_24, Irrepx, 24)
C         Call checksum("List-24",Work,Length_24)
C         Call Dzero(Work,Length_25)
C         Call Putall(Work, Length_25, Irrepx, 25)
C         Call checksum("List-25",Work,Length_25)
C         Call Dzero(Work,Length_26)
C         Call Putall(Work, Length_26, Irrepx, 26)
C         Call checksum("List-26",Work,Length_26)
C         Call Dzero(Work,Length_17)
C         Call Putall(Work, Length_17, Irrepx, 121)
C         Call Dzero(Work,Length_18)
C         Call Putall(Work, Length_18, Irrepx, 122)

C         Call Dzero(Work,Length_17)
C         Call Putall(Work, Length_17, Irrepx, 17)
C         Call Dzero(Work,Length_18)
C         Call Putall(Work, Length_18, Irrepx, 18)
C         Call Dzero(Work,Length_21)
C         Call Putall(Work, Length_21, Irrepx, 21)
C         Call Dzero(Work,Length_22)
C         Call Putall(Work, Length_22, Irrepx, 22)
C
c         Call Dzero(Work,Length_19)
c         Call Putall(Work, Length_19, Irrepx, 19)
c         Call Dzero(Work,Length_20)
c         Call Putall(Work, Length_20, Irrepx, 20)

C         Call Dzero(Work,Length_14)
C         Call Putall(Work, Length_14, Irrepx, 14)
C         Call Dzero(Work,Length_15)
C         Call Putall(Work, Length_15, Irrepx, 15)
C         Call Dzero(Work,Length_16)
C         Call Putall(Work, Length_16, Irrepx, 16)
          
C
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
         Write(6,"(a)") " Lists 23,24,25 and 26 after scalling 14-22"
         Call Getall(Work, Length_23, Irrepx, 23)
         Call checksum("List-23",Work,Length_23)
         Call Getall(Work, Length_24, Irrepx, 24)
         Call checksum("List-24",Work,Length_24)
         Call Getall(Work, Length_25, Irrepx, 25)
         Call checksum("List-25",Work,Length_25)
         Call Getall(Work, Length_26, Irrepx, 26)
         Call checksum("List-26",Work,Length_26)
         Call Getall(Work, Length_121, Irrepx, 121)
         Call checksum("List-121",Work,Length_121)
         Call Getall(Work, Length_122, Irrepx, 122)
         Call checksum("List-122",Work,Length_122)

      Else

C For RHF Gamma(Ij|Ab)

         Length_14 = IDSYMSZ(IRREPX,ISYTYP(1,14),ISYTYP(2,14))
         Length_16 = IDSYMSZ(IRREPX,ISYTYP(1,16),ISYTYP(2,16))
         Length_18 = IDSYMSZ(IRREPX,ISYTYP(1,18),ISYTYP(2,18))
         Length_19 = IDSYMSZ(IRREPX,ISYTYP(1,19),ISYTYP(2,19))
         Length_21 = IDSYMSZ(IRREPX,ISYTYP(1,21),ISYTYP(2,21))

C For RHF, scale the lists 14,16,18,19 and 21. Make sure that
C 18 is not used (IA,JB) contributions 

C#ifdef _NOSKIP
         Call Getall(Work, Length_14, Irrepx, 14)
         Call Dscal(Length_14, Scale_factor, Work, 1)
         Call Putall(Work, Length_14, Irrepx, 14)
C#endif 

C#ifdef _NOSKIP
         Call Getall(Work, Length_16, Irrepx, 16)
         Call Dscal(Length_16, Scale_factor, Work, 1)
         Call Putall(Work, Length_16, Irrepx, 16)
C#endif 

C#ifdef _NOSKIP
         Call Getall(Work, Length_18, Irrepx, 18)
         Call Dscal(Length_18, Scale_factor, Work, 1)
         Call Putall(Work, Length_18, Irrepx, 18)
C#endif 

C#ifdef _NOSKIP
         Call Getall(Work, Length_19, Irrepx, 19)
         Call Dscal(Length_19, Scale_factor, Work, 1)
         Call Putall(Work, Length_19, Irrepx, 19)
C#endif 

C#ifdef _NOSKIP
         Call Getall(Work, Length_21, Irrepx, 21)
         Call Dscal(Length_21, Scale_factor, Work, 1)
         Call Putall(Work, Length_21, Irrepx, 21)
C#endif 
         Write(6,"(a)") " Lists 23 and 25 after scalling 14-22"
         Call Getall(Work, Length_23, Irrepx, 23)
         Call checksum("List-23",Work,Length_23)
         Call Getall(Work, Length_25, Irrepx, 25)
         Call checksum("List-25",Work,Length_25)
      Endif 
C
C For DCC methods, we have Coulomb only <pp|hh> lists 
C
      If (Ispar .and. Coulomb) Then

      Write(6,"(a,a,F7.5)") "<pp|hh> integrals for DCC are ",
     &                      "incremented by ", Delta*scale 
         If (Uhf) Then
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
         Else
             Length_119 = IDSYMSZ(IRREPX,ISYTYP(1,119),ISYTYP(2,119))
C#ifdef _NOSKIP
         Call Getall(Work, Length_18, Irrepx, 18)
         Call Dscal(Length_18, Scale_factor, Work, 1)
         Call Putall(Work, Length_18, Irrepx, 18)
C#endif 
C#ifdef _NOSKIP
         Call Getall(Work, Length_21, Irrepx, 21)
         Call Dscal(Length_21, Scale_factor, Work, 1)
         Call Putall(Work, Length_21, Irrepx, 21)
C#endif 
C#ifdef _NOSKIP
             Call Getall(Work, Length_119, Irrepx, 119)
             Call Dscal(Length_119, Scale_factor, Work, 1)
             Call Putall(Work, Length_119, Irrepx, 119)
C#endif 
         Endif 
      Endif 
#endif 

C#ifdef _SCALE_IBJA

      Write(6,"(a,F7.5)") "<ph||ph> integrals are incremented by ", 
     &                    delta*scale 
      If (Uhf) Then
         Length_23  = IDSYMSZ(IRREPX,ISYTYP(1,23),ISYTYP(2,23))
         Length_24  = IDSYMSZ(IRREPX,ISYTYP(1,24),ISYTYP(2,24))
         Length_25  = IDSYMSZ(IRREPX,ISYTYP(1,25),ISYTYP(2,25))
         Length_26  = IDSYMSZ(IRREPX,ISYTYP(1,26),ISYTYP(2,26))
         Length_121 = IDSYMSZ(IRREPX,ISYTYP(1,17),ISYTYP(2,17))
         Length_122 = IDSYMSZ(IRREPX,ISYTYP(1,18),ISYTYP(2,18))

C         Call Getall(Work, Length_23, Irrepx, 23)
C         Call Dscal(Length_23, Scale_factor, Work, 1)
C         Call Putall(Work, Length_23, Irrepx, 23)
C
C         Call Getall(Work, Length_24, Irrepx, 24)
C         Call Dscal(Length_24, Scale_factor, Work, 1)
C         Call Putall(Work, Length_24, Irrepx, 24)
C          
C         Call Getall(Work, Length_25, Irrepx, 25)
C         Call Dscal(Length_25, Scale_factor, Work, 1)
C         Call Putall(Work, Length_25, Irrepx, 25)
C
C         Call Getall(Work, Length_26, Irrepx, 26)
C         Call Dscal(Length_26, Scale_factor, Work, 1)
C         Call Putall(Work, Length_26, Irrepx, 26)

         Call Getall(Work, Length_121, Irrepx, 121)
         Call Dscal(Length_121, Scale_factor, Work, 1)
         Call Putall(Work, Length_121, Irrepx, 121)

         Call Getall(Work, Length_122, Irrepx, 122)
         Call Dscal(Length_122, Scale_factor, Work, 1)
         Call Putall(Work, Length_122, Irrepx, 122)

      Else

C RHF needs all list 122,147 (instead of 25 (used in G(ij,ab)
C 23 and 25 scalled. This yield the sum of three G(ia,jb) 
C contributions.

        Length_122 = IDSYMSZ(IRREPX,ISYTYP(1,18),ISYTYP(2,18))
        Length_147 = IDSYMSZ(IRREPX,ISYTYP(1,25),ISYTYP(2,25))
        Length_23  = IDSYMSZ(IRREPX,ISYTYP(1,23),ISYTYP(2,23))
        Length_25  = IDSYMSZ(IRREPX,ISYTYP(1,25),ISYTYP(2,25))

         Call Getall(Work, Length_122, Irrepx, 122)
         Call Dscal(Length_122, Scale_factor, Work, 1)
         Call Putall(Work, Length_122, Irrepx, 122)

         Call Getall(Work, Length_147, Irrepx, 147)
         Call Dscal(Length_147, Scale_factor, Work, 1)
         Call Putall(Work, Length_147, Irrepx, 147)

         Call Getall(Work, Length_23, Irrepx, 23)
         Call Dscal(Length_23, Scale_factor, Work, 1)
         Call Putall(Work, Length_23, Irrepx, 23)

         Call Getall(Work, Length_25, Irrepx, 25)
         Call Dscal(Length_25, Scale_factor, Work, 1)
         Call Putall(Work, Length_25, Irrepx, 25)

      Endif 
           
C#endif 
#ifdef _SCALE_IJKA

      Write(6,"(a,F7.5)") "<hh||hp> integrals are incremented by ", 
     &                    delta*scale 

      If (Uhf) Then
         Length_7 = IDSYMSZ(IRREPX,ISYTYP(1,7),ISYTYP(2,7))
         Length_8 = IDSYMSZ(IRREPX,ISYTYP(1,8),ISYTYP(2,8))
         Length_9 = IDSYMSZ(IRREPX,ISYTYP(1,9),ISYTYP(2,9))
         Length_10 = IDSYMSZ(IRREPX,ISYTYP(1,10),ISYTYP(2,10))

         Call Getall(Work, Length_7, Irrepx, 7)
         Call Dscal(Length_7, Scale_factor, Work, 1)
         Call Putall(Work, Length_7, Irrepx, 7)

         Call Getall(Work, Length_8, Irrepx, 8)
         Call Dscal(Length_8, Scale_factor, Work, 1)
         Call Putall(Work, Length_8, Irrepx, 8)
          
         Call Getall(Work, Length_9, Irrepx, 9)
         Call Dscal(Length_9, Scale_factor, Work, 1)
         Call Putall(Work, Length_9, Irrepx, 9)

         Call Getall(Work, Length_10, Irrepx, 10)
         Call Dscal(Length_10, Scale_factor, Work, 1)
         Call Putall(Work, Length_10, Irrepx, 10)
      Else
         Length_7 = IDSYMSZ(IRREPX,ISYTYP(1,7),ISYTYP(2,7))
         Length_10 = IDSYMSZ(IRREPX,ISYTYP(1,10),ISYTYP(2,10))

         Call Getall(Work, Length_7, Irrepx, 7)
         Call Dscal(Length_7, Scale_factor, Work, 1)
         Call Putall(Work, Length_7, Irrepx, 7)

         Call Getall(Work, Length_10, Irrepx, 10)
         Call Dscal(Length_10, Scale_factor, Work, 1)
         Call Putall(Work, Length_10, Irrepx, 10)
      Endif

#endif 

#ifdef _SCALE_ABCI

      Write(6,"(a,F7.5)") "<pp||ph> integrals are incremented by ", 
     &                    delta*scale 

      If (Uhf) Then
         Length_27 = IDSYMSZ(IRREPX,ISYTYP(1,27),ISYTYP(2,27))
         Length_28 = IDSYMSZ(IRREPX,ISYTYP(1,28),ISYTYP(2,28))
         Length_29 = IDSYMSZ(IRREPX,ISYTYP(1,29),ISYTYP(2,29))
         Length_30 = IDSYMSZ(IRREPX,ISYTYP(1,30),ISYTYP(2,30))

         Call Getall(Work, Length_27, Irrepx, 27)
         Call Dscal(Length_27, Scale_factor, Work, 1)
         Call Putall(Work, Length_7, Irrepx, 27)

         Call Getall(Work, Length_8, Irrepx, 28)
         Call Dscal(Length_28, Scale_factor, Work, 1)
         Call Putall(Work, Length_28, Irrepx, 28)

         Call Getall(Work, Length_29, Irrepx, 29)
         Call Dscal(Length_29, Scale_factor, Work, 1)
         Call Putall(Work, Length_29, Irrepx, 29)

         Call Getall(Work, Length_30, Irrepx, 30)
         Call Dscal(Length_30, Scale_factor, Work, 1)
         Call Putall(Work, Length_30, Irrepx, 30)
      Else
         Write(6,*) ISYTYP(1,27), ISYTYP(2,27)
CSSS         Length_27 = IDSYMSZ(IRREPX,ISYTYP(1,27),ISYTYP(2,27))
         Length_30 = IDSYMSZ(IRREPX,ISYTYP(1,30),ISYTYP(2,30))

CSSS         Call Getall(Work, Length_27, Irrepx, 27)
CSSS         Call Dscal(Length_27, Scale_factor, Work, 1)
CSSS         Call Putall(Work, Length_27, Irrepx, 27)

         Call Getall(Work, Length_30, Irrepx, 30)
         Call Dscal(Length_30, Scale_factor, Work, 1)
         Call Putall(Work, Length_30, Irrepx, 30)

      Endif

#endif
#ifdef _SCALE_IJKL

      Write(6,"(a,F7.5)") "<hh||hh> integrals are incremented by ", 
     &                    delta*scale 
      If (Uhf) Then
         Length_11 = IDSYMSZ(IRREPX,ISYTYP(1,11),ISYTYP(2,11))
         Length_12 = IDSYMSZ(IRREPX,ISYTYP(1,12),ISYTYP(2,12))
         Length_13 = IDSYMSZ(IRREPX,ISYTYP(1,13),ISYTYP(2,13))

#ifdef _NOSKIP      
         Call Getall(Work, Length_11, Irrepx, 11)
         Call Dscal(Length_11, Scale_factor, Work, 1)
         Call Putall(Work, Length_11, Irrepx, 11)
#endif 
#ifdef _NOSKIP      
         Call Getall(Work, Length_12, Irrepx, 12)
         Call Dscal(Length_12, Scale_factor, Work, 1)
         Call Putall(Work, Length_12, Irrepx, 12)
#endif 
C#ifdef _NOSKIP      
         Call Getall(Work, Length_13, Irrepx, 13)
         Call Dscal(Length_13, Scale_factor, Work, 1)
         Call Putall(Work, Length_13, Irrepx, 13)
C#endif 

      Else
         Length_13 = IDSYMSZ(IRREPX,ISYTYP(1,13),ISYTYP(2,13))

         Call Getall(Work, Length_13, Irrepx, 13)
         Call Dscal(Length_13, Scale_factor, Work, 1)
         Call Putall(Work, Length_13, Irrepx, 13)
      Endif

#endif

#ifdef _SCALE_ABCD

      Write(6,"(a,F7.5)") "<pp||pp> integrals are incremented by ", 
     &                    delta*scale 

      If (Uhf) Then
         Length_231 = IDSYMSZ(IRREPX,ISYTYP(1,231),ISYTYP(2,231))
         Length_232 = IDSYMSZ(IRREPX,ISYTYP(1,232),ISYTYP(2,232))
         Length_233 = IDSYMSZ(IRREPX,ISYTYP(1,233),ISYTYP(2,233))

         Call Getall(Work, Length_231, Irrepx, 231)
         Call Dscal(Length_231, Scale_factor, Work, 1)
         Call Putall(Work, Length_231, Irrepx, 231)

         Call Getall(Work, Length_232, Irrepx, 232)
         Call Dscal(Length_232, Scale_factor, Work, 1)
         Call Putall(Work, Length_232, Irrepx, 232)

         Call Getall(Work, Length_233, Irrepx, 233)
         Call Dscal(Length_233, Scale_factor, Work, 1)
         Call Putall(Work, Length_233, Irrepx, 233)

      Else
         Length_233 = IDSYMSZ(IRREPX,ISYTYP(1,233),ISYTYP(2,233))

         Call Getall(Work, Length_233, Irrepx, 233)
         Call Dscal(Length_233, Scale_factor, Work, 1)
         Call Putall(Work, Length_233, Irrepx, 233)

      Endif

#endif
      Return 
C
C One electron terms 
C
      Call Getrec(20,"JOBARC","NBASTOT ", 1, Norbs)
      LUNIT   = 10
      ILNBUF  = 600
      LDIM = NORBS*(NORBS+1)/2
      I000 = 1
      I010 = I000 + Ldim
      I020 = I010 + Norbs*Norbs

      CALL GFNAME('IIII    ',FNAME,ILENGTH)
      OPEN(LUNIT,FILE=FNAME(1:ILENGTH),FORM='UNFORMATTED',
     &           ACCESS='SEQUENTIAL')
      CALL LOCATE(LUNIT,'ONEHAMIL')
      CALL ZERO(Work(I010),LDIM)
      NUT = ILNBUF
      DO WHILE (NUT.EQ.ILNBUF)
         READ(LUNIT) BUF, IBUF, NUT
         DO INT = 1, NUT
            Work((I000-1)+(IBUF(INT))) = BUF(INT)
         END DO
      END DO
      CALL EXPND2(Work(I000),Work(I010),Norbs)

CSSS      Write(6,*) "The one electron ints in AO basis"
CSSS      call output(Work(I010),1,norbs,1,norbs,norbs,norbs,1)
#ifdef _DFULL
      Write(6,"(a,F7.5)") "<p|h|q> integrals are incremented by ", 
     &                    delta*scale 
      Call Dscal(Norbs*Norbs, Scale_factor, Work(I010), 1)
      Call Putrec(20,"JOBARC","ONEINTAO",Norbs*Norbs,Work(I010))
#endif 
    
      I030 = I020 +  Norbs*Norbs
      I040 = I030 +  Norbs*Norbs
      I050 = I040 +  2*Norbs*Norbs
      If(I050 .GT. Maxcor) Call Insmem("@-Scale_ints",I050,Maxcor)

      Call Ao2mo2(Work(I010),Work(I020),Work(I030),Work(I040),
     &                Norbs,Norbs,1)
      If (UHF) Call Ao2mo2(Work(I010),Work(I020),Work(I030),
     &                     Work(I040),Norbs,Norbs,2)

CSSS      Write(6,*) "The one electron ints in MO basis"
CSSS      call output(Work(I020),1,norbs,1,norbs,norbs,norbs,1)

#ifdef _OO_BLOCK
      Write(6,"(a,F7.5)") "<h|h|h> integrals are incremented by ", 
     &                    delta*scale 

      Call OO_block(Work(I020),Work(I000),Norbs,Pop,Vrt,Nocco,
     &              Nvrto,Nirrep,1)
      Call Dscal(Norbs*Norbs, Scale_factor, Work(I000), 1)

CSSS      Write(6,*) "The OO one electron ints in MO basis"
CSSS      call output(Work(I000),1,nocco(1),1,Nocco(1),norbs,norbs,1)
      
      Call Mo2ao2(Work(I000),Work(I010),Work(I030),Work(I040),
     &            Norbs,1)

CSSS      Write(6,*) "The OO one electron ints in AO basis"
CSSS      call output(Work(I010),1,norbs,1,norbs,norbs,norbs,1)

      Call Putrec(20,"JOBARC","ONEINTSA",Norbs*Norbs,Work(I010))

      If (Uhf) Then
         Call OO_block(Work(I020),Work(I000),Norbs,Pop,Vrt,Nocco,
     &                 Nvrto,Nirrep,2)
         Call Dscal(Norbs*Norbs, Scale_factor, Work(I000), 1)
         Call Mo2ao2(Work(I000),Work(I010),Work(I030),Work(I040),
     &               Norbs,2)
         Call Putrec(20,"JOBARC","ONEINTSB",Norbs*Norbs,Work(I010))
      Endif 
#endif 

#ifdef _VV_BLOCK
      Write(6,"(a,F7.5)") "<p|h|p> integrals are incremented by ", 
     &                    delta*scale 
      Call VV_block(Work(I020),Work(I000),Norbs,Pop,Vrt,Nocco,
     &                 Nvrto,Nirrep,1)
      Call Dscal(Norbs*Norbs, Scale_factor, Work(I000), 1)

CSSS      Write(6,*) "The VV one electron ints in MO basis"
CSSS      call output(Work(I000),1,norbs,1,norbs,norbs,norbs,1)

      Call Mo2ao2(Work(I000),Work(I010),Work(I030),Work(I040),
     &                Norbs,1)
      Call Putrec(20,"JOBARC","ONEINTSA",Norbs*Norbs,Work(I010))

      If (Uhf) Then
         Call VV_block(Work(I020),Work(I000),Norbs,Pop,Vrt,Nocco,
     &                 Nvrto,Nirrep,2)
         Call Dscal(Norbs*Norbs, Scale_factor, Work(I000), 1)
         Call Mo2ao2(Work(I000),Work(I010),Work(I030),Work(I040),
     &                Norbs,2)
         Call Putrec(20,"JOBARC","ONEINTSB",Norbs*Norbs,Work(I010))
      Endif 
#endif 

#ifdef _OV_BLOCK
      Write(6,"(a,F7.5)") "<h|h|p> integrals are incremented by ", 
     &                    delta*scale 
      Call OV_block(Work(I020),Work(I000),Norbs,Pop,Vrt,Nocco,
     &                 Nvrto,Nirrep,1)
      Call Dscal(Norbs*Norbs, Scale_factor, Work(I000), 1)

CSSS      Write(6,*) "The OV one electron ints in MO basis"
CSSS      call output(Work(I000),1,norbs,1,norbs,norbs,norbs,1)

      Call Mo2ao2(Work(I000),Work(I010),Work(I030),Work(I040),
     &                Norbs,1)
      Call Putrec(20,"JOBARC","ONEINTSA",Norbs*Norbs,Work(I010))

      If (Uhf) Then
         Call OV_block(Work(I020),Work(I000),Norbs,Pop,Vrt,Nocco,
     &                 Nvrto,Nirrep,2)
         Call Dscal(Norbs*Norbs, Scale_factor, Work(I000), 1)
         Call Mo2ao2(Work(I000),Work(I010),Work(I030),Work(I040),
     &                Norbs,2)
         Call Putrec(20,"JOBARC","ONEINTSB",Norbs*Norbs,Work(I010))
      Endif
#endif
      Return
      End
     
      Subroutine OO_Block(Hfull,Hblock,Norbs,Pop,Vrt,Nocco,Nvrto,
     &                    Nirrep, Ispin)

      Implicit Integer (A-Z)

      Double Precision Hfull(Norbs,Norbs), Hblock(Norbs,Norbs)
      Integer Pop(8,2),Vrt(8,2),Nocco(2),Nvrto(2)

      Ioff  = 0
      Call Dzero(Hblock,Norbs*Norbs)
      Do Irrep = 1, Nirrep
         Do J = 1, Pop(Irrep,Ispin)
            Do I = 1, Pop(Irrep,Ispin)
               Hblock(I+Ioff,J+Joff) = Hfull(I+Ioff,J+Ioff) 
            Enddo 
         Enddo 
         Ioff = Ioff + Pop(Irrep,Ispin)
      Enddo 
      Return
      End

      Subroutine VV_Block(Hfull,Hblock,Norbs,Pop,Vrt,Nocco,Nvrto,
     &                    Nirrep,Ispin)

      Implicit Integer (A-Z)

      Double Precision Hfull(Norbs,Norbs), Hblock(Norbs,Norbs)
      Integer Pop(8,2),Vrt(8,2),Nocco(2),Nvrto(2)

      Ioff   = Nocco(Ispin)
      Call Dzero(Hblock,Norbs*Norbs)
      Do Irrep = 1, Nirrep
         Do J = 1, Vrt(Irrep,Ispin)
            Do I = 1, Vrt(Irrep,Ispin)
               Hblock(I+Ioff,J+Ioff) = Hfull(I+Ioff,J+Ioff)
            Enddo
         Enddo
         Ioff = Ioff + Vrt(Irrep,Ispin)
      Enddo

      Return
      End

      Subroutine OV_Block(Hfull,Hblock,Norbs,Pop,Vrt,Nocco,Nvrto,
     &                    Nirrep,Ispin)

      Implicit Integer (A-Z)

      Double Precision Hfull(Norbs,Norbs), Hblock(Norbs,Norbs)
      Integer Pop(8,2),Vrt(8,2),Nocco(2),Nvrto(2)

      Ioff  = 0
      Aoff  = Nocco(Ispin)
      Call Dzero(Hblock,Norbs*Norbs)
      Do Irrep = 1, Nirrep
         Do I = 1, Pop(Irrep,Ispin)
            Do A = 1, Vrt(Irrep,Ispin)
               Hblock(A+Aoff,I+Ioff) = Hfull(A+Aoff,I+Ioff)
               Hblock(I+Ioff,A+Aoff) = Hfull(A+Aoff,I+Ioff)
            Enddo
         Enddo
         Ioff = Ioff + Pop(Irrep,Ispin)
         Aoff = Ioff + Vrt(Irrep,Ispin)
      Enddo

      Return
      End


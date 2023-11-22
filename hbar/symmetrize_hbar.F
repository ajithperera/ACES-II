      Subroutine Symmetrize_Hbar(Work,Maxcor,Iuhf)
 
      Implicit Integer(A-Z)
   
      Double Precision Work(Maxcor) 
      Logical UHF  
        
#include "machsp.com"
#include "sympop.com"
#include "sym.com"
#include "info.com"
#include "syminf.com"

      UHF    = (IUhf .EQ. 1)
      IRREPX = 1
C
C First symmetrize the Hbar(A,E), Hbar(M,I). The Hbar(M,E) is 
C symmetric for when T1 is zero.
C
      IRREPX =1
      HHA_LENGTH = IRPDPD(IRREPX,21)
      I000 = 1
      I010 = I000 + HHA_LENGTH
      Call Getlst(Work(I000),1,1,1,1,91)
      Ioff = 0
      Do Irrep = 1, Nirrep
         Nocca = Pop(Irrep,1)
         If (Nocca .NE. 0) Then
            Call Symmet2(Work(I000+Ioff),Nocca)
            Ioff = Ioff + Nocca * Nocca
         Endif 
      Enddo
      Call Putlst(Work(I000),1,1,1,1,91)

      PPA_LENGTH = IRPDPD(IRREPX,19)
      I000 = 1
      I010 = I000 + PPA_LENGTH
      CALL Getlst(Work(I000),1,1,1,1,92)
      Ioff = 0
      Do Irrep = 1, Nirrep
         Nvrta = Vrt(Irrep,1)
         If (Nvrta .NE. 0) Then
            Call Symmet2(Work(I000+Ioff),NVRTA)
            Ioff = Ioff + Nvrta * Nvrta
         Endif 
      Enddo
      CALL Putlst(Work(I000),1,1,1,1,92)

      IF (UHF) Then

         HHB_LENGTH = IRPDPD(IRREPX,22)
         I000 = 1
         I010 = I000 + HHB_LENGTH
         Call Getlst(Work(I000),1,1,1,2,91)
         Ioff = 0
         Do Irrep = 1, Nirrep
            Noccb = Pop(Irrep,2)
            If (Noccb .NE. 0) Then
               Call Symmet2(Work(I000+Ioff),Noccb)
               Ioff = Ioff + Noccb * Noccb
            Endif 
         Enddo
         Call Putlst(Work(I000),1,1,1,2,91)
 
         PPB_LENGTH = IRPDPD(IRREPX,20)
         I000 = 1
         I010 = I000 + PPB_LENGTH
         CALL Getlst(Work(I000),1,1,1,2,92)
         Ioff = 0
         Do Irrep = 1, Nirrep
            Nvrtb = Vrt(Irrep,2)
            If (Nvrtb .NE. 0) Then
               Call Symmet2(Work(I000+Ioff),NVRTB)
              Ioff = Ioff + Nvrtb * Nvrtb
            Endif 
         Enddo
         Call Putlst(Work(I000),1,1,1,2,92)
      ENDIF
      
C Hbar(AI,BC) and Hbar(AB,CI)

      IF (UHF) Then

         AAAA_LENGTH_ABCI = IDSYMSZ(IRREPX,ISYTYP(1,127),ISYTYP(2,127))
         BBBB_LENGTH_ABCI = IDSYMSZ(IRREPX,ISYTYP(1,128),ISYTYP(2,128))
         ABAB_LENGTH_ABCI = IDSYMSZ(IRREPX,ISYTYP(1,129),ISYTYP(2,129))
         ABBA_LENGTH_ABCI = IDSYMSZ(IRREPX,ISYTYP(1,130),ISYTYP(2,130))

         MEM_REQUIRED  = 1
         MEM_REQUIRED  = MAX(MEM_REQUIRED, AAAA_LENGTH_ABCI) 
         MEM_REQUIRED  = MAX(MEM_REQUIRED, BBBB_LENGTH_ABCI) 
         MEM_REQUIRED  = MAX(MEM_REQUIRED, ABAB_LENGTH_ABCI) 
         MEM_REQUIRED  = MAX(MEM_REQUIRED, ABBA_LENGTH_ABCI) 
         MEM_REQUIRED  = 3 * MEM_REQUIRED

         IF (MEM_REQUIRED .GT. MAXCOR) Call Insmem("@-Symmetrize_Hbar",
     &                                      Mem_Required,Maxcor)
         I000 = 1
         I010 = I000 + AAAA_LENGTH_ABCI
         I020 = I010 + AAAA_LENGTH_ABCI
         I030 = I020 + AAAA_LENGTH_ABCI

         Call Getall(Work(I000), AAAA_LENGTH_ABCI, IRREPX, 27)
         Call Getall(Work(I010), AAAA_LENGTH_ABCI, IRREPX, 127)
         Call Dcopy(AAAA_LENGTH_ABCI, Work(I000), 1, Work(I020), 1)
         Call Daxpy(AAAA_LENGTH_ABCI,1.0D0,Work(I010),1,Work(I000), 1)
         Call Daxpy(AAAA_LENGTH_ABCI,1.0D0,Work(I000),1,Work(I010), 1)
         Call Dscal(AAAA_LENGTH_ABCI,0.5D0,Work(I000),1)
         Call Dscal(AAAA_LENGTH_ABCI,0.5D0,Work(I010),1)
         Call Putall(Work(I000), AAAA_LENGTH_ABCI,IRREPX, 27)
         Call Putall(Work(I010), AAAA_LENGTH_ABCI,IRREPX, 127)

         I010 = I000 + BBBB_LENGTH_ABCI
         I020 = I010 + BBBB_LENGTH_ABCI
         I030 = I020 + BBBB_LENGTH_ABCI

         Call Getall(Work(I000), BBBB_LENGTH_ABCI, IRREPX, 28)
         Call Getall(Work(I010), BBBB_LENGTH_ABCI, IRREPX, 128)
         Call Dcopy(BBBB_LENGTH_ABCI, Work(I000), 1, Work(I020), 1)
         Call Daxpy(BBBB_LENGTH_ABCI,1.0D0,Work(I010),1,Work(I000), 1)
         Call Daxpy(BBBB_LENGTH_ABCI,1.0D0,Work(I000),1,Work(I010), 1)
         Call Dscal(BBBB_LENGTH_ABCI,0.5D0,Work(I000),1)
         Call Dscal(BBBB_LENGTH_ABCI,0.5D0,Work(I010),1)
         Call Putall(Work(I000), BBBB_LENGTH_ABCI,IRREPX, 28)
         Call Putall(Work(I010), BBBB_LENGTH_ABCI,IRREPX, 128)

         I010 = I000 + ABAB_LENGTH_ABCI
         I020 = I010 + ABAB_LENGTH_ABCI
         I030 = I020 + ABAB_LENGTH_ABCI

         Call Getall(Work(I000), ABAB_LENGTH_ABCI, IRREPX, 29)
         Call Getall(Work(I010), ABAB_LENGTH_ABCI, IRREPX, 129)
         Call Dcopy(ABAB_LENGTH_ABCI, Work(I000), 1, Work(I020), 1)
         Call Daxpy(ABAB_LENGTH_ABCI,1.0D0,Work(I010),1,Work(I000), 1)
         Call Daxpy(ABAB_LENGTH_ABCI,1.0D0,Work(I000),1,Work(I010), 1)
         Call Dscal(ABBA_LENGTH_ABCI,0.5D0,Work(I000),1)
         Call Dscal(ABBA_LENGTH_ABCI,0.5D0,Work(I010),1)
         Call Putall(Work(I000), ABAB_LENGTH_ABCI,IRREPX, 29)
         Call Putall(Work(I010), ABAB_LENGTH_ABCI,IRREPX, 129)

         I010 = I000 + ABBA_LENGTH_ABCI
         I020 = I010 + ABBA_LENGTH_ABCI
         I030 = I020 + ABBA_LENGTH_ABCI

         Call Getall(Work(I000), ABBA_LENGTH_ABCI, IRREPX, 30) 
         Call Getall(Work(I010), ABBA_LENGTH_ABCI, IRREPX, 130)
         Call Dcopy(ABBA_LENGTH_ABCI, Work(I000), 1, Work(I020), 1)
         Call Daxpy(ABBA_LENGTH_ABCI,1.0D0,Work(I010),1,Work(I000), 1)
         Call Daxpy(ABBA_LENGTH_ABCI,1.0D0,Work(I020),1,Work(I010), 1)
         Call Dscal(ABBA_LENGTH_ABCI,0.5D0,Work(I000),1)
         Call Dscal(ABBA_LENGTH_ABCI,0.5D0,Work(I010),1)
         Call Putall(Work(I000), ABBA_LENGTH_ABCI,IRREPX, 30)
         Call Putall(Work(I010), ABBA_LENGTH_ABCI,IRREPX, 130)

      Else

         ABBA_LENGTH_ABCI = IDSYMSZ(IRREPX,ISYTYP(1,30),ISYTYP(2,30))
         ABBA_LENGTH_ABCI = IDSYMSZ(IRREPX,ISYTYP(1,130),ISYTYP(2,130))

         I000 = 1
         I010 = I000 + ABBA_LENGTH_ABCI
         I020 = I010 + ABBA_LENGTH_ABCI
         I030 = I020 + ABBA_LENGTH_ABCI

         MEM_REQUIRED = 1
         MEM_REQUIRED = MAX(MEM_REQUIRED, ABBA_LENGTH_ABCI) 
         MEM_REQUIRED = 3 * MEM_REQUIRED

         IF (MEM_REQUIRED .GT. MAXCOR) Call Insmem("@-Symmetrize_Hbar",
     &                                      Mem_Required,Maxcor)

         Call Getall(Work(I000), ABBA_LENGTH_AIBC, IRREPX, 30) 
         Call Getall(Work(I010), ABBA_LENGTH_ABCI, IRREPX, 130)
         Call Dcopy(ABBA_LENGTH_ABCI, Work(I000), 1, Work(I020), 1)
         Call Daxpy(ABBA_LENGTH_ABCI,1.0D0,Work(I010),1,Work(I000), 1)
         Call Daxpy(ABBA_LENGTH_ABCI,1.0D0,Work(I020),1,Work(I010), 1)
         Call Dscal(ABBA_LENGTH_ABCI,0.5D0,Work(I000),1)
         Call Dscal(ABBA_LENGTH_ABCI,0.5D0,Work(I010),1)
         Call Putall(Work(I000), ABBA_LENGTH_ABCI,IRREPX, 30)
         Call Putall(Work(I010), ABBA_LENGTH_ABCI,IRREPX, 130)

         write(6,"(6(1X,F10.6))")(work(i000-1+i),i=1,ABBA_LENGTH_ABCI)
         Write(*,*)
         write(6,"(6(1X,F10.6))")(work(i010-1+i),i=1,ABBA_LENGTH_ABCI)

      Endif 

C Hbar(IJ,KA) and Hbar(IA,JK)

      IF (UHF) Then

         AAAA_LENGTH_IJKA = IDSYMSZ(IRREPX,ISYTYP(1,107),ISYTYP(2,107))
         BBBB_LENGTH_IJKA = IDSYMSZ(IRREPX,ISYTYP(1,108),ISYTYP(2,108))
         ABAB_LENGTH_IJKA = IDSYMSZ(IRREPX,ISYTYP(1,109),ISYTYP(2,109))
         ABBA_LENGTH_IJKA = IDSYMSZ(IRREPX,ISYTYP(1,110),ISYTYP(2,110))

         MEM_REQUIRED = 1
         MEM_REQUIRED = MAX(MEM_REQUIRED, AAAA_LENGTH_ABCI) 
         MEM_REQUIRED = MAX(MEM_REQUIRED, BBBB_LENGTH_ABCI) 
         MEM_REQUIRED = MAX(MEM_REQUIRED, ABAB_LENGTH_ABCI) 
         MEM_REQUIRED = MAX(MEM_REQUIRED, ABBA_LENGTH_ABCI) 
         MEM_REQUIED  = 3 * MEM_REQUIED

         IF (MEM_REQUIRED .GT. MAXCOR) Call Insmem("@-Symmetrize_Hbar",
     &                                      Mem_Required,Maxcor)

         I000 = 1
         I010 = I000 + AAAA_LENGTH_IJKA
         I020 = I010 + AAAA_LENGTH_IJKA
         I030 = I020 + AAAA_LENGTH_IJKA

         Call Getall(Work(I000), AAAA_LENGTH_IJKA, IRREPX, 7)
         Call Getall(Work(I010), AAAA_LENGTH_IJKA, IRREPX, 107)
         Call Dcopy(AAAA_LENGTH_IJKA, Work(I000), 1, Work(I020), 1)
         Call Daxpy(AAAA_LENGTH_IJKA,1.0D0,Work(I010),1,Work(I000), 1)
         Call Daxpy(AAAA_LENGTH_IJKA,1.0D0,Work(I020),1,Work(I010), 1)
         Call Dscal(AAAA_LENGTH_IJKA,0.5D0,Work(I000),1)
         Call Dscal(AAAA_LENGTH_IJKA,0.5D0,Work(I010),1)
         Call Putall(Work(I000), AAAA_LENGTH_IJKA,IRREPX, 7)
         Call Putall(Work(I010), AAAA_LENGTH_IJKA,IRREPX, 107)
     
         I010 = I000 + BBBB_LENGTH_IJKA
         I020 = I010 + BBBB_LENGTH_IJKA
         I030 = I020 + BBBB_LENGTH_IJKA

         Call Getall(Work(I000), BBBB_LENGTH_IJKA, IRREPX, 8)
         Call Getall(Work(I010), BBBB_LENGTH_IJKA, IRREPX, 108)
         Call Dcopy(BBBB_LENGTH_IJKA, Work(I000), 1, Work(I020), 1)
         Call Daxpy(BBBB_LENGTH_IJKA,1.0D0,Work(I010),1,Work(I000), 1)
         Call Daxpy(BBBB_LENGTH_IJKA,1.0D0,Work(I020),1,Work(I010), 1)
         Call Dscal(BBBB_LENGTH_IJKA,0.5D0,Work(I000),1)
         Call Dscal(BBBB_LENGTH_IJKA,0.5D0,Work(I010),1)
         Call Putall(Work(I000), BBBB_LENGTH_IJKA,IRREPX, 8)
         Call Putall(Work(I010), BBBB_LENGTH_IJKA,IRREPX, 108)

         I010 = I000 + ABAB_LENGTH_IJKA
         I020 = I010 + ABAB_LENGTH_IJKA
         I030 = I020 + ABAB_LENGTH_IJKA

         Call Getall(Work(I000), ABAB_LENGTH_IJKA, IRREPX, 9)
         Call Getall(Work(I010), ABAB_LENGTH_IJKA, IRREPX, 109)
         Call Dcopy(ABAB_LENGTH_IJKA, Work(I000), 1, Work(I020), 1)
         Call Daxpy(ABAB_LENGTH_IJKA,1.00D0,Work(I010),1,Work(I000), 1)
         Call Daxpy(ABAB_LENGTH_IJKA,1.00D0,Work(I020),1,Work(I010), 1)
         Call Dscal(BBBB_LENGTH_IJKA,0.5D0,Work(I000),1)
         Call Dscal(BBBB_LENGTH_IJKA,0.5D0,Work(I010),1)
         Call Putall(Work(I000), ABAB_LENGTH_IJKA,IRREPX, 9)
         Call Putall(Work(I010), ABAB_LENGTH_IJKA,IRREPX, 109)

         I010 = I000 + ABBA_LENGTH_IJKA
         I020 = I010 + ABBA_LENGTH_IJKA
         I030 = I020 + ABBA_LENGTH_IJKA

         Call Getall(Work, ABBA_LENGTH_IJKA, IRREPX, 10)
         Call Getall(Work(I010), ABBA_LENGTH_IJKA, IRREPX, 110)
         Call Dcopy(ABBA_LENGTH_IJKA, Work(I000), 1, Work(I020), 1)
         Call Daxpy(ABBA_LENGTH_IJKA,1.0D0,Work(I010),1,Work(I000), 1)
         Call Daxpy(ABBA_LENGTH_IJKA,1.00D0,Work(I020),1,Work(I010), 1)
         Call Dscal(BBBB_LENGTH_IJKA,0.5D0,Work(I000),1)
         Call Dscal(BBBB_LENGTH_IJKA,0.5D0,Work(I010),1)
         Call Putall(Work(I000), ABAB_LENGTH_IJKA,IRREPX, 10)
         Call Putall(Work(I010), ABAB_LENGTH_IJKA,IRREPX, 110)

      Else

         ABBA_LENGTH_IJKA = IDSYMSZ(IRREPX,ISYTYP(1,10),ISYTYP(2,10))

         MEM_REQUIERD = 1 
         MEM_REQUIRED = MAX(MEM_REQUIRED, ABAB_LENGTH_IJKA) 
         MEM_REQUIRED = 3 * MEM_REQUIRED

         IF (MEM_REQUIRED .GT. MAXCOR) Call Insmem("@-Symmetrize_Hbar",
     &                                      Mem_Required,Maxcor)
         I000 = 1
         I010 = I000 + ABBA_LENGTH_IJKA
         I020 = I010 + ABBA_LENGTH_IJKA
         I030 = I020 + ABBA_LENGTH_IJKA

         Call Getall(Work(I000), ABBA_LENGTH_IJKA, IRREPX, 10)
         Call Getall(Work(I010), ABBA_LENGTH_IJKA, IRREPX, 110)
         Call Dcopy(ABBA_LENGTH_IJKA, Work(I000), 1, Work(I020), 1)
         Call Daxpy(ABBA_LENGTH_IJKA,1.0D0,Work(I010),1,Work(I000), 1)
         Call Daxpy(ABBA_LENGTH_IJKA,1.0D0,Work(I020),1,Work(I010), 1)
         Call Dscal(ABBA_LENGTH_IJKA,0.5D0,Work(I000),1)
         Call Dscal(ABBA_LENGTH_IJKA,0.5D0,Work(I010),1)
         Call Putall(Work(I000), ABBA_LENGTH_IJKA,IRREPX, 10)
         Call Putall(Work(I010), ABBA_LENGTH_IJKA,IRREPX, 110)

         write(6,"(6(1X,F10.6))")(work(i000-1+i),i=1,ABBA_LENGTH_IJKA)
         Write(*,*)
         write(6,"(6(1X,F10.6))")(work(i010-1+i),i=1,ABBA_LENGTH_IJKA)
      Endif
C
C Symmetrize the Hbar(MB,EJ) 
C
      NDIM_AA  = 0
      NDIM_BB  = 0
      NDIM_AB  = 0
      NDIM_ABA = 0
      NDIM_ABB = 0

      IF (UHF) THEN
          DO IRREP = 1, NIRREP 
             NDIM_AA  = NDIM_AA + IRPDPD(IRREP,ISYTYP(1,54)) 
             NDIM_BB  = NDIM_BB + IRPDPD(IRREP,ISYTYP(1,55)) 
             NDIM_ABA = NDIM_AB + IRPDPD(IRREP,ISYTYP(1,58)) 
             NDIM_ABB = NDIM_AB + IRPDPD(IRREP,ISYTYP(1,59))
          ENDDO 
      ELSE
          DO IRREP = 1, NIRREP 
             NDIM_AA = NDIM_AA + IRPDPD(IRREP,ISYTYP(1,54)) 
             NDIM_AB = NDIM_AB + IRPDPD(IRREP,ISYTYP(1,58))
          ENDDO 
      ENDIF 

      IF (UHF) Then

         AAAA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,54),ISYTYP(2,54))
         BBBB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,55),ISYTYP(2,55))
         AABB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,56),ISYTYP(2,56))
         BBAA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,57),ISYTYP(2,57))
         ABAB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,58),ISYTYP(2,58))
         BABA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,59),ISYTYP(2,59))

         MEM_REQUIRED = 1
         MEM_REQUIRED = MAX(MEM_REQUIRED, AAAA_LENGTH_MBEJ) 
         MEM_REQUIRED = MAX(MEM_REQUIRED, BBBB_LENGTH_MBEJ) 
         MEM_REQUIRED = MAX(MEM_REQUIRED, AABB_LENGTH_MBEJ) 
         MEM_REQUIRED = MAX(MEM_REQUIRED, BBAA_LENGTH_MBEJ) 
         MEM_REQUIRED = MAX(MEM_REQUIRED, ABAB_LENGTH_MBEJ) 
         MEM_REQUIRED = MAX(MEM_REQUIRED, BABA_LENGTH_MBEJ) 

         IF (MEM_REQUIRED .GT. MAXCOR) Call Insmem("@-Symmetrize_Hbar",
     &                                      Mem_Required,Maxcor)

         I000 = 1
         I010 = I000 + AAAA_LENGTH_MBEJ
         Call Getall(Work, AAAA_LENGTH_MBEJ, IRREPX, 54)
         Ioff = 0
         Do I = 1, Nirrep
            NDIM_AA = IRPDPD(IRREP,ISYTYP(1,54))
            Call Symmet2(Work(I000+Ioff),NDIM_AA)
            Ioff = Ioff + NDIM_AA * NDIM_AA 
         Enddo 
         Call Putall(Work(I000), AAAA_LENGTH_MBEJ, IRREPX, 54)

         I000 = 1
         I010 = I000 + BBBB_LENGTH_MBEJ
         Ioff = 0
         Do I = 1, Nirrep
            NDIM_BB = IRPDPD(IRREP,ISYTYP(1,54))
            Call Symmet2(Work(I000+Ioff),NDIM_BB)
            Ioff = Ioff + NDIM_BB * NDIM_BB 
         Enddo 
         Call Putall(Work(I000), BBBB_LENGTH_MBEJ, IRREPX, 55)

         I000 = 1
         I010 = I000 + AABB_LENGTH_MBEJ
         I020 = I010 + AABB_LENGTH_MBEJ
         I030 = I020 + BBAA_LENGTH_MBEJ

         Call Getall(Work(I000), AABB_LENGTH_MBEJ, IRREPX, 56)
         Call Getall(Work(I020), BBAA_LENGTH_MBEJ, IRREPX, 57)
         Ioff = 0
         Do I = 1, Nirrep
            NRDIM = IRPDPD(IRREP,ISYTYP(1,54))
            NCDIM = IRPDPD(IRREP,ISYTYP(1,55))
            Call Transp(Work(I000+Ioff), Work(I010+Ioff), NRDIM, NCDIM)
            Ioff = Ioff + NRDIM * NCDIM 
         Enddo 
         Call Daxpy(BBAA_LENGTH_MBEJ,1.0D0,Work(I010),1,Work(I020),1)
         Call Putall(Work(I020), BBAA_LENGTH_MBEJ, IRREPX, 57)

         I000 = 1
         I010 = I000 + BBAA_LENGTH_MBEJ
         I020 = I010 + BBAA_LENGTH_MBEJ
         I030 = I020 + AABB_LENGTH_MBEJ

         Call Getall(Work(I000), BBAA_LENGTH_MBEJ, IRREPX, 57)
         Call Getall(Work(I020), AABB_LENGTH_MBEJ, IRREPX, 56)
         Ioff = 0
         Do I = 1, Nirrep
            NRDIM = IRPDPD(IRREP,ISYTYP(1,55))
            NCDIM = IRPDPD(IRREP,ISYTYP(1,54))
            Call Transp(Work(I000+Ioff), Work(I010+Ioff), NRDIM, NCDIM)
            Ioff = Ioff + NRDIM * NCDIM 
         Enddo 
         Call Daxpy(BBAA_LENGTH_MBEJ,1.0D0,Work(I010),1,Work(I020),1)
         Call Putall(Work(I020), BBAA_LENGTH_MBEJ, IRREPX, 57)

         I000 = 1
         I010 = I000 + ABAB_LENGTH_MBEJ
         Call Getall(Work(I000), ABAB_LENGTH_MBEJ, IRREPX, 58)
         Ioff = 0
         Do I = 1, Nirrep
            NDIM_ABA = IRPDPD(IRREP,ISYTYP(1,58))
            Call Symmet2(Work(I000+Ioff),NDIM_ABA)
            Ioff = Ioff + NDIM_ABA * NDIM_ABA
         Enddo 
         Call Putall(Work(I000), ABAB_LENGTH_MBEJ, IRREPX, 58)

         I000 = 1
         I010 = I000 + BABA_LENGTH_MBEJ
         Call Getall(Work(I000), BABA_LENGTH_MBEJ, IRREPX, 59)
         Ioff = 0
         Do I = 1, Nirrep
            NDIM_ABB = IRPDPD(IRREP,ISYTYP(1,59))
            Call Symmet2(Work(I000+Ioff),NDIM_ABB)
            Ioff = Ioff + NDIM_ABB * NDIM_ABB
         Enddo 
         Call Putall(Work(I000), BABA_LENGTH_MBEJ, IRREPX, 59)

      Else

         AAAA_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,54),ISYTYP(2,54))
         AABB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,56),ISYTYP(2,56))
         ABAB_LENGTH_MBEJ = IDSYMSZ(IRREPX,ISYTYP(1,58),ISYTYP(2,58))

         MEM_REQUIRED = 1
         MEM_REQUIRED = MAX(MEM_REQUIRED, AAAA_LENGTH_MBEJ) 
         MEM_REQUIRED = MAX(MEM_REQUIRED, ABAB_LENGTH_MBEJ) 
         MEM_REQUIRED = MAX(MEM_REQUIRED, AABB_LENGTH_MBEJ) 

         IF (MEM_REQUIRED .GT. MAXCOR) Call Insmem("@-Symmetrize_Hbar",
     &                                      Mem_Required,Maxcor)
         I000 = 1
         I010 = I000 + AAAA_LENGTH_MBEJ
         Call Getall(Work(I000), AAAA_LENGTH_MBEJ, IRREPX, 54)

         Ioff = 0
         Do IRREP = 1, Nirrep
            NDIM_AA = IRPDPD(IRREP,ISYTYP(1,54))
      call output(work(i000),1,ndim_aa,1,ndim_aa,ndim_aa,ndim_aa,1)
            Call Symmet2(Work(I000+Ioff),NDIM_AA)
      call output(work(i000),1,ndim_aa,1,ndim_aa,ndim_aa,ndim_aa,1)
            Ioff = Ioff + NDIM_AA * NDIM_AA
         Enddo
         Call Putall(Work(I000), AAAA_LENGTH_MBEJ, IRREPX, 54)


         I000 = 1
         I010 = I000 + AABB_LENGTH_MBEJ
         Call Getall(Work(I000), AABB_LENGTH_MBEJ, IRREPX, 56)
         Ioff = 0
         Do IRREP = 1, Nirrep
            NDIM_AA = IRPDPD(IRREP,ISYTYP(1,56))
      call output(work(i000),1,ndim_aa,1,ndim_aa,ndim_aa,ndim_aa,1)
            Call Symmet2(Work(I000+Ioff),NDIM_AA)
      call output(work(i000),1,ndim_aa,1,ndim_aa,ndim_aa,ndim_aa,1)
            Ioff = Ioff + NDIM_AA * NDIM_AA
         Enddo
         Call Putall(Work(I000), AABB_LENGTH_MBEJ, IRREPX, 56)

         I000 = 1
         I010 = I000 + ABAB_LENGTH_MBEJ
         Call Getall(Work(I000), ABAB_LENGTH_MBEJ, IRREPX, 58)
         Ioff = 0
         Do IRREP = 1, Nirrep
            NDIM_ABA = IRPDPD(IRREP,ISYTYP(1,58))
      call output(work(i000),1,ndim_aba,1,ndim_aba,ndim_aba,ndim_aba,1)
            Call Symmet2(Work(I000+Ioff),NDIM_ABA)
      call output(work(i000),1,ndim_aba,1,ndim_aba,ndim_aba,ndim_aba,1)
            Ioff = Ioff + NDIM_ABA * NDIM_ABA
         Enddo
         Call Putall(Work(I000), ABAB_LENGTH_MBEJ, IRREPX, 58)
      Endif

      Return
      End



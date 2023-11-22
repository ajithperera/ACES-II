      Subroutine test_hbar(
C     +                    Hbar_mia,Hbar_mib,Hbar_aea,Hbar_aeb,
C     +                    Hbar_mea,Hbar_meb,
C     +                    Hbar_ijkl1,
C     +                     Hbar_ijkl2,Hbar_ijkl3,Hbar_abef1,
C     +                     Hbar_abef2,Hbar_abef3,Hbar_mbej1, 
C     +                     Hbar_mbej2,Hbar_mbej3,Hbar_mbej4,
C     +                     Hbar_mbej5,Hbar_mbej6,
     +                     Hbar_abci1,
     +                     Hbar_abci2,Hbar_abci3,Hbar_abci4,
     +                     Hbar_aibc1,Hbar_aibc2,Hbar_aibc3,
     +                     Hbar_aibc4,
C     +                     Hbar_ijka1,Hbar_ijka2,
C     +                     Hbar_ijka3,Hbar_ijka4,Hbar_iajk1,
C     +                     Hbar_iajk2,Hbar_iajk3,Hbar_iajk4,
     +                     Nocc_a,Nocc_b,
     +                     Nvrt_a,Nvrt_b)


      Implicit Double Precision(A-H,O-Z)

C      Dimension Hbar_abef1(Nvrt_a,Nvrt_a,Nvrt_a,Nvrt_a)
C      Dimension Hbar_abef2(Nvrt_b,Nvrt_b,Nvrt_b,Nvrt_b)
C      Dimension Hbar_abef3(Nvrt_a,Nvrt_b,Nvrt_a,Nvrt_b)
C
C      Dimension Hbar_ijkl1(Nocc_a,Nocc_a,Nocc_a,Nocc_a)
C      Dimension Hbar_ijkl2(Nocc_b,Nocc_b,Nocc_b,Nocc_b)
C      Dimension Hbar_ijkl3(Nocc_a,Nocc_b,Nocc_a,Nocc_b)
C
C      Dimension Hbar_ijab1(Nocc_a,Nocc_a,Nvrt_a,Nvrt_a)
C      Dimension Hbar_ijab2(Nocc_b,Nocc_b,Nvrt_b,Nvrt_b)
C      Dimension Hbar_ijab3(Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)
C
C      Dimension Hbar_abij1(Nocc_a,Nocc_a,Nvrt_a,Nvrt_a)
C      Dimension Hbar_abij2(Nocc_b,Nocc_b,Nvrt_b,Nvrt_b)
C      Dimension Hbar_abij3(Nocc_a,Nocc_b,Nvrt_a,Nvrt_b)
C
C      Dimension Hbar_ijka1(Nocc_a,Nocc_a,Nocc_a,Nvrt_a)
C      Dimension Hbar_ijka2(Nocc_b,Nocc_b,Nocc_b,Nvrt_b)
C      Dimension Hbar_ijka3(Nocc_a,Nocc_b,Nocc_a,Nvrt_b)
C      Dimension Hbar_ijka4(Nocc_b,Nocc_a,Nocc_b,Nvrt_a)
C
C      Dimension Hbar_iajk1(Nocc_a,Nvrt_a,Nocc_a,Nocc_a)
C      Dimension Hbar_iajk2(Nocc_b,Nvrt_b,Nocc_b,Nocc_b)
C      Dimension Hbar_iajk3(Nocc_a,Nvrt_b,Nocc_a,Nocc_b)
C      Dimension Hbar_iajk4(Nocc_b,Nvrt_a,Nocc_b,Nocc_a)

      Dimension Hbar_abci1(Nvrt_a,Nvrt_a,Nvrt_a,Nocc_a)
      Dimension Hbar_abci2(Nvrt_b,Nvrt_b,Nvrt_b,Nocc_b)
      Dimension Hbar_abci3(Nvrt_a,Nvrt_b,Nvrt_a,Nocc_b)
      Dimension Hbar_abci4(Nvrt_b,Nvrt_a,Nvrt_b,Nocc_a)

      Dimension Hbar_aibc1(Nvrt_a,Nocc_a,Nvrt_a,Nvrt_a)
      Dimension Hbar_aibc2(Nvrt_b,Nocc_b,Nvrt_a,Nvrt_b)
      Dimension Hbar_aibc3(Nvrt_a,Nocc_b,Nvrt_a,Nvrt_b)
      Dimension Hbar_aibc4(Nvrt_b,Nocc_a,Nvrt_b,Nvrt_a)

C      Dimension Hbar_mbej1(Nocc_a,Nvrt_a,Nvrt_a,Nocc_a)
C      Dimension Hbar_mbej2(Nocc_b,Nvrt_b,Nvrt_b,Nocc_b)
C      Dimension Hbar_mbej3(Nocc_a,Nvrt_b,Nvrt_a,Nocc_b)
C      Dimension Hbar_mbej4(Nocc_b,Nvrt_a,Nvrt_b,Nocc_a)
C      Dimension Hbar_mbej5(Nocc_a,Nvrt_b,Nocc_a,Nvrt_b)
C      Dimension Hbar_mbej6(Nocc_b,Nvrt_a,Nocc_b,Nvrt_a)

C      Dimension Hbar_mia(Nocc_a,Nocc_a)
C      Dimension Hbar_mib(Nocc_b,Nocc_b)
C      Dimension Hbar_aea(Nvrt_a,Nvrt_a)
C      Dimension Hbar_aeb(Nvrt_b,Nvrt_b)
C      Dimension Hbar_mea(Nocc_a,Nvrt_a)
C      Dimension Hbar_meb(Nocc_b,Nvrt_b)

      Data  Ione /1/


#ifdef _DEBUG_LVL0
      L_aaaa =Nvrt_a*Nocc_a*Nvrt_a*Nvrt_a
      L_bbbb =Nvrt_b*Nocc_b*Nvrt_a*Nvrt_b
      L_abab =Nvrt_a*Nocc_b*Nvrt_a*Nvrt_b
      L_baba =Nvrt_b*Nocc_a*Nvrt_b*Nvrt_a
      Write(6,*) L_aaaa,L_bbbb,L_abab,L_baba
      call checksum("Hbar_abci1:",Hbar_abci1,L_aaaa)
      call checksum("Hbar_abci2:",Hbar_abci2,L_bbbb)
      call checksum("Hbar_abci3:",Hbar_abci3,L_abab)
      call checksum("Hbar_abci4:",Hbar_abci4,L_baba)
#endif

      Return
      End 

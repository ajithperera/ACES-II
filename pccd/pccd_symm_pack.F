      Subroutine Pccd_symm_pack(F,Opij,Opab,Opai,Opia,Nocc,Nvrt,Nbas)

      Implicit Double Precision(A-H,O-Z)

#include "sym.com"
#include "syminf.com"
#include "info.com"

      Dimension OPij(Nocc*Nocc)
      Dimension OPab(Nvrt*Nvrt)
      Dimension OPai(Nvrt*Nocc)
      Dimension OPia(Nocc*Nvrt)
      Dimension F(Nbas,Nbas)
      Dimension Ioffo(8)
      Dimension Ioffv(8)

      ISPIN   = 1
      IRREPX  = 1
      ITHRUIJ = 0
      ITHRUAB = 0
      ITHRUAI = 0
      ITHRUIA = 0
      IOFFO(1)= 0
      IOFFV(1)=NOCCO(ISPIN)

      DO IRREP=1,NIRREP-1
         IOFFO(IRREP+1)=IOFFO(IRREP)+POP(IRREP,ISPIN)
         IOFFV(IRREP+1)=IOFFV(IRREP)+VRT(IRREP,ISPIN)
      ENDDO 

      CALL DZERO(OPIJ,NOCC*NOCC)
      CALL DZERO(OPAB,NVRT*NVRT)
      CALL DZERO(OPAI,NVRT*NOCC)
      CALL DZERO(OPIA,NVRT*NOCC)

C OO part 

      DO IRREPR=1,NIRREP
        IRREPL=DIRPRD(IRREPR,IRREPX)
        DO J=1,POP(IRREPR,ISPIN)
         DO I=1,POP(IRREPL,ISPIN)
          ITHRUIJ=ITHRUIJ+1
          OPIJ(ITHRUIJ)=F(I+IOFFO(IRREPL),J+IOFFO(IRREPR))
         ENDDO
        ENDDO
      ENDDO 

C VV part 

      DO IRREPR=1,NIRREP
        IRREPL=DIRPRD(IRREPR,IRREPX)
        DO  J=1,VRT(IRREPR,ISPIN)
         DO I=1,VRT(IRREPL,ISPIN)
          ITHRUAB=ITHRUAB+1
          OPAB(ITHRUAB)=F(I+IOFFV(IRREPL),J+IOFFV(IRREPR))
         ENDDO
        ENDDO
      ENDDO 

C VO part 

      DO IRREPR=1,NIRREP
        IRREPL=DIRPRD(IRREPR,IRREPX)
        DO I=1,POP(IRREPR,ISPIN)
         DO A=1,VRT(IRREPL,ISPIN)
          ITHRUAI=ITHRUAI+1
          OPAI(ITHRUAI)=F(A+IOFFV(IRREPL),I+IOFFO(IRREPR))
         ENDDO
        ENDDO
      ENDDO 

C OV part 

      DO  IRREPR=1,NIRREP
        IRREPL=DIRPRD(IRREPR,IRREPX)
        DO  I=1,POP(IRREPR,ISPIN)
         DO A=1,VRT(IRREPL,ISPIN)
          ITHRUIA=ITHRUIA+1
          OPIA(ITHRUIA)=F(I+IOFFO(IRREPL),A+IOFFV(IRREPR))
         ENDDO
        ENDDO
      ENDDO 

      RETURN
      END 

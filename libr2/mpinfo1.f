      SUBROUTINE MPINFO1(IOFFT2,IOFFINT,ITOTSIZ,ISPIN)
C
C COMPUTE ABSOLUTE STARTING ADDRESSES TO EACH IRREP OF T2 AND
C W(XX,XX) AO INTEGRAL VECTORS.
C
      IMPLICIT INTEGER (A-Z)
      DIMENSION IOFFT2(8),IOFFINT(8)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMINF/NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON/AOSYM/IAOPOP(8),IOFFAO(8),IOFFV(8,2),IOFFO(8,2),
     &             IRPDPDAO(8),IRPDPDAOS(8),ISTART(8,8),ISTARTMO(8,3)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/INFO/NOCCO(2),NVRTO(2)
C
      IOFFT2(1)=1
      IOFFINT(1)=1
      LISTT=43+ISPIN
      DO 10 IRREP=2,NIRREP
       IOFFT2(IRREP)=IOFFT2(IRREP-1)+
c     &               IRPDPD(IRREP-1,ISYTYP(1,LISTT))*
     &               IRPDPDAO(IRREP-1)*
     &               IRPDPD(IRREP-1,ISYTYP(2,LISTT))
       IOFFINT(IRREP)=IOFFINT(IRREP-1)+
     &                IRPDPDAO(IRREP-1)*IRPDPDAO(IRREP-1) 
10    CONTINUE
      ITOTSIZ=IOFFINT(NIRREP)+IRPDPDAO(NIRREP)*IRPDPDAO(NIRREP)-1
      RETURN
      END        
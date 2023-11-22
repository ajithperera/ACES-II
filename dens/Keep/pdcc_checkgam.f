










      SUBROUTINE PDCC_CHECKGAM(ICORE,LISTC,LISTW,LISTG,FACT,
     &                         ISPIN)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER DIRPRD,DISSYW
      DIMENSION ICORE(*)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
      COMMON/ADD/SUM

c sym.com : begin
      integer      pop(8,2), vrt(8,2), nt(2), nfmi(2), nfea(2)
      common /sym/ pop,      vrt,      nt,    nfmi,    nfea
c sym.com : end
c info.com : begin
      integer       nocco(2), nvrto(2)
      common /info/ nocco,    nvrto
c info.com : end

      logical ispar,coulomb
      double precision paralpha, parbeta, pargamma
      double precision pardelta, Parepsilon
      double precision Fae_scale,Fmi_scale,Wmnij_scale,Wmbej_scale
      double precision Gae_scale,Gmi_scale
      common/parcc_real/ paralpha,parbeta,pargamma,pardelta,Parepsilon
      common/parcc_log/ ispar,coulomb
      common/parcc_scale/Fae_scale,Fmi_scale,Wmnij_scale,Wmbej_scale,
     &                   Gae_scale,Gmi_scale 


CSSS      RETURN
      E=0.0D+0
      E1=0.0D+0
      E2=0.0D+0
      sum = 0.0D+0
      
CSSS#ifdef _NOSKIP
      LEN_C = IDSYMSZ(1,ISYTYP(1,LISTC),ISYTYP(2,LISTC))
      LEN_WC= IDSYMSZ(1,ISYTYP(1,LISTW),ISYTYP(2,LISTW))
      Len_14 = IDSYMSZ(1,ISYTYP(1,LISTW),ISYTYP(2,LISTW))

      I000   = 1
      I010   = I000 + NIRREP 
      I020   = I010 + NIRREP 

      If (Ispin .EQ. 1) Then
         CALL GETREC(20,'JOBARC','SVAVA2X ',NIRREP,ICORE(I000))
         CALL GETREC(20,'JOBARC','SOAOA2X ',NIRREP,ICORE(I010))
      Else
         CALL GETREC(20,'JOBARC','SVBVB2X ',NIRREP,ICORE(I000))
         CALL GETREC(20,'JOBARC','SOBOB2X ',NIRREP,ICORE(I010))
      Endif 
      LEN_C = IDSYMSZ(1,ISYTYP(1,LISTC),ISYTYP(2,LISTC))
      LEN_WC= IDSYMSZ(1,ISYTYP(1,LISTW),ISYTYP(2,LISTW))
      LEN_G= IDSYMSZ(1,ISYTYP(1,LISTG),ISYTYP(2,LISTG))

      IF (ISPIN .EQ. 1) THEN
         LEN_W = IDSYMSZ(1,19,21)
      ELSE 
         LEN_W = IDSYMSZ(1,20,22)
      ENDIF

      LEN1 = 0
      LEN2 = 0
      LEN3 = 0
C
      DO IRREPR = 1, NIRREP
         IRREPL = DIRPRD(IRREPR,1)
         IF (ISPIN .EQ. 1) THEN
             LEN1 = LEN1 + IRPDPD(IRREPL,19)
             LEN2 = LEN2 + IRPDPD(IRREPR,21)
             LEN3 = LEN3 + IRPDPD(IRREPR,3)
CSSS             LEN3 = LEN3 + IRPDPD(IRREPR,7)
         ELSE 
             LEN1 = LEN1 + IRPDPD(IRREPL,20)
             LEN2 = LEN2 + IRPDPD(IRREPR,22)
             LEN3 = LEN3 + IRPDPD(IRREPR,4)
CSSS             LEN3 = LEN3 + IRPDPD(IRREPR,8)
         ENDIF 
      ENDDO 

      NSCRSZ=NVRTO(ISPIN)*NVRTO(ISPIN)+NOCCO(ISPIN)*NOCCO(ISPIN)+
     &       NVRTO(ISPIN)*NOCCO(ISPIN)
      
      I030   = I020 + LEN_WC
      I040   = I030 + LEN_W
      I050   = I040 + NSCRSZ
C
      CALL GETALL(ICORE(I020), LEN_C, 1, LISTC) 
CSSS      call checksum("W(AI,BJ)",ICORE(I020),LEN_C)
C
C Switch AIBJ to (AB,IJ)
C
CSSS      Write(6,*) "<AI|BJ>"
CSSS      call print(Icore(i020),len_c)
      CALL SSTGEN(ICORE(I020),ICORE(I030),LEN_C,VRT(1,ISPIN),
     &            POP(1,ISPIN),VRT(1,ISPIN),POP(1,ISPIN),
     &            ICORE(I040),1,"1324") 

      CALL DCOPY(LEN_W, ICORE(I030), 1, ICORE(I020), 1)
      call checksum("W(AB,IJ)",ICORE(I020),LEN_W)
CSSS      Write(6,*) "<AB|IJ>"
CSSS      call print(Icore(i020),len_w)


CSSS#ifdef _NOSKIP
      I030   = I020 + LEN_W
      I040   = I030 + LEN_W
      I050   = I040 + LEN_W
      CALL TRANSP(ICORE(I020),ICORE(I030),LEN2,LEN1)

      IOFF1 = I030 
      IOFF2 = I020 
      
      DO IRREPR = 1, NIRREP
         IRREPL = DIRPRD(IRREPR,1)
         NUMDISC=IRPDPD(IRREPR,18+ISPIN)

         NDSSIZC=IRPDPD(IRREPL,20+ISPIN)
CSSS         NDSSIZW=IRPDPD(IRREPL,6+ISPIN)
         NDSSIZW=IRPDPD(IRREPL,2+ISPIN)

         CALL SQSYM(IRREPL,POP(1,ISPIN),NDSSIZW,NDSSIZC,NUMDISC,
     &              ICORE(IOFF2),ICORE(IOFF1))
      ENDDO 
     
      CALL TRANSP(ICORE(I020),ICORE(I030),LEN1,LEN3)

CSSS      call checksum("W(AB,IJ)",ICORE(I030),LEN_W)

      IOFF1 = I030 
      IOFF2 = I040

      DO IRREPR = 1, NIRREP
         IRREPL = DIRPRD(IRREPR,1)

CSSS         NUMDISC=IRPDPD(IRREPR,6+ISPIN)
         NUMDISC=IRPDPD(IRREPR,2+ISPIN)
         NDSSIZC=IRPDPD(IRREPL,18+ISPIN)
CSSS         NDSSIZW=IRPDPD(IRREPL,4+ISPIN)
         NDSSIZW=IRPDPD(IRREPL,0+ISPIN)


         CALL SQSYM(IRREPL,VRT(1,ISPIN),NDSSIZW,NDSSIZC,NUMDISC,
     &             ICORE(IOFF2),ICORE(IOFF1))

      ENDDO 

      CALL PUTALL(ICORE(I040), LEN_WC, 1, LISTW)
CSSS      write(6,*) "<A<B|I<J>"
CSSS      call print(Icore(i040), len_wc) 
      call checksum("W(A<B,I<J)",ICORE(I040),LEN_WC)
CSSS#endif 
CSSS#endif 

      DO 1000 IRREP=1,NIRREP
      NUMSYW=IRPDPD(IRREP,ISYTYP(2,LISTW))
      DISSYW=IRPDPD(IRREP,ISYTYP(1,LISTW))
      IOFFW=1
      IOFFW2=1+NUMSYW*DISSYW*IINTFP
      CALL GETLST(ICORE(IOFFW),1,NUMSYW,1,IRREP,LISTW)
      CALL GETLST(ICORE(IOFFW2),1,NUMSYW,2,IRREP,LISTG)
CSSS      write(6,*) "G(Ij,ab)"
CSSS      call print(Icore(ioffw2), NUMSYW*DISSYW) 
      E=E+SDOT(NUMSYW*DISSYW,ICORE(IOFFW),1,ICORE(IOFFW2),1)
      E1=E1+SDOT(NUMSYW*DISSYW,ICORE(IOFFW),1,ICORE(IOFFW),1)
      E2=E2+SDOT(NUMSYW*DISSYW,ICORE(IOFFW2),1,ICORE(IOFFW2),1)
1000  CONTINUE
      sum=sum+FACT*e
      write(6,500)e,e1,e2
      write(6,501)sum
      return
500   format(' Energy contribution : ',3(F13.10,1X))
501   format(' Cumulative energy   : ',F20.10)
      end

      subroutine print(w,Len)

      Implicit double precision (a-h,o-z)
      Dimension W(Len)

      Write(6,"(5(1x,F15.8))") (W(i), i=1,Len)

      Return
      End 

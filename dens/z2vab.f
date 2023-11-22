      SUBROUTINE Z2VAB(G,W,T,AIVV,FACT,ISPIN,POP1,POP2,VRT1,VRT2,
     &                 DISSYT,NUMSYT,LISTG,LISTW,IRREP,TMP,IUHF,
     &                 LARGE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL LARGE
      INTEGER DISSYT,DIRPRD,POP1,POP2,VRT1,VRT2
      DIMENSION G(NUMSYT,1),W(NUMSYT,2),T(DISSYT,1),AIVV(1),TMP(1)
      DIMENSION POP1(8),POP2(8),VRT1(8),VRT2(8) 
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
      common /dropgeo/ ndrgeo
C
      data ONE,ONEM,HALF,TWO /1.0D0,-1.0D0,0.5D0,2.0D+0/
C
      CALL GETLST(T,1,NUMSYT,1,IRREP,LISTG)
      CALL TRANSP(T,G,NUMSYT,DISSYT)
C
      CALL GETLST(T,1,NUMSYT,2,IRREP,LISTW)
      CALL TRANSP(T,W,NUMSYT,DISSYT)
C
c       call checksum("Z2VAB listw =",W, numsyt*dissyt)
c       call checksum("Z2VAB listg =",G, numsyt*dissyt)
c      print *,"numsyt, dissyt = ",numsyt, dissyt
C
C  HOWEVER, IN RHF TAKE ADVANTAGE IF THE SPIN SYMMETRY
C
      IF(IUHF.EQ.1.AND.ISPIN.EQ.1) THEN
        CALL SYMTR1(IRREP,VRT1,VRT2,NUMSYT,G,TMP,TMP(1+NUMSYT),
     &              TMP(1+2*NUMSYT))
        CALL SYMTR1(IRREP,VRT1,VRT2,NUMSYT,W,TMP,TMP(1+NUMSYT),
     &              TMP(1+2*NUMSYT))
      ENDIF
C
C JOFF OFFSET IN THE VIRTUAL-VIRTUAL BLOCK OF I AND W
C IOFF OFFSET IN AIVV
C
      IOFF=1
      JOFF=1
      DO 90 IRREPI=1,NIRREP
C          
C GET NUMBER OF VIRTUAL ORBITALS FOR IRREPJ
C
       NVRTI=VRT1(IRREPI)
C
C DETERMINE IRREPI WHOSE DIRECT PRODUCT WITH IRREPJ GIVES IRREP
C
       IRREPJ=DIRPRD(IRREP,IRREPI)
C
C GET NUMBER OF VIRTUAL ORBITALS FOR KRREP
C
       NVRTJ=VRT2(IRREPJ)
C
c       Print*, NUMSYT*NVRTJ 
       IF(MIN(NVRTJ,NVRTI).NE.0) THEN
C
          CALL XGEMM('T','N',NVRTI,NVRTI,NUMSYT*NVRTJ,FACT,
     &               W(1,JOFF),NVRTJ*NUMSYT,G(1,JOFF),
     &               NVRTJ*NUMSYT,ONE,AIVV(IOFF),NVRTI)
       ENDIF

c         print *,"sizes: ",nvrti, numsyt
c         do i=1,nvrti*nvrti
c          print *,"I(",i,") = ",aivv(ioff+i)
c         enddo
C
       JOFF=JOFF+NVRTJ*NVRTI
       IOFF=IOFF+NVRTI*NVRTI
90    CONTINUE
C
      RETURN
      END

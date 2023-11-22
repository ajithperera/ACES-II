      SUBROUTINE PCCD_FORM_FD(DOO,DVV,DVO,AIOO,AIVV,AIOV,XAI,
     &                        ICORE,MAXCOR,IUHF,NONHF_TERMS_EXIST)
C
C   THIS ROUTINE CALCULATES THE FOLLOWING CONTRIBUTIONS TO THE
C   INTERMEDIATES
C
C    I(P,Q) <= - 1/2 SUM R [ f(P,R) D(Q,R) + f(R,P) D(R,Q) ]
C
C   THERE ARE ACTUALLY THE FOLLOWING TERMS TO CONSIDER
C
C    I(I,A) <= - SUM M f(I,M) D(A,M)     (NON-HF)
C
C              - SUM E f(I,E) D(A,E)     (NON-HF)
C
C    - I(A,I) <=  SUM E f(A,E) D(E,I)    (NON-HF)
C
C                  + SUM M f(A,M) D(M,I)  (NON-HF)
C
C    I(I,J) <= - SUM M f(I,M) D(J,M)
C
C              - SUM E f(I,E) D(J,E)     (NON-HF)
C
C    I(A,B) <= - SUM M f(A,M) D(B,M)     (NON-HF)
C
C              - SUM E f(A,E) D(B,E)
C
C  FOR QRHF WE HAVE TO ADD HERE SOME TERMS WHICH TAKE CARE
C  OF THE RELAXATION EFFECTS IN THE REFERENCE STATE. THIS
C  TERM IS (CONSIDERING THE STRUCTURE OF F) ESSENTIALLY THE
C  SAME AS THE CONTRIBUTION OF THE ENERGY WEIGHTED DENSITY
C  MATRIX. IN HF THIS HOWEVER CONTRIBUTES ONLY TO THE OVERLAP
C  TERMS AND NOT TO TERMS WHICH INVOLVE ORBITAL RELAXATIONS.
C  NOTE THAT THE QRHF REFERENCEE IS NOT STATIONARY WITH
C  RESPECT TO ORBITAL VARIATIONS
C
C  ACTUALLY WE HAVE TO ADD
C
C   I(I,J) <= - f(I,J)      (FOR HF I(I,I) <= -f(I,I))
C
C   I(A,I) <= - f(A,I)       (NON-HF)
C
C  FOR CANONICAL ORBITALS, THE CONTRIBUTION REDUCES TO
C
C    I(I,I) <= -f(I,I) D(I,I)
C
C    I(A,B) <= -f(A,A) D(A,A)
C
CEND
C
C CODED OCT/90 JG
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION ICORE(MAXCOR),DOO(1),DVV(1),DVO(1),AIOO(1),
     &          AIVV(1),AIOV(1),XAI(1)
      LOGICAL QRHF,ROHF,NONHF_TERMS_EXIST
      INTEGER DIRPRD,POP,VRT
      DIMENSION I0FAI(2),I0FAE(2),I0FMI(2),I0OO(2),I0VV(2),I0VO(2)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/INFO/NOCCO(2),NVRTO(2)
      COMMON/SYMINF/NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
      COMMON /SYM/  POP(8,2),VRT(8,2),NT(2),NF1(2),NF2(2)
C
      DATA ONE,ONEM /1.0D+0, -1.0D+0/
C
      MXCOR=MAXCOR
      ISHIFT = 0
C
C   SKIP THE FIRST PART FOR HF REFERENCE FUNCTIONS
C
      Write(6,"(a)") "---------Printing from pCCD_form_df-------"
      Write(6,"(a,3l)") "Logicals: QRHF ROHF NONHF_TERMS_EXIST:", 
     +                QRHF,ROHF,NONHF_TERMS_EXIST

      Call Dzero(Aioo,nf1(1))
      Call Dzero(Aivv,nf2(1))
      Call Dzero(aiov,Nt(1)) 
      Call Dzero(xai,Nt(1)) 

      IF(QRHF.OR.ROHF .OR. NONHF_TERMS_EXIST) THEN
C
C  GET THE FOCK MATRICES
C
      I0FAI(1)=MXCOR+1-IINTFP*NT(1)   
      I0FAE(1)=I0FAI(1)-IINTFP*NF2(1)
      I0FMI(1)=I0FAE(1)-IINTFP*NF1(1)
      MXCOR=MXCOR-IINTFP*(NT(1)+NF1(1)+NF2(1))

      CALL GETLST(ICORE(I0FMI(1)),1,1,1,5,91+ISHIFT)
      CALL GETLST(ICORE(I0FAE(1)),1,1,1,5,92+ISHIFT)
      CALL GETLST(ICORE(I0FAI(1)),1,1,1,5,93+ISHIFT)
      IF(IUHF.EQ.0) THEN
       I0FAI(2)=I0FAI(1)
       I0FAE(2)=I0FAE(1)
       I0FMI(2)=I0FMI(1)
      ELSE
       I0FAI(2)=I0FMI(1)-IINTFP*NT(2)
       I0FAE(2)=I0FAI(2)-IINTFP*NF2(2)
       I0FMI(2)=I0FAE(2)-IINTFP*NF1(2)
       MXCOR=MXCOR-IINTFP*(NT(2)+NF1(2)+NF2(2))
       CALL GETLST(ICORE(I0FMI(2)),1,1,1,6,91+ISHIFT)
       CALL GETLST(ICORE(I0FAE(2)),1,1,1,6,92+ISHIFT)
       CALL GETLST(ICORE(I0FAI(2)),1,1,1,6,93+ISHIFT)    
      ENDIF

#ifdef _DEBUG_LVLM 
      call checksum("@-dfini-91",ICORE(I0FMI(1)),NF1(1))
      call checksum("@-dfini-92",ICORE(I0FAE(1)),NF2(1))
      call checksum("@-dfini-dv",Dvv,NF2(1))
      call checksum("@-dfini-93",ICORE(I0FAI(1)),NT(1))
      call checksum("@-dfini-do",Dvo,NT(1))
#endif 
C
      I0OO(1)=1
      I0VV(1)=1
      I0VO(1)=1
      I0OO(2)=1+NF1(1)*IUHF
      I0VV(2)=1+NF2(1)*IUHF
      I0VO(2)=1+NT(1)*IUHF

C
C   LOOP OVER ALL SPINS
C
      DO 1000 ISPIN=1,IUHF+1
C
C   FIRST CALCULATE FIRST CONTRIBUTION TO I(I,A)
C
C    I(I,A) <= - SUM M f(I,M) D(A,M)
C
       IOFFI=I0VO(ISPIN)
       IOFFD=I0VO(ISPIN)
       IOFFF=I0FMI(ISPIN)
C
       DO 100 IRREP=1,NIRREP
C
        NOCC=POP(IRREP,ISPIN)
        NVRT=VRT(IRREP,ISPIN)
C
        IF(MIN(NOCC,NVRT).NE.0) THEN
C
#ifdef _DEBUG_LVLM
        Print*, "Dvo"
        Call output(Dvo(ioffd),1,nvrt,1,nocc,nvrt,nocc,1)
        Print*, "foo"
        Call output(icore(iofff),1,nocc,1,nocc,nocc,nocc,1)
#endif 
         CALL XGEMM('N','N',NVRT,NOCC,NOCC,ONEM,DVO(IOFFD),NVRT,
     &              ICORE(IOFFF),NOCC,ONE,AIOV(IOFFI),NVRT)
        ENDIF

#ifdef _DEBUG_LVL0
        Print*, "Aiov-1"
        Call output(Aiov(ioffi),1,nocc,1,nvrt,nocc,nvrt,1)
#endif 
C
C  UPDATE OFFSETS
C
        IOFFI=IOFFI+NOCC*NVRT
        IOFFD=IOFFD+NOCC*NVRT
        IOFFF=IOFFF+NOCC*NOCC*IINTFP
C
100    CONTINUE 
C
C    CALCULATE SECOND CONTRIBUTION TO I(I,A)
C
C    I(I,A) <= - SUM E f(E,I) D(A,E)
C
       IOFFI=I0VO(ISPIN)
       IOFFD=I0VV(ISPIN)
       IOFFF=I0FAI(ISPIN)
C
       DO 200 IRREP=1,NIRREP
C
        NOCC=POP(IRREP,ISPIN)
        NVRT=VRT(IRREP,ISPIN)
C
        IF(MIN(NOCC,NVRT).NE.0) THEN

#ifdef _DEBUG_LVLM
        Print*, "Dvv"
        Call output(Dvv(ioffd),1,nvrt,1,nvrt,nvrt,nvrt,1)
        Print*, "fvo"
        Call output(icore(iofff),1,nvrt,1,nvrt,nvrt,nvrt,1)
#endif 
C
         CALL XGEMM('N','N',NVRT,NOCC,NVRT,ONEM,DVV(IOFFD),NVRT,
     &              ICORE(IOFFF),NVRT,ONE,AIOV(IOFFI),NVRT)
        ENDIF

#ifdef _DEBUG_LVL0
        Print*, "Aiov-2"
        Call output(Aiov(ioffi),1,nocc,1,nvrt,nocc,nvrt,1)
#endif 
C
C  UPDATE OFFSETS
C
        IOFFI=IOFFI+NOCC*NVRT
        IOFFD=IOFFD+NVRT*NVRT
        IOFFF=IOFFF+NOCC*NVRT*IINTFP
C
200    CONTINUE 
C
C    CALCULATE FIRST CONTRIBUTION TO - I(A,I)
C
C    - I(A,I) <= - SUM E f(A,E) D(E,I)
C
       IOFFI=I0VO(ISPIN)
       IOFFD=I0VO(ISPIN)
       IOFFF=I0FAE(ISPIN)
C
       DO 300 IRREP=1,NIRREP
C
        NOCC=POP(IRREP,ISPIN)
        NVRT=VRT(IRREP,ISPIN)
C
        IF(MIN(NOCC,NVRT).NE.0) THEN
C
         CALL XGEMM('N','N',NVRT,NOCC,NVRT,ONE,ICORE(IOFFF),NVRT,
     &              DVO(IOFFD),NVRT,ONE,XAI(IOFFI),NVRT)
        ENDIF
#ifdef _DEBUG_LVL0
       Print*,"Xai-1"
       Call output(Xai(ioffi),1,nvrt,1,nocc,nvrt,nocc,1)
#endif 
C
C  UPDATE OFFSETS
C
        IOFFI=IOFFI+NOCC*NVRT
        IOFFD=IOFFD+NOCC*NVRT
        IOFFF=IOFFF+NVRT*NVRT*IINTFP
C
300    CONTINUE 
C
C    CALCULATE SECOND CONTRIBUTION TO - I(A,I)
C
C    - I(A,I) <=  SUM E f(A,M) D(M,I)
C
       IOFFI=I0VO(ISPIN)
       IOFFD=I0OO(ISPIN)
       IOFFF=I0FAI(ISPIN)
C
       DO 400 IRREP=1,NIRREP
C
        NOCC=POP(IRREP,ISPIN)
        NVRT=VRT(IRREP,ISPIN)
C
        IF(MIN(NOCC,NVRT).NE.0) THEN
C
         CALL XGEMM('N','N',NVRT,NOCC,NOCC,ONE,ICORE(IOFFF),NVRT,
     &              DOO(IOFFD),NOCC,ONE,XAI(IOFFI),NVRT)
        ENDIF
#ifdef _DEBUG_LVL0
       Print*,"Xai-2"
       Call output(Xai(ioffi),1,nvrt,1,nocc,nvrt,nocc,1)
#endif
C
C  UPDATE OFFSETS
C
        IOFFI=IOFFI+NOCC*NVRT
        IOFFD=IOFFD+NOCC*NOCC
        IOFFF=IOFFF+NVRT*NOCC*IINTFP
C
400    CONTINUE 

       CALL SAXPY(NT(ISPIN),ONE,ICORE(I0FAI(ISPIN)),1,
     &            XAI(I0VO(ISPIN)),1)
#ifdef _DEBUG_LVL0
       Print*,"Xai-3"
       Call output(Xai(I0VO(ISPIN)),1,nvrt,1,nocc,nvrt,nocc,1)
#endif

C    ADD HERE -f(I,J) TO I(I,J)
C
       CALL SAXPY(NF1(ISPIN),ONEM,ICORE(I0FMI(ISPIN)),1,
     &            AIOO(I0OO(ISPIN)),1)
C
C    CALCULATE FIRST CONTRIBUTION TO  I(I,J)
C
C     I(I,J) <= - SUM E f(I,E) D(E,J)
C
       IOFFI=I0OO(ISPIN)
       IOFFD=I0VO(ISPIN)
       IOFFF=I0FAI(ISPIN)
C
       DO 500 IRREP=1,NIRREP
C
        NOCC=POP(IRREP,ISPIN)
        NVRT=VRT(IRREP,ISPIN)
C
        IF(MIN(NOCC,NVRT).NE.0) THEN
C
         CALL XGEMM('T','N',NOCC,NOCC,NVRT,ONEM,ICORE(IOFFF),NVRT,
     &              DVO(IOFFD),NVRT,ONE,AIOO(IOFFI),NOCC)
        ENDIF
C
C  UPDATE OFFSETS
C
        IOFFI=IOFFI+NOCC*NOCC
        IOFFD=IOFFD+NOCC*NVRT
        IOFFF=IOFFF+NVRT*NOCC*IINTFP
C
500    CONTINUE 
C
C    CALCULATE SECOND CONTRIBUTION TO  I(I,J)
C
C     I(I,J) <= -  SUM M f(I,M) D(M,J)
C
       IOFFI=I0OO(ISPIN)
       IOFFD=I0OO(ISPIN)
       IOFFF=I0FMI(ISPIN)
C
       DO 600 IRREP=1,NIRREP
C
        NOCC=POP(IRREP,ISPIN)
C
        IF(NOCC.NE.0) THEN
C
         CALL XGEMM('N','N',NOCC,NOCC,NOCC,ONEM,ICORE(IOFFF),NOCC,
     &              DOO(IOFFD),NOCC,ONE,AIOO(IOFFI),NOCC)
        ENDIF
C
C  UPDATE OFFSETS
C
        IOFFI=IOFFI+NOCC*NOCC
        IOFFD=IOFFD+NOCC*NOCC
        IOFFF=IOFFF+NOCC*NOCC*IINTFP
C
600    CONTINUE 

#ifdef _DEBUG_LVL0
       Print*,"IOO"
       Call output(AIOO(1),1,nocc,1,nocc,nocc,nocc,1)
#endif 
C
C    CALCULATE FIRST CONTRIBUTION TO  I(A,B)
C
C     I(A,B) <= - SUM E f(A,E) D(E,B)
C
       IOFFI=I0VV(ISPIN)
       IOFFD=I0VV(ISPIN)
       IOFFF=I0FAE(ISPIN)
C
       DO 700 IRREP=1,NIRREP
C
        NVRT=VRT(IRREP,ISPIN)
C
        IF(NVRT.NE.0) THEN
C
         CALL XGEMM('N','N',NVRT,NVRT,NVRT,ONEM,ICORE(IOFFF),NVRT,
     &              DVV(IOFFD),NVRT,ONE,AIVV(IOFFI),NVRT)
        ENDIF
C
C  UPDATE OFFSETS
C
        IOFFI=IOFFI+NVRT*NVRT
        IOFFD=IOFFD+NVRT*NVRT
        IOFFF=IOFFF+NVRT*NVRT*IINTFP
C
700    CONTINUE 
C
C    CALCULATE SECOND CONTRIBUTION TO  I(A,B)
C
C     I(A,B) <= -  SUM M f(A,M) D(M,B)
C
       IOFFI=I0VV(ISPIN)
       IOFFD=I0VO(ISPIN)
       IOFFF=I0FAI(ISPIN)
C
       DO 800 IRREP=1,NIRREP
C
        NOCC=POP(IRREP,ISPIN)
        NVRT=VRT(IRREP,ISPIN)
C
        IF(MIN(NOCC,NVRT).NE.0) THEN
C
         CALL XGEMM('N','T',NVRT,NVRT,NOCC,ONEM,ICORE(IOFFF),NVRT,
     &              DVO(IOFFD),NVRT,ONE,AIVV(IOFFI),NVRT)
        ENDIF
C
C  UPDATE OFFSETS
C
        IOFFI=IOFFI+NVRT*NVRT
        IOFFD=IOFFD+NOCC*NVRT
        IOFFF=IOFFF+NOCC*NVRT*IINTFP
C
800    CONTINUE 
#ifdef _DEBUG_LVL0
       Print*,"IVV"
       Call output(AIvv(1),1,nvrt,1,nvrt,nvrt,nvrt,1)
#endif 
1000  CONTINUE
C
      ELSE
#ifdef _DEBUG_LVL0
      Write(6,*) "---------HF-----------"
#endif 
C
C  GET OFFSETS FOR HF CASE 
C
C  NOTE I0FMI IS NOW THE OFFSET FOR THE DIAGONAL ELEMENTS
       I0FMI(1)=MXCOR+1-IINTFP*(NOCCO(1)+NVRTO(1))
       I0OO(1)=1
       I0VV(1)=1
       IF(IUHF.EQ.0) THEN
        I0FMI(2)=I0FMI(1)
        I0OO(2)=1
        I0VV(2)=1
       ELSE
        I0FMI(2)=I0FMI(1)-IINTFP*(NOCCO(2)+NVRTO(2))
        I0OO(2)=1+NF1(1)*IUHF
        I0VV(2)=1+NF2(1)*IUHF
       ENDIF
      ENDIF

C      Call Dzero(Aioo,nf1(1))
C      Call Dzero(Aivv,nf2(1))
C      Call Dzero(aiov,Nt(1)) 
C      Call Dzero(xai,Nt(1)) 

#ifdef _DEBUG_LVL0
       Print*, "-------At Entry (should be zero)-----------"
       Print*, "Aiov"
       Call output(Aiov(1),1,nocc,1,nvrt,nocc,nvrt,1)
       Print*,"Xai"
       Call output(Xai(1),1,nvrt,1,nocc,nvrt,nocc,1)
       Print*,"Aioo"
       Call output(AIoo(1),1,nocc,1,nocc,nocc,nocc,1)
       Print*,"AIVV"
       Call output(AIvv(1),1,nvrt,1,nvrt,nvrt,nvrt,1)
#endif 
C
C  DEAL NOW WITH THE DIAGONAL PARTS OF F
C
      CALL GETREC(20,'JOBARC','SCFEVALA',IINTFP*(NOCCO(1)+NVRTO(1)),
     &            ICORE(I0FMI(1)))
      nbas=(NOCCO(1)+NVRTO(1))
      IF(IUHF.EQ.1) THEN
       CALL GETREC(20,'JOBARC','SCFEVALB',IINTFP*(NOCCO(2)+NVRTO(2)),
     &             ICORE(I0FMI(2)))
      ENDIF
C
C  FIRST CONSIDER THE NON-HF TERMS
C 
      DO 2000 ISPIN=1,IUHF+1
C
C  NEGATE THE OCCUPIED DIAGONAL PARTS OF F
C
      CALL VMINUS(ICORE(I0FMI(ISPIN)),NOCCO(ISPIN))
C
      IF(QRHF.OR.ROHF .OR. NONHF_TERMS_EXIST) THEN
C
       INDOCC=-1
       INDVRT=NOCCO(ISPIN)-1
       IOFFVO=I0VO(ISPIN)
       IOFFOO=I0OO(ISPIN)
C
       DO 900 IRREP=1,NIRREP
C
        NOCC=POP(IRREP,ISPIN)
        NVRT=VRT(IRREP,ISPIN)
C
        DO 910 I=1,NOCC
         INDOCC=INDOCC+1
         CALL SAXPY(NVRT,ICORE(I0FMI(ISPIN)+IINTFP*INDOCC),
     &              DVO(IOFFVO+(I-1)*NVRT),1,AIOV(IOFFVO+(I-1)*NVRT),1)
910     CONTINUE
C
        DO 920 IA=1,NVRT
         INDVRT=INDVRT+1
         CALL SAXPY(NOCC,ICORE(I0FMI(ISPIN)+IINTFP*INDVRT),
     &              DVO(IOFFVO+IA-1),NVRT,XAI(IOFFVO+IA-1),NVRT)
920     CONTINUE
C
        IOFFVO=IOFFVO+NOCC*NVRT  
        IOFFOO=IOFFOO+NOCC*NOCC    
900    CONTINUE
C
      ENDIF
#ifdef _DEBUG_LVL0
       Print*, "-------After adding first block-----------"
       Print*, "Aiov"
       Call output(Aiov(1),1,nocc,1,nvrt,nocc,nvrt,1)
       Print*,"Xai"
       Call output(Xai(1),1,nvrt,1,nocc,nvrt,nocc,1)
       Print*,"Aioo"
       Call output(AIoo(1),1,nocc,1,nocc,nocc,nocc,1)
       Print*,"AIVV"
       Call output(AIvv(1),1,nvrt,1,nvrt,nvrt,nvrt,1)
#endif 
C
C
C  NEGATE THE DIAGONAL ELEMENTS OF THE VIRTUAL-VIRTUAl BLOCK OF F
C
      CALL VMINUS(ICORE(I0FMI(ISPIN)+IINTFP*NOCCO(ISPIN)),NVRTO(ISPIN)) 
C
C  NOW THE HF CONTRIBUTIONS
C
      INDOCC=-1
      INDVRT=NOCCO(ISPIN)-1
      IOFFOO=I0OO(ISPIN)
      IOFFVV=I0VV(ISPIN)
C
       DO 950 IRREP=1,NIRREP
C
        NOCC=POP(IRREP,ISPIN)
        NVRT=VRT(IRREP,ISPIN)
C
C  SUBTRACT HERE ALSO THE DIAGONAL ELEMENTS OF F FROM I(I,J)
C  (THIS IS ACTUALLY A CONTRIBUTION OF THE REFERENCE FUNCTION)
C
        CALL SAXPY(NOCC,ONE,ICORE(I0FMI(ISPIN)+IINTFP*(INDOCC+1)),1,
     &             AIOO(IOFFOO),NOCC+1)
C
        DO 960 I=1,NOCC
         INDOCC=INDOCC+1
         CALL SAXPY(NOCC,ICORE(I0FMI(ISPIN)+IINTFP*INDOCC),
     &              DOO(IOFFOO+I-1),NOCC,AIOO(IOFFOO+I-1),NOCC)
960     CONTINUE
C
        DO 970 IA=1,NVRT
         INDVRT=INDVRT+1
         CALL SAXPY(NVRT,ICORE(I0FMI(ISPIN)+IINTFP*INDVRT),
     &              DVV(IOFFVV+IA-1),NVRT,AIVV(IOFFVV+IA-1),NVRT)
970     CONTINUE
C
        IOFFVV=IOFFVV+NVRT*NVRT      
        IOFFOO=IOFFOO+NOCC*NOCC
950    CONTINUE

#ifdef _DEBUG_LVL0
       Print*, "-------After adding the second block-----------"
       Print*, "Aiov"
       Call output(Aiov(1),1,nocc,1,nvrt,nocc,nvrt,1)
       Print*,"Xai"
       Call output(Xai(1),1,nvrt,1,nocc,nvrt,nocc,1)
       Print*,"Aioo"
       Call output(AIoo(1),1,nocc,1,nocc,nocc,nocc,1)
       Print*,"AIVV"
       Call output(AIvv(1),1,nvrt,1,nvrt,nvrt,nvrt,1)
#endif 
C
2000  CONTINUE
C
C  ALL DONE , RETURN
C 
      RETURN
      END

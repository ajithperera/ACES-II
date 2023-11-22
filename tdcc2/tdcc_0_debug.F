C
      SUBROUTINE OPRBAR0(MIJ,MAB,MAI,ICORE,MAXCOR,IUHF,IRREPX)
C
C THIS ROUTINE CALCULATES THE FOLLOWING QUANTITIES FOR THE GENERAL
C ONE-PARTICLE OPERATOR, M
C
C       <0|(1+LAMBDA)exp(-T) M exp(T)|0>  [GROUND STATE EXPECTATION VALUE]
C
C       <0|exp(-T) M exp(T)|0>            
C
      IMPLICIT INTEGER (A-Z)
      DOUBLE PRECISION ONE,ONEM,ZILCH,SNRM2,TWO,SDOT,FACT
      DOUBLE PRECISION MIJ(*),MAB(*),MAI(*),MEXPECT,MBAR00
      DIMENSION ICORE(MAXCOR)
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMINF/NSTART,NIRREP,IRREPS(255,2),DIRPRD(8,8)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NFMI(2),NFEA(2)
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/SYMLOC/ISYMOFF(8,8,25)
      COMMON/FLAGS/ IFLAGS(100)
C
      DATA ONE,ONEM,ZILCH,TWO /1.0D0,-1.0D0,0.0D0,2.0D0/
C
      Mu_0  = ZILCH
      Mut_0 = ZILCH
C
      IF (IRREPX .EQ. 1) THEN
C
         IOFFOO = 1
         IOFFVO = 1
         IOFFVV = 1
C
         DO 10 ISPIN=1,1+IUHF
C
            I000=1
            CALL GETLST(ICORE(I000),1,1,1,ISPIN,90)
C
C Virtual-occupied block
C   _
C   M = Sum M(i,a)*t(i,a)  = -Sum M(a,i)*t(a,i) 
C       a,b                   a,b

            Mu_0 = Mu_0 +SDOT(NT(ISPIN),ICORE(I000),1,
     &                          MAI(IOFFVO),1)
C
C Occupied-virtual block
C   _
C   M = Sum M(a,i)*L(i,a)  = Sum M(a,i)*L(i,a) 
C       a,b                  a,b
C     
            CALL GETLST(ICORE(I000),1,1,1,ISPIN,190)
            Mut_0 =Mut_0 +SDOT(NT(ISPIN),ICORE(I000),1,
     &                         MAI(IOFFVO),1)
C     
            IOFFOO = IOFFOO+NFMI(1)
            IOFFVV = IOFFVV+NFEA(1)
            IOFFVO = IOFFVO+NT(1)
C
 10      CONTINUE
C
         FACT = DFLOAT(2-IUHF)
C
         Mu_0  = Mu_0*FACT
         Mut_0 = Mut_0*FACT
C
      ENDIF
C
      RETURN
      END

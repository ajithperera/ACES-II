C
      SUBROUTINE TDEE_FORM_MUBAR_00(MIJ,MAB,MAI,ICORE,MAXCOR,IUHF,
     &                              IRREPX,MEXPECT,MBAR00)
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
      COMMON/LISTDENS/LDENS
C
      DATA ONE,ONEM,ZILCH,TWO /1.0D0,-1.0D0,0.0D0,2.0D0/
C
      MEXPECT = ZILCH
      MBAR00  = ZILCH
      IANTI   = 0
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
C                         _
C Occupied-occupied block M = Sum M(i,j)*D(j,i)
C                             i,j 
            CALL GETLST(ICORE(I000),1,1,1,ISPIN,ldens)
            MEXPECT = MEXPECT+SDOT(NFMI(ISPIN),ICORE(I000),1,
     &                MIJ(IOFFOO),1)
C                       _
C Virtual-virtual block M = Sum M(a,b)*D(b,a)
C                             a,b 
            CALL GETLST(ICORE(I000),1,1,1,2+ISPIN,ldens)
            MEXPECT = MEXPECT+SDOT(NFEA(ISPIN),ICORE(I000),1,
     &                MAB(IOFFVV),1)
C
            CALL GETLST(ICORE(I000),1,1,1,4+ISPIN,ldens)
C
C Virtual-occupied block
C   _
C   M = Sum M(i,a)*D(a,i)  = -Sum M(a,i)*D(a,i) 
C       a,b                   a,b
C   _
C   M = -Sum M(a,i)*D(a,i) when the perturbation is antisymmetric.
C        a,b   
C           
            IF (IANTI .EQ. 1) THEN
               MEXPECT = MEXPECT-SDOT(NT(ISPIN),ICORE(I000),1,
     &                   MAI(IOFFVO),1)
            ELSE
               MEXPECT = MEXPECT+SDOT(NT(ISPIN),ICORE(I000),1,
     &                   MAI(IOFFVO),1)
            ENDIF
C
            CALL GETLST(ICORE(I000),1,1,1,ISPIN,90)
C
C Virtual-occupied block
C   _
C   M = Sum M(i,a)*t(i,a)  = -Sum M(a,i)*t(a,i) 
C       a,b                   a,b
C   _
C   M = -Sum M(a,i)*t(a,i) when the perturbation is antisymmetric.
C        a,b   
            IF (IANTI .EQ. 1) THEN
               MBAR00 = MBAR00-SDOT(NT(ISPIN),ICORE(I000),1,
     &                  MAI(IOFFVO),1)
               MEXPECT = MEXPECT-SDOT(NT(ISPIN),ICORE(I000),1,
     &                   MAI(IOFFVO),1)
            ELSE
               MBAR00 = MBAR00+SDOT(NT(ISPIN),ICORE(I000),1,
     &                  MAI(IOFFVO),1)
               MEXPECT = MEXPECT+SDOT(NT(ISPIN),ICORE(I000),1,
     &                   MAI(IOFFVO),1)
            ENDIF
C
C Occupied-virtual block
C   _
C   M = Sum M(a,i)*L(i,a)  = Sum M(a,i)*L(i,a) 
C       a,b                  a,b
C     
            CALL GETLST(ICORE(I000),1,1,1,ISPIN,190)
            MEXPECT = MEXPECT+SDOT(NT(ISPIN),ICORE(I000),1,
     &                MAI(IOFFVO),1)
C     
            IOFFOO = IOFFOO+NFMI(1)
            IOFFVV = IOFFVV+NFEA(1)
            IOFFVO = IOFFVO+NT(1)
C
 10      CONTINUE
C
         FACT = DFLOAT(2-IUHF)
C
         MEXPECT = MEXPECT*FACT
         MBAR00  = MBAR00*FACT
C
         IF (IFLAGS(1) .GE. 20) THEN
            WRITE(6,1000)MEXPECT
 1000       FORMAT(T3,'Ground state expectation value of operator ',
     &             F15.8,' a.u.')
            WRITE(6,1001) MBAR00
 1001       FORMAT(T3,'<0|exp(-T) M exp(T)|0> = ',F15.8,' a.u.')
         ENDIF
C
      ENDIF

#ifdef _DEBUG_LVL0
        Write(6,*)
        Write(6,"(a,2F12.6)") "The mu_0/mu_bar alpha+Beta",
     +                   mexpect,mbar00
#endif
C
      RETURN
      END

      SUBROUTINE INTGRT_ONEDEN(NCNTR,NATOMS,NFCT,NTANGM,IANGX,IANGY,
     $   IANGZ,TOTWT,XCDNT,YCDNT,ZCDNT,RSQRD,ALPHA,NANGMOM,NMOMFCT,
     $   FCTDEN,FCTGRX,FCTGRY,FCTGRZ,IPOPF,IUATMS,IFACT,IUHF,NRADPT,
     $   TOTELE,XALPHA,BECKE,dlda,
     &   ELYP,CUTOFF,NAPTS,NCNTO,NMOA,NMOB,
     $   PCOEFFA,PCOEFFB,VALMO,XGVAL,YGVAL,ZGVAL,ANGXVL,ANGYVL,
     $   ANGZVL,MAXANG,RRR,PGRP)
C     
C     This routine does a numerical integration over spherical 
C     coordinates.  The radial quadrature formula is taken from P. 
C     M. W. Gill, B. G. Johnson, J. A. Pople and M. J. Frisch, 
C     Chem. Phys. Lett. 197, 499 (1992).
C     The angular quadrature formula is taken from V. I. Lebedev, 
C     Zh. Vychisl. Mat. Fiz. 15, 48 (1975) and 16, 293 (1976), 
C     (english translation in U.S.S.R. Comput. Math and Math Phys).
C     
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     
      CHARACTER*8 PGRP
C
      COMMON /PAR/ PI,PIX4
      COMMON /IPAR/ LUOUT
      COMMON /CNST/ ALPH,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,X12,X13,X14,X15
      COMMON /FRAC/ F1,F2,F3,F4,F5,F6,F7,F8,F9,F10,F11,F12,F13,F14,F15
C     
      DIMENSION IANGX(56),IANGY(56),IANGZ(56),
     &   TOTWT(NRADPT,194),IPOPF(IUATMS),
     &   XCDNT(NATOMS,NRADPT,194),YCDNT(NATOMS,NRADPT,194),
     &   ZCDNT(NATOMS,NRADPT,194),RSQRD(NATOMS,NRADPT,194),
     &   ALPHA(NFCT),NANGMOM(NATOMS),NMOMFCT(NATOMS*NTANGM),
     &   FCTDEN(NFCT),FCTGRX(NFCT),FCTGRY(NFCT),FCTGRZ(NFCT),
     $   PCOEFFA(NMOA*NFCT),PCOEFFB(NMOB*NFCT),VALMO(NMOA),
     $   XGVAL(NMOA),YGVAL(NMOA),ZGVAL(NMOA),ANGXVL(8),
     $   ANGYVL(8),ANGZVL(8),RRR(NATOMS,NRADPT,194)
C
C     EULER-MACLAURIN SCHEME, RADIAL INTEGRATION
C     GAUSS-MARKOV QUADRATURE, ANGULAR INTEGRATION
C     
C     Symmetry multiplication factor
      FACT=REAL(IFACT)
C     
C     Initialize integral values
C
      TOTELE=0.0D+00
      XALPHA=0.0D+00
      BECKE=0.0D+00
      dlda=0.0d+00
      ELYP=0.0D+00
C
C     Functional value
      FUNCT=0.D+00
      CUTOFF=1.0D-10

#ifdef _DEBUG_LVLM1
        Write(6,*)
        Write(6,"(a)") "The weight array in one-den"
        Write(6,"(4F15.5)") ((TOTWT(j, i),j=1,NRADPT), i=1,194)
#endif

C     
C     Loop over angular points
C
      DO 10 IANGPT=1,NAPTS
     
C     Loop over radial points
         DO 20 IRADPT=1,NRADPT
C
C     Set function counter
            IFCT=1
C     Loop over atoms
            DO 30 IATM=1,NATOMS
C
C     ANGULAR MOMENTUM
               ANGXVL(1)=1.D+00
               ANGYVL(1)=1.D+00
               ANGZVL(1)=1.D+00
               ANGXVL(2)=1.D+00
               ANGYVL(2)=1.D+00
               ANGZVL(2)=1.D+00
               ANGXVL(3)=XCDNT(IATM,IRADPT,IANGPT)
               ANGYVL(3)=YCDNT(IATM,IRADPT,IANGPT)
               ANGZVL(3)=ZCDNT(IATM,IRADPT,IANGPT)
#ifdef _DEBUG_LVLM
              If(IRADPT .le. 10 .and. IANGPT .le. 10)
     &        Write(6,"(3(1x,F15.7))")  XCDNT(IATM, IRADPT, IANGPT),
     &                            YCDNT(IATM, IRADPT, IANGPT),
     &                            ZCDNT(IATM, IRADPT, IANGPT)
              write(6,"(F15.7)") RSQRD(IATM,IRADPT,IANGPT)
              write(6,*) "one_den:END"
#endif 
C
               DO 35 IANG=3,MAXANG+2
                  ANGXVL(IANG+1)=ANGXVL(IANG)*XCDNT(IATM,IRADPT,IANGPT)
                  ANGYVL(IANG+1)=ANGYVL(IANG)*YCDNT(IATM,IRADPT,IANGPT)
                  ANGZVL(IANG+1)=ANGZVL(IANG)*ZCDNT(IATM,IRADPT,IANGPT)
  35           CONTINUE
C
C     Loop over angular momentum for each atom 
               DO 40 IANGMOM=1,NANGMOM(IATM) 
C     Angular intermediates for the density and gradient
C     
C     Density
                  ANGDEN=ANGXVL(IANGX(IANGMOM)+2)
     $               *ANGYVL(IANGY(IANGMOM)+2)
     $               *ANGZVL(IANGZ(IANGMOM)+2)
C     Gradient minus one component
                  ANGXM=IANGX(IANGMOM)
     $               *ANGXVL(IANGX(IANGMOM)+1)
     $               *ANGYVL(IANGY(IANGMOM)+2)
     $               *ANGZVL(IANGZ(IANGMOM)+2)
                  ANGYM=IANGY(IANGMOM)
     $               *ANGXVL(IANGX(IANGMOM)+2)
     $               *ANGYVL(IANGY(IANGMOM)+1)
     $               *ANGZVL(IANGZ(IANGMOM)+2)
                  ANGZM=IANGZ(IANGMOM)
     $               *ANGXVL(IANGX(IANGMOM)+2)
     $               *ANGYVL(IANGY(IANGMOM)+2)
     $               *ANGZVL(IANGZ(IANGMOM)+1)
C     Gradient plus one component
                  ANGXP=ANGXVL(IANGX(IANGMOM)+3)
     $               *ANGYVL(IANGY(IANGMOM)+2)
     $               *ANGZVL(IANGZ(IANGMOM)+2)
                  ANGYP=ANGXVL(IANGX(IANGMOM)+2)
     $               *ANGYVL(IANGY(IANGMOM)+3)
     $               *ANGZVL(IANGZ(IANGMOM)+2)
                  ANGZP=ANGXVL(IANGX(IANGMOM)+2)
     $               *ANGYVL(IANGY(IANGMOM)+2)
     $               *ANGZVL(IANGZ(IANGMOM)+3)
C     
C     Loop over number of functions for each angular momentum
                  DO 50 IMOMFCT=1,NMOMFCT((IATM-1)*NTANGM+IANGMOM)
                     TEXP=DEXP(-ALPHA(IFCT)*RSQRD(IATM,IRADPT,IANGPT))
                     TWOALP=2.D+00*ALPHA(IFCT)
                     FCTDEN(IFCT)=ANGDEN*TEXP
                     FCTGRX(IFCT)=(ANGXM-TWOALP*ANGXP)*TEXP
                     FCTGRY(IFCT)=(ANGYM-TWOALP*ANGYP)*TEXP
                     FCTGRZ(IFCT)=(ANGZM-TWOALP*ANGZP)*TEXP
                     IFCT=IFCT+1
 50               CONTINUE
 40            CONTINUE
 30         CONTINUE
C
C     Alpha spin
C
C     Contract functions with alpha coefficients and add to get 
C     alpha MO values at this point.
            CALL XGEMM('T','N',NMOA,1,NFCT,1.0D+00,PCOEFFA,NFCT,
     $         FCTDEN,NFCT,0.0D+00,VALMO,NMOA)
C     Form alpha density at this point
            CALL XGEMM('T','N',1,1,NMOA,1.0D+00,VALMO,NMOA,VALMO,
     $         NMOA,0.0D+00,ROA,1)
C
C     Contract gradients with alpha coefficients and add to get 
C     alpha MO gradient values at this point.
            CALL XGEMM('T','N',NMOA,1,NFCT,1.0D+00,PCOEFFA,NFCT,
     $         FCTGRX,NFCT,0.0D+00,XGVAL,NMOA)
            CALL XGEMM('T','N',NMOA,1,NFCT,1.0D+00,PCOEFFA,NFCT,
     $         FCTGRY,NFCT,0.0D+00,YGVAL,NMOA)
            CALL XGEMM('T','N',NMOA,1,NFCT,1.0D+00,PCOEFFA,NFCT,
     $         FCTGRZ,NFCT,0.0D+00,ZGVAL,NMOA)
C     Form alpha x, y, and z gradient components at this point
            CALL XGEMM('T','N',1,1,NMOA,2.0D+00,VALMO,NMOA,XGVAL,
     $         NMOA,0.0D+00,GRADXA,1)
            CALL XGEMM('T','N',1,1,NMOA,2.0D+00,VALMO,NMOA,YGVAL,
     $         NMOA,0.0D+00,GRADYA,1)
            CALL XGEMM('T','N',1,1,NMOA,2.0D+00,VALMO,NMOA,ZGVAL,
     $         NMOA,0.0D+00,GRADZA,1)
C     
C     Form alpha gradient invarient (grad dot grad)
            GRDAA=GRADXA*GRADXA+GRADYA*GRADYA+GRADZA*GRADZA
C     
            IF(IUHF.EQ.1)THEN
C
C     Beta Spin
C
C     Contract functions with beta coefficients and add to get beta
C     MO values at this point.
               CALL XGEMM('T','N',NMOB,1,NFCT,1.0D+00,PCOEFFB,NFCT,
     $            FCTDEN,NFCT,0.0D+00,VALMO,NMOB)
C     Form beta density at this point
               CALL XGEMM('T','N',1,1,NMOB,1.0D+00,VALMO,NMOB,
     $            VALMO,NMOB,0.0D+00,ROB,1)
C     
C     Contract gradients with coefficients and add to get MO values at
C     this point.
               CALL XGEMM('T','N',NMOB,1,NFCT,1.0D+00,PCOEFFB,NFCT,
     $            FCTGRX,NFCT,0.0D+00,XGVAL,NMOB)
               CALL XGEMM('T','N',NMOB,1,NFCT,1.0D+00,PCOEFFB,NFCT,
     $            FCTGRY,NFCT,0.0D+00,YGVAL,NMOB)
               CALL XGEMM('T','N',NMOB,1,NFCT,1.0D+00,PCOEFFB,NFCT,
     $            FCTGRZ,NFCT,0.0D+00,ZGVAL,NMOB)
C     Form gradients at this point
               CALL XGEMM('T','N',1,1,NMOB,2.0D+00,VALMO,NMOB,
     $            XGVAL,NMOB,0.0D+00,GRADXB,1)
               CALL XGEMM('T','N',1,1,NMOB,2.0D+00,VALMO,NMOB,
     $            YGVAL,NMOB,0.0D+00,GRADYB,1)
               CALL XGEMM('T','N',1,1,NMOB,2.0D+00,VALMO,NMOB,
     $            ZGVAL,NMOB,0.0D+00,GRADZB,1)
C
C     Form beta gradient invarient
               GRDBB=GRADXB*GRADXB+GRADYB*GRADYB+GRADZB*GRADZB
C     Form alpha-beta gradient invarient
               GRDAB=GRADXA*GRADXB+GRADYA*GRADYB+GRADZA*GRADZB
C
            ELSE
               ROB=ROA
               GRDAB=GRDAA
               GRDBB=GRDAA
            ENDIF
C
C      Determine the value of the functionals
C
            XALPHA=XALPHA+TOTWT(IRADPT,IANGPT)*
     $         XALPHF(ROA,ROB,GRDAA,GRDAB,GRDBB)
            BECKE=BECKE+TOTWT(IRADPT,IANGPT)*
     $         BECKEF(ROA,ROB,GRDAA,GRDAB,GRDBB)
            dlda=dlda+totwt(iradpt,iangpt)*
     $         dldaf(roa,rob,grdaa,grdab,grdbb)
            ELYP=ELYP+TOTWT(IRADPT,IANGPT)*
     $         ELYPF(ROA,ROB,GRDAA,GRDAB,GRDBB)
C
C      Determine total density
C
            TOTELE=TOTELE+TOTWT(IRADPT,IANGPT)*(ROA+ROB)
C
C     If the radial density has fallen below the density threshold skip the 
C     remaining radial points and go to the next angular point

CSSS            IF (IANGPT .EQ. 1) Write(6,*) "The density =", (ROA+ROB)
           
            IF(ROA+ROB.LT.CUTOFF) GOTO 10
C
C        Next radial point
 20      CONTINUE
C
C     Next angular point
 10   CONTINUE
C     
C     Determine the number of equivalent atoms for this integration center
      IF(PGRP.NE.'C1')THEN
        FACT=FACT*REAL(IPOPF(NCNTO))
      ENDIF
C
C     Determine total electrons for this integration region
      TOTELE=FACT*TOTELE
C 
C     XALPHA exchange functional value for this integration region
      XALPHA=FACT*XALPHA
C
C     Becke exchange functional value for this integration region
      BECKE=FACT*BECKE
      dlda=fact*dlda
C
C     LYP correlation functional value for this integration region
      ELYP=FACT*ELYP
C     
      RETURN
      END

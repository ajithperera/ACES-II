      SUBROUTINE HDIVIDEDD(SCR, MAXCOR, IUHF, IRREPX, NSIZEC, ROOT,
     &   IPOWER)
C
C  THE DIAGONAL OF (ROOT -HBAR^-1)**IPOWER  *  DOUBLE COMPONENTS  OF INCREMENT VECTOR
c IS CALCULATED, AND PUT IN PLACE OF THE INPUT VECTOR
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL SHIFTROOT
      DIMENSION SCR(MAXCOR)
C
      COMMON /SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),ID(18)
      COMMON/LISTPROJ/LISTH0, ICOLPR1, ICOLPR2
      COMMON/PROJECT/IPROJECT, IPATTERN, NCALC, ICALC, IWINDOW(8)
      COMMON/SHFTROOT/SHIFTROOT
C
      ONE = 1.0D0
      thresh = 8.0d-2
C
      I000 = 1
      I010 = I000 + NSIZEC
      I020 = I010 + NSIZEC
C
      N1AA = IRPDPD(IRREPX,16)
      IF (IUHF .EQ. 0) THEN
        N1BB = 0
      ELSE
        N1BB = IRPDPD(IRREPX,17)
      ENDIF
      N1 = N1AA + N1BB
c
      CALL GETLST(SCR(I010),1,1,1,1,472)
c
      if (shiftroot .and. ipower.eq. 1) then
C
C MAKE SURE THAT ROOT DOES NOT LEAD TO SINGULAR DENOMINATOR.
C
C  DETERMINE RANGE AROUND ROOT THAT DOES NOT HAVE SINGULARITIES
C
        ITIMES = 0
  100   ITIMES = ITIMES + 1
        RPLUS = 10.0D0
        RMIN = 10.0D0
        DO I = N1, NSIZEC -1
          DENOM =  SCR(I010+I) - ROOT
          IF (DENOM .GE. 0.0D0) THEN
            RPLUS = DMIN1(RPLUS,DENOM)
          ELSE
            RMIN = DMIN1(RMIN, -DENOM)
          ENDIF
        ENDDO
c
c      write(6,*) ' @hdividedd, rmin, rplus, root'
c      write(6,*) rmin, rplus, root
C
C enclosing singularities located at root-rmin and root+rplus
C
        IF ((RMIN + RPLUS) .GE. 2.1D0*THRESH) THEN
          IF (RMIN .LT. THRESH) THEN
            write(6,*) ' Root is located in a singular region',
     &         root
            ROOT = ROOT-RMIN+1.04d0*THRESH
          ENDIF
          IF (RPLUS .LT. THRESH) THEN
            write(6,*) ' Root is located in a singular region',
     &         root
            ROOT = ROOT+RPLUS-1.04D0*THRESH
          ENDIF
        ELSE
C
C  COMPLICATED ROOT. SHIFT ROOT TO LY BEYOND ROOT + RPLUS
C  AND REPEAT PROCEDURE
C
          write(6,*) ' Root is located in a singular region', root
          if (itimes .ge. 3) then
            write(6,*) ' This result is going to be suspect / garbage'
          else
            ROOT = ROOT + RPLUS + 1.04D0*THRESH
            GO TO 100
          endif
C
        ENDIF
      endif
C          
      DO I = 1, NSIZEC
        SCR(I010+I-1) = ONE / (ROOT-SCR(I010+I-1))**IPOWER
      ENDDO
C
      IF (IPROJECT .GE. 2) THEN
C
C PROJECT OUT UNINTERESTING COMPONENTS.
C
        CALL GETLST(SCR(I000),ICOLPR2,1,1,1,LISTH0)       
        CALL VECPRD(SCR(I000),SCR(I010),SCR(I010),NSIZEC)
      ENDIF
C
      CALL LOADVEC1(IRREPX,SCR,MAXCOR,IUHF,490,0,460,NSIZEC,
     &   .FALSE.)
C
      CALL VECPRD(SCR(I000+N1),SCR(I010+N1),SCR(I000+N1),NSIZEC-N1)
C
C      DO I = N1, NSIZEC - 1
C        SCR(I000 + I) =  SCR(I000+I) / (ROOT - SCR(I010+I))**IPOWER
C      ENDDO
C
C  DO NOT CHANGE THE C1 COEFFICENTS
C
C      CALL GETLST(SCR(I000), 1,1,1,1,490)
C      IF (IUHF .NE. 0) CALL GETLST(SCR(I000+N1AA),1,1,1,2,490)
C
      CALL UPDATES(IRREPX,SCR(I000),444,0,490,IUHF)
C
      RETURN
      END
      

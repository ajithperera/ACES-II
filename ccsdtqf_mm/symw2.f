      SUBROUTINE SYMW2m(nu,W,W2,w3,NROW1,NCOL1,NROW2,NCOL2,ISPIN1,
     1ISPIN2,IOFF1,IOFF2,IRPA,IRPC,IRPD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)



c maxbasfn.par : begin

c MAXBASFN := the maximum number of (Cartesian) basis functions

c This parameter is the same as MXCBF. Do NOT change this without changing
c mxcbf.par as well.

      INTEGER MAXBASFN

      PARAMETER (MAXBASFN=1000)



c maxbasfn.par : end


      INTEGER A,C,D,POP,VRT,absvrt,absocc
c      DIMENSION W2(NROW2,NCOL2),W(NROW1,NCOL1)
      DIMENSION W2(NROW2,NCOL2),W(NROW1,NCOL1),w3(nu,nu,nu)
      COMMON /SYM/    POP(8,2),VRT(8,2),NTAA,NTBB,NF1AA,NF2AA,
     1                NF1BB,NF2BB
      common/actorb/absvrt(maxbasfn,8,2),absocc(maxbasfn,8,2)
c      write(6,*)'entering symw2'
c      write(6,*)'irpa,irpc,irpd:',irpa,irpc,irpd
C
      DO 30     D=1,VRT(IRPD,ISPIN2)
      DO 20     C=1,VRT(IRPC,ISPIN1)
      DO 10     A=1,VRT(IRPA,ISPIN1)
C
        W(IOFF1 + (D-1)*VRT(IRPC,ISPIN1) + C,A) =
     1 W2(IOFF2 + (C-1)*VRT(IRPA,ISPIN1) + A,D)
C
   10 CONTINUE
   20 CONTINUE
   30 CONTINUE
      DO 130 D=1,VRT(IRPD,ISPIN2)
         DO 120 C=1,VRT(IRPC,ISPIN1)
            DO 110 A=1,VRT(IRPA,ISPIN1)
               ia=absvrt(a,irpa,1)
               ic=absvrt(c,irpc,1)
               id=absvrt(d,irpd,2)
               w3(ia,ic,id)=W2(IOFF2+(C-1)*VRT(IRPA,ISPIN1) + A,D)
 110        CONTINUE
 120     CONTINUE
 130  CONTINUE
      RETURN
      END
















































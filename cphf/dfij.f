 
          
      SUBROUTINE DFIJ(IRREPX,IOFFX,NX,DFIJA,
     &                UAIA,FAIA,SCR,
     &                UIJA,FIJA,EVALA)
C
C THIS ROUTINE CALCULATES THE FINAL CONTRIBUTION TO THE
C OCCUPIED-OCCUPIED BLOCK OF THE TOTAL DERIVATIVES OF
C THE FOCK MATRIX
C
C  d f(I,J)         x        x                            x
C  -------- = f(I,J)  + SUM U (R,I)  f(R,J) + SUM f(I,R) U (R,J)
C    d x                 R                     R
C
C                              x                           x
C HOWEVER, THE TERMS SUM -1/2 S (M,I) f(M,J) + SUM f(I,M) -1 /2 S (M,J)
C                     M                    M
C HAVE BEEN ALREADY ADDED (SEE FORMFD).
C
C THE REMAINING TERMS ARE SIMPLY:
C
C                   x                           x
C              SUM U (E,I) f(E,J) + SUM f(I,E) U (E,J)
C               E                    E 
C
C NOTE THAT THIS VERSION USES PERTURBED STANDARD ORBITALS 
C (MEANING ALPHA = BETA), THE ALPHA D F(IJ)/D X HAS THUS TO
C BE MODIFIED IN THE FOLLOWING WAY
C
C  d f(1,j) / d x :    U(M,1)^chi f(M,J) + 1/2 S(M,1)^chi f(M,J)
C
C  d f(I,1) / d x :    f(I,M) U(M,1)^chi + 1/2 f(I,M) S(M,1)^chi
C
C ALL d f(I,J) / d x
C
C   U(M,I)^chi f(M,J) + 1/2 S(M,I)^chi f(M,J)
C
C   + f(I,M) U(M,J)^chi + 1/2 f(I,M) S(M,J)^chi
C
C 
CEND
C
C CODED MAY/91 JG
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER DIRPRD,POP,VRT
      LOGICAL FIELD,GEOM,SPIN,THIRD,MAGN
C
      DIMENSION DFIJA(100),UAIA(100),FAIA(100),SCR(100)
      DIMENSION UIJA(100),FIJA(100),EVALA(100)
C
      COMMON/PERT/NTPERT,NPERT(8),IIPERT(8),IXPERT,IYPERT,IZPERT,
     &            IYZPERT,IXZPERT,IXYPERT,ITRANSX,ITRANSY,ITRANSZ,
     &            NUCIND
      COMMON/PERT2/FIELD,GEOM,THIRD,MAGN,SPIN
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),NJUNK(18)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NF1(2),NF2(2)
      COMMON/SYMINF/NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
      COMMON/BASSYM/NBAS(8),NBASIS,NBASSQ,NBASTT
C
      DATA ONE,HALF,TWOM /1.D0,0.5D0,-2.D0/
C
      CALL GETLST(FIJA,1,1,1,3,91)
C
      DO 10 ISPIN=1,2
C
       CALL GETLST(FAIA,1,1,1,2+ISPIN,93)
C
       DO 12 IPERT=IOFFX+1,NX
        CALL GETLST(UAIA,IPERT,1,1,IRREPX,181+ISPIN)
        CALL ZERO(DFIJA,IRPDPD(IRREPX,20+ISPIN))
C
        IOFFF=1
        IOFFD=1
C
        DO 13 IRREPR=1,NIRREP
C
         IRREPL=DIRPRD(IRREPX,IRREPR)
C
         IOFFU=1
         DO 89 IRREP=1,IRREPL-1
          IOFFU=IOFFU+POP(IRREP,ISPIN)*VRT(DIRPRD(IRREP,IRREPX),ISPIN)
89       CONTINUE
C 
         NOCCR=POP(IRREPR,ISPIN)
         NOCCL=POP(IRREPL,ISPIN)
         NVRTR=VRT(IRREPR,ISPIN)
C
         CALL XGEMM('T','N',NOCCL,NOCCR,NVRTR,ONE,UAIA(IOFFU),
     &              NVRTR,FAIA(IOFFF),NVRTR,ONE,DFIJA(IOFFD),NOCCL)
C
         IOFFF=IOFFF+NVRTR*NOCCR
         IOFFD=IOFFD+NOCCR*NOCCL
C
13      CONTINUE
C
        CALL SYMMET6(IRREPX,POP(1,ISPIN),DFIJA,SCR,
     &               IRPDPD(IRREPX,20+ISPIN))
C
        CALL GETLST(SCR,IPERT,1,1,IRREPX,175+ISPIN)
        CALL SAXPY(IRPDPD(IRREPX,20+ISPIN),ONE,SCR,1,DFIJA,1)
C
        IF(ISPIN.EQ.1) THEN
         CALL GETREC(20,'JOBARC','SCFEVLA0',IINTFP*NBASIS,
     &               EVALA)
         CALL GETLST(UIJA,IPERT,1,1,IRREPX,198)
         IF(IPERT.LE.NPERT(IRREPX).AND.(.NOT.FIELD)) THEN
          CALL GETLST(SCR,IPERT,1,1,IRREPX,170)
          CALL SAXPY(IRPDPD(IRREPX,21),HALF,SCR,1,UIJA,1) 
         ENDIF
         CALL SSCAL(IRPDPD(IRREPX,21),TWOM,UIJA,1)
         CALL SXEVAL1(IRREPX,DFIJA,EVALA,UIJA,POP(1,1),0)
         CALL SXF(IRREPX,DFIJA,FIJA,UIJA,POP(1,1))
        ENDIF
C
        CALL PUTLST(DFIJA,IPERT,1,1,IRREPX,175+ISPIN)
c        call checksum('dfij',scr,irpdpd(irrepx,20+ispin))
C
12     CONTINUE
10    CONTINUE
C
      RETURN
      END

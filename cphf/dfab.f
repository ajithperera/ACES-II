      SUBROUTINE DFAB(IRREPX,IOFFX,NX,DFABA,
     &                UAIA,SAIA,FAIA,SCR,
     &                UABB,FABB,EVALB)
C
C THIS ROUTINE CALCULATES THE FINAL CONTRIBUTION TO THE
C VIRTUAL-VIRTUAL BLOCK OF THE TOTAL DERIVATIVES OF
C THE FOCK MATRIX
C
C  d f(A,B)         x        x                            x
C  -------- = f(A,B)  + SUM U (R,A)  f(R,B) + SUM f(A,R) U (R,B)
C    d x                 R                     R
C
C                                x                               x
C HOWEVER, THE TERMS  - 1/2 SUM S (E,A) f(E,B) - 1/2 SUM f(A,E) S (E,B)
C                            E                        E
C
C HAVE BEEN ALREADY ADDED (SEE FORMFD).
C THE REMAINING TERMS ARE SIMPLY :
C
C                   x                           x
C              SUM U (M,A) f(M,B) + SUM f(A,M) U (M,B)
C               M                    M 
C
C HOWEVER IN ORDER TO USE STANDARD PERTURBED ORBITALS (ALPHA=BETA)
C THE BETA FOCK D f(ab)/d X HAS TO BE MODIFIED
C
C   SUM  (U(e,a) + 1/2 S(e,a)) f(e,b) + SUM f(a,e) (U(e,b) + 1/2 S(e,b))
C    E                                   E
C 
CEND
C
C CODED MAY/91 JG
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      LOGICAL GEOM,FIELD,THIRD,MAGN,SPIN
      INTEGER DIRPRD,POP,VRT
C
      DIMENSION DFABA(100),UAIA(100),FAIA(100),SCR(100)
      DIMENSION UABB(100),FABB(100),EVALB(100) 
C
      COMMON/MACHSP/IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON/PERT/NTPERT,NPERT(8),IIPERT(8),IXPERT,IYPERT,IZPERT,
     &            IYZPERT,IXZPERT,IXYPERT,ITRANSX,ITRANSY,ITRANSZ,
     &            NUCIND
      COMMON/PERT2/FIELD,GEOM,THIRD,MAGN,SPIN
      COMMON/SYMPOP/IRPDPD(8,22),ISYTYP(2,500),NJUNK(18)
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NF1(2),NF2(2)
      COMMON/SYMINF/NSTART,NIRREP,IRREPA(255,2),DIRPRD(8,8)
      COMMON/BASSYM/NBAS(8),NBASIS,NBASSQ,NBASTT
      COMMON/INFO/NOCCO(2),NVRTO(2)
C
      DATA ONE,ONEM,HALF,TWOM /1.D0,-1.D0,0.5D0,-2.D0/
C
      CALL GETLST(FABB,1,1,1,4,92)
C
      DO 10 ISPIN=1,2
C
       CALL GETLST(FAIA,1,1,1,2+ISPIN,93)
C
       DO 12 IPERT=IOFFX+1,NX
        CALL GETLST(UAIA,IPERT,1,1,IRREPX,181+ISPIN)
        CALL VMINUS(UAIA,IRPDPD(IRREPX,8+ISPIN))
        IF((.NOT.FIELD).AND.(IPERT.LE.NPERT(IRREPX))) THEN
         CALL GETLST(SAIA,IPERT,1,1,IRREPX,173+ISPIN) 
         CALL SAXPY(IRPDPD(IRREPX,8+ISPIN),ONEM,SAIA,1,UAIA,1)
        ENDIF
        CALL ZERO(DFABA,IRPDPD(IRREPX,18+ISPIN))
C
        IOFFF=1
        IOFFD=1
        IOFFU=1
C
        DO 13 IRREPR=1,NIRREP
C
         IRREPL=DIRPRD(IRREPX,IRREPR)
C
         NOCCR=POP(IRREPR,ISPIN)
         NVRTL=VRT(IRREPL,ISPIN)
         NVRTR=VRT(IRREPR,ISPIN)
C
         CALL XGEMM('N','T',NVRTL,NVRTR,NOCCR,ONE,UAIA(IOFFU),
     &              NVRTL,FAIA(IOFFF),NVRTR,ONE,DFABA(IOFFD),NVRTL)
C
         IOFFF=IOFFF+NVRTR*NOCCR
         IOFFU=IOFFU+NVRTL*NOCCR
         IOFFD=IOFFD+NVRTR*NVRTL
C
13      CONTINUE
C
        CALL SYMMET6(IRREPX,VRT(1,ISPIN),DFABA,SCR,
     &               IRPDPD(IRREPX,18+ISPIN))
C
        CALL GETLST(SCR,IPERT,1,1,IRREPX,177+ISPIN)
        CALL SAXPY(IRPDPD(IRREPX,18+ISPIN),ONE,SCR,1,DFABA,1)
C
        IF(ISPIN.EQ.2) THEN
         CALL GETREC(20,'JOBARC','SCFEVLB0',IINTFP*NBASIS,EVALB)
         CALL GETLST(UABB,IPERT,1,1,IRREPX,199)
         IF(IPERT.LE.NPERT(IRREPX).AND.(.NOT.FIELD)) THEN
          CALL GETLST(SCR,IPERT,1,1,IRREPX,173)
          CALL SAXPY(IRPDPD(IRREPX,20),HALF,SCR,1,UABB,1)
         ENDIF
         CALL SSCAL(IRPDPD(IRREPX,20),TWOM,UABB,1)
         CALL SXEVAL1(IRREPX,DFABA,EVALB(1+NOCCO(2)),UABB,VRT(1,2),0)
         CALL SXF(IRREPX,DFABA,FABB,UABB,VRT(1,2))
        ENDIF
C
        CALL PUTLST(DFABA,IPERT,1,1,IRREPX,177+ISPIN)
c        call checksum('dfab',scr,irpdpd(irrepx,18+ispin))
C
12     CONTINUE
C
10    CONTINUE
      RETURN
      END

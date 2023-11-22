C
      SUBROUTINE MKDBLAB2(W, T, Z, MAXSIZE, DISSYG, NUMSYG, DISSYT,
     &                    NUMSYT, LISTG, LISTT, IRREPTL, IRREPTR, 
     &                    IRREPGL, IRREPGR, IRREPQL, IRREPQR,
     &                    IRREPX, POP1, POP2, VRT1, VRT2, LISTZ, IUHF)
C
C This routine computes the ABAB hole-hole ladder contribution.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      INTEGER DISSYT,DISSYG,DIRPRD,POP1,POP2,VRT1,VRT2
      INTEGER DISREAD,DISMAX,DISLEFT
      DIMENSION W(DISSYG),T(DISSYT,NUMSYT)
      DIMENSION Z(DISSYT,NUMSYG)
      DIMENSION POP1(8),POP2(8),VRT1(8),VRT2(8)
C
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      COMMON/SYMPOP/ IRPDPD(8,22),ISYTYP(2,500),NTOT(18)
      COMMON/FLAGS/ IFLAGS(100)
      COMMON/FILES/ LUOUT, MOINTS
C   
C Common blocks used in the quadratic term
C 
      COMMON /QADINTI/ INTIMI, INTIAE, INTIME
      COMMON /QADINTG/ INGMNAA, INGMNBB, INGMNAB, INGABAA,
     &                 INGABBB, INGABAB, INGMCAA, INGMCBB,
     &                 INGMCABAB, INGMCBABA, INGMCBAAB,
     &                 INGMCABBA
      COMMON /APERTT1/ IAPRT1AA             
      COMMON /BPERTT1/ IBPRT1AA
      COMMON /APERTT2/ IAPRT2AA1, IAPRT2BB1, IAPRT2AB1, IAPRT2AB2, 
     &                 IAPRT2AB3, IAPRT2AB4, IAPRT2AA2, IAPRT2BB2,
     &                 IAPRT2AB5 
C
      DATA ONE, AZERO /1.0D0, 0.D0/
C
C Pick up the relevent T2 and G pieces
C
      CALL GETLST(T, 1, NUMSYT, 1, IRREPTR, LISTT)
C
C Check if we can do it in-core 
C     
      IF (MAXSIZE .GE. NUMSYG*DISSYG) THEN
C
         CALL ZERO(Z, NUMSYG*DISSYT)
C
         CALL GETLST(W, 1, NUMSYG, 2, IRREPGR, LISTG)
C     
C Compute the ladder contribution
C     
         CALL XGEMM('N', 'N', DISSYT, NUMSYG, DISSYG, ONE, T, DISSYT,
     &               W, DISSYG, AZERO, Z, DISSYT)
C     
      ELSE
C     
C We have to do it out of core 
C     
         CALL ZERO(Z, NUMSYG*DISSYT)
C     
         DISMAX  = MAXSIZE/DISSYG
         IOFFSET = 1
         DISLEFT = NUMSYG
C
 10      CONTINUE
         DISREAD = MIN(DISMAX, DISLEFT)
C
         IF(DISREAD .EQ. 0) THEN
            WRITE(6, *) ' @-MKDBLAB2-F, no out-of-core algorithm.'
            CALL ERREX
         ENDIF
C
         DISLEFT = DISLEFT - DISREAD
C
         CALL GETLST(W, IOFFSET, DISREAD, 2, IRREPGR, LISTG)
C     
C Compute now the ladder contribution
C     
         CALL XGEMM('N', 'N', DISSYT, DISREAD, DISSYG, ONE, T, DISSYT,
     &               W, DISSYG, ONE, Z(1,IOFFSET), DISSYT)
C
         IOFFSET = IOFFSET + DISREAD
         IF(DISLEFT .NE. 0) GO TO 10
C
      ENDIF
C
      IF (IFLAGS(1) .GE. 20) THEN

         NSIZE = NUMSYG*DISSYT
C
         CALL HEADER('Checksum @-MNIJIND per sym. block', 0, LUOUT)
            
         WRITE(LUOUT, *) ' ABAB =  ', SDOT(NSIZE, Z, 1, Z, 1)

      ENDIF
C     
C Augment the T2 increment in the disk.
C     
      CALL SUMSYM2(Z, T, NUMSYG*DISSYT, 1, IRREPQR, LISTZ)
C
      RETURN
      END

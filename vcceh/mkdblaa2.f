C
      SUBROUTINE MKDBLAA2(W, T, Z, MAXSIZE, DISSYG, NUMSYG, DISSYT,
     &                    NUMSYT, LISTG, LISTT, IRREPTL, IRREPTR,
     &                    IRREPGL, IRREPGR, IRREPQL, IRREPQR, 
     &                    IRREPX, POP, VRT, LISTZ, ISPIN)
C
C This routine computes the AAAA or BBBB hole-hole ladder contributions.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DISSYT,DISSYG,DISLEFT,DISREAD,DISMAX,DIRPRD,POP,VRT
      DIMENSION W(DISSYG,1),T(DISSYT,NUMSYT)
      DIMENSION Z(DISSYT,NUMSYG)
      DIMENSION POP(8),VRT(8) 
      CHARACTER*8 SPCASE(2)
C
      COMMON /SYMINF/ NSTART,NIRREP,IRREPA(255),IRREPB(255),
     &                DIRPRD(8,8)
      COMMON /FILES/ LUOUT, MOINTS
      COMMON /FLAGS/ IFLAGS(100)
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
      DATA SPCASE /'AAAA =  ', 'BBBB =  '/
      DATA ONE,AZERO /1.0D0, 0.0D0/
C
C Pick up the relevent T2 and W pieces. Here we do not need to make
C TAU amplitudes.

      CALL GETLST(T, 1, NUMSYT, 1, IRREPTR, LISTT)
C     
C Can we do it in in-core
C     
      IF(MAXSIZE .GE. DISSYG*NUMSYG) THEN
C
         CALL ZERO(Z, NUMSYG*DISSYT)
C     
C Incore
C     
         CALL GETLST(W, 1, NUMSYG, 2, IRREPGR, LISTG)
C     
C Compute the ladder contribution.
C     
         CALL XGEMM('N', 'N', DISSYT, NUMSYG, DISSYG, ONE, T, DISSYT,
     &               W, DISSYG, AZERO, Z, DISSYT)
C
      ELSE
C     
C We have to do it out-of core 
C     
         CALL ZERO(Z, NUMSYG*DISSYT)
C     
C Maximum number of distribution which can be held in core 
C     
         DISMAX  = MAXSIZE/DISSYG
         IOFFSET = 1
         DISLEFT = NUMSYG
C
 10      CONTINUE
         DISREAD = MIN(DISMAX, DISLEFT)
C
         IF(DISREAD .EQ. 0) THEN
            WRITE(6, *) ' @-MKDBLAA2-F, no out-of-core algorithm.'
            CALL ERREX
         ENDIF
C
         DISLEFT = DISLEFT - DISREAD
         CALL GETLST(W, IOFFSET, DISREAD, 2, IRREPGR, LISTG)
C     
C Now compute the ladder contribution
C     
         CALL XGEMM('N', 'N', DISSYT, DISREAD, DISSYG, ONE, T, DISSYT,
     &               W, DISSYG, ONE, Z(1,IOFFSET), DISSYT)
C         
         IOFFSET = IOFFSET + DISREAD
         IF(DISLEFT .NE. 0) GO TO 10
      ENDIF
C     
C Update the doubles contribution in the disk
C     
         IF (IFLAGS(1) .GE. 20) THEN
C
            NSIZE = NUMSYG*DISSYT
            CALL HEADER('Checksum @-GMNIJIND per sym. block', 0, LUOUT)
C            
            WRITE(LUOUT, *) SPCASE(ISPIN), SDOT(NSIZE, Z, 1, Z, 1)
         ENDIF
C
      CALL SUMSYM2(Z, T, NUMSYG*DISSYT, 1, IRREPQR, LISTZ)
C
C      CALL CHKSUM(Z, T, DISSYT*NUMSYG, IRREPQR, LISTZ, SPCASE(ISPIN),
C     &            'MKDBLAA2')
C     
      RETURN
      END

C
      SUBROUTINE MKDBLAB5(T1T1, Q, H, T1A, T1B, POP1, POP2, VRT1, VRT2,
     &                    DISSYH, NUMSYH, DISSYQ, NUMSYQ, DISTMP,
     &                    NUMTMP, LISTH, LISTQ, IRREPHBR, IRREPTR,
     &                    IRREPQR, IRREPX, MAXSIZE, ITYPE, ISPIN,
     &                    IUHF)    
C   
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      INTEGER DISSYH, DISSYQ, DISTMP, VRT1, VRT2, POP1, POP2,
     &        DIRPRD, DISMAX, DISLEFT, DISREAD
      DIMENSION T1A(1), T1B(1), T1T1(DISTMP, NUMTMP),
     &          Q(DISSYQ, NUMSYQ), H(DISSYH, 1)
      DIMENSION POP1(8), POP2(8), VRT1(8), VRT2(8)
C 
      COMMON /SYMINF/ NSTART, NIRREP, IRREPA(255), IRREPB(255),
     &                DIRPRD(8,8)
      COMMON /MACHSP/ IINTLN, IFLTLN, IINTFP, IALONE, IBITWD
      COMMON /FILES/ LUOUT, MOINTS
      COMMON /SYMPOP/ IRPDPD(8,22), ISYTYP(2,500), NTOT(18)
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
      DATA ONE, HALF, HALFM, ZILCH / 1.00D+00, 0.5D+00, -0.5D+00,
     &                              0.00D+00/
C 
C Build the T(A,M)*T(b,n) amplitude. This is only term which contribute
C in ABAB spin case.
C
      CALL ZERO (Q, DISSYQ*NUMSYQ)
C     
      CALL MKT1TAU(T1T1, T1A, T1B, DISTMP, NUMTMP, POP1, POP2, 
     &             VRT1, VRT2, IRREPTR, IRREPX, IRREPX, ISPIN,
     &             ONE)
C
C Check if we can do it in core.
C
      IF (MAXSIZE .GE. DISSYH*NUMSYH) THEN
C     
C Load the Hbar(Mn,Ij) elements to the memory.
C     
         CALL GETLST(H, 1, NUMSYH, 1, IRREPHBR, LISTH)
C
C Carry out the multiplication. The result is ordered as Q(Ab,Ij)      
C 
         IF (ITYPE .EQ. 1) THEN
C
            CALL XGEMM('N', 'N', DISTMP, NUMSYH, DISSYH, ONE, T1T1, 
     &                  DISTMP, H, DISSYH, ZILCH, Q, DISSYQ)
C     
         ELSE IF (ITYPE .EQ. 6) THEN
C     
            CALL XGEMM('N', 'N', DISSYH, NUMTMP, NUMSYH, ONE, H,
     &                  DISSYH, T1T1, DISTMP, ZILCH, Q, DISSYQ)
C
         ENDIF
C
         IF (IUHF .EQ. 0 .AND. ISYTYP(1, 233) .EQ. 5 .AND. ITYPE
     &       .EQ. 6) THEN
C
            CALL SYMEXP4(IRREPQR, VRT1, POP1, DISSYQ, DISSYH, NUMSYQ,
     &                   Q, Q, T1T1)
         ENDIF
C
      ELSE
C
C We have to do it out of core
C
         DISMAX  = MAXSIZE/DISSYH
         IOFFSET = 1
         DISLEFT = NUMSYH
C
 10      CONTINUE
C
         DISREAD = MIN(DISMAX, DISLEFT)
C
         IF (DISREAD .EQ. 0) THEN
            WRITE(LUOUT, *) ' @-MKDBLAA5',' No out of core algorithm.'
            CALL ERREX
         ENDIF 
C     
         DISLEFT = DISLEFT - DISREAD
C     
         CALL GETLST(H, IOFFSET, DISREAD, 2, IRREPHBR, LISTH)
C     
         IF (ITYPE .EQ. 1) THEN
C     
            CALL XGEMM('N', 'N', DISTMP, DISREAD, DISSYH, ONE, T1T1, 
     &                  DISTMP, H, DISSYH, ONE, Q(1, IOFFSET),
     &                  DISSYQ)
C     
         ELSE IF (ITYPE .EQ. 6) THEN
C     
            CALL XGEMM('N', 'N', DISSYH, NUMTMP, DISREAD, ONE, H,
     &                  DISSYH, T1T1(IOFFSET, 1), DISTMP, ONE, Q,
     &                  DISSYQ)
C
         ENDIF
C
         IOFFSET = IOFFSET + DISREAD
C
         IF(DISLEFT .NE. 0) GOTO 10
C
      ENDIF         
C
      IF (IUHF .EQ. 0 .AND. ISYTYP(1, 233) .EQ. 5 .AND. ITYPE
     &    .EQ. 6) THEN
C
         CALL SYMEXP4(IRREPQR, VRT1, POP1, DISSYQ, DISSYH, NUMSYQ,
     &                Q, Q, T1T1)
      ENDIF
C
         IF (IFLAGS(1) .GE. 20) THEN
C
            NSIZE = NUMSYQ*DISSYQ
            CALL HEADER('Checksum @-T1T1IND2 per sym. block', 0, LUOUT)
            
            WRITE(LUOUT, *) ' ABAB =  ', SDOT(NSIZE, Q, 1, Q, 1)

         ENDIF
C
C Update the Q(Ab,Ij) distributions in the disk.
C
      CALL SUMSYM2(Q, T1T1, DISSYQ*NUMSYQ, 1, IRREPQR, LISTQ)
C
      RETURN
      END

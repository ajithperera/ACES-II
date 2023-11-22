C
      SUBROUTINE DRVSYM (MAXCENT, NMAX, NSYMOP, ISYTYP, MULNUC, IPTCNT,
     &                   IDEGEN, NUCDEG, NUCPRE, NUCNUM, COOO, INDEX, 
     &                   IRREPERT, CSTRA, SYMTRN, SYMINVRT)
C
C  Calculate symmetry-adapted nuclear coordinates, in the case
C  of magnetic perturbations calculate symmetry adapted nuclear
C  magnetic moments.
C
C  Variable directory
C  NMAX   : Number of symmetry unique centers
C  NSYMOP : The total number of unique symmetry operations
C  ISYTYP : The transformation properties of X, Y, Z with 
C           respect to symmetry operations
C  MULNUC : Information about symmetry equivalent nuclei
C  MAXLOP : MAXLOP + 1 is number of irreducible representations
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION ISYTYP(3), MULNUC(MAXCENT), IPTCNT(6*MAXCENT, 0:7),
     &          IDEGEN(MAXCENT), NUCDEG(MAXCENT), NUCPRE(MAXCENT),
     &          NUCNUM(MAXCENT, 8), COOO(MAXCENT, 3), NCRREP(0:7),
     &          CSTRA(6*MAXCENT, 6*MAXCENT), ISYTYPO(3), PT(0:7),
     &          FMULT(0:7), MULT(0:7), KA(8, 3),
     &          INDEX(6*MAXCENT), SYMTRN(6*MAXCENT, 6*MAXCENT),
     &          IRREPERT(6*MAXCENT), SYMINVRT(6*MAXCENT, 6*MAXCENT)
C
      LOGICAL JFC, JPSO, JSD, NUCLEI
C
#include "mxcbf.par"
#ifdef _NAMN_THROUGH_CHRTMP
C SG 8/25/98 NAMN changed to be character instead of integer
      CHARACTER*4 NAMN(MXCBF)
      INTEGER DIM_CHRTMP
      PARAMETER (DIM_CHRTMP=MXCBF*4/10)
      CHARACTER*(DIM_CHRTMP) CHRTMP
      CHARACTER*1 CHRNUM
#else
      DIMENSION NAMN(MXCBF)
#endif /* _NAMN_THROUGH_CHRTMP */
C
      CHARACTER*4 KWO(20)
      DATA KWO /'s   ','x   ','y   ','z   ',
     &          'xx  ','xy  ','xz  ','yy  ','yz  ','zz  ',
     &          'xxx ','xxy ','xxz ','xyy ','xyz ','xzz ',
     &          'yyy ','yyz ','yzz ','zzz '/
C
      COMMON /NMR/JFC, JPSO, JSD, NUCLEI
      COMMON /MACHSP/ IINTLN, IFLTLN, IINTFP, IALONE, IBITWD
      COMMON /FILES/ LUOUT, MOINTS
      COMMON /FLAGS/ IFLAGS(100)
C     
C Evaluate the variable MAXLOP
C     
      IBTAND(I, J) = AND(I, J)
      IBTOR(I, J)  = OR(I, J)
      IBTSHL(I, J) = ISHFT(I, J)
      IBTSHR(I, J) = ISHFT(I, -J)
      IBTXOR(I, J) = XOR(I, J)
C      
      MAXLOP = 2**NSYMOP - 1 
      IPRINT = IFLAGS(1) 
C
      DO 5 I = 1, 8
         IS = 0
         DO 6  K = 1, 3
            KA(I, K) = (I - 1)/(2**(K - 1))- 2*((I - 1)/(2**K))
            IS = IS + KA(I, K)
 6       CONTINUE
         PT(I - 1)    = (-1.0)**IS
         MULT(I - 1)  = 2**MAX(0, NSYMOP - IS)
         FMULT(I - 1) = FLOAT(MULT(I-1))                         
 5    CONTINUE
C
      DO 10 IREP = 0, 7
         DO 20 ICENTA = 1, 6*MAXCENT
            IPTCNT(ICENTA, IREP) = 0
 20      CONTINUE
 10   CONTINUE
C
      CALL ZERO(CSTRA, 36*MAXCENT*MAXCENT)
C
      IOFF = 0
C
C NUCLEAR COORDINATES
C
      IF (NUCLEI) THEN
C
         DO 30 IREP = 0, MAXLOP
            ICENTA = 0
            DO 40 ICENT = 1, NMAX
               MULC = MULNUC(ICENT)
               DO 50 IDIRN = 1, 3
                  ICENTA = ICENTA + 1
                  IF (IBTAND(MULC, IBTXOR(IREP, ISYTYP(IDIRN))) 
     &                .EQ. 0) THEN
                     IOFF = IOFF + 1
                     IPTCNT(ICENTA, IREP) = IOFF
                  ENDIF
 50            CONTINUE
 40         CONTINUE
 30      CONTINUE
C
      ELSE IF(JFC) THEN
C     
         DO 60 IREP = 0, MAXLOP
            ICENTA = 0
            DO 70 ICENT = 1, NMAX
               MULC = MULNUC(ICENT)
               ICENTA = ICENTA + 1
               IF(IBTAND(MULC, IREP) .EQ. 0) THEN
                  IOFF = IOFF + 1
                  IPTCNT(ICENTA, IREP) = IOFF
               ENDIF
 70         CONTINUE
 60      CONTINUE
C
      ELSE IF(JSD) THEN
C     
         DO 80 IREP = 0, MAXLOP
            ICENTA = 0
            DO 90 ICENT = 1, NMAX
               MULC = MULNUC(ICENT)
               DO 100 IDIRN = 1, 3
                  DO 120 JDIRN = IDIRN, 3
                     ICENTA = ICENTA + 1
                     ISYTIJ = IBTXOR(ISYTYP(IDIRN), ISYTYP(JDIRN))
                     IF (IBTAND(MULC, IBTXOR(IREP, ISYTIJ)) .EQ. 0)
     &                  THEN
                        IOFF = IOFF + 1
                        IPTCNT(ICENTA, IREP) = IOFF
                     ENDIF
 120              CONTINUE
 100           CONTINUE
 90         CONTINUE
 80      CONTINUE
C     
      ELSE IF(JPSO) THEN
C
C NUCLEAR MAGNETIC MOMENTS
C
         DO 130 IREP = 0, MAXLOP
            ICENTA = 0
            DO 140 ICENT = 1, NMAX
               MULC = MULNUC(ICENT)
               DO 150 IDIRN = 1, 3
                  IF(IDIRN .EQ. 1) THEN
                     ISYMTYP =IBTXOR(ISYTYP(2), ISYTYP(3))
                  ELSE IF(IDIRN .EQ. 2) THEN
                     ISYMTYP = IBTXOR(ISYTYP(1), ISYTYP(3))
                  ELSE IF(IDIRN.EQ.3) THEN
                     ISYMTYP = IBTXOR(ISYTYP(1), ISYTYP(2))
                  ENDIF
                  ICENTA = ICENTA + 1
                  IF (IBTAND(MULC,IBTXOR(IREP,ISYMTYP)) .EQ. 0) THEN
                     IOFF = IOFF + 1
                     IPTCNT(ICENTA, IREP) = IOFF
                  ENDIF
 150           CONTINUE
 140        CONTINUE
 130     CONTINUE
C
      ENDIF
C
#ifdef _NAMN_THROUGH_CHRTMP
c   o NAMN is divided into 10 'columns'
      DO J = 0, 9
c      o CHRTMP will hold column 1+J
         WRITE (CHRNUM,'(I1)') J
         CALL GETCREC (20,'JOBARC','ATMSYMB'//CHRNUM,DIM_CHRTMP,CHRTMP)
         DO I = 0, DIM_CHRTMP-4, 4
            NAMN(1+(DIM_CHRTMP*J+I)/4) = CHRTMP(1+I:4+I)
         END DO
      END DO
#else
      CALL GETREC (20, 'JOBARC', 'ATOMSYMB', 300, NAMN)
#endif /* _NAMN_THROUGH_CHRTMP */
C
      IF (JFC) THEN
         CALL HEADER('Symmetry-adapted nuclei',-1, 6)
      ELSE IF(JSD) THEN
        CALL HEADER('Symmetry-adapted spin-dipolar integrals',-1, 6) 
      ELSE
         CALL HEADER('Symmetry-adapted PSO integrals',-1, 6)
      ENDIF
C
       ICENTA = 0
C
       IF(JFC) THEN
          WRITE(LUOUT, '(11X,8I4)') (IREP, IREP = 0,MAXLOP)
          DO 1000 ICENT =1, NMAX
             WRITE(LUOUT, '(4X,A4,3X,8I4)')
     &       NAMN(ICENT), (IPTCNT(ICENTA + 1, IREP), IREP = 0, MAXLOP)
             ICENTA=ICENTA+1
 1000     CONTINUE
       ELSE IF(JSD) THEN
          WRITE(LUOUT, '(12X,8I4)') (IREP, IREP = 0,MAXLOP)
          DO 1010 ICENT = 1, NMAX
             WRITE(LUOUT,'(4X,A4,2X,A2,8I4)')
     &       NAMN(ICENT),KWO(5),(IPTCNT(ICENTA+1,IREP),IREP=0,MAXLOP)
             WRITE(LUOUT,'(4X,A4,2X,A2,8I4)')
     &       NAMN(ICENT),KWO(6),(IPTCNT(ICENTA+2,IREP),IREP=0,MAXLOP)
             WRITE(LUOUT,'(4X,A4,2X,A2,8I4)')
     &       NAMN(ICENT),KWO(7),(IPTCNT(ICENTA+3,IREP),IREP=0,MAXLOP)
             WRITE(LUOUT,'(4X,A4,2X,A2,8I4)')
     &       NAMN(ICENT),KWO(8),(IPTCNT(ICENTA+4,IREP),IREP=0,MAXLOP)
             WRITE(LUOUT,'(4X,A4,2X,A2,8I4)')
     &       NAMN(ICENT),KWO(9),(IPTCNT(ICENTA+5,IREP),IREP=0,MAXLOP)
             WRITE(LUOUT,'(4X,A4,2X,A2,8I4)')
     &       NAMN(ICENT),KWO(10),(IPTCNT(ICENTA+6,IREP),IREP=0,MAXLOP)
             ICENTA = ICENTA + 6
 1010     CONTINUE
       ELSE
          WRITE(LUOUT,'(11X,8I4)') (IREP, IREP = 0,MAXLOP)
          DO 1020 ICENT = 1, NMAX
             WRITE(LUOUT,'(4X,A4,2X,A1,8I4)')
     &       NAMN(ICENT), KWO(2),(IPTCNT(ICENTA+1,IREP),IREP=0,MAXLOP)
             WRITE(LUOUT,'(4X,A4,2X,A1,8I4)')
     &       NAMN(ICENT), KWO(3),(IPTCNT(ICENTA+2,IREP),IREP=0,MAXLOP)
             WRITE(LUOUT,'(4X,A4,2X,A1,8I4)')
     &       NAMN(ICENT), KWO(4),(IPTCNT(ICENTA+3,IREP),IREP=0,MAXLOP)
             ICENTA = ICENTA + 3
 1020     CONTINUE
      ENDIF
C
C Determine the number of symmetry equivalent atoms (degeneracy)
C      
      DO 160 NN = 1, NMAX
         IVARB = MULNUC(NN)
         NDEG = 0
         DO 170 LA = 0, MAXLOP
            IF (IBTAND(IVARB, LA) .EQ. 0) NDEG = NDEG + 1
 170     CONTINUE
         IDEGEN(NN) = NDEG
 160  CONTINUE
C     
      NUCPRE(1) = 0
      NUCDEP    = 0
C
      DO 180 JATOM = 1, NMAX
         IDEG = IDEGEN(JATOM)
         NUCDEG(JATOM) = IDEG
         NUCPRE(JATOM + 1) = NUCPRE(JATOM) + IDEG
         ICENT = NUCPRE(JATOM)
         MULX  = MULNUC(JATOM)
         DO 190 ISYMOP = 0, MAXLOP
            IF (IBTAND(MULX, ISYMOP) .EQ. 0) THEN
               NUCDEP = NUCDEP + 1
               ICENT = ICENT + 1
               NUCNUM(JATOM, ISYMOP + 1) = ICENT
            ELSE
               NUCNUM(JATOM, ISYMOP + 1) = 0
            ENDIF
 190     CONTINUE
 180  CONTINUE
C
      NATOMS = NUCDEP
C
C Calculate the value of MAXLOT. It is the really the number
C of irreducible representation minus one. The following piece
C of code would do it but this should be simply the MAXLOP.
C
      MAXLOT = 1
      MAXLO  = 0
C     
      IF (NSYMOP .GT. 0) THEN
         II = 1
         DO 210 I = 1, NSYMOP
            IF(IBTAND(II, MAXLO) .EQ. 0) MAXLOT = 2*II
            II = 2*II
 210     CONTINUE
      ENDIF
C
      MAXLOT = MAXLOT - 1
C
      IF (JFC .OR. JSD .OR. NUCLEI) THEN
         ISYTYPO(1)=ISYTYP(1)
         ISYTYPO(2)=ISYTYP(2)
         ISYTYPO(3)=ISYTYP(3)
      ELSE 
         ISYTYPO(1) = IBTXOR(ISYTYP(2), ISYTYP(3))
         ISYTYPO(2) = IBTXOR(ISYTYP(3), ISYTYP(1))
         ISYTYPO(3) = IBTXOR(ISYTYP(1), ISYTYP(2))
      ENDIF
C
      IF (JFC) NCOOR = NATOMS
      IF (JPSO) NCOOR = 3*NATOMS
      IF (JSD) NCOOR = 6*NATOMS
C
      IF (JFC) THEN
         DO 220 IREP = 0, MAXLOP
            NCRREP(IREP) = 0
            DO 230 ICENT = 1, NMAX
               MULCNT = MULNUC(ICENT)
               ISCOOR = IPTCNT(ICENT, IREP)
               IF (ISCOOR .GT. 0) THEN
                  NCRREP(IREP) = NCRREP(IREP) + 1
                  IVAR = IBTXOR(IREP, 0)
                  DO 240 ISYMOP = 0, MAXLOT
                     IF (IBTAND(ISYMOP, MULCNT) .EQ. 0) THEN
                        ICCOOR = NUCNUM(ICENT, ISYMOP + 1)
                        CSTRA(ISCOOR, ICCOOR)=PT(IBTAND(ISYMOP, IVAR))
                     ENDIF
 240              CONTINUE
               ENDIF
 230        CONTINUE
 220     CONTINUE
      ELSE IF (JSD) THEN
         DO 250 IREP = 0, MAXLOP
            NCRREP(IREP) = 0
            DO 260 ICENT = 1, NMAX
               MULCNT = MULNUC(ICENT)
               KCOOR = 0
               DO 270 ICOOR = 1, 3
                  DO 280 JCOOR = ICOOR, 3
                     KCOOR = KCOOR + 1
                     ISCOOR = IPTCNT(6*(ICENT - 1) + KCOOR, IREP)
                     IF (ISCOOR .GT. 0) THEN
                        NCRREP(IREP) = NCRREP(IREP) + 1
                        IVAR = IBTXOR(IREP, IBTXOR(ISYTYPO(ICOOR),
     &                                ISYTYPO(JCOOR)))
                        DO 290 ISYMOP = 0, MAXLOT
                           IF (IBTAND(ISYMOP, MULCNT) .EQ. 0) THEN
                              ICCOOR = 6*(NUCNUM(ICENT, ISYMOP + 1) 
     &                                 - 1) + KCOOR
                              CSTRA(ISCOOR, ICCOOR) = PT(IBTAND(ISYMOP,
     &                                                IVAR))
                           ENDIF
 290                    CONTINUE
                     ENDIF
 280              CONTINUE
 270           CONTINUE
 260        CONTINUE
 250     CONTINUE
      ELSE
         DO 300 IREP = 0, MAXLOP
            NCRREP(IREP) = 0
            DO 310 ICENT = 1, NMAX
               MULCNT = MULNUC(ICENT)
               DO 320 ICOOR = 1, 3
                  ISCOOR = IPTCNT(3*(ICENT - 1) + ICOOR, IREP)
                  IF (ISCOOR .GT. 0) THEN
                     NCRREP(IREP) = NCRREP(IREP) + 1
                     IVAR = IBTXOR(IREP, ISYTYPO(ICOOR))
                     DO 330 ISYMOP = 0, MAXLOT
                        IF (IBTAND(ISYMOP, MULCNT) .EQ. 0) THEN
                           ICCOOR = 3*(NUCNUM(ICENT, ISYMOP + 1) 
     &                              - 1) + ICOOR
                           CSTRA(ISCOOR, ICCOOR) = PT(IBTAND(ISYMOP,
     &                                                IVAR))
                        ENDIF
 330                 CONTINUE
                  ENDIF
 320           CONTINUE
 310        CONTINUE
 300     CONTINUE
      ENDIF
C
C Everything is done, print if required 
C 
      WRITE (LUOUT,'(/,1X,A,8I3,/)')
     &      ' Number of coordinates in each symmetry: ',
     &       (NCRREP(I),I=0,MAXLOP)
C
      IF (IPRINT .GE. 20) THEN
         CALL HEADER('Transformation from Cartesian to symmetry '
     &               //'coordinates',-1, 6)
         DO 340 ISYM = 1, NCOOR
            WRITE (LUOUT,'(9F6.2)') (CSTRA(ISYM,J),J=1,NCOOR)
 340     CONTINUE
      END IF
C
C For our case we need to rearrange the CSTRA matrix
C
      KK = 0
      IF (JFC) THEN
         NPERT = NMAX
      ELSE IF (JPSO) THEN
         NPERT = 3*NMAX
      ELSE
         NPERT = 6*NMAX
      ENDIF
C 
      DO 350 I = 1, NPERT
         DO 360 J = 0, MAXLOP
            IF (IPTCNT(I, J) .GT. 0) THEN
               KK = KK + 1
               INDEX(KK) = IPTCNT(I, J)
               IRREPERT(KK) = J + 1 
            ENDIF
 360     CONTINUE
 350  CONTINUE

      DO 370 I = 1, NCOOR
         DO 380 J = 1, NCOOR
            SYMTRN(I, J) = CSTRA(INDEX(I), J)
 380        CONTINUE
 370     CONTINUE
C
      IF (IPRINT .GE. 20) THEN
         WRITE (LUOUT,'(/,1X,A,8I3,/)')
     &         ' Number of coordinates in each symmetry: ',
     &        (NCRREP(I),I=0,MAXLOP)
         CALL HEADER('Transformation from Cartesian to symmetry '
     &               //'coordinates', -1, 6)
         DO 390 ISYM = 1, NCOOR
            WRITE (LUOUT,'(9F6.2)') (SYMTRN(ISYM,J),J=1,NCOOR)
 390     CONTINUE
      END IF
C
C Inverse transformation matrix
C
      DO 400 I = 1, NCOOR
         DO 410 J = 1, NCOOR
            SYMINVRT(I, J) = CSTRA(J, I)
 410     CONTINUE
 400  CONTINUE
C
      DO 420 ICOL = 1, NCOOR
         SUM = 0.0D+00
         DO 430 IROW = 1, NCOOR
            SUM = SUM + ABS(SYMINVRT(IROW, ICOL))
 430     CONTINUE
C     
         DO 440 IROW = 1, NCOOR
            SYMINVRT(IROW, ICOL) = SYMINVRT(IROW, ICOL)/SUM
 440     CONTINUE
 420  CONTINUE
C
      IF (IPRINT .GE. 20) THEN
         CALL HEADER('Inverse transformation matrix', -1, 6)
         DO 450 ISYM = 1, NCOOR
            WRITE (LUOUT,'(9F6.2)') (SYMINVRT(ISYM,J),J=1,NCOOR)
 450     CONTINUE
      ENDIF
C
      RETURN
      END

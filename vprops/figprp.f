      SUBROUTINE FIGPRP(WHAT,WHERE,NMAX,NATOMS)
C
C THIS ROUTINE IS THE CRAPS-VPROPS INTERFACE AND TELLS
C VPROPS WHICH PROPERTIES SHOULD BE COMPUTED.
C
C Added 08/93 and 03/94 by Ajith
C Rewritten 8/27/97 by SG
C
      IMPLICIT NONE
C
      INTEGER NMAX, NATOMS
      INTEGER WHAT(NMAX),WHERE(NMAX)
C
      INTEGER IFLAGS, IFLAGS2
      COMMON /FLAGS/ IFLAGS(100)
      COMMON /FLAGS2/ IFLAGS2(500)
C
      INTEGER I
C
      CALL IZERO(WHERE, NMAX)
      CALL IZERO(WHAT, NMAX)
      DO 5 I = 1,(9+9*NATOMS)
        WHAT(I) = 27
 5    CONTINUE
C
      IF (IFLAGS2(142) .EQ. 1) THEN
C
C CALCULATE MOST OF THE REGULARLY USED INTEGRALS. THE ONES THAT ARE
C NOT COMPUTED: Thid moment(15), velocity(16), overlap(10), hyperfine 
C splitting(13), Diamagnetic sucep (11), dipole velocity(22), 
C electric field(2), planar density(7) and linear density(8),
C    
        WHAT(1)=4
        WHAT(2)=5
        WHAT(3)=6
        WHAT(4)=12
        WHAT(5)=14
        WHAT(6)=22
        WHAT(8)=24
        WHAT(9)=26
C      
        DO 10 I=1,NATOMS
C
          WHAT(9+I)= 9
          WHAT(9+NATOMS+I)= 3
          WHAT(9+2*NATOMS+I)= 1
C
C Added 08/93 03/94 Ajith
C
C          WHAT(6+3*NATOMS+I)= 2
C
          WHAT(9+4*NATOMS+I)= 17
          WHAT(9+5*NATOMS+I)= 18
          WHAT(9+6*NATOMS+I)= 19
          WHAT(9+7*NATOMS+I)= 20
          WHAT(9+8*NATOMS+I)= 21
          WHAT(9+9*NATOMS+I)= 23
C
          WHERE(9+I)=I
          WHERE(9+NATOMS+I)= I
          WHERE(9+2*NATOMS+I)= I
C
C Added 08/93 03/94 Ajith
C
C         WHERE(7+3*NATOMS+I)= I
C
          WHERE(9+4*NATOMS+I)= I
          WHERE(9+5*NATOMS+I)= I
          WHERE(9+6*NATOMS+I)= I
          WHERE(9+7*NATOMS+I)= I         
          WHERE(9+8*NATOMS+I)= I         
          WHERE(9+9*NATOMS+I)= I         
C
 10     CONTINUE
C
      ELSE

C The core excitation calculation need quadrupole,octupole and
C angular momentum integrals. Ajith Perera,02/2020.
C The transition moments between states also need the same moments
C Ajith Perera, 11/2023

        IF (IFLAGS2(119) .GT. 0 .OR. IFLAGS2(185) .GT. 0) Then
           WHAT(1) = 4
           WHAT(2) = 5
           WHAT(3) = 6
           WHAT(4) = 12
           WHAT(5) = 14
           WHAT(6) = 24
           WHAT(7) = 26
           WHAT(8) = 22 
        ELSE

C Figure out which integrals are needed We always need 
C second moment and dipole moment, velocity and angular 
C momentum integrals
C
           WHAT(1) = 4
           WHAT(2) = 14
           WHAT(3) = 22
           WHAT(4) = 24
        ENDIF 
C
        IF (IFLAGS(18) .EQ. 5) THEN
          WHAT(3) = 22
          DO 100 I = 1,NATOMS
            WHAT(3 + NATOMS + I) = 20
            WHAT(3 + 2*NATOMS + I) = 21
            WHERE(3 + NATOMS + I) = I
            WHERE(3 + 2*NATOMS + I) = I
 100      CONTINUE
C
        ELSE IF (IFLAGS(18) .EQ. 8) THEN
          DO 110 I = 1,NATOMS
            WHAT(2 + NATOMS + I) = 17
            WHAT(2 + 2*NATOMS + I) = 18
            WHERE(2 + NATOMS + I) = I
            WHERE(2 + 2*NATOMS + I) = I
 110      CONTINUE
C
        ELSE IF (IFLAGS(18) .EQ. 9) THEN
          DO 120 I = 1,NATOMS
            WHAT(2 + NATOMS + I) = 9
            WHERE(2 + NATOMS + I) = I
 120      CONTINUE
C
        ELSE IF (IFLAGS(18) .EQ. 10) THEN
          DO 130 I = 1,NATOMS
            WHAT(2 + NATOMS + I) = 19
            WHERE(2 + NATOMS + I) = I
 130      CONTINUE
C
        ELSE IF (IFLAGS(18) .EQ. 11) THEN
          WHAT(2) = 5
          WHAT(3) = 6
          WHAT(4) = 14
C
        ELSE IF (IFLAGS(18) .EQ. 13) THEN
          DO 140 I = 1,NATOMS
            WHAT(2 + NATOMS + I) = 9
            WHAT(2 + 2*NATOMS + I) = 17
            WHAT(2 + 3*NATOMS + I) = 18
            WHAT(2 + 4*NATOMS + I) = 19
            WHERE(2 + NATOMS + I) = I
            WHERE(2 + 2*NATOMS + I) = I
            WHERE(2 + 3*NATOMS + I) = I
            WHERE(2 + 4*NATOMS + I) = I
 140      CONTINUE
C
        ELSE IF (IFLAGS(18) .NE. 0 .OR. IFLAGS2(103) .EQ. 1) THEN
          WHAT(2) = 5
          WHAT(3) = 6
          WHAT(4) = 12
          WHAT(5) = 14
          WHAT(6) = 24
          DO 150 I = 1,NATOMS
            WHAT(6 + NATOMS + I) = 9
            WHAT(6 + 2*NATOMS + I) = 3
            WHAT(6 + 3*NATOMS + I) = 1
            WHAT(6 + 4*NATOMS + I) = 23
            WHERE(6 + NATOMS + I) = I
            WHERE(6 + 2*NATOMS + I) = I
            WHERE(6 + 3*NATOMS + I) = I
            WHERE(6 + 4*NATOMS + I) = I
 150      CONTINUE
C Moved out of this block. The logic now say if DKH do these 
C integrals regardless of what the property key-word says.
C Ajith Perera, 12/2018 (without this change Young could not do
C DKH time-dependent EOM calcs).
        
CSSS       ELSE IF (IFLAGS2 (167) .GT. 0) THEN ! Watson, h_IFLAGS2_dkh_order
CSSS          DO 160 I = 1, NATOMS
CSSS            WHAT(5 + 4*NATOMS + I) = 23
CSSS            WHERE(5 + 4*NATOMS + I) = I
CSSS160       CONTINUE

       ENDIF
C
      ENDIF

      IF (IFLAGS2 (167) .GT. 0) THEN ! Watson, h_IFLAGS2_dkh_order
          DO 160 I = 1, NATOMS
            WHAT(5 + 4*NATOMS + I) = 23
            WHERE(5 + 4*NATOMS + I) = I
160       CONTINUE
      ENDIF
C
      RETURN
      END

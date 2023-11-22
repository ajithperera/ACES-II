      SUBROUTINE INSING(METHOD,IUHF)
C
      IMPLICIT INTEGER(A-Z)
      LOGICAL NONHF
      COMMON/NHFREF/NONHF
      COMMON/SYM/POP(8,2),VRT(8,2),NT(2),NF1(2),NF2(2)
C
c      if (calc="sdq-mbpt(4)" or calc="mbpt(4)")
      IF (METHOD.EQ.3.OR.METHOD.EQ.4) THEN
         DO ISPIN = 1, IUHF+1
            CALL UPDMOI(1,NF1(ISPIN),ISPIN,91,0,0)
            CALL UPDMOI(1,NF2(ISPIN),ISPIN,92,0,0)
            CALL UPDMOI(1,NT(ISPIN), ISPIN,93,0,0)
cYAU - initialize F(M,I), F(E,A), and F(A,I) to zero
            call aces_list_memset(ispin,91,0)
            call aces_list_memset(ispin,92,0)
            call aces_list_memset(ispin,93,0)
         END DO
      END IF
C
      IF (METHOD.GE.5.OR.NONHF) THEN
         DO ISPIN = 1, IUHF+1
            CALL UPDMOI(1,NT(ISPIN), ISPIN+2,90,0,0)
            CALL UPDMOI(1,NF1(ISPIN),ISPIN,  91,0,0)
            CALL UPDMOI(1,NF2(ISPIN),ISPIN,  92,0,0)
            CALL UPDMOI(1,NT(ISPIN), ISPIN,  93,0,0)
cYAU - initialize F(M,I), F(E,A), and F(A,I) to zero
            call aces_list_memset(ispin,91,0)
            call aces_list_memset(ispin,92,0)
            call aces_list_memset(ispin,93,0)
         END DO
      END IF
C
      RETURN
      END

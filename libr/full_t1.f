
      SUBROUTINE FULL_T1(T1_EXPAND, T1_DROP, NVRT_ACT, NOCC_ACT,
     &                   NVRT, NOCC, NVRT_DRP, NOCC_DRP)

      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      DIMENSION T1_EXPAND(NVRT,NOCC), T1_DROP(NVRT_ACT,NOCC_ACT)

      CALL ZERO(T1_EXPAND,NVRT*NOCC)
      DO IORBS = 1, NVRT_ACT
         DO JORBS = 1, NOCC_ACT
cSSS            IF_INDEX = IORBS + NVRT_DRP
            IF_INDEX = IORBS
            JF_INDEX = JORBS + NOCC_DRP
            T1_EXPAND(IF_INDEX,JF_INDEX) = T1_DROP(IORBS,JORBS)
         END DO
      END DO

CSSS      call output(T1_DROP,1,nvrt_act,1,nocc_act,nvrt_act,
CSSS     &            nocc_act,1)
CSSS      write(6,*) "expanded"
CSSS      call output(T1_EXPAND,1,Nvrt,1,nocc,nvrt,nocc,1)
      RETURN
      END


      Subroutine Set_lmn_delta(LMNVAL,LVAL)

      Implicit none
      integer MAXANG
      parameter (Maxang = 7)
      integer LMNval(3, MAXANG*(MAXANG+1)*(MAXANG+2)/2)
      integer ii,lval,m,left,l,n,i,i1
C
      II = 0
      DO L = LVAL,0,-1
         LEFT = LVAL - L
         DO M = LEFT,0,-1
            N = LEFT - M
            II = II + 1
            LMNVAL(1,II) = L
            LMNVAL(2,II) = M
            LMNVAL(3,II) = N
          ENDDO
      ENDDO
C     
#ifdef _DEBUG_LVL0     
      Write(6,"(a)") "The lmn values"
      Do i=1, 3
         Write(6,"(10I3)") (lmnval(i,
     &                      i1),i1=1,MAXANG*(MAXANG+1)*(MAXANG+2)/2)
         Write(6,*)
      Enddo
#endif
c     
      Return
      End

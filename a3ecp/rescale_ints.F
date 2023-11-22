      Subroutine rescale_ints(SOI, ACES_SOI, Nbasis, ERD_index,
     &                        ERD_scale)

      Implicit double precision (a-h,o-z)

      Dimension SOI(Nbasis,Nbasis,Nbasis,Nbasis), ERD_scale(Nbasis)
      Dimension ACES_SOI(Nbasis,Nbasis,Nbasis,Nbasis)
      Integer ERD_index(Nbasis)

#ifdef  _DEBUG_LVL0
        Write(*,"(a)") "Reorder Index array"
        Write(*,"(6I5)") (ERD_Index(i), i=1, Nbasis)
        Write(*,"(a)") "The scale values"
        Write(*,"(1x,F15.7)") (ERD_scale(i), i=1, Nbasis)
#endif

      Do I = 1, Nbasis
         DO J = 1, Nbasis
            DO K = 1, Nbasis
               DO L= 1, Nbasis

                  II = ERD_index(I)
                  JJ = ERD_index(J)
                  KK = ERD_index(K)
                  LL = ERD_index(L)
                    
                  C1 = ERD_scale(I)
                  C2 = ERD_scale(J)
                  C3 = ERD_scale(K)
                  C4 = ERD_scale(L)
                  
                  SCALE = C1*C2*C3*C4
                  
                  ACES_SOI(II,JJ,KK,LL) = SOI(I,J,K,L)*SCALE

               ENDDO
           ENDDO
         ENDDO
      ENDDO

      RETURN
      END







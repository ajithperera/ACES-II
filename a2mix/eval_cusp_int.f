      SUBROUTINE EVAL_CUSP_INT(IMOMN, JMOMN, KMOMN, LMOMN, NIPRIM,
     &                         NJPRIM, NKPRIM, NLPRIM, NTOTPRIM,
     &                         NTOTCRF,ALPHA, PCOEF, COORD_IX,
     &                         COORD_IY, COORD_IZ, COORD_JX, COORD_JY,
     &                         COORD_JZ, COORD_KX, COORD_KY, COORD_KZ, 
     &                         COORD_LX, COORD_LY, COORD_LZ, UVALUE, 
     &                         CUSPINT)
C
      DIMENSION ALPHA(NTOTPRIM), PCOEF(NTOTPRIM, NTOTCRF)
C
      DO IPRIM = 1, NIPRIM
         DO JPRIM = 1, NJPRIM
            DO KPRIM = 1, NKPRIM
               DO LPRIM = 1, NLPRIM

                  ALPHAA = ALPHA(IPRIM)
                  ALPHAB = ALPHA(JPRIM)
                  ALPHAC = ALPHA(KPRIM)
                  ALPHAD = ALPHA(LPRIM)
C
                  CALL ABCDPU(ALPHAA, ALPHAB, ALPHAC, ALPHAD, IMOM, 
     &                        JMOM, KMOM, LMOM, UVALUE)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
      RETURN
      END

      SUBROUTINE A2BUILD_TENSOR(WORK, MAXCOR, PRDUTINT, TENSOR_ATR,
     &                          NMBR_OF_PERTS, NAOBFNS, NATOMS, 
     &                          IBEGIN_P_DENS, IBEGIN_P_OPRTS,
     &                          IBEGIN_AO_OVRLP, IBEGIN_MO_OVRLP, 
     &                          IBEGIN_MO_VECTS, IBEGIN_MEM_LEFT,
     &                          JPERT)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      DIMENSION WORK(MAXCOR), PRDUTINT(NAOBFNS, NAOBFNS),
     &          TENSOR_ATR(NMBR_OF_PERTS)
C
      DATA IZR0, IONE /0, 1/
C
#ifdef _TRASHON
      IOFF_4TMP    = IBEGIN_MEM_LEFT
      IOFF_4PRDINT = IOFF_4TMP + NAOBFNS*NAOBFNS
      IEND         = IOFF_4PRDINT + NAOBFNS*NAOBFNS

      IF (IEND .GE. MAXCOR) CALL INSMEM("A2BUILD_TENSOR", IEND,
     &                                   MAXCOR)
C
      CALL XGEMM('T', 'N', NAOBFNS, NAOBFNS, NAOBFNS, 1.0D0,
     &            WORK(IBEGIN_MO_VECTS), NAOBFNS, PRDUTINT,
     &            NAOBFNS, 0.0D0,WORK(IOFF_4TMP),NAOBFNS)

      CALL XGEMM('N', 'N', NAOBFNS, NAOBFNS, NAOBFNS, 1.0D0,
     &            WORK(IOFF_4TMP), NAOBFNS, WORK(IBEGIN_MO_VECTS),
     &            NAOBFNS, 0.0D0,WORK(IOFF_4PRDINT),NAOBFNS)
C
#ifdef _DEBUG_LVLM1
      CALL HEADER ('The  MO vectors', 1, 6)
      call output(WORK(IBEGIN_MO_VECTS),1,NAOBFNS,1,NAOBFNS,NAOBFNS,
      CALL HEADER ('The AO Product at R after MO Scaling', 1, 6)
      call output(WORK(IOFF_4PRDINT),1,NAOBFNS,1,NAOBFNS,NAOBFNS,
     &            NAOBFNS, 1)
#endif
#endif
C
#ifdef _DEBUG_LVLM1
      CALL HEADER ('The AO Product at R before MO Scaling', 1, 6)
      call output(PRDUTINT,1,NAOBFNS,1,NAOBFNS,NAOBFNS,NAOBFNS, 1)
#endif

      DO IPERT = 1, NMBR_OF_PERTS
C
         I_TOT_PDEN_N = IBEGIN_P_DENS + (IPERT-1)*
     &                  NAOBFNS*NAOBFNS

         TOT_PDEN_N = DDOT(NAOBFNS*NAOBFNS, PRDUTINT, 1,
     &                     WORK(I_TOT_PDEN_N), 1)
C
CSSS         TENSOR_ATR(IPERT, JPERT) = TOT_PDEN_N
#ifdef _DEBUG_LVLM1
      Write(6,*)
      write(6,'(a,a,I2)') "Computing the perturbed total/spin", 
     &                    " density for perturbation: ", ipert
      call output(WORK(I_TOT_PDEN_N),1,NAOBFNS,1,NAOBFNS,NAOBFNS, 
     &            NAOBFNS, 1)
      Write(6,*)
      Write(*,'((1x,i5),F10.5)') I_TOT_PDEN_N, TOT_PDEN_N
      Write(6,*)
#endif
C
CSSS         DO JPERT = 1, NMBR_OF_PERTS
C
            I_TOT_POPRT_N = IBEGIN_P_OPRTS + (IPERT-1)*
     &                      NAOBFNS*NAOBFNS
C
            TOT_OPRT_N = DDOT(NAOBFNS*NAOBFNS, PRDUTINT, 1,
     &                        WORK(I_TOT_POPRT_N), 1)

            TOT_AO_OVRLP = DDOT(NAOBFNS*NAOBFNS, PRDUTINT, 1,
     &                        WORK(IBEGIN_AO_OVRLP), 1)

            TOT_MO_OVRLP = DDOT(NAOBFNS*NAOBFNS, PRDUTINT, 1,
     &                        WORK(IBEGIN_MO_OVRLP), 1)
C#ifdef _DEBUG_LVL0
      Write(6,*)
      write(6,'(a,a)') "Computing the perturbed total/spin",
     &                    " density and others for props: "
      Write(6,*) 
      Write(*,'((2I5, 1x),F10.5)')ipert,jpert,TOT_PDEN_N
      Write(6,*)
C#endif
C
            TENSOR_ATR(IPERT) = TOT_PDEN_N
C
#ifdef _DEBUG_LVL0
            write(6,*)
            Write(6,"(2I5,F10.5)") ipert, jpert,TENSOR_ATR(ipert)
#endif
CSSS         ENDDO
C
      ENDDO
C
      RETURN
      END 
C  

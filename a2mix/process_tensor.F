      SUBROUTINE PROCESS_TENSOR(CC_ARRAY, CC_TENSOR, GRID_POINTS,
     &                          NTPERT, IPOINT, IJPAIR, MAXCENT,
     &                          MAXGRD, MAXPAIRS, NATOMS, 
     &                          NPOINT_PAIR)

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      DIMENSION CC_ARRAY(MAXCENT,MAXGRD*MAXPAIRS), 
     &          NPOINT_PAIR(MAXPAIRS), GRID_POINTS(3,MAXPAIRS*MAXGRD),
     &          CC_TENSOR(MAXPAIRS,MAXGRD*MAXPAIRS)
C
      CHARACTER*12 FILENAME
      LOGICAL FILE_EXIST

#ifdef _DEBUG_LVLM1
      Write(6, *)
      Write(6,"(a,I3)") "The total # of points: ", IPOINT
      Write(6,"(a,5(1x,I2))") "No. points for IJPAIR: ", 
     &         (NPOINT_pAIR(I), I=1, IJPAIR)
C
      Write(6,*)
      Write(6,*) "The grid points"
      Write(6,*)
      ioff = 0
      Do i=1, IJPAIR
         Write(6,*)
         DO J = 1, NPOINT_pAIR(I)
            ioff = ioff + 1
            Write(6, "(3F10.5)"), (GRID_POINTS(k, ioff), k=1, 3)
         Enddo
      Enddo
C
      DO I =1, IJPAIR
         Write(6,*) "New Pair"
         Index = (I-1)*NPOINT_PAIR(I)
         DO K = 1, NPOINT_PAIR(I)
            Index = INdex + 1
            Write(6,*) Index
            CALL OUTPUT(CC_ARRAY(1,Index), 1, NTPERT, 1, 1,
     &                  NTPERT, 1, 1)
         ENDDO
       ENDDO
#endif
C
C Build the coupling tensor from the coupling array.
C
      IJPAIR = 0
      IKEEP  = 0
      DO IATOM = 1, NATOMS
         DO JATOM = 1, NATOMS
C
            IF (IATOM .LT. JATOM)  THEN
                IJPAIR = IJPAIR + 1

                IMV   = (IJPAIR - 1)*NPOINT_PAIR(IJPAIR)
                IFX   =  IJPAIR*NPOINT_PAIR(IJPAIR)
C
                DO NPOINT = 1, NPOINT_PAIR(IJPAIR)
                     IMV = IMV + 1
                   IKEEP = IKEEP + 1
                   CC_TENSOR(IJPAIR, IKEEP) = CC_ARRAY(JATOM,IMV)
     &                                      + CC_ARRAY(IATOM,IFX)
#ifdef _DEBUG_LVLM1
                   Write(6, "(I2)") IATOM, JATOM, IMV, IFX
                   Write(6,"(3F10.5)") CC_ARRAY(JATOM,IMV),
     &                                 CC_ARRAY(IATOM,IFX),
     &                                 CC_TENSOR(IJPAIR, IKEEP)
#endif
                ENDDO
C
            ENDIF
         ENDDO
      ENDDO
#ifdef _DEBUG_LVLM1
      Index = 0
      DO NPAIRS = 1, IJPAIR
         Write(6, "(a)") "The NMR FC-CC for the pair"
         DO K = 1, NPOINT_PAIR(NPAIRS)
           Index = INdex + 1
           Write(6, "(F10.5)") CC_TENSOR(NPAIRS, Index)
         ENDDO
      ENDDO
#endif
C
#ifdef _PRINT_OPTION1
      ICOUNT = 0
      IFLTAG = 9
C
      DO NPAIRS = 1, IJPAIR
C
         IFLTAG = IFLTAG + 1
         WRITE(FILENAME,"(a,I2)") "PDEN_VALUE",IFLTAG
         IUNIT = 70 + NPAIR
         INQUIRE(FILE=FILENAME,EXIST=FILE_EXIST)
         IF (FILE_EXIST) THEN
            OPEN (UNIT=IUNIT, FILE=FILENAME, FORM="FORMATTED",
     &            STATUS="OLD")
            CLOSE(UNIT=IUNIT, STATUS="DELETE")
         ENDIF 
         OPEN(UNIT=IUNIT, FILE=FILENAME, FORM="FORMATTED", 
     &        STATUS="NEW")

         DO NPOINTS = 1, NPOINT_PAIR(NPAIRS)
C
            ICOUNT = ICOUNT + 1
C
            XCOORD = GRID_POINTS(1, ICOUNT)
            YCOORD = GRID_POINTS(2, ICOUNT)
            ZCOORD = GRID_POINTs(3, ICOUNT)
C
            R2 = XCOORD**2 + YCOORD**2 + ZCOORD**2
C
            WRITE(IUNIT,'(5F20.10)') R2, (CC_ARRAY(j, ICOUNT), j=1, 
     &            NTPERT), CC_TENSOR(NPAIRS, ICOUNT)
C
         ENDDO
      ENDDO
#endif
#ifdef _PRINT_OPTION2
C
      IJPAIR = 0
      IKEEP  = 0
      IFLTAG = 9
      DO IATOM = 1, NATOMS
         DO JATOM = 1, NATOMS
C
            IF (IATOM .LT. JATOM)  THEN
                IJPAIR = IJPAIR + 1
                IFLTAG = IFLTAG + 1
                WRITE(FILENAME,"(a,I2)") "PDEN_VALUE",IFLTAG
                IUNIT = 70 + IJPAIR
                INQUIRE(FILE=FILENAME,EXIST=FILE_EXIST)
                IF (FILE_EXIST) THEN 
                   OPEN (UNIT=IUNIT, FILE=FILENAME, FORM="FORMATTED",
     &             STATUS="OLD")
                   CLOSE(UNIT=IUNIT, STATUS="DELETE")
                ENDIF
                OPEN(UNIT=IUNIT, FILE=FILENAME, FORM="FORMATTED",
     &               STATUS="NEW")

                IMV   = (IJPAIR - 1)*NPOINT_PAIR(IJPAIR)
                IFX   =  IJPAIR*NPOINT_PAIR(IJPAIR)
C
                DO NPOINT = 1, NPOINT_PAIR(IJPAIR)
                      IMV = IMV + 1
                    IKEEP = IKEEP + 1
                   XCOORD = GRID_POINTS(1, IKEEP)
                   YCOORD = GRID_POINTS(2, IKEEP)
                   ZCOORD = GRID_POINTs(3, IKEEP)
C
                      R = DSQRT(XCOORD**2 + YCOORD**2 + ZCOORD**2)
#ifdef _DEBUG_LVLM1
             Write(6,"(a,4F10.5)") "Perturbation is centered on ",
     &                                XCOORD, YCOORD, ZCOORD,
#endif
C
                      WRITE(IUNIT,'(4F10.4)') R, 
     &                                         CC_ARRAY(JATOM,IMV),
     &                                         CC_ARRAY(IATOM,IFX),
     &                                         CC_TENSOR(IJPAIR,IKEEP)

                ENDDO
C
            ENDIF
         ENDDO
      ENDDO
#endif
C
      RETURN
      END

      SUBROUTINE FIXDUMMY(VMLORDLABEL,LABEL_NODUMMY,TMP,ZMTORDLABEL,
     &                    IMAP,ICHARGE,NATOMS,NREALS)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER ATMCHRG

      CHARACTER*5 ZMTORDLABEL(NATOMS)
      CHARACTER*5 LABEL_NODUMMY(NREALS)
      CHARACTER*5 VMLORDLABEL(NREALS)
      CHARACTER*5 TMP(NATOMS)
      CHARACTER*5 ANAME(36)
      DIMENSION   IMAP(NATOMS)
      DIMENSION   ATMMASS(NATOMS)
      DIMENSION   ATMCHRG(NATOMS)
      DIMENSION   ICHARGE(NREALS)
      LOGICAL     FOUND 
      CHARACTER*5 LABEL_REF 

      DATA ANAME/ 'H    ', 'HE   ', 'LI   ', 'BE  ', 'B    ',
     &            'C    ', 'N    ', 'O    ', 'F   ', 'NE   ',
     &            'NA   ', 'MG   ', 'AL   ', 'SI  ', 'P    ',
     &            'S    ', 'CL   ', 'AR   ', 'K   ', 'CA   ',
     &            'SC   ', 'TI   ', 'V    ', 'CR  ', 'MN   ',
     &            'FE   ', 'CO   ', 'NI   ', 'CU  ', 'ZN   ',
     &            'GA   ', 'GE   ', 'AS   ', 'SE  ', 'BR   ',
     &            'KR   '/

C Eliminate dummy atoms from the zmtordlabel. The resulting array is
C the input order of atoms (i.e. if we were using Cartesian coordinates 
C inputs this would be how they appear in the ZMAT), 

      ICOUNT = 0
      CALL GETCREC(20,"JOBARC",'ZSYM',5*NATOMS,TMP)
      DO IATOM = 1, NATOMS
         IF (TMP(IATOM) .NE. "X    ") THEN
            ICOUNT = ICOUNT + 1
            ZMTORDLABEL(ICOUNT) = TMP(IATOM)
         ENDIF
      ENDDO 

      DO IATOM = 1, NREALS
         LABEL_NODUMMY(IATOM) = ZMTORDLABEL(IATOM)
         VMLORDLABEL(IATOM)   = ANAME(ICHARGE(IATOM))
      ENDDO 

#ifdef _DEBUG_LVLM
      Print*, (LABEL_NODUMMY(i),i=1,natoms)
      Print*, (VMLORDLABEL(i),i=1,natoms)
#endif 

C The Label array present the internal order without the dummy atoms and
C the imap carries the input-to-internal order, but it takes into 
C consideration of dummy atoms. The following logic reconstruct the imap
C array to eliminate dummies.

      DO ICOUNT = 1, NREALS
         LABEL_REF(1:5) = VMLORDLABEL(ICOUNT)(1:2)
         FOUND  = .FALSE.
         JCOUNT = 0
         DO WHILE (.NOT. FOUND)
            JCOUNT = JCOUNT + 1
            IF (LABEL_NODUMMY(JCOUNT) .EQ. LABEL_REF) THEN
               IMAP(ICOUNT)  = JCOUNT
               LABEL_NODUMMY(JCOUNT) = "DLTD"
               FOUND = .TRUE.
            ENDIF  
         ENDDO 
      ENDDO 

#ifdef _DEBUG_LVL0
      Write(6,"(a)") " Modified IMAP array after fixing for dummy atoms"
      write(6,"(6(1x,I2))") (IMAP(I),I=1,NREALS)
#endif 

      RETURN
      END

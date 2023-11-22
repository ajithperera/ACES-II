
#ifdef _UNICOS
c This is very dangerous. I.E., changing the bindings for each arch.
      SUBROUTINE INCENT(VLIST,ATOMSYM,NOCMX,NON,NOADC,INU,CENTR,IOUTU)
#else
      SUBROUTINE INCENT(VLIST,NOCMX,NON,NOADC,INU,CENTR,IOUTU)
#endif /* _UNICOS */

C     This routine fills the arrays VLIST and CENTR. This routine had
C     to be modified somewhat in order to compile under f90 on crunch.
C     VLIST is a 64-bit quantity, but parts of it are characters. Prev-
C     iously, column 5 of VLIST contained 2H character variables (atomic
C     symbols), and columns 9 and 10 contained 2H blanks. These columns
C     are now filled by CHARACTER*8 variables in subroutine INCENT_CH.
C
C     Information about contents of VLIST :
C
C     Columns  1-3 : x,y,z coordinates
C     Column     4 : charge of center
C     Column     5 : atomic symbols
C     Column     6 : atomic masses
C     Column     7 : SPIN array elements
C     Column     8 : QUADN array elements
C     Columns 9,10 : Blanks
C
C     CENTR is treated as an integer array, but really it is a character.
C     On entry it contains things like 'H #1'. In PUTLN the # is removed
C     and replaced by an integer (K).
C
C     This routine will not work properly if NOADC is GT 0.
C
CJDW  1/15/98.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C SG 8/25/98 CENTR changed to character*4
      CHARACTER*4 CENTR(NOCMX)
#ifdef _UNICOS
      CHARACTER*8 ATOMSYM(NOCMX)
#endif /* _UNICOS */
      COMMON /COMINF/ TOTMSS, COM(3), RMOI(3), PATRAN(3,3)
      DIMENSION VLIST(NOCMX,10)
      DIMENSION ASYMB(103), AMASS(103), SPIN(103), QUADN(103)
c      DATA (ASYMB(I),I = 1,103)
c     1/2HH , 2HHe, 2HLi, 2HBe, 2HB , 2HC , 2HN , 2HO , 2HF , 2HNe,
c     2 2HNa, 2HMg, 2HAl, 2HSi, 2HP , 2HS , 2HCl, 2HAr, 2HK , 2HCa,
c     3 2HSc, 2HTi, 2HV , 2HCr, 2HMn, 2HFe, 2HCo, 2HNi, 2HCu, 2HZn,
c     4 2HGa, 2HGe, 2HAs, 2HSe, 2HBr, 2HKr, 2HRb, 2HSr, 2HY , 2HZr,
c     5 2HNb, 2HMo, 2HTc, 2HRu, 2HRh, 2HPd, 2HAg, 2HCd, 2HIn, 2HSn,
c     6 2HSb, 2HTe, 2HI , 2HXe, 2HCs, 2HBa, 2HLa, 2HCe, 2HPr, 2HNd,
c     7 2HPm, 2HSm, 2HEu, 2HGd, 2HTb, 2HDy, 2HHo, 2HEr, 2HTm, 2HYb,
c     8 2HLu, 2HHf, 2HTa, 2HW , 2HRe, 2HOs, 2HIr, 2HPt, 2HAu, 2HHg,
c     9 2HTl, 2HPb, 2HBi, 2HPo, 2HAt, 2HRn, 2HFr, 2HRa, 2HAc, 2HTh,
c     O 2HPa, 2HU , 2HNp, 2HPu, 2HAm, 2HCm, 2HBk, 2HCf, 2HEs, 2HFm,
c     1 2HMd, 2HNo, 2HLr /
      DATA (AMASS(I),I = 1,92)
     1/ 1.0078246, 4.002601, 7.01600, 9.01218, 11.009307,
     2  12.000000, 14.0030738, 15.9949141, 18.9984022, 19.992441,
     3  22.9898, 23.98504, 26.98153, 27.976929, 30.973764,
     4  31.9720727, 34.9688531, 39.962386, 38.96371, 39.96259,
     5  44.95592, 48., 50.9440, 51.9405, 54.9380,
     6  55.9349, 58.9332, 57.9353, 62.9296, 63.9291,
     7  68.9257, 73.9219, 74.9216, 79.9165, 78.91839,
     8  83.91151, 84.9117, 87.9056, 88.9059, 89.9043,
     9  92.9060, 97.9055, 98., 101.9037, 102.9048,
     O  107.90389, 106.90509, 113.9036, 114.9041, 120.,
     1  120.9038, 129.9067, 126.90466, 131.90416, 132.9051,
     2  137.9050, 138.9061, 139.9053, 140.9074, 141.9075,
     3  145., 151.9195, 152.9209, 157.9241, 159.9250,
     4  163.9288, 164.9303, 165.9304, 168.9344, 173.9390,
     5  174.9409, 179.9468, 180.9480, 183.9510, 186.9560,
     6  192., 192.9633, 194.9648, 196.9666, 201.970625,
     7  204.9745, 207.9766, 208.9804, 209., 210.,
     8  222., 223., 226., 227., 232.,
     9  231., 238. /
      DATA (AMASS(I),I = 93,103)
     1/ 237., 244., 243.,
     2  247., 247., 251., 252., 257.,
     3  258., 259., 260./
      DATA (SPIN(I),I = 1,103)
     1/ 103*0./
      DATA (QUADN(I),I = 1,103)
     1/ 103*0./
C
      IF(NOADC .GT. 0)THEN
       WRITE(6,*) '  @INCENT-F, Invalid value of NOADC ',NOADC
       CALL ERREX
      ENDIF
C
C....    READ CENTRE SPECIFICATIONS
C
c      WRITE(6,600) NON,NOADC
600   FORMAT(/' THE NUMBER OF CENTERS=',I5,' ADDITIONAL CENTERS=',I5)
c      WRITE(IOUTU) ILBL
c      WRITE(IOUTU)NON
c      WRITE(6,104)
C      WRITE(6,106)
      NOC=NON+NOADC
      DO 1 I=1,NON
      READ(INU   ) CENTR(I),K,(VLIST(I,J),J=1,4)
      CALL PUTLN(CENTR(I),K)
      IATOM = VLIST(I,4) + 0.5
      IF (IATOM .GT. 0) THEN
cold    VLIST(I,5) = ASYMB(IATOM)
        IF (VLIST(I,6) .LT. 0.) VLIST(I,6) = AMASS(IATOM)
        VLIST(I,7) = SPIN(IATOM)
        VLIST(I,8) = QUADN(IATOM)
cold    VLIST(I,9) = 8H
cold    VLIST(I,10) = 8H
      ELSE
cold    VLIST(I,5) = 2H
        VLIST(I,6) = 0.
        VLIST(I,7) = 0.
        VLIST(I,8) = 0.
cold    VLIST(I,9) = 8H
cold    VLIST(I,10) = 8H
      ENDIF
#ifdef _UNICOS
      CALL INCENT_CH(ATOMSYM,NOCMX,I,IATOM)
#else
      CALL INCENT_CH(VLIST,NOCMX,I,IATOM)
#endif /* _UNICOS */
C
C....    WRITE ON IOUTU NOW PERFORMED IN INPUT
C
C      WRITE(IOUTU) CENTR(I), (VLIST(I,J),J=1,4)
C      WRITE(6,107) CENTR(I), (VLIST(I,J), J=1,4)
1     CONTINUE
      IF(NOADC.LT.1) GO TO 10
c      WRITE(6,114)
c      WRITE(6,115)
      NFIRST = NON + 1
      DO 11 I=NFIRST,NOC
      READ(INU     ) CENTR(I),(VLIST(I,J),J=1,3)
      VLIST(I,4) = 0.0
cold  VLIST(I,5) = 2H
      VLIST(I,6) = 0.
      VLIST(I,7) = 0.
      VLIST(I,8) = 0.
C      WRITE(6,1020) CENTR(I),(VLIST(I,J),J=1,3)
11    CONTINUE
10    CONTINUE
      CALL MOI(VLIST,NOC)
114   FORMAT(//,1X,18HADDITIONAL CENTRES,//)
115   FORMAT(//,1X,6HCENTRE,18X,12HCO-ORDINATES//)
106   FORMAT(/ ,1X,6HCENTRE,18X,12HCO-ORDINATES,17X,6HCHARGE,//)
107   FORMAT(1X,A6,6X,3F12.8,6X,F4.1)
104   FORMAT(/ ,1X,15HNUCLEAR CENTRES)
1020  FORMAT(A6,4X,4F12.0)
      RETURN
      END

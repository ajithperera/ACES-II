      Subroutine Filabl(Lun, Label, Ierr)
C
C Description
C
C   Look for a record label in property integral 
C   program and place the file position pointer 
C   at that label.
C
C Variables
C
C   Lun   = Unit number corresponds to vprop integral file (Input)
C   Label = Record label (input)
C   Ierr  = Error Handler 
C            = 0 Success
C            = 1 unable to find the record label
C Declarations
C
      Implicit Double Precision (A-H, O-Z)
      Character*8 Star, Blank, Check(4), Label
      Parameter (Star  ='********')
      Parameter (Blank ='        ')
C
C Common Block informations
C
      Common/Machsp/Iintln,Ifltln,Iintfp,Ialone,Ibitwd
      Common/Files/Luout,Moints
C
      Ierr = 0
 100  Read(UNIT=Lun, END=1000, ERR=2000, IOSTAT=IOS) Check(1)
      If (Check(1) .ne. Star) Goto 100
C
C Now check the rest to make sure ...
      Backspace (Lun)
      Read(UNIT=Lun, END=1000, ERR=2000, IOSTAT=IOS) Check
      If (Check(4) .eq. Label) Return
      Goto 100
C
C Handle Errors
C
 1000 Ierr = 1
      Return
C     
 2000 Continue

C Get here via an error reading the file
C
      Write(luout,*) ' @Filabl-I, Read error on file attched to
     &                unit ', Lun, '.'

      Return
      End

      subroutine Opnfil(Lun, Filnam, Format, Iopned, Ierr)

C Description
C   Open files attach to a particular unit number
C 
C Variables
C   Lun     : Unit number corresponding to Vprop integral file (input)
C   Filnam  : File name for property integrals
C   Format  : Structure of the file (input)
C   Iopned  : Check to see the file is being used(input) 
C              = 0 then close it and continue
C              = 1 Print a error message and abort
C   Ierr    : Return code
C              = 0 Succesful completion
C              = 1 unable to open/read file
C
C   Ioerr   : Holds return from IOSTAT
C
C Declarations
      Implicit Double Precision(A-H, O-Z)
      Character*(*) Filnam, Format
      Logical Isopen
C
C Common block informations
C
      Common/Machsp/Iintln,Ifltln,Iintfp,Ialone,Ibitwd
      Common/Files/Luout,Moints
     
C See whether particular file exists
C
      Inquire(UNIT=Lun, OPENED=Isopen)
      If (Isopen) then 
         Ierr = 2 + Iopned
         If (Iopned .eq. 1) Return
         Close(UNIT=Lun, IOSTAT=Ioerr)
         If (Ioerr .ne. 0) then
            Ierr = 5
            Write(Luout,*) '@Opnfil-I I/O error on Unit ', Lun,'.'
            Return
         Endif
      Endif
         
C It is not opend, good deal
C
      Open(UNIT=Lun, FILE=Filnam, FORM=Format, Access='SEQUENTIAL',
     &     STATUS='UNKNOWN',IOSTAT=Ioerr) 
C
      If (Ioerr .ne. 0) then
         Ierr = 7
         Write(Luout,*) '@OPnfil-I I/O error on Unit ', Lun,'.'
         Return
      Endif
C
      Rewind (Lun)
C
      Return
      End


      

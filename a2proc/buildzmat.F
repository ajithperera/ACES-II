#include "iachar.h"
      Subroutine Buildzmat(Plusminus,Coords,Label,Nreals,
     +                     Imode,Ipoints)
      Implicit Double Precision (a-h,o-z)

#include "linelen.par"
#include "io_units.par"
#include "fnamelen.par"
   
      Integer Fnblnk,Linblnk,F_strpbrk
      Character*7 file_name0
      Character*7 file_name1
      Character*1 blank
      Character*5 Label(Nreals)
      Character*5 Plusminus
      Dimension Coords(3*Nreals)
      Character*2 V1
      Character*2 V2
      Character*1 V3
      Character*1 czPercent,czAsterisk,czHash,czSpace,czFirstNLChar
      Character*1 czTab
      Character*1 Achar 
      Character*10 String
      Character*(Linelen) Wrk,Sztmp,Wrk2
      Character*(linelen) Trmblk
      Logical Btmp,Bopened,Bdone
      Logical First,Last,bDelLastChar

      czAsterisk = achar(_IACHAR_ASTERISK)
      czPercent  = achar(_IACHAR_PERCENT)
      czSpace    = achar(_IACHAR_SPACE)
      czHash     = achar(_IACHAR_HASH)
      czTab      = achar(_IACHAR_TAB)

      iunit = 10
      blank = ""
 
      If (Imode .Gt. 99) Then
         Write(6,"(a)") " The maximum number of modes can not be",
     +                  "  higher than 100."
         Write(6,"(a)") " This artificial limitation can be easily",
     +                  " removed if needed."
         Call Errex 
      Endif 
                        
      Write(V1,"(i0)") Ipoints 
      Write(V2,"(i0)") Imode
      If (Plusminus .Eq. "Plus ") V3 = "p"
      If (Plusminus .Eq. "Minus") V3 = "m"

      String = V3//"zmat_"//Trim(V2)//Trim(V1)

      Inquire(file=ZFil,opened=bOpened)
      If (.Not.bopened) then
         Open(Unit=LUZ,file=ZFil,form='FORMATTED',status='OLD')
      End if

      rewind(luz)
      Open(Unit=Iunit,File=String,Form="Formatted")

C Skip to the *NEW_INPUT namelist
      bTmp = .true.
      Do while (bTmp)
         If (wrk(i:(i+9)).eq.'*NEW_INPUT') Then
            bTmp = .false.
         Else
            i = 0
            do while (i.eq.0)
               read(luz,'(a)',err=8000,end=5400) wrk
               i = fnblnk(wrk)
               if (i.ne.0) then
                  if (wrk(i:i).ne.czAsterisk) i = 0
               end if
            end do
         End if
      End do

C count the number of lines with keywords on them
c There are three ways to terminate the NEW_INPUT namelist:
c - there is an unpaired close parenthesis, ')'
c - an asterisk is the first non-blank character on the line
c - the end-of-file is reached

      czFirstNLChar = czSpace
      bDelLastChar  = .True.
      First         = .True. 
      Last          = .False. 

      ilines = 1
      Call Count_lines(Ilines,Wrk)

      do i = 1, iLines
         backspace(luz)
      end do

      Do iLine = 1, 2
         read(luz,'(a)') wrk

         iLastC = index(wrk,czHash)
         if (iLastC.eq.1) then
            wrk = ' '
         Else
            If (iLastC.eq.0) then
               wrk2 = trmblk(wrk)
            Else
               wrk2 = trmblk(wrk(1:iLastC-1))
            End if
            If (iLine.eq.1) then
c            o throw out the name of the namelist
               wrk  = wrk2(11:)
               wrk2 = trmblk(wrk)
            End if
         Endif
         If (Iline .Eq. 2) Write(Iunit,"(a5)") Wrk
      Enddo

      Do I = 1, Nreals 
         Ioff = 3*(I-1) + 1
         Write(Iunit,"(a5,3(1x,F15.10))") Label(i),(Coords(Ioff+J-1),
     +         J =1, 3)
      Enddo

      Write(iunit,"(a1)") czSpace

      Do iLine = 3, iLines
         read(luz,'(a)') wrk

         iLastC = index(wrk,czHash)
         if (iLastC.eq.1) then
            wrk = ' '
         Else
            If (iLastC.eq.0) then
               wrk2 = trmblk(wrk)
            Else
               wrk2 = trmblk(wrk(1:iLastC-1))
            End if
         Endif 

         I = Index(wrk2,czSpace)
         If (Iline .Eq.3) Then 
            Write(Iunit,"(a)") "*ACES2("//Wrk2(1:I)
         Else
            If (Wrk2(2:2) .Eq. ":" .OR. Wrk2(3:3) .Eq. ":") Then
                If (First) Write(iunit,"(a1)") czSpace
                First = .False. 
            Endif 
            If (.Not. First) Then
               Write(Iunit,"(a)") Wrk2(1:I-1)
            Else
               Write(Iunit,"(a)") Wrk2(1:I-1)//")"
            Endif 
         Endif 
      Enddo

      Write(iunit,"(a1)") czSpace

      close(iunit)
      Return

 8000 Write(6,"(a)") '@-buildzmat: There was a problem reading ZMAT.'
      Call Errex

 5400 Write(6,"(2a)") '@-buildzmat: ZMAT is missing the NEW_INPUT',
     +                ' namelist.'
      close(unit=luz,status='KEEP')
      call errex

      Return
      End

#line 7 "ishell.F"
      function ishell(string)
      character*(*) string
      integer system
      external system
      character*80 cmd
      if ( len(string) .lt. 80 ) then
         cmd = string // char(0)
         ishell = system(cmd)
      else
         write(*,*)
     &      '@ISHELL: The command buffer is only 80 characters long.'
         write(*,*)
     &      '         calling errex while trying to execute: ',
     &      '(string between arrows)'
         write(*,'(3a)') '-->',string,'<--'
         call errex
      end if
      return
      end

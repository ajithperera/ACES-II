
#include "flags.h"

      subroutine parse_cli
      implicit none

      integer iLast
      integer i, f_iargc, num_args
      character arg*80
      logical bDoNext, ACESIII

      character*1 achar
      integer  f_strpbrk, c_atol
      external f_strpbrk, c_atol

#include "parallel_aces.com"

c ----------------------------------------------------------------------

      num_args = f_iargc()

c   o loop over arguments
      bDoNext = .true.
      do i = 1, num_args

         if (bDoNext) then
c         o process this argument

            arg=' '
            call f_getarg(i,arg)

c         o MPI_COMM_SIZE
            if (arg(1:6).eq.'-procs') then
               if (i+1.le.num_args) then
                  call f_getarg(i+1,arg)
                  iLast = f_strpbrk(arg,' ')
                  arg(iLast:iLast) = achar(0)
                  nprocs = c_atol(arg)
                  bDoNext = .false.
                  print *, 'nprocs = ',nprocs
               else
                  print *, '@PARSE_CLI: -procs is missing an argument'
                  call aces_exit(1)
               end if
            end if

c         o MPI_COMM_RANK
            if (arg(1:5).eq.'-rank') then
               if (i+1.le.num_args) then
                  call f_getarg(i+1,arg)
                  iLast = f_strpbrk(arg,' ')
                  arg(iLast:iLast) = achar(0)
                  irank = c_atol(arg)
                  bDoNext = .false.
                  print *, 'irank = ',irank
               else
                  print *, '@PARSE_CLI: -rank is missing an argument'
                  call aces_exit(1)
               end if
            end if

c        else if (.not.bDoNext) then
         else

c         o resume processing arguments
            bDoNext = .true.

c        end if (bDoNext)
         end if

c     end do i = 1, num_args
      end do
C
C This is useful for ACESIII. Test for this flag to terminate after
C gradients are done (in ACES II, indiviudal MEs are called in 
C sequence for gradient calculation), Ajith Perera, 04/2012.
C
      ACESIII = .FALSE. 
      If (Nprocs .Ne. 0 .And. Irank .Ne. 0 .And. ACESIII) 
     &    Call Putrec(20,"JOBARC", "NFDMANUL", 1, 1)

c ----------------------------------------------------------------------

c   o check consistency

c   o MPI_COMM_SIZE
      if (nprocs.lt.1) then
         print *, '@PARSE_CLI: resetting number of processes to 1'
         nprocs = 1
      end if

c   o MPI_COMM_RANK
      if (irank.lt.0) then
         print *, '@PARSE_CLI: resetting process rank to 0'
         irank = 0
      end if
      if (irank.ge.nprocs) then
         print *, '@PARSE_CLI: resetting process rank to ',nprocs-1
         irank = nprocs-1
      end if

c ----------------------------------------------------------------------

      return
c     end subroutine parse_cli
      end


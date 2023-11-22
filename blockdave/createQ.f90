      subroutine createQ(colsV, Nrows,Ncols)
      integer,intent(in):: Nrows,Ncols
      double precision,intent(inout) :: colsV(Nrows,Ncols)
      double precision:: tau(Ncols),work(Ncols)
      integer::lwork,info

!double precision:: Q(Nrows,Nrows),H(Nrows,Nrows)
!double precision:: Hscr(Nrows,Nrows)!,R(Ncols,Ncols)
!integer::i,j,offset
real::start,finish
! colsV - the matrix needing to be QR factorized
! Nrows - rows of colsV
! Ncols - columns of colsV
! H - Current iterations Householder reflector
! Hscr - scatch space to hold the matrix multiplication of H1*H2*...*Hm
! Q - final form of the orthogonal eigenvectors of colsV (Nrows x Nrows)

! Purpose:: Computes a reduced QR factorization of colsV and stored eigenvectors
!           back into colsV(Nrows, Ncols) 

! First, perform MLK QR decomposition
! Then, call supporting routine to build the Q 
! From Household reflectors.
! Def. of LAPACK routine:::
! https://www.netlib.org/lapack/explore-html/dd/d9a/group__double_g_ecomputational_ga3766ea903391b5cf9008132f7440ec7b.html
! Def. of dorgqr:
! https://software.intel.com/content/www/us/en/develop/documentation/onemkl-developer-reference-fortran/top/lapack-routines/lapack-least-squares-and-eigenvalue-problem-routines/lapack-least-squares-and-eigenvalue-problem-computational-routines/orthogonal-factorizations-lapack-computational-routines/orgqr.html#orgqr


      tau=0.0d0
      work=0.0d0
      lwork=Ncols
      print*,'Ncols inside QR:', Ncols
      call cpu_time(start)
      call dgeqrf(Nrows,Ncols,colsV,Nrows,tau,work,lwork,info)
      call dorgqr(Nrows,Ncols,Ncols,colsV,Nrows,tau,work,lwork,info) 

      call cpu_time(finish)
      print*,'Time spent in LAPACK dgeqrf routine:', finish-start

! ** THE FOLLOWING IS OBSOLETE: SLOW WAY TO COMPUTE Q FROM
!    HOUSEHOLDER REFLECTORS
!
!print*,'cols V:', colsV
! create H1
!Q=0.0d0
!H=0.0d0
!Hscr=0.0d0
!offset=Nrows-2+1
!call cpu_time(start)
!call createH(colsV(2:Nrows,1),offset,tau(1),1,Nrows,H)
!  call cpu_time(finish)
!  print*,'Time spent in 1 call to createH routine:', finish-start
!!print*,'H1'
!!print*,H
!Hscr=H
!! create H2,...,Hm, while matrix multiplying to create the final
!! orthogonal Q
!call cpu_time(start)
!do i=2,Nrows
!  offset=Nrows-(i+1)+1
!  H=0.0d0
!  !print*,'Hscr'
!  !print*,Hscr
!  if (i.ne.Nrows) then
!    call createH(ColsV(i+1:Nrows,i),offset,tau(i),i,Nrows,H)
!! Handles the overcount of ColsV; for square matrix this H_nrows returns I
!  else 
!    call createH(ColsV(i:Nrows,i),offset,tau(i),i,Nrows,H)
!  endif
!  !print*,'H inside createQ:'
!  !print*,H
!  call dgemm('N','N',Nrows,Nrows,Nrows,1.0d0,Hscr,Nrows,H,Nrows,0.0d0,Q,Nrows)
!  !print*,'RESULTS OF MATRIX MULT:'
!  !print*,Q
!  Hscr=Q
!
!enddo
!  call cpu_time(finish)
!  print*,'Time spent in lOOP of QR routine:', finish-start

!colsV=Q(1:Nrows,1:Ncols)
   
   Return 
   End

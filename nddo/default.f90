subroutine default(ikey)
! sets keyword defaults
use tables
use control
interface
subroutine atoi(input,j)
character,intent(in)::input
integer,intent(inout)::j
end subroutine atoi
end interface
integer,intent(in)::ikey
integer::iexp
character(20)::ch

scftol=1D-6
METHOD='AM1'
DEORTHO=.false.
GRADIENT=.false.
TRIAL=.false.
XYZ=.false.
OPTIMIZE=.false.
OPTTOL=1.0d0
parameters=.false.
core='standard'
CENTERS=.false.
newton=.false.
restricted=.true.
MULT=1
sys_charge=0.0d0
DFTHF=.false.
LAPACK=.FALSE.
SPARKLES=.FALSE.
densityin=.false.
densityout=.false.
FREQUENCY=.FALSE.

do i=1,ikey
ch=' '
ch=keyword(i)


  if(ch=='AM1')then
  METHOD='AM1'
  core='GAUSSIAN'
  elseif(ch=='PM3')then
  METHOD='PM3'
  core='GAUSSIAN'
  elseif(ch=='MNDO')then
  METHOD='MNDO'
  elseif(ch=='MNDOD')then
  METHOD='MNDOD'

  elseif(ch(1:4)=='CONV')then
  call atoi(ch(6:),iexp)
  scftol=10**(-float(iexp))

 
 elseif(ch=='DEORTHO')then
  DEORTHO=.true.
 elseif(ch=='GRADIENT')then
  GRADIENT=.true.
 elseif(ch=='TRIAL')then
  TRIAL=.true.
 elseif(ch=='XYZ')then
  XYZ=.true.
  elseif(ch=='OPTIMIZE')then
  OPTIMIZE=.true.

elseif(ch(1:6)=='OPTTOL')then
  call atoi(ch(8:),iexp)
  OPTTOL=10**(-float(iexp))
elseif(ch(1:6)=='CHARGE')then
  if(ch(8:8)=='-')then
     call atoi(ch(9:),iexp)
     SYS_CHARGE=-float(iexp)
  else
     call atoi(ch(8:),iexp)
     SYS_CHARGE=float(iexp)
     end if

elseif(ch(1:4)=='MULT')then
  call atoi(ch(6:),iexp)
  mult=int(iexp)

elseif(ch=='PARAMETERS')then
  parameters=.true.
elseif(ch=='GAUSSIAN')then
  core='GAUSSIAN'

elseif(ch=='3CENTER')then
  CENTERS=.true.
elseif(ch=='NEWTON')then
  newton=.true.
elseif(ch=='UHF')then
  restricted=.false.
elseif(ch=='HFDFT')then
  DFTHF=.true.
elseif(ch=='LAPACK')then
  LAPACK=.TRUE.
elseif(ch=='SPARKLES')then
  SPARKLES=.TRUE.
elseif(ch=='DENSITYIN')then
  densityin=.TRUE.
elseif(ch=='DENSITYOUT')then
  densityout=.TRUE.
elseif(ch=='FREQUENCY')then
  FREQUENCY=.TRUE.
 end if



end do
end subroutine default




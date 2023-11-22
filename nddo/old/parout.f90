subroutine parout
use tables
use constants
use control
implicit double precision (a-h,o-z)
write(*,*)''
write(*,*)'PARAMETERS USED IN THIS CALCULATION:'
WRITE(*,*)''
do i=1,itype
write(*,*)'PARAMETERS FOR ',PERIODIC(NTYPE(I)),':'

write(*,10)'USS  ',USS(I), ' S ORBITAL ENERGY (eV)'
write(*,10)'UPP  ',UPP(I), ' P ORBITAL ENERGY (eV)'
if(NBAS(I)==9)write(*,10)'UDD  ',UDD(I), ' D ORBITAL ENERGY (eV)'
write(*,10)'BETAS  ',BETAS(I),' S BONDING PARAMETER (eV)'
write(*,10)'BETAP  ',BETAP(I),' P BONDING PARAMETER (eV)'
if(NBAS(I)==9)write(*,10)'BETAD  ',BETAD(I),' D BONDING PARAMETER (eV)'
write(*,10)'ZS  ',ZS(I),' S ORBITAL EXPONENT (au)'
write(*,10)'ZP  ',ZP(I),' P ORBITAL EXPONENT (au)'
if(NBAS(I)==9)write(*,10)'ZD  ',ZETAD(I),' D ORBITAL EXPONENT (au)'
write(*,10)'GSS  ',GSS(I),' (SS,SS) (eV)'
write(*,10)'GSP  ',GSP(I),' (SS,PP) (eV)'
write(*,10)'GPP  ',GPP(I),' (PP,PP) (eV)'
write(*,10)'GP2  ',GP2(I),' (PP,P*P*) (eV)'
write(*,10)'HSP  ',HSP(I),' (SP,SP) (eV)'
write(*,10)'HPP  ',HALF*(GPP(I)-GP2(I)),' (PP*,PP*) (eV)'
write(*,11)'NBAS  ',NBAS(I),' NUMBER OF BASIS ORBITALS'
write(*,10)'ALPHA  ',ALPHA(I),' ALPHA CORE PARAMETER (1/Ang)'
if(core=='GAUSSIAN' .or.method=='AM1'  .or. method=='PM3')then
write(*,10)'A1  ',G1(1,I),' PRE-EXPONENTIAL -  CORE TERM'
write(*,10)'B1  ',G2(1,I),' GAUSSIAN EXPONENT - CORE TERM'
write(*,10)'C1  ',G3(1,I),' GAUSSIAN WIDTH - CORE TERM'
write(*,10)'A2  ',G1(2,I),' PRE-EXPONENTIAL -  CORE TERM'
write(*,10)'B2  ',G2(2,I),' GAUSSIAN EXPONENT - CORE TERM'
write(*,10)'C2  ',G3(2,I),' GAUSSIAN WIDTH - CORE TERM'
write(*,10)'A3  ',G1(3,I),' PRE-EXPONENTIAL -  CORE TERM'
write(*,10)'B3  ',G2(3,I),' GAUSSIAN EXPONENT - CORE TERM'
write(*,10)'C3  ',G3(3,I),' GAUSSIAN WIDTH - CORE TERM'
write(*,10)'A4  ',G1(4,I),' PRE-EXPONENTIAL -  CORE TERM'
write(*,10)'B4  ',G2(4,I),' GAUSSIAN EXPONENT - CORE TERM'
write(*,10)'C4  ',G3(4,I),' GAUSSIAN WIDTH - CORE TERM'
end if
write(*,10)'-----------------------------------'
end do
10 format(A10,F20.5,A40)
11 format(A10,I20,A40)
end subroutine parout

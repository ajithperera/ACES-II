        subroutine Read_rotinfo(Type,UorL,Eps,Iop1,Iop2,Iop3,
     +                          Iop4)

        Implicit Double Precision(A-H,O-Z)
c
c read in orbital rotation parameters, Fermi level is designated as 
C o and occupied is counted as -1,-2... and the virtual are labelled
C as 1, 2....
c
        Character*2 Type 
        Character*1 UorL
	open(unit=1,file='porbrot.info',
     &     form='formatted',access='sequential',status='old')

        read(1, "(a)") Type
        read(1, "(a)") UorL
	read(1,*) Eps
        read(1,*) Iop1
        read(1,*) Iop2
        read(1,*) Iop3
        read(1,*) Iop4
        close(1)
        return
        end

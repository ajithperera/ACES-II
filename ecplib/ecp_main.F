#include "flags.h"

      Subroutine ECP_MAIN(Iecp,Natoms,Namat,Iqmstr,Jcostr,Nucstr,
     &                    Nrcstr)
C      
      Implicit Double Precision(A-H, O-Z)

#include "ecp.par"
#include "ecp.com"
#include "flags.com"
#include "flags2.com"

      Character*2 mtype(Mxatms)
      Character*80 Cnick(Mxatms)
      Character*80 Ecpnam(Mxatms)
      Character*4 namat(Mxatms)
      Logical Nharm, Grads
C
      Dimension iqmstr(Mxatms),jcostr(Mxatms,Maxang), 
     &          nucstr(Mxatms,Maxang,Mxprim),
     &          nrcstr(Mxatms,Maxang,Mxprim),IGenby(Mxatms),
     &          Coord(3,Mxatms)
C
      Dimension  Cint(Maxmem)
C
      COMMON /INDX/ PC(512),DSTRT(8,MXCBF),NTAP,LU2,NRSS,NUCZ,ITAG,
     & MAXLOP,MAXLOT,KMAX,NMAX,KHKT(7),MULT(8),ISYTYP(3),ITYPE(7,28),
     & AND(8,8),OR(8,8),EOR(8,8),NPARSU(8),NPAR(8),MULNUC(Mxatms),
     & NHKT(MXTNSH),MUL(MXTNSH),NUCO(MXTNSH),NRCO(MXTNSH),JSTRT(MXTNSH),
     & NSTRT(MXTNSH),MST(MXTNSH),JRS(MXTNSH)
C
      COMMON /DAT/  EXPA(MXTNPR),CONT(MXTNCC),CENT(3,MXTNSH),
     &              CORD(Mxatms,3),CHARGE(Mxatms),FMULT(8),TLA, TLC
C
C ECP integrals. Except for simplifications, I have kept the original
C form of all the input/output ECP processing routines. This was 
C originaly written by Christien Huber, Ajith Perera, 12/2001.
C
C - Natoms is the number of symmetry unique atoms
C - Ipseux is 0 when no ECP definitions; > 0 otherwise.
C - Cnick is the name of then ECP 
C - iqmstr highest ang. mom. qun. num.(AMQN) for each atom.
C - jcostr has the AQMN for each shell (each atom). Normally the # of shells
C   is the same as the Max AMQN but some cases there may be mutiple shells
C   with the same AQMN, so jcostr tells you how many.
C - nucstr is the number of primitives per each shell (per atom,per AQMN).
C - nrcstr is the number of contracted functions per each shell
C   (per atom,per AQMN)
C - nhram is logical that tells you whether this is sperical or 
C   cartesian basis calculation (key-word driven)
C - I must add that the names choosen for those 4 variables are dumb!!
C   
      nharm = (Iflags(h_IFLAGS_spherical) .EQ. 1)
      Itol  = Iflags(h_IFLAGS2_intgrl_tol)
      Thres = Dble(10.0D+00*(-Itol))
 
      if (iecp.ne.0) then
c-----------------------------------------------------------------------
c get atomsymbols and put them to mtype-array
c-----------------------------------------------------------------------
         ierr=0
         call chrges(CHARGE,mtype,natoms,ierr)
c-----------------------------------------------------------------------
c read ecp data and prepare data for evulation of the c integrals
c-----------------------------------------------------------------------

         if (ierr.eq.0) then
             ipr=0
             call ecpdef(natoms,mtype,Ecpnam,cnick)

             if (ierr.eq.0) then
                ipr=0
                call ecppar(ipr,cnick,CHARGE,natoms,.FALSE.,0)
C
c-----------------------------------------------------------------------
c  Summation over the charges and storage of this value
c-----------------------------------------------------------------------

                chgsum=0
                do i=1,natoms
C
C account for symmetry 
C
                    chgsum=chgsum+charge(i)*fmult(mulnuc(i))
                Enddo

                nproton=idint(chgsum)
c
                call putrec(20,'JOBARC','NMPROTON',1,nproton)

c-----------------------------------------------------------------------
c prepare data for ecp integral evaluation
c-----------------------------------------------------------------------
                Call ecp_init
                Call ecp_basis_init(namat,iqmstr,jcostr,nucstr,
     &                              nrcstr, NHARM)
          
c-----------------------------------------------------------------------
             else
                write(6,"(a)") 'ERROR while reading ECP data!'
                write(6,"(a)") 'Please check yourECPDATA file.'
                Call Errex
             endif

         else
             write(6,"(a)") 'E R R O R in subroutine chrges!'
             Call Errex
         endif
c-----------------------------------------------------------------------

      Endif

      Call Dzero(Cint, Maxmem)
      Call Getrec(20, 'JOBARC', 'NATOMS  ', 1, Ntotatoms) 
      Call Getrec(20, 'JOBARC', 'COORD   ', 3*Ntotatoms, Coord)

      Call Symeqv(Ntotatoms, IGenby)
#ifdef _DEBUG_LVL0
      Write(6,"(a)") "The symmetry equivalent atoms"
      Write(6,"(6(1x,I3))") (IGenby(i), i=1, Ntotatoms)
      Write(6,*)
      Do j=1, Ntotatoms
      Write(6, "(3(1x,F10.6))") (Coord(i,j), i=1, 3)
      Enddo
      Write(6,*) 
#endif
      Call ECP_int_driver(Natoms, Ntotatoms, IGenby, Coord, Cint, 
     &                    .False.)
      
      Return
      End

#include "flags.h"

      Subroutine ECP_set4_grdint(Iecp,Natoms, Namat, Iqmstr,
     &                           Jcostr, Nucstr,Nrcstr)
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
     &          nrcstr(Mxatms,Maxang,Mxprim)
C
      Dimension  Cint(Maxmem)
C
      COMMON /PRIMIT/ PRIEXP(MXPRIM), PRICCF(MXPRIM,MXCONT),
     &                PRICRX(MXPRIM), PRICRY(MXPRIM), PRICRZ(MXPRIM)

      COMMON /NUCLEI/ CHARGE(MXATMS), CORD(MXATMS,3),
     &                DCORD(MXATMS,3),DCORGD(MXATMS,3),
     &                DOPERT(0:3*MXATMS)

      COMMON /SYMMET/ FMULT(0:7), PT(0:7),
     &                MAXLOP, MAXLOT, MULT(0:7), ISYTYP(3),
     &                ITYPE(4,10), NPARSU(8), NPAR(8), NAOS(8),
     &                NPARNU(8,8), IPTSYM(MXCORB,0:7),
     &                IPTCNT(3*MXATMS,0:7), NCRREP(0:7),
     &                IPTCOR(MXATMS*3), NAXREP(0:7), IPTAX(3),
     &                IPTXYZ(3,0:7)
C
C ECP deriavtive integrals. Except for simplifications, I have kept the 
C original  form of all the input/output ECP processing routines. This was 
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
             iunit=55
             call ecpdef(natoms,mtype,Ecpnam,cnick)

             if (ierr.eq.0) then
                ipr=0
                call ecppar(ipr,cnick,CHARGE,natoms,.FALSE.,0)
C
c-----------------------------------------------------------------------
c prepare data for ecp integral evaluation
c-----------------------------------------------------------------------
                Call ecp_init
                Call ecp_basis_init_4grd(namat,iqmstr,jcostr,nucstr,
     &                                   nrcstr,.FALSE., Natoms)
          
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

      Return
      End

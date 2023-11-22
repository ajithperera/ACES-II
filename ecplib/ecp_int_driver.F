      subroutine Ecp_int_driver(Natoms, Ntotatoms, IGenby, Coord, Cint,
     &                          Grads)
c----------------------------------------------------------------------
c     Ecp_int_driver drives the calculation of 1e-matrices in the basis of
c     cartesian atomic orbitals using subroutine ecp_int.
c----------------------------------------------------------------------
c     tol    = integral cutoff (ints to be neglected, if value expected
c              to be < exp(-tol)(i)
c     Cint   = 1e-integral array ('batch') of dimension nfij
c       
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
C-----------------------------------------------------------------------
#include "ecp.par"
#include "ecp.com"

CSSS      parameter (nftmax = (ndilmx+1)*(ndilmx+2)*(ndilmx*3)/6)
c     symmetry parameter:
C-----------------------------------------------------------------------
C     nftmax  total number of reducible gaussians up to {ndilmx+1}
C     SPECIAL parameter needed in gradient calculations
C-----------------------------------------------------------------------

      logical iandj,forget,Zero_int, Grads
C 
C Ecpints need to be (2*Maxang+1)**2 times the square of the 
C maximum number of contracted functions per shell (currently
C set at 20 for ECP integrals). The current Mxcbf is 1000. 

      dimension Cint(Maxmem), Ecpint_4Shell(Maxints_4shell),
     &          Ecpint(Mxcbf*(Mxcbf+1)/2), IGenby(Ntotatoms),
     &          Coord(3,Ntotatoms)

CSSS      common /infoa / xyz(3,Mxatms),charg(Mxatms),wmass(Mxatms),natoms
C
CSSS      common /symshe/ mulsh(Mxtnsh)
CSSS      common /modez / zetm(Mxtnsh)
CSSS      common /pairij/ ipq(ndi4+1)
CSSS      common /typinf/ ltmax,ltdrv,nftnft,nftdrv
CSSS      common /powers/ jx(nftmax),jy(nftmax),jz(nftmax),
CSSS     &                ix(nftmax),iy(nftmax),iz(nftmax)
CSSS      common/xa/xand(3,8,4*(Maxang-1)-3)
CSSS      common/symsh2/mij,kij,hkij,kb
C
      integer and,eor,or,dstrt
C
c-----------------------------------------------------------------------
c     common from AcesII (created in Readin)
c-----------------------------------------------------------------------
      COMMON /DAT/ EXPA(MXTNPR),CONT(MXTNCC),CENT(3,MXTNSH),
     &             CORD(Mxatms,3),CHARGE(Mxatms),FMULT(8),TLA,TLC

      common/indx/pc(512),dstrt(8,MXCBF),ntap,lu2,nrss,nucz,itag,
     &            maxlop,maxlot,kmax,nmax,khkt(7),mult(8),isytyp(3),
     &            itype(7,28),and(8,8),or(8,8),eor(8,8),nparsu(8),
     &            npar(8),mulnuc(Mxatms),nhkt(MXTNSH),mul(MXTNSH),
     &            nuco(MXTNSH),nrco(MXTNSH),jstrt(MXTNSH),
     &            nstrt(MXTNSH),mst(MXTNSH),jrs(MXTNSH)
c-----------------------------------------------------------------------
#ifdef _DEBUG_LVL0
CSS      Nshell = 2
      Write(6,"(a,I3)") "The Debug value for Nshell: ", Nshell
CSS      nprims(2) = 1
#endif

      Iloc   = 0
      Indecp = 0
      Do Icent = 1, Ntotatoms
      do ishell = 1, Nshell
CSSS      do ishell = 6, 6
        iatom  = katom(ishell)

        If (IGenby(Icent) .EQ. Iatom) Then

        xa = Coord(1,icent)
        ya = Coord(2,icent)
        za = Coord(3,icent)
        la = ktype(ishell) - 1

        kprimi=kprim(ishell)-1
        idegen=ndegen(ishell)
        numcoi=numcon(ishell)
        imin  =kmini(ishell)
        imax  =kmaxi(ishell)

C begin loop jshell 

        Jloc = 0
        Do Jcent = 1, Icent
        do jshell=1, Nshell
CSSS        do jshell=5, 5

          Jdegen=ndegen(Jshell)
          jatom=katom(jshell)
          lb = ktype(jshell) - 1
          jmin = kmini(Jshell)
          jmax = kmaxi(Jshell)
C
#ifdef _DEBUG_LVLM
      Write(6,"(a)") "Icent, Jcent, Jatom, Ishell, Jshell: " 
      write(6,"(6(1x,i4))") Icent, Jcent, Jatom, Ishell, Jshell, Jloc
      Write(6,*)
#endif
          If (IGenby(Jcent) .EQ. Jatom) Then
          If (.NOT. (ICent .EQ. Jcent .AND. JShell .GT. IShell)) Then
C
          xb=Coord(1,jcent)
          yb=Coord(2,jcent)
          zb=Coord(3,jcent)

          kprimj=kprim(jshell)-1
          numcoj=numcon(jshell)

#ifdef _DEBUG_LVLM
      Nprims(1) = 1
      Write(6,"(a,I3)") "The Debug value for Nprims: ", Nprims(1)
#endif
          Call Dzero(Ecpint_4shell, Maxints_4shell)
C
          Do iprim = 1, nprims(ishell)
CSSS                  Do iprim = 3, 3

              exp1 = expnt(kprim(ishell)+iprim-1)
              Indnpri = Kprimi + Iprim
                     
              Do jprim = 1, nprims(jshell)
CSSS                     Do jprim = 3, 3
                 exp2 = expnt(kprim(jshell)+jprim-1)
                 Indnprj = Kprimj + Jprim
                 Call Dzero(Cint, Maxmem)
C
#ifdef _DEBUG_LVLM
      Write(6, "(a, 4(1x,I2))") "The Ishell, Jshell, La, Lb: ", 
     &                           Ishell, Jshell, La, Lb
      Write(6, "(a, 4(1x,I2))") "The prim. pair; Iprim, Jprim",
     &                           Iprim, Jprim
#endif
#ifdef _DEBUG_LVLM
      Write(6, "(a, 3(1x,F10.6))") "The XYZ of A : ", Xa, ya, Za
      Write(6, "(a, 3(1x,F10.6))") "The XYZ of B : ", Xb, yb, Zb
      Write(6, "(a, 2(1x,F10.6))") "The Primtive pair: ", Exp1, Exp2
      Write(6,*)
      Write(6,*) "Initilize debug data"
#endif 
#ifdef _DEBUG_LVLM
      Write(6,"(a)") "Initilize debug xyz of A and B"
      xa = 1.0
      ya = 1.0
      za = 1.0
      xb = -1.0
      yb = -1.0
      zb = -1.0
      exp1 = 1.5
      exp2 = 1.5 
      La = 2
      Lb = 2
      Write(6, "(a, 3(1x,F10.6))") "The XYZ of A : ", Xa, ya, Za
      Write(6, "(a, 3(1x,F10.6))") "The XYZ of B : ", Xb, yb, Zb
      Write(6, "(a, 2(1x,F10.6))") "The Primtive pair: ", Exp1, Exp2
      Write(6,*)
#endif
                 Call Ecp_int_4prim(Xa, Ya, Za, Xb, Yb, Zb, 
     &                              La, Lb, Coord, Exp1, Exp2,
     &                              Natoms, Ntotatoms, IGenby,
     &                              Cint, Int, Zero_int, 
     &                              Grads)

                 If (.NOT. Zero_int)
     &              Call ecp_int_4shell(Cint, Ecpint_4shell, 
     &                                  La, Lb, Numcoi, Numcoj, 
     &                                  Indnpri, Indnprj, Iprim,
     &                                  Jprim, Jnt)
             Enddo
          Enddo        
C 
#ifdef _DEBUG_LVL0
      Write(6,*)
      Write(6,"(a,4(1x,I2))")
     &"The Contracted integral for shell pair", Ishell, JShell, Idegen,
     & Jdegen

      If (Ishell .NE. Jshell) Then
          Do I = 1, Idegen
             Write(6,"(4(1x,F20.13))")(Ecpint_4shell((I-1)*Jdegen+J),
     &                                 J=1,Jdegen)
          Enddo
      Else
          Do I = 1, Idegen
             Write(6,"(4(1x,F20.13))")(Ecpint_4shell((I-1)*Jdegen+J),
     &                                 J=1,I)
          Enddo
      Endif
#endif  
                  Call Ecp_int(Ecpint, Ecpint_4shell, La, Lb, Iloc, 
     &                         Jloc, Numcoi, Numcoj, Imin, Imax, 
     &                         Jmin, Jmax, Indecp)
C
          Jloc = Jloc + (Lb+1)*(Lb+2)*Numcoj/2
C
        Endif
        Endif
C
        Enddo
        Enddo

        Iloc = Iloc + (La+1)*(La+2)*Numcoi/2
C
      Endif 
      Enddo 
      Enddo
C 
#ifdef _DEBUG_LVL0
      Write(6,*)  Iloc*(Iloc+1)/2
      w=0.0d0
      Write(6,"(a)") "@-Ecpint_Driver, The ECP integrals"
      Write(6, "(6(1x,F10.7))") (Ecpint(I), I=1, Iloc*(Iloc+1)/2)
      do i =1, iloc*(iloc+1)/2
      W = w+ecpint(i)*ecpint(i)
      enddo
      Write(6, "(a,(1x,F20.13))") "The int. check sum = ", w
#endif
C
      Call Putrec(20, "JOBARC", "ECP1INTS", iloc*(iloc+1)/2, Ecpint)
C
      return
      end

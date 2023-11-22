      Subroutine Oprots_driver(Work,Maxcor,Iuhf)

      Implicit Double Precision(A-H,O-Z)
      Character*2 Type
      Character*1 UorL
      Logical One_pair 
      Dimension Work(Maxcor)
      Dimension Ioffbfn(8),Nbfirr(8)

#include "syminf.com"
#include "sym.com"
#include "machsp.com"
#include "info.com"

      Data Ione,Inull /1,0/

C This is written to debug the orbital rotation codes. The algorithim
C takes  pair of orbitals (OO,VV,OV and VO) and apply a unitary rotation (U)
C or mixing (L) of oribitals with  +eps and -eps 
C p' = cos (eps) p + Sin (eps) q
C q' = Sin(eps) p  - Cos (eps) q

      Call read_rotinfo(Type,UorL,Eps,Ip1,Ip2,Ip3,Ip4)

      If (Isym .Gt. Nirrep) Then 
         write(6,"(a)") " The number of irreps", i2, "exceed", Ieight 
         Call Errex
      Endif 
      If (Type .EQ. "OO") Then
      If (Abs(Ip1) .Gt. Nocco(1) .OR. Abs(Ip2) .Gt. Nocco(1)) Then 
      write(6,"(a,i4,a,I4,a,i4)") " The selected occupied orbitals", 
     +     Abs(Ip1)," or",Abs(Ip2)," exceeded the number of occupied:", 
     +     Nocco(1)
      Call Errex
      Endif 
      If (Ip1 .Gt. Inull .or. Ip2 .Gt. Inull .or.
     +    Ip3 .Gt. Inull .or. Ip4 .Gt. Inull) Then
         Write(6,"(2a)") " The occupied indices must be negative since",
     +                   " they are counted down from the Fermi level."
      Call Errex
      Endif 
      If (Abs(Ip3) .Gt. Nocco(1) .OR. Abs(Ip4) .Gt. Nocco(1)) Then 
      write(6,"(a,i4,a,I4,a,i4)") " The selected occupied orbitals", 
     +     Abs(Ip3)," or",Abs(Ip4)," exceeded the number of occupied:", 
     +     Nocco(1)
      Call Errex
      Endif 
      Endif 

      If (Type .EQ. "VV") Then
      If (Ip1 .Gt. Nvrto(1) .OR. Ip2 .Gt. Nvrto(1)) Then 
      write(6,"(a,i4,a,I4,a,i4)") " The selected virtual orbital", Ip1, 
     +      " or",Ip2," exceeded the number of virtuals:", Nvrto(1)
         Call Errex
      Endif 
      If (Ip2 .Gt. Nvrto(1) .OR. Ip3 .Gt. Nvrto(1)) Then 
      write(6,"(a,i4,a,I4,a,i4)") " The selected virtual orbital", Ip2, 
     +      " or",Ip3," exceeded the number of virtuals:", Nvrto(1)
         Call Errex
      Endif 
      If (Ip1 .Lt. Inull .or. Ip2 .Lt. Inull .or.
     +    Ip3 .Lt. Inull .or. Ip4 .Lt. Inull) Then
         Write(6,"(2a)") " The virtual indices must be non negative",
     +           " since they are counted up from the Fermi level."
      Call Errex
      Endif 
      Endif 

      If (Type .EQ. "OV") Then
      If (Abs(Ip1) .Gt. Nocco(1) .OR. Ip2 .Gt. Nvrto(1)) Then 
      write(6,"(2a,i4,a,I4,2a,2i4)") " The selected occupied or virtual",
     +      " orbital", Abs(Ip1)," or",Ip2," exceeded the number of"
      write(6,"(a,2(1x,i4))") " occupied or virtuals:", Nocco(1),
     +                          Nvrto(1)
      Call Errex
      Endif 
      If (Abs(Ip3) .Gt. Nocco(1) .OR. Ip4 .Gt. Nvrto(1)) Then 
      write(6,"(2a,i4,a,I4,2a)") " The selected occupied or virtual",
     +      " orbital", Abs(Ip2)," or",Ip3," exceeded the number of"
      write(6,"(a,2(1x,i4))") " occupied or virtuals:", Nocco(1),
     +                          Nvrto(1)
      Call Errex
      Endif 
      If (Ip1 .Gt. Inull .or. Ip2 .Lt. Inull .or.
     +    Ip3 .Gt. Inull .or. Ip4 .Lt. Inull) Then
         Write(6,"(2a)") " The occupied indices must be negative",
     +                   " and virtual indice must be postive."
         write(6,"(2a)") " since they are counted down or up from",
     +                   " the Fermi level."
      Call Errex
      Endif
      Endif 

      If (Type .EQ. "VO") Then
      If (Ip1 .Gt. Nvrto(1) .OR. Abs(Ip2) .Gt. Nocco(1)) Then 
      write(6,"(2a,i4,a,I4,2a,2i4)") " The selected virtual or occupied"
     +      " orbitals", Ip1," or",Abs(Ip2)," exceeded the number of",
     +      " virtuals or occupied:", Nvrto(1),Nocco(1)
         Call Errex
      Endif 
      If (Ip3 .Gt. Nvrto(1) .OR. Abs(Ip4) .Gt. Nocco(1)) Then 
      write(6,"(2a,i4,a,I4,2a,2i4)") " The selected virtual or occupied"
     +      " orbitals", Ip3," or",Abs(Ip4)," exceeded the number of",
     +      " virtuals or occupied:", Nvrto(1),Nocco(1)
         Call Errex
      Endif 
      If (Ip1 .Lt. Inull .or. Ip2 .Gt. Inull .or.
     +    Ip3 .Lt. Inull .or. Ip4 .Gt. Inull) Then
         Write(6,"(2a)") " The occupied indices must be negative",
     +                   " and virtual indice must be postive."
         write(6,"(2a)") " since they are counted down or up from",
     +                   " Fermi level."
      Call Errex
      Endif
      Endif 

C Type is the nature of the orbitals that are rotated (i.e. "OO", "VV",
C "VO" and "OV") and Ip1 and Ip2 are the MO orbital indices. These indices
C simply refere to the energy ordering of the MOS. 

      Nbasis_a   = Nocco(1) + Nvrto(1)
      Nbasis_b   = Nocco(2) + Nvrto(2)
      Nbasis     = Max(Nbasis_a,Nbasis_b)
      Call Getrec(20,"JOBARC","NBASTOT",Ione,Naobfns)
      
      Icmo_a   = Ione
      Icmo_b   = Icmo_a + Nbasis*Naobfns
      Imap_a   = Icmo_b + Nbasis*Naobfns
      Imap_b   = Imap_a + Nbasis 
      Isym     = Imap_b + Nbasis 
      Iend     = Isym   + Nbasis 

      Call Getrec(20,"JOBARC","SCFEVCA0",Nbasis*Naobfns,Work(Icmo_a))
      If (Iuhf .Eq. 1) Then
         Call Getrec(20,"JOBARC","SCFEVCB0",Nbasis*Naobfns,Work(Icmo_b))
      Endif 
#ifdef _DEBUG_LVL0
      call output(Work(Icmo_a),1,Nbasis,1,Nbasis,Nbasis,Nbasis,1)
#endif
      Call Getrec(20,"JOBARC","ORBMAP_A",Nbasis,Work(Imap_a))
      If (Iuhf .Eq. 1) Then
         Call Getrec(20,"JOBARC","ORBMAP_B",Nbasis,Work(Imap_b))
      Endif 
      Call Getrec(20,"JOBARC","IRPORDMO",Nbasis,Work(Imap_a))
      Call Getrec(20,"JOBARC","ORBSYM  ",Nbasis,Work(Isym))
      Call Getrec(20,"JOBARC","NUMBASIR",Nirrep,Nbfirr)

      Ioffbfn(1) = 1

      Do Irrep =2, Nirrep
         Ioffbfn(Irrep)=Nbfirr(Irrep-1) + Ioffbfn(Irrep-1)
      Enddo

      Call Get_orbsym(Type,Ip1,Ip2,Ip3,Ip4,Work(Imap_a),Work(Isym),
     +                Nbasis,Ip1sym,Ip2sym,Ip3sym,Ip4sym,Index1,
     +                Index2,Index3,Index4,One_pair)

#ifdef _DEBUG_LVL0
      Write(6,"(a)") " The sym-energy ordering map"
      Write(6,"(10(1x,i4))") (Work(Imap_a+i-1),i=1,Nbasis)
      Write(6,"(a)") " The orbital symmetries"
      Write(6,"(10(1x,i4))") (Work(Isym+i-1),i=1,Nbasis)
      Write(6,"(a,2a)") " Type                : ", Type
      Write(6,"(a)") " Symmetry offsets"
      Write(6,"(8(1x,i4))") (Ioffbfn(i),i=1,Nirrep)
      Write(6,"(a,2(1x,i4))") " Pair indices  : ", Ip1, Ip2
      Write(6,"(a,2(1x,i4))") " Pair symmetry : ", Ip1sym, Ip2sym
      if (.not. one_pair) then
      Write(6,"(a,2(1x,i4))") " Pair indices  : ", Ip1, Ip2
      Write(6,"(a,2(1x,i4))") " Pair symmetry : ", Ip1sym, Ip2sym
      endif 
      Write(6,"(a,F10.5)")    " Increment(eps): ", eps
#endif
      Write(6,"(a,2(1x,I3))")" Rotated orbital indices: ",Index1,Index2
      If (.not. one_pair) Then
      Write(6,"(a,2(1x,I3))")" Rotated orbital indices: ",Index3,Index4
      Endif 

      If (UorL .EQ. "U") Then
      Call Rotate_pairs(Work(Icmo_a),Eps,Index1,Index2,Ip1sym,Ip2sym,
     +                   Nbasis,Type)
      If (.Not. One_pair) Then
      Call Rotate_pairs(Work(Icmo_a),Eps,Index3,Index4,Ip3sym,Ip4sym,
     +                   Nbasis,Type)
      Endif 
      Elseif (UorL .EQ. "L") Then
 
      Call Mix_pairs(Work(Icmo_a),Eps,Index1,Index2,Ip1sym,Ip2sym,
     +                   Nbasis,Type)
      If (.Not. One_pair) Then
      Call Mix_pairs(Work(Icmo_a),Eps,Index3,Index4,Ip3sym,Ip4sym,
     +                   Nbasis,Type)
      Endif 
      Endif 

C Write the rotated pairs back into the JOBARC file.

      Call Putrec(20,"JOBARC","SCFEVCA0",Nbasis*Naobfns,Work(Icmo_a))
      If (Iuhf .Eq. 1) Then
         Call Putrec(20,"JOBARC","SCFEVCB0",Nbasis*Naobfns,Work(Icmo_b))
      Endif 
C  
      Return
      End




 

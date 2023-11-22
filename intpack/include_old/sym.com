
#ifndef _SYM_COM_
#define _SYM_COM_

c This common block contains information related to symmetry and lists.  A
c list is any four index quantity.  Since these quantities are so large,
c they are stored on disk using symmetry information to greatly reduce
c their size.  The four index quantities are treated as two-dimensional
c arrays, each dimension being a two index quantity.
c
c The lists are not stored in a single file.  Instead, they are split
c into several files (currently 5), each of which holds a certain number
c of lists (currently 100).  Each list has a set number of sublists
c (currently 10).  Because of the usefulness of lists, (currently) 2
c additional scratch list files are used.  These may be used within a
c module and will be automatically cleaned up at the end.
c
c A typical use of the sublist is to store each of the symmetry irreps in
c a separate sublist.  For Abelian symmetry, this if fine.  Some of the
c lists have subtypes not related to irreps, so this is certainly not the
c only way in which subtypes may be (or are) used.
c
c The actual physical description of the lists (data used ONLY in the
c getlst and putlst routines) are stored in the common block /lists/.
c This should NEVER appear in a standard routine.  The variables are
c described here however.
c
c maxirrep  : The highest number of symmetry irreps allowed.  For Abelian
c             symmetry, this is 8.
c num2comb  : This is the number of two-index combinations used in Aces3
c             lists.  These combinations include:
c                [ 1]  A<B    (alpha)    [12]  a,I  (ba)       [23]  a,B  (ba)
c                [ 2]  a<b    (beta)     [13]  A,b  (ab)       [24]  i,J  (ba)
c                [ 3]  I<J    (alpha)    [14]  I,j  (ab)       [25]  i,A  (ba)
c                [ 4]  i<j    (beta)     [15]  A,b  (ab)
c                [ 5]  A<=B   (alpha)    [16]  I,A  (alpha)
c                [ 6]  a<=b   (beta)     [17]  i,a  (beta)
c                [ 7]  I<=J   (alpha)    [18]  I,a  (ab)
c                [ 8]  i<=j   (beta)     [19]  A,B  (alpha)
c                [ 9]  A,I    (alpha)    [20]  a,b  (beta)
c                [10]  a,i    (beta)     [21]  I,J  (alpha)
c                [11]  A,i    (ab)       [22]  i,j  (beta)
c             where uppercase refers to alpha spin, lowercase to beta spin,
c             a and b are virtual orbitals, i and j are occupied orbitals.
c             Note that 13 and 15 are identical due to an early mistake in
c             Aces2.  Not all possible combinations are given, but these
c             are the only ones used.  Only the first 22 are used in the
c             lists.
c max2comb  : The number of 2 index combinations included in the list
c             above.
c
c numlistfile : The total number of list files (currently 5).
c scrlistfile : The number of scratch list files (currently 2).
c totlistfile : The total of numlistfile and scrlistfile.
c listperfile : The number of lists in each file (currently 100).
c numlist     : The maximum number of lists supported in the standard files.
c scrlist     : The maximum number of lists supported in the scratch files.
c totlist     : The maximum number of lists supported in the all files.
c numsublis   : The maximum number of sublists supported.
c
c nirrep    : The order of the computational point group (i.e. number of
c             symmetry irreps in the point group).
c dirprd    : The direct product table for this point group.  Given two
c             irreps (i and j), the irrep of "i x j" is given as dirprd(i,j).
c irpdpd    : The total population of all two index combinations by irrep.
c             irpdpd(irrep,comb) refers to the population of a two index
c             quantity of type comb (see table above) in the irrep.  Since
c             it is inconveniant to always be looking up the two index type
c             in the table, isytyp is defined.
c isytyp    : Each list has a row and column two index combination.  Each
c             of these combinations is one of the combinations in the list
c             above.  isytyp(1,ilist) returns the type of combination of
c             the row (left-hand) two index combination for ilist.
c             isytyp(2,ilist) returns the type for the column (right-hand)
c             combination of ilist.  Each type is the integer between 1
c             and num2comb.
c pop       : Number of occupied orbitals of each spin within each irrep.
c vrt       : Number of virtual orbitals of each spin within each irrep.
c nt        : Length of both alpha-alpha and beta-beta virtual-occupied
c             vectors which are overall totally symmetric.
c nfmi      : Length of both alpha-alpha and beta-beta occupied-occupied
c             vectors which are overall totally symmetric.
c nfea      : Length of both alpha-alpha and beta-beta virtual-virtual
c             vectors which are overall totally symmetric.
c isymoff   : This contains offsets in a 2-index vector of a given type
c             (all 25 types given above are included) which contains all
c             symmetry blocks.  The offset isymoff(irr1,irr2,type) looks
c             in a 2-index vector of the given type.  Inside that vector
c             are all combinations of symmetry irreps.  This returns the
c             offset of (irr1,irr2) in that list.  If type is between 1
c             and 22, this can be useful when working with lists.  As
c             stored here, the offset starts at 1 (probably not as useful
c             as if it started from 0) so you must subtract 1 from the
c             answer.
c occup(i,2): Orbital occupancy for each irrep for alpha/beta spin-orbitals.
c totocc(2) : The total orbital occupancy for alpha/beta spin-orbitals.
c totocca/b : Same information as totocc (we need scalars since these will
c           : be used as array bounds
c
c Many matrices are block diagonal, one block per irrep.  These blocks
c are also symmetric sometimes and can be stored in upper triangular form.
c The following arrays give information useful in storing these arrays.
c
c numbasir  : The number of SOs in each irrep.
c irrorboff : The offset of the first orbital in each irrep.  For example,
c             if irrep 1,2,3,4 contain 13,4,7,2 orbitals respectively, then
c             irrorboff contains 1,14,18,25.  The nirrep+1 element contains
c             the last orbital + 1.
c irrtrioff : When stored in triangular form as one long vector, this gives
c             the offset for the first element in each irrep.
c irrtrilen : This gives the length of each block.
c irrtritot : The sum of the elements in irrtrilen.
c irrsqroff : Similar to irrtrioff but for storing in square form.
c irrsqrlen :
c irrsqrtot :
c maxirrtri : The maximum length in irrtrilen.
c maxirrsqr : The maximum length in irrsqrlen.
c
c Some related variables appear in lists.com.

      integer maxirrep,num2comb,max2comb,numlist,numsublis,
     &    numlistfile,scrlistfile,totlistfile,scrlist,totlist,
     &    listperfile
      parameter (maxirrep=8)
      parameter (num2comb=22)
      parameter (max2comb=25)
      parameter (numlistfile=5)
      parameter (scrlistfile=2)
      parameter (totlistfile=numlistfile+scrlistfile)
      parameter (listperfile=100)
      parameter (numlist=numlistfile*listperfile)
      parameter (scrlist=scrlistfile*listperfile)
      parameter (totlist=numlist+scrlist)
      parameter (numsublis=10)

      integer nirrep,dirprd(maxirrep,maxirrep),
     &    irpdpd(maxirrep,num2comb),isytyp(2,totlist),
     &    pop(maxirrep,2),vrt(maxirrep,2),nt(2),nfmi(2),nfea(2),
     &    isymoff(maxirrep,maxirrep,max2comb),numbasir(maxirrep),
     &    occup(maxirrep,2),totocc(2),totocca,totoccb,
     &    irrorboff(maxirrep+1),maxirrtri,maxirrsqr,
     &    irrtrioff(maxirrep),irrtrilen(maxirrep),irrtritot,
     &    irrsqroff(maxirrep),irrsqrlen(maxirrep),irrsqrtot

      common /sym/ nirrep,occup,dirprd,irpdpd,isytyp,isymoff,
     &    numbasir,totocc,totocca,totoccb,
     &    pop,vrt,nt,nfmi,nfea,
     &    irrorboff,maxirrtri,maxirrsqr,irrtrioff,irrtrilen,
     &    irrsqroff,irrsqrlen,irrtritot,irrsqrtot

#endif /* _SYM_COM_ */


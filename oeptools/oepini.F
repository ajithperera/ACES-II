      subroutine oepini (
     &     luzmt, ncache, maxlen, maxpair,
     &     cache, list1, list2, str1,loslater)
c
c This routine parse the OEP options supplied in ZMAT via a separate ASV
c list and initializes the OEP configuration common block.
c 
      implicit none
c
c     Arguments
c     
      integer
     &     luzmt, ncache, maxlen, maxpair
c     
      character*(50)
     &     cache, list1 (maxpair), list2 (maxpair), str1
c
c     The OEP configuration common block
       logical loslater
#include <machsp.com>      
#include "oep.com"      
c
c     Local variables
c     
      logical report
      integer npair, asv_strlen
      character*1 asv_null
c
c     Messages to the user
c     
 5910 format (
     &     t3,'@OEPINI-I: Updating OEP1ON to ',i2,' (default: ',i2,').'/
     &     t3,'The first-order OEP method will be turned on when',/,
     &     t3,'the density difference is less than 10**(',i3,') a.u.'/)
 5920 format (
     &     t3,'@OEPINI-I: Updating OEP2ON to ',i2,' (default: ',i2,').'/
     &     t3,'The second-order OEP method will be turned on when',/,
     &     t3,'the density difference is less than 10**(',i3,') a.u.'/)
 5921 format (
     &     t3,'@OEPINI-I: Updating R12ON to ',i2,' (default: ',i2,').'/
     &     t3,'The second-order R12 method will be turned on when',/,
     &     t3,'the density difference is less than 10**(',i3,') a.u.'/)     
 5930 format (
     &     t3,'@OEPINI-I: Updating OEP2OFF to ',i2,' (default: ',i2,
     &     ').',/,t3,'The second-order OEP will not be updated if',/,
     &     t3,'the density difference is less than 10**(',i3,') a.u.'/)
 5935 format (
     &     t3,'@OEPINI-I: Updating OEP2ALWAYS to ',i2,' (default: ',i2,
     &     ').',/,t3,'The second-order OEP will be updated on every',/,
     &     t3,'iteration if the density difference is less than 10**(',
     &     i3,') a.u.'/)
 5940 format (
     &     t3,'@OEPINI-I: Updating OEP2UPDATE to ',i2,' (default: ',i2,
     &     t3,').',/,t3,'The second-order OEP will be updated every',/,
     &     t3,i2,' iterations.',/)
 5950 format (
     &     t3,'@OEPINI-I: Updating HYBRID to ',i3,'. The hybrid ',/,
     &     t3,'exchange operator will be constructed using ',i3,'% of',/
     &     t3,'the nonlocal and ',i3,'% of the local exact exchange',/,
     &     t3,'operators.',/)
 5955 format (
     &     t3,'@OEPINI-I: Update DFF to FALSE. The contribution of',/,
     &     t3,'the DFF diagram to the second-order potential will',/,
     &     t3,'be neglected.',/)
 5956 format (
     &     t3,'@OEPINI-I: Update DFW to FALSE. The contribution of',/,
     &     t3,'the DFW and DWF diagrams to the second-order potential',/,
     &     t3,'will be neglected.',/)
 5957 format (
     &     t3,'@OEPINI-I: Update DWW to FALSE. The contribution of',/,
     &     t3,'the DWW diagram to the second-order potential will',/,
     &     t3,'be neglected.',/)
 5958 format (
     &     t3,'@OEPINI-I: Update FDF to FALSE. The contribution of',/,
     &     t3,'the FDF diagram to the second-order potential will',/,
     &     t3,'be neglected.',/)
 5959 format (
     &     t3,'@OEPINI-I: Update WDW to FALSE. The contribution of',/,
     &     t3,'the WDW diagram to the second-order potential will',/,
     &     t3,'be neglected.',/)
 5960 format (
     &     t3,'@OEPINI-I: ',a,' will be used',/,
     &     t3,'to construct the reference Hamiltonian.',/)
 5961 format (
     &     t3,'@OEPINI-F: ',a,': unknown reference Hamiltonian.',/)
 5962 format (
     &     t3,'@OEPINI-I: Updating SVDTHR to ',i2,' (default: ',i2,').'/
     &     t3,'Eigenvalues smaller than 10**(',i3') will be discarded',/
     &     t3,'while inverting the response function.',/)
 5964 format (
     &     t3,'@OEPINI-I: ',a,' potential will be used to correct',/
     &     t3,'the asymptotic behavior of the exchange OEP.',/)
 5965 format (
     &     t3,'@OEPINI-I: Update OEPR12 to TRUE. The contribution of',/,
     &     t3,'MBPT(2)-F12 will be included',/)
 5970 format (
     &     t3,'@OEPINI-I: ',/,
     &     t3,70('-'),/,
     &     t3,13x,'OEP ASV',4x,'TYPE',1x,13x,'CURRENT',1x,13x,'DEFAULT'/
     &     t3,70('-'))
 5973 format (
     &     t3,a20,1x,a7,2(1x,a20))
 5976 format (
     &     t3,a20,1x,a7,2(1x,i20))
 5979 format (
     &     t3,70('-'),/)
c
c     Update the OEP control variables: <> Rewind the file. <> Parse the
c     OEP ASV list. <> Update the variables. <> Reports all the OEP
c     control ASV if requested.
c     
      asv_null=achar(0)
      call asv_parse_list (
     &     luzmt, ncache, maxlen, maxpair,
     &     cache, 'OEP'//asv_null, npair, list1, list2)
      call asv_update_boolean (
     &     maxlen, npair, list1, list2,
     &     'OEPR12'//asv_null, str1, .false., BOEPR12)
      if (OEP2.and.BOEPR12) write(6,5965)
      call asv_update_integer (
     &     maxlen, npair, list1, list2,
     &     'OEP1ON'//asv_null, str1, -20, 20, OEP1ONTHRDEF, OEP1ONTHR)
      if (OEP1ONTHR .ne. OEP1ONTHRDEF)
     &     write (6, 5910) OEP1ONTHR, OEP1ONTHRDEF, -OEP1ONTHR
      call asv_update_integer (
     &     maxlen, npair, list1, list2,
     &     'OEP2ON'//asv_null, str1,  -20, 20, OEP2ONTHRDEF,OEP2ONTHR)
      if (OEP2 .and. OEP2ONTHR .ne. OEP2ONTHRDEF)
     &     write (6, 5920) OEP2ONTHR, OEP2ONTHRDEF, -OEP2ONTHR
C
      call asv_update_integer (
     &     maxlen, npair, list1, list2,
     &     'R12ON'//asv_null, str1,  -20, 20, R12ONTHRDEF,R12ONTHR)      
C
      if (BOEPR12 .and. OEP2) then
        if (R12ONTHRDEF.ne.R12ONTHR) then
          write (6, 5921) R12ONTHR, R12ONTHRDEF, -R12ONTHR 
        else if (R12ONTHRDEF.ne.OEP2ONTHR) then
          R12ONTHR=OEP2ONTHR
          write (6, 5921) R12ONTHR, R12ONTHRDEF, -R12ONTHR
        endif
      endif
c     
      call asv_update_integer (
     &     maxlen, npair, list1, list2,
     &     'OEP2OFF'//asv_null, str1, -20, 20, OEP2OFFTHRDEF,OEP2OFFTHR)
      if (OEP2 .and. OEP2OFFTHR .ne. OEP2OFFTHRDEF)
     &     write (6, 5930) OEP2OFFTHR, OEP2OFFTHRDEF, -OEP2OFFTHR
c     
      call asv_update_integer (
     &     maxlen, npair, list1, list2,
     &     'OEP2ALWAYS'//asv_null, str1, -20, 20,OEP2ALWAYSTHRDEF,
     &      OEP2ALWAYSTHR)
      if (OEP2 .and. OEP2ALWAYSTHR .ne. OEP2ALWAYSTHRDEF)
     &     write (6, 5930) OEP2ALWAYSTHR,OEP2ALWAYSTHRDEF,-OEP2ALWAYSTHR
c     
      call asv_update_integer (
     &     maxlen, npair, list1, list2,
     &     'OEP2UPDATE'//asv_null, str1, 0, 999, OEP2UPDATEDEF, 
     &     OEP2UPDATE)
      if (OEP2 .and. OEP2UPDATE .ne. OEP2UPDATEDEF)
     &     write (6, 5940) OEP2UPDATE, OEP2UPDATEDEF, OEP2UPDATE
c     
      call asv_update_integer (
     &     maxlen, npair, list1, list2,
     &     'HYBRID'//asv_null, str1, 0, 999, OEP1HYBRIDDEF,
     &     OEP1HYBRIDINT)
      if (OEP1HYBRIDINT .ne. OEP1HYBRIDDEF)
     &     write (6,5950) OEP1HYBRIDINT,OEP1HYBRIDINT,100-OEP1HYBRIDINT
      OEP1HYBRID = dble (OEP1HYBRIDINT) / 100d0
c     
      call asv_update_boolean (
     &     maxlen, npair, list1, list2,
     &     'DFF'//asv_null, str1, .true., OEPDFF)
      if (OEP2.and..not.OEPDFF) write (6, 5955)
c     
      call asv_update_boolean (
     &     maxlen, npair, list1, list2,
     &     'DFW'//asv_null, str1, .true., OEPDFW)
      if (OEP2.and..not.OEPDFW) write (6, 5956)
c     
      call asv_update_boolean (
     &     maxlen, npair, list1, list2,
     &     'DWW'//asv_null, str1, .true., OEPDWW)
      if (OEP2.and..not.OEPDWW) write (6, 5957)
c     
      call asv_update_boolean (
     &     maxlen, npair, list1, list2,
     &     'FDF'//asv_null, str1, .true., OEPFDF)
      if (OEP2.and..not.OEPFDF) write (6, 5958)
c     
      call asv_update_boolean (
     &     maxlen, npair, list1, list2,
     &     'WDW'//asv_null, str1, .true., OEPWDW)
      if (OEP2.and..not.OEPWDW) write (6, 5959)
c     
      call asv_update_string (
     &     maxlen, npair, list1, list2,
     &     'REFERENCE'//asv_null, str1, 
     &     'OEP,DIAG,SEMI,KS,DF,SC'//asv_null, OEPH0DEF)
      if (str1 (1:4).eq. 'OEP'//asv_null .or. str1(1:3) .eq. 'KS') then
         OEPH0(1:5)='OEP '//asv_null
         write (6, 5960)
     &        'The OEP Hamiltonian eigenvalues'
      elseif (str1 (1:5) .eq. 'DIAG'//asv_null .or. 
     &        str1 (1:3) .eq. 'DF') then
         OEPH0(1:5)='DIAG'//asv_null
         write (6, 5960)
     &        'The diagonal Fock elements'
      elseif (str1 (1:5) .eq. 'SEMI'//asv_null .or. 
     &        str1 (1:3) .eq. 'SC') then
         OEPH0(1:5)='SEMI'//asv_null
         write (6, 5960)
     &        'The diagonal Fock elements and semicanonical orbitals'
      else
         write (6, 5961)
         call errex
      endif
      call putrec (
     &     20, 'JOBARC', 'OEPH0   ', (50-mod(50,iintln))/iintln, OEPH0)
c     
      call asv_update_boolean (
     &     maxlen, npair, list1, list2,
     &     'PRINTASV'//asv_null, str1, .false., report)
c     
      call asv_update_integer (
     &     maxlen, npair, list1, list2,
     &     'SVDTHR'//asv_null, str1, 0, 999, OEPSVDTHRDEF, OEPSVDTHR)
      if (OEPSVDTHR .ne. OEPSVDTHRDEF)
     &     write (6, 5962) OEPSVDTHR, OEPSVDTHRDEF, OEPSVDTHR
c     
      call asv_update_string (
     &     maxlen, npair, list1, list2,
     &     'ASYMPTOTIC'//asv_null, OEPAXP, 'NONE,FA,SLATER'//asv_null, 
     &     OEPAXPDEF)
      if     (OEPAXP (1:3) .eq. 'FA'//asv_null) then
         write (6, 5964) 'The scaled Coulomb'
      elseif (OEPAXP (1:7) .eq. 'SLATER'//asv_null) then
         write (6, 5964) 'The Slater'
              LOslater=.true.
      endif
      call putrec (
     &     20, 'JOBARC', 'OEPAXP  ', (50-mod(50,iintln))/iintln, OEPAXP)
c     
      if (report) then
         write (6, 5970)
         write (6, 5976) 'OEP1ON' , 'integer', OEP1ONTHR, OEP1ONTHRDEF
         write (6, 5976) 'HYBRID','integer',OEP1HYBRIDINT,OEP1HYBRIDDEF
         if (OEP2) then
            write (6, 5976)
     &           'OEP2ON',  'integer', OEP2ONTHR, OEP2ONTHRDEF
            write (6, 5976)
     &           'OEP2OFF', 'integer', OEP2OFFTHR, OEP2OFFTHRDEF
            write (6, 5976)
     &           'OEP2ALWAYS', 'integer',OEP2ALWAYSTHR,OEP2ALWAYSTHRDEF
            write (6, 5976)
     &           'OEP2UPDATE','integer',OEP2UPDATE,OEP2UPDATEDEF
         endif
         if (BOEPR12) 
     &     write (6, 5976) 'R12ON',  'integer', R12ONTHR, R12ONTHRDEF
         str1 (1:5) = 'FALSE'
         if (OEPDFF) str1 (1:5) = ' TRUE'
         write (6, 5973) 'DFF', 'boolean', str1 (1:5), 'TRUE'
         str1 (1:5) = 'FALSE'
         if (OEPDFW) str1 (1:5) = ' TRUE'
         write (6, 5973) 'DFW', 'boolean', str1 (1:5), 'TRUE'
         str1 (1:5) = 'FALSE'
         if (OEPDWW) str1 (1:5) = ' TRUE'
         write (6, 5973) 'DWW', 'boolean', str1 (1:5), 'TRUE'
         str1 (1:5) = 'FALSE'
         if (OEPFDF) str1 (1:5) = ' TRUE'
         write (6, 5973) 'FDF', 'boolean', str1 (1:5), 'TRUE'
         str1 (1:5) = 'FALSE'
         if (OEPWDW) str1 (1:5) = ' TRUE'
         write (6, 5973) 'WDW', 'boolean', str1 (1:5), 'TRUE'
         str1 (1:5) = 'FALSE'
c         if (OEPR12)  str1 (1:5) = ' TRUE'
c         write (6, 5973) 'OEPR12', 'boolean', str1 (1:5), 'FALSE'
c         str1 (1:5) = 'FALSE'
         write (6, 5973)
     &        'OEPH0', 'string',
     &        OEPH0 (1:asv_strlen (maxlen, OEPH0)-1),
     &        OEPH0DEF (1:asv_strlen (maxlen, OEPH0DEF)-1)
         write (6, 5976) 'SVDTHR' , 'integer', OEPSVDTHR, OEPSVDTHRDEF
         write (6, 5973)
     &        'ASYMPTOTIC', 'string',
     &        OEPAXP (1:asv_strlen (maxlen, OEPAXP)-1),
     &        OEPAXPDEF (1:asv_strlen (maxlen, OEPAXPDEF)-1)
         write (6, 5973) 'PRINTASV', 'boolean', ' TRUE', 'FALSE'
         write (6, 5979)
      endif
c
      return
      end

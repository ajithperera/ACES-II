C
      subroutine dkh__contract(zcore,maxcor,coreham,dkh_ints,coeffs,
     +                         nshells,nprims,ncnfns,npfps,npcoef,
     +                         npcoef_uc,ncfps,Ivangmom,Spherical)

      Implicit Double Precision (a-h, o-z)

#include "ecplims.par"

      Dimension zcore(maxcor),coreham(nprims*nprims),coeffs(npcoef)
      Dimension dkh_ints(ncnfns*ncnfns)
      Integer npfps(nshells),ncfps(nshells)
      Integer ivangmom(Nshells)
      Integer coef_start,zcor_start,prim_off(max_shells)
      Integer cont_off(max_shells),Ccoef_off(max_shells)
      Integer Tot_prm_fns,Tot_cnt_fns,Tot_cnt_ish,Tot_cnt_jsh
      Integer Icnt_block
      
      logical spherical 

      prim_off(1)  = 1
      cont_off(1)  = 1
      ccoef_off(1) = 1

      Do ishell = 2, Nshells
         ncont_shell = ncfps(ishell-1)
         nprim_shell = npfps(ishell-1)

         if (spherical) then
            ndegen = 2*ivangmom(ishell-1) + 1
         else
            ndegen = (ivangmom(ishell-1)+1) * (ivangmom(ishell-1)+2)/2
         endif

         prim_off(ishell) = prim_off(ishell-1) + nprim_shell * ndegen
         cont_off(ishell) = cont_off(ishell-1) + ncont_shell * ndegen   
         Ccoef_off(ishell)= Ccoef_off(ishell-1) + 
     +                      ncont_shell * nprim_shell
      Enddo

#ifdef _DEBUG_LVLM
      write(6,*) "The primitive DKH integrals"
      call output(coreham,1,nprims,1,nprims,nprims,nprims,1)
      write(6,*) "The contraction coefs."
      write(6,"(6(1x,f10.6))") (coeffs(i),i=1,npcoef)
      write(6,*) "The number of primitive functions"
      write(6,"(6(1x,i4))") (npfps(i),i=1,nshells)
      write(6,*) "The number of contracted functions"
      write(6,"(6(1x,i4))") (ncfps(i),i=1,nshells)
      write(6,"(a,l)") "Spherical?",Spherical
      Write(6,"(a,1x,6i4)") "prim. off-sets:",(prim_off(i),i=1,Nshells)
      Write(6,"(a,1x,6i4)") "contr. off-sets:",(cont_off(i),i=1,Nshells)
      Write(6,"(a,1x,6i4)") "Ccoef. off-sets:",(Ccoef_off(i),
     +                       i=1,Nshells)
      ioff = 0
      ncount = 0
#endif 


      Do Jshell = 1, Nshells
     
         Do Ishell = 1, Nshells 

            Iprim_4shell = Npfps(Ishell)
            Jprim_4shell = Npfps(Jshell)
            Icntr_4shell = Ncfps(Ishell)
            Jcntr_4shell = Ncfps(Jshell)

            If (Spherical) Then
               Idegen = 2*ivangmom(ishell) + 1
               Jdegen = 2*ivangmom(jshell) + 1
            Else 
               Idegen = (ivangmom(ishell)+1) * (ivangmom(ishell)+2)/2
               Jdegen = (ivangmom(jshell)+1) * (ivangmom(jshell)+2)/2
            Endif
           
            Iprimfns = Iprim_4shell * Idegen
            Jprimfns = Jprim_4shell * Jdegen
            
            Icntrfns = Icntr_4shell * Idegen
            Jcntrfns = Jcntr_4shell * Jdegen

            Tot_prm_fns = Iprimfns * Jprimfns 
            Tot_cnt_fns = Icntrfns * Jcntrfns 
            Tot_cnt_ish = Iprim_4shell * Icntr_4shell 
            Tot_cnt_jsh = Jprim_4shell * Jcntr_4shell 
            Tot_bfns    = Ncnfns * Ncnfns
C

            Istart      = 1
            Ipnt_block  = Istart 
            Icn_block   = Ipnt_block + Tot_prm_fns 
            Jcn_block   = Icn_block  + Tot_cnt_ish
            Ict_block   = Jcn_block  + Tot_cnt_jsh
            Inext       = Ict_block  + Tot_cnt_fns 

#ifdef _DEBUG_LVLM
      write(6,*)
      write(6,"(a,2(1x,i4))")"Ishell,Jshell",Ishell,Jshell
      write(6,"(a,2(1x,i4))")"Iprim_4shell,Iprim_4shell",
     +                        Iprim_4shell,Jprim_4shell
      write(6,"(a,2(1x,i4))")"Icntr_4shell,Icntr_4shell",
     +                        Icntr_4shell,Jcntr_4shell
      write(6,"(a,2(1x,i4))")"Idegen,Jdegen",Idegen,Jdegen
      write(6,"(a,2(1x,i4))")"Iprimfns,Jprimfns",Iprimfns,Jprimfns
      write(6,"(a,6(1x,i4))")"Icntrfns,Jcntrfns",Icntrfns,Jcntrfns
      write(6,"(a,6(1x,i4))")"Icntrfns,Jcntrfns",Icntrfns,Jcntrfns
      write(6,"(a,6(1x,i4))")"Tot_prm_fns,Tot_cnt_fns",
     +                        Tot_prm_fns,Tot_cnt_fns
      Write(6,"(a,6(1x,i4))")"Tot_cnt_ish.Tot_cnt_jsh",Tot_cnt_ish,
     +                  Tot_cnt_jsh
#endif 

            Call get_ints_block(coreham,zcore(Ipnt_block),ishell,
     +                          jshell,Nprims,prim_off,cont_off,
     +                          Iprimfns,jprimfns,Nshells)

            Call get_ccoefs_block(Coeffs,zcore(Icn_block),ishell,
     +                            Npcoef,ccoef_off,
     +                            Tot_cnt_ish,Nshells)

            Call get_ccoefs_block(Coeffs,zcore(Jcn_block),jshell,
     +                            Npcoef,Ccoef_off,
     +                            Tot_cnt_jsh,NShells)

#ifdef _DEBUG_LVLM
      write(6,*) 
      write(6,"(a)") "The DKH integrals for shell-block"
      ncount = ncount + 1
      ioff = ioff + Tot_prm_fns*(ncount-1)
      write(6,"(6(1x,E12.5))") (zcore(Ipnt_block+i-1),i=1,Tot_prm_fns)
      write(6,*) 
      Write(6,"(a,2i4)") "The contraction coefs. of shell:",ishell
      write(6,"(6(1x,E12.5))") (zcore(Icn_block+i-1),i=1, 
     +                          Tot_cnt_ish)
      Write(6,"(a,2i4)") "The contraction coefs. of shell:",jshell
      write(6,"(6(1x,E12.5))") (zcore(Icn_block+i-1),i=1, 
     +                          Tot_cnt_jsh)
#endif 
            Call contract(zcore(Ipnt_block),zcore(Icn_block),
     +                    zcore(Jcn_block),Iprimfns,jprimfns,
     +                    Icntrfns,Jcntrfns,Iprim_4shell,
     +                    Jprim_4shell,Icntr_4shell,Jcntr_4shell,
     +                    Tot_cnt_fns,Tot_prm_fns,Idegen,Jdegen,
     +                    zcore(Ict_block))

            Call form_dkh_ints(zcore(Ict_block),Dkh_ints,
     +                         Ishell,jshell,Icntrfns,Jcntrfns,
     +                         Ncnfns,prim_off,cont_off,NShells)

         Enddo
      Enddo

#ifdef _DEBUG_LVLM
      Write(6,*)
      Write(6,"(a,a)") "@-dkh__contract: The contracted DKH",
     +                 " 1-el Hamiltonian"
      call output(dkh_ints,1,Ncnfns,1,Ncnfns,Ncnfns,Ncnfns,1)
CSS      call dkh__sq_print(Ncnfns,dkh_ints,"")
#endif 

      Return 
      End

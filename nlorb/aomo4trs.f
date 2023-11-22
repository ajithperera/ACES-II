      subroutine aomo4trs(aobas,coef,mobas,modim1,modim2,
     & modim3,modim4,nbas,nocc,nvir,oldmo1,oldmo2)

      integer modim1, modim2, modim3, modim4, nbas, nocc, nvir,
     & iii, jjj, aaa, bbb, mu, nu, lam, sig, ilev, jlev, alev, blev

      double precision aobas(nbas,nbas,nbas,nbas),
     & coef(nbas,nbas), mobas(modim1,modim2,modim3,modim4),
     & oldmo1(nbas,nbas,nbas,nbas), oldmo2(nbas,nbas,nbas,nbas)

      if (modim3.eq.nocc) then
         ilev = 0
      else if (modim3.eq.nvir) then
         ilev = nocc
      end if

      if (modim4.eq.nocc) then
         jlev = 0
      else if (modim4.eq.nvir) then
         jlev = nocc
      end if

      if (modim1.eq.nocc) then
         alev = 0
      else if (modim1.eq.nvir) then
         alev = nocc
      end if

      if (modim2.eq.nocc) then
         blev = 0
      else if (modim2.eq.nvir) then
         blev = nocc
      end if

      call dzero4arr(nbas,nbas,nbas,nbas,oldmo1)

      do bbb = 1, modim2
         do mu = 1, nbas
            do nu = 1, nbas
               do sig = 1, nbas
                  do lam = 1, nbas
       oldmo1(mu,nu,sig,bbb+blev) = oldmo1(mu,nu,sig,bbb+blev) +
     & coef(lam,bbb+blev)*aobas(mu,nu,sig,lam)
                  end do
               end do
            end do
         end do
      end do

      call dzero4arr(nbas,nbas,nbas,nbas,oldmo2)

      do jjj = 1, modim4
         do bbb = 1, modim2
            do mu = 1, nbas
               do nu = 1, nbas
                  do sig = 1, nbas
       oldmo2(mu,nu,jjj+jlev,bbb+blev) = oldmo2(mu,nu,jjj+jlev,
     & bbb+blev) + coef(sig,jjj+jlev)*oldmo1(mu,nu,sig,bbb+blev)
                  end do
               end do
            end do
         end do
      end do

      call dzero4arr(nbas,nbas,nbas,nbas,oldmo1)
      
      do aaa = 1, modim1
         do jjj = 1, modim4
            do bbb = 1, modim2
               do mu = 1, nbas
                  do nu = 1, nbas
       oldmo1(mu,aaa+alev,jjj+jlev,bbb+blev) = oldmo1(mu,
     & aaa+alev,jjj+jlev,bbb+blev) + coef(nu,aaa+alev)*
     & oldmo2(mu,nu,jjj+jlev,bbb+blev)
                  end do
               end do
            end do
         end do
      end do

      call dzero4arr(nbas,nbas,nbas,nbas,oldmo2)

      do iii = 1, modim3
         do aaa = 1, modim1
            do jjj = 1, modim4
               do bbb = 1, modim2
                  do mu = 1, nbas
       oldmo2(iii+ilev,aaa+alev,jjj+jlev,bbb+blev) = oldmo2(iii+ilev,
     & aaa+alev,jjj+jlev,bbb+blev) + coef(mu,iii+ilev)*
     & oldmo1(mu,aaa+alev,jjj+jlev,bbb+blev)
                  end do
               end do
            end do
         end do
      end do

      do iii = 1, modim3
         do aaa = 1, modim1
            do jjj = 1, modim4
               do bbb = 1, modim2
      mobas(aaa,bbb,iii,jjj) = oldmo2(iii+ilev,
     & aaa+alev,jjj+jlev,bbb+blev)
               end do
            end do
         end do
      end do

      return
      end


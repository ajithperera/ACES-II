C  This subroutine performs the ZMP procedure
C given a trial density (tmp_4) and all the standard HF SCF input
C data (onehao, trnfor, twoint, norbs) and gives the density for 
C a series of lambdas (alternatively we could make it so that lambda
C is included as an input)

      Subroutine getzmp(onehao,twoint,trnfor,norbs,
     &                  fock,dens,tmp_1,tmp_2,tmp_3,tmp_4)
      Implicit Double Precision (a-h, o-z)
      Integer numocc,numprot,indtot1
      Dimension onehao(norbs,norbs), 
     & twoint(norbs*norbs*norbs*norbs),trnfor(norbs,norbs),
     & fock(norbs,norbs),dens(norbs,norbs),tmp_1(norbs,norbs),
     & tmp_2(norbs,norbs),tmp_3(norbs,norbs),tmp_4(norbs,norbs)
      COMMON /MACHSP/ IINTLN,IFLTLN,IINTFP,IALONE,IBITWD
      ione=1
      Call zero(dens,norbs*norbs*IINTFP)  
      Call getrec(20,'JOBARC','NMPROTON',ione,numprot)
      Call getrec(20,'JOBARC','NUCREP ',ione,anucrep)
      numocc=numprot/2
      write(luout,*)'the number of orbitals =',norbs
      write(luout,*)'the number of occupied orbitals =',numocc
C    tmp_4=P_o
C    dens=P_lambda
C  Get the relaxed density from CC calculation
C alternatively this could be input as tmp_4
      Call getrec(20,'JOBARC','RELDENSA',norbs*norbs*iintfp,tmp_4)
C Pick the maximum iteration (maxit) and lambda
       maxit=200
       lambda=13
 
C########### Begin SCF Iteration ############
      Do 50 iteration=1,maxit
         Call zero(tmp_1,norbs*norbs)
         Call zero(tmp_2,norbs*norbs)
C**** Calculate F from H and two-e ints ****
       Do 20 imu=1,norbs
         Do 20 inu=imu,norbs
           Do 10 ilam=1,norbs
             Do 10 isig=1,norbs
             Call mysort(imu,inu,isig,ilam,indtot1)
      tmp_1(imu,inu)=
     &   twoint(indtot1)*(dens(ilam,isig)*(1.d0-1.d0/dble(2*numocc)+
     &   dble(lambda))-dble(lambda)*tmp_4(ilam,isig))+tmp_1(imu,inu)
 10   Continue
         tmp_1(inu,imu)=tmp_1(imu,inu)
         fock(imu,inu)=onehao(imu,inu)+tmp_1(imu,inu)
         fock(inu,imu)=fock(imu,inu)
 20    Continue
C***** F Calculated
C calculate F' from X^T(trnfor, 't'), F and X(trnfor, 'n')
      Call xgemm('t','n',norbs,norbs,norbs,1.D0,trnfor,norbs,
     &           fock,norbs,0.D0,tmp_2,norbs)
C tmp_2=X^T*F
      Call xgemm('n','n',norbs,norbs,norbs,1.D0,tmp_2,norbs,
     &           trnfor,norbs,0.D0,tmp_1,norbs)
C tmp_1=F'
C calculate C' and E from F'
      Call eig(tmp_1,tmp_2,norbs,norbs,1)
C calculate C= XC'
      Call xgemm('n','n',norbs,norbs,norbs,1.D0,trnfor,norbs,
     &           tmp_2,norbs,0.D0,tmp_1,norbs)
C tmp_1=C
C calculate P of occupied orbitals from C
      Call dgemm('n','t',norbs,norbs,numocc,2.D0,tmp_1,norbs,tmp_1,
     &           norbs,0.D0,dens,norbs)
      Do 111 imu=1,norbs
      Do 111 inu=1,norbs
        dens(imu,inu)=(dens(imu,inu)+tmp_3(imu,inu))/2.d0
 111    Continue

C ******** begin Check for convergence   ********
      bigdiff=0.d0
      Do 30 imu=1,norbs
         Do 30 inu=1,norbs
          diff=abs(dens(imu,inu)-tmp_3(imu,inu))
          tmp_3(imu,inu)=dens(imu,inu)
         If (diff.gt.bigdiff) Then
            bigdiff=diff
         End If

 30      Continue
C         write(luout,*)'bigdiff of iteration',iteration,' is',
C     &   bigdiff
      If (bigdiff.lt.(1.d-6)) Then
       write(luout,*)'SCF has converged'
       GoTo 60
      End If
      If (iteration.gt.maxit) Then
       write(luout,*)'SCF has not converged'
       GoTo 60
      End If
C ******* end Check for convergence  **********
C############ end SCF Iteration  #############

 50   Continue
 60   Continue
      
C  Now we want to calculate the coefficient matrix of the vxc

      Do 90 imu=1,norbs
        Do 90 inu=1,norbs
          tmp_2(imu,inu)=(-1.d0/dble(2*numocc)+dble(lambda))
     &    *dens(imu,inu)-dble(lambda)*tmp_4(imu,inu)
 90   Continue
       write(luout,*)'The Vxc matrix'
       call output(tmp_2,1,norbs,1,norbs,norbs,norbs,0)
      Return
      End 












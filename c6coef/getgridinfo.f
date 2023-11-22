      subroutine getgridinfo(numpoints,wghts,freqs,grdvl,omega0)
      implicit none

      integer numpoints

      double precision wghts(numpoints), freqs(numpoints),
     & grdvl(numpoints)

      double precision omega0

      if (numpoints.eq.4) then

      grdvl(1) = -0.86113D0 
      grdvl(2) = -0.33998D0
      grdvl(3) =  0.33998D0
      grdvl(4) =  0.86113D0 

      wghts(1) =  0.34785D0
      wghts(2) =  0.65214D0
      wghts(3) =  0.65214D0
      wghts(4) =  0.34785D0
      
      freqs(1) = omega0*(1+0.86113D0)/(1-0.86113D0)
      freqs(2) = omega0*(1+0.33998D0)/(1-0.33998D0)
      freqs(3) = omega0*(1-0.33998D0)/(1+0.33998D0)
      freqs(4) = omega0*(1-0.86113D0)/(1+0.86113D0)

      end if
      
      if (numpoints.eq.8) then
      
      grdvl(1) = -0.96028D0
      grdvl(2) = -0.79666D0
      grdvl(3) = -0.52553D0
      grdvl(4) = -0.18343D0
      grdvl(5) =  0.18343D0
      grdvl(6) =  0.52553D0
      grdvl(7) =  0.79666D0
      grdvl(8) =  0.96028D0
      
      wghts(1) =  0.10122D0
      wghts(2) =  0.22238D0
      wghts(3) =  0.31370D0
      wghts(4) =  0.36268D0
      wghts(5) =  0.36268D0
      wghts(6) =  0.31370D0
      wghts(7) =  0.22238D0
      wghts(8) =  0.10122D0
      
      freqs(1) = omega0*(1+0.96028D0)/(1-0.96028D0)
      freqs(2) = omega0*(1+0.79666D0)/(1-0.79666D0)
      freqs(3) = omega0*(1+0.52553D0)/(1-0.52553D0)
      freqs(4) = omega0*(1+0.18343D0)/(1-0.18343D0)
      freqs(5) = omega0*(1-0.18343D0)/(1+0.18343D0)
      freqs(6) = omega0*(1-0.52553D0)/(1+0.52553D0)
      freqs(7) = omega0*(1-0.79666D0)/(1+0.79666D0)
      freqs(8) = omega0*(1-0.96028D0)/(1+0.96028D0)
     
      end if
      
      if (numpoints.eq.12) then

      grdvl(1) = -0.98156D0
      grdvl(2) = -0.90411D0
      grdvl(3) = -0.76990D0
      grdvl(4) = -0.58731D0
      grdvl(5) = -0.36783D0
      grdvl(6) = -0.12523D0
      grdvl(7) =  0.12523D0
      grdvl(8) =  0.36783D0
      grdvl(9) =  0.58731D0
      grdvl(10)=  0.76990D0
      grdvl(11)=  0.90411D0
      grdvl(12)=  0.98156D0

      wghts(1) =  0.04717D0
      wghts(2) =  0.10693D0
      wghts(3) =  0.16007D0
      wghts(4) =  0.20316D0
      wghts(5) =  0.23349D0
      wghts(6) =  0.24914D0
      wghts(7) =  0.24914D0
      wghts(8) =  0.23349D0
      wghts(9) =  0.20316D0
      wghts(10)=  0.16007D0
      wghts(11)=  0.10693D0
      wghts(12)=  0.04717D0
      
      freqs(1) = omega0*(1+0.98156D0)/(1-0.98156D0)
      freqs(2) = omega0*(1+0.90411D0)/(1-0.90411D0)
      freqs(3) = omega0*(1+0.76990D0)/(1-0.76990D0)
      freqs(4) = omega0*(1+0.58731D0)/(1-0.58731D0)
      freqs(5) = omega0*(1+0.36783D0)/(1-0.36783D0)
      freqs(6) = omega0*(1+0.12523D0)/(1-0.12523D0)
      freqs(7) = omega0*(1-0.12523D0)/(1+0.12523D0)
      freqs(8) = omega0*(1-0.36783D0)/(1+0.36783D0)
      freqs(9) = omega0*(1-0.58731D0)/(1+0.58731D0)
      freqs(10)= omega0*(1-0.76990D0)/(1+0.76990D0)
      freqs(11)= omega0*(1-0.90411D0)/(1+0.90411D0)
      freqs(12)= omega0*(1-0.98156D0)/(1+0.98156D0)

      end if
      
      if (numpoints.eq.16) then
      
      grdvl(1) = -0.98940D0
      grdvl(2) = -0.94457D0
      grdvl(3) = -0.86563D0
      grdvl(4) = -0.75540D0
      grdvl(5) = -0.61787D0
      grdvl(6) = -0.45801D0
      grdvl(7) = -0.28160D0
      grdvl(8) = -0.09501D0
      grdvl(9) =  0.09501D0
      grdvl(10)=  0.28160D0
      grdvl(11)=  0.45801D0
      grdvl(12)=  0.61787D0
      grdvl(13)=  0.75540D0
      grdvl(14)=  0.86563D0
      grdvl(15)=  0.94457D0
      grdvl(16)=  0.98940D0
       
      wghts(1) =  0.02715D0
      wghts(2) =  0.06225D0
      wghts(3) =  0.09515D0
      wghts(4) =  0.12462D0
      wghts(5) =  0.14959D0
      wghts(6) =  0.16915D0
      wghts(7) =  0.18260D0
      wghts(8) =  0.18945D0
      wghts(9) =  0.18945D0
      wghts(10)=  0.18260D0
      wghts(11)=  0.16915D0
      wghts(12)=  0.14959D0
      wghts(13)=  0.12462D0
      wghts(14)=  0.09515D0
      wghts(15)=  0.06225D0
      wghts(16)=  0.02715D0
       
      freqs(1) = omega0*(1+0.98940D0)/(1-0.98940D0)
      freqs(2) = omega0*(1+0.94457D0)/(1-0.94457D0)
      freqs(3) = omega0*(1+0.86563D0)/(1-0.86563D0)
      freqs(4) = omega0*(1+0.75540D0)/(1-0.75540D0)
      freqs(5) = omega0*(1+0.61787D0)/(1-0.61787D0)
      freqs(6) = omega0*(1+0.45801D0)/(1-0.45801D0)
      freqs(7) = omega0*(1+0.28160D0)/(1-0.28160D0)
      freqs(8) = omega0*(1+0.09501D0)/(1-0.09501D0)
      freqs(9) = omega0*(1-0.09501D0)/(1+0.09501D0)
      freqs(10)= omega0*(1-0.28160D0)/(1+0.28160D0)
      freqs(11)= omega0*(1-0.45801D0)/(1+0.45801D0)
      freqs(12)= omega0*(1-0.61787D0)/(1+0.61787D0)
      freqs(13)= omega0*(1-0.75540D0)/(1+0.75540D0)
      freqs(14)= omega0*(1-0.86563D0)/(1+0.86563D0)
      freqs(15)= omega0*(1-0.94457D0)/(1+0.94457D0)
      freqs(16)= omega0*(1-0.98940D0)/(1+0.98940D0)

      end if

      return
      end


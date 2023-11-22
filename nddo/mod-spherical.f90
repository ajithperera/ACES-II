 module spherical
! this is the transformation matrix from cartesian gaussians to spherical harmonics
double precision,dimension(9,10)::trans_columns=reshape( source =(/ &
1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0, &
0.0d0,1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0, &
0.0d0,0.0d0,1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0, &
0.0d0,0.0d0,0.0d0,1.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0, &
0.0d0,0.0d0,0.0d0,0.0d0,0.866025403d0,-0.5d0,0.0d0,0.0d0,0.0d0, &
0.0d0,0.0d0,0.0d0,0.0d0,-0.866025403d0,-0.5d0,0.0d0,0.0d0,0.0d0, &
0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,1.0d0,0.0d0,0.0d0,0.0d0, &
0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,1.0d0,0.0d0,0.0d0, &
0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,1.0d0,0.0d0, &
0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,1.0d0 /), shape=(/9,10/) )


end module spherical

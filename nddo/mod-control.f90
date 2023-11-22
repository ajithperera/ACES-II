module control
! this globally accessed module holds all the available keywords
double precision::scftol,OPTTOL,mult,sys_charge
character(20),dimension(:),allocatable::keyword
logical::RESTRICTED,DEORTHO,GRADIENT,TRIAL,keep,XYZ,OPTIMIZE,parameters,CENTERS &
,newton,DFTHF,LAPACK,SPARKLES,densityin,densityout,FREQUENCY,SAVE_TREE
character(10)::core
character(6)::method
end module control

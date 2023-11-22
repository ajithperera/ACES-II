      SUBROUTINE PCCD_DODIIS(ERR,AUGERR,TMP,NPHYS,NDIM)                                                                      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                                                               
      DIMENSION ERR(NPHYS,NPHYS),AUGERR(NPHYS+1,NPHYS+2)                                                                
      DIMENSION TMP(*),DET(2)                                                                                           
      CALL ZERO(AUGERR,(NPHYS+1)*(NPHYS+2))                                                                             
      CALL BLKCPY(ERR,NPHYS,NPHYS,AUGERR,NPHYS+1,NPHYS+1,2,2)                                                           
      DO 10 I=2,NDIM+1                                                                                                  
       AUGERR(I,1)=-1.0D0                                                                                               
       AUGERR(1,I)=-1.0D0                                                                                               
10    CONTINUE                                                                                                          
      AUGERR(1,NDIM+2)=-1.0d0                                                                                           
      CALL ZERO(AUGERR(2,NDIM+2),NDIM)                                                                                  
      CALL MINV(AUGERR,NDIM+1,NPHYS+1,TMP,DET,0.0D0,1,0)                                                                
      CALL SCOPY(NDIM,AUGERR(2,NDIM+2),1,TMP,1)                                                                         
      RETURN                                                                                                            
      END                                                                                                               

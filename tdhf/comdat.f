C   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX    A  00030
C  XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX  A  00040
C   ******************************************************************  A  00050
C                                                                       A  00060
C     ===== FILE INFORMATION =====                                      A  00070
C                                                                       A  00080
C     FT03F001        P******** FIEL >>>>> E1                           A  00090
C     FT08F001        2-ELECTRON A.O. INTEGRAL FIEL                     A  00100
C     FT10F001        DIRECT ACCESS FILE                                A  00110
C                              FROM GAMESS PROGRAM                      A  00120
C      REC #  = IODA(1) :  ON OR OFF                                    A  00130
C               IODA(12) :  OVERLAP INTEGRAL                            A  00140
C               IODA(15) : EIGEN FUNCTIONS(ORBITAL COEFFICIENT)         A  00150
C               IODA(16) : DENSITY MATRIX                               A  00160
C               IODA(17) : EIGEN VALUES                                 A  00170
C                                                                       A  00180
C     ======  MARK ======                                               A  00190
C      NSIZ1 : Dimension of Orbital Space 
C      NSIZ2  : NSIZ1*NSIZ1         
C      NSIZ3 : ((NSIZ1+1)*NSIZ1)/2   
C              ((NBASIS+1)*NBASIS)/2
C                                                                       A  00230
C      INP : SWITCH FOR DIFFERENT INPUTS 
C                                                                       A  00250
CEND====================================================================A  00260
      BLOCK DATA COMDAT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ITITLE(18),LCOMP(21,4)       
      COMMON/INFT/ITITLE,NCOMP,LCOMP,MSIZVO
      COMMON/ITOPR/ IDCSHG,IOKE,IDCOR,IIDRI,ITHG
      COMMON/SWPPP/INP,IAMO,IFAMO,IINDO,IORTH                
      COMMON/IPRNT/ IOPEV,IOPDA,IOPU,IOPFE,IOPPR,IWRPA
      COMMON/ILINEA/ISALPH,IDALPH
      COMMON/INBETA/ISBETA,ISHG,IEOPE,IOR
C     COMMON/TDHFIN/FREQ,NTDHF,IFIL1,IPROP 
      COMMON/MFREQ/NFREQ
      COMMON/THRES/ TOLPER,DEGEN,EPSI
      COMMON/THRE1/ NITER,MAXIT
      COMMON /LRDUCE/ CONVI
      COMMON /LRDUC1/ KMAX
      COMMON/CONST/ EBETA,EGAMMA 
C     MSIZVO( = 100); the maximum size of the matrix in the non-iterative 
C     solution. This is overridden in Main and MSIZVO=NSIZVO, if LSIZVO
C     is not specified.
C     KMAX( = 20); maximum size of reduced equation
C     NCOMP( = 21 ); identical commponent
      DATA NCOMP/21/,KMAX/20/,MSIZVO/100/,IOPEV/0/,IOPDA/0/,
     X IOPU/0/,IOPFE/0/,IOPPR/0/,IWRPA/0/,NITER/20/,MAXIT/50/,
     X IDCSHG/0/,IOKE/0/,IDCOR/0/,IIDRI/0/,ITHG/0/,INP/0/,CONVI/1.D-10/
     X ,NFREQ/1/,TOLPER/5.D-5/,EPSI/0./,IAMO/0/,IFAMO/0/,IORTH/0/
     X ,IINDO/0/,DEGEN/1.D-8/,ISALPH/0/,IDALPH/0/,ISBETA/0/,ISHG/0/
     X ,IEOPE/0/,IOR/0/
      DATA LCOMP/1,2,3,1,2,1,2,1,2,2,3,2,3,2,3,1,3,1,3,1,3,
     X           1,2,3,1,2,2,1,2,1,2,3,3,2,3,2,1,3,3,1,3,1,
     X           1,2,3,2,1,2,1,1,2,3,2,3,2,2,3,3,1,3,1,1,3,
     X           1,2,3,2,1,1,2,2,1,3,2,2,3,3,2,3,1,1,3,3,1/
      END

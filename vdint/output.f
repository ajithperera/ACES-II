CLIC                                                                            
C                                                                               
C...   Copyright (c) 1992 - 2010 by the authors of CFOUR (see below).           
C...   All Rights Reserved.                                                     
C...                                                                            
C...         The source code in this file is part of                            
C...                                                                            
C...   CFOUR, Coupled Cluster techniques for Computational Chemistry,           
C...   a quantum-chemical program package by                                    
C...                                                                            
C...   J.F. Stanton, J. Gauss, M.E. Harding, P.G. Szalay                        
C...                                                                            
C...   with contributions from                                                  
C...                                                                            
C...   A.A. Auer, R.J. Bartlett, U. Benedikt, C. Berger, D.E. Bernholdt,        
C...   Y.J. Bomble, O. Christiansen, M. Heckert, O. Heun, C. Huber, T.-C. Jagau,
C...   D. Jonsson, J. Juselius, K. Klein, W.J. Lauderdale, D. A. Matthews,      
C...   T. Metzroth, D.P. O'Neill, D.R. Price, E. Prochnow, K. Ruud,             
C...   F. Schiffmann, S. Stopkowicz, A. Tajti, J. Vazquez, F. Wang,             
C...   J.D. Watts                                                               
C...   and the integral packages                                                
C...   MOLECULE (J. Almlof and P.R. Taylor),                                    
C...   PROPS (P.R. Taylor),                                                     
C...   ABACUS (T. Helgaker, H.J. Aa. Jensen, P. Jorgensen, and J. Olsen),       
C...   and ECP routines by A. V. Mitin and C. van Wuellen.                      
C...                                                                            
C...   see http://www.cfour.de for the current version.                         
C...                                                                            
C...   This source code is provided under a written licence and may be          
C...   used, copied, transmitted, or stored only in accord with that            
C...   written licence.                                                         
C...                                                                            
C...   In particular, no part of the source code or compiled modules may        
C...   be distributed outside the research group of the licence holder.         
C...   This means also that persons (e.g. post-docs) leaving the research       
C...   group of the licence holder may not take any part of CFOUR,              
C...   including modified files, with him/her, unless that person has           
C...   obtained his/her own licence.                                            
C...                                                                            
C...   For questions concerning this copyright write to:                        
C...                          info@cfour.de                                     
C...                                                                            
CLIC                                                                            
                                                                                
      SUBROUTINE OUTPUT (AMATRX,ROWLOW,ROWHI,COLLOW,COLHI,ROWDIM,COLDIM,                                                
     *                   NCTL,LUOUT)                                                                                    
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                                                               
      INTEGER   ROWLOW,ROWHI,COLLOW,COLHI,ROWDIM,COLDIM,BEGIN,KCOL                                                      
      DIMENSION AMATRX(ROWDIM,COLDIM)                                                                                   
      CHARACTER*1 ASA(3), BLANK, CTL                                                                                    
      CHARACTER   PFMT*20, COLUMN*8                                                                                     
      PARAMETER (ZERO=0.D00, KCOLP=4, KCOLN=6)                                                                          
      PARAMETER (FFMIN=1.D-3, FFMAX = 1.D3)                                                                             
      DATA COLUMN/'Column  '/, BLANK/' '/, ASA/' ', '0', '-'/                                                           
      IF (ROWHI.LT.ROWLOW) GO TO 3                                                                                      
      IF (COLHI.LT.COLLOW) GO TO 3                                                                                      
      AMAX = ZERO                                                                                                       
      DO 10 J = COLLOW,COLHI                                                                                            
         DO 10 I = ROWLOW,ROWHI                                                                                         
            AMAX = MAX( AMAX, ABS(AMATRX(I,J)) )                                                                        
   10 CONTINUE                                                                                                          
      IF (AMAX .EQ. ZERO) THEN                                                                                          
         WRITE (LUOUT,'(/T6,A)') 'Zero matrix.'                                                                         
         GO TO 3                                                                                                        
      END IF                                                                                                            
      IF (FFMIN .LE. AMAX .AND. AMAX .LE. FFMAX) THEN                                                                   
         PFMT = '(A1,I7,2X,8F15.8)'                                                                                     
      ELSE                                                                                                              
         PFMT = '(A1,I7,2X,1P8D15.6)'                                                                                   
      END IF                                                                                                            
      IF (NCTL .LT. 0) THEN                                                                                             
         KCOL = KCOLN                                                                                                   
      ELSE                                                                                                              
         KCOL = KCOLP                                                                                                   
      END IF                                                                                                            
      MCTL = ABS(NCTL)                                                                                                  
      IF ((MCTL.LE.3).AND.(MCTL.GT.0)) THEN                                                                             
         CTL = ASA(MCTL)                                                                                                
      ELSE                                                                                                              
         CTL = BLANK                                                                                                    
      END IF                                                                                                            
      LAST = MIN(COLHI,COLLOW+KCOL-1)                                                                                   
      DO 2 BEGIN = COLLOW,COLHI,KCOL                                                                                    
         WRITE (LUOUT,1000) (COLUMN,I,I = BEGIN,LAST)                                                                   
         DO 1 K = ROWLOW,ROWHI                                                                                          
            DO 4 I = BEGIN,LAST                                                                                         
               IF (AMATRX(K,I).NE.ZERO) GO TO 5                                                                         
    4       CONTINUE                                                                                                    
         GO TO 5
    5       WRITE (LUOUT,PFMT) CTL,K,(AMATRX(K,I), I = BEGIN,LAST)                                                      
    1    CONTINUE                                                                                                       
    2 LAST = MIN(LAST+KCOL,COLHI)                                                                                       
    3 RETURN                                                                                                            
 1000 FORMAT (/12X,6(3X,A6,I4,2X),(3X,A6,I4))                                                                           
      END                                                                                                               

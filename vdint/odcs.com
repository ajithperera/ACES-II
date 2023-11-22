#ifndef _ODCS_COM_
#define _ODCS_COM_

      COMMON/ODCS/ODC10X(mxqn*mxqn*(2*mxqn+4)),
     &            ODC20X(mxqn*mxqn*(2*mxqn+4)),
     &            ODC30X(mxqn*mxqn*(2*mxqn+4)), 
     &            ODC40X(mxqn*mxqn*(2*mxqn+4)),
     &            ODC00X(mxqn*(mxqn+5)*(2*mxqn+4)),
     &            ODC10Y(mxqn*mxqn*(2*mxqn+4)),
     &            ODC20Y(mxqn*mxqn*(2*mxqn+4)), 
     &            ODC30Y(mxqn*mxqn*(2*mxqn+4)), 
     &            ODC40Y(mxqn*mxqn*(2*mxqn+4)),
     &            ODC00Y(mxqn*(mxqn+5)*(2*mxqn+4)),
     &            ODC10Z(mxqn*mxqn*(2*mxqn+4)),
     &            ODC20Z(mxqn*mxqn*(2*mxqn+4)),
     &            ODC30Z(mxqn*mxqn*(2*mxqn+4)), 
     &            ODC40Z(mxqn*mxqn*(2*mxqn+4)),
     &            ODC00Z(mxqn*(mxqn+5)*(2*mxqn+4)) 

#endif  /* _ODCS_COM */



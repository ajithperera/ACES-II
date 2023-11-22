/* ecp_printstats.cpp: print out timing and related information at end of job
 *
 * Written by Tom Grimes, 26-March 2010
 *
 */

#include "ecpints.h"

extern "C" {

void ecp_printstats_()
{
   // initialize output data connection
   FILE *outf = fdopen(output_fd,"a");
   if (outf == NULL) { printf("Error connecting to output file!\n"); exit(1); };

   // print the timer stats
   ecp_data.time.print("Overall timing information for ECPs:",outf);
   fprintf(outf,"\n");
   fflush(outf);
   return;
};

};

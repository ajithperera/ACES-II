/* timer.cpp - a class that implements a simple function timer
 *
 * Written by Tom Grimes, 25-March 2010
 *
 */

#include "ecpints.h"

// initialize a key (might be unnecessary)
// this can also reset an element
void timer::init(string key)
{
   data[key]=(mytime) { tbb::tick_count(), false, 0., 0., 0., 0 };
   return;
};

// punch the clock, assuming the key exists
void timer::punch(string key)
{
   tbb::tick_count curr_t = tbb::tick_count::now();

   if (data.find(key)==data.end()) data[key]=(mytime) { tbb::tick_count(), false, 0., 0., 0., 0 };

   mytime *tptr = &data[key];
   if (tptr->running)
   {
      tptr->running=false;
      double tdiff=(curr_t-tptr->mark).seconds();
      if (tdiff>tptr->max) tptr->max=tdiff;
      if ((tdiff<tptr->min)||(tptr->N==0)) tptr->min=tdiff;
      tptr->cum+=tdiff;
      tptr->N++;
   } else {
      tptr->mark=curr_t;
      tptr->running=true;
   };
   return;
};

// print out the stats
void timer::print(string header, FILE* outf)
{
   // close out running timers
   for (map<string,mytime>::iterator iter=data.begin(); iter!=data.end(); iter++)
   {
      if (iter->second.running) this->punch(iter->first);
   };

   // print header
   fprintf(outf,"  %s\n",header.c_str());
   fprintf(outf,"  Timing in seconds, resolution approx. 1e-6 s\n",header.c_str());
   fprintf(outf,"                       cumulative      max         min         avg     N\n");
   for (map<string,mytime>::iterator iter=data.begin(); iter!=data.end(); iter++)
   {
      fprintf(outf,"  %-20s % 5.3e  % 5.3e  % 5.3e  % 5.3e  %i\n",iter->first.c_str(),iter->second.cum,iter->second.max,
              iter->second.min,iter->second.cum/iter->second.N,iter->second.N);
   };
   fflush(outf);

   return;
};

// fortran access routines
extern "C" {

// create a timer and pass back the address
// int means 32-bit dependence!
void get_timer_(int &ptr)
{
   ptr = reinterpret_cast<int>(new timer());
   return;
};

// punch a given entry
void punch_timer_(int &ptr, const char *key, unsigned int keylen)
{
   reinterpret_cast<timer*>(ptr)->punch(string(key,keylen));
   return;
};

// print and release a timer
// note that the fd is stolen from other ECP code (c_open_output.cpp)
void destroy_timer_(int &ptr, const char *header, unsigned int headerlen)
{
   FILE *outf = fdopen(output_fd,"a");
   if (outf != NULL)
   {
      reinterpret_cast<timer*>(ptr)->print(string(header,headerlen),outf);
      fprintf(outf,"\n");
   } else {
      reinterpret_cast<timer*>(ptr)->print(string(header,headerlen));
      printf("\n");
   };

   delete reinterpret_cast<timer*>(ptr);
   return;
};

};


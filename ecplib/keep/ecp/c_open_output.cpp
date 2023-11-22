/* c_open_output.cpp: saves the file descriptor for the output when opened */

#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/file.h>
#include <errno.h>

extern "C" {

int output_fd = -1;

int c_open_output_ ( char *file_name, /* access read: name of the file to open (null terminated) */
                     int *open_flags, /* access read: READ/WRITE, see file.h or open(2) */
                     int *create_mode, /* access read: set if the file is to be created */
                     int *unit_num, /* access read: logical unit number to be opened */
                     int filenam_len ) /* access read: number of characters in file_name */
{

   /*
    *
    * ** The returned value is the following:
    * ** value >= 0 is a valid file descriptor
    * ** value < 0 is an error
    * */

   int return_value;

   /*
   printf(" %s: Opening FILENAME = %s\n", __FILE__, file_name);
   printf(" %s: open_flags = 0x%8.8x\n", __FILE__, *open_flags);

   if ( *open_flags & O_CREAT )
      printf(" %s: the file is being created, create_mode = 0x%8.8x\n", __FILE__, *create_mode);

   printf(" %s: open() ", __FILE__);
   */

   return_value = open(file_name, *open_flags, *create_mode);
   output_fd=return_value;

   /*
   if (return_value < 0) { printf("FAILED.\n"); }
   else { printf("SUCCEEDED.\n"); };
   */

   return (return_value);

};

}; // extern C

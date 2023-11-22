
#ifndef _CALLSTACK_COM_
#define _CALLSTACK_COM_

c This contains variables used in the callstack.  See the callstack_ops
c file for a description of the callstack.
c
c callstack_num : the number of routines on the stack
c callstack_len : the length of each element on the stack
c callstack     : the strings containing the callstack
c callstack_ptr : the first unused element in the callstack

c callstack.com
      integer callstack_num,callstack_len
      parameter (callstack_num=100)
      parameter (callstack_len=20)
      character*(callstack_len) callstack_curr
      common /callstack_curr/ callstack_curr

      character*(callstack_len) callstack(callstack_num)
      integer callstack_ptr

      common /callstackc/ callstack
      common /callstackp/  callstack_ptr
      save /callstackc/
      save /callstackp/

#endif /* _CALLSTACK_COM_ */


/** \file  ListExample.cpp
    \brief Test and demonstration program.
           Copyright 2005 by Erik Schlögl

		   Redistribution and use in source and binary forms, with or without modification, are permitted provided
		   that the following conditions are met:
		
           -# Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

           -# Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

           -# The name of the copyright holder may not be used to endorse or promote products derived from this software without specific prior written permission.
		
           THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
		   INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
		   PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY DIRECT, INDIRECT,
		   INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
		   OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
		   CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
		   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
		   POSSIBILITY OF SUCH DAMAGE.
*/

#include <iostream>
#include <list>

int main(int argc,char *argv[])
{
  int i;
  // Declare list.
  std::list<int> my_list;
  // Fill with integers 1 to 10.
  for (i=1;i<=10;i++) my_list.push_back(i);
  // Declare list iterator.
  std::list<int>::iterator my_iter;
  // Print list.
  for (my_iter=my_list.begin();my_iter!=my_list.end();my_iter++) 
	  std::cout << *my_iter << ' ' << std::flush;
  std::cout << std::endl;
  // Insert the first program argument as the third element in the list.
  int j = 0;
  if (argc>1) j = atoi(argv[1]);
  my_iter=my_list.begin();
  for (i=0;i<2;i++) my_iter++;
  my_list.insert(my_iter,j);
  // Print list.
  for (my_iter=my_list.begin();my_iter!=my_list.end();my_iter++) 
	  std::cout << *my_iter << ' ' << std::flush;
  std::cout << std::endl;
  // Erase the fifth element in the list.
  my_iter=my_list.begin();
  for (i=0;i<4;i++) my_iter++;
  my_list.erase(my_iter);
  // Print list.
  for (my_iter=my_list.begin();my_iter!=my_list.end();my_iter++) 
	  std::cout << *my_iter << ' ' << std::flush;
  std::cout << std::endl;
}

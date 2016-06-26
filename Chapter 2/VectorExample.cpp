/** \file  VectorExample.cpp
    \brief Test and demonstration program.
           Copyright 2013 by Erik Schlögl

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
#include <vector>

int main(int argc,char *argv[])
{
  int i;
  // Instantiate vector
  std::vector<int> my_vector(10);
  // Fill with integers 1 to 10.
  for(i=1;i<=10;i++) my_vector[i-1] = i;
  // and print all the elements in a vector by
  for(i=0;i<my_vector.size();i++)
    std::cout << my_vector[i] << ' ' << std::flush;  
  std::cout << std::endl;
  // Declare vector iterator.
  std::vector<int>::iterator my_iter;
  // Print vector.
  for (my_iter=my_vector.begin();my_iter!=my_vector.end();my_iter++) 
	  std::cout << *my_iter << ' ' << std::flush;
  std::cout << std::endl;
  // Insert the first program argument as the third element in the vector.
  int j = 0;
  if (argc>1) j = atoi(argv[1]);
  my_iter=my_vector.begin();
  for (i=0;i<2;i++) my_iter++;
  my_vector.insert(my_iter,j);
  // Print vector.
  for (my_iter=my_vector.begin();my_iter!=my_vector.end();my_iter++) 
	  std::cout << *my_iter << ' ' << std::flush;
  std::cout << std::endl;
  // Erase the fifth element in the vector.
  my_iter=my_vector.begin();
  for (i=0;i<4;i++) my_iter++;
  my_vector.erase(my_iter);
  // Print vector.
  for (my_iter=my_vector.begin();my_iter!=my_vector.end();my_iter++) 
	  std::cout << *my_iter << ' ' << std::flush;
  std::cout << std::endl;
}

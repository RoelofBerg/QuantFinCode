/** \file  for_eachExample.cpp
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
#include <map>
#include <algorithm>

template <typename T>
class number_of_occurences {
private:
  std::map<T,int> count_occurences;
public:
  inline void operator()(const T& x) { count_occurences[x]++; };
  inline operator std::map<T,int>&() { return count_occurences; };
};

int main(int argc,char *argv[])
{
  int i;
  // Declare list.
  std::list<int> my_list;
  // Fill with integers 1 to 10.
  for (i=1;i<=10;i++) my_list.push_back(i);
  // append another 5 and another 3
  my_list.push_back(5);
  my_list.push_back(3);
  // instantiate counting functor
  number_of_occurences<int> f;
  // count occurences of each integer using the for_each algorithm template
  number_of_occurences<int> g = for_each(my_list.begin(),my_list.end(),f);
  std::map<int,int>& occurences = g;
  // Print number of occurences of each integer
  std::map<int,int>::iterator my_iter;
  for (my_iter=occurences.begin();my_iter!=occurences.end();my_iter++)
      std::cout << "The integer " << my_iter->first << " occurs " << my_iter->second << " times in the list." << std::endl;
}

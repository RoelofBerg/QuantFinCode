/** \file  MapExample.cpp
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
#include <string>
#include <map>

int main(int argc,char *argv[])
{
  // Create an empty map
  std::map<std::string,double> my_map;
  // Add stocks and their prices to the map
  my_map["Mended Hill Properties"] = 25.2;
  my_map["Yoyodyne Industries"] = 105.1;
  my_map["Peak Petroleum Corporation"] = 803.5;
  my_map["Titanic Shipping Ltd."] = 3.15;
  // Elements of the map can be accessed through their keys
  std::cout << "The share price of Yoyodyne Industries is " << my_map["Yoyodyne Industries"] << std::endl;
  // Print elements in order
  std::map<std::string,double>::iterator my_iter;
  for (my_iter=my_map.begin();my_iter!=my_map.end();my_iter++)
	  std::cout << "The share price of " << my_iter->first << " is " << my_iter->second << std::endl;
}

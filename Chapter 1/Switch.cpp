/** \file  Switch.cpp
    \brief Example program for a "switch" statement.
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

#include <cstdlib>
#include <iostream>

using namespace std;

int main(int argc, char *argv[])
{
    int c;
    cout << "Enter an integer between 1 and 7: " << endl;
    cin >> c;
    switch (c) {
      case 1: cout << "Sunday" << endl;
              break; 
      case 2: cout << "Monday" << endl;
              break; 
      case 3: cout << "Tuesday" << endl;
              break; 
      case 4: cout << "Wednesday" << endl;
              break; 
      case 5: cout << "Thursday" << endl;
              break; 
      case 6: cout << "Friday" << endl;
              break; 
      case 7: cout << "Saturday" << endl;
              break; 
      default: cout << "Invalid input" << endl;
              break;                            }
    return EXIT_SUCCESS;
}

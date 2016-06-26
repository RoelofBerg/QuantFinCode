/** \file  Static.cpp
    \brief Example program for function calls.
           Copyright 2013 by Erik Schl�gl

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

int static_variable_function(int x)
{
	static int previous_value;
	int result = previous_value;
    previous_value = x;
    return result;
}

int main(int argc, char *argv[])
{
    int y;
    int result;
    
    y = 1;
    result = static_variable_function(y);
    cout << "After first function call: " << endl;
    cout << "y is " << y << endl;
    cout << "The function returned " << result << endl << endl;
    
    y = 2;
    result = static_variable_function(y);
    cout << "After second function call: " << endl;
    cout << "y is " << y << endl;
    cout << "The function returned " << result << endl << endl;

    y = 3;
    result = static_variable_function(y);
    cout << "After third function call: " << endl;
    cout << "y is " << y << endl;
    cout << "The function returned " << result << endl << endl;

    return EXIT_SUCCESS;
}

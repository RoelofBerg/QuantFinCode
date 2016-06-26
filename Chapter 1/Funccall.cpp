/** \file  Funccall.cpp
    \brief Example program for function calls.
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

int add_one_by_value(int x)
{
    x = x + 1;
    return x;
}

int add_one_by_reference(int& x)
{
    x = x + 1;
    return x;
}

int add_one_by_pointer(int* x)
{
    *x = *x + 1;
    return *x;
}

int main(int argc, char *argv[])
{
    int y;
    int result;
    
    y = 1;
    result = add_one_by_value(y);
    cout << "After function call by value: " << endl;
    cout << "y is " << y << endl;
    cout << "The function returned " << result << endl << endl;
    
    y = 1;
    result = add_one_by_reference(y);
    cout << "After function call by reference: " << endl;
    cout << "y is " << y << endl;
    cout << "The function returned " << result << endl << endl;
    
    y = 1;
    result = add_one_by_pointer(&y);
    cout << "After function call by pointer: " << endl;
    cout << "y is " << y << endl;
    cout << "The function returned " << result << endl << endl;

    return EXIT_SUCCESS;
}

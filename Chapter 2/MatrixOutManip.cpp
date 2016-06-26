/** \file  MatrixIdx.cpp
    \brief Example program demonstrating usage of overloaded
           operator () for subscripting.
           Copyright 2005 by Erik Schloegl 
     */

#include <iostream>
#include <iomanip>
#include <cmath>
#include "Matrix.h"

using namespace std;

int main(int argc, char *argv[])
{
    int i,j;
    Matrix A(3,4),B(3,4);
    // assign values to matrix elements
    for (i=0;i<3;i++) {
      for (j=0;j<4;j++) {
		A(i,j) = 0.1*i*j; 
		B(i,j) = std::sqrt(1.0 - A(i,j)); }}
	cout << "Matrix B:\n" << std::setprecision(12) << B;

    return EXIT_SUCCESS;
}

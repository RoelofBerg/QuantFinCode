/** \file  QFQuasiRandom.cpp
    \brief C++ source file implementing classes for quasi-random number generators.
           Copyright 2010, 2013 by Erik Schlögl (except as noted below)

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

    -----------------------------------------------------------------------------
    Licence pertaining to sobol.cc and the accompanying sets of direction numbers
    -----------------------------------------------------------------------------
    Copyright (c) 2008, Frances Y. Kuo and Stephen Joe
    All rights reserved.
    
    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:
    
        * Redistributions of source code must retain the above copyright
          notice, this list of conditions and the following disclaimer.
    
        * Redistributions in binary form must reproduce the above copyright
          notice, this list of conditions and the following disclaimer in the
          documentation and/or other materials provided with the distribution.
    
        * Neither the names of the copyright holders nor the names of the
          University of New South Wales and the University of Waikato
          and its contributors may be used to endorse or promote products derived
          from this software without specific prior written permission.
    
    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ``AS IS'' AND ANY
    EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
    WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS BE LIABLE FOR ANY
    DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
    (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
    LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
    ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
		   */

#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "QFQuasiRandom.hpp"

using namespace quantfin;

const char quantfin::default_direction_number_file[] = "new-joe-kuo-6.21201.dat";

Sobol::Sobol(int dimension,size_t maximum_number_of_points,const char* direction_number_file)
: maxN(maximum_number_of_points),current_point(dimension),currentX(dimension),dim(dimension)
{
  initialise(direction_number_file);
}

void Sobol::initialise(const char* direction_number_file)
{
  unsigned i,j,k;
  std::ifstream infile(direction_number_file,std::ios::in);
  if (!infile) throw std::runtime_error("Input file containing direction numbers cannot be found!");
  if (maxN>std::pow((double)2.0,(int)32)) throw std::runtime_error("Maximum number of quasirandom points cannot be greater than 2^32");
  // Clear first line in input file.
  char buffer[1000];
  infile.getline(buffer,1000,'\n');
  // L = max number of bits needed 
  unsigned L = (unsigned)std::ceil(std::log((double)maxN)/std::log(2.0)); 
  scaling  = std::pow(2.0,32);
  currentX = 0;
  current_point = 0.0;
  // For the first dimension, compute direction numbers V[1] to V[L], scaled by pow(2,32)
  boost::shared_ptr<Array<unsigned,2> > tmp(new Array<unsigned,2>(L+1,dim));
  V = tmp;
  for (i=1;i<=L;i++) (*V)((int)i,(int)0) = 1 << (32-i); // all m's = 1
  // ----- Compute the remaining dimensions -----
  for (j=1;j<=dim-1;j++) {
    // Read in parameters from file 
    unsigned d, s;
    unsigned a;
    infile >> d >> s >> a;
    Array<unsigned,1> m(s+1);
    for (i=1;i<=s;i++) infile >> m(i);
    // Compute direction numbers V[1] to V[L], scaled by pow(2,32)
    if (L <= s) {
      for (i=1;i<=L;i++) (*V)((int)i,(int)j) = m(i) << (32-i); }
    else {
      for (i=1;i<=s;i++) (*V)((int)i,(int)j) = m(i) << (32-i); 
      for (i=s+1;i<=L;i++) {
	    (*V)((int)i,(int)j) = (*V)((int)(i-s),(int)j) ^ ((*V)((int)(i-s),(int)j) >> s); 
		for (k=1;k<=s-1;k++) (*V)((int)i,(int)j) ^= (((a >> (s-1-k)) & 1) * (*V)((int)(i-k),(int)j)); }}}
  currentC = 1;
  currenti = 0;
  // discard initial point (zero)
  random();
}

Array<double,1>& Sobol::random()
{
  int j;
  if (currenti) {
	for (j=0;j<currentX.extent(blitz::firstDim);j++) currentX(j) = currentX(j) ^ (*V)((int)currentC,j);
    // ----- Create point -----
	for (j=0;j<currentX.extent(blitz::firstDim);j++) current_point(j) = ((double)(currentX(j)))/scaling; }
  // Update i and C
  currentC = 1;
  unsigned value = currenti;
  while (value & 1) {
    value >>= 1;
    currentC++; }
  currenti++;
  return current_point;
}

Array<double,1>& Sobol::random(const Array<double,1>& shift)
{
  int j;
  Array<double,1>& result = random();
  result += shift;
  for (j=0;j<result.extent(blitz::firstDim);j++) if (result(j)>1.0) result(j) -= 1.0;
  return result;
}

void Sobol::reset()
{
  currentX = 0;
  current_point = 0.0;
  currentC = 1;
  currenti = 0;
  // discard initial point (zero)
  random();
}

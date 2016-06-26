/** \file  CSV2ArrayExample.cpp
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
#include <fstream>
#include <cstdlib>
#include "CSV2Array.hpp"

using namespace quantfin;
 
/** Test and demonstration for reading CSV files.

    Command-line arguments:
	  -# option parameter CSV file name.
	  -# name of CSV file with numerical values (only).
  */
int main(int argc,char* argv[]) 
{
  using std::cout;
  using std::endl;
  using std::flush;
  using std::cerr;

  int i;
  double mat,strike,vol,S,r;
  try {
	// echo command line
    for (i=0;i<argc;i++) cout << argv[i] << ' ';
	cout << endl;
	if (argc<2) throw std::logic_error("Missing option parameter CSV file name");
	// Read data from file and initialise option parameters
    std::ifstream is_inputs(argv[1]);
    if (!is_inputs.is_open()) throw std::logic_error("Failed to open option parameter CSV file");
    blitz::Array<std::string,2> inputs_matrix(CSV2Array<std::string>(is_inputs,make_string_strip_spaces));
    std::map<std::string,std::string> inputs_map;
    for (i=0;i<inputs_matrix.rows();i++) inputs_map[inputs_matrix(i,0)] = inputs_matrix(i,1);
    std::string str("Option expiry");
    if (inputs_map.count(str)) mat = std::atof(inputs_map[str].data());
    else throw std::logic_error("Missing option expiry"); 
	str = "Strike";
    if (inputs_map.count(str)) strike = std::atof(inputs_map[str].data());
    else throw std::logic_error("Missing strike"); 
	str = "Volatility";
    if (inputs_map.count(str)) strike = std::atof(inputs_map[str].data());
    else throw std::logic_error("Missing volatility"); 
	str = "Initial underlying price";
    if (inputs_map.count(str)) strike = std::atof(inputs_map[str].data());
    else throw std::logic_error("Missing initial underlying price"); 
	str = "Interest rate";
    if (inputs_map.count(str)) strike = std::atof(inputs_map[str].data());
    else throw std::logic_error("Missing interest rate"); 
	if (argc<3) throw std::logic_error("Missing name of CSV file containing numerical values");
	// Read numerical data from file 
    std::ifstream is_num(argv[2]);
    if (!is_inputs.is_open()) throw std::logic_error("Failed to open CSV file containing numerical values");
    blitz::Array<double,2> numerical_matrix(CSV2Array<double>(is_num,std::atof));
	cout << "The numerical data is: " << numerical_matrix << endl;

  } // end of try block

  catch (std::logic_error xcpt) {
    std::cerr << xcpt.what() << endl; }
  catch (std::runtime_error xcpt) {
    std::cerr << xcpt.what() << endl; }
  catch (...) {
    std::cerr << "Other exception caught" << endl; }
  
  return 0;
}

//// 
//// Autogenerated by xlw 
//// Do not edit this file, it will be overwritten 
//// by InterfaceGenerator 
////

#include "xlw/MyContainers.h"
#include <xlw/CellMatrix.h>
#include "ESXLWProject.hpp"
#include <xlw/xlw.h>
#include <xlw/XlFunctionRegistration.h>
#include <stdexcept>
#include <xlw/XlOpenClose.h>
#include <xlw/HiResTimer.h>
using namespace xlw;

namespace {
const char* LibraryName = "EriksTestLibrary";
};


// registrations start here


namespace
{
  XLRegistration::XLFunctionRegistrationHelper
registerEriksEmptyArgFunction("xlEriksEmptyArgFunction",
"EriksEmptyArgFunction",
" tests empty args ",
LibraryName,
0,
0
,false
,false
,""
,""
,false
,false
,false
);
}



extern "C"
{
LPXLFOPER EXCEL_EXPORT
xlEriksEmptyArgFunction(
)
{
EXCEL_BEGIN;

	if (XlfExcel::Instance().IsCalledByFuncWiz())
		return XlfOper(true);

std::string result(
	EriksEmptyArgFunction());
return XlfOper(result);
EXCEL_END
}
}



//////////////////////////

namespace
{
XLRegistration::Arg
ErikEchoArrayArgs[]=
{
{ "Echoee"," argument to be echoed ","XLF_OPER"}
};
  XLRegistration::XLFunctionRegistrationHelper
registerErikEchoArray("xlErikEchoArray",
"ErikEchoArray",
" echoes an array ",
LibraryName,
ErikEchoArrayArgs,
1
,false
,false
,""
,""
,false
,false
,false
);
}



extern "C"
{
LPXLFOPER EXCEL_EXPORT
xlErikEchoArray(
LPXLFOPER Echoeea)
{
EXCEL_BEGIN;

	if (XlfExcel::Instance().IsCalledByFuncWiz())
		return XlfOper(true);

XlfOper Echoeeb(
	(Echoeea));
MyArray Echoee(
	Echoeeb.AsArray("Echoee"));

MyArray result(
	ErikEchoArray(
		Echoee)
	);
return XlfOper(result);
EXCEL_END
}
}



//////////////////////////

namespace
{
XLRegistration::Arg
ESLogLinearArgs[]=
{
{ "given_times"," timeline for supplied discount factors ","XLF_OPER"},
{ "given_discount_factors"," supplied discount factors ","XLF_OPER"},
{ "target_times"," times for which discount factors are desired ","XLF_OPER"}
};
  XLRegistration::XLFunctionRegistrationHelper
registerESLogLinear("xlESLogLinear",
"ESLogLinear",
" log-linear interpolation given set of discount factors ",
LibraryName,
ESLogLinearArgs,
3
,false
,false
,""
,""
,false
,false
,false
);
}



extern "C"
{
LPXLFOPER EXCEL_EXPORT
xlESLogLinear(
LPXLFOPER given_timesa,
LPXLFOPER given_discount_factorsa,
LPXLFOPER target_timesa)
{
EXCEL_BEGIN;

	if (XlfExcel::Instance().IsCalledByFuncWiz())
		return XlfOper(true);

XlfOper given_timesb(
	(given_timesa));
MyArray given_times(
	given_timesb.AsArray("given_times"));

XlfOper given_discount_factorsb(
	(given_discount_factorsa));
MyArray given_discount_factors(
	given_discount_factorsb.AsArray("given_discount_factors"));

XlfOper target_timesb(
	(target_timesa));
MyArray target_times(
	target_timesb.AsArray("target_times"));

MyArray result(
	ESLogLinear(
		given_times,
		given_discount_factors,
		target_times)
	);
return XlfOper(result);
EXCEL_END
}
}



//////////////////////////

namespace
{
XLRegistration::Arg
ESLogLinearMatrixArgs[]=
{
{ "given"," timeline for supplied discount factors (first column) and supplied discount factors (second column) ","XLF_OPER"},
{ "target_times"," times for which discount factors are desired ","XLF_OPER"}
};
  XLRegistration::XLFunctionRegistrationHelper
registerESLogLinearMatrix("xlESLogLinearMatrix",
"ESLogLinearMatrix",
" log-linear interpolation given set of discount factors ",
LibraryName,
ESLogLinearMatrixArgs,
2
,false
,false
,""
,""
,false
,false
,false
);
}



extern "C"
{
LPXLFOPER EXCEL_EXPORT
xlESLogLinearMatrix(
LPXLFOPER givena,
LPXLFOPER target_timesa)
{
EXCEL_BEGIN;

	if (XlfExcel::Instance().IsCalledByFuncWiz())
		return XlfOper(true);

XlfOper givenb(
	(givena));
MyMatrix given(
	givenb.AsMatrix("given"));

XlfOper target_timesb(
	(target_timesa));
MyArray target_times(
	target_timesb.AsArray("target_times"));

MyArray result(
	ESLogLinearMatrix(
		given,
		target_times)
	);
return XlfOper(result);
EXCEL_END
}
}



//////////////////////////

namespace
{
XLRegistration::Arg
ESBlackScholesArgs[]=
{
{ "S"," current price of underlying asset ","B"},
{ "K"," exercise price ","B"},
{ "sigma"," volatility of underlying asset ","B"},
{ "r"," risk-free interest rate (continuously compounded) ","B"},
{ "T"," time to maturity (in years) ","B"},
{ "callput"," 1 for call, -1 for put ","B"}
};
  XLRegistration::XLFunctionRegistrationHelper
registerESBlackScholes("xlESBlackScholes",
"ESBlackScholes",
" Black/Scholes European option price ",
LibraryName,
ESBlackScholesArgs,
6
,false
,false
,""
,""
,false
,false
,false
);
}



extern "C"
{
LPXLFOPER EXCEL_EXPORT
xlESBlackScholes(
double S,
double K,
double sigma,
double r,
double T,
double callputa)
{
EXCEL_BEGIN;

	if (XlfExcel::Instance().IsCalledByFuncWiz())
		return XlfOper(true);






int callput(
	static_cast<int>(callputa));

double result(
	ESBlackScholes(
		S,
		K,
		sigma,
		r,
		T,
		callput)
	);
return XlfOper(result);
EXCEL_END
}
}



//////////////////////////

namespace
{
XLRegistration::Arg
ESCoxRossRubinsteinArgs[]=
{
{ "S"," current price of underlying asset ","B"},
{ "K"," exercise price ","B"},
{ "sigma"," volatility of underlying asset ","B"},
{ "r"," risk-free interest rate (continuously compounded) ","B"},
{ "T"," time to maturity (in years) ","B"},
{ "callput"," 1 for call, -1 for put ","B"},
{ "N"," number of steps in the binomial lattice ","B"}
};
  XLRegistration::XLFunctionRegistrationHelper
registerESCoxRossRubinstein("xlESCoxRossRubinstein",
"ESCoxRossRubinstein",
" Cox/Ross/Rubinstein European option price ",
LibraryName,
ESCoxRossRubinsteinArgs,
7
,false
,false
,""
,""
,false
,false
,false
);
}



extern "C"
{
LPXLFOPER EXCEL_EXPORT
xlESCoxRossRubinstein(
double S,
double K,
double sigma,
double r,
double T,
double callputa,
double Na)
{
EXCEL_BEGIN;

	if (XlfExcel::Instance().IsCalledByFuncWiz())
		return XlfOper(true);






int callput(
	static_cast<int>(callputa));

int N(
	static_cast<int>(Na));

double result(
	ESCoxRossRubinstein(
		S,
		K,
		sigma,
		r,
		T,
		callput,
		N)
	);
return XlfOper(result);
EXCEL_END
}
}



//////////////////////////

namespace
{
XLRegistration::Arg
ESCoxRossRubinsteinAmericanArgs[]=
{
{ "S"," current price of underlying asset ","B"},
{ "K"," exercise price ","B"},
{ "sigma"," volatility of underlying asset ","B"},
{ "r"," risk-free interest rate (continuously compounded) ","B"},
{ "T"," time to maturity (in years) ","B"},
{ "N"," number of steps in the binomial lattice ","B"}
};
  XLRegistration::XLFunctionRegistrationHelper
registerESCoxRossRubinsteinAmerican("xlESCoxRossRubinsteinAmerican",
"ESCoxRossRubinsteinAmerican",
" Cox/Ross/Rubinstein American put option price ",
LibraryName,
ESCoxRossRubinsteinAmericanArgs,
6
,false
,false
,""
,""
,false
,false
,false
);
}



extern "C"
{
LPXLFOPER EXCEL_EXPORT
xlESCoxRossRubinsteinAmerican(
double S,
double K,
double sigma,
double r,
double T,
double Na)
{
EXCEL_BEGIN;

	if (XlfExcel::Instance().IsCalledByFuncWiz())
		return XlfOper(true);






int N(
	static_cast<int>(Na));

double result(
	ESCoxRossRubinsteinAmerican(
		S,
		K,
		sigma,
		r,
		T,
		N)
	);
return XlfOper(result);
EXCEL_END
}
}



//////////////////////////

//////////////////////////
// Methods that will get registered to execute in AutoOpen
//////////////////////////

//////////////////////////
// Methods that will get registered to execute in AutoClose
//////////////////////////

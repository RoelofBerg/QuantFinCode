#include "TestObject.hpp"
#include "TSObject.hpp"
#include <oh/enumerations/typefactory.hpp>
#include <ohxl/objecthandlerxl.hpp>
#include <ohxl/register/register_all.hpp>
#include <ohxl/functions/export.hpp>
#include <ohxl/utilities/xlutilities.hpp>
#include <ohxl/objectwrapperxl.hpp>
#include <ohxl/callingrange.hpp>

#ifdef BOOST_MSVC
#  define BOOST_LIB_DIAGNOSTIC
#  include <oh/auto_link.hpp>
#  undef BOOST_LIB_DIAGNOSTIC
#endif
#include <sstream>

DLLEXPORT int xlAutoOpen() {

    // Instantiate the ObjectHandler Repository
    static ObjectHandler::RepositoryXL repositoryXL;
    // Instantiate the Enumerated Type Registry
    static ObjectHandler::EnumTypeRegistry enumTypeRegistry;
    // Instantiate the Enumerated Class Registry
    static ObjectHandler::EnumClassRegistry enumClassRegistry;
    // Instantiate the Enumerated Pair Registry
    static ObjectHandler::EnumPairRegistry enumPairRegistry;
	//Instantiate the Processor Factory
	static ObjectHandler::ProcessorFactory processorFactory;
    // Instantiate the Serialization Factory
    //static AccountExample::SerializationFactory factory;

    static XLOPER xDll;

    try {
		ObjectHandler::logSetFile("./example.log");

        Excel(xlGetName, &xDll, 0);

        ObjectHandler::Configuration::instance().init();

        registerOhFunctions(xDll);

		/*
        Parameter codes for function registration:
		  C		char*
		  N		long*
		  E		double*
		  L		bool*
		  P		OPER*
		*/

        Excel(xlfRegister, 0, 7, &xDll,
            TempStrNoSize("\x0E""ohESTestObject"),          // function code name
            TempStrNoSize("\x05""CCNE#"),              // parameter codes: char* (char*,long*,double*)
            TempStrNoSize("\x0E""ohESTestObject"),          // function display name
            TempStrNoSize("\x17""ObjectID,integer,double"), // comma-delimited list of parameters
            TempStrNoSize("\x01""1"),                   // function type (0 = hidden function, 1 = worksheet function, 2 = command macro)
            TempStrNoSize("\x07""Example"));            // function category

        Excel(xlfRegister, 0, 7, &xDll,
            TempStrNoSize("\x19""ohESTestObjectQueryDouble"),    // function code name
            TempStrNoSize("\x03""EC#"),                // parameter codes
            TempStrNoSize("\x19""ohESTestObjectQueryDouble"),    // function display name
            TempStrNoSize("\x08""ObjectID"),    // comma-delimited list of parameters
            TempStrNoSize("\x01""1"),                   // function type (0 = hidden function, 1 = worksheet function, 2 = command macro)
            TempStrNoSize("\x07""Example"));            // function category

        Excel(xlfRegister, 0, 7, &xDll,
            TempStrNoSize("\x1A""ohESTestObjectQueryInteger"),    // function code name
            TempStrNoSize("\x03""NC#"),                // parameter codes
            TempStrNoSize("\x1A""ohESTestObjectQueryInteger"),    // function display name
            TempStrNoSize("\x08""ObjectID"),    // comma-delimited list of parameters
            TempStrNoSize("\x01""1"),                   // function type (0 = hidden function, 1 = worksheet function, 2 = command macro)
            TempStrNoSize("\x07""Example"));            // function category

        Excel(xlfRegister, 0, 7, &xDll,
            TempStrNoSize("\x0B""TSLogLinear"),    // function code name
            TempStrNoSize("\x05""CCPP#"),                // parameter codes
            TempStrNoSize("\x0B""TSLogLinear"),    // function display name
            TempStrNoSize("\x22""ObjectID,timeline,discount factors"),    // comma-delimited list of parameters
            TempStrNoSize("\x01""1"),                   // function type (0 = hidden function, 1 = worksheet function, 2 = command macro)
            TempStrNoSize("\x07""Example"));            // function category

        Excel(xlfRegister, 0, 7, &xDll,
            TempStrNoSize("\x0D""TSLogLinearDF"),    // function code name
            TempStrNoSize("\x04""ECE#"),                // parameter codes
            TempStrNoSize("\x0D""TSLogLinearDF"),    // function display name
            TempStrNoSize("\x11""ObjectID,maturity"),    // comma-delimited list of parameters
            TempStrNoSize("\x01""1"),                   // function type (0 = hidden function, 1 = worksheet function, 2 = command macro)
            TempStrNoSize("\x07""Example"));            // function category

        Excel(xlFree, 0, 1, &xDll);
        return 1;

    } catch (const std::exception &e) {

        std::ostringstream err;
        err << "Error loading ExampleXllStatic: " << e.what();
        Excel(xlcAlert, 0, 1, TempStrStl(err.str()));
        Excel(xlFree, 0, 1, &xDll);
        return 0;

    } catch (...) {

        Excel(xlFree, 0, 1, &xDll);
        return 0;

    }
}

DLLEXPORT int xlAutoClose() {

    static XLOPER xDll;

    try {

        Excel(xlGetName, &xDll, 0);

        unregisterOhFunctions(xDll);

        Excel(xlFree, 0, 1, &xDll);
        return 1;

    } catch (const std::exception &e) {

        std::ostringstream err;
        err << "Error unloading ExampleXllStatic: " << e.what();
        Excel(xlcAlert, 0, 1, TempStrStl(err.str()));
        Excel(xlFree, 0, 1, &xDll);
        return 0;

    } catch (...) {

        Excel(xlFree, 0, 1, &xDll);
        return 0;

    }

}

DLLEXPORT void xlAutoFree(XLOPER *px) {

    freeOper(px);

}

// function to instantiate ESTestObject in the repository
DLLEXPORT char *ohESTestObject(
        char *objectID,
        long *input_int,
        double *input_double) {

    boost::shared_ptr<ObjectHandler::FunctionCall> functionCall;

    try {

        functionCall = boost::shared_ptr<ObjectHandler::FunctionCall>
            (new ObjectHandler::FunctionCall("ohESTestObject"));

        boost::shared_ptr<ObjectHandler::ValueObject> valueObject(
            new ESTestValueObject(objectID,*input_int,*input_double));

        boost::shared_ptr<ObjectHandler::Object> object(
            new ESTestObject(valueObject,*input_int,*input_double));

        std::string returnValue = 
            ObjectHandler::RepositoryXL::instance().storeObject(objectID, object);

        static char ret[XL_MAX_STR_LEN];
        ObjectHandler::stringToChar(returnValue, ret);
        return ret;

    } catch (const std::exception &e) {

        ObjectHandler::RepositoryXL::instance().logError(e.what(), functionCall);
        return 0;

    }
}

DLLEXPORT double *ohESTestObjectQueryDouble(char *objectID, OPER *trigger) {

    boost::shared_ptr<ObjectHandler::FunctionCall> functionCall;

    try {

        functionCall = boost::shared_ptr<ObjectHandler::FunctionCall>
            (new ObjectHandler::FunctionCall("ohESTestObjectQueryDouble"));

        OH_GET_OBJECT(testObject,     // name of variable to hold a shared_ptr to object
			          objectID,       // object ID
					  ESTestObject    // C++ class of object
					  )

        static double ret;
        ret = testObject->double_value();
        return &ret;

    } catch (const std::exception &e) {

        ObjectHandler::RepositoryXL::instance().logError(e.what(), functionCall);
        return 0;

    }
}

DLLEXPORT long *ohESTestObjectQueryInteger(char *objectID, OPER *trigger) {

    boost::shared_ptr<ObjectHandler::FunctionCall> functionCall;

    try {

        functionCall = boost::shared_ptr<ObjectHandler::FunctionCall>
            (new ObjectHandler::FunctionCall("ohESTestObjectQueryInteger"));

        OH_GET_OBJECT(testObject,     // name of variable to hold a shared_ptr to object
			          objectID,       // object ID
					  ESTestObject    // C++ class of object
					  )

        static long ret;
        ret = testObject->integer_value();
        return &ret;

    } catch (const std::exception &e) {

        ObjectHandler::RepositoryXL::instance().logError(e.what(), functionCall);
        return 0;

    }
}

// function to create term structure object with loglinear interpolation of ZCB prices
DLLEXPORT char *TSLogLinear(
        char *ObjectId,
        OPER *xT,
        OPER *xB) {
	int i;

    // declare a shared pointer to the Function Call object

    boost::shared_ptr<ObjectHandler::FunctionCall> functionCall;

    try {

        // instantiate the Function Call object
    
        functionCall = boost::shared_ptr<ObjectHandler::FunctionCall>(
            new ObjectHandler::FunctionCall("TSLogLinear"));

        // convert input datatypes to C++ datatypes

        std::vector<double> vecT = ObjectHandler::operToVector<double>(*xT,"timeline");          // second argument is only relevant for outputting error msg if conversion fails
        std::vector<double> vecB = ObjectHandler::operToVector<double>(*xB,"discount factors");
		if (vecT.size()!=vecB.size()) throw std::logic_error("Timeline and discount factors must have the same number of entries");
		blitz::Array<double,1> T(vecT.size()),B(vecB.size());
        for (i=0;i<vecT.size();i++) T(i) = vecT[i];
        for (i=0;i<vecB.size();i++) B(i) = vecB[i];
		boost::shared_ptr<quantfin::TSLogLinear> ts(new quantfin::TSLogLinear(T,B));

        // Strip the Excel cell update counter suffix from Object IDs
        
        std::string ObjectIdStrip = ObjectHandler::CallingRange::getStub(ObjectId);

        // Construct the Value Object

        boost::shared_ptr<ObjectHandler::ValueObject> valueObject(new ESTSValueObject(ObjectId));

        // Construct the Object
        
        boost::shared_ptr<ObjectHandler::Object> object(new ESTSObject(valueObject,ts));

        // Store the Object in the Repository

        std::string returnValue =
            ObjectHandler::RepositoryXL::instance().storeObject(ObjectIdStrip, object);
            //ObjectHandler::RepositoryXL::instance().storeObject(ObjectIdStrip, object, *Overwrite);

        // Convert and return the return value

        static char ret[XL_MAX_STR_LEN];
        ObjectHandler::stringToChar(returnValue, ret);
        return ret;

    } catch (const std::exception &e) {
        ObjectHandler::RepositoryXL::instance().logError(e.what(), functionCall);
        return 0;
    }
}

// function to query discount factor for maturity t from existing TSLogLinear object
DLLEXPORT double *TSLogLinearDF(char *objectID,double* t) {

    boost::shared_ptr<ObjectHandler::FunctionCall> functionCall;

    try {

        functionCall = boost::shared_ptr<ObjectHandler::FunctionCall>
            (new ObjectHandler::FunctionCall("TSLogLinearDF"));

        OH_GET_OBJECT(ts_object,     // name of variable to hold a shared_ptr to object
			          objectID,       // object ID
					  ESTSObject    // C++ class of object
					  )

        static double ret;
        ret = ts_object->discount_factor(*t);
        return &ret;

    } catch (const std::exception &e) {

        ObjectHandler::RepositoryXL::instance().logError(e.what(), functionCall);
        return 0;

    }
}


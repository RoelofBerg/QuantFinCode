#include <oh/libraryobject.hpp>
#include <oh/valueobject.hpp>

class ESTestVal {
private:
  long ival;
  double dval;
public: 
  ESTestVal(const long& integer_input,const double& double_input) : ival(integer_input),dval(double_input) { };
  const double &double_value() { return dval; };
  const long &integer_value() { return ival; };
};

class ESTestValueObject : public ObjectHandler::ValueObject {
private:
  static const char* mPropertyNames[];
  long ival;
  double dval;
public:
  ESTestValueObject() {};
  ESTestValueObject(const std::string& objectID,const long& integer_input,const double& double_input)
            : ObjectHandler::ValueObject(objectID,"ESTestObject",true),ival(integer_input),dval(double_input) { };
  const std::set<std::string>& getSystemPropertyNames() const;
  std::vector<std::string> getPropertyNamesVector() const;
  ObjectHandler::property_t getSystemProperty(const std::string& name) const;
  void setSystemProperty(const std::string& name, const ObjectHandler::property_t& value);
};

class ESTestObject : public ObjectHandler::LibraryObject<ESTestVal> {
public:
  ESTestObject(const boost::shared_ptr<ObjectHandler::ValueObject>& properties,const long& integer_input,const double& double_input)
            : ObjectHandler::LibraryObject<ESTestVal>(properties,true) 
  {
    libraryObject_ = boost::shared_ptr<ESTestVal>(new ESTestVal(integer_input,double_input));
  }
  const double &double_value()  { // pass through to underlying LibraryObject
    return libraryObject_->double_value();
  }
  const long &integer_value()  { // pass through to underlying LibraryObject
    return libraryObject_->integer_value();
  }
};

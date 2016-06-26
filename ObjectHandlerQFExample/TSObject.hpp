#include <oh/libraryobject.hpp>
#include <oh/valueobject.hpp>
#include "TSBootstrap.hpp"

class ESTSValueObject : public ObjectHandler::ValueObject {
private:
  static const char* mPropertyNames[];
public:
  ESTSValueObject() {};
  ESTSValueObject(const std::string& objectID) : ObjectHandler::ValueObject(objectID,"ESTSObject",true) { };
  const std::set<std::string>& getSystemPropertyNames() const;
  std::vector<std::string> getPropertyNamesVector() const;
  ObjectHandler::property_t getSystemProperty(const std::string& name) const;
  void setSystemProperty(const std::string& name, const ObjectHandler::property_t& value);
};

class ESTSObject : public ObjectHandler::LibraryObject<quantfin::TSLogLinear> {
public:
  ESTSObject(const boost::shared_ptr<ObjectHandler::ValueObject>& properties,boost::shared_ptr<quantfin::TSLogLinear> xts)
            : ObjectHandler::LibraryObject<quantfin::TSLogLinear>(properties,true) 
  {
    libraryObject_ = xts;
  }
  double discount_factor(double t)  { // pass through to underlying LibraryObject
    return libraryObject_->operator()(t);
  }
};

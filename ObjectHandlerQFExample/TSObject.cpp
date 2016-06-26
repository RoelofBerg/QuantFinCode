#include <oh/exception.hpp>
#include <iostream>
#include <string>
#include <boost/algorithm/string/case_conv.hpp>
#include "TSObject.hpp"

const char* ESTSValueObject::mPropertyNames[] = {
        "Permanent"};

const std::set<std::string>& ESTSValueObject::getSystemPropertyNames() const {
        static std::set<std::string> ret;
        if(ret.empty())
            ret = std::set<std::string>(mPropertyNames,
                mPropertyNames + sizeof(mPropertyNames)/sizeof(const char*));
        return ret;
    }

std::vector<std::string> ESTSValueObject::getPropertyNamesVector() const {
        std::vector<std::string> ret(mPropertyNames,
            mPropertyNames + sizeof(mPropertyNames)/sizeof(const char*));
        for (std::map<std::string, ObjectHandler::property_t>::const_iterator i
            = userProperties.begin(); i != userProperties.end(); ++i)
            ret.push_back(i->first);
        return ret;
    }

ObjectHandler::property_t ESTSValueObject::getSystemProperty(const std::string& name) const {
        std::string nameUpper = boost::algorithm::to_upper_copy(name);
        if (strcmp(nameUpper.c_str(), "OBJECTID")==0)
            return objectId_;
        else if (strcmp(nameUpper.c_str(), "CLASSNAME")==0)
            return className_;
        else if (strcmp(nameUpper.c_str(), "PERMANENT")==0)
            //return (long)permanent_;
            return (bool)permanent_;
        else
            OH_FAIL("Error: attempt to retrieve non-existent Property: '" + name + "'");
    }

void ESTSValueObject::setSystemProperty(const std::string& name, const ObjectHandler::property_t& value) {
        std::string nameUpper = boost::algorithm::to_upper_copy(name);
        if (strcmp(nameUpper.c_str(), "OBJECTID")==0)
            objectId_ = boost::get<std::string>(value);
        else if (strcmp(nameUpper.c_str(), "CLASSNAME")==0)
            className_ = boost::get<std::string>(value);
        else if (strcmp(nameUpper.c_str(), "PERMANENT")==0)
            permanent_ = boost::get<bool>(value);
        else
            OH_FAIL("Error: attempt to retrieve non-existent Property: '" + name + "'");
    }

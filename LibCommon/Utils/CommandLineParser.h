//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//
//    .--------------------------------------------------.
//    |  This file is part of NTCodeBase                 |
//    |  Created 2018 by NT (https://ttnghia.github.io)  |
//    '--------------------------------------------------'
//                            \o/
//                             |
//                            / |
//
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

#pragma once

#include <LibCommon/CommonSetup.h>
#include <algorithm>

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
namespace NTCodeBase {
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
class CommandLineParser : public std::map<String, String> {
public:
    CommandLineParser() = default;
    CommandLineParser(int argc, char** argv) { init(argc, argv); }

    void init(int argc, char** argv) {
        this->insert(std::make_pair("Program", argv[0]));

        for(int arg = 1; arg < argc; ++arg) {
            String argstr(argv[arg]);
            if(argstr.find('=') == String::npos) {
                fprintf(stderr, "Invalid parameter input: variable assignment must be '-name=value' or 'name=value' without spacing\n");
                exit(EXIT_FAILURE);
            }

            const String tag(argstr.substr(0, argstr.find('=')));
            const String val(argstr.substr(tag.size() + 1));
            if(this->find(tag) == this->end()) {
                this->insert(std::make_pair(tag, val));
            } else {
                (*this)[tag] = val;
            }
        }
    }

    String getProgramName() { return (*this)[String("Program")]; }

    bool hasParam(const String& tag) const {
        return this->find(tag) != this->end();
    }

    bool getString(const String& tag, String& value) const {
        auto it = this->find(tag);
        if(it == this->end()) {
            return false;
        }
        value = it->second;
        return true;
    }

    String getStringValue(const String& tag) const {
        auto it = this->find(tag);
        __NT_REQUIRE(it != this->end());
        return it->second;
    }

    bool getInt(const String& tag, int& value) const {
        auto it = this->find(tag);
        if(it == this->end()) {
            return false;
        }
        value = std::stoi(it->second);
        return true;
    }

    int getIntValue(const String& tag) const {
        auto it = this->find(tag);
        __NT_REQUIRE(it != this->end());
        value = std::stoi(it->second);
    }

    template<class Real_t>
    bool getReal(const String& tag, Real_t& value) const {
        auto it = this->find(tag);
        if(it == this->end()) {
            return false;
        }
        value = static_cast<Real_t>(std::stof(it->second));
        return true;
    }

    template<class Real_t>
    Real_t getRealValue(const String& tag) const {
        auto it = this->find(tag);
        __NT_REQUIRE(it != this->end());
        return static_cast<Real_t>(std::stof(it->second));
    }
};
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
} // end namespace NTCodeBase

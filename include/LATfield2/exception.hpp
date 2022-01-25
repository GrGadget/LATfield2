#pragma once

#include <string>
#include <stdexcept>

namespace LATfield2
{
    class bad_dimensions : public std::runtime_error
    {
        public: 
        explicit bad_dimensions(const std::string& what_arg):
            std::runtime_error(what_arg)
        {}
        explicit bad_dimensions( const char* what_arg):
            std::runtime_error(what_arg)
        {}
    };
}

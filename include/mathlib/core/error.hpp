#pragma once
#include <stdexcept>
#include <string>

namespace mathlib::core {

    struct dimension_error : std::logic_error {
        explicit dimension_error(const std::string& msg)
            : std::logic_error("mathlib dimension error: " + msg) {}
    };

    struct domain_error : std::domain_error {
        explicit domain_error(const std::string& msg)
            : std::domain_error("mathlib domain error: " + msg) {}
    };

} // namespace mathlib::core
#pragma once

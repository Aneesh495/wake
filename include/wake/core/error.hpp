#pragma once

#include <stdexcept>
#include <string>

namespace wake {

class WakeException : public std::runtime_error {
public:
    explicit WakeException(const std::string& message)
        : std::runtime_error("wake: " + message) {}
};

class NumericalDivergenceException : public WakeException {
public:
    explicit NumericalDivergenceException(const std::string& details)
        : WakeException("Numerical divergence detected: " + details) {}
};

class InvalidConfigurationException : public WakeException {
public:
    explicit InvalidConfigurationException(const std::string& details)
        : WakeException("Invalid configuration: " + details) {}
};

} // namespace wake

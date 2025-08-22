#pragma once
#include <nlohmann/json.hpp>
#include <functional>
#include <vector>
#include <string>
#include <stdexcept>
#include <cmath>
#include <map>

class KinematicCorrections {
public:
    using CorrFunc = std::function<double(const std::vector<double>&, double)>;

    // Construct from JSON file path
    KinematicCorrections(const std::string& filename);
    // Construct from JSON object (can be empty → no corrections)
    KinematicCorrections(const nlohmann::json& j);

    // Apply corrections
    double deltaP(double p, double theta, bool isFD) const;
    double deltaTheta(double p, double theta, bool isFD) const;
    double deltaPhi(double p, double theta, bool isFD) const;

private:
    bool enabled_ = false;
    
    struct CorrEntry {
        CorrFunc func;
        nlohmann::json coeffsJson;
    };

    std::map<std::string, CorrEntry> corrections_;

    // Lambda factory
    CorrFunc makeProfileFunc(const std::string& form) const;

    // Evaluate polynomial in theta: coeff[0] + coeff[1]*theta + ...
    double evalPoly(const std::vector<double>& coeffs, double theta) const;

    // Convert JSON coefficients to numeric vector for this theta
    std::vector<double> getCoeffsFromTheta(const nlohmann::json& coeffsJson, double theta) const;
};

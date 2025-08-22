// SFCuts.h
#pragma once
#include "nlohmann/json.hpp"
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <unordered_map>

class SFCuts {
public:
    // Construct from JSON file path
    SFCuts(const std::string& filename);
    // Construct from JSON object
    SFCuts(const nlohmann::json& j);

    double mu_p(int sec, double p) const;
    double sigma_p(int sec, double p) const;
    bool pass(int sec, double sf, double p) const;

private:
    bool enabled_ = false;
    struct SectorCoeffs {
        std::vector<double> mu_coeffs;
        std::vector<double> sigma_coeffs;
    };
    std::unordered_map<int, SectorCoeffs> coeffsMap;

    double evalPoly(const std::vector<double>& coeffs, double p) const;
};

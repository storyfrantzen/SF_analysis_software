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
    struct SectorCoeffs {
        std::vector<double> mu_coeffs;
        std::vector<double> sigma_coeffs;
    };

    // Construct from JSON file path
    SFCuts(const std::string& filename);
    // Construct from JSON object
    SFCuts(const nlohmann::json& j);

    double mu_p(int sec, double p) const;
    double sigma_p(int sec, double p) const;
    bool passSigmaCut(int sec, double sf, double p, float numSigma = 3.5) const;
    bool passTriangleCut(double ePCAL, double eECIN, double p, float yScale = 1, float xScale = 1, 
                         float hypotenuse = 0.2, float HTCC_THRESHOLD = 4.5) const;

private:
    bool enabled_ = false;
    std::unordered_map<int, SectorCoeffs> coeffsMap;
    double evalPoly(const std::vector<double>& coeffs, double p) const;
};

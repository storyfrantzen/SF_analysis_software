#include "SFCuts.h"
#include <fstream>
#include <cmath>
#include <stdexcept>

using json = nlohmann::json;

// --- Construct from JSON file ---
SFCuts::SFCuts(const std::string& filename) {
    std::ifstream f(filename);
    if (!f.is_open()) {
        enabled_ = false;  // file missing → treat as no-op
        return;
    }
    json j;
    f >> j;
    *this = SFCuts(j);  // delegate to JSON constructor
}

// --- Construct from JSON object ---
SFCuts::SFCuts(const json& j) {
    if (j.is_null() || j.empty()) {
        enabled_ = false;
        return;
    }

    for (int sec = 1; sec <= 6; ++sec) {
        std::string key = "sector_" + std::to_string(sec);
        if (!j.contains(key)) continue;

        SectorCoeffs sc;
        sc.mu_coeffs = j[key]["mu_coeffs"].get<std::vector<double>>();
        sc.sigma_coeffs = j[key]["sigma_coeffs"].get<std::vector<double>>();
        coeffsMap[sec] = sc;
    }

    if (coeffsMap.empty()) enabled_ = false;
    else enabled_ = true;
}

// --- Evaluate polynomial ---
double SFCuts::evalPoly(const std::vector<double>& coeffs, double p) const {
    double val = 0.0;
    for (size_t i = 0; i < coeffs.size(); i++) {
        val += coeffs[i] * std::pow(p, coeffs.size() - 1 - i);
    }
    return val;
}

double SFCuts::mu_p(int sec, double p) const {
    return evalPoly(coeffsMap.at(sec).mu_coeffs, p);
}

double SFCuts::sigma_p(int sec, double p) const {
    return evalPoly(coeffsMap.at(sec).sigma_coeffs, p);
}

// Note: currently, triangle cut is enabled IFF sigma cut is enabled.
bool SFCuts::passTriangleCut(double ePCAL, double eECIN, double p, float yScale, float xScale, float hypotenuse, float HTCC_THRESHOLD) const {
    if (!enabled_) return true;
    if (p < HTCC_THRESHOLD) return true;
    double y = ePCAL / p;
    double x = eECIN / p;
    return (y * yScale + x * xScale) > hypotenuse;
}

bool SFCuts::passSigmaCut(int sec, double sf, double p, float numSigma) const {
    if (!enabled_) return true;
    double mu = mu_p(sec, p);
    double sigma = sigma_p(sec, p);
    return (sf > mu - numSigma * sigma && sf < mu + numSigma * sigma);
}

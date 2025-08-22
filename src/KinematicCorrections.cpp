#include <fstream>
#include <iostream>
#include "KinematicCorrections.h"
#include "nlohmann/json.hpp"

// Construct from JSON file
KinematicCorrections::KinematicCorrections(const std::string& filename) {
    std::ifstream f(filename);
    if (!f.is_open()) {
        std::cerr << "KinematicCorrections file not found: " << filename 
                  << " → corrections disabled.\n";
        enabled_ = false;
        return;
    }
    nlohmann::json j;
    f >> j;
    *this = KinematicCorrections(j);
}

// Construct from JSON object
KinematicCorrections::KinematicCorrections(const nlohmann::json& j) {
    if (j.is_null() || j.empty()) {
        enabled_ = false;
        return;
    }
    enabled_ = true;
    for (auto it = j.begin(); it != j.end(); ++it) {
        CorrEntry entry;
        entry.func = makeProfileFunc(it.value()["form"].get<std::string>());
        entry.coeffsJson = it.value()["coeffs"];
        corrections_[it.key()] = entry;
    }
}

KinematicCorrections::CorrFunc KinematicCorrections::makeProfileFunc(const std::string& form) const {
    if (form == "[0] + [1]/p + [2]/(p^2)")
        return [](const std::vector<double>& c, double p) { return c[0] + c[1]/p + c[2]/(p*p); };
    if (form == "[0] + [1]/p")
        return [](const std::vector<double>& c, double p) { return c[0] + c[1]/p; };
    if (form == "[0] + [1]*p + [2]*p^2")
        return [](const std::vector<double>& c, double p) { return c[0] + c[1]*p + c[2]*p*p; };
    if (form == "[0] + [1]*p")
        return [](const std::vector<double>& c, double p) { return c[0] + c[1]*p; };
    throw std::runtime_error("Unsupported correction form: " + form);
}

double KinematicCorrections::evalPoly(const std::vector<double>& coeffs, double theta) const {
    double val = 0.0, powTheta = 1.0;
    for (double c : coeffs) {
        val += c * powTheta;
        powTheta *= theta;
    }
    return val;
}

std::vector<double> KinematicCorrections::getCoeffsFromTheta(const nlohmann::json& coeffsJson, double theta) const {
    std::vector<double> coeffs;
    for (auto it = coeffsJson.begin(); it != coeffsJson.end(); ++it) {
        coeffs.push_back(evalPoly(it.value().get<std::vector<double>>(), theta));
    }
    return coeffs;
}

double KinematicCorrections::deltaP(double p, double theta_rad, bool isFD) const {
    if (!enabled_) return 0.0;  // no correction
    const std::string key = isFD ? "p_delta_p_FD" : "p_delta_p_CD";
    auto it = corrections_.find(key);
    if (it == corrections_.end()) return 0.0;  // No correction
    double theta_deg = theta_rad * 180.0 / M_PI;
    std::vector<double> coeffs = getCoeffsFromTheta(it->second.coeffsJson, theta_deg);
    return it->second.func(coeffs, p);; 
}

double KinematicCorrections::deltaTheta(double p, double theta_rad, bool isFD) const {
    if (!enabled_) return 0.0;  // no correction
    const std::string key = isFD ? "p_delta_theta_FD" : "p_delta_theta_CD";
    auto it = corrections_.find(key);
    if (it == corrections_.end()) return 0.0;  
    double theta_deg = theta_rad * 180.0 / M_PI;
    std::vector<double> coeffs = getCoeffsFromTheta(it->second.coeffsJson, theta_deg);
    double delta_theta_deg = it->second.func(coeffs, p);
    return delta_theta_deg * M_PI / 180.0;
}

double KinematicCorrections::deltaPhi(double p, double theta_rad, bool isFD) const {
    if (!enabled_) return 0.0;  // no correction
    const std::string key = isFD ? "p_delta_phi_FD" : "p_delta_phi_CD";
    auto it = corrections_.find(key);
    if (it == corrections_.end()) return 0.0;  
    double theta_deg = theta_rad * 180.0 / M_PI;
    std::vector<double> coeffs = getCoeffsFromTheta(it->second.coeffsJson, theta_deg);
    double delta_phi_deg = it->second.func(coeffs, p);
    return delta_phi_deg * M_PI / 180.0;
}

#pragma once
#include <iostream>     // for std::cout, std::cerr
#include <iomanip>      // for std::put_time
#include <ctime>        // for std::time, std::localtime
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <algorithm>

class ConfigLoader {
public:
    ConfigLoader(const std::string& filename) {
        std::ifstream infile(filename);
        std::string line;
        while (std::getline(infile, line)) {
            // Trim leading/trailing whitespace first
            trim(line);

            // Explicitly skip full-line comments or empty lines
            if (line.empty() || line[0] == '#') continue;
            
            size_t comment_pos = line.find('#');
            if (comment_pos != std::string::npos)
                line = line.substr(0, comment_pos);

            // Trim whitespace
            trim(line);
            if (line.empty()) continue;

            size_t eq_pos = line.find('=');
            if (eq_pos != std::string::npos) {
                std::string key = line.substr(0, eq_pos);
                std::string value = line.substr(eq_pos + 1);
                trim(key);
                trim(value);
                params_[key] = value;
            }

        }
    }

    // Generic getter functions defined here for simplicity. Get value as string, double, int, bool
    std::string get(const std::string& key) const { return params_.at(key); }
    double getDouble(const std::string& key) const { return std::stod(params_.at(key));}
    int getInt(const std::string& key) const { return std::stoi(params_.at(key)); }
    bool getBool(const std::string& key) const {
        std::string v = params_.at(key);
        std::transform(v.begin(), v.end(), v.begin(), ::tolower);
        return (v == "true" || v == "1");
    }
    bool contains(const std::string& key) const { return params_.find(key) != params_.end(); }

    // Helper function which outputs .txt file containing a list of key:val config tags parsed successfully
    void writeSummary() const {
        auto now = std::time(nullptr);
        std::ostringstream output_fs;
        output_fs << "config_summary_"
                << std::put_time(std::localtime(&now), "%Y%m%d_%H%M%S")
                << ".txt";

        std::ofstream out(output_fs.str());
        if (!out) {
            std::cerr << "Failed to open config summary file for writing.\n";
            return;
        }

        out << "# Summary of configuration tags parsed by ConfigLoader :\n";
        for (const auto& [key, value] : params_) {
            out << key << " = " << value << "\n";
        }

        std::cout << "Config summary written to " << output_fs.str() << std::endl;
    }

    // Specialized helper functions defined in ConfigLoader.cpp:

    std::vector<std::string> getQARequirements() const;
    std::vector<std::string> getFiducialCuts() const;
    std::vector<std::string> getTopology() const;

    std::string getOutputFile() const;
    std::string getChannel() const;
    
    int getTorus() const;
    double getEbeam() const;


private:
    std::unordered_map<std::string, std::string> params_;

    static void trim(std::string& s) {
        size_t start = s.find_first_not_of(" \t");
        size_t end   = s.find_last_not_of(" \t");
        if (start == std::string::npos || end == std::string::npos)
            s = "";
        else
            s = s.substr(start, end - start + 1);
    }
};

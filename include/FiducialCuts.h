// FiducialCuts.h
#ifndef FIDUCIALCUTS_H
#define FIDUCIALCUTS_H

#include <array>

#include <clas12reader.h>
#include <region_particle.h>

class FiducialCuts {
public:

    FiducialCuts();

    void addCut(const std::string& tag);

    bool inFTHole(float x, float y) const;
    bool inExcludedECALRegion(int detector, int sector, char coord, double value) const;
    std::pair<double, double> rotateToSector1Frame(double x, double y, int sector);

    bool passesDC(clas12::region_particle* p, const int& torus);
    bool passesFT(clas12::region_particle* p);
    bool passesECAL(clas12::region_particle* p);
    bool passesCVT(clas12::region_particle* p);

    struct Hole { float x, y, radius; };

private:

    std::vector<Hole> ftHoles_; // List of standard FT holes

    // Key: detector_sector_coordinate (e.g., "PCal_1_w"), value: list of [min, max] intervals to exclude
    std::map<std::string, std::vector<std::pair<float, float>>> rgaExclusionMap_;
    std::map<std::string, std::vector<std::pair<float, float>>> rgkExclusionMap_;

    // Each row = one parameter (Psplit, Tleft, Tright, ...)
    // Each column = sector 1 to 6 (index 0 to 5)
    std::array<std::array<double, 6>, 9> rgkParams_;

    // Fiducial flags here:
    bool DCedgeCut_      = false;
    bool FTstandardCut_  = false;
    bool ECALrgaCut_     = false;
    bool ECALrgkCut_     = false;
    bool ECALrgkEdgeCut_ = false;
    bool CVTrgaCut_      = false;
};

#endif

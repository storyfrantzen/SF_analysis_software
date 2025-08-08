#include "FiducialCuts.h"

#include <cmath>  // For cos, sin, M_PI
#include <algorithm>  // for std::transform

using namespace clas12;

FiducialCuts::FiducialCuts() {

    // Format: (hole_x, hole_y, hole_radius)
    ftHoles_ = {
        {-8.42f, 9.89f, 1.60f},
        {-9.89f, -5.33f, 1.60f},
        {-6.15f, -13.0f, 2.30f},
        {3.70f, -6.50f, 2.00f}
    };

    // Note: map below does NOT include SP19 PCAL S2 v cuts, unless cut == ECALrgaS19Cut
    rgaExclusionMap_ = {
        {"1_1_w",  {{72.0, 94.5}, {220.5, 234.0}}},
        {"1_2_v",  {{99.0, 117.0}}},
        {"1_3_w",  {{346.5, 378.0}}},
        {"1_4_w",  {{0.0, 13.5}}},
        {"1_4_v",  {{229.5, 243.0}}},
        {"1_6_w",  {{166.5, 193.5}}},
        {"4_1_v",  {{67.5, 94.5}}},
        {"4_5_v",  {{0.0, 23.5}}},
        {"7_1_v",  {{0.0, 40.5}}},
        {"7_5_u",  {{193.5, 216.0}}}
    };

    // Note: map below is done mostly by eye and not rigorously. Want to do this systematically if possible.
    rgkExclusionMap_ = {
        {"1_1_w", {{72, 94.5}, {207, 234}}},
        //{"1_1_u", {{}}},
        {"1_2_v", {{99, 112.5}}},
        {"1_3_v", {{346.5, 382.5}}},
        {"1_4_v", {{229.5, 243}}},
        {"1_6_w", {{171, 193.5}}},
        {"2_1_v", {{67.5, 99}}},
        
    };

    rgkEdgeParams_ = {{
        { 87,        82,        85,        77,        78,        82},       // Psplit
        { 58.7356,   62.8204,   62.2296,   53.7756,   58.2888,   54.5822},  // Tleft
        { 58.7477,   51.2589,   59.2357,   56.2415,   60.8219,   49.8914},  // Tright
        { 0.582053,  0.544976,  0.549788,  0.56899,   0.56414,   0.57343},  // Sleft
        {-0.591876, -0.562926, -0.562246, -0.563726, -0.568902, -0.550729}, // Sright
        { 64.9348,   64.7541,   67.832,    55.9324,   55.9225,   60.0997},  // Rleft
        { 65.424,    54.6992,   63.6628,   57.8931,   56.5367,   56.4641},  // Rright
        { 0.745578,  0.606081,  0.729202,  0.627239,  0.503674,  0.717899}, // Qleft
        {-0.775022, -0.633863, -0.678901, -0.612458, -0.455319, -0.692481}  // Qright
    }};

}

void FiducialCuts::addCut(const std::string& tag) {
    std::string cut = tag;

    if (cut == "DCEdges_RGA") DCEdges_RGA_ = true;

    else if (cut == "FT_RGA") FT_RGA_ = true;

    else if (cut == "ECAL_RGA") ECAL_RGA_ = true; 

    else if (cut == "ECAL_RGAS19") {
        rgaExclusionMap_["1_2_v"].emplace_back(31.5, 49.5);
        ECAL_RGA_ = true;
    }

    else if (cut == "ECAL_RGK") ECAL_RGK_ = true;

    else if (cut == "ECALEdges_RGK") ECALEdges_RGK_ = true;

    else if (cut == "CVT_RGA") CVT_RGA_ = true;

    else {
        std::cerr << "[FiducialCuts] Warning: Unrecognized cut tag '" << tag << ", " << cut << "'\n";
    }
}

bool FiducialCuts::inFTHole(float x, float y) const {
    for (const auto& hole : ftHoles_) {
        float dx = x - hole.x;
        float dy = y - hole.y;
        float dist = std::sqrt(dx * dx + dy * dy);
        if (dist < hole.radius) return true;
    }
    return false;
}

// Note: currently compatible with only 1 exclusion map at a time
bool FiducialCuts::inExcludedECALRegion(int detector, int sector, char coord, double value) const {
    std::string key = std::to_string(detector) + "_" + std::to_string(sector) + "_" + coord;
    std::map<std::string, std::vector<std::pair<float, float>>> exclusionMap;
    if (ECAL_RGA_) exclusionMap = rgaExclusionMap_;
    else if (ECAL_RGK_) exclusionMap = rgkExclusionMap_;
    auto it = exclusionMap.find(key);
    if (it != exclusionMap.end()) {
        for (const auto& range : it->second) {
            if (value >= range.first && value <= range.second) {
                //std::cout << "CUT ! key = " << key << ", and value = " << value << std::endl;
                return true;
            }
        }
    }
    return false;
}

// Rotate point (x, y) by -60 * (sector - 1) degrees to align with sector 1
std::pair<double, double> FiducialCuts::rotateToSector1Frame(double x, double y, int sector) {
    // sector 1 needs no rotation
    if (sector == 1) return {x, y};

    double angle_deg = -60.0 * (sector - 1); // negative for counter-clockwise to clockwise
    double angle_rad = angle_deg * M_PI / 180.0;

    double x_rot =  x * cos(angle_rad) - y * sin(angle_rad);
    double y_rot =  x * sin(angle_rad) + y * cos(angle_rad);

    return {x_rot, y_rot};
}

bool FiducialCuts::passesDC(clas12::region_particle* p, const int& torus) {

    if (!DCEdges_RGA_ ) return true;

    int charge = p->par()->getCharge();
    // Might need to check, but this ensures particle-independent cut:
    bool outbending = torus * charge < 0 ? true : false;

    int edge1  = p->traj(DC, 6)  ? p->traj(DC, 6)->getEdge()  : 0;
    int edge2  = p->traj(DC, 18) ? p->traj(DC, 18)->getEdge() : 0;
    int edge3  = p->traj(DC, 36) ? p->traj(DC, 36)->getEdge() : 0;

    int pid    = p->par()->getPid();

    // ELECTRON EDGE VALUES:
    if (std::abs(pid) == 11) {
        if (edge3 <= 10) return false;
        if (outbending)  return edge1 > 3 && edge2 > 3;
        else             return p->getTheta() > 10 ? (edge1 > 3 && edge2 > 3) : (edge1 > 10 && edge2 > 10);
    }

    // PROTON EDGE VALUES:
    else if (std::abs(pid) == 2212) {
        if (edge3 <= 5)  return false;
        if (!outbending) return edge1 > 3 && edge2 > 3;
    }

    return true;
}

bool FiducialCuts::passesFT(clas12::region_particle* p) {
    if (!FT_RGA_) return true;

    double x = p->ft(FTCAL)->getX();
    double y = p->ft(FTCAL)->getY();
    double r = sqrt(x * x + y * y);
    // Reject if outside min and max radii:
    if (r < 8.5 || r > 15.5) return false;

    // Remove circles of low efficiency, i.e., if |r - r_hole| > r_allowed:
    // Note for later: can add config tag for generic hole centered at (x, y) with radius r
    return (!inFTHole(x, y));
}

bool FiducialCuts::passesECAL(clas12::region_particle* p) {
    if (ECAL_RGA_ || ECAL_RGK_) {

        int sector = p->cal(1) ? p->cal(1)->getSector() : -1;
        if (sector < 0) return false;

        int layers[] = {1, 4, 7};
        double u[3], v[3], w[3];

        for (int i = 0; i < 3; ++i) {
            // NOTE: not sure if the line below rejects possible events 
            if (!p->cal(layers[i])) return false;

            u[i] = p->cal(layers[i])->getLu();
            v[i] = p->cal(layers[i])->getLv();
            w[i] = p->cal(layers[i])->getLw();

            // FOR NOW, HARDCODE LOOSE THRESHOLD > 9 , MEDIUM > 13.5 , TIGHT > 18: 
            if (v[0] < 9.0 || w[0] < 9.0) return false; // NEED TO CHECK IF THIS IS WHAT THE SLIDES MEAN!

            if (inExcludedECALRegion(layers[i], sector, 'u', u[i])) { 
                //std::cout << "Excluded: " << "Layer: " << layers[i] << ", Sector: " << sector << ", " << "u: " << u[i] << std::endl;
                return false;
            } 
            if (inExcludedECALRegion(layers[i], sector, 'v', v[i])) {
                //std::cout << "Excluded: " << "Layer: " << layers[i] << ", Sector: " << sector << ", " << "v: " << v[i] << std::endl;
                return false;
            }
            if (inExcludedECALRegion(layers[i], sector, 'w', w[i])) {
                //std::cout << "Excluded: " << "Layer: " << layers[i] << ", Sector: " << sector << ", " << "w: " << w[i] << std::endl;
                return false;
            };
        }
    }
    if (ECALEdges_RGK_) {

        int sector = p->cal(1) ? p->cal(1)->getSector() : 999;
        if (sector < 0) return false;

        double x = p->cal(1)->getX();
        double y = p->cal(1)->getY();

        auto [x_local, y_local] = rotateToSector1Frame(x, y, sector);

        int    pSplit = rgkEdgeParams_[0][sector - 1];
        double tLeft  = rgkEdgeParams_[1][sector - 1];
        double tRight = rgkEdgeParams_[2][sector - 1];
        double sLeft  = rgkEdgeParams_[3][sector - 1];
        double sRight = rgkEdgeParams_[4][sector - 1];
        double rLeft  = rgkEdgeParams_[5][sector - 1];
        double rRight = rgkEdgeParams_[6][sector - 1];
        double qLeft  = rgkEdgeParams_[7][sector - 1];
        double qRight = rgkEdgeParams_[8][sector - 1];

        bool condition1 = x_local > pSplit && (y_local < sLeft * (x_local - tLeft)) && (y_local > sRight * (x_local - tRight));
        bool condition2 = x_local < pSplit && (y_local < qLeft * (x_local - rLeft)) && (y_local > qRight * (x_local - rRight));

        return (condition1 || condition2);
    } 
    return true;
}

bool FiducialCuts::passesCVT(clas12::region_particle* p) {
    //NEED TO CHECK...
    if (CVT_RGA_) {
        // Define CVT layer IDs to check
        std::vector<int> layers = {1, 3, 5, 7, 12};
        for (int layer : layers) {
            double edge = NAN;
            if (p->traj(CVT, layer)) edge = p->traj(CVT, layer)->getEdge();
            if (edge <= 0 || std::isnan(edge)) return false;
        }

        // Get phi in degrees
        double x_cvt = p->traj(CVT, 1) ? p->traj(CVT, 1)->getX() : NAN;
        double y_cvt = p->traj(CVT, 1) ? p->traj(CVT, 1)->getY() : NAN;
        double z_cvt = p->traj(CVT, 1) ? p->traj(CVT, 1)->getZ() : NAN;

        double r_cvt = std::sqrt(x_cvt * x_cvt + y_cvt * y_cvt + z_cvt * z_cvt);
        double phi = std::atan2(y_cvt, x_cvt) * 180.0 / M_PI;

        if (std::isnan(phi)) return false;

        if (phi < 0) phi += 360;
        // RGA phi dead zones:
        if ((25 <= phi && phi <= 40) || (143 <= phi && phi <= 158) || (265 <= phi && phi <= 280)) return false;
    }
    return true; // CVT_RGA_ not enabled, automatically pass
}

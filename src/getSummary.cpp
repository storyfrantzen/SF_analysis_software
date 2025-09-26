
#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include <TSystem.h>

int main(int argc, char** argv) {
    // Load compiled dictionary. Note: NOT necessary, just suppresses warning messages.
    gSystem->Load("/work/clas12/storyf/SF_analysis_software/build/install/lib/libBranchVarsDict.so");

    if (argc < 2) {
        std::cerr << "Usage: getSummary <rootfile>" << std::endl;
        return 1;
    }
    TFile f(argv[1], "READ");
    if (f.IsZombie()) return 1;

    TTree* t = (TTree*)f.Get("Summary");
    if (!t) {
        std::cerr << "No Summary tree in file." << std::endl;
        return 1;
    }

    double charge; int events, fills;
    t->SetBranchAddress("TotalCharge", &charge);
    t->SetBranchAddress("EventsProcessed", &events);
    t->SetBranchAddress("Fills", &fills);
    t->GetEntry(0);

    std::cout << "Events processed: " << events
              << "\nAccumulated beam charge: " << charge << " nC" 
              << "\nFills made: " << fills << std::endl;
    return 0;
}

#ifndef PHYSICALCONSTANTS_H
#define PHYSICALCONSTANTS_H

#include <TDatabasePDG.h>

inline const TDatabasePDG* db = TDatabasePDG::Instance();

inline const double ELECTRON_MASS = db->GetParticle(11)->Mass();    // GeV
inline const double PROTON_MASS   = db->GetParticle(2212)->Mass();  // GeV
inline const double PI0_MASS      = db->GetParticle(111)->Mass();   // GeV

#endif

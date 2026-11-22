// SPDX-License-Identifier: LGPL-3.0-or-later
// SPDX-FileCopyrightText: Copyright CERN for the benefit of the SHiP
// Collaboration

#include "UpstreamTaggerScoringPlaneHit.h"

#include <iostream>

#include "UpstreamTaggerScoringPoint.h"
#include "FairLogger.h"
#include "TRandom.h"

UpstreamTaggerScoringPlaneHit::UpstreamTaggerScoringPlaneHit() : SHiP::DetectorHit() {}

UpstreamTaggerScoringPlaneHit::UpstreamTaggerScoringPlaneHit(const UpstreamTaggerScoringPoint* p,
                                         Double_t t0, Double_t pos_res,
                                         Double_t time_res,
                                         TRandom* rng)
    : SHiP::DetectorHit() {
  if (!p) {
    LOG(error) << "UpstreamTaggerScoringPlaneHit: null UpstreamTaggerScoringPoint pointer";
    return;
  }
  SetDetectorID(p->GetDetectorID());
  TRandom& random = rng ? *rng : *gRandom;
  fX = pos_res > 0. ? random.Gaus(p->GetX(), pos_res) : p->GetX();
  fY = pos_res > 0. ? random.Gaus(p->GetY(), pos_res) : p->GetY();
  fZ = p->GetZ();
  const Double_t t = p->GetTime() + t0;
  fTime = time_res > 0. ? random.Gaus(t, time_res) : t;
  fPx = p->GetPx();
  fPy = p->GetPy();
  fPz = p->GetPz();
}

void UpstreamTaggerScoringPlaneHit::Print() const {
  std::cout << "-I- UpstreamTaggerScoringPlaneHit: detector " << fDetectorID << std::endl;
  std::cout << "    Position: (" << fX << ", " << fY << ", " << fZ << ") cm"
            << "  time " << fTime << " ns" << std::endl;
  std::cout << "    Momentum: (" << fPx << ", " << fPy << ", " << fPz << ") GeV/c"
            << std::endl;
}

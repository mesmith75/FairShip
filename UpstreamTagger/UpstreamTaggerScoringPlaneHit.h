// SPDX-License-Identifier: LGPL-3.0-or-later
// SPDX-FileCopyrightText: Copyright CERN for the benefit of the SHiP
// Collaboration

#ifndef UPSTREAMTAGGERSCORINGPLANE_UPSTREAMTAGGERSCORINGPLANEHIT_H_
#define UPSTREAMTAGGERSCORINGPLANE_UPSTREAMTAGGERSCORINGPLANEHIT_H_

#include "DetectorHit.h"
#include "TVector3.h"

class UpstreamTaggerScoringPoint;
class TRandom;

/**
 * @brief Digitised UBT scoring-plane hit.
 *
 * Position and time are smeared from the MC point; the momentum is the
 * MC-truth momentum at the plane, and the detector ID is that of the MC
 * point. The MC truth link is kept in the separate
 * digiUpstreamTaggerScoringPlane2MC branch.
 */
class UpstreamTaggerScoringPlaneHit : public SHiP::DetectorHit {
 public:
  UpstreamTaggerScoringPlaneHit();

  /**
   * @param p        MC point
   * @param t0       event time offset (ns)
   * @param pos_res  Gaussian x/y resolution (cm), 0 = no smearing
   * @param time_res Gaussian time resolution (ns), 0 = no smearing
   * @param rng      generator for the smearing (gRandom if null); nothing is
   *                 drawn without smearing
   */
  UpstreamTaggerScoringPlaneHit(const UpstreamTaggerScoringPoint* p, Double_t t0,
                      Double_t pos_res, Double_t time_res,
                      TRandom* rng = nullptr);

  ~UpstreamTaggerScoringPlaneHit() override = default;
  UpstreamTaggerScoringPlaneHit(const UpstreamTaggerScoringPlaneHit& hit) = default;
  UpstreamTaggerScoringPlaneHit& operator=(const UpstreamTaggerScoringPlaneHit& hit) = default;

  Double_t GetX() const { return fX; }
  Double_t GetY() const { return fY; }
  Double_t GetZ() const { return fZ; }
  TVector3 GetXYZ() const { return TVector3(fX, fY, fZ); }
  Double_t GetTime() const { return fTime; }

  Double_t GetPx() const { return fPx; }
  Double_t GetPy() const { return fPy; }
  Double_t GetPz() const { return fPz; }
  TVector3 GetMom() const { return TVector3(fPx, fPy, fPz); }
  Double_t GetP() const { return GetMom().Mag(); }


  using SHiP::DetectorHit::Print;
  void Print() const;

 private:
  Double_t fX = 0.;     ///< smeared x (cm)
  Double_t fY = 0.;     ///< smeared y (cm)
  Double_t fZ = 0.;     ///< z of the MC point (cm)
  Double_t fTime = 0.;  ///< smeared time incl. t0 (ns)
  Double_t fPx = 0.;    ///< MC-truth momentum (GeV/c)
  Double_t fPy = 0.;
  Double_t fPz = 0.;

  ClassDefOverride(UpstreamTaggerScoringPlaneHit, 1);
};

#endif  // UPSTREAMTAGGERSCORINGPLANE_UPSTREAMTAGGERSCORINGPLANEHIT_H_

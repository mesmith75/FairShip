# SPDX-License-Identifier: LGPL-3.0-or-later
# SPDX-FileCopyrightText: Copyright CERN for the benefit of the SHiP Collaboration

"""Digitisation of the UBT scoring plane.

Each accepted UpstreamTaggerScoringPoint becomes one UpstreamTaggerScoringPlaneHit
(UpstreamTagger/UpstreamTaggerScoringPlaneHit.h), with Gaussian-smeared x, y and
time and the MC-truth momentum. The resolutions, the forward-only cut and the
ignored PDG codes come from ShipGeo.UpstreamTaggerScoringPlane.
"""

import logging

import global_variables
import ROOT
from BaseDetector import BaseDetector

logger = logging.getLogger(__name__)


class UpstreamTaggerScoringPlaneDetector(BaseDetector):
    """One UpstreamTaggerScoringPlaneHit per accepted UpstreamTaggerScoringPoint.

    Output branches: ``Digi_UpstreamTaggerScoringPlaneHits`` and ``digiUpstreamTaggerScoringPlane2MC``
    (per hit, the list with the MC track index).
    """

    def __init__(self, name, intree, outtree=None) -> None:
        super().__init__(name, intree, mcBranchName="digiUpstreamTaggerScoringPlane2MC", outtree=outtree)
        cfg = global_variables.ShipGeo.UpstreamTaggerScoringPlane
        self.pos_res = cfg.get("PositionResolution", 0.0)
        self.time_res = cfg.get("TimeResolution", 0.0)
        self.forward_only = cfg.get("ForwardOnly", True)
        self.ignored = {abs(int(p)) for p in cfg.get("IgnoredPdg", [12, 14, 16])}
        # Own generator, seeded like shipVeto.py's: drawing the smearing from it
        # leaves the random sequences of the other detectors unchanged.
        self.random = ROOT.TRandom3()
        self.random.SetSeed(13)

    def digitize(self) -> None:
        t0 = self.intree.t0
        for point in self.intree.UpstreamTaggerScoringPoint:
            if abs(point.PdgCode()) in self.ignored:
                continue
            if self.forward_only and point.GetPz() <= 0:
                continue
            self.det.push_back(ROOT.UpstreamTaggerScoringPlaneHit(point, t0, self.pos_res, self.time_res, self.random))
            if self.MCdet is not None:  # always set (mcBranchName given); narrows the type
                links = ROOT.std.vector("int")()
                links.push_back(point.GetTrackID())
                self.MCdet.push_back(links)

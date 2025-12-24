// SPDX-License-Identifier: LGPL-3.0-or-later
// SPDX-FileCopyrightText: Copyright CERN for the benefit of the SHiP Collaboration

// RPC Timing Detector
// 17/12/2019
// celso.franco@cern.ch

#include "UpstreamTagger.h"
#include "UpstreamTaggerPoint.h"
#include "UpstreamTaggerHit.h"

#include "FairVolume.h"
#include "FairGeoVolume.h"
#include "FairGeoNode.h"
#include "FairRootManager.h"
#include "FairGeoLoader.h"
#include "FairGeoInterface.h"
#include "FairGeoMedia.h"
#include "FairGeoBuilder.h"
#include "FairRun.h"
#include "FairRuntimeDb.h"
#include "ShipDetectorList.h"
#include "ShipStack.h"

#include "TClonesArray.h"
#include "TVirtualMC.h"
#include "TGeoManager.h"
#include "TGeoBBox.h"
#include "TGeoCompositeShape.h"
#include "TGeoTube.h"
#include "TGeoMaterial.h"
#include "TGeoMedium.h"
#include "TMath.h"
#include "TParticle.h"
#include "TVector3.h"

#include <ROOT/TSeq.hxx>
#include <iostream>
#include <sstream>
using std::cout;
using std::endl;
using ROOT::TSeq;
using ShipUnit::m;
using ShipUnit::cm;


UpstreamTagger::UpstreamTagger()
  : FairDetector("UpstreamTagger", kTRUE, kUpstreamTagger),
    fTrackID(-1),
    fVolumeID(-1),
    fPos(),
    fMom(),
    fTime(-1.),
    fLength(-1.),
    fELoss(-1),
    //
    det_zPos(0),

    //
    UpstreamTagger_fulldet(0),
    scoringPlaneUBText(0),
  //
    fUpstreamTaggerPointCollection(new TClonesArray("UpstreamTaggerPoint"))
{
}

//UpstreamTagger::UpstreamTagger(const char* name, Bool_t active)
UpstreamTagger::UpstreamTagger(std::string medium)
  : FairDetector("UpstreamTagger", kTRUE, kUpstreamTagger),
    fMedium(medium),
    fTrackID(-1),
    fVolumeID(-1),
    fPos(),
    fMom(),
    fTime(-1.),
    fLength(-1.),
    fELoss(-1),
    //
    det_zPos(0),

    //
    UpstreamTagger_fulldet(0),
    scoringPlaneUBText(0), // Initialize new scoring plane to nullptr
        //
    fUpstreamTaggerPointCollection(new TClonesArray("UpstreamTaggerPoint"))
{
}


void UpstreamTagger::Initialize()
{
  FairDetector::Initialize();
}


UpstreamTagger::~UpstreamTagger()
{
    std::cout<<"destroyer droid"<<std::endl;
  if (fUpstreamTaggerPointCollection) {
    fUpstreamTaggerPointCollection->Delete();
    delete fUpstreamTaggerPointCollection;
  }
}



Int_t UpstreamTagger::InitMedium(const char* name)
{

   static FairGeoLoader *geoLoad=FairGeoLoader::Instance();
   static FairGeoInterface *geoFace=geoLoad->getGeoInterface();
   static FairGeoMedia *media=geoFace->getMedia();
   static FairGeoBuilder *geoBuild=geoLoad->getGeoBuilder();

   FairGeoMedium *ShipMedium=media->getMedium(name);

   if (!ShipMedium)
   {
     Fatal("InitMedium","Material %s not defined in media file.", name);
     return -1111;
   }
   TGeoMedium* medium=gGeoManager->GetMedium(name);
   if (medium!=NULL)
     return ShipMedium->getMediumIndex();

   return geoBuild->createMedium(ShipMedium);

  return 0;
}



Bool_t  UpstreamTagger::ProcessHits(FairVolume* vol)
{
  std::cout<<"processing a hit"<<std::endl;
  /** This method is called from the MC stepping */
  //Set parameters at entrance of volume. Reset ELoss.
  if ( gMC->IsTrackEntering() ) {
    fELoss  = 0.;
    fTime   = gMC->TrackTime() * 1.0e09;
    fLength = gMC->TrackLength();
    gMC->TrackPosition(fPos);
    gMC->TrackMomentum(fMom);
  }
  std::cout<<"gello"<<std::endl;
  // Sum energy loss for all steps in the active volume
  fELoss += gMC->Edep();
  std::cout<<"eLoss: "<<fELoss<<std::endl;
  fELoss += 0.01;

  // Create vetoPoint at exit of active volume
  if ( gMC->IsTrackExiting()    ||
       gMC->IsTrackStop()       ||
       gMC->IsTrackDisappeared()   ) {
    if (fELoss == 0. ) { return kFALSE; }
    std::cout<<"going well"<<std::endl;
    fTrackID  = gMC->GetStack()->GetCurrentTrackNumber();

    Int_t uniqueId;
    gMC->CurrentVolID(uniqueId);
    if (uniqueId>1000000) //Solid scintillator case
    {
      Int_t vcpy;
      gMC->CurrentVolOffID(1, vcpy);
      if (vcpy==5) uniqueId+=4; //Copy of half
    }

    TParticle* p = gMC->GetStack()->GetCurrentTrack();
    Int_t pdgCode = p->GetPdgCode();
    TLorentzVector Pos;
    gMC->TrackPosition(Pos);
    TLorentzVector Mom;
    gMC->TrackMomentum(Mom);
    Double_t xmean = (fPos.X()+Pos.X())/2. ;
    Double_t ymean = (fPos.Y()+Pos.Y())/2. ;
    Double_t zmean = (fPos.Z()+Pos.Z())/2. ;

    AddHit(fTrackID, uniqueId, TVector3(xmean, ymean,  zmean),
           TVector3(fMom.Px(), fMom.Py(), fMom.Pz()), fTime, fLength,
           fELoss,pdgCode,TVector3(Pos.X(), Pos.Y(), Pos.Z()),
	   TVector3(Mom.Px(), Mom.Py(), Mom.Pz()) );

    // Increment number of veto det points in TParticle
    ShipStack* stack = dynamic_cast<ShipStack*>(gMC->GetStack());
    stack->AddPoint(kUpstreamTagger);
  }

  return kTRUE;
}



void UpstreamTagger::EndOfEvent()
{
  fUpstreamTaggerPointCollection->Clear();
}



void UpstreamTagger::Register()
{

  /** This will create a branch in the output tree called
      UpstreamTaggerPoint, setting the last parameter to kFALSE means:
      this collection will not be written to the file, it will exist
      only during the simulation.
  */

  FairRootManager::Instance()->Register("UpstreamTaggerPoint", "UpstreamTagger",
                                        fUpstreamTaggerPointCollection, kTRUE);
}



TClonesArray* UpstreamTagger::GetCollection(Int_t iColl) const
{
  if (iColl == 0) { return fUpstreamTaggerPointCollection; }
  else { return NULL; }
}



void UpstreamTagger::Reset()
{
  fUpstreamTaggerPointCollection->Clear();
}


void UpstreamTagger::SetzPositions(Double_t z1, Double_t z2, Double_t z3, Double_t z4)
{
     f_T1_z = z1;                                               //!  z-position of tracking station 1
     f_T2_z = z2;                                               //!  z-position of tracking station 2
     f_T3_z = z3;                                               //!  z-position of tracking station 3
     f_T4_z = z4;                                               //!  z-position of tracking station 4
}

void UpstreamTagger::SetApertureArea(Double_t width, Double_t height)
{
     f_aperture_width = width;                                  //!  Aperture width (x)
     f_aperture_height = height;                                //!  Aperture height (y)
}

void UpstreamTagger::SetStrawDiameter(Double_t outer_straw_diameter, Double_t wall_thickness)
{
     f_outer_straw_diameter = outer_straw_diameter;             //!  Outer straw diameter
     f_inner_straw_diameter =
       outer_straw_diameter - 2 * wall_thickness;               //!  Inner straw diameter
}

void UpstreamTagger::SetStrawPitch(Double_t straw_pitch, Double_t layer_offset)
{
     f_straw_pitch = straw_pitch;                               //!  Distance (y) between straws in a layer
     f_offset_layer = layer_offset;                             //!  Offset (y) of straws between layers
}

void UpstreamTagger::SetDeltazLayer(Double_t delta_z_layer)
{
     f_delta_z_layer = delta_z_layer;                           //!  Distance (z) between layers
}

void UpstreamTagger::SetStereoAngle(Double_t stereo_angle)
{
     f_view_angle = stereo_angle;                                //!  Stereo view angle
}

void UpstreamTagger::SetWireThickness(Double_t wire_thickness)
{
     f_wire_thickness = wire_thickness;                         //!  Sense wire thickness
}

void UpstreamTagger::SetDeltazView(Double_t delta_z_view)
{
     f_delta_z_view = delta_z_view;                             //!  Distance (z) between stereo views
}

void UpstreamTagger::SetFrameMaterial(TString frame_material)
{
     f_frame_material = frame_material;                         //!  Structure frame material
}

void UpstreamTagger::SetStationEnvelope(Double_t x, Double_t y, Double_t z)
{
     f_station_width = x;                                       //!  Station envelope width (x)
     f_station_height = y;                                      //!  Station envelope height (y)
     f_station_length = z;                                      //!  Station envelope length (z)
}

void UpstreamTagger::ConstructGeometry()
{

 /** If you are using the standard ASCII input for the geometry
      just copy this and use it for your detector, otherwise you can
      implement here you own way of constructing the geometry. */

    std::cout<<"Making a geometry"<<std::endl;

    TGeoVolume *top               = gGeoManager->GetTopVolume();
    InitMedium("air");
    TGeoMedium *air               = gGeoManager->GetMedium("air");
    InitMedium("ShipSens");
    TGeoMedium *Se                = gGeoManager->GetMedium("ShipSens");
    InitMedium("aluminium");
    TGeoMedium *Al                = gGeoManager->GetMedium("aluminium");
    InitMedium("mylar");
    TGeoMedium *mylar             = gGeoManager->GetMedium("mylar");
    InitMedium("STTmix8020_1bar");
    TGeoMedium *sttmix8020_1bar   = gGeoManager->GetMedium("STTmix8020_1bar");
    InitMedium("tungsten");
    TGeoMedium *tungsten          = gGeoManager->GetMedium("tungsten");
    InitMedium(f_frame_material);
    TGeoMedium *FrameMatPtr       = gGeoManager->GetMedium(f_frame_material);
    InitMedium(fMedium.c_str());
    TGeoMedium* med = gGeoManager->GetMedium(fMedium.c_str());

    gGeoManager->SetVisLevel(4);
    gGeoManager->SetTopVisible();

    // Epsilon to avoid overlapping volumes
    Double_t eps = 0.01;
    // Straw (half) length
    Double_t straw_length = f_aperture_width + 2. * eps;
    // Width of frame: standard HEA 500 I-beam width
    Double_t frame_width =  49.;
    // Offset due to floor space limitation
    Double_t floor_offset = 14.;

    Double_t rmin, rmax, T_station_z;

    // Arguments of boxes are half-lengths
    TGeoBBox* detbox1 = new TGeoBBox(
    "ubt_detbox1", (f_station_width + frame_width), (f_station_height + frame_width - floor_offset / 2.), f_station_length);
    TGeoBBox* detbox2 = new TGeoBBox(
	"ubt_detbox2",
    (straw_length + eps),
	(f_station_height + TMath::Tan(f_view_angle * TMath::Pi() / 180.0) * straw_length * 2 + f_offset_layer / TMath::Cos(f_view_angle * TMath::Pi() / 180.0) + eps),
	f_station_length + eps);

    std::cout<<"frame height: "<<f_station_height + frame_width - floor_offset / 2.<<" - frame width: "<<f_station_width + frame_width<<std::endl;
    std::cout<<"f_station_width: "<<f_station_width<<" - f_station_height: "<<f_station_height<<" - frame_width: "<<frame_width<<std::endl;
    TGeoTranslation* move_up = new TGeoTranslation("move_up", 0, floor_offset / 2., 0);
    move_up->RegisterYourself();

    // Composite shape to create frame
    TGeoCompositeShape* detcomp1 = new TGeoCompositeShape("ubt_detcomp1", "ubt_detbox1:move_up - ubt_detbox2");

    // Volume: straw
    rmin = f_inner_straw_diameter / 2.;
    rmax = f_outer_straw_diameter / 2.;
    // Third argument is half-length of tube
    TGeoTube *straw_tube = new TGeoTube("ubt_straw", rmin, rmax, straw_length);
    TGeoVolume *straw = new TGeoVolume("ubt_straw", straw_tube, mylar);
    straw->SetLineColor(4);
    straw->SetVisibility(kTRUE);

    // Volume: gas
    rmin = f_wire_thickness / 2. + eps;
    rmax = f_inner_straw_diameter / 2. - eps;
    TGeoTube *gas_tube = new TGeoTube("gas", rmin, rmax, straw_length - 2. * eps);
    TGeoVolume *gas = new TGeoVolume("gas", gas_tube, sttmix8020_1bar);
    gas->SetLineColor(5);
    // Only the gas is sensitive
    AddSensitiveVolume(gas);

    // Volume: wire
    rmin = 0.;
    rmax = f_wire_thickness / 2.;
    TGeoTube *wire_tube = new TGeoTube("wire", rmin, rmax, straw_length - 4. * eps);
    TGeoVolume *wire = new TGeoVolume("wire", wire_tube, tungsten);
    wire->SetLineColor(6);

    // Tracking stations
    // statnb = station number; vnb = view number; lnb = layer number; snb = straw number

    // Station box to contain all components
    TGeoBBox* statbox = new TGeoBBox("ubt_statbox", 0.01*f_station_width, 0.01*f_station_height - floor_offset / 2., 0.01*f_station_length);

    f_frame_material.ToLower();

    for (Int_t statnb = 1; statnb < 2; statnb++) {
        // Tracking station loop
        TString nmstation = "UBT";
        std::stringstream ss;
        ss << statnb;
        nmstation = nmstation + ss.str();
        switch (statnb) {
            case 1:
	        T_station_z = f_T1_z;
                break;
            default:
                T_station_z = f_T1_z;
        }

	TGeoVolume* vol = new TGeoVolume(nmstation, statbox, med);
	// z-translate the station to its (absolute) position
	top->AddNode(vol, statnb, new TGeoTranslation(0, floor_offset / 2., T_station_z));

	TGeoVolume* statframe = new TGeoVolume(nmstation + "_frame", detcomp1, FrameMatPtr);
	vol->AddNode(statframe, statnb * 1e6, new TGeoTranslation(0, -floor_offset / 2., 0));
	statframe->SetLineColor(kOrange);

        for (Int_t vnb = 0; vnb < 4; vnb++) {
            // View loop
            TString nmview;
            Double_t angle;
	        Double_t stereo_growth;
    	    Double_t stereo_pitch;
	        Double_t offset_layer;
	        Int_t straws_per_layer;

            switch (vnb) {
                case 0:
                    angle = 0.;
                    nmview = nmstation + "_y1";
                    break;
                case 1:
                    angle = f_view_angle;
                    nmview = nmstation + "_u";
                    break;
                case 2:
                    angle = -f_view_angle;
                    nmview = nmstation + "_v";
                    break;
                case 3:
                    angle = 0.;
                    nmview = nmstation + "_y2";
                    break;
                default:
                    angle = 0.;
                    nmview = nmstation + "_y1";
            }

	    // Adjustments in the stereo views
	    // stereo_growth: extension of stereo views beyond aperture
	    // stereo_pitch: straw pitch in stereo views
	    // offset_layer: layer offset in stereo views
	    // straws_per_layer: number of straws in one layer with stereo extension
	    // If angle == 0., all numbers return the case of non-stereo views.
	    stereo_growth = TMath::Tan(TMath::Abs(angle) * TMath::Pi() / 180.0) * straw_length;
	    stereo_pitch = f_straw_pitch / TMath::Cos(TMath::Abs(angle) * TMath::Pi() / 180.0);
	    offset_layer = f_offset_layer / TMath::Cos(TMath::Abs(angle) * TMath::Pi() / 180.0);
	    straws_per_layer = std::ceil(2 * (f_aperture_height + stereo_growth) / stereo_pitch);

            for (Int_t lnb = 0; lnb < 2; lnb++) {
                // Layer loop
                TString nmlayer = nmview + "_layer_";
                nmlayer += lnb;
                TGeoBBox* layer = new TGeoBBox(
                    "layer box", straw_length + eps / 4, f_aperture_height + stereo_growth * 2 + offset_layer + eps / 4, f_outer_straw_diameter / 2. + eps / 4);
                TGeoVolume* layerbox = new TGeoVolume(nmlayer, layer, med);

                // The layer box sits in the viewframe.
                // Hence, z-translate the layer w.r.t. the view
                vol->AddNode(layerbox, statnb * 1e6 + vnb * 1e5 + lnb * 1e4, new TGeoTranslation(0, -floor_offset / 2., (vnb - 3. / 2.) * f_delta_z_view + (lnb - 1. / 2.) * f_delta_z_layer));

                TGeoRotation r6s;
                TGeoTranslation t6s;
                for (Int_t snb = 1; snb <= straws_per_layer; snb++) {
                    // Straw loop
		    // y-translate the straw to its position
		    t6s.SetTranslation(0, f_aperture_height + stereo_growth - (snb - 1. / 2.) * stereo_pitch + lnb * offset_layer, 0);
		    // Rotate the straw with stereo angle
                    r6s.SetAngles(90 + angle, 90, 0);
                    TGeoCombiTrans c6s(t6s, r6s);
                    TGeoHMatrix* h6s = new TGeoHMatrix(c6s);
                    layerbox->AddNode(straw, statnb * 1e6 + vnb * 1e5 + lnb * 1e4 + 1e3 + snb, h6s);
                    layerbox->AddNode(gas, statnb * 1e6 + vnb * 1e5 + lnb * 1e4 + 2e3 + snb, h6s);
                    layerbox->AddNode(wire, statnb * 1e6 + vnb * 1e5 + lnb * 1e4 + 3e3 + snb, h6s);
                    // End of straw loop
                }

		if (lnb == 1) {
		   // Add one more straw at the bottom of the second layer to cover aperture entirely
		   t6s.SetTranslation(0, f_aperture_height + stereo_growth - (straws_per_layer - 1. / 2.) * stereo_pitch - lnb * offset_layer, 0);
                   r6s.SetAngles(90 + angle, 90, 0);
                   TGeoCombiTrans c6s(t6s, r6s);
                   TGeoHMatrix* h6s = new TGeoHMatrix(c6s);
                   layerbox->AddNode(straw, statnb * 1e6 + vnb * 1e5 + lnb * 1e4 + 1e3 + straws_per_layer + 1, h6s);
                   layerbox->AddNode(gas, statnb * 1e6 + vnb * 1e5 + lnb * 1e4 + 2e3 + straws_per_layer + 1, h6s);
                   layerbox->AddNode(wire, statnb * 1e6 + vnb * 1e5 + lnb * 1e4 + 3e3 + straws_per_layer + 1, h6s);
		}
                // End of layer loop
            }
            // End of view loop
        }
    // A layer of plastic scintillator detector
    InitMedium("polyvinyltoluene");
    TGeoMedium *Vacuum_box =gGeoManager->GetMedium("polyvinyltoluene");
    UpstreamTagger_plastic = gGeoManager->MakeBox("Upstream_Tagger_Plastic", Vacuum_box, 150, 150, 1);
    UpstreamTagger_plastic->SetLineColor(kGreen);
    vol->AddNode(UpstreamTagger_plastic, 1, new TGeoTranslation(0.0, 0.0, -50));
    std::cout<<"geometry constructed"<<std::endl;
    // End of tracking station loop
    }
/*
  TGeoVolume *top = gGeoManager->GetTopVolume();

  //////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////

  InitMedium("STTmix8020_1bar");
  TGeoMedium *Vacuum_box =gGeoManager->GetMedium("STTmix8020_1bar");
  ///////////////////////////////////////////////////////////////////

  // Adding UBT Extension
  if (!Vacuum_box) {
      Fatal("ConstructGeometry", "Medium 'vacuum' not found.");
  }

  UpstreamTagger_fulldet = gGeoManager->MakeBox("Upstream_Tagger", Vacuum_box, xbox_fulldet/2.0, ybox_fulldet/2.0, zbox_fulldet/2.0);
  UpstreamTagger_fulldet->SetLineColor(kGreen);

  UpstreamTagger_fulldet_1 = gGeoManager->MakeBox("Upstream_Tagger_firstLayer", Vacuum_box, xbox_fulldet/2.0, ybox_fulldet/2.0, zbox_fulldet/2.0);
  UpstreamTagger_fulldet_1->SetLineColor(kGreen);
  UpstreamTagger_fulldet_3 = gGeoManager->MakeBox("Upstream_Tagger_lastLayer", Vacuum_box, xbox_fulldet/2.0, ybox_fulldet/2.0, zbox_fulldet/2.0);
  UpstreamTagger_fulldet_3->SetLineColor(kGreen);

  top->AddNode(UpstreamTagger_fulldet, 1, new TGeoTranslation(0.0, 0.0, det_zPos));
  top->AddNode(UpstreamTagger_fulldet_1, 1, new TGeoTranslation(0.0, 0.0, det_zPos-4.*cm));
  top->AddNode(UpstreamTagger_fulldet_3, 1, new TGeoTranslation(0.0, 0.0, det_zPos+4.*cm));

  AddSensitiveVolume(UpstreamTagger_fulldet);
  AddSensitiveVolume(UpstreamTagger_fulldet_1);
  AddSensitiveVolume(UpstreamTagger_fulldet_3);

  cout << " Z Position (Upstream Tagger1) " << det_zPos << endl;
  //////////////////////////////////////////////////////////////////

  return;
*/
}



UpstreamTaggerPoint* UpstreamTagger::AddHit(Int_t trackID, Int_t detID,
			TVector3 pos, TVector3 mom,
			Double_t time, Double_t length,
			Double_t eLoss, Int_t pdgCode,TVector3 Lpos, TVector3 Lmom)
{
  TClonesArray& clref = *fUpstreamTaggerPointCollection;
  Int_t size = clref.GetEntriesFast();

  return new(clref[size]) UpstreamTaggerPoint(trackID, detID, pos, mom,
		         time, length, eLoss, pdgCode,Lpos,Lmom);
}

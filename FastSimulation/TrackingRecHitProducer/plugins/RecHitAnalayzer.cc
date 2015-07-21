// -*- C++ -*-
//
// Package:    FastSimulation/RecHitAnalayzer
// Class:      RecHitAnalayzer
// 
/**\class RecHitAnalayzer RecHitAnalayzer.cc FastSimulation/RecHitAnalayzer/plugins/RecHitAnalayzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Abideh Jafari
//         Created:  Sat, 18 Jul 2015 06:38:24 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/DetId/interface/DetId.h"

//RecHits
#include "DataFormats/TrackerRecHit2D/interface/SiTrackerGSRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiTrackerGSRecHit2DCollection.h"
#include "DataFormats/Common/interface/OwnVector.h" 
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

// PSimHits
#include "SimDataFormats/CrossingFrame/interface/MixCollection.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"

//ROOT
#include "TFile.h"
#include "TH1.h"

//C++
#include <vector>
#include <iostream>
//
// class declaration
//

class RecHitAnalayzer : public edm::EDAnalyzer {
   public:
      explicit RecHitAnalayzer(const edm::ParameterSet&);
      ~RecHitAnalayzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
      
      typedef std::vector<PSimHit> PSimHitCollection;


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;


      // ----------member data ---------------------------
      edm::ParameterSet iPset;
      edm::InputTag track_label, simhit_label;
      int verbose;
      TFile * fout;
      TH1D * DX;
      TH1D * DY;
      TH1D * DZ;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
RecHitAnalayzer::RecHitAnalayzer(const edm::ParameterSet& iConfig): iPset(iConfig), 
	track_label( iPset.getParameter<edm::InputTag>("track_label")),
	simhit_label( iPset.getParameter<edm::InputTag>("simhit_label")),
	verbose( iPset.getParameter<int>("verbose"))

{

}


RecHitAnalayzer::~RecHitAnalayzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each iEvent  ------------
void
RecHitAnalayzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;

   Handle<reco::TrackCollection> track_handle;
   iEvent.getByLabel(track_label, track_handle);
   
   edm::Handle<PSimHitCollection> simhit_handle; 
   iEvent.getByLabel(simhit_label, simhit_handle);  
   int iTrack = 0;
   reco::TrackCollection::const_iterator itTrk = track_handle->begin();
   reco::TrackCollection::const_iterator trkEnd = track_handle->end();
   for(; itTrk != trkEnd; itTrk++){
	trackingRecHit_iterator itRecHit = itTrk->recHitsBegin();
    	trackingRecHit_iterator rhEnd = itTrk->recHitsEnd();
	int iRechit = 0;
	for(; itRecHit != rhEnd; itRecHit++){
		if(!(*itRecHit)->isValid()) continue;
		float xRec = (*itRecHit)->localPosition().x();
                float yRec = (*itRecHit)->localPosition().y();
                float zRec = (*itRecHit)->localPosition().z();
		if(verbose > 3)
			std::cout<<"Track "<< iTrack<< " and RecHit "<<iRechit<<": "<<xRec <<"\t" << yRec <<"\t"<< zRec <<std::endl;
		if(verbose > 0)
			std::cout<<"Looking for detId ";
		DetId detId = (*itRecHit)->geographicalId().rawId();
		if(verbose > 0)
			std::cout<<detId<<endl;
		const PSimHit* simHit = NULL;
		double minDR = 99999999;
		PSimHitCollection::const_iterator itSimHit = simhit_handle->begin();
		PSimHitCollection::const_iterator itSimHitEnd = simhit_handle->begin();
		for (; itSimHit!= itSimHitEnd; itSimHitEnd++) {		
			if(verbose > 0)
				std::cout<<"comparing detID with unitId: "<<detId << " vs. "<<itSimHit->detUnitId()<<endl;
			if (detId == itSimHit->detUnitId()){
				float xSim = itSimHit->localPosition().x();
				float ySim = itSimHit->localPosition().y();
				float zSim = itSimHit->localPosition().z();
				double DR = sqrt(pow((xRec-xSim),2)+pow((yRec-ySim),2)+pow((zRec-zSim),2));
				if(DR < minDR){
					minDR = DR;
					//simHit = new PSimHit(*itSimHit);
					simHit = &(*itSimHit);
				}
			}
		}    
		if(verbose > 0)
			std::cout<<"After simHit loop, th eclosest simHit is: "<<simHit<<endl;
		if(simHit == NULL)
			continue;
		if(verbose > 0)
			std::cout<<"simHit coordinates are: ";
		float xSim = simHit->localPosition().x();
		float ySim = simHit->localPosition().y();
		float zSim = simHit->localPosition().z();
		if(verbose > 0)
			std::cout<<xSim<<"\t"<<ySim<<"\t"<<zSim<<endl;
		DX->Fill(xRec-xSim);
		DY->Fill(yRec-ySim);			
		DZ->Fill(zRec-zSim);
	}
   }
}


// ------------ method called once each job just before starting iEvent loop  ------------
void 
RecHitAnalayzer::beginJob()
{
   DX = new TH1D("DX", "#Delta x",2000,-10,10);
   DY = new TH1D("DY", "#Delta y",2000,-10,10);
   DZ = new TH1D("DZ", "#Delta z",2000,-10,10);
}

// ------------ method called once each job just after ending the iEvent loop  ------------
void 
RecHitAnalayzer::endJob() 
{
   fout = new TFile("rechitanalyzer.root","recreate");
   fout->cd();
   DX->Write();
   DY->Write();
   DZ->Write();
   fout->Close();
}

// ------------ method called when starting to processes a run  ------------
/*
void 
RecHitAnalayzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
RecHitAnalayzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
RecHitAnalayzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
RecHitAnalayzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
RecHitAnalayzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(RecHitAnalayzer);

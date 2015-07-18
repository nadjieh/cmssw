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

// PSimHits
#include "SimDataFormats/CrossingFrame/interface/MixCollection.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"

//ROOT
#include "TFile.h"
#include "TH1.h"

//C++
#include <vector>
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

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
      edm::ParameterSet iPset;
      edm::InputTag RecHits_label;
      TFile * fout;
      TH1D * DX;
      TH1D * DY;
      TH1D * DZ;
      // PSimHits
      std::vector<edm::InputTag> trackerContainers;
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
RecHitAnalayzer::RecHitAnalayzer(const edm::ParameterSet& iConfig): iPset(iConfig), RecHits_label( iPset.getParameter<edm::InputTag>("RecHits_label") )

{
  trackerContainers.clear();
  trackerContainers = iPset.getParameter<std::vector<edm::InputTag> >("ROUList");

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
   Handle<SiTrackerFullGSRecHit2DCollection> RecHits_handle;
   iEvent.getByLabel(RecHits_label, RecHits_handle);
   
   edm::Handle<CrossingFrame<PSimHit> > cf_simhit; 
   std::vector<const CrossingFrame<PSimHit> *> cf_simhitvec;
   for(uint32_t i=0; i<trackerContainers.size(); i++){
		iEvent.getByLabel(trackerContainers[i], cf_simhit);
    	cf_simhitvec.push_back(cf_simhit.product());
   }
   std::auto_ptr<MixCollection<PSimHit> > allTrackerHits(new MixCollection<PSimHit>(cf_simhitvec)); 
   
   // loop on RecHits
   unsigned int iRecHit = 0;
   const std::vector<DetId> theDetIds = RecHits_handle->ids();
   // loop over detunits
   for ( std::vector<DetId>::const_iterator iDetId = theDetIds.begin(); iDetId != theDetIds.end(); iDetId++ ) {
   	unsigned int detid = (*iDetId).rawId();
   	if(detid!=999999999){ // valid detector
      	SiTrackerFullGSRecHit2DCollection::range RecHit_range = RecHits_handle->get((*iDetId));
      	SiTrackerFullGSRecHit2DCollection::const_iterator RecHit_rangeIteratorBegin = RecHit_range.first;
      	SiTrackerFullGSRecHit2DCollection::const_iterator RecHit_rangeIteratorEnd   = RecHit_range.second;
      	SiTrackerFullGSRecHit2DCollection::const_iterator iterRecHit = RecHit_rangeIteratorBegin;
      	// loop over RecHits of the same detector
      	for(iterRecHit = RecHit_rangeIteratorBegin; iterRecHit != RecHit_rangeIteratorEnd; ++iterRecHit) {
      		iRecHit++;
			// search the associated original PSimHit
			PSimHit* simHit = NULL;
			int simHitNumber = (*iterRecHit).simhitId();
			int simHitCounter = -1;
			for (MixCollection<PSimHit>::iterator isim=(*allTrackerHits).begin(); isim!= (*allTrackerHits).end(); isim++) {
	  			simHitCounter++;
	  			if(simHitCounter == simHitNumber) {
	    			simHit = const_cast<PSimHit*>(&(*isim));
	    			break;
	  			}
			}    
	
			float xRec = (*iterRecHit).localPosition().x();
			float yRec = (*iterRecHit).localPosition().y();
			float zRec = (*iterRecHit).localPosition().z();
			float xSim = simHit->localPosition().x();
			float ySim = simHit->localPosition().y();
			float zSim = simHit->localPosition().z();
			DX->Fill(xRec-xSim);
			DY->Fill(yRec-ySim);			
			DZ->Fill(zRec-zSim);
		}
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

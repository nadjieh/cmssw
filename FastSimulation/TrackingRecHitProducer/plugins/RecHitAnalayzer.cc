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
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"

//RecHits
#include "DataFormats/TrackerRecHit2D/interface/SiTrackerGSRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiTrackerGSRecHit2DCollection.h"
#include "DataFormats/TrackerRecHit2D/interface/SiTrackerGSMatchedRecHit2D.h"
#include "DataFormats/Common/interface/OwnVector.h" 
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

// PSimHits
#include "SimDataFormats/CrossingFrame/interface/MixCollection.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"

//ROOT
#include "TFile.h"
#include "TH1.h"
#include "TString.h"
#include "TDirectory.h"

//C++
#include <vector>
#include <iostream>
#include <sstream>
#include <string>
using namespace std;
//
// class declaration
//
class AuxilaryPlots{
public:
      AuxilaryPlots(TString name_): name(name_){
	DX = new TH1D("DX", "#Delta x",200,-0.1,0.1);
	DY = new TH1D("DYWide", "#Delta y",2000,-10,10);
	DYLR = new TH1D("DY", "#Delta y",200,-0.1,0.1);
	DZ = new TH1D("DZ", "#Delta z",200,-0.01,0.01);
	recErrXX = new TH1D("recXX", "xx",200,-0.1,0.1);
	recErrXY = new TH1D("recXY", "xy",200,-0.1,0.1);
	recErrYY = new TH1D("recYY", "yy",200,-0.1,0.1);
      }
      void Fill(double dx, double dy, double dz){
	DX->Fill(dx);
        DY->Fill(dy);
        DYLR->Fill(dy);
        DZ->Fill(dz);
      }
      void FillRecHitErr(double xx, double xy, double yy){
	recErrXX->Fill(sqrt(xx));
        recErrXY->Fill(sqrt(xy));
        recErrYY->Fill(sqrt(yy));
      }
      void Write(TDirectory * f){
	f->cd();
	TDirectory * dir = f->mkdir(name); 
	dir->cd();
        DX->Write();
        DY->Write();
        DYLR->Write();
        DZ->Write();
        dir->mkdir("recHitErr")->cd();
        recErrXX->Write();
        recErrXY->Write();
        recErrYY->Write();
	dir->cd();
	f->cd();
      }
      TString name;
private:
      TH1D * DX;
      TH1D * DY;
      TH1D * DZ;
      TH1D * DYLR;
      TH1D * recErrXX;
      TH1D * recErrXY;
      TH1D * recErrYY;
};
class SubDetLayersHistograms{
public:
	SubDetLayersHistograms(TString name_, unsigned int n): name(name_),nLayer(n){
		stringstream s;
		for(unsigned int i = 0; i < nLayer; i++){
			s.str("");
			s << name <<"_"<<i+1;
			layerPlots.push_back(new AuxilaryPlots(s.str().c_str()));		
		}
	}
	~SubDetLayersHistograms(){}
	void Fill(int ID, double dx, double dy, double dz){
		if(ID > (int)layerPlots.size()){
			cout<<"INVALID LAYER FOR "<<name<<endl;	
			return;
		}
		layerPlots[ID-1]->Fill(dx,dy,dz);
	}
	void FillRecHitErr(int ID, double xx, double xy, double yy){
	        if(ID > (int)layerPlots.size()){
                        cout<<"INVALID LAYER FOR "<<name<<endl;
                        return;
                }
                layerPlots[ID-1]->FillRecHitErr(xx,xy,yy);
        }
	void Write(TDirectory * f){
		f->cd();
		TDirectory * dir = f->mkdir(name);
		for(unsigned int i = 0; i < nLayer; i++){
			layerPlots[i]->Write(dir);
		}
	}
private:
	TString name;
	unsigned int nLayer;
	std::vector<AuxilaryPlots*> layerPlots;
};

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
      void FillAuxilaryPlots(DetId, unsigned int, const TrackerTopology*, double, double, double, double, double, double);

      // ----------member data ---------------------------
      edm::ParameterSet iPset;
      edm::InputTag track_label, simhit_label;
      int verbose;
      bool isFS;
      std::string fname;
      TFile * fout;
      SubDetLayersHistograms * PixleBarrel;
      SubDetLayersHistograms * PixleFwd;
      SubDetLayersHistograms * TIB;
      SubDetLayersHistograms * TID;
      SubDetLayersHistograms * TOB;
      SubDetLayersHistograms * TEC;
      TH1D * NSimRecSameDetId;
      TH1D * Rfs;
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
	verbose( iPset.getParameter<int>("verbose")),
	isFS( iPset.getParameter<bool>("isFastSimOnly")),
	fname( iPset.getParameter<string>("outfile"))

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

   edm::ESHandle<TrackerTopology> tTopoHand;
   iSetup.get<TrackerTopologyRcd>().get(tTopoHand);
   const TrackerTopology *tTopo=tTopoHand.product();

   int iTrack = 0;
   reco::TrackCollection::const_iterator itTrk = track_handle->begin();
   reco::TrackCollection::const_iterator trkEnd = track_handle->end();
   for(; itTrk != trkEnd; itTrk++){
	trackingRecHit_iterator itRecHit = itTrk->recHitsBegin();
    	trackingRecHit_iterator rhEnd = itTrk->recHitsEnd();
	int iRechit = 0;
	for(; itRecHit != rhEnd; itRecHit++){
		if(!(*itRecHit)->isValid()) continue;
		TrackingRecHit * tempHit = (*itRecHit)->clone();
		float xRec = tempHit->localPosition().x();
                float yRec = tempHit->localPosition().y();
                float zRec = tempHit->localPosition().z();
		float xxRecErr = tempHit->localPositionError().xx();
                float xyRecErr = tempHit->localPositionError().xy();
                float yyRecErr = tempHit->localPositionError().yy();
		if(verbose > 3)
			std::cout<<"Track "<< iTrack<< " and RecHit "<<iRechit<<": "<<xRec <<"\t" << yRec <<"\t"<< zRec <<std::endl;
		if(verbose > 3)
			std::cout<<"Looking for detId ";
		DetId detId = tempHit->geographicalId().rawId();
		delete tempHit;
		
//		if(verbose > 3)
//			std::cout<<detId<<endl;
		const PSimHit* simHit = NULL;
		double minDR = 99999999;
		PSimHitCollection::const_iterator itSimHit = simhit_handle->begin();
		PSimHitCollection::const_iterator itSimHitEnd = simhit_handle->end();
		if(verbose > 3)
                        std::cout<<simhit_handle->size()<<endl;
		int nSimRecSameDetId = 0;
		for (; itSimHit!= itSimHitEnd; itSimHit++) {		
//			if(verbose > 6)
//				std::cout<<"comparing detID with unitId: "<<detId << " vs. "<<itSimHit->detUnitId()<<endl;
			if (detId == itSimHit->detUnitId()){
				nSimRecSameDetId++;
				float xSim = itSimHit->localPosition().x();
				float ySim = itSimHit->localPosition().y();
				float zSim = itSimHit->localPosition().z();
				if(verbose > 0)
					cout<<"SimHit Coordinates: "<<xSim <<"\t"<<ySim<<"\t"<<zSim<<endl;
				double DR = sqrt(pow((xRec-xSim),2)+pow((yRec-ySim),2)+pow((zRec-zSim),2));
				if(verbose > 0)
					cout<< "DR is "<<DR<<endl;
				if(DR < minDR){
					if(verbose > 0)
	                                        cout<< "keep new simHit :-) "<<endl;
					minDR = DR;
					//simHit = new PSimHit(*itSimHit);
					simHit = &(*itSimHit);
				}
			}
		}    
		if(verbose > 0)
			std::cout<<"After simHit loop, the closest simHit is: "<<simHit<<endl;
		if(isFS){
			const SiTrackerGSMatchedRecHit2D * rechit = (const SiTrackerGSMatchedRecHit2D*) (*itRecHit);
			if(verbose > 0){
				cout<<"From SiTrackerGSMatched:\n";
				if(rechit->isMatched()) {
					cout<<"\tThe RecHit is Matched with SimHit ======= "<<endl;
					cout<<"position: "<<rechit->localPosition().x()<<"\t"<<rechit->localPosition().y()<<"\t"<<rechit->localPosition().z()<<endl;
				} else  cout<<"\tThe RecHit is NOT Matched with SimHit"<<endl;
			}
			if(rechit->isMatched()){
				Rfs->Fill(sqrt(pow(rechit->localPosition().x(),2)+pow(rechit->localPosition().y(),2)+pow(rechit->localPosition().z(),2)));
			}
		}
		if(simHit == NULL)
			continue;
		if(verbose > 0)
			std::cout<<"simHit coordinates are: ";
		float xSim = simHit->localPosition().x();
		float ySim = simHit->localPosition().y();
		float zSim = simHit->localPosition().z();
		if(verbose > 0)
			std::cout<<xSim<<"\t"<<ySim<<"\t"<<zSim<<endl;
		double dx = xRec-xSim;
                double dy = yRec-ySim;
                double dz = zRec-zSim;

		//Layer definitions and ID
        	unsigned int subDetId = DetId(simHit->detUnitId()).subdetId();				
		detId = simHit->detUnitId();
	        FillAuxilaryPlots(detId, subDetId, tTopo,dx,dy,dz,xxRecErr,xyRecErr,yyRecErr);
		NSimRecSameDetId->Fill(nSimRecSameDetId);
		if(verbose > 3)
			cout<<"================================= END OF RECHITS ==============================="<<endl;
	}
   }
   if(verbose > 3)
	cout<<"================================= END OF TRACKS ==============================="<<endl;
}


// ------------ method called once each job just before starting iEvent loop  ------------
void 
RecHitAnalayzer::beginJob()
{
   PixleBarrel = new SubDetLayersHistograms("PixleBarrel",3);
   PixleFwd = new SubDetLayersHistograms("PixleFwd",2);
   TIB = new SubDetLayersHistograms("TIB", 4);
   TID = new SubDetLayersHistograms("TID",3);
   TOB = new SubDetLayersHistograms("TOB",6);
   TEC = new SubDetLayersHistograms("TEC",7);
   NSimRecSameDetId = new TH1D("NSimRecSameDetId","N (SimHit,RecHit)_{same detid}", 11, -5.5, 5.5);
   Rfs = new TH1D("Rfs","Matched SiTrack2D position (R)", 2000, -10, 10);
}

// ------------ method called once each job just after ending the iEvent loop  ------------
void 
RecHitAnalayzer::endJob() 
{
   fout = new TFile(fname.c_str(),"recreate");
   fout->cd();
   PixleBarrel->Write(fout);
   PixleFwd->Write(fout);
   TIB->Write(fout);
   TID->Write(fout);
   TOB->Write(fout);
   TEC->Write(fout);
   NSimRecSameDetId->Write();
   Rfs->Write();
   fout->Close();
}

void
RecHitAnalayzer::FillAuxilaryPlots(DetId detId, unsigned int subDetId, const TrackerTopology * tt, double dx, double dy, double dz
				, double recXX, double recXY, double recYY){
   using namespace std;
   if(subDetId > 6)
	return;
   switch(subDetId){
	case 1:
	    {
                if(verbose > 3 )
			cout << "\tPixel Barrel Layer " << tt->pxbLayer(detId)<<endl;
		PixleBarrel->Fill((int)tt->pxbLayer(detId),dx,dy,dz);
		PixleBarrel->FillRecHitErr((int)tt->pxbLayer(detId),recXX,recXY,recYY);
		break;
	    }
        case 2:
	    {
                if(verbose > 3 )
                	cout << "\tPixel Forward Disk " << tt->pxfDisk(detId)<<endl;
		PixleFwd->Fill((int)tt->pxfDisk(detId),dx,dy,dz);
		PixleFwd->FillRecHitErr((int)tt->pxfDisk(detId),recXX,recXY,recYY);
		break;
	    }
	case 3:
	    {
                if(verbose > 3 )
                	cout << "\tTIB Layer " << tt->tibLayer(detId)<<endl;
		TIB->Fill((int)tt->tibLayer(detId),dx,dy,dz);
		TIB->FillRecHitErr((int)tt->tibLayer(detId),recXX,recXY,recYY);
		break;
	    }
        case 4:
	    {
                if(verbose > 3 )
                	cout << "\tTID " << tt->tidRing(detId)<<endl;
		TID->Fill((int)tt->tidRing(detId), dx,dy,dz);
		TID->FillRecHitErr((int)tt->tidRing(detId),recXX,recXY,recYY);
		break;
	    }
	case 5:
	    {
                if(verbose > 3 )
             		cout << "\tTOB Layer " << tt->tobLayer(detId)<<endl;
		TOB->Fill((int)tt->tobLayer(detId), dx,dy,dz);
		TOB->FillRecHitErr((int)tt->tobLayer(detId),recXX,recXY,recYY);
		break;
	    }
        case 6:
	    {
                if(verbose > 3 )
                	cout << "\tTEC " << tt->tecRing(detId)<<endl;
		TEC->Fill((int)tt->tecRing(detId), dx,dy,dz);
		TEC->FillRecHitErr((int)tt->tecRing(detId),recXX,recXY,recYY);
		break;
	    }
	default:
	    {
   	        cout<< "\t INVALID SUBDETECTOR =========================="<<endl;
	        return;
		break;
	    }
    }    

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

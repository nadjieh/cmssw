#ifndef TOPDILEPTONOFFLINEDQM
#define TOPDILEPTONOFFLINEDQM

#include <string>
#include <vector>

#include "FWCore/Framework/interface/Event.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DQM/Physics/interface/TopDQMHelpers.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Framework/interface/EDConsumerBase.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
/**
   \class   MonitorEnsemble TopDQMHelpers.h "DQM/Physics/interface/TopDQMHelpers.h"

   \brief   Helper class to define histograms for monitoring of muon/electron/jet/met quantities.

   Helper class to contain histograms for the monitoring of muon/electron/jet/met quantities. 
   This class can be instantiated several times after several event selection steps. It can 
   be used to fill histograms in three different granularity levels according to STANDARD 
   (<10 histograms), VERBOSE(<20 histograms), DEBUG(<30 histgorams). Note that for the sake 
   of simplicity and to force the analyst to keep the number of histograms to be monitored 
   small the MonitorEnsemble class contains the histograms for all objects at once. It should 
   not contain much more than 10 histograms though in the STANDARD configuration, as these 
   histograms will be monitored at each SelectionStep. Monitoring of histograms after selec-
   tion steps within the same object collection needs to be implemented within the Monitor-
   Ensemble. It will not be covered by the SelectionStep class.
*/

namespace TopDiLeptonOffline {

  class MonitorEnsemble {
  public:
    /// different verbosity levels
    enum Level{ STANDARD, VERBOSE, DEBUG };
    /// make clear which LorentzVector to use
    /// for jet, electrons and muon buffering
    typedef reco::LeafCandidate::LorentzVector LorentzVector;
    /// different decay channels
    enum DecayChannel{ NONE, DIMUON, DIELEC, ELECMU };
    
  public:
    /// default contructor
//    MonitorEnsemble(const char* label, const edm::ParameterSet& cfg);
    MonitorEnsemble(const char* label, const edm::ParameterSet& cfg, edm::ConsumesCollector && iC );
    /// default destructor
    ~MonitorEnsemble(){};
    
    /// book histograms in subdirectory _directory_
    void book(std::string directory);
    /// fill monitor histograms with electronId and jetCorrections
    void fill(const edm::Event& event, const edm::EventSetup& setup);

  private:
    /// deduce monitorPath from label, the label is expected
    /// to be of type 'selectionPath:monitorPath'
    std::string monitorPath(const std::string& label) const { return label.substr(label.find(':')+1); };  
    /// deduce selectionPath from label, the label is 
    /// expected to be of type 'selectionPath:monitorPath' 
    std::string selectionPath(const std::string& label) const { return label.substr(0, label.find(':')); };  
    /// determine dileptonic decay channel 
    DecayChannel decayChannel(const std::vector<const reco::PFCandidate*>& muons, const std::vector<const reco::GsfElectron*>& elecs) const;

    /// set labels for event logging histograms
    void loggerBinLabels(std::string hist);
    /// set configurable labels for trigger monitoring histograms
    void triggerBinLabels(std::string channel, const std::vector<std::string> labels);
    /// fill trigger monitoring histograms
    void fill(const edm::Event& event, const edm::TriggerResults& triggerTable, std::string channel, const std::vector<std::string> labels) const;

    /// check if histogram was booked
    bool booked(const std::string histName) const { return hists_.find(histName.c_str())!=hists_.end(); };
    /// fill histogram if it had been booked before
    void fill(const std::string histName, double value) const { if(booked(histName.c_str())) hists_.find(histName.c_str())->second->Fill(value); };
    /// fill histogram if it had been booked before (2-dim version)
    void fill(const std::string histName, double xValue, double yValue) const { if(booked(histName.c_str())) hists_.find(histName.c_str())->second->Fill(xValue, yValue); };
    /// fill histogram if it had been booked before (2-dim version)
    void fill(const std::string histName, double xValue, double yValue, double zValue) const { if(booked(histName.c_str())) hists_.find(histName.c_str())->second->Fill(xValue, yValue, zValue); };

  private:
    /// verbosity level for booking
    Level verbosity_;
    /// instance label 
    std::string label_;
    /// input sources for monitoring
    //edm::InputTag elecs_, muons_, jets_; 
    edm::EDGetTokenT<edm::View<reco::Jet> >  jets_; 
    edm::EDGetTokenT<edm::View<reco::PFCandidate> > muons_;
    edm::EDGetTokenT<edm::View<reco::GsfElectron> > elecs_;

    /// considers a vector of METs
    //std::vector<edm::InputTag> mets_;
    std::vector<edm::EDGetTokenT<edm::View<reco::MET> > > mets_;

    /// trigger table
    //edm::InputTag triggerTable_;
    edm::EDGetTokenT<edm::TriggerResults> triggerTable_;
    /// trigger paths for monitoring, expected 
    /// to be of form signalPath:MonitorPath
    std::vector<std::string> elecMuPaths_;
    /// trigger paths for di muon channel
    std::vector<std::string> diMuonPaths_;

    /// electronId label
    //edm::InputTag electronId_;
    edm::EDGetTokenT<edm::ValueMap<float> > electronId_;
    /// electronId pattern we expect the following pattern:
    ///  0: fails
    ///  1: passes electron ID only
    ///  2: passes electron Isolation only
    ///  3: passes electron ID and Isolation only
    ///  4: passes conversion rejection
    ///  5: passes conversion rejection and ID
    ///  6: passes conversion rejection and Isolation
    ///  7: passes the whole selection
    /// As described on https://twiki.cern.ch/twiki/bin/view/CMS/SimpleCutBasedEleID
    int eidPattern_;
    /// extra isolation criterion on electron
    StringCutObjectSelector<reco::GsfElectron>* elecIso_;
    /// extra selection on electrons
    StringCutObjectSelector<reco::GsfElectron>* elecSelect_;

    /// extra isolation criterion on muon
    StringCutObjectSelector<reco::PFCandidate, true>* muonIso_;
    
    /// extra selection on muons
    StringCutObjectSelector<reco::PFCandidate, true>* muonSelect_;

    /// jetCorrector
    std::string jetCorrector_;
    /// jetID as an extra selection type 
    //edm::InputTag jetIDLabel_;
    edm::EDGetTokenT<reco::JetIDValueMap> jetIDLabel_;
    /// extra jetID selection on calo jets
    StringCutObjectSelector<reco::JetID>* jetIDSelect_;
    /// extra selection on jets (here given as std::string as it depends
    /// on the the jet type, which selections are valid and which not)
    std::string jetSelect_;
    /// mass window upper and lower edge
    double lowerEdge_, upperEdge_;

    /// number of logged interesting events
    int elecMuLogged_, diMuonLogged_, diElecLogged_;
    /// storage manager
    DQMStore* store_;
    /// histogram container  
    std::map<std::string,MonitorElement*> hists_;
  };

  inline void 
  MonitorEnsemble::loggerBinLabels(std::string hist)
  {
    // set axes titles for selected events
    hists_[hist.c_str()]->getTH1()->SetOption("TEXT");
    hists_[hist.c_str()]->setBinLabel( 1 , "Run"             , 1);
    hists_[hist.c_str()]->setBinLabel( 2 , "Block"           , 1);
    hists_[hist.c_str()]->setBinLabel( 3 , "Event"           , 1);
    hists_[hist.c_str()]->setBinLabel( 6 , "pt_{L2L3}(jet1)" , 1);
    hists_[hist.c_str()]->setBinLabel( 7 , "pt_{L2L3}(jet2)" , 1);
    hists_[hist.c_str()]->setBinLabel( 8 , "MET_{Calo}"      , 1);
    hists_[hist.c_str()]->setAxisTitle("logged evts"         , 2);

    if(hist=="diMuonLogger_"){
      hists_[hist.c_str()]->setBinLabel( 4 , "pt(muon)" , 1);
      hists_[hist.c_str()]->setBinLabel( 5 , "pt(muon)" , 1);
    }
    if(hist=="diElecLogger_"){
      hists_[hist.c_str()]->setBinLabel( 4 , "pt(elec)" , 1);
      hists_[hist.c_str()]->setBinLabel( 5 , "pt(elec)" , 1);
    }
    if(hist=="elecMuLogger_"){
      hists_[hist.c_str()]->setBinLabel( 4 , "pt(elec)" , 1);
      hists_[hist.c_str()]->setBinLabel( 5 , "pt(muon)" , 1);
    }
  }

  inline void 
  MonitorEnsemble::triggerBinLabels(std::string channel, const std::vector<std::string> labels)
  {
    for(unsigned int idx=0; idx<labels.size(); ++idx){
      hists_[(channel+"Mon_").c_str()]->setBinLabel( idx+1, "["+monitorPath(labels[idx])+"]", 1);
      hists_[(channel+"Eff_").c_str()]->setBinLabel( idx+1, "["+selectionPath(labels[idx])+"]|["+monitorPath(labels[idx])+"]", 1);
    }
  }

  inline void 
  MonitorEnsemble::fill(const edm::Event& event, const edm::TriggerResults& triggerTable, std::string channel, const std::vector<std::string> labels) const
  {
    for(unsigned int idx=0; idx<labels.size(); ++idx){
      if( accept(event, triggerTable, monitorPath(labels[idx])) ){
	fill((channel+"Mon_").c_str(), idx+0.5 );
	// take care to fill triggerMon_ before evts is being called
	int evts = hists_.find((channel+"Mon_").c_str())->second->getBinContent(idx+1);
	double value = hists_.find((channel+"Eff_").c_str())->second->getBinContent(idx+1);
	fill((channel+"Eff_").c_str(), idx+0.5, 1./evts*(accept(event, triggerTable, selectionPath(labels[idx]))-value));
      }
    }
  }
  
  inline MonitorEnsemble::DecayChannel
  MonitorEnsemble::decayChannel(const std::vector<const reco::PFCandidate*>&  muons, const std::vector<const reco::GsfElectron*>& elecs) const   
  {
    DecayChannel type=NONE;
    if( muons.size()>1 ){ type=DIMUON; } else if( elecs.size()>1 ){ type=DIELEC; } else if( !elecs.empty() && !muons.empty() ){ type=ELECMU; }
    return type;
  } 
  
}

#include <utility>

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Common/interface/TriggerResults.h"
  
/**
   \class   TopDiLeptonOfflineDQM TopDiLeptonOfflineDQM.h "DQM/Physics/plugins/TopDiLeptonOfflineDQM.h"

   \brief   Module to apply a monitored selection of top like events in the di-leptonic channel

   Plugin to apply a monitored selection of top like events with some minimal flexibility 
   in the number and definition of the selection steps. To achieve this flexibility it 
   employes the SelectionStep class. The MonitorEnsemble class is used to provide a well 
   defined set of histograms to be monitored after each selection step. The SelectionStep 
   class provides a flexible and intuitive selection via the StringCutParser. SelectionStep 
   and MonitorEnsemble classes are interleaved. The monitoring starts after a preselection 
   step (which is not monitored in the context of this module) with an instance of the 
   MonitorEnsemble class. The following objects are supported for selection:

    - jets  : of type reco::Jet
    - elecs : of type reco::GsfElectron
    - muons : of type reco::PFCandidate
    - met   : of type reco::MET

   These types have to be present as prefix of the selection step paramter _label_ separated 
   from the rest of the label by a ':' (e.g. in the form "jets:step0"). The class expects 
   selection labels of this type. They will be disentangled by the private helper functions 
   _objectType_ and _seletionStep_ as declared below.
*/

/// define MonitorEnsembple to be used
//using TopDiLeptonOffline::MonitorEnsemble;

class TopDiLeptonOfflineDQM : public edm::EDAnalyzer  {
 public: 
  /// default constructor
  TopDiLeptonOfflineDQM(const edm::ParameterSet& cfg);
  /// default destructor
  ~TopDiLeptonOfflineDQM(){ 
    if( beamspotSelect_ ) delete beamspotSelect_; 
    if( vertexSelect_ ) delete vertexSelect_;

    if( MuonStep) delete MuonStep;

    if( ElectronStep) delete ElectronStep;

    if( PvStep) delete PvStep;

    if( METStep) delete METStep;

    for(unsigned int i = 0; i < JetSteps.size(); i++)

	if( JetSteps[i]) delete JetSteps[i];

   for(unsigned int i = 0; i < CaloJetSteps.size(); i++)

	if( CaloJetSteps[i]) delete CaloJetSteps[i];

  for(unsigned int i = 0; i < PFJetSteps.size(); i++)

	if( PFJetSteps[i]) delete PFJetSteps[i];
  }
  
  /// do this during the event loop
  virtual void analyze(const edm::Event& event, const edm::EventSetup& setup);
    
 private:
  /// deduce object type from ParameterSet label, the label
  /// is expected to be of type 'objectType:selectionStep'
  std::string objectType(const std::string& label) { return label.substr(0, label.find(':')); };  
  /// deduce selection step from ParameterSet label, the 
  /// label is expected to be of type 'objectType:selectionStep' 
  std::string selectionStep(const std::string& label) { return label.substr(label.find(':')+1); };  

 private:
  /// trigger table
  //edm::InputTag triggerTable_;
  edm::EDGetTokenT<edm::TriggerResults> triggerTable_;
  /// trigger paths
  std::vector<std::string> triggerPaths_;
  /// primary vertex 
  //edm::InputTag vertex_;
  edm::EDGetTokenT<std::vector<reco::Vertex> > vertex_;
  /// string cut selector
  StringCutObjectSelector<reco::Vertex>* vertexSelect_;
  /// beamspot 
  //edm::InputTag beamspot_;
  edm::EDGetTokenT<reco::BeamSpot> beamspot_;
  /// string cut selector
  StringCutObjectSelector<reco::BeamSpot>* beamspotSelect_;

  /// needed to guarantee the selection order as defined by the order of
  /// ParameterSets in the _selection_ vector as defined in the config
  std::vector<std::string> selectionOrder_;
  /// this is the heart component of the plugin; std::string keeps a label 
  /// the selection step for later identification, edm::ParameterSet keeps
  /// the configuration of the selection for the SelectionStep class, 
  /// MonitoringEnsemble keeps an instance of the MonitorEnsemble class to 
  /// be filled _after_ each selection step
  std::map<std::string, std::pair<edm::ParameterSet, TopDiLeptonOffline::MonitorEnsemble*> > selection_;
  SelectionStep<reco::PFCandidate> * MuonStep;
  SelectionStep<reco::GsfElectron> * ElectronStep;
  SelectionStep<reco::Vertex> * PvStep;
  SelectionStep<reco::MET> * METStep;
  std::vector<SelectionStep<reco::Jet> * > JetSteps;
  std::vector<SelectionStep<reco::CaloJet> * > CaloJetSteps;
  std::vector<SelectionStep<reco::PFJet> * > PFJetSteps;
};

#endif

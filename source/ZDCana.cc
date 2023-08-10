//____________________________________________________________________________..
//
//
//____________________________________________________________________________..

#include "ZDCana.h"

// Fun4all stuff
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllHistoManager.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>
#include <ffaobjects/EventHeader.h>

#include <phool/PHCompositeNode.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4Shower.h>
#include <g4main/PHG4HitDefs.h>

// ROOT stuff
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TTree.h>
#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <HepMC/GenEvent.h>
#include <HepMC/GenVertex.h> // for GenVertex, GenVertex::part...
#pragma GCC diagnostic pop

#include <HepMC/GenParticle.h> // for GenParticle
#include <HepMC/GenRanges.h>
#include <HepMC/HeavyIon.h>      // for HeavyIon
#include <HepMC/IteratorRange.h> // for children, descendants
#include <HepMC/SimpleVector.h>

//____________________________________________________________________________..
ZDCana::ZDCana(const std::string &name = "ZDCana", const std::string &outName = "ZDCresult") : SubsysReco(name), G4Hits(nullptr), HEPMCinfo(nullptr), m_ZDC_Hit_j(), m_ZDC_Hit_k(), m_ZDC_Hit_Edep(), m_ZDC_Hit_Evis(), m_SMD_Hit_Evis(), m_SMD_Hit_Edep(), m_SMD_Hit_j(), m_SMD_Hit_k(), m_SMD_Hit_x(), m_SMD_Hit_y(), m_SMD_Hit_z(), m_BBC_Hit_Edep(), m_BBC_Hit_Evis(), m_e(), m_px(), m_py(), m_pz(), m_pid(), m_bphi(), m_b()

                                                                                               ,
                                                                                               Outfile(outName)
{
  std::cout << "ZDCana::ZDCana(const std::string &name) Calling ctor" << std::endl;
}

//____________________________________________________________________________..
ZDCana::~ZDCana()
{
  std::cout << "ZDCana::~ZDCana() Calling dtor" << std::endl;
}

//____________________________________________________________________________..
int ZDCana::Init(PHCompositeNode *topNode)
{

  std::cout << "ZDCana::Init(PHCompositeNode *topNode) Initializing" << std::endl;
  out = new TFile(Outfile.c_str(), "RECREATE");
  G4Hits = new TTree("TreeG4Hit", "Tree for G4Hits");
  G4Hits->Branch("ZDC_Hit_j", &m_ZDC_Hit_j);
  G4Hits->Branch("ZDC_Hit_k", &m_ZDC_Hit_k);
  G4Hits->Branch("ZDC_Hit_Edep", &m_ZDC_Hit_Edep);
  G4Hits->Branch("ZDC_Hit_Evis", &m_ZDC_Hit_Evis);
  G4Hits->Branch("ZDC_primary_px", &m_ZDC_primary_px);
  G4Hits->Branch("ZDC_primary_py", &m_ZDC_primary_py);
  G4Hits->Branch("ZDC_primary_pz", &m_ZDC_primary_pz);
  G4Hits->Branch("ZDC_primary_pid", &m_ZDC_primary_pid);
  G4Hits->Branch("ZDC_showerid", &m_ZDC_showerid);
  G4Hits->Branch("SMD_Hit_Evis", &m_SMD_Hit_Evis);
  G4Hits->Branch("SMD_Hit_Edep", &m_SMD_Hit_Edep);
  G4Hits->Branch("SMD_Hit_j", &m_SMD_Hit_j);
  G4Hits->Branch("SMD_Hit_k", &m_SMD_Hit_k);
  G4Hits->Branch("SMD_Hit_x", &m_SMD_Hit_x);
  G4Hits->Branch("SMD_Hit_y", &m_SMD_Hit_y);
  G4Hits->Branch("SMD_Hit_z", &m_SMD_Hit_z);
  G4Hits->Branch("SMD_showerid", &m_SMD_showerid);
  G4Hits->Branch("BBC_Hit_Edep", &m_BBC_Hit_Edep);
  G4Hits->Branch("BBC_Hit_Evis", &m_BBC_Hit_Evis);

  HEPMCinfo = new TTree("TreeHEPMC", "Tree of HEPMC particles");
  HEPMCinfo->Branch("e", &m_e);
  HEPMCinfo->Branch("px", &m_px);
  HEPMCinfo->Branch("py", &m_py);
  HEPMCinfo->Branch("pz", &m_pz);
  HEPMCinfo->Branch("pid", &m_pid);
  HEPMCinfo->Branch("bphi", &m_bphi);
  HEPMCinfo->Branch("b", &m_b);

  G4Truth = new TTree("TreeG4Truth", "Tree of G4Truth particles and shower info");
  G4Truth->Branch("truth_px", &m_truth_px);
  G4Truth->Branch("truth_py", &m_truth_py);
  G4Truth->Branch("truth_pz", &m_truth_pz);
  G4Truth->Branch("truth_e", &m_truth_e);
  G4Truth->Branch("truth_pid", &m_truth_pid);
  G4Truth->Branch("truth_showerE", &m_truth_showerE);
  G4Truth->Branch("showerid", &m_showerid);
  G4Truth->Branch("e00", &m_e00);
  G4Truth->Branch("e01", &m_e01);
  G4Truth->Branch("e02", &m_e02);
  G4Truth->Branch("e10", &m_e10);
  G4Truth->Branch("e11", &m_e11);
  G4Truth->Branch("e12", &m_e12);
  G4Truth->Branch("particle_px", &m_particle_px);
  G4Truth->Branch("particle_py", &m_particle_py);
  G4Truth->Branch("particle_pz", &m_particle_pz);
  G4Truth->Branch("particle_e", &m_particle_e);
  G4Truth->Branch("particle_pid", &m_particle_pid);


  Fun4AllServer *se = Fun4AllServer::instance();
  se->Print("NODETREE");
  hm = new Fun4AllHistoManager("MYHISTOS");

  se->registerHistoManager(hm);

  se->registerHisto(G4Hits->GetName(), G4Hits);
  se->registerHisto(HEPMCinfo->GetName(), HEPMCinfo);

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int ZDCana::InitRun(PHCompositeNode *topNode)
{
  std::cout << "ZDCana::InitRun(PHCompositeNode *topNode) Initializing for Run XXX" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int ZDCana::process_event(PHCompositeNode *topNode)
{
  std::cout << "ZDCana::process_event(PHCompositeNode *topNode) Processing Event" << std::endl;
  // HEPMC stuff
  if (m_HepMCFlag)
  {
    PHHepMCGenEventMap *genevtmap = findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");
    for (PHHepMCGenEventMap::Iter iter = genevtmap->begin(); iter != genevtmap->end(); ++iter)
    {
      PHHepMCGenEvent *genevt = iter->second;
      HepMC::GenEvent *event = genevt->getEvent();
      if (!event)
      {
        std::cout << PHWHERE << " no evt pointer under HEPMC Node found" << std::endl;
      }
      else
      {
        HepMC::HeavyIon *hi = event->heavy_ion();
        if (!hi)
        {
          std::cout << PHWHERE << ": Fermi Motion Afterburner needs the Heavy Ion Event Info, GenEvent::heavy_ion() returns NULL" << std::endl;
          exit(1);
        }
        float b = hi->impact_parameter();
        float bphi = hi->event_plane_angle();
        m_b.push_back(b);
        m_bphi.push_back(bphi);
        for (HepMC::GenEvent::particle_const_iterator p = event->particles_begin(), prev = event->particles_end(); p != event->particles_end(); prev = p, ++p)
        {
          int id = (*p)->pdg_id();
          HepMC::GenParticle *n = (*p);
          m_px.push_back(n->momentum().px());
          m_py.push_back(n->momentum().py());
          m_pz.push_back(n->momentum().pz());
          m_e.push_back(n->momentum().e());
          m_pid.push_back(id);
        }
      }
    }
  }
  // truth container here
  PHG4TruthInfoContainer *truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");

  // get primary particle range
  PHG4TruthInfoContainer::Range range = truthinfo->GetPrimaryParticleRange();
  // loop over primary particles
  for (PHG4TruthInfoContainer::ConstIterator iter = range.first; iter != range.second; ++iter)
  {
    // get primary particle
    PHG4Particle *primary = iter->second;
    // get momentum
    TLorentzVector mom(primary->get_px(), primary->get_py(), primary->get_pz(), primary->get_e());
    // get pid
    int pid = primary->get_pid();
    // print all info
    m_particle_px.push_back(mom.Px());
    m_particle_py.push_back(mom.Py());
    m_particle_pz.push_back(mom.Pz());
    m_particle_e.push_back(mom.E());
    m_particle_pid.push_back(pid);

  }
  //----------------------------------------------------------------------------------------

  PHG4HitContainer *hits_ZDC = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_ZDC");
  if (!hits_ZDC)
    std::cout << "ZDCana::process_event(PHCompositeNode *topNode) No ZDC hits found" << std::endl;
  // PHG4HitContainer *hits_BBC = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_BBC");





  // get shower
  PHG4TruthInfoContainer::ShowerRange showerrange = truthinfo->GetPrimaryShowerRange();
  // loop over showers
  for (PHG4TruthInfoContainer::ShowerIterator showeriter = showerrange.first; showeriter != showerrange.second; ++showeriter)
  {
    float ZDC_E[2][3] = {{0}};
    PHG4Shower *shower = showeriter->second;
    // get particle
    PHG4Particle *particle = truthinfo->GetParticle(shower->get_parent_particle_id());
    // get momentum
    TLorentzVector mom(particle->get_px(), particle->get_py(), particle->get_pz(), particle->get_e());
    // get pid
    int pid = particle->get_pid();
    // get light yield
    float totalE = shower->get_edep();
    // get shower id
    int showerid = shower->get_id();

    // Get the g4hit_ids map from the shower object
    // const PHG4Shower::HitIdMap &hit_ids = shower->g4hit_ids();

    PHG4Shower::HitIdIter shower_hit_iter;
    PHG4Shower::HitIdIter shower_hit_begin = shower->begin_g4hit_id();
    PHG4Shower::HitIdIter shower_hit_end = shower->end_g4hit_id();

    for (shower_hit_iter = shower_hit_begin; shower_hit_iter != shower_hit_end; ++shower_hit_iter)
    {
      // int volume_id = shower_hit_iter->first;                               // the volume ID
      std::set<PHG4HitDefs::keytype> &hit_id_set = shower_hit_iter->second; // the set of hit IDs

      for (PHG4HitDefs::keytype hit_id : hit_id_set)
      {
        PHG4Hit *hit = hits_ZDC->findHit(hit_id);

        if (hit)
        {
          float light_yield = hit->get_light_yield();
          float edep = hit->get_edep();
          int idx_j = hit->get_index_j();
          int idx_k = hit->get_index_k();
          float x = hit->get_x(0);
          float y = hit->get_y(0);
          float z = hit->get_z(0);
          int showerid = hit->get_shower_id();

          if (idx_j < 2)
          {

            ZDC_E[idx_j][idx_k] += light_yield;
            m_ZDC_Hit_Evis.push_back(light_yield);
            m_ZDC_Hit_Edep.push_back(edep);
            m_ZDC_Hit_j.push_back(idx_j);
            m_ZDC_Hit_k.push_back(idx_k);
            m_ZDC_primary_px.push_back(particle->get_px());
            m_ZDC_primary_py.push_back(particle->get_py());
            m_ZDC_primary_pz.push_back(particle->get_pz());
            m_ZDC_primary_pid.push_back(particle->get_pid());
            m_ZDC_showerid.push_back(showerid);
          }
          else
          {
            m_SMD_Hit_Evis.push_back(light_yield);
            m_SMD_Hit_Edep.push_back(edep);
            m_SMD_Hit_j.push_back(idx_j);
            m_SMD_Hit_k.push_back(idx_k);
            m_SMD_Hit_x.push_back(x);
            m_SMD_Hit_y.push_back(y);
            m_SMD_Hit_z.push_back(z);
            m_SMD_showerid.push_back(showerid);
          }
         
        }
        else
        {
          std::cout << "ZDCana::process_event(PHCompositeNode *topNode) No ZDC hits found" << std::endl;
        }
      }
    }

    // print all info
    if (totalE > 0)
    {
      m_truth_px.push_back(mom.Px());
      m_truth_py.push_back(mom.Py());
      m_truth_pz.push_back(mom.Pz());
      m_truth_e.push_back(mom.E());
      m_truth_pid.push_back(pid);
      m_truth_showerE.push_back(totalE);
      m_showerid.push_back(showerid);
      m_e00.push_back(ZDC_E[0][0]);
      m_e01.push_back(ZDC_E[0][1]);
      m_e02.push_back(ZDC_E[0][2]);
      m_e10.push_back(ZDC_E[1][0]);
      m_e11.push_back(ZDC_E[1][1]);
      m_e12.push_back(ZDC_E[1][2]);

      //print it all
      //std::cout<<"shower id: "<<showerid<<" pid: "<<pid<<" px: "<<mom.Px()<<" py: "<<mom.Py()<<" pz: "<<mom.Pz()<<" E: "<<mom.E()<<" totalE: "<<totalE<<std::endl;
      //std::cout<<"E00:" <<ZDC_E[0][0]<<" E01: "<<ZDC_E[0][1]<<" E02: "<<ZDC_E[0][2]<<" E10: "<<ZDC_E[1][0]<<" E11: "<<ZDC_E[1][1]<<" E12: "<<ZDC_E[1][2]<<std::endl;
    }
  }
  /*
    if (hits_ZDC)
    {

      PHG4HitContainer::ConstRange hit_range = hits_ZDC->getHits();
      for (PHG4HitContainer::ConstIterator hit_iter = hit_range.first; hit_iter != hit_range.second; hit_iter++)
      {
        float light_yield = hit_iter->second->get_light_yield();
        float edep = hit_iter->second->get_edep();
        int idx_j = hit_iter->second->get_index_j();
        int idx_k = hit_iter->second->get_index_k();
        float x = hit_iter->second->get_x(0);
        float y = hit_iter->second->get_y(0);
        float z = hit_iter->second->get_z(0);
        int showerid = hit_iter->second->get_shower_id();
        if (idx_j < 2)
        {
          m_ZDC_Hit_Evis.push_back(light_yield);
          m_ZDC_Hit_Edep.push_back(edep);
          m_ZDC_Hit_j.push_back(idx_j);
          m_ZDC_Hit_k.push_back(idx_k);
          PHG4Shower *shower = truthinfo->GetShower(showerid);
          PHG4Particle *particle = truthinfo->GetParticle(shower->get_parent_particle_id());
          m_ZDC_primary_px.push_back(particle->get_px());
          m_ZDC_primary_py.push_back(particle->get_py());
          m_ZDC_primary_pz.push_back(particle->get_pz());
          m_ZDC_primary_pid.push_back(particle->get_pid());
        }
        else
        {
          m_SMD_Hit_Evis.push_back(light_yield);
          m_SMD_Hit_Edep.push_back(edep);
          m_SMD_Hit_j.push_back(idx_j);
          m_SMD_Hit_k.push_back(idx_k);
          m_SMD_Hit_x.push_back(x);
          m_SMD_Hit_y.push_back(y);
          m_SMD_Hit_z.push_back(z);
        }
        if (light_yield > 10 || light_yield < 0)
          std::cout << idx_j << " " << idx_j << " " << light_yield << " " << edep << std::endl;
      }
    }
    if (hits_BBC)
    {

      PHG4HitContainer::ConstRange hit_range = hits_BBC->getHits();
      for (PHG4HitContainer::ConstIterator hit_iter = hit_range.first; hit_iter != hit_range.second; hit_iter++)
      {
        float light_yield = hit_iter->second->get_light_yield();
        float edep = hit_iter->second->get_edep();

        m_BBC_Hit_Evis.push_back(light_yield);
        m_BBC_Hit_Edep.push_back(edep);
      }
    }
  */
  HEPMCinfo->Fill();
  G4Hits->Fill();
  G4Truth->Fill();
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int ZDCana::ResetEvent(PHCompositeNode *topNode)
{
  std::cout << "ZDCana::ResetEvent(PHCompositeNode *topNode) Resetting internal structures, prepare for next event" << std::endl;

  m_ZDC_Hit_j.clear();
  m_ZDC_Hit_k.clear();
  m_ZDC_Hit_Edep.clear();
  m_ZDC_Hit_Evis.clear();
  m_ZDC_primary_px.clear();
  m_ZDC_primary_py.clear();
  m_ZDC_primary_pz.clear();
  m_ZDC_primary_pid.clear();
  m_ZDC_showerid.clear();
  m_SMD_Hit_Evis.clear();
  m_SMD_Hit_Edep.clear();
  m_SMD_Hit_j.clear();
  m_SMD_Hit_k.clear();
  m_SMD_Hit_x.clear();
  m_SMD_Hit_y.clear();
  m_SMD_Hit_z.clear();
  m_SMD_showerid.clear();
  m_BBC_Hit_Edep.clear();
  m_BBC_Hit_Evis.clear();
  m_e.clear();
  m_px.clear();
  m_py.clear();
  m_pz.clear();
  m_pid.clear();
  m_bphi.clear();
  m_b.clear();
  m_truth_px.clear();
  m_truth_py.clear();
  m_truth_pz.clear();
  m_truth_e.clear();
  m_truth_pid.clear();
  m_truth_showerE.clear();
  m_showerid.clear();
  m_e00.clear();
  m_e01.clear();
  m_e02.clear();
  m_e10.clear();
  m_e11.clear();
  m_e12.clear();
  m_particle_px.clear();
  m_particle_py.clear();
  m_particle_pz.clear();
  m_particle_e.clear();
  m_particle_pid.clear();

  

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int ZDCana::EndRun(const int runnumber)
{
  std::cout << "ZDCana::EndRun(const int runnumber) Ending Run for Run " << runnumber << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int ZDCana::End(PHCompositeNode *topNode)
{
  std::cout << "ZDCana::End(PHCompositeNode *topNode) This is the End..." << std::endl;

  out->cd();

  G4Hits->Write();
  HEPMCinfo->Write();
  G4Truth->Write();
  out->Close();
  delete out;
  hm->dumpHistos(Outfile.c_str(), "UPDATE");
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int ZDCana::Reset(PHCompositeNode *topNode)
{
  std::cout << "ZDCana::Reset(PHCompositeNode *topNode) being Reset" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void ZDCana::Print(const std::string &what) const
{
  std::cout << "ZDCana::Print(const std::string &what) const Printing info for " << what << std::endl;
}

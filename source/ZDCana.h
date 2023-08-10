// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef ZDCANA_H
#define ZDCANA_H

#include <fun4all/SubsysReco.h>
//ROOT stuff
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TTree.h>


#include <string>

namespace HepMC
{
class GenEvent;
}


class PHCompositeNode;
class Fun4AllHistoManager;

class ZDCana : public SubsysReco
{
public:

  ZDCana(const std::string &name , const std::string &outName );

  ~ZDCana() override;

  /** Called during initialization.
      Typically this is where you can book histograms, and e.g.
      register them to Fun4AllServer (so they can be output to file
      using Fun4AllServer::dumpHistos() method).
   */
  int Init(PHCompositeNode *topNode) override;

  /** Called for first event when run number is known.
      Typically this is where you may want to fetch data from
      database, because you know the run number. A place
      to book histograms which have to know the run number.
   */
  int InitRun(PHCompositeNode *topNode) override;

  /** Called for each event.
      This is where you do the real work.
   */
  int process_event(PHCompositeNode *topNode) override;

  /// Clean up internals after each event.
  int ResetEvent(PHCompositeNode *topNode) override;

  /// Called at the end of each run.
  int EndRun(const int runnumber) override;

  /// Called at the end of all processing.
  int End(PHCompositeNode *topNode) override;

  /// Reset
  int Reset(PHCompositeNode * /*topNode*/) override;

  void Print(const std::string &what = "ALL") const override;

private:
  TTree *G4Hits;
  TTree *HEPMCinfo;
  TTree *G4Truth;
  


  std::vector<int> m_ZDC_Hit_j;
  std::vector<int> m_ZDC_Hit_k;
  std::vector<float> m_ZDC_Hit_Edep;
  std::vector<float> m_ZDC_Hit_Evis;
  std::vector<float> m_ZDC_primary_px;
  std::vector<float> m_ZDC_primary_py;
  std::vector<float> m_ZDC_primary_pz;
  std::vector<int> m_ZDC_primary_pid;
  std::vector<int> m_ZDC_showerid;
  std::vector<float> m_SMD_Hit_Evis;
  std::vector<float> m_SMD_Hit_Edep;
  std::vector<int> m_SMD_Hit_j;
  std::vector<int> m_SMD_Hit_k;
  std::vector<float> m_SMD_Hit_x;
  std::vector<float> m_SMD_Hit_y;
  std::vector<float> m_SMD_Hit_z;
  std::vector<int> m_SMD_showerid;
  std::vector<float> m_BBC_Hit_Edep;
  std::vector<float> m_BBC_Hit_Evis;
  
  

  std::vector<float> m_e;
  std::vector<float> m_px;
  std::vector<float> m_py;
  std::vector<float> m_pz;
  std::vector<float> m_pid;
  std::vector<float> m_bphi;
  std::vector<float> m_b;

  std::vector<float> m_truth_px;
  std::vector<float> m_truth_py;
  std::vector<float> m_truth_pz;
  std::vector<float> m_truth_e;
  std::vector<float> m_truth_pid;
  std::vector<float> m_truth_showerE;
  std::vector<float> m_showerid;
  std::vector<float> m_e00;
  std::vector<float> m_e01;
  std::vector<float> m_e02;
  std::vector<float> m_e10;
  std::vector<float> m_e11;
  std::vector<float> m_e12; 
  std::vector<float> m_particle_px;
  std::vector<float> m_particle_py;
  std::vector<float> m_particle_pz;
  std::vector<float> m_particle_e;
  std::vector<float> m_particle_pid;

  



  std::string Outfile;
  TFile *out;
  Fun4AllHistoManager *hm = nullptr;
  bool m_HepMCFlag = true;
};

#endif // ZDCANA_H

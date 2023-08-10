// read TTree and calculate the energy response at different location
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TMath.h"
#include "TString.h"
#include "TTree.h"
#include "TVector3.h"

void ZDCEScan()
{
  // read TTree TreeHEPMC and TreeG4Truth
  TChain *TreeHEPMC = new TChain("TreeHEPMC");
  TChain *TreeG4Truth = new TChain("TreeG4Truth");

  TreeHEPMC->Add("/sphenix/user/shuhangli/macros/detectors/sPHENIX/rootfiles/cross/zdcall.root");
  TreeG4Truth->Add("/sphenix/user/shuhangli/macros/detectors/sPHENIX/rootfiles/cross/zdcall.root");

  // get branch

  std::vector<float> *m_e = 0;
  std::vector<float> *m_px = 0;
  std::vector<float> *m_py = 0;
  std::vector<float> *m_pz = 0;
  std::vector<int> *m_pid = 0;
  std::vector<float> *m_bphi = 0;
  std::vector<float> *m_b = 0;

  std::vector<float> *m_truth_px = 0;
  std::vector<float> *m_truth_py = 0;
  std::vector<float> *m_truth_pz = 0;
  std::vector<float> *m_truth_e = 0;
  std::vector<int> *m_truth_pid = 0;
  std::vector<float> *m_truth_showerE = 0;
  std::vector<int> *m_showerid = 0;
  std::vector<float> *m_e00 = 0;
  std::vector<float> *m_e01 = 0;
  std::vector<float> *m_e02 = 0;
  std::vector<float> *m_e10 = 0;
  std::vector<float> *m_e11 = 0;
  std::vector<float> *m_e12 = 0;
  std::vector<float> *m_particle_px = 0;
  std::vector<float> *m_particle_py = 0;
  std::vector<float> *m_particle_pz = 0;
  std::vector<float> *m_particle_e = 0;
  std::vector<float> *m_particle_pid = 0;

  TreeHEPMC->SetBranchAddress("e", &m_e);
  TreeHEPMC->SetBranchAddress("px", &m_px);
  TreeHEPMC->SetBranchAddress("py", &m_py);
  TreeHEPMC->SetBranchAddress("pz", &m_pz);
  TreeHEPMC->SetBranchAddress("pid", &m_pid);
  TreeHEPMC->SetBranchAddress("bphi", &m_bphi);
  TreeHEPMC->SetBranchAddress("b", &m_b);

  TreeG4Truth->SetBranchAddress("truth_px", &m_truth_px);
  TreeG4Truth->SetBranchAddress("truth_py", &m_truth_py);
  TreeG4Truth->SetBranchAddress("truth_pz", &m_truth_pz);
  TreeG4Truth->SetBranchAddress("truth_e", &m_truth_e);
  TreeG4Truth->SetBranchAddress("truth_pid", &m_truth_pid);
  TreeG4Truth->SetBranchAddress("truth_showerE", &m_truth_showerE);
  TreeG4Truth->SetBranchAddress("showerid", &m_showerid);
  TreeG4Truth->SetBranchAddress("e00", &m_e00);
  TreeG4Truth->SetBranchAddress("e01", &m_e01);
  TreeG4Truth->SetBranchAddress("e02", &m_e02);
  TreeG4Truth->SetBranchAddress("e10", &m_e10);
  TreeG4Truth->SetBranchAddress("e11", &m_e11);
  TreeG4Truth->SetBranchAddress("e12", &m_e12);
  TreeG4Truth->SetBranchAddress("particle_px", &m_particle_px);
  TreeG4Truth->SetBranchAddress("particle_py", &m_particle_py);
  TreeG4Truth->SetBranchAddress("particle_pz", &m_particle_pz);
  TreeG4Truth->SetBranchAddress("particle_e", &m_particle_e);
  TreeG4Truth->SetBranchAddress("particle_pid", &m_particle_pid);

  // create output file
  TFile *out = new TFile("ZDCEScancrossnew.root", "RECREATE");

  // create histograms
  TH2F *hitmap = new TH2F("hitmap", "hitmap", 100, -10, 10, 100, -10, 10);

  // make x y angle bins
  const int nxbins = 13;
  const int nybins = 13;
  TH1F *htotalE[nxbins][nybins];
  TH1F *htowerE[nxbins][nybins][3];

  for (int i = 0; i < nxbins; i++)
  {
    for (int j = 0; j < nybins; j++)
    {
      TString hname = Form("htotalE_%d_%d", i, j);
      TString htitle = Form("htotalE_%d_%d", i, j);
      htotalE[i][j] = new TH1F(hname, htitle, 200, 0, 2000);
      for (int k = 0; k < 3; k++)
      {
        TString hname2 = Form("htowerE_%d_%d_%d", i, j, k);
        TString htitle2 = Form("htowerE_%d_%d_%d", i, j, k);
        htowerE[i][j][k] = new TH1F(hname2, htitle2, 150, 0, 1500);
        // set title
        htowerE[i][j][k]->GetXaxis()->SetTitle("E (PE)");
        htowerE[i][j][k]->GetYaxis()->SetTitle("Arb. Units");
      }
    }
  }
  TH2F *h2bins = new TH2F("h2bins", "h2bins", nxbins, -3.25, 3.25, nybins, -3.25, 3.25);
  TH2F *h2resolution = new TH2F("h2resolution", "ZDC resolution", nxbins, -3.25, 3.25, nybins, -3.25, 3.25);
  h2resolution->GetXaxis()->SetTitle("x angle (mrad)");
  h2resolution->GetYaxis()->SetTitle("y angle (mrad)");
  TH2F *h2response = new TH2F("h2response", "ZDC response", nxbins, -3.25, 3.25, nybins, -3.25, 3.25);
  h2response->GetXaxis()->SetTitle("x angle (mrad)");
  h2response->GetYaxis()->SetTitle("y angle (mrad)");

  // total energy vs b
  TProfile *hpEvsb = new TProfile("hpEvsb", "E vs b", 100, 0, 20);
  hpEvsb->GetXaxis()->SetTitle("b (fm)");

  TProfile *hpEvsbt0 = new TProfile("hpEvsbt0", "E vs b for tower 0", 100, 0, 20);
  hpEvsbt0->GetXaxis()->SetTitle("b (fm)");
  TProfile *hpEvsbt1 = new TProfile("hpEvsbt1", "E vs b for tower 1", 100, 0, 20);
  hpEvsbt1->GetXaxis()->SetTitle("b (fm)");
  TProfile *hpEvsbt2 = new TProfile("hpEvsbt2", "E vs b for tower 2", 100, 0, 20);
  hpEvsbt2->GetXaxis()->SetTitle("b (fm)");

  TH2F *hneutron = new TH2F("hneutron", "neutron", 100, -5, 5, 100, -5, 5);
  hneutron->GetXaxis()->SetTitle("x (mrad)");
  hneutron->GetYaxis()->SetTitle("y (mrad)");

  TH1F *hangle = new TH1F("hangle", "hangle", 100, 0, 15);
  hangle->GetXaxis()->SetTitle("angle (mrad)");
  hangle->GetYaxis()->SetTitle("Arb. Units");

  TH1F *hZDCE_n = new TH1F("hZDCE_n", "hZDCE_n", 900, 0, 90000);
  hZDCE_n->GetXaxis()->SetTitle("ZDC Energy (PE)");
  hZDCE_n->GetYaxis()->SetTitle("Arb. Units");
  TH1F *hZDCE0_n = new TH1F("hZDCE0_n", "hZDCE0_n", 900, 0, 90000);
  hZDCE0_n->GetXaxis()->SetTitle("ZDC Energy (PE)");
  hZDCE0_n->GetYaxis()->SetTitle("Arb. Units");
  TH1F *hZDCE1_n = new TH1F("hZDCE1_n", "hZDCE1_n", 900, 0, 90000);
  hZDCE1_n->GetXaxis()->SetTitle("ZDC Energy (PE)");
  hZDCE1_n->GetYaxis()->SetTitle("Arb. Units");
  TH1F *hZDCE2_n = new TH1F("hZDCE2_n", "hZDCE2_n", 900, 0, 90000);
  hZDCE2_n->GetXaxis()->SetTitle("ZDC Energy (PE)");
  hZDCE2_n->GetYaxis()->SetTitle("Arb. Units");

  TH1F *hZDCE_s = new TH1F("hZDCE_s", "hZDCE_s", 900, 0, 90000);
  hZDCE_s->GetXaxis()->SetTitle("ZDC Energy (PE)");
  hZDCE_s->GetYaxis()->SetTitle("Arb. Units");
  TH1F *hZDCE0_s = new TH1F("hZDCE0_s", "hZDCE0_s", 900, 0, 90000);
  hZDCE0_s->GetXaxis()->SetTitle("ZDC Energy (PE)");
  hZDCE0_s->GetYaxis()->SetTitle("Arb. Units");
  TH1F *hZDCE1_s = new TH1F("hZDCE1_s", "hZDCE1_s", 900, 0, 90000);
  hZDCE1_s->GetXaxis()->SetTitle("ZDC Energy (PE)");
  hZDCE1_s->GetYaxis()->SetTitle("Arb. Units");
  TH1F *hZDCE2_s = new TH1F("hZDCE2_s", "hZDCE2_s", 900, 0, 90000);
  hZDCE2_s->GetXaxis()->SetTitle("ZDC Energy (PE)");
  hZDCE2_s->GetYaxis()->SetTitle("Arb. Units");

  // loop over events
  for (int i = 0; i < TreeHEPMC->GetEntries(); i++)
  {
    // print every 1000 events
    if (i % 1000 == 0) std::cout << "Event " << i << std::endl;
    TreeHEPMC->GetEntry(i);
    TreeG4Truth->GetEntry(i);

    float b = m_b->at(0);
    float totalZDCE = 0;
    float totalZDCE0 = 0;
    float totalZDCE1 = 0;
    float totalZDCE2 = 0;
    float totalZDCE_s = 0;
    float totalZDCE0_s = 0;
    float totalZDCE1_s = 0;
    float totalZDCE2_s = 0;
    // loop over truth showers:
    for (int j = 0; j < m_truth_px->size(); j++)
    {
      float px = m_truth_px->at(j);
      float py = m_truth_py->at(j);
      float pz = m_truth_pz->at(j);
      float e = m_truth_e->at(j);
      float pid = m_truth_pid->at(j);
      float edep = m_truth_showerE->at(j);
      // if not nuetron, skip
      if (pid != 2112) continue;

      if (abs(pz) < 99) continue;
      // if (edep < 100) continue;
      if (pz > 0)
      {
        float totalE = m_e00->at(j) + m_e01->at(j) + m_e02->at(j);
        totalZDCE += totalE;
        totalZDCE0 += m_e00->at(j);
        totalZDCE1 += m_e01->at(j);
        totalZDCE2 += m_e02->at(j);

        // calculate angle in mrad
        float xa = px / pz * 1000;
        float ya = py / pz * 1000;
        hitmap->Fill(xa, ya);
        // find bin in hbins
        h2bins->Fill(xa, ya);
        int xbin = h2bins->GetXaxis()->FindBin(xa) - 1;
        int ybin = h2bins->GetYaxis()->FindBin(ya) - 1;
        if (xbin >= nxbins || ybin >= nybins) continue;
        if (xbin < 0 || ybin < 0) continue;
        htotalE[xbin][ybin]->Fill(totalE);
        htowerE[xbin][ybin][0]->Fill(m_e00->at(j));
        htowerE[xbin][ybin][1]->Fill(m_e01->at(j));
        htowerE[xbin][ybin][2]->Fill(m_e02->at(j));
      }
      else
      {
        float totalE = m_e10->at(j) + m_e11->at(j) + m_e12->at(j);
        totalZDCE_s += totalE;
        totalZDCE0_s += m_e10->at(j);
        totalZDCE1_s += m_e11->at(j);
        totalZDCE2_s += m_e12->at(j);
      }
    }
    hpEvsb->Fill(b, totalZDCE);
    hpEvsbt0->Fill(b, totalZDCE0);
    hpEvsbt1->Fill(b, totalZDCE1);
    hpEvsbt2->Fill(b, totalZDCE2);
    hZDCE_n->Fill(totalZDCE);
    hZDCE0_n->Fill(totalZDCE0);
    hZDCE1_n->Fill(totalZDCE1);
    hZDCE2_n->Fill(totalZDCE2);
    hZDCE_s->Fill(totalZDCE_s);
    hZDCE0_s->Fill(totalZDCE0_s);
    hZDCE1_s->Fill(totalZDCE1_s);
    hZDCE2_s->Fill(totalZDCE2_s);
    

    // loop over truth particles
    for (int j = 0; j < m_particle_px->size(); j++)
    {
      float px = m_particle_px->at(j);
      float py = m_particle_py->at(j);
      float pz = m_particle_pz->at(j);
      float e = m_particle_e->at(j);
      float pid = m_particle_pid->at(j);
      // if not nuetron, skip
      if (pid != 2112) continue;
      if (pz < 0) continue;
      if (abs(pz) < 99) continue;
      // calculate angle in mrad
      float xa = px / pz * 1000;
      float ya = py / pz * 1000;
      float angle = sqrt(xa * xa + ya * ya);
      hangle->Fill(angle);
      hneutron->Fill(xa, ya);
    }
  }
  // loop over energy histogram and calculate energy resolution
  for (int i = 0; i < nxbins; i++)
  {
    for (int j = 0; j < nybins; j++)
    {
      float mean = htotalE[i][j]->GetMean();
      float rms = htotalE[i][j]->GetRMS();
      float resolution = rms / mean;
      float meanerror = htotalE[i][j]->GetMeanError();
      float rmserror = htotalE[i][j]->GetRMSError();
      // calculate error on resolution
      float relativeerror = sqrt(meanerror * meanerror / mean / mean + rmserror * rmserror / rms / rms);
      float error = resolution * relativeerror;
      h2resolution->SetBinContent(i + 1, j + 1, resolution);
      h2resolution->SetBinError(i + 1, j + 1, error);
      h2response->SetBinContent(i + 1, j + 1, mean);
      h2response->SetBinError(i + 1, j + 1, meanerror);
    }
  }
  // scale the response by the bin at 0,0
  float scale = h2response->GetBinContent(7, 7);
  h2response->Scale(1. / scale);
  out->Write();
  out->Close();
}
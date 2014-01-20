/*
 * analysis.cpp
 *
 *  Created on: 11 Sep 2013
 *      Author: Nicola Mori
 */

/* Analysis program for the example run. */

TString dictName("/path/to/libGGSReader.so");
TString fileName("RunOutput.root");
float cube1Size = 5.; // in cm
float primaryEnergy = 10.; // in GeV

void analysis() {

  // Load the library with dictionaries
  gROOT->ProcessLine(TString(".L ") + dictName);

  // Create the reader and open the data file
  GGSTRootReader *reader = new GGSTRootReader();
  reader->OpenFile(fileName);

  // Retrieve the hit reader
  GGSTHitsReader *hReader = reader->GetReader<GGSTHitsReader>();
  // Set which hit branches are to be read
  hReader->SetDetector("cube1Logic", true); // Read also particle hits
  hReader->SetDetector("cube2Logic");

  // Retrieve the MC truth reader
  GGSTMCTruthReader *mcReader = reader->GetReader<GGSTMCTruthReader>();

  // Prepare the histograms
  TH1I *eDepHisto = new TH1I("eDepHisto", "Total energy deposit (events inside acceptance)", 100, primaryEnergy * 0.8,
      primaryEnergy * 1.1); // Total energy release
  TGraph *impactGraph = new TGraph;
  impactGraph->SetTitle("Trajectory points at Z=0 (Red: hits on cube1)");
  TGraph *missGraph = new TGraph;
  missGraph->SetTitle("Trajectory points at Z=0 (Red: hits on cube1)");

  // Event loop
  for (int iEv = 0; iEv < reader->GetEntries(); iEv++) {
    reader->GetEntry(iEv);

    // Compute trajectory points at Z=0 by extrapolating the MC trajectory
    float *pos = mcReader->GetPrimaryParticle(0)->pos;
    float *mom = mcReader->GetPrimaryParticle(0)->mom;
    float deltaZ = pos[2]; // Upstream face is at Z=0
    float pointX = pos[0] + mom[0] / mom[2] * deltaZ;
    float pointY = pos[1] + mom[1] / mom[2] * deltaZ;
    // check if the point is on the upstream face of cube1
    if (fabs(pointX) < cube1Size / 2. && fabs(pointY) < cube1Size / 2.) {
      impactGraph->SetPoint(impactGraph->GetN(), pointX, pointY);
      // Compute total energy release for events inside acceptance (i.e. that hit the upstream face of the cube)
      double eDep = 0;
      if (hReader->GetNHits("cube1Logic") > 0)
        eDep += hReader->GetHit("cube1Logic", 0)->eDep; // Each detector has only one element, so only one hit at maximum
      if (hReader->GetNHits("cube2Logic") > 0)
        eDep += hReader->GetHit("cube2Logic", 0)->eDep;
      eDepHisto->Fill(eDep);
    }
    else {
      missGraph->SetPoint(missGraph->GetN(), pointX, pointY);
    }
  }

  new TCanvas;
  eDepHisto->Draw();

  new TCanvas;
  if (missGraph->GetN() > 0) {
    missGraph->SetMarkerStyle(7);
    missGraph->Draw("ap");
  }
  impactGraph->SetMarkerStyle(7);
  impactGraph->SetMarkerColor(kRed);
  if (missGraph->GetN() > 0)
    impactGraph->Draw("p");
  else
    impactGraph->Draw("ap");

}


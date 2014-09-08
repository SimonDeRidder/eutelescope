#include "TH1F.h"
#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"
#include "TProfile.h"
#include "TGraphPainter.h"
#include "TCanvas.h"

using namespace TMath;

const int startNumber = 46;
const int endNumber = 69;
const int dutSensorID = 20;
const Double_t pi = 3.14159265359;
const int anglebinSize = 20;
const double interval = 0.2;

void make_angle_sizeY()
{
	Double_t tempAngle;
	Int_t tempSize;
	UInt_t nEntries = 0;

	TCanvas canvas;
	TProfile outputHist("outputHist", "Measured cluster size in y vs fitted track incidence angle w.r.t. y", 							anglebinSize, -interval, interval);
	outputHist.GetXaxis()->SetTitle("Angle");
	outputHist.GetXaxis()->CenterTitle();
	outputHist.GetYaxis()->SetTitle("Average cluster size");
	outputHist.GetYaxis()->CenterTitle();
	gStyle->SetPadTickY(1);
	gStyle->SetOptStat(0);

	for (int i=0; i<=endNumber-startNumber; i++)
	{
		if(startNumber+i==55)// || startNumber+i==52)
		{//55 doesn't exist, 52 has huge errors
			continue;
		}
		TString inputFileName = "output/histograms/run0035";
		inputFileName.Append(TString::Itoa(startNumber+i,10));
		inputFileName.Append("-NTuple.root");
		TFile* inputFile = new TFile(inputFileName,"READ");
		TTree* inputTree;
		inputTree = (TTree) inputFile->Get("EUFit");
		inputTree->SetBranchAddress("dutTx",&tempAngle);
		inputTree->SetBranchAddress("dutClusterSizeY",&tempSize);
		nEntries = inputTree->GetEntries();
		if(nEntries==0)
		{
			continue;
		}
		for(int j=0; j<nEntries; j++)
		{
			inputTree->GetEntry(j);
			outputHist.Fill(((tempAngle/pi)*180),tempSize);
		}
	}

	outputHist.Draw();
	canvas.SaveAs("angle_sizeY.png");
	//inputFile->Close();
  return;
}

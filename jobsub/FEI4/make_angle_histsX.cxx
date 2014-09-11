#include "TH1F.h"
#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TGraphPainter.h"
#include "TCanvas.h"

using namespace TMath;
using namespace std;

const int startNumber = 46;
const int endNumber = 69;
const int dutSensorID = 20;
const Double_t pi = 3.14159265359;
const int nBins = 30;

void make_angle_histsX()
{
	TVectorD angles;
	Double_t tempAngle;
	UInt_t nEntries = 0;
	int index = 0;
	Double_t minAngle;
	Double_t maxAngle;
	Double_t angleSum;
	Double_t angleMean;
	for (int i=0; i<=endNumber-startNumber; i++)
	{
		if(startNumber+i==55)
		{//55 doesn't exist
			continue;
		}
		TString inputFileName = "output/histograms/run0035";
		TString runNr = TString::Itoa(startNumber+i,10);
		inputFileName.Append(runNr);
		inputFileName.Append("-NTuple.root");
		TFile* inputFile = new TFile(inputFileName,"READ");
		TTree* inputTree;
		inputTree = (TTree) inputFile->Get("EUFit");
		inputTree->SetBranchAddress("dutTx",&tempAngle);
		nEntries = inputTree->GetEntries();
		if(nEntries==0)
		{
			continue;
		}

		angles.ResizeTo(nEntries);
		inputTree->GetEntry(0);
		angles[0] = (tempAngle/pi)*180;
		minAngle = angles[0];
		maxAngle = angles[0];
		angleSum = 0.;
		angleMean = 0.;
		for(int j=0; j<nEntries; j++)
		{
			inputTree->GetEntry(j);
			angles[j] = (tempAngle/pi)*180;
			angleSum += (tempAngle/pi)*180;
			if(angles[j]<minAngle)
			{
				minAngle = angles[j];
			}
			if(angles[j]>maxAngle)
			{
				maxAngle = angles[j];
			}
		}
		angleMean = angleSum/nEntries;
		TCanvas canvas;
		TString histName = "Angular distribution in x for ";
		histName.Append(TString::Itoa(angleMean+0.5,10));
		histName.Append(" degrees");
		double dist  = (maxAngle - minAngle)/(nBins-2);
		TH1D* anglehist = new TH1D("ToT_Hist", histName, nBins, minAngle-dist, maxAngle+dist);
		anglehist->SetNameTitle("ToT_Hist", histName);
		for(int j=0; j<nEntries; j++)
		{
			anglehist->Fill(angles[j]);
		}
		anglehist->GetXaxis()->SetTitle("Angle");
		anglehist->GetXaxis()->CenterTitle();
		anglehist->Draw();

		TString outName = "angleDist/angleX_";
		outName.Append(TString::Itoa(angleMean+0.5,10));
		outName.Append(".png");
		canvas.SaveAs(outName);
		//inputFile->Close();
	}
  return;
}

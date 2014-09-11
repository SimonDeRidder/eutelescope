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

void make_angle_charge_hists()
{
	TVectorD charge;
	Double_t tempAngle;
	Double_t tempCharge;
	UInt_t nEntries = 0;
	int chargeMax = 0;
	Double_t angleSum;
	Double_t angleMean;
	int index = 0;
	for (int i=0; i<=endNumber-startNumber; i++)
	{
		if(startNumber+i==55)
		{//55 doesn't exist
			continue;
		}
		TString inputFileName = "output/histograms/run0035";
		inputFileName.Append(TString::Itoa(startNumber+i,10));
		inputFileName.Append("-NTuple.root");
		TFile* inputFile = new TFile(inputFileName,"READ");
		TTree* inputTree;
		inputTree = (TTree) inputFile->Get("EUFit");
		inputTree->SetBranchAddress("dutTx",&tempAngle);
		inputTree->SetBranchAddress("dutQ",&tempCharge);
		nEntries = inputTree->GetEntries();
		if(nEntries==0)
		{
			continue;
		}
		charge.ResizeTo(nEntries);
		chargeMax = 0;
		angleSum = 0;
		for(int j=0; j<nEntries; j++)
		{
			inputTree->GetEntry(j);
			charge[j] = tempCharge;
			angleSum += (tempAngle/pi)*180;
			if (tempCharge>chargeMax)
			{
				chargeMax = tempCharge;
			}
		}
		angleMean = angleSum/nEntries;
		TCanvas canvas;
		TString histName = "Measured total cluster charge disribution for incidence angle of ";
		stringstream ss;
		ss << angleMean;
		histName.Append(ss.str());
		histName.Append(" degrees, eta=");
		double tant = tan(pi*angleMean/360);
		if (tant<=0)
		{
			histName.Append("0");
		} else
		{
			stringstream s;
			ss << -log(tant);
			histName.Append(ss.str());
		}
		TH1D anglehist("ToT_Hist", histName, 100, 0, 100);
		for(int j=0; j<nEntries; j++)
		{
			anglehist.Fill(charge[j]);
		}
		anglehist.GetXaxis()->SetTitle("Total cluster charge");
		anglehist.GetXaxis()->CenterTitle();
		anglehist.Draw();

		TString outName = "ToT_";
		outName.Append(TString::Itoa(angleMean+0.5,10));
		outName.Append(".png");
		canvas.SaveAs(outName);
		//inputFile->Close();
	}
  return;
}

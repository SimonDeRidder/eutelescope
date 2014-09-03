#include "TH1F.h"
#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TGraphPainter.h"
#include "TCanvas.h"

using namespace TMath;

const int startNumber = 46;
const int endNumber = 69;
const int dutSensorID = 20;
const Double_t pi = 3.14159265359;

void make_angle_charge()
{
	TVectorD angle;
	TVectorD charge;
	Double_t tempAngle;
	Double_t tempCharge;
	UInt_t nEntries = 0;
	TVectorD angleList;
	TVectorD angleErrorList;
	TVectorD chargeList;
	TVectorD chargeErrorList;
	Double_t angleSum;
	Double_t chargeSum;
	Double_t angleResSum;
	Double_t chargeResSum;
	Double_t angleMean;
	Double_t chargeMean;
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
		angle.ResizeTo(nEntries);
		charge.ResizeTo(nEntries);
		for(int j=0; j<nEntries; j++)
		{
			inputTree->GetEntry(j);
			tempAngle = (tempAngle/pi)*180;
			angle[j] = tempAngle;
			charge[j] = tempCharge;
		}
		angleSum = 0;
		chargeSum = 0;
		angleResSum = 0;
		chargeResSum = 0;
		angleList.ResizeTo(index+1);
		angleErrorList.ResizeTo(index+1);
		chargeList.ResizeTo(index+1);
		chargeErrorList.ResizeTo(index+1);
		if(nEntries==1)
		{
			angleList[index] = angle[0];
			chargeList[index] = charge[0];
			angleErrorList[index] = 0;
			chargeErrorList[ndex] = charge[0];;
		} else
		{
			for(int j=0; j<nEntries; j++)
			{
				angleSum += angle[j];
				chargeSum += charge[j];
			}
			angleMean = angleSum/nEntries;
			chargeMean = chargeSum/nEntries;
			for(int j=0; j<nEntries; j++)
			{
				angleResSum += (angle[j]-angleMean)*(angle[j]-angleMean);
				chargeResSum += (charge[j]-chargeMean)*(charge[j]-chargeMean);
			}
			angleList[index] = angleMean;
			chargeList[index] = chargeMean;
			angleErrorList[index] = sqrt(angleResSum/(nEntries*(nEntries-1)));
			chargeErrorList[index] = sqrt(chargeResSum/(nEntries*(nEntries-1)));
		}
		index++;
	}
	TCanvas canvas;
	TGraphErrors outputGraphCharge(angleList,chargeList,angleErrorList,chargeErrorList);
	outputGraphCharge.SetTitle("measured cluster charge vs fitted track incidence angle");
	outputGraphCharge.Draw("AP");
	canvas.SaveAs("angle_charge.png");
	inputFile->Close();
	~angle;
	~charge;
	~angleList;
	~angleErrorList;
	~chargeList;
	~chargeErrorList;
	~inputFileName;
  return;
}

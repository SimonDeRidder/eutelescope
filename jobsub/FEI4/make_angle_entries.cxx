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

void make_angle_entries()
{
	TVectorD angle = TVectorD(0);
	Double_t tempAngle;
	UInt_t nEntries = 0;
	TVectorD angleList;
	TVectorD angleErrorList;
	TVectorD entryList;
	TVectorD entryErrorList;
	Double_t angleSum;
	Double_t angleResSum;
	Double_t angleMean;
	UInt_t off;
	int index = 0;
	for (int i=0; i<=endNumber-startNumber; i++)
	{
		if(startNumber+i==55 || startNumber+i==52)
		{//55 doesn't exist, 52 has huge errors
			off++;
			continue;
		}
		TString inputFileName = "output/histograms/run0035";
		inputFileName.Append(TString::Itoa(startNumber+i,10));
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
		angle.ResizeTo(nEntries);
		for(int j=0; j<nEntries; j++)
		{
			inputTree->GetEntry(j);
			tempAngle = (tempAngle/pi)*180;
			angle[j] = tempAngle;
		}
		angleSum = 0;
		angleResSum = 0;
		angleList.ResizeTo(index+1);
		entryList.ResizeTo(index+1);
		angleErrorList.ResizeTo(index+1);
		entryErrorList.ResizeTo(index+1);
		if(nEntries==1)
		{
			angleList[index] = angle[0];
			entryList[index] = 1;
			angleErrorList[index] = 0;
			entryErrorList[index] = 0;
		} else
		{
			for(int j=0; j<nEntries; j++)
			{
				angleSum += angle[j];
			}
			angleMean = angleSum/nEntries;
			for(int j=0; j<nEntries; j++)
			{
				angleResSum += (angle[j]-angleMean)*(angle[j]-angleMean);
			}
			angleList[index] = angleMean;
			entryList[index] = nEntries;
			angleErrorList[index] = sqrt(angleResSum/(nEntries*(nEntries-1)));
			entryErrorList[index] = 0;
		}
		index ++;
	}
	TCanvas canvas;
	TGraphErrors outputGraphEntries(angleList,entryList,angleErrorList,entryErrorList);
	outputGraphEntries.SetTitle("accepted number of tracks with DUT hits vs incidence angle");
	outputGraphEntries.GetHistogram()->GetXaxis()->SetTitle("Angle");
	outputGraphEntries.GetHistogram()->GetXaxis()->CenterTitle();
	outputGraphEntries.GetHistogram()->GetYaxis()->SetTitle("Number of accepted hits");
	outputGraphEntries.GetHistogram()->GetYaxis()->CenterTitle();
	canvas.SetLogy();
	outputGraphEntries.Draw("AP");
	canvas.SaveAs("angle_entries.png");
	inputFile->Close();
	~angle;
	~angleList;
	~angleErrorList;
	~entryList;
	~entryErrorList;
	~inputFileName;
  return;
}

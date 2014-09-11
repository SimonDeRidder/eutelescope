#include "TH1F.h"
#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TGraphPainter.h"
#include "TCanvas.h"

using namespace TMath;
using namespace std;

const int startNumber = 46;
const int endNumber = 69;
const int dutSensorID = 20;
const Double_t pi = 3.14159265359;
const Int_t numCols =80;


void make_angle_TOT_profile()
{
	double tempAngle;
	UInt_t nEntries = 0;
	double tempPoeX;
	double tempPosX;
	vector<int>* tempTot;
	vector<int>* tempPixX;

	int start;
	int step;

	int index;
	int pixColSum[numCols];
	int pixColCount[numCols];

	TCanvas canvas;
	TProfile2D outputGraph("TOT_Graph","measured TOT along hit for different fitted track incidence angles", numCols, 0, numCols, 96, 0, 96);
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
		inputTree = (TTree*) inputFile->Get("EUFit");
		inputTree->SetBranchAddress("dutTx",&tempAngle);
		inputTree->SetBranchAddress("dutPixTOT",&tempTot);
		inputTree->SetBranchAddress("dutPixX",&tempPixX);
		inputTree->SetBranchAddress("dutPOEX",&tempPoeX);
		inputTree->SetBranchAddress("dutXLocal",&tempPosX);
		nEntries = inputTree->GetEntries();
		if(nEntries==0)
		{
			cout << "WARNING: zero entry" << endl;
			continue;
		}
		if(nEntries==1)
		{
			cout << "WARNING: only 1 entry" << endl;
			continue;
		}
		for(int j=0; j<nEntries; j++)
		{
			inputTree->GetEntry(j);

			start = tempPixX->at(0);//should exist because of loop break when nEntries==0
			if (tempPoeX>=tempPosX)
			{
				for (int m=1; m<tempPixX->size(); m++)
				{
					if(tempPixX->at(m)>start)
					{
						start = tempPixX->at(m);
					}
				}
				step = -1;
			} else
			{
				for (int m=1; m<tempPixX->size(); m++)
				{
					if(tempPixX->at(m)<start)
					{
						start = tempPixX->at(m);
					}
				}
				step = 1;
			}
			for (int k=0; k<numCols; k++)
			{
				pixColSum[k] = 0;
				pixColCount[k] = 0;
			}
			for (int m=0; m<tempPixX->size(); m++)
			{
				index = (tempPixX->at(m)-start)/step;
				if (index<numCols && index>=0)
				{
					pixColSum[index] += tempTot->at(m);
					pixColCount[index]++;
				}
			}
			for (int k=0; k<numCols; k++)
			{
				if(pixColCount[k]!=0)
				{
					outputGraph.Fill(k,(tempAngle/pi)*180,(Double_t) pixColSum[k]);
				}
			}
		}
		inputFile->Close();
	}
	//canvas.DrawFrame(0.,0.,10.,96.);
	//canvas.SetTitle("measured TOT along hit for different fitted track incidence angles");

	outputGraph.GetXaxis()->SetTitle("N'th row along cluster");
	outputGraph.GetXaxis()->CenterTitle();
	outputGraph.GetYaxis()->SetTitle("Angle");
	outputGraph.GetYaxis()->CenterTitle();
	outputGraph.GetZaxis()->SetTitle("TOT");
	outputGraph.GetZaxis()->CenterTitle();
	gStyle->SetOptStat(0);
	outputGraph.Draw("LEGO2Z");
	canvas.SetTheta(45);
	canvas.SetPhi(-30);
	canvas.SaveAs("angle_tot_profile.png");
	canvas.SetTheta(89.999);
	canvas.SetPhi(0.00001);
	canvas.SaveAs("angle_tot_profile_above.png");
	delete inputFile;
	return;
}

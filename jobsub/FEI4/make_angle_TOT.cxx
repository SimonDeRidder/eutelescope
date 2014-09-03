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


void make_angle_TOT()
{
	//x-values for the pixel columns
	Double_t* xValues = new Double_t[numCols];
	for (int i=1; i<=numCols; i++)
	{
		xValues[i-1] = (Double_t) i;
	}

	double tempAngle;
	TVectorD angle;
	UInt_t nEntries = 0;
	Double_t angleSum;
	Double_t angleMean;
	UInt_t off = 0;
	Int_t listSize = 0;
	double tempPoeX;
	double tempPosX;
	vector<int>* tempTot;
	vector<int>* tempPixX;
	vector<Double_t> totSum(numCols);
	TVectorD angleList;
	TVectorD xList;
	TVectorD totList;
	for (int i=0; i<=endNumber-startNumber; i++)
	{
		if(startNumber+i==55)
		{//55 doesn't exist
			off++;
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
			off++;
			cout << "WARNING: zero entry" << endl;
			continue;
		}
		if(nEntries==1)
		{
			off++;
			cout << "WARNING: only 1 entry" << endl;
			continue;
		}
		angle.ResizeTo(nEntries);
		angleSum = 0;
		for (int k=0; k<numCols; k++)
		{
			totSum[k] = 0.;
		}
		for(int j=0; j<nEntries; j++)
		{
			inputTree->GetEntry(j);
			tempAngle = (tempAngle/pi)*180;
			angleSum += tempAngle;
			int start = tempPixX->at(0);//should exist because of loop break when nEntries==0
			int step;
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
			int index;
			int pixColSum[numCols];
			int pixColCount[numCols];
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
					totSum[k] += ((Double_t) pixColSum[k]);///pixColCount[k];
				}
			}
		}
		angleMean = angleSum/nEntries;
		angleList.ResizeTo(listSize+numCols);
		xList.ResizeTo(listSize+numCols);
		totList.ResizeTo(listSize+numCols);
		for(int k=0; k<numCols; k++)
		{
			angleList[listSize+k] = angleMean;
			xList[listSize+k] = xValues[k];
			totList[listSize+k] = totSum[k]/nEntries;//mean TOT per column
		}
		listSize+=numCols;
		inputFile->Close();
	}
	TCanvas canvas;
	canvas.DrawFrame(0.,0.,10.,96.);
	canvas.SetTitle("measured TOT along hit for different fitted track incidence angles");
	TGraph2D* graph2d = new TGraph2D("TOT_Graph","measured TOT along hit for different fitted track incidence angles",
								listSize, xList.GetMatrixArray(), angleList.GetMatrixArray(), totList.GetMatrixArray());
	graph2d->SetNpx(numCols);
	graph2d->SetNpy(96);
	graph2d->GetHistogram()->GetXaxis()->SetTitle("N'th row along cluster");
	graph2d->GetXaxis()->CenterTitle();
	graph2d->GetHistogram()->GetYaxis()->SetTitle("Angle");
	graph2d->GetYaxis()->CenterTitle();
	graph2d->GetHistogram()->GetZaxis()->SetTitle("TOT");
	graph2d->GetHistogram()->GetZaxis()->CenterTitle();
	graph2d->Draw("LEGO2Z");
	canvas.SetTheta(45);
	canvas.SetPhi(-30);
	canvas.SaveAs("angle_tot.png");
	canvas.SetTheta(89.999);
	canvas.SetPhi(0.00001);
	canvas.SaveAs("angle_tot_above.png");
	delete inputFile;
	return;
}

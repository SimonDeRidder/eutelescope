#include "TH1F.h"
#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraphPainter.h"
#include "TGraphErrors.h"
#include "TCanvas.h"

using namespace TMath;

const int startNumber = 46;
const int endNumber = 69;
const int dutSensorID = 20;
const Double_t pi = 3.14159265359;

void make_angle_res()
{
	TCanvas canvas;
	TVectorD angle;
	TVectorD res;
	Double_t tempAngle;
	Double_t tempRes;
	UInt_t nEntries = 0;
	TVectorD angleList;
	TVectorD angleErrorList;
	TVectorD resList;
	TVectorD resErrorList;
	Double_t angleSum;
	Double_t resSum;
	Double_t angleResSum;
	Double_t resResSum;
	Double_t angleMean;
	Double_t resMean;
	int index = 0;
	for (int i=0; i<=endNumber-startNumber; i++)
	{
		if(startNumber+i==55 || startNumber+i==52)
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
		TString resName = "resX_";
		resName.Append(TString::Itoa(dutSensorID,10));
		inputTree->SetBranchAddress(resName,&tempRes);
		nEntries = inputTree->GetEntries();
		if(nEntries==0)
		{
			off++;
			continue;
		}
		angle.ResizeTo(nEntries);
		res.ResizeTo(nEntries);
		for(int j=0; j<nEntries; j++)
		{
			inputTree->GetEntry(j);
			tempAngle = (tempAngle/pi)*180;
			angle[j] = tempAngle;
			res[j] = tempRes;
		}
		angleSum = 0;
		resSum = 0;
		angleResSum = 0;
		resResSum = 0;
		angleList.ResizeTo(index+1);
		angleErrorList.ResizeTo(index+1);
		resList.ResizeTo(index+1);
		resErrorList.ResizeTo(index+1);
		if(nEntries==1)
		{
			angleList[index] = angle[0];
			resList[index] = res[0];
			angleErrorList[index] = 0;
			resErrorList[index] = res[0];
		} else
		{
			for(int j=0; j<nEntries; j++)
			{
				angleSum += angle[j];
				resSum += res[j];
			}
			angleMean = angleSum/nEntries;
			resMean = resSum/nEntries;
			for(int j=0; j<nEntries; j++)
			{
				angleResSum += (angle[j]-angleMean)*(angle[j]-angleMean);
				resResSum += (res[j]-resMean)*(res[j]-resMean);
			}
			angleList[index] = angleMean;
			resList[index] = resMean;
			angleErrorList[index] = sqrt(angleResSum/(nEntries*(nEntries-1)));
			resErrorList[index] = sqrt(resResSum/(nEntries*(nEntries-1)));
		}
		index++;
	}
	TGraphErrors* outputGraphRes = new TGraphErrors(angleList,resList,angleErrorList,resErrorList);
	outputGraphRes->SetTitle("residual between fitted and measured hit vs fitted track incidence angle");
	outputGraphRes->GetHistogram()->GetXaxis()->SetTitle("Angle");
	outputGraphRes->GetHistogram()->GetXaxis()->CenterTitle();
	outputGraphRes->GetHistogram()->GetYaxis()->SetTitle("residuals");
	outputGraphRes->GetHistogram()->GetYaxis()->CenterTitle();
	outputGraphRes->Draw("AP");
	canvas.SaveAs("angle_res.png");
	inputFile->Close();
	~angle;
	~res;
	~angleList;
	~angleErrorList;
	~resList;
	~resErrorList;
	~inputFileName;
	~resName;
  return;
}

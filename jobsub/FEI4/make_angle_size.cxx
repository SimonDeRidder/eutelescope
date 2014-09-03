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

void make_angle_size()
{
	TVectorD angle;
	vector<int> size;
	Double_t tempAngle;
	Int_t tempSize;
	UInt_t nEntries = 0;
	TVectorD angleList;
	TVectorD angleErrorList;
	TVectorD sizeList;
	TVectorD sizeErrorList;
	Double_t angleSum;
	Double_t sizeSum;
	Double_t angleResSum;
	Double_t sizeResSum;
	Double_t angleMean;
	Double_t sizeMean;
	int index = 0;
	for (int i=0; i<=endNumber-startNumber; i++)
	{
		if(startNumber+i==55 || startNumber+i==52)
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
		inputTree->SetBranchAddress("dutClusterSizeX",&tempSize);
		nEntries = inputTree->GetEntries();
		if(nEntries==0)
		{
			continue;
		}
		angle.ResizeTo(nEntries);
		size.resize(nEntries);
		for(int j=0; j<nEntries; j++)
		{
			inputTree->GetEntry(j);
			tempAngle = (tempAngle/pi)*180;
			angle[j] = tempAngle;
			size[j] = tempSize;
		}
		angleSum = 0;
		sizeSum = 0;
		angleResSum = 0;
		sizeResSum = 0;
		angleList.ResizeTo(index+1);
		angleErrorList.ResizeTo(index+1);
		sizeList.ResizeTo(index+1);
		sizeErrorList.ResizeTo(index+1);
		if(nEntries==1)
		{
			angleList[index] = angle[0];
			sizeList[index] = size[0];
			angleErrorList[index] = 0;
			sizeErrorList[index] = size[0];;
		} else
		{
			for(int j=0; j<nEntries; j++)
			{
				angleSum += angle[j];
				sizeSum += size[j];
			}
			angleMean = angleSum/nEntries;
			sizeMean = sizeSum/nEntries;
			for(int j=0; j<nEntries; j++)
			{
				angleResSum += (angle[j]-angleMean)*(angle[j]-angleMean);
				sizeResSum += (size[j]-sizeMean)*(size[j]-sizeMean);
			}
			angleList[index] = angleMean;
			sizeList[index] = sizeMean;
			angleErrorList[index] = sqrt(angleResSum/(nEntries*(nEntries-1)));
			sizeErrorList[index] = sqrt(sizeResSum/(nEntries*(nEntries-1)));
		}
		index++;
	}
	//fitting function
	TF1* sizeFit = new TF1("sizeFit", "min([0] + abs([1]*tan((((x-[2])/180.)*pi))),80)", -2.,2.);
	sizeFit->SetParName(0, "Background");
	sizeFit->SetParName(1, "Factor");
	sizeFit->SetParName(2, "Phase");
	sizeFit->SetParameter(0, 0.997);
	sizeFit->SetParameter(1, 0.621);
	sizeFit->SetParameter(2, 0.1);
	sizeFit->SetLineWidth(1);
	sizeFit->SetLineColor(4);

	//create image
	TCanvas canvas;
	TGraphErrors outputGraphSize(angleList,sizeList,angleErrorList,sizeErrorList);
	outputGraphSize.Fit("sizeFit");
	TF1* fitResult = outputGraphSize.GetFunction("sizeFit");
	Double_t chi2 = fitResult->GetChisquare();
	Double_t ndf = fitResult->GetNDF();
	cout << "chi2/ndf= " << chi2 << "/" << ndf << " = " << chi2/ndf << endl;
	outputGraphSize.SetTitle("measured cluster size in x vs fitted track incidence angle w.r.t. x");
	outputGraphSize.GetHistogram()->GetXaxis()->SetTitle("Angle");
	outputGraphSize.GetHistogram()->GetXaxis()->CenterTitle();
	outputGraphSize.GetHistogram()->GetYaxis()->SetTitle("average cluster size");
	outputGraphSize.GetHistogram()->GetYaxis()->CenterTitle();
	outputGraphSize.Draw("AP");

	TText t;
	t.SetTextSize(0.03);
	t.SetTextAlign(22);
	Double_t xmin = 0;
	Double_t y = -4.1;
	t.DrawText(xmin, y, Form("%d", 0.));
	for (int i=1; i<4; i++)
	{
		t.DrawText(xmin+(90-((atan(exp(-i))*360)/pi)), y, Form("%d", i));
	}
	t.DrawText(xmin+90, y, "inf");
	t.DrawText(xmin+90+(atan(exp(-3))*360)/pi, y, Form("%d", 3));

	canvas.SaveAs("angle_size.png");
	//inputFile->Close();
  return;
}

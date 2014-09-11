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

void make_angle_ToTTotal()
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
	//fitting function
	TF1* chargeFit = new TF1("chargeFit", "[0] + abs([1]/cos(((x-[2])/180.)*pi))", -2.,2.);
	chargeFit->SetParName(0, "Background");
	chargeFit->SetParName(1, "Factor");
	chargeFit->SetParName(2, "Phase");
	chargeFit->SetParameter(0, 0.78);
	chargeFit->SetParameter(1, 7.4);
	chargeFit->SetParameter(2, -0.00001);
	chargeFit->SetLineWidth(1);
	chargeFit->SetLineColor(4);

	//create image
	TCanvas canvas;
	TGraphErrors outputGraphCharge(angleList,chargeList,angleErrorList,chargeErrorList);
	/*outputGraphCharge.Fit("chargeFit");
	TF1* fitResult = outputGraphCharge.GetFunction("chargeFit");
	Double_t chi2 = fitResult->GetChisquare();
	Double_t ndf = fitResult->GetNDF();
	cout << "chi2/ndf= " << chi2 << "/" << ndf << " = " << chi2/ndf << endl;*/
	outputGraphCharge.SetTitle("Measured total cluster charge vs fitted track incidence angle");
	outputGraphCharge.GetHistogram()->GetXaxis()->SetTitle("Angle");
	outputGraphCharge.GetHistogram()->GetXaxis()->CenterTitle();
	outputGraphCharge.GetHistogram()->GetYaxis()->SetTitle("Mean cluster charge (sum of ToT)");
	outputGraphCharge.GetHistogram()->GetYaxis()->CenterTitle();
	outputGraphCharge.Draw("AP");

	TText t;
	t.SetTextSize(0.03);
	t.SetTextAlign(22);
	Double_t xmin = 0;
	Double_t y = -40.;
	t.DrawText(xmin, y, Form("%d", 0.));
	for (int i=1; i<4; i++)
	{
		t.DrawText(xmin+(90-((atan(exp(-i))*360)/pi)), y, Form("%d", i));
	}
	t.DrawText(xmin+90, y, "inf");
	t.DrawText(xmin+90+(atan(exp(-3))*360)/pi, y, Form("%d", 3));

	canvas.SaveAs("angle_ToTTotal.png");
	inputFile->Close();
  return;
}

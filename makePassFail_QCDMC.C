#include <iostream>
#include <string>
#include "TCanvas.h"
#include "TH1.h"
#include "TFile.h"
#include <sstream>
#include "TStyle.h"
#include "assert.h"
#include "TF1.h"
#include "TMath.h"
#include "TLegend.h"
#include "TPaveStats.h"
#include "TPaveText.h"
#include <iomanip>
#include "IOColors.hh"
#include "TGraphErrors.h"
#include "TMatrixDSym.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"

using namespace std;

const static float fDATA_LUMI = 10000; //pb-1
const static float passFailRatio_fitrange_xmin = 50.0;
const static float passFailRatio_fitrange_xmax = 150.0;
std::vector<float> incMHTBins;
const static int nMHTbins = 5;
//const static int nMHTbins = 2;
//const static float arrHTbins[nHTbins] = {0,500,800,1000,1200,1400,8000};
const static float arrMHTbins[nMHTbins] = {200,300,450,600,8000};
//const static float arrMHTbins[nMHTbins] = {500,800};
//float CONST_C = 0.0217;
float CONST_C = 0.01;
const static bool bFLATSAMPLE = 0;
const static string sQCD_FLAT_FILE_NAME("QCD_Flat_8TeV.root");
const static float fQCD_FLAT_LUMI = 9998154/2.99913994E10; //pb-1
const Int_t nBins = 1; //mean + 6-systs
TFile *files[nBins];  

struct Predictions_t
{
		float incl_mean[nMHTbins];
		float incl_statErr[nMHTbins];
		float incl_fitErr[nMHTbins];
		float excl_mean[nMHTbins-1];
		float excl_statErr[nMHTbins-1];
		float excl_fitErr[nMHTbins-1];
		float incl_signal_mean[nMHTbins];
		float incl_signal_statErr[nMHTbins];
		float excl_signal_mean[nMHTbins-1];
		float excl_signal_statErr[nMHTbins-1];
};

void DumpHist(const TH1* hist)
{
	cout << ">>>> HIST INFO OF :" << hist->GetName() << endl;
	for (int bin = 1; bin < hist->GetNbinsX(); ++bin)
	{
		cout << "bin/val/err=" << bin << "\t" << hist->GetBinContent(bin) << "\t" << hist->GetBinError(bin) << endl;
	}
}


double GetFitFunctionError(TF1* f1, TFitResultPtr fitResPtr, const double val)
{
	TMatrixDSym covMat_ratio1pol2 = fitResPtr->GetCorrelationMatrix();
	double parErrp0 = fitResPtr->ParError(0);
	double parErrp1 = fitResPtr->ParError(1);

	const double err = sqrt( val*val*parErrp1*parErrp1*covMat_ratio1pol2(1, 1)*covMat_ratio1pol2(1, 1) + parErrp0*parErrp0*covMat_ratio1pol2(0, 0)*covMat_ratio1pol2(0, 0) + 2*val*parErrp0*parErrp1*covMat_ratio1pol2(0, 1) );
	double centralVal = f1->Eval(val);

	return err;
}
vector<float> 
GetVarBinVector (float down, float up5, float up10, float up20, float up50, 
					float width1, float width2, float width3, float width4)
{
	vector<float> result;
	float point = down;
	const unsigned nregion = 4;
	const float up [] = {up5, up10, up20, up50};
	const float step [] = {width1, width2, width3, width4};		//1j pet

	result.push_back (point);
	for (unsigned region = 0; region != nregion; ++ region)
	{
		while (point + step[region] < up[region] + 1)
		{
			point += step[region];
			result.push_back (point);
	 	};
	};
	
	if (point + 1 < up[nregion-1])
		result.push_back (up[nregion-1]);
	
	return result;
};

void NormalizeBinWidth(TH1* hist)
{
	assert (hist != NULL && "requirement failed"); //spec
	assert (hist->GetDimension() == 1 && "CommonTools::NormalizeBinWidth:: hist is not 1-D");

	for (unsigned bin = 1; bin <= unsigned (hist->GetNbinsX()); ++ bin)
	{
		const float width = hist->GetXaxis()->GetBinWidth (bin);
		hist->SetBinContent (bin, hist->GetBinContent (bin) / width);
		hist->SetBinError (bin, hist->GetBinError (bin) / width);
	};

}

auto_ptr<TH1> FillVarBinHist (const vector<float>& bins, TH1 *input)
{
	assert (input != NULL && "requirement failed"); //spec
	assert (input->GetDimension() == 1 && "requirement failed"); //spec
	const unsigned nbin = unsigned (input->GetNbinsX ());

	auto_ptr<TH1> result (new TH1F (input->GetName(), input->GetTitle(),
						 bins.size() - 1, &bins[0]));

	result->SetBinContent (0, input->GetBinContent (0));
	result->SetBinError (0, input->GetBinError (0));
	result->SetBinContent (bins.size(), input->GetBinContent (nbin));
	result->SetBinError (bins.size(), input->GetBinError (nbin));
	for (unsigned bin = 1; bin != nbin; ++ bin)
	{
		const float low = input->GetXaxis()->GetBinLowEdge (bin);
		const float high = input->GetXaxis()->GetBinUpEdge (bin);
		const unsigned bin1 = result->FindBin (0.99 * low + 0.01 * high);
		const unsigned bin2 = result->FindBin (0.01 * low + 0.99 * high);
		assert (bin1 == bin2 && "requirement failed: histogram bin edges don't match"); //spec

		const float va = input->GetBinContent (bin);
		const float ea = input->GetBinError (bin);
		const float vb = result->GetBinContent (bin2);
		const float eb = result->GetBinError (bin2);
		const float vc = va + vb;
		const float ec = sqrt (ea * ea + eb * eb);
		result->SetBinContent (bin2, vc);
		result->SetBinError (bin2, ec);
	};

	assert (result.get() != NULL && "postcondition failed"); //spec
	return result;
}

TH1* MakeVariableBinHist (TH1 *hist, const float xmin, const float xpoint1,
							const float xpoint2, const float xpoint3, const float xpoint4,
				 			const float width1, const float width2, const float width3, 
							const float width4)
{
	assert (hist != NULL && "CommonTools::MakeVariableBinHist:: hist is NULL!");
	assert (hist->GetDimension() == 1 && "CommonTools::MakeVariableBinHist:: hist is not 1-D");
	
  	auto_ptr<TH1> result = FillVarBinHist ( 
											GetVarBinVector(xmin, xpoint1, xpoint2, 
															xpoint3, xpoint4, width1, 
															width2, width3, width4)
											, hist);
  	NormalizeBinWidth(result.get());
  	return result.release ();
};

void PrintExclResults(const Predictions_t& gaus, const Predictions_t& exp,
							const Predictions_t& gaus_cplus, Predictions_t& exp_cplus)
{
	cout << setprecision(3) 
					<< setw(10) << " HT "
					<< setw(15) << "MHT"
					<< setw(15) << " mean " 
					<< setw(20) << "statFitErr"
					<< setw(30) << "Signal+/-stat"
					<< endl;
	for (int i=0; i< nMHTbins -1 ; ++i)
	{
		cout << setprecision(4) 
					/*//<< setw(5) << HTmin << "-" << HTmax << " & "
					<< setw(10) << incMHTBins.at(i) << "-" << incMHTBins.at(i+1) << " & "
					//<< setw(15) << gaus.excl_mean[i]   << "& $\\pm$ " << gaus.excl_statErr[i] << " & ("
					<< "\t" << gaus.excl_mean[i]   << "& $\\pm$ " << gaus.excl_statErr[i] << " & ("
					//<< setw(15) << gaus_cplus.excl_mean[i] << ") & "
					<< "\t" << gaus_cplus.excl_mean[i] << ") & "
					//<< setw(15) << exp.excl_mean[i]   << "& $\\pm$ " << exp.excl_statErr[i] << " & ("
					<< "\t" << exp.excl_mean[i]   << "& $\\pm$ " << exp.excl_statErr[i] << " & ("
					//<< setw(15) << exp_cplus.excl_mean[i] << ") & "
					<< "\t" << exp_cplus.excl_mean[i] << ") & "
					//<< setw(15) << gaus.excl_signal_mean[i] << "& $\\pm$ " << gaus.excl_signal_statErr[i]
					<< "\t" << gaus.excl_signal_mean[i] << "& $\\pm$ " << gaus.excl_signal_statErr[i]
					<< endl;
					*/

					<< setw(10) << incMHTBins.at(i) << "-" << incMHTBins.at(i+1) << " & "<< fixed 
					<< "\t" << gaus.excl_mean[i]   << "& $\\pm$ " << gaus.excl_statErr[i] << " & ("
					<< "\t" << gaus_cplus.excl_mean[i] << ") & "
					<< "\t" << exp.excl_mean[i]   << "& $\\pm$ " << exp.excl_statErr[i] << " & ("
					<< "\t" << exp_cplus.excl_mean[i] << ") & "
					<< "\t" << gaus.excl_signal_mean[i] << "& $\\pm$ " << gaus.excl_signal_statErr[i]
					<< endl;
	}
}
void PrintExclResults(const Predictions_t& res, const float HTmin=0, const float HTmax=0)
{
	cout << setprecision(3) 
					<< setw(10) << " HT "
					<< setw(15) << "MHT"
					<< setw(15) << " mean " 
					<< setw(20) << "statFitErr"
					<< setw(30) << "Signal+/-stat"
					<< endl;
	for (int i=0; i< nMHTbins -1 ; ++i)
	{
			const float mean          = res.excl_mean[i]; 
			const float statAndFitErr = sqrt(pow(res.excl_statErr[i],2)+pow(res.excl_fitErr[i],2));
			const float signal        = res.excl_signal_mean[i];
			const float signalErr     = res.excl_signal_statErr[i];

		cout << setprecision(4) 
					<< setw(5) << HTmin << "-" << HTmax << " & "
					<< setw(10) << incMHTBins.at(i) << "-" << incMHTBins.at(i+1) << " & "
					<< setw(15) << mean   << "& $\\pm$ " << statAndFitErr << " & "
					<< setw(15) << signal << "& $\\pm$ " << signalErr
					<< endl;
	}
}
void PrintInclResults(const Predictions_t& res, const float HTmin=0)
{
	cout << setprecision(3) 
					<< setw(10) << " HT "
					<< setw(15) << "MHT"
					<< setw(15) << " mean " 
					<< setw(20) << "statFitErr"
					<< setw(30) << "Signal+/-stat"
					<< endl;
	for (int i=0; i< nMHTbins ; ++i)
	{
			const float mean          = res.incl_mean[i]; 
			const float statAndFitErr = sqrt(pow(res.incl_statErr[i],2)+pow(res.incl_fitErr[i],2));
			const float signal        = res.incl_signal_mean[i];
			const float signalErr     = res.incl_signal_statErr[i];

		cout << setprecision(4) 
					<< setw(10) << HTmin << " & "
					<< setw(10) << incMHTBins.at(i) << " & "
					<< setw(15) << mean   << " & $\\pm$ " << statAndFitErr << " & "
					<< setw(15) << signal << " & $\\pm$ " << signalErr
					<< endl;
	}
	cout << "out of " << __FUNCTION__ << endl;
}


//cross section for each sample in pt order
const float xSec[] = {
	1759.549,  //pt 300-470
	113.8791,  //470-600
	26.9921,  //600-800
	3.550036, //8000-1000
	0.737844, //1000-1400
	0.03352235, //1400-1800
	0.001829005 //1800
};

const float nEvents[] = {
	5927300, // 300-470  #numbers from DBS, PREP page numbers are approximate
	3994848, // 470-600
	3992760, // 600-800
	3998563, //800-1000
	1964088, //1000-1400
	2000062, //1400-1800
	977586 //1800
};

void OpenFiles()
{
	//files[0] = new TFile ("Mean/qcd_allht.root");
	files[0] = new TFile ("qcd_all_HTranged.root");
	/*files[1] = new TFile ("Syst1/qcd_allht.root");
	files[2] = new TFile ("Syst2/qcd_allht.root");
	files[3] = new TFile ("Syst3/qcd_allht.root");
	files[4] = new TFile ("Syst4/qcd_allht.root");
	files[5] = new TFile ("Syst5/qcd_allht.root");
	files[6] = new TFile ("Syst6/qcd_allht.root");
*/
	for (int i=0; i<nBins; ++i)
	{
		if (files[i]->IsZombie())
		{
			cout << "QCD file # " << i << " not found!" <<  endl;
			assert (false);
		}
	}
}

TH1* GetHist(const std::string histname, const float scaleTo=1.0)
{
	TH1 *res_hist = 0;

	if (! bFLATSAMPLE)
	{
		TH1 *hists[nBins];

		for (int i=0; i<nBins; ++i)
		{
			hists[i] = 0; 
			if (files[i]->IsZombie())
			{
				cout << files[i]->GetName() << " not found!" <<  endl;
				assert (false);
			} else
			{
				hists[i] = dynamic_cast<TH1*> (files[i]->Get(histname.c_str()));
				if (hists[i] == 0 )
				{
					cout << "hist_pass " << histname << " not found in " << files[i]->GetName() << "!" << endl;
					assert (false);
				} else 
				{
					hists[i]->Sumw2();

					//temp: no scaling for lumi wgted samples
					//const float scale = scaleTo/ ( nEvents[i] / xSec[i] );
					//hists[i]->Scale(scale);

					if (i == 0) 
					{ 
						res_hist = dynamic_cast<TH1*> (hists[i]->Clone("histcopy")); 
						res_hist->SetDirectory(0);
					} else { res_hist->Add(hists[i]); }
				}
			}
		}

	} else
	{
		TFile f(sQCD_FLAT_FILE_NAME.c_str());
		if (f.IsZombie())
		{
			cout << __FUNCTION__ << ":" << __LINE__ << ":file " << sQCD_FLAT_FILE_NAME << " not found!"<< endl;
			assert(false);
		}
		TH1* hist = dynamic_cast<TH1*> (f.Get(histname.c_str()));
		if (hist == NULL)
		{
			cout << __FUNCTION__ << ":" << __LINE__ << ": " << histname << " not found!"<< endl;
			assert(false);
		}
		res_hist = dynamic_cast<TH1*> (hist->Clone());
		res_hist->SetDirectory(0);
	}

	return res_hist;
}

double expFitFunc_Cup(double *x, double *par)
{
	double fitval=0.0;
	//if(x[0]>0.0) fitval=par[0] * exp(par[1] * x[0]) + 0.03;
	//if(x[0]>0.0) fitval=par[0] * exp(par[1] * x[0]) + 0.0217;
	//if(x[0]>0.0) fitval=par[0] * exp(par[1] * x[0]) + (CONST_C * 2);  //100% error for C
	if(x[0]>0.0) fitval=par[0] * exp(par[1] * x[0]) + par[2];  //100% error for C
	else fitval=-1.0E6;
	return fitval;
}

double gausFitFunc_Cup(double *x, double *par)
{
	double arg = par[0] * exp(par[1] * x[0]);
	//double fitval = 1.0 / TMath::Erf (arg) - 1 + (CONST_C * 2); //100% ERROR FOR C
	double fitval = 1.0 / TMath::Erf (arg) - 1 + par[2]; //100% ERROR FOR C
	return fitval;
}


double expFitFunc(double *x, double *par)
{
	double fitval=0.0;
	//if(x[0]>0.0) fitval=par[0] * exp(par[1] * x[0]) + 0.03;
	//if(x[0]>0.0) fitval=par[0] * exp(par[1] * x[0]) + 0.0217;
	//if(x[0]>0.0) fitval=par[0] * exp(par[1] * x[0]) + CONST_C;
	if(x[0]>0.0) fitval=par[0] * exp(par[1] * x[0]) + par[2];
	else fitval=-1.0E6;
	return fitval;
}

double gausFitFunc(double *x, double *par)
{
	double arg = par[0] * exp(par[1] * x[0]);
	//double fitval = 1.0 / TMath::Erf (arg) - par[2];
	//double fitval = 1.0 / TMath::Erf (arg) - 1 + 0.03;
	//double fitval = 1.0 / TMath::Erf (arg) - 1 + 0.02;
	//double fitval = 1.0 / TMath::Erf (arg) - 1 + 0.0217;
	//double fitval = 1.0 / TMath::Erf (arg) - 1 + CONST_C;
	double fitval = 1.0 / TMath::Erf (arg) - 1 + par[2];
	return fitval;
}

Predictions_t GetPredictions(const TH1* hist, TF1* f1, const TH1* signalHist)
{
	assert(hist != NULL && "GetPredictions:: hist not found!");
	assert(f1 != NULL && "GetPredictions:: func not found!");
	assert(signalHist != NULL && "GetPredictions:: signalHist not found!");
	//new TCanvas(); gPad->SetLogy(); hist->SetStats(1); hist->DrawCopy(); //gPad->Print("control.eps");
	
	assert (hist->GetNbinsX() == signalHist->GetNbinsX() && "GetPredictions:: controlHist and signalHist have different binning!!!");

	//std::cout << "INITIAL INFO FOR HIST:"; hist->Print();
	int bin1 = 0;
	for (int bin =0; bin<hist->GetNbinsX(); ++bin) 
	{ 
		if (hist->GetBinLowEdge(bin)<incMHTBins.at(1)) continue;
		else { bin1 = bin; break; }
	}
	int bin2 = hist->GetNbinsX()+1;
	double err =0;
	double integral = hist->IntegralAndError(bin1, bin2, err);
//	std::cout << "Intergral for MHT> "<< incMHTBins.at(0) << ":" << integral << "+/-" << err << std::endl;
/*	double sig_sum = 0;
	double sig_err2 =0;
	for (int bin =0; bin<signalHist->GetNbinsX()+1; ++bin) 
	{
		cout << std::setw(10) << signalHist->GetBinLowEdge(bin) << ", "  << std::setw(10) 
							 << signalHist->GetXaxis()->GetBinUpEdge(bin)<< "]" 
							 << std::setw(10) << signalHist->GetBinContent(bin) << endl;

		if (signalHist->GetBinContent(bin)>0 && signalHist->GetBinCenter(bin)>= arrMHTbins[0] && signalHist->GetBinCenter(bin)< arrMHTbins[1])
		{
			sig_sum += signalHist->GetBinContent(bin);
			sig_err2 += pow(signalHist->GetBinError(bin),2);
		}
	}
	std::cout << red <<"Signal hist  [" << arrMHTbins[0] << "-" << arrMHTbins[1] << "]:= " << sig_sum << " +/- " << sqrt(sig_err2) << clearatt << std::endl;
*/
	double sumGaus[incMHTBins.size()];
	double Gaus_StatErr[incMHTBins.size()];
	double Gaus_FitErr[incMHTBins.size()];
	double sumGaus_excl[incMHTBins.size()-1];
	double Gaus_StatErr_excl[incMHTBins.size()-1];
	double Gaus_FitErr_excl[incMHTBins.size()-1];
	double signal_mean_excl[incMHTBins.size()-1];
	double signal_statErr_excl[incMHTBins.size()-1];
	double signal_mean_incl[incMHTBins.size()];
	double signal_statErr_incl[incMHTBins.size()];

	Predictions_t results;

	for (unsigned mhtBin = 0; mhtBin < incMHTBins.size(); ++mhtBin)
	{
		sumGaus[mhtBin]      = 0;
		Gaus_StatErr[mhtBin] = 0;
		Gaus_FitErr[mhtBin]  = 0;
		signal_mean_incl[mhtBin]  = 0;
		signal_statErr_incl[mhtBin] = 0;
		if (mhtBin+1<incMHTBins.size()) 
		{
			sumGaus_excl[mhtBin]      = 0;
			Gaus_StatErr_excl[mhtBin] = 0;
			Gaus_FitErr_excl[mhtBin]  = 0;
			signal_mean_excl[mhtBin]  = 0;
			signal_statErr_excl[mhtBin] = 0;
		}

		for (int bin = 0; bin <= hist->GetNbinsX()+1; ++bin)
		{
			if (hist->GetBinContent(bin)>0)
			{
							/* cout << std::setprecision(4) << std::setw(5) << bin << std::setw(3) << "[" 
							 << std::setw(10) << hist->GetBinLowEdge(bin) << ", "  << std::setw(10) 
							 << hist->GetXaxis()->GetBinUpEdge(bin)<< "]" 
							 << std::setw(10) << hist->GetBinContent(bin) 
							 << std::setw(10) << hist->GetBinContent(bin) * f1->Eval(hist->GetBinCenter(bin))
							 << endl;
							 */				//inclusive bin stuff

				const float binCenter = hist->GetBinCenter(bin);
				const float binVal    = hist->GetBinContent(bin);
				const float binErr    = hist->GetBinError(bin);
				const float funcVal   = f1->Eval(binCenter);
				const float res       = (binVal * funcVal);
				const float statErr2  = pow(binErr * funcVal, 2);
			 	//const float fitErr    = binVal * GetFitFunctionError(f1, binCenter);
			 	const float fitErr    = 99999999.99;
				const float binSig    = signalHist->GetBinContent(bin); 
				const float binSigStatErr2 = pow(signalHist->GetBinError(bin),2); 

			//cout << setprecision(3) << setw(15) << "MHT > " << incMHTBins.at(mhtBin) 
			//<< setw(20) << binSig  << "&$\\pm$" << binSigStatErr2 << std::endl;

				if (hist->GetBinCenter(bin) > incMHTBins.at(mhtBin))
				{
					sumGaus[mhtBin]      += res;
					Gaus_StatErr[mhtBin] += statErr2;
					Gaus_FitErr[mhtBin]  += fitErr;
					signal_mean_incl[mhtBin]  += binSig;
					signal_statErr_incl[mhtBin]  += binSigStatErr2;
				}

				//exclsuive bin stuff
				if (mhtBin+1<incMHTBins.size())
				{
					if (binCenter >= incMHTBins.at(mhtBin) 
							&& binCenter < incMHTBins.at(mhtBin+1))
					{
						//cout << green << "binCenter/MhtBin " << binCenter << "/" << incMHTBins.at(mhtBin) << clearatt<< endl;
						sumGaus_excl[mhtBin]      += res;
						Gaus_StatErr_excl[mhtBin] += statErr2;
						Gaus_FitErr_excl[mhtBin]  += fitErr;
						signal_mean_excl[mhtBin]  += binSig;
						signal_statErr_excl[mhtBin]  += binSigStatErr2;
						//cout << __LINE__ << ":" << mhtBin << " = " << binSig  << " [ " << signal_mean_excl[mhtBin]  << " ]" << endl; 
					}
				}
			}
		}

		Gaus_StatErr[mhtBin] = sqrt(Gaus_StatErr[mhtBin]);
		Gaus_StatErr_excl[mhtBin] = sqrt(Gaus_StatErr_excl[mhtBin]);
		signal_statErr_incl[mhtBin] = sqrt(signal_statErr_incl[mhtBin]);
		signal_statErr_excl[mhtBin] = sqrt(signal_statErr_excl[mhtBin]);

		results.incl_mean[mhtBin] = sumGaus[mhtBin];
		results.incl_statErr[mhtBin] = Gaus_StatErr[mhtBin];
		results.incl_fitErr[mhtBin] = Gaus_FitErr[mhtBin];
		results.incl_signal_mean[mhtBin] = signal_mean_incl[mhtBin];
		results.incl_signal_statErr[mhtBin] = signal_statErr_incl[mhtBin];

		
		if (mhtBin+1<incMHTBins.size())
		{
			std::stringstream excl_pred;
			excl_pred << setprecision(3) << setw(10) << incMHTBins.at(mhtBin) << "<MHT<" << incMHTBins.at(mhtBin+1)
				<< setw(20) << sumGaus_excl[mhtBin]  << "&$\\pm$" << Gaus_StatErr_excl[mhtBin] << " &$\\pm$ " << Gaus_FitErr_excl[mhtBin];
				//cout << ">>> " <<  excl_pred.str() << endl;
			results.excl_mean[mhtBin]    = sumGaus_excl[mhtBin];
			results.excl_statErr[mhtBin] = Gaus_StatErr_excl[mhtBin];
			results.excl_fitErr[mhtBin]  = Gaus_FitErr_excl[mhtBin];
			results.excl_signal_mean[mhtBin] = signal_mean_excl[mhtBin];
			results.excl_signal_statErr[mhtBin] = signal_statErr_excl[mhtBin];
		}
	}
	
	return results;
}


double GetAvgVal(TH1* hist, const float& xmin)
{
	assert(hist != NULL && "GetAvgVal:: hist is null!");
	assert(hist->GetDimension() == 1 && "GetAvgVal:: hist is not 1D!");
	
	double sum = 0, sumw = 0, N = 0;
	double wgt_sum = 0;
	
	for (int bin=1; bin<= hist->GetNbinsX(); ++bin)
	{
		if (hist->GetBinCenter(bin)< xmin) continue;
		
		const double w = hist->GetBinError(bin);
		const double v = hist->GetBinContent(bin);
		cout << "Adding bin " << bin << "[" << hist->GetBinLowEdge(bin) << "] = " << v << "+/-" << w <<endl; 
		sum     += v; 
		if (w>0) wgt_sum += (v / (w * w) ); 
		if (w>0) sumw    += (1.0/ (w * w) );
		++N;
	}

	const double avg = sum / N;
	const double wgtavg = sumw>0? wgt_sum/ sumw : 0.111111;

	cout << __FUNCTION__ << ": avg/wgtavg  = " << avg << " / " << wgtavg << endl;

	return wgtavg;

}
void makePassFail_QCDMC(const string numeHistName, const string denoHistName, 
			const string title, const string HTbinlabel, 
			const string signalHistName, const string controlHistName,
			const pair<unsigned, unsigned>& jetBin,
			const float fitrange_xmin = 50, const float fitrange_xmax = 150) 
{
	incMHTBins.clear();
	for (int i=0; i<nMHTbins; ++i) incMHTBins.push_back(arrMHTbins[i]);
	const bool logScale = true;
	const float scaleTo = fDATA_LUMI; // pb-1

	const int nSystematics = 6;
	vector<TH1*> vHist_ratio;
	//0= mean
	//>0 all systematic variations

/*	for (int i=0; i< nSystematics; ++i)
	{

		TH1* hist = GetRatioHist(i, numeHistName, denoHistName, scaleTo);

	}
*/



	TH1 *h_num = GetHist(numeHistName, scaleTo);
	TH1 *h_den = GetHist(denoHistName, scaleTo);

	//const float xmin = 50, xpoint1 = 160, xpoint2 = 300, xpoint3 = 400, xpoint4 = 1000, width1 = 10, width2 = 20, width3 = 50, width4 = 300;
	const float xmin = 50, xpoint1 = 150, xpoint2 = 200, xpoint3 = 240, xpoint4 = 640, width1 = 10, width2 = 20, width3 = 50, width4 = 400;
	TH1 *Hist_pass = MakeVariableBinHist (h_num, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4);
	TH1 *Hist_fail = MakeVariableBinHist (h_den, xmin, xpoint1, xpoint2, xpoint3, xpoint4, width1, width2, width3, width4);

	Hist_pass->SetMarkerColor(kBlack);
	Hist_pass->SetMarkerStyle(20);
	Hist_pass->SetLineColor(kBlack);
	/*new TCanvas();
	gPad->SetLogy();
	Hist_pass->SetLineColor(kRed);
	Hist_pass->SetMarkerColor(kRed);
	Hist_pass->DrawCopy();
	Hist_fail->DrawCopy("same");
	*/
	//gPad->Print("range.eps");
	cout << __LINE__ << ": Integrals= " << Hist_fail->Integral() << "/" << Hist_pass->Integral() << endl;
	if (Hist_fail->GetEntries()<1 || Hist_pass->GetEntries()<1) 
	{
		cout << __LINE__ << ": not enough intries to make the plots!!! " << Hist_fail->Integral() << "/" << Hist_pass->Integral() << endl;
		return;
	}

//	DumpHist(Hist_pass);
//	DumpHist(Hist_fail);

	Hist_pass->Divide(Hist_fail);

//	cout << ">>>>>>> AFTER DIVIDE <<<<< " << endl;
//	DumpHist(Hist_pass);
	const int maxbin = Hist_pass->GetMaximumBin();
	const double max = Hist_pass->GetBinContent(maxbin);
	//Hist_pass->GetYaxis()->SetRangeUser(-0.05, max+0.05);


	//fit range
	//const float fitrange_xmin = 50, fitrange_xmax = 120;
	stringstream newtitle;
	//newtitle << "Fit Range " << fitrange_xmin << "--" << fitrange_xmax << title;
	//newtitle << title << " (fit range = " <<  fitrange_xmin << "--" << fitrange_xmax << "), c = " << CONST_C << ";MHT [GeV];Ratio (r);";
	//newtitle << title << " (fit range = " <<  fitrange_xmin << "-" << fitrange_xmax << ");MHT [GeV];Ratio (r);";
	newtitle <<  title << " (from MG MC)" << ";MHT [GeV];Ratio (r);";


	gStyle->SetOptStat(0);	
	TCanvas *fitc = new TCanvas();
	//gPad->SetGridy();
	gPad->SetTickx();
	if (logScale) gPad->SetLogy();
	Hist_pass->SetLineColor(9);
	Hist_pass->SetTitle(newtitle.str().c_str());
	Hist_pass->SetLineWidth(2);
	Hist_pass->SetLabelSize(0.05,"XY");
	Hist_pass->SetTitleSize(0.05,"XY");

	gStyle->SetOptFit(1);
	//Hist_pass->SetStats(0);
	Hist_pass->Draw();
	//return;


	//use mean value of the last several bin as the 'C'
	const float mean_c_initial = GetAvgVal(Hist_pass, 300);

	//do fittings exp and gaus
	const float C_UPLIMIT = mean_c_initial+0.0001; //this only to get this values on the stat box
	const float C_LOLIMIT = mean_c_initial-0.0001;


	Hist_pass->SetMinimum(C_LOLIMIT/10);


	TF1 *expFit=new TF1("fit_1",expFitFunc,fitrange_xmin,fitrange_xmax,3);
	expFit->SetParameter(0,0.09);
	expFit->SetParameter(1,-0.0002);
	expFit->FixParameter(2,mean_c_initial);
	TFitResultPtr expFitResPtr = Hist_pass->Fit(expFit,"E0S","",fitrange_xmin, fitrange_xmax);
	gPad->Modified();
	gPad->Update();
	TPaveStats *e_stats = (TPaveStats*) Hist_pass->FindObject("stats");
	e_stats->SetTextColor(kRed);
	TPaveStats *exp_stats = (TPaveStats*) e_stats->Clone("exp_stats");

	TF1 *expFit2=new TF1("fit_2",expFitFunc,50,1000.0,3);
	expFit2->SetParameter(0,expFit->GetParameter(0));
	expFit2->SetParameter(1,expFit->GetParameter(1));
	expFit2->SetParameter(2,expFit->GetParameter(2));
	expFit2->SetLineColor(kRed+1);
	expFit2->SetLineWidth(2);

	TF1 *gausFit=new TF1("fit_3",gausFitFunc,fitrange_xmin, fitrange_xmax,3);
	gausFit->SetParameter(0,0.09); 
	gausFit->SetParameter(1,-0.0002);
	gausFit->FixParameter(2,mean_c_initial);
	gausFit->SetLineColor(kBlack);
	TFitResultPtr gausFitResPtr = Hist_pass->Fit(gausFit,"E0S","", fitrange_xmin, fitrange_xmax);
	gPad->Modified();
	gPad->Update();
	TPaveStats *gaus_stats = (TPaveStats*) Hist_pass->FindObject("stats");
	gaus_stats->SetTextColor(kBlack);
	//TPaveStats *gaus_stats = (TPaveStats*) g_stats->Clone("gaus_stats");

	//TF1 *gausFit2=new TF1("fit_4",gausFitFunc,50,1000.0,2);
	TF1 *gausFit2=new TF1("fit_4",gausFitFunc,50,1000.0,3);
	gausFit2->SetParameter(0,gausFit->GetParameter(0)); 
	gausFit2->SetParameter(1,gausFit->GetParameter(1));
	gausFit2->SetParameter(2,gausFit->GetParameter(2));
	gausFit2->SetLineColor(kBlack);
	gausFit2->SetLineWidth(2);

	gausFit2->Draw("same");
	//expFit2->Draw("same");
	//TLegend *leg  = new TLegend(0.7,0.8,0.9,0.9);
	TLegend *leg  = new TLegend(0.7,0.7,0.9,0.9);
	leg->AddEntry(gausFit2,"Gaussian");
	//leg->AddEntry(expFit2,"Exponential");
	leg->SetTextFont(42);
	leg->Draw();

	const float sb_xmin=0.2, sb_xmax=0.45, sb_ymin=0.7, sb_ymax=0.9;
/*	gaus_stats->SetX1NDC(sb_xmin);
	gaus_stats->SetX2NDC(sb_xmax);
	gaus_stats->SetY1NDC(sb_ymin);
	gaus_stats->SetY2NDC(sb_ymax);
*/	gaus_stats->SetX1NDC(sb_xmax);
	gaus_stats->SetX2NDC(sb_xmax+0.25);
	gaus_stats->SetY1NDC(sb_ymin);
	gaus_stats->SetY2NDC(sb_ymax);
	gaus_stats->Draw("same");
	/*exp_stats->SetX1NDC(sb_xmax);
	exp_stats->SetX2NDC(sb_xmax+0.25);
	exp_stats->SetY1NDC(sb_ymin);
	exp_stats->SetY2NDC(sb_ymax);
	exp_stats->Draw("same");
*/
	//fit quality
	stringstream gaus_fit_res, exp_fit_res;
	const float gausFit_goodness =  gausFit->GetChisquare()/gausFit->GetNDF();
	gaus_fit_res << "#chi^{2}/ndof = " << gausFit->GetChisquare()
			 << "/"<< gausFit->GetNDF() << " = " << gausFit_goodness;
	const float expFit_goodness =  expFit->GetChisquare()/expFit->GetNDF();
	exp_fit_res << "#chi^{2}/ndof = " << expFit->GetChisquare()
			 << "/"<< expFit->GetNDF() << " = " << expFit_goodness;

	cout << gaus_fit_res.str() << endl;
	cout << exp_fit_res.str() << endl;

	TPaveText *pt1 = new TPaveText(0.5,0.8,0.7,0.9);
	pt1->AddText(gaus_fit_res.str().c_str());
	pt1->AddText(exp_fit_res.str().c_str());
//	pt1->Draw("same");

	//stringstream epsname;
	//epsname << "factnomht/HTbin" << HTbin << "_fitrange_" << fitrange_xmin << "to" << fitrange_xmax << ".eps";
	//epsname << "factnomht_" << HTbinlabel << ".eps";
	//gPad->Print(epsname.str().c_str());
	gPad->Print("ratios.eps");
	delete fitc;
	//return;

	 //make a predicion
/*	TH1 *sigHist      = GetHist(signalHistName, scaleTo);
	TH1 *controlHist  = GetHist(controlHistName, scaleTo);

	
	if (sigHist->GetNbinsX() != controlHist->GetNbinsX())
	{
		cout << red << __LINE__ << ": Mismatch in number of sigHist/controlHist bins: " << sigHist->GetNbinsX() << "/" << controlHist->GetNbinsX() << clearatt << endl; 
		assert(false);
	}
*/
//	new TCanvas();
//	sigHist->SetLineColor(kRed);
//	sigHist->Rebin(50);
//	controlHist->Rebin(50);
//	sigHist->Draw();
//	controlHist->Draw("same");


	//collect results
/*	Predictions_t pred_gaus       = GetPredictions(controlHist, gausFit2, sigHist); 
	Predictions_t pred_gaus_sigma = GetPredictions(controlHist, gausFit_sigma, sigHist); 
	Predictions_t pred_exp        = GetPredictions(controlHist, expFit2, sigHist); 
	Predictions_t pred_exp_sigma  = GetPredictions(controlHist, expFit_sigma, sigHist); 

//	PrintExclResults(pred_gaus);
//	PrintExclResults(pred_exp);
	PrintExclResults(pred_gaus, pred_exp, pred_gaus_sigma, pred_exp_sigma);
*/
	/*****************************************/
	// TEMP HACK TO GET RESULTS for all HT bins 
	// using these inclusive fits
	/*****************************************/
/*	vector<pair<float, float> > htBins_temp;
	pair<float, float> htbin1(500,800);	
	pair<float, float> htbin2(800,1000);	
	pair<float, float> htbin3(1000,1250);	
	pair<float, float> htbin4(1250,1500);	
	pair<float, float> htbin5(1500,8000);	

	htBins_temp.push_back(htbin1);
	htBins_temp.push_back(htbin2);
	htBins_temp.push_back(htbin3);
	htBins_temp.push_back(htbin4);
	htBins_temp.push_back(htbin5);

	TFile file("qcd_all.root");
	if (file.IsZombie()) 
	{ 
		cout << __FUNCTION__ << ":" << __LINE__ 
				<< "File to get exlcusive HT prediction is not found!" << endl;assert (false);
	}

	for (unsigned htbin=0; htbin<htBins_temp.size(); ++htbin)
	{
		stringstream folder, control_hist_name, signal_hist_name;
		folder << "Hist/Njet" << jetBin.first << "to" << jetBin.second 
			<< "HT"   << htBins_temp.at(htbin).first << "to" << htBins_temp.at(htbin).second 
			<< "MHT0to8000";
		control_hist_name << folder.str() << "/smeared_failFineBin1"; 
		signal_hist_name << folder.str() << "/smeared_signalFineBin"; 
		TH1* control_hist = (TH1*) (file.Get(control_hist_name.str().c_str())); 
		if (control_hist == NULL) { 
			cout << __LINE__ << ": control hist " << control_hist_name.str() << " not found for htbin " << htBins_temp.at(htbin).first 
				<< "-" << htBins_temp.at(htbin).second << endl; 
			assert(false);		
		} 
		TH1* signal_hist = (TH1*) (file.Get(signal_hist_name.str().c_str())); 
		if (signal_hist == NULL) { 
			cout << __LINE__ << ": signal hist not found for htbin " << htBins_temp.at(htbin).first << "-" << htBins_temp.at(htbin).second << endl; 
			assert(false);
		} 

		control_hist->Print();
		signal_hist->Print();

		//temp hack to get # of bins the same for easy debugging
		//signal_hist->Rebin(2);

		Predictions_t pred_gaus       = GetPredictions(control_hist, gausFit2, signal_hist); 
//		Predictions_t pred_gaus_sigma = GetPredictions(control_hist, gausFit_sigma, signal_hist); 
		Predictions_t pred_exp        = GetPredictions(control_hist, expFit2, signal_hist); 
//		Predictions_t pred_exp_sigma  = GetPredictions(control_hist, expFit_sigma, signal_hist); 

		cout << ">>>>>>>>>>>>>>> PREDICTIONS FOR " << jetBin.first << "-" << jetBin.second 
					<< ", HT=" << htBins_temp.at(htbin).first << "-" << htBins_temp.at(htbin).second << endl;
//		PrintExclResults(pred_gaus, pred_exp, pred_gaus_sigma, pred_exp_sigma);
	}
	file.Close();
*/


	//draw systematic band around the fits using TGraphErrors
	//These total syst error values comes from final predictions from : abcdMethod_MC
	
	//now get the corresponding y errors for each jet bins
	vector<pair<float,float> > totsyst;
	if (jetBin.first==3 && jetBin.second == 5)
	{
		//gaus /exp mean value
		totsyst.push_back(make_pair(50, 0.066229));
		totsyst.push_back(make_pair(70, 0.0484432));
		totsyst.push_back(make_pair(90, 0.0467324));
		totsyst.push_back(make_pair(110, 0.0397972));
		totsyst.push_back(make_pair(130, 0.021121));
		totsyst.push_back(make_pair(150, 0.00611814));
		totsyst.push_back(make_pair(190, 0.00241428));
		totsyst.push_back(make_pair(230, 0.000976787));
		totsyst.push_back(make_pair(270, 0.00248154));
		totsyst.push_back(make_pair(370, 0.00273476));
		totsyst.push_back(make_pair(520, 0.00519916));
		totsyst.push_back(make_pair(670, 0.0186649));
		totsyst.push_back(make_pair(820, 0.0221335));
		totsyst.push_back(make_pair(820, 0.0221335));
		totsyst.push_back(make_pair(820, 0.0221335));
		totsyst.push_back(make_pair(820, 0.0221335));
		totsyst.push_back(make_pair(820, 0.0221335));
		totsyst.push_back(make_pair(820, 0.0221335));
		totsyst.push_back(make_pair(820, 0.0221335));
		totsyst.push_back(make_pair(820, 0.0221335));
	} else if (jetBin.first==6 && jetBin.second == 7)
	{
		totsyst.push_back(make_pair(50, 0.0328614));
		totsyst.push_back(make_pair(70, 0.0211245));
		totsyst.push_back(make_pair(90, 0.0560358));
		totsyst.push_back(make_pair(110, 0.0740068));
		totsyst.push_back(make_pair(130, 0.0569941));
		totsyst.push_back(make_pair(150, 0.0242288));
		totsyst.push_back(make_pair(190, 0.0115296));
		totsyst.push_back(make_pair(230, 0.00623851));
		totsyst.push_back(make_pair(270, 0.0367438));
		totsyst.push_back(make_pair(370, 0.025785));
		totsyst.push_back(make_pair(520, 0.044696));
		totsyst.push_back(make_pair(670, 0.118435));
		totsyst.push_back(make_pair(820, 0.124019));
		totsyst.push_back(make_pair(820, 0.124019));
		totsyst.push_back(make_pair(820, 0.124019));
		totsyst.push_back(make_pair(820, 0.124019));
		totsyst.push_back(make_pair(820, 0.124019));
		totsyst.push_back(make_pair(820, 0.124019));
		totsyst.push_back(make_pair(820, 0.124019));
		totsyst.push_back(make_pair(820, 0.124019));

	} else if (jetBin.first>=8)
	{
		totsyst.push_back(make_pair(50, 0.0683995));
		totsyst.push_back(make_pair(70, 0.0748122));
		totsyst.push_back(make_pair(90, 0.0754261));
		totsyst.push_back(make_pair(110, 0.0852469));
		totsyst.push_back(make_pair(130, 0.0852212));
		totsyst.push_back(make_pair(150, 0.0415323));
		totsyst.push_back(make_pair(190, 0.0265027));
		totsyst.push_back(make_pair(230, 0.0681772));
		totsyst.push_back(make_pair(270, 0.047416));
		totsyst.push_back(make_pair(370, 0.0882638));
		totsyst.push_back(make_pair(520, 0.152387));
		totsyst.push_back(make_pair(670, 0.283193));
		totsyst.push_back(make_pair(820, 0.493378));
		totsyst.push_back(make_pair(820, 0.493378));
		totsyst.push_back(make_pair(820, 0.493378));
		totsyst.push_back(make_pair(820, 0.493378));
		totsyst.push_back(make_pair(820, 0.493378));
		totsyst.push_back(make_pair(820, 0.493378));
		totsyst.push_back(make_pair(820, 0.493378));
		totsyst.push_back(make_pair(820, 0.493378));
	}



const float xbin_edges[] = {50,60,70,80,90,100,110,120,130,140,150,160,180,200,220,240,260,280,300,400,1000};
const Int_t nbin_edges = 21; 

	if ( totsyst.size() != (nbin_edges-1) )
	{
		cout << __FUNCTION__ << ":" << __LINE__ << ": totsyst vector size (" << totsyst.size() 
				<< ") does not match with number of bins of (" << nbin_edges -1 << ")" << endl;
				assert(false);
	}

Float_t xvals[nbin_edges -1 ]; //common for both fits
Float_t xerr[nbin_edges-1]; //common for both fits
Float_t yvals_gaus[nbin_edges - 1];
Float_t yvals_exp[nbin_edges - 1];
Float_t yerr_gaus[nbin_edges - 1];
Float_t yerr_exp[nbin_edges - 1];
	for (int i = 1; i < nbin_edges; ++i)
	{
		const double mid_x = (xbin_edges[i] + xbin_edges[i-1])/2.0;
		xerr[i-1]        = (xbin_edges[i] - xbin_edges[i-1]) /2.0;
		xvals[i-1]       = mid_x;
		yvals_gaus[i-1]  = gausFit2->Eval(mid_x);
		yvals_exp[i-1]   = expFit2->Eval(mid_x);
		yerr_exp[i - 1]  = totsyst.at(i-1).second;
		const double gaus_fiterr    = GetFitFunctionError(gausFit2, gausFitResPtr, mid_x);
		yerr_gaus[i - 1] =  sqrt( pow(totsyst.at(i-1).second,2) + pow(gaus_fiterr, 2));

	}
	
	
	TCanvas *tgc = new TCanvas();
	gPad->SetLogy();
	tgc->cd();
//	tgc->SetGrid();
	gPad->SetName("TG");
	//TGraphErrors *tg = new TGraphErrors( nbin_edges-1, xvals, yvals_gaus);
	TGraphErrors *tg1 = new TGraphErrors( nbin_edges-1, xvals, yvals_gaus, xerr, yerr_gaus);
	//TGraphErrors *tg_exp = new TGraphErrors( nbin_edges-1, xvals, yvals_exp);
	TGraphErrors *tg2_exp = new TGraphErrors( nbin_edges-1, xvals, yvals_exp, xerr, yerr_gaus);

	stringstream tg_title;
	if (jetBin.second>100) tg_title << "#geq " << jetBin.first << " Jets;MHT [GeV];Ratio (r);";
	else tg_title << jetBin.first << "-" << jetBin.second << " Jets;MHT [GeV];Ratio (r);";
	tg1->SetTitle(tg_title.str().c_str());
	tg1->GetXaxis()->SetLabelSize(0.05);
	tg1->GetYaxis()->SetLabelSize(0.05);
	tg1->GetXaxis()->SetTitleSize(0.05);
	tg1->GetYaxis()->SetTitleSize(0.05);
	cout << "y offset = " << tg1->GetYaxis()->GetTitleOffset() << endl;
	cout << "x offset = " << tg1->GetXaxis()->GetTitleOffset() << endl;
	tg1->GetYaxis()->SetTitleOffset(0.8);
	tg1->GetXaxis()->SetTitleOffset(0.9);

	//tg->GetXaxis()->SetRangeUser(50,700); 
	tg1->GetXaxis()->SetRangeUser(50,700); 
	//tg_exp->GetXaxis()->SetRangeUser(50,700); 
	tg2_exp->GetXaxis()->SetRangeUser(50,700); 
	tg2_exp->GetYaxis()->SetRangeUser(0.01,5); 
	tg1->GetYaxis()->SetRangeUser(0.01,5); 
	tg1->SetLineColor(kBlue);
	//tg->SetMarkerColor(kBlack);
	//tg->SetMarkerStyle(kCircle);
	tg1->SetFillColor(kGreen);
	tg1->SetFillStyle(3001);

	tg2_exp->SetFillStyle(3013);
	tg2_exp->SetFillColor(kRed);

	tg1->DrawClone("E3AL");
//	tg2_exp->DrawClone("E3Lsame");

	//tg->DrawClone("PEsame");
	Hist_pass->DrawClone("same");
	gausFit2->Draw("same");
//	expFit2->Draw("same");
	leg->Draw();


	stringstream tgeps;
	tgeps << "tg_gaus_" << jetBin.first << "to" << jetBin.second << ".eps";
	gPad->Print(tgeps.str().c_str());
//	tg1->Print();

	delete tgc;
}

void makePassFail_QCDMC()
{

	//jetbins
	//
	vector<pair<unsigned, unsigned> > jetBins;
	vector<pair<float, float> > htBins, mhtBins;
	pair<float, float> jetbin1(2,2);	
	pair<float, float> jetbin2(3,5);	
	pair<float, float> jetbin3(6,7);	
	pair<float, float> jetbin4(8,1000);	

	pair<float, float> htbin1(500,800);	
	pair<float, float> htbin2(800,1000);	
	pair<float, float> htbin3(1000,1250);	
	pair<float, float> htbin4(1250,1500);	
	pair<float, float> htbin5(1500,8000);	


	pair<float, float> mhtbin1(0,8000);	
	
	//jetBins.push_back(jetbin1);
	jetBins.push_back(jetbin2);
	jetBins.push_back(jetbin3);
	jetBins.push_back(jetbin4);
	
	htBins.push_back(htbin1);
	htBins.push_back(htbin2);
	htBins.push_back(htbin3);
	htBins.push_back(htbin4);

	mhtBins.push_back(mhtbin1);

	vector<string> dphibins;
//	dphibins.push_back("0.15");
	dphibins.push_back("0.20");
//	dphibins.push_back("0.25");
//	dphibins.push_back("0.30");
//	dphibins.push_back("0.35");
//	dphibins.push_back("0.40");

	OpenFiles();
	TCanvas *c = new TCanvas("print");	
	c->Draw();
	c->Print("ratios.eps[");
	
/*	for (unsigned d = 0; d < dphibins.size(); ++d)
	{	
		for (unsigned jetbin = 0; jetbin < jetBins.size(); ++jetbin)
		{	
			for (unsigned htbin = 0; htbin < htBins.size(); ++htbin)
			{
				for (unsigned mhtbin = 0; mhtbin < mhtBins.size(); ++mhtbin)
				{
					stringstream htrange;
					htrange << htBins.at(htbin).first << "<HT<" << htBins.at(htbin).second << " GeV, "
						<< mhtBins.at(mhtbin).first << "<MHT<" << mhtBins.at(mhtbin).second << " GeV";

					stringstream folder;
					folder << "Hist/Njet" << jetBins.at(jetbin).first << "to" << jetBins.at(jetbin).second 
						<< "HT"   << htBins.at(htbin).first << "to" << htBins.at(htbin).second 
						<< "MHT"  << mhtBins.at(mhtbin).first << "to" << mhtBins.at(mhtbin).second;

					stringstream njetlabel, title, numeHistName, denoHistName, signalHistName, controlHistName;
					njetlabel << jetBins.at(jetbin).first << "-" << jetBins.at(jetbin).second << " Jets";
					title << "Njet " << njetlabel.str() << ", " << htrange.str() << ", #Delta #Phi _{min}<" << dphibins.at(d);
					numeHistName << folder.str() << "/smear_signal";
					//denoHistName << folder.str() << "/smeared_fail1";
					denoHistName << folder.str() << "/smeared_fail" << d;
					signalHistName << folder.str() << "/smear_signalFineBin";
					controlHistName << folder.str() << "/smeared_failFineBin" << d;
					//makePassFail_QCDMC(numeHistName.str(), denoHistName.str(), 
					//		title.str(), htrange.str(), signalHistName.str(), controlHistName.str()); 
					makePassFail_QCDMC(numeHistName.str(), denoHistName.str(), 
							title.str(), htrange, signalHistName.str(), controlHistName.str(), jetBins.at(jetbin), 50, 150); 
				}

			}
		}
	}
*/	


	for (unsigned d = 0; d < dphibins.size(); ++d)
	{	
		for (unsigned jetbin = 0; jetbin < jetBins.size(); ++jetbin)
		{	
			stringstream folder;
			folder << "Hist/Njet" << jetBins.at(jetbin).first << "to" << jetBins.at(jetbin).second; 
			string htrange("HT>500 GeV");
			cout << folder.str() << endl;

			stringstream njetlabel, title, numeHistName, denoHistName, signalHistName, controlHistName;
			njetlabel << jetBins.at(jetbin).first << "-" << jetBins.at(jetbin).second << " Jets";
			//title << "Jets " << njetlabel.str() << ", " << htrange << ", #Delta #Phi _{min}<" << dphibins.at(d);
			if (jetBins.at(jetbin).second>100) title << "#geq" << jetBins.at(jetbin).first << " Jets";
			else title << njetlabel.str();
			numeHistName << folder.str() << "/smeared_signal";
			denoHistName << folder.str() << "/smeared_fail" << d;
			signalHistName << folder.str() << "/smeared_signalFineBin";
			controlHistName << folder.str() << "/smeared_failFineBin" << d;
			makePassFail_QCDMC(numeHistName.str(), denoHistName.str(), 
					title.str(), htrange, signalHistName.str(), controlHistName.str(), jetBins.at(jetbin), 50, 150); 
		}
	}


	c->Print("ratios.eps]");
}

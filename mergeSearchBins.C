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

using namespace std;

void mergeSearchBins()
{

	vector<string> histnames;
	histnames.push_back("smeared_signal");
	histnames.push_back("smeared_signalFineBin");
	histnames.push_back("smeared_fail0");
	histnames.push_back("smeared_fail1");
//	histnames.push_back("smeared_fail2");
	histnames.push_back("smeared_failFineBin0");
	histnames.push_back("smeared_failFineBin1");
//	histnames.push_back("smeared_failFineBin2");
	histnames.push_back("smeared_ht");
	histnames.push_back("reco_mht");
	histnames.push_back("reco_ht");
	histnames.push_back("recojet1_pt");
	histnames.push_back("recojet2_pt");

	histnames.push_back("reco_signal");
	histnames.push_back("reco_signalFineBin");
	histnames.push_back("reco_fail0");
	histnames.push_back("reco_fail1");
//	histnames.push_back("reco_fail2");
	histnames.push_back("reco_failFineBin0");
	histnames.push_back("reco_failFineBin1");
//	histnames.push_back("reco_failFineBin2");

	//for debugging
	histnames.push_back("gen_ht");


	//jetbins
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
	
//	jetBins.push_back(jetbin1);
	jetBins.push_back(jetbin2);
	jetBins.push_back(jetbin3);
	jetBins.push_back(jetbin4);
	
	htBins.push_back(htbin1);
	htBins.push_back(htbin2);
	htBins.push_back(htbin3);
	htBins.push_back(htbin4);
	htBins.push_back(htbin5);

	mhtBins.push_back(mhtbin1);

	vector<string> dphibins;
	dphibins.push_back("0.15");
	dphibins.push_back("0.20");
//	dphibins.push_back("0.25");
//	dphibins.push_back("0.30");
//	dphibins.push_back("0.35");
//	dphibins.push_back("0.40");

	//TFile fold("qcd_all.root");
//	TFile fold("TotBG_MCpuWgted_Lumi_12187pb.root");
	TFile fold("MG_QCD_1945pb.root");
	if (fold.IsZombie())
	{
		cout << "Input root file not found!" << endl;
		assert(false);
	}

	TFile fnew("MG_QCD_1945pb_HTranged.root","RECREATE");
	//TFile fnew("seed_ht_19p4fb.root","RECREATE");
	//TFile fnew("TotBG_MCpuWgted_Lumi_12187pb_HTranged.root","RECREATE");
	fnew.mkdir("Hist");
	
	for (unsigned jetbin = 0; jetbin < jetBins.size(); ++jetbin)
	{	
		stringstream jetbinfolder;
		jetbinfolder << "Njet" << jetBins.at(jetbin).first << "to" << jetBins.at(jetbin).second;

		fnew.cd("Hist");
		gDirectory->mkdir(jetbinfolder.str().c_str());
		gDirectory->cd(jetbinfolder.str().c_str());
		gDirectory->pwd();
		gDirectory->ls();

		for (unsigned h = 0; h < histnames.size(); ++h)
		{	
			TH1* hist_temp = 0;

			for (unsigned htbin = 0; htbin < htBins.size(); ++htbin)
			{
				for (unsigned mhtbin = 0; mhtbin < mhtBins.size(); ++mhtbin)
				{
					stringstream folder;
					folder << "Hist/" << jetbinfolder.str() 
						<< "HT"   << htBins.at(htbin).first << "to" << htBins.at(htbin).second 
						<< "MHT"  << mhtBins.at(mhtbin).first << "to" << mhtBins.at(mhtbin).second;

					stringstream histname;
					histname << folder.str() << "/" << histnames.at(h);

					TH1* temp = dynamic_cast<TH1*> (fold.Get(histname.str().c_str()));
					if (temp == NULL){
						cout << "hist not found = " << histname.str() << endl; 
						assert(false);
					}
					temp->SetDirectory(0);
					if (htbin ==0 && mhtbin==0) hist_temp = (TH1*) temp->Clone();
					else hist_temp->Add(temp);

				}
			}

			//if (h==0 || h == 1)
			if (hist_temp != 0) hist_temp->Write();
		}
	}
}

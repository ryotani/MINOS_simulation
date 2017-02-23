#include <TROOT.h>
#include <TChain.h>
#include <TTree.h>
#include <TBranch.h>
#include <TRandom.h>
#include <TFile.h>
#include <TMath.h>
#include <TF1.h>
#include <iostream>
#include <TH1F.h>
#include <TCanvas.h>
#include <TRint.h>
#include <TGraph.h>
#include <TObject.h>
#include "../lib/ExN03DataoRings.hh"
#include "../lib/ExN03Datai.hh"
#include "../lib/ExN03Setup.hh"
#include <stdio.h>
#define PI 3.14159 
double Qradt[20000];
bool PadTouch[20000];
TRandom Rand;
double conv_fit(double x[], double p[]);
double polya(double x[], double p[]);
double conv(double x[], double p[]);
void CalcTimeElectrons(ExN03DataoRings *datao, ExN03Setup *setup, int &j, TH1D *h_conv, TF1 *f_conv, TF1 *fit);

int main(int argc, char** argv)
{
   
   	TRint *theApp = new TRint("Rint",&argc,argv,0,0);  
	
	//TFile *MyFilei = new TFile(("/home/irfulx168/mnt/acorsi/MINOS_Bfield/result/SAMURAI/simu_Li11_250MeV.root"),"READ");
	TFile *MyFilei = new TFile(("/home/local1/workspace/result/simu_250MeV.root"),"READ");
	TTree *MyTreei = (TTree *) MyFilei->Get("tpc");
	ExN03Datai *tpc = new ExN03Datai;
	MyTreei->SetBranchAddress("tpc", &tpc);
	MyTreei->GetBranch("tpc")->SetFile("simu_Li11_250MeV.root");
	Int_t nevent = MyTreei->GetEntries();

	TFile *MyFileo = new TFile(("/home/local1/workspace/result/SAMURAI/simu_250MeV.root"),"RECREATE");
	TTree *MyTreeo = new TTree("MyTreeo","Realistic TPC events");
	ExN03DataoRings *tpco = new ExN03DataoRings;
	MyTreeo->Branch("tpco","ExN03DataoRings",&tpco);
	ExN03Setup *setup = new ExN03Setup;
	
	int NRings = tpco->NRings;
	int NSegments = tpco->NSegments;
	int NPads = NRings*NSegments;

	for(int k=0;k<NPads;k++)
	{
		PadTouch[k]=false;
	}
        TH1D *h_conv; 
	TF1 *f_conv; 
	TF1 *fit; 

	for (Int_t j=0;j<nevent;j++) 
	{
    		//if(j%10==0) 
		cout << "\n---> Begin of event: " << j << endl;
		int nb_pads=0;  	
  		MyTreei->GetEvent(j);
  		tpco->datai = *tpc;
  		if(tpco->datai.event == 5) CalcTimeElectrons(tpco, setup, nb_pads, h_conv, f_conv, fit); //add if(tpc

  		if(nb_pads>0) 
		{
			MyTreeo->Fill();
		}
		//cout<< "even details: "<<(*tpc).x_tpc.size() << " " << (*tpc).Et_tpc_tot << endl;
	}
	cout << "NbEvents = " << nevent << endl;

	MyFileo->Write();
	MyFileo->Close();
	theApp->Run(kTRUE);
	delete theApp;
	return 0;
}

void CalcTimeElectrons(ExN03DataoRings *datao, ExN03Setup *setup, int &nb_pads, TH1D *h_conv, TF1 *f_conv, TF1 *fit)
{
	int NPads = datao->NRings*datao->NSegments;
	TF1 *fonc_gain = new TF1("fonc_gain", polya, 0., 10000., 2);
	fonc_gain->FixParameter(0, datao->Gain);
	fonc_gain->FixParameter(1, datao->Theta);
	int rmm,thetamm;
	double x = 0.,y = 0.,z = 0.,r = 0.,dr = 0.;
	vector<int> pad;
	vector<double> time;
	datao->x_pad.clear();
	datao->y_pad.clear();
	datao->t_pad.clear();
	datao->q_pad.clear();
	datao->Et_pad=0;
	datao->nb_pads=0;
	double thr = datao->Threshold*datao->NoiseRMS;
	float Btheta=0.00;//10Gauss on Y:0.003;
	float Bphi=PI/2;
	for(int i=0;i<(int)(datao->datai.z_tpc.size());i++)
	{
		double sT=datao->sigT*TMath::Sqrt(datao->datai.z_tpc[i]/10.);
		double sL=datao->sigL*TMath::Sqrt(datao->datai.z_tpc[i]/10.);
		dr=(datao->LastRing-datao->FirstRing)/datao->NRings;
		double emm = 1./datao->Ionis*datao->datai.e_tpc[i];

		int nel=(int)(Rand.Gaus(emm,TMath::Sqrt(emm)));

		for(int k=0;k<nel;k++)
		{
			double theta=0.;
			x=Rand.Gaus(datao->datai.x_tpc[i],sT) + datao->datai.z_tpc[i]*tan(Btheta)*cos(Bphi);
			y=Rand.Gaus(datao->datai.y_tpc[i],sT) + datao->datai.z_tpc[i]*tan(Btheta)*sin(Bphi);
			z=Rand.Gaus(datao->datai.z_tpc[i],sL)/cos(Btheta);
//cout<<"test "<<datao->datai.z_tpc[i]*tan(Btheta)<<" "<<datao->datai.z_tpc[i]<<"; "<<Rand.Gaus(datao->datai.z_tpc[i],sL)/cos(Btheta)<< " "<<Rand.Gaus(datao->datai.z_tpc[i],sL)<<endl;
			double t_tpc=z/(datao->driftV);
			r=TMath::Sqrt(x*x+y*y);

			if(r>datao->FirstRing && r<datao->LastRing)
			{	
				rmm=(int)((r-datao->FirstRing)/dr);

				if(x>0&&y>0) theta=180./TMath::Pi()*TMath::ATan(y/x);
				if(x<0&&y>0) theta=180.-180./TMath::Pi()*TMath::ATan(-1.*y/x);
				if(x<0&&y<0) theta=180.+180./TMath::Pi()*TMath::ATan(y/x);
				if(x>0&&y<0) theta=360.-180./TMath::Pi()*TMath::ATan(-1.*y/x);
			
				if(theta < 0) theta = 360. + theta;
				thetamm=(int)(theta/(360./datao->NSegments));

				int npads = 0;
				for(int j=0; j<rmm; j++)
				{
					npads = npads + datao->NSegments;
				}
				int num_pad = npads + thetamm; 
				pad.push_back(num_pad);
				time.push_back(t_tpc);
				PadTouch[npads + thetamm] = true;
			}
		}
	}

	const int max_time= (int) (setup->ChamberLength / datao->driftV + 6.*datao->ShapingTime);
	const int samples = (int) ((setup->ChamberLength / datao->driftV + 6.*datao->ShapingTime)/datao->TimeBinSize);
	nb_pads=0;
	double noise[max_time];
	double tab_noise[max_time];
		bool tst = false;
	for(int i=0; i<NPads; i++)
	{
		/*TObject *obj = NULL;
		if((obj = gROOT->FindObject("noise")) != NULL) delete obj;  
		if((obj = gROOT->FindObject("After shaping")) != NULL) delete obj;  */
		for(int ii=0; ii<max_time; ii++)
		{
			noise[ii] = Rand.Gaus(0., datao->NoiseRMS);
		}
		for(int ii=0; ii<samples; ii++)
		{
			tab_noise[ii] = noise[(int) (ii*datao->TimeBinSize)];
		}
		double max_noise = tab_noise[TMath::LocMax(samples, tab_noise)];
		double t_max_noise = TMath::LocMax(samples, tab_noise)*datao->TimeBinSize;
		char s[10];
		sprintf(s, "%d", i);
		int nb_param=0;

		for(unsigned int j=0; j<pad.size(); j++)
		{
			if(pad[j] == i) nb_param+=2;
		}
		double p2[nb_param+2];
		p2[0] = nb_param+2;
		p2[1] = datao->ShapingTime;
		int ii=0;
		int kk=0;
		double telec[1000000];

		for(unsigned int j=0; j<pad.size(); j++)
		{
			if(pad[j] == i)
			{
				double val_gain = fonc_gain->GetRandom();
				double charge = val_gain;
				p2[2*ii+2+1] = time[j];
				p2[2*ii+2] = charge;
				ii++;
				

				telec[kk] = time[j];
				kk++;

		
			}
		}
		// temps barycentre
		double t0 = 0.;
		for(int ii=0; ii<kk; ii++)
		{
			t0 = t0 + telec[ii];
		}
		t0 = t0 / kk;
		
		//TH1D *
		h_conv = new TH1D("After shaping", (string(s) + " after shaping").c_str(), samples, 0., max_time);
		//TF1 *
		f_conv = new TF1("F_conv", conv, 0.,setup->ChamberLength / datao->driftV + 6.*datao->ShapingTime, int(p2[0]));
		for(int ii=0; ii<p2[0]; ii++)
		{
			f_conv->FixParameter(ii, p2[ii]);
		}
		
		if(PadTouch[i])
		{ 		
			for(int jj=0; jj<samples; jj++)
			{
				h_conv->Fill(jj*datao->TimeBinSize,(f_conv->Eval(jj*datao->TimeBinSize)) + noise[(int) (jj*datao->TimeBinSize)]);
			}
		}
		double t_rec=0.;
		double max, time0;
		if(PadTouch[i])
		{ 		
			max = h_conv->GetBinContent(h_conv->GetMaximumBin());
			time0 = h_conv->GetMaximumBin()*datao->TimeBinSize-datao->TimeBinSize/2.;
		}
		else
		{
			max = max_noise;
			time0 = t_max_noise;
		}
		//TF1 *
		fit = new TF1("FIT", conv_fit, 0., setup->ChamberLength / datao->driftV + 6.*datao->ShapingTime, 4);
		fit->SetParLimits(0, 0., 100000000);
		fit->SetParameter(0, max);
		fit->SetParLimits(1, time0-datao->ShapingTime*2., time0);
		fit->SetParameter(1, time0-datao->ShapingTime);
		fit->FixParameter(2, p2[1]);
		fit->SetParameter(3, 10.);
		fit->SetNpx(100000);
		if(PadTouch[i]) h_conv->Fit(fit, "NOQ"); //Q=quiet
		fit->SetLineColor(2);
		t_rec = fit->GetParameter(1);
		double integral = fit->GetMaximum(0., max_time);
			
		int npads = 0;
		int val_ring = 0, val_seg = 0;
		for(int k=0; k<datao->NRings; k++)
		{
			double v = (double)(i) / (datao->NSegments+npads);	
			if(v>1) npads += datao->NSegments;
			else
			{
				val_ring = k;
				k = datao->NRings;
				val_seg = i - npads;
			}
		}
		if(PadTouch[i])
		{
			nb_pads++;
			/*if(!tst){
			TCanvas *c = new TCanvas("can","",800,800);
			h_conv->DrawCopy();
			fit->DrawCopy("same");
			tst=true;
			} */
		}
		double rr=datao->FirstRing+((double)(val_ring)+0.5)*dr;
		double cost=TMath::Cos((double)(val_seg+0.5)/datao->NSegments*2.*TMath::Pi());
		double sint=TMath::Sin((double)(val_seg+0.5)/datao->NSegments*2.*TMath::Pi());
		if(max > thr) 
		{
			datao->x_pad.push_back(rr*cost);
			datao->y_pad.push_back(rr*sint);
			datao->t_pad.push_back(t_rec);
			datao->q_pad.push_back(integral);
			datao->Et_pad+=integral;
			datao->dt.push_back(t_rec-t0);
		}
		delete h_conv, f_conv, fit;

	}
	for(int i=0; i<NPads; i++)
	{
		PadTouch[i] = false;
		Qradt[i] = 0.;
	}
	datao->nb_pads=nb_pads;
}


double conv(double x[], double p[])
{
	double val=0.;
	for(int i=0; i<(p[0]-2)/2; i++)
	{
		if(!(x[0]<p[2*i+2+1] || x[0]>1000000.)) val += p[2*i+2] * 22.68113723 * exp(-3.*(x[0]-p[2*i+2+1])/p[1]) * sin((x[0]-p[2*i+2+1])/p[1]) * pow((x[0]-p[2*i+2+1])/p[1], 3);
	}
	return(val);

}
double polya(double x[], double p[])
{
	double val = (pow(p[1]+1, p[1]+1) / TMath::Gamma(p[1]+1)) *
	pow(x[0]/p[0], p[1]) * exp(-(p[1]+1) * x[0] / p[0]);
	return val;
}
double conv_fit(double x[], double p[])
{
	double val=0.;
	if(!(x[0]<p[1] || x[0]>1000000.)) val += p[0] * 22.68113723 * exp(-3.*(x[0]-p[1])/p[2]) * sin((x[0]-p[1])/p[2]) * pow((x[0]-p[1])/p[2], 3) + p[3];
	else val += p[3];
	return(val);
}

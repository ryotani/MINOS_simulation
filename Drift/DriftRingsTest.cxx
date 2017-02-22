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
#include <TMinuit.h>
#include "/home/user/MINOS_Bfield/lib/ExN03DataoRings.hh"
#include "/home/user/MINOS_Bfield/lib/ExN03Datai.hh"
#include "/home/user/MINOS_Bfield/lib/ExN03Setup.hh"
#include <stdio.h>
#include <time.h>

bool PadTouch[20000];
TRandom Rand;
double conv_fit(double x[], double p[]);
double polya(double x[], double p[]);
double conv(double x[], double p[]);
void CalcTimeElectrons(ExN03DataoRings *datao, ExN03Setup *setup, int &j);
int main(int argc, char** argv)
{
        time_t start, stop;
   	time(&start);
   	TRint *theApp = new TRint("Rint",&argc,argv,0,0);  
   
	TFile *MyFilei = new TFile(("/home/user/MINOS_Bfield/result/(p,2p)/simuB4T_theta55_energy60.root"),"READ");
	TTree *MyTreei = (TTree *) MyFilei->Get("tpc");
	ExN03Datai *tpc = new ExN03Datai;
	MyTreei->SetBranchAddress("tpc", &tpc);
	MyTreei->GetBranch("tpc")->SetFile("simuB4T_theta55_energy60.root");
	Int_t nevent = MyTreei->GetEntries();

	TFile *MyFileo = new TFile(("/home/user/MINOS_Bfield/Drift/result/(p,2p)/simuB4T_theta55_energy60.root"),"RECREATE");
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
//nevent=20;
	int nb_pads=0;  	
	for (Int_t j=0;j<nevent;j++) 
	{
    		if(j%10==0) cout << "\n---> Begin of event: " << j << endl;
  		MyTreei->GetEvent(j);
  		tpco->datai = *tpc;
  		CalcTimeElectrons(tpco, setup, nb_pads); 

  		if(nb_pads>0) 
		{
			MyTreeo->Fill();
		}
	}
	cout << "NbEvents = " << nevent << endl;

	MyFileo->Write();
	MyFileo->Close();
   	time(&stop);
   	cerr << "Running time: " << stop-start << " s" << endl;
	theApp->Run(kTRUE);
	delete theApp;
	return 0;
}

void CalcTimeElectrons(ExN03DataoRings *datao,ExN03Setup *setup, int &nb_pads)
{
	int NPads = datao->NRings*datao->NSegments;
	TF1 *fonc_gain = new TF1("fonc_gain", polya, 0., 10000., 2);
	fonc_gain->FixParameter(0, datao->Gain);
	fonc_gain->FixParameter(1, datao->Theta);
	int rmm=0,thetamm=0,nel=0,npads=0,num_pad=0;
	double x = 0.,y = 0.,z = 0.,r = 0.,dr = 0., sT=0., sL=0., emm=0., theta=0., t_tpc=0.;
	vector<int> pad;
	vector<double> time;
	//vector<bool> PadTouch;
	
	datao->x_pad.clear();
	datao->y_pad.clear();
	datao->t_pad.clear();
	datao->q_pad.clear();
	datao->Et_pad=0;
	double thr = datao->Threshold*datao->NoiseRMS;
	for(int i=0;i<(int)(datao->datai.z_tpc.size());i++)
	{
		sT=datao->sigT*TMath::Sqrt(datao->datai.z_tpc[i]/10.);
		sL=datao->sigL*TMath::Sqrt(datao->datai.z_tpc[i]/10.);
		dr=(datao->LastRing-datao->FirstRing)/datao->NRings;
		emm = 1./datao->Ionis*datao->datai.e_tpc[i];

		nel=(int)(Rand.Gaus(emm,TMath::Sqrt(emm)));
		for(int k=0;k<nel;k++)
		{
			theta=0.;
			x=Rand.Gaus(datao->datai.x_tpc[i],sT);
			y=Rand.Gaus(datao->datai.y_tpc[i],sT);
			z=Rand.Gaus(datao->datai.z_tpc[i],sL);

			t_tpc=z/(datao->driftV);
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

				npads = 0;
				for(int j=0; j<rmm; j++)
				{
					npads = npads + datao->NSegments;
				}
				num_pad = npads + thetamm; 
				pad.push_back(num_pad);
				time.push_back(t_tpc);
				PadTouch[npads + thetamm] = true;
			}
		}
	}

	const int max_time= (int) (setup->ChamberLength / datao->driftV + 6.*datao->ShapingTime);
	const int max_time_noise = (int) (setup->ChamberLength / datao->driftV);
	const int samples = (int) ((setup->ChamberLength / datao->driftV + 6.*datao->ShapingTime)/datao->TimeBinSize);
	nb_pads=0;
	vector<double> noise;
	double val_noise=0., max_noise=0.;
	for(int ii=0; ii<max_time; ii++)
	{
		val_noise = Rand.Gaus(0., datao->NoiseRMS);
		noise.push_back(val_noise);
		if (val_noise>max_noise) max_noise = val_noise;
	}
	int nb_param=0,val_ring=0,val_seg=0;
	double charge=0.,t0=0.,t_rec=0.,time0=0.,integral=0.,rr=0.,cost=0.,sint=0.,v=0.,max=0.; 
	
	TF1 *fit = new TF1("FIT", conv_fit, 0., setup->ChamberLength / datao->driftV + 6.*datao->ShapingTime, 4);
	fit->SetParLimits(0, 0., 100000000);
	fit->FixParameter(2, datao->ShapingTime);
	vector<double> p2;
	vector<double> telec;
	
	for(int i=0; i<NPads; i++)
	{
		nb_param=0;
		for(unsigned int j=0; j<pad.size(); j++)
		{
			if(pad[j] == i) nb_param+=2;
		}

		for(unsigned int j=0; j<pad.size(); j++)
		{
			if(pad[j] == i)
			{
				charge = fonc_gain->GetRandom();
				p2.push_back(charge);
				p2.push_back(time[j]);				
				telec.push_back(time[j]);

		
			}
		}
		// temps barycentre
		t0 = 0.;
		for(unsigned int ii=0; ii<telec.size(); ii++)
		{
			t0 = t0 + telec[ii];
		}
		t0 = t0 / telec.size();
		
		TH1D *h_conv = new TH1D("After shaping", "", samples, 0., max_time);
		TF1 *f_conv = new TF1("F_conv", conv, 0.,setup->ChamberLength / datao->driftV + 6.*datao->ShapingTime, p2.size()+2);
		for(unsigned int ii=2; ii<p2.size(); ii++)
		{
			f_conv->FixParameter(ii, p2[ii]);
		}
		f_conv->FixParameter(0, p2.size()+2);
		f_conv->FixParameter(1, datao->ShapingTime);
		if(PadTouch[i])
		{ 		
			for(int jj=0; jj<samples; jj++)
			{
				h_conv->Fill(jj*datao->TimeBinSize,(f_conv->Eval(jj*datao->TimeBinSize)) + noise[Rand.Integer(max_time)]);
			}
			max = h_conv->GetBinContent(h_conv->GetMaximumBin());
			time0 = h_conv->GetMaximumBin()*datao->TimeBinSize-datao->TimeBinSize/2.;
			fit->SetParameter(0, max);
			fit->SetParLimits(1, time0-datao->ShapingTime*2., time0);
			fit->SetParameter(1, time0-datao->ShapingTime);
			fit->SetParameter(3, 0.);
			h_conv->Fit(fit, "NOQ");
			t_rec = fit->GetParameter(1);
			integral = fit->GetParameter(0);
			nb_pads++; 
		}
		else
		{
			integral = max_noise;
			t_rec = Rand.Uniform(max_time_noise);
		}
		npads = 0;
		for(int k=0; k<datao->NRings; k++)
		{
			v = (double)(i) / (datao->NSegments+npads);	
			if(v>1) npads += datao->NSegments;
			else
			{
				val_ring = k;
				k = datao->NRings;
				val_seg = i - npads;
			}
		}
		rr=datao->FirstRing+((double)(val_ring)+0.5)*dr;
		cost=TMath::Cos((double)(val_seg+0.5)/datao->NSegments*2.*TMath::Pi());
		sint=TMath::Sin((double)(val_seg+0.5)/datao->NSegments*2.*TMath::Pi());
		if(integral > thr) 
		{
			datao->x_pad.push_back(rr*cost);
			datao->y_pad.push_back(rr*sint);
			datao->t_pad.push_back(t_rec);
			datao->q_pad.push_back(integral);
			datao->Et_pad+=integral;
		}
		delete h_conv;
		delete f_conv;
		p2.clear();
		telec.clear();
	}
	for(int i=0; i<NPads; i++)
	{
		PadTouch[i] = false;
	}
	delete fit;
	delete fonc_gain;
	pad.clear();
	time.clear();
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
	if(!(x[0]<p[1])) val += p[0] * 22.68113723 * exp(-3.*(x[0]-p[1])/p[2]) * sin((x[0]-p[1])/p[2]) * pow((x[0]-p[1])/p[2], 3) + p[3];
	else val += p[3];
	return(val);
}

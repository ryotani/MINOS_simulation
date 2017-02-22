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
#include "lib/TMinosTPCSimData.h"
#include "lib/ExN03DataoRings.hh"
#include "lib/ExN03Datai.hh"
#include <stdio.h>
#include <TClonesArray.h>


void convtree(char *filename)
{
	
    TFile *f = new TFile(filename,"READ");//input from simulation
		TTree *MyTreeo = (TTree *) f->Get("MyTreeo");		

		ExN03DataoRings *tpc = new ExN03DataoRings;
		MyTreeo->SetBranchAddress("tpco", &tpc);
		MyTreeo->GetBranch("tpco");
		
		
		Int_t nevent = MyTreeo->GetEntries();
		cout<<"----------- MyTreeo has "<<nevent<<" entries"<<endl;

    		// ROOT output file (locate number of run, time & date)    
    		TFile *fout = new TFile("results/test_out_noB.root","RECREATE");
    		TTree *tree = new TTree("tree","Converted tree");

		TClonesArray TPCData;
    		TPCData.SetClass("TMinosTPCSimData");
    		tree->Branch("TPCData",&TPCData);
		double x0, y0, z0;
		tree->Branch("x0",&x0);
		tree->Branch("y0",&y0);
		tree->Branch("z0",&z0);
		double xB, yB, zB, thetaB,phiB;
		tree->Branch("xB",&xB);
		tree->Branch("yB",&yB);
		tree->Branch("zB",&zB);
		tree->Branch("thetaB",&thetaB);
		tree->Branch("phiB",&phiB);
		double vdrift;
		tree->Branch("vdrift",&vdrift);


		for (Int_t j=0;j<nevent;j++) //loop over events
		{
			
   			TPCData.Clear();
			TMinosTPCSimData *minosdataTPC;

			MyTreeo->GetEvent(j); //cout<<"entry "<<j<<" "<<tpc->x_pad.size()<<endl;
			
			for(unsigned int i=0; i<tpc->x_pad.size(); i++)
			{
				minosdataTPC = (TMinosTPCSimData*)TPCData.ConstructedAt(i);
				minosdataTPC->Set(tpc->x_pad[i], tpc->y_pad[i], tpc->t_pad[i], tpc->dt[i], tpc->q_pad[i]);
				
			}
			//cout<<tpc->datai.x0.size()<<endl;
			if(tpc->datai.x0.size()>1)
			{
			xB = tpc->datai.x0[0];
			yB = tpc->datai.y0[0];
			zB = tpc->datai.z0[0];
			thetaB = tpc->datai.theta0[0];
			phiB = tpc->datai.phi0[0];
	
			x0 = tpc->datai.x0[1];
			y0 = tpc->datai.y0[1];
			z0 = tpc->datai.z0[1];
			}
			vdrift = tpc->driftV;

			tree->Fill();
			
		}
				
		MyTreeo->Delete();
		f->Close();
		
		if(tree->GetEntries()!=0) {fout->Write();}
//		tree->Delete();
//		fout->Close();

}

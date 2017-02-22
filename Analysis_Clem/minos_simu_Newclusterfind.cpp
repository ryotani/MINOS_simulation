// minos_simu_clusterfind_proj.cpp
// ->Makefile_simu_cluster_proj
//
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Cluster finding !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//
// Output file: Ryr_mth_d-h_min_s-run.root
//
//

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! FOR SEVERAL FEMINOS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//                                                        (w/ new DAQ)
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! DO NOT USE THIS FILE IF FEMINOS WITH AGET CHIPS IS USED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//run with: ./minos_simu_clusterfind_proj ../Drift/simu***.root

#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <TString.h>
// root headers
#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TRint.h>
#include <TChain.h>
#include <TBranch.h>
#include <TRandom.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TF2.h>
#include <TH3.h>
#include <TF3.h>
#include <TRint.h>
#include <time.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TMath.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TMarker3DBox.h>
#include <TPolyLine3D.h>
#include <TClonesArray.h>
#include <arpa/inet.h>

// home-made headers
#include "lib/TMinosTPCSimData.h"
#include "lib/TMinosClust.h"


#define FILEDEF "test.acq"

double PI = TMath::Pi();
using namespace std;
using namespace ROOT::Math;


class Comp
{
public:
    Comp(vector<double>& inVec): _V(inVec) {}
    bool operator()(int i, int j) {return (_V.at(i)>_V.at(j));print(i,j);}
    void print(int i, int j){std::cout<<"Cmp "<<_V.at(i)<<" "<<_V.at(j)<<std::endl;}
private:
    vector<double> _V;
};

//int Obertelli_filter(TCanvas *c1, int &DiagTaken, int entry, int iteration, vector<double> *x, vector<double> *y,vector<double> *z,vector<double> *q, vector<double> *x_out,vector<double> *y_out,vector<double> *z_out, vector<double> *q_out, vector<int> *ringbool)
int Obertelli_filter(int &DiagTaken, int entry, int iteration, vector<double> *x,vector<double> *y,vector<double> *z,vector<double> *q, vector<double> *x_out,vector<double> *y_out,vector<double> *z_out, vector<double> *q_out, vector<int> *ringbool)
{
	double bint1=2.;
	double bint2=2.;
	int maxt = 360.;
	int mint = 0.;
	int nt1=(maxt-mint)/bint1;
	int nt2=(maxt-mint)/bint1;
    
	double Rint = 53;
	double Rext = 53 + 72;
	int filter_result = 0;
	
	TH2F *hp_xy = new TH2F("hp_xy",Form("hp_xy_Evt%d_cluster%d",entry, iteration),nt1,mint,maxt,nt2,mint,maxt);
	TH2F *hpDiag_xy = new TH2F("hpDiag_xy",Form("hpDiag_xy_Evt%d_cluster%d",entry, iteration),nt1,mint,maxt,nt2,mint,maxt);
	//TH2F *hc_xy = new TH2F("hc_xy",Form("hc_xy_Evt%d_cluster%d",entry, iteration),100,-85,85,100,-85,85);
	//TH2F *hcnew_xy = new TH2F("hcnew_xy",Form("hcnew_xy_Evt%d_cluster%d",entry, iteration),100,-85,85,100,-85,85);

	double max_xy;
	TF1* line_xy = new TF1("line_xy","[0] + [1]*x",-85,85);

	vector<double> xTemp, yTemp, zTemp, qTemp;

	double theta1, theta2, xt1, yt1, xt2, yt2;
	double line0=0., line1=0.;
	double delta=0., AA=0., BB=0., CC=0.;
	double maxtheta1=0., maxtheta2=0., xmax1=0., ymax1=0., xmax2=0., ymax2=0.;
	double par0=0., par1=0.;
	double r_mm=0.;
	int ringsum=0;
	bool maxfound = false;
	
	for(unsigned int i=0;i<x->size();i++)
	{

		xTemp.push_back(x->at(i));
		yTemp.push_back(y->at(i));
		zTemp.push_back(z->at(i));
		qTemp.push_back(q->at(i));

		//Fill coordinate space histograms for plots
		//hc_xy->Fill(x->at(i),y->at(i),q->at(i));
 
		//Loop of indices
		for(int j=0;j<nt1;j++)
		{

			theta1 = (j+0.5)*bint1 + mint;
			xt1 = Rint * TMath::Cos(theta1*PI/180.);
			yt1 = Rint * TMath::Sin(theta1*PI/180.);
			line1 = (yt1 - y->at(i))/(xt1 - x->at(i));
			line0 = yt1 - xt1 * line1;
			AA = 1 + line1*line1;
			BB = 2*line0*line1;
			CC = line0*line0 - Rext*Rext;

			delta = BB*BB - 4*AA*CC;
			if(delta>=0)
			{
				xt2 = (-BB - sqrt(delta))/(2*AA);
				yt2 = line0 + line1*xt2;
                    		if(xt2<=0)	      theta2=  180 - asin(yt2/Rext)*180/PI;
                    		else if(xt2>0)
				{
					if(yt2>0)	      theta2=  asin(yt2/Rext)*180/PI;
                    			else if(yt2<=0)	      theta2=  360 + asin(yt2/Rext)*180/PI;
				}
				

				if( (xt1*x->at(i) + yt1*y->at(i))>=0 && (xt2*x->at(i) + yt2*y->at(i))>=0  && (xt1*xt2+yt1*yt2)>=0)
				{
					hp_xy->Fill(theta1,theta2);
					if(abs(theta1-theta2)<=10) hpDiag_xy->Fill(theta1,theta2);
				}
				else
				{
					if(delta!=0)
					{
						xt2 = (-BB + sqrt(delta))/(2*AA);
						yt2 = line0 + line1*xt2;
                    				if(xt2<=0)	      theta2=  180 - asin(yt2/Rext)*180/PI;
                    				else if(xt2>0)
						{
							if(yt2>0)	      theta2=  asin(yt2/Rext)*180/PI;
                    					else if(yt2<=0)	      theta2=  360 + asin(yt2/Rext)*180/PI;
						}
						if( (xt1*x->at(i) + yt1*y->at(i))>=0 && (xt2*x->at(i) + yt2*y->at(i))>=0  && (xt1*xt2+yt1*yt2)>=0)
						{
							hp_xy->Fill(theta1,theta2);
							if(abs(theta1-theta2)<=10) hpDiag_xy->Fill(theta1,theta2);
						}
					}
				}

			}
			
		}
	}

	x->clear();
	y->clear();
	q->clear();
	z->clear();

	if(hpDiag_xy->GetMaximum()>=10) 
	{
		max_xy = hpDiag_xy->GetMaximum();
		DiagTaken = 1;
//		cout << "Max taken in diag... with value=" << max_xy << endl;
	}
	else 
	{
		max_xy = hp_xy->GetMaximum();
		DiagTaken = 0;
	}
	
	for(int ii=0; ii<nt1; ii++)
	{
		if(maxfound ==true) break;
		for(int jj=0; jj<nt2; jj++)
		{
			if(hp_xy->GetBinContent(ii+1, jj+1) == max_xy)
			{
				maxtheta1 = (ii+0.5)*bint1 + mint;
				maxtheta2 = (jj+0.5)*bint2 + mint;
				maxfound = true;
				cout << "xy: theta max are " << maxtheta1 << " , " << maxtheta2 << endl;

			}
			if(maxfound ==true) break;
		}
	}

	xmax1 = Rint * TMath::Cos(maxtheta1*PI/180.);
	ymax1 = Rint * TMath::Sin(maxtheta1*PI/180.);
	xmax2 = Rext * TMath::Cos(maxtheta2*PI/180.);
	ymax2 = Rext * TMath::Sin(maxtheta2*PI/180.);

	// xy PEAK
	par1 = (ymax2-ymax1)/(xmax2-xmax1);
	par0 = (ymax1 - xmax1*par1);
/*
	cout<<"xmax1 "<<xmax1<<" ymax1 "<<ymax1<<" xmax2 "<<xmax2<<" ymax2 "<<ymax2<<endl;
        line_xy->SetParameter(0,par0);
        line_xy->SetParameter(1,par1);
        hc_xy->GetListOfFunctions()->Add(line_xy);
	line_xy->SetLineWidth(1);
*/


	//Selection of x,y points IN the maxmean+/-1 found in Obertelli transform of xy plane
	for(unsigned int i=0;i<xTemp.size();i++)
	{
		if( (abs(par1*xTemp[i]-yTemp[i]+par0)/sqrt(1+par1*par1))<= 6 && ((xmax1*xTemp[i] + ymax1*yTemp[i]) >= 0) && ((xmax2*xTemp[i] + ymax2*yTemp[i]) >= 0) && ((xmax1*xmax2 + ymax1*ymax2) >= 0))
		{
//			cout << "Taken points= " << xTemp[i] << " , " << yTemp[i] << " , " << zTemp[i] << endl;
			//hcnew_xy->Fill(xTemp[i],yTemp[i],qTemp[i]);
			x_out->push_back(xTemp[i]);
			y_out->push_back(yTemp[i]);
			z_out->push_back(zTemp[i]);
			q_out->push_back(qTemp[i]);
			//cerr << "Taken points= " << xTemp[i] << " , " << yTemp[i] << " , cluster = " << iteration +1 << endl;
			filter_result++;
			r_mm = sqrt(xTemp[i]*xTemp[i]+yTemp[i]*yTemp[i]);
			if(r_mm<(52+5*3)) ringsum++;
		}
		else
		{	
			x->push_back(xTemp[i]);
			y->push_back(yTemp[i]);
			z->push_back(zTemp[i]);
			q->push_back(qTemp[i]);

		}
		
	}

	for(int ip=0; ip<filter_result; ip++)
	{
		if(ringsum>2) ringbool->push_back(1);
		else ringbool->push_back(0);
	}
/*
	c1->Divide(4,1);
	// Coordinate space
	c1->cd(1);
	hc_xy->Draw("colz");
	// Hough space
	c1->cd(2);
	hp_xy->Draw("colz"); 
	c1->cd(3);
	hpDiag_xy->Draw("colz"); 
	// Coordinate space : New plot
	c1->cd(4);
	hcnew_xy->Draw("colz");

	c1->Update();
*/
	delete hp_xy;
	delete hpDiag_xy;

	return filter_result;

}


int main(int argc, char** argv)
{
    // Definition of run configuration (electronics, TPC)
    double thresh = 4000.;

    // Variables to calculate the elapsed time of the process
    time_t start,stop;


    // Filename
    char* filename;
    filename=FILEDEF;
    if(argc < 2)
    {
        cerr << "Missing file argument" << endl;
    }
    
    filename = argv[1];
    
    // Definition of input file (put before the def of root file that will come out as filename changes)
    cout << "Input File = " << filename << endl;
    
    ifstream file(filename,ifstream::in);
    if(!file.is_open())
    {
        cout<<" File "<< filename <<" does not exist..."<<endl;
        exit(1);
    }
    
    TFile *f = new TFile(filename);
    TClonesArray *TPCData = new TClonesArray("TMinosTPCSimData");
    TTree *tree = (TTree*)f->Get("tree");
    tree->GetBranch("TPCData");
    tree->SetBranchAddress("TPCData",&TPCData);
    double x0, y0, z0;
    tree->GetBranch("x0");
    tree->SetBranchAddress("x0",&x0);
    tree->GetBranch("y0");
    tree->SetBranchAddress("y0",&y0);
    tree->GetBranch("z0");
    tree->SetBranchAddress("z0",&z0);
    double xB, yB, zB, thetaB, phiB;
    tree->GetBranch("xB");
    tree->SetBranchAddress("xB",&xB);
    tree->GetBranch("yB");
    tree->SetBranchAddress("yB",&yB);
    tree->GetBranch("zB");
    tree->SetBranchAddress("zB",&zB);
    tree->GetBranch("thetaB");
    tree->SetBranchAddress("thetaB",&thetaB);
    tree->GetBranch("phiB");
    tree->SetBranchAddress("phiB",&phiB);
    double vdrift;
    tree->GetBranch("vdrift");
    tree->SetBranchAddress("vdrift",&vdrift);

    int nentries = tree->GetEntries(); // number of events
    cout << "nentries = " << nentries << endl;
    
    // ROOT output file (locate number of run, time & date)
    string pchroot= filename;
    string pchtxt= filename;
    string pch2 = ".root";
    pchroot.replace(pchroot.find(pch2),pch2.length(),Form("_Obert-thr%d.root",int(thresh)));
    pchtxt.replace(pchtxt.find(pch2),pch2.length(),Form("_Obert-thr%d.txt",int(thresh)));
 
    const char *OutFile = pchroot.c_str();
    const char *OutTextFile = pchtxt.c_str();
    freopen(OutTextFile,"w",stdout);
    

    cerr << "OutFile: " << OutFile << " " << OutTextFile << endl;    

    TFile *fout = new TFile(OutFile,"RECREATE");
    TTree *tree_out = new TTree("tree_out","analyzed tree");
    TClonesArray data_out;
    data_out.SetClass("TMinosClust");
    tree_out->Branch("data_out",&data_out);
    int evtOrig;
    tree_out->Branch("evtOrig",&evtOrig,"evtOrig/I");
    int trackNbr;
    tree_out->Branch("trackNbr",&trackNbr,"trackNbr/I");
    double x_int, y_int, z_int;
    tree_out->Branch("x0",&x_int,"x0/D");
    tree_out->Branch("y0",&y_int,"y0/D");
    tree_out->Branch("z0",&z_int,"z0/D");
    double x_beam, y_beam, z_beam, theta_beam, phi_beam;
    tree_out->Branch("xB",&x_beam,"xB/D");
    tree_out->Branch("yB",&y_beam,"yB/D");
    tree_out->Branch("zB",&z_beam,"zB/D");
    tree_out->Branch("thetaB",&theta_beam,"thetaB/D");
    tree_out->Branch("phiB",&phi_beam,"phiB/D");
    
    time(&start);

    // Definition of variables    
    double x_mm,y_mm,z_mm,r_mm,q_pad,t_pad;
    vector<double> Xpad, Ypad, Qpad, Zpad, radius, index;
    vector<double> XpadNew, YpadNew, ZpadNew, QpadNew;
    vector<int> clusternbr;
    vector<int> clusterpads;
    vector<int> clusterringbool;
    int Iteration=0;
    int filter_result=0;

    //vector<TCanvas*> Filter_canvas;
    int max_diag;

    int filled=0;
    int goodevts=0;

    
    // Start of loop over events    
//    for(int i=0; i<50; i++)
    for(int i=0; i<nentries; i++) // loop over events
    {
        
        ///////////////////////////////////////////////////First reading=>put in vectors/////////////////////////////////////////////////////
        
        // Clear event
	//tpco->Clear();
        data_out.Clear();
        Xpad.clear();
        Ypad.clear();
        Zpad.clear();
        Qpad.clear();
        XpadNew.clear();
        YpadNew.clear();
        ZpadNew.clear();
        QpadNew.clear();
        clusternbr.clear();
        clusterpads.clear();
        clusterringbool.clear();
        radius.clear();
        index.clear();
        
        filled = 0;
	trackNbr=0;
	Iteration=0;
	filter_result=0;

        //TH2F *padpattern0_clone = (TH2F*)padpattern0->Clone(Form("padpattern0_nev%d",i));
        
        tree->GetEntry(i);
        if(i%1000==0)cerr <<  "-----ENTRY " << i << " / " << nentries << "-----" << endl;
        cout << "-----ENTRY " << i << " / " << nentries <<"-----" << endl;

        evtOrig = i;
	x_int = x0;
	y_int = y0;
	z_int = z0;

	x_beam = xB;
	y_beam = yB;
	z_beam = zB;
	theta_beam = thetaB;
	phi_beam = phiB;
        
        TMinosTPCSimData *minosTPCData;
        TMinosClust *minosdata_out;
        
        //build vectors containing x,y,z,radius and index of entry for all events
        for (int hitcount=0; hitcount<TPCData->GetEntriesFast(); hitcount++)
        {
        	minosTPCData = (TMinosTPCSimData*)TPCData->At(hitcount);
	    	x_mm = minosTPCData->x_mm;
	    	y_mm = minosTPCData->y_mm;
	    	t_pad = minosTPCData->t_pad;
		z_mm = (t_pad*vdrift); //time bin * vdrift(mm/ns)
	    	q_pad = minosTPCData->q_pad;

            	if((q_pad-250)>=thresh)
		{
                    r_mm=sqrt(x_mm*x_mm+y_mm*y_mm);
                    
                    radius.push_back(r_mm);
                    Xpad.push_back(x_mm);
                    Ypad.push_back(y_mm);
		    Zpad.push_back(z_mm);
                    Qpad.push_back(q_pad-250);
                    index.push_back(filled);
                    filled++;
                }
            
        }//END of entries in tclonesarray for the entry
        
        
        
	///////////////////////////////////////////////// Find clusters :: HOUGH ////////////////////////////////////////////////////
	if(filled>0)
	{
		while(Xpad.size()>=10 && Iteration<20)
		{
			filter_result = 0;
			Iteration++;
  			//Filter_canvas.push_back(new TCanvas(Form("Event%d_cluster%d", i, Iteration), Form("Event%d_cluster%d", i, Iteration)));
//			filter_result = Obertelli_filter(Filter_canvas.back(),max_diag, i, Iteration, &Xpad, &Ypad, &Zpad, &Qpad, &XpadNew, &YpadNew, &ZpadNew, &QpadNew, &clusterringbool);
			filter_result = Obertelli_filter(max_diag, i, Iteration,&Xpad, &Ypad, &Zpad, &Qpad, &XpadNew, &YpadNew, &ZpadNew, &QpadNew, &clusterringbool);
			cout << "Evt" << i << " ::: Obertelli Cluster n°" << Iteration << " of size=" << filter_result << ", max taken in Diag=" << max_diag << ", new evt size=" << XpadNew.size() << " clusteringbool "<<clusterringbool.back()<<endl; 
			//Filter_canvas.back()->Write();

			for(int ik=0; ik<filter_result; ik++)
			{
				clusterpads.push_back(filter_result);
				clusternbr.push_back(Iteration);
			}

			if(filter_result>10 && clusterringbool.back()==1) trackNbr++;
		}
        
		if(trackNbr==1 || trackNbr==2 || trackNbr==3 || trackNbr==4) goodevts++;
			
		for(unsigned int il=0; il<XpadNew.size(); il++)
		{
				minosdata_out = (TMinosClust*)data_out.ConstructedAt(il);
                		minosdata_out->Set(XpadNew[il], YpadNew[il], ZpadNew[il] , QpadNew[il], clusternbr[il], clusterpads[il], clusterringbool[il]);
				//cerr << "Taken points= " << XpadNew[il] << " , " << YpadNew[il] << " , cluster = " << clusternbr[il] << endl;
		}
            
	    	tree_out->Fill();

	    	cout << "trackNbr=" << trackNbr << endl;
	   	//cerr << "same size ? " << Phi.size() << " " << Xpad.size() << " " << Ypad.size() << " " << clusternbr.size() << " " << clusterpads.size() << endl;
            


        }//end of if filled
        
        
    }//end of LOOP on entries
    
    tree_out->Print();
    fout->Write();
    f->Close();
    fout->Close();
    
    time(&stop);
    cerr << "Evts with 1/2/3/4p tracks : " << goodevts << " / " << nentries << endl;
    printf("Elapsed time: %.1f seconds\n",difftime(stop,start));
    
    return 0;
    
}



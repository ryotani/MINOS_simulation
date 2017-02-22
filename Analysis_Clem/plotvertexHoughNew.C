// plotvertexHoughNew.C
// 
// Created 13/11/15
//
//  root[0]> .L plotvertexHoughNew.C++
//  root[1]> plotvertexHoughNew("filepath/filename.root")
//

// To reconstruct the vertex position -- modified version of vertex calculation 
//                                       w/ rejection of clusters w/ bad chi2 fits 
//                                       & creation of a vertexdiff.txt info file 
//                                               gives events w/ bad chi2

#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TGraph2D.h>
#include <TH2Poly.h>
#include <TClonesArray.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <TString.h>
#include <TROOT.h>
#include <TRint.h>
#include <TTree.h>
#include <TFile.h>
#include <TCutG.h>
#include <TChain.h>
#include <TH1F.h>
#include <TH2.h>
#include <TH2F.h>
#include <TF1.h>
#include <TH3.h>
#include <TH3F.h>
#include <TF3.h>
#include <TMarker3DBox.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TRandom2.h>
#include <TStyle.h>
#include <TF2.h>
#include <TH1.h>
#include <TMinuit.h>
#include <TVirtualFitter.h>
#include <TFitter.h>
#include <TPolyLine3D.h>
#include <TH2Poly.h>
#include "/Users/acorsi/programs/root/include/Math/Vector3D.h"
#include <TLine.h>
#include <TArc.h>
#include <sstream>
#include <TPolyMarker.h>

// Home-made libraries
#include "lib/TMinosClust.h"
#include "lib/TMinosResult.h"


double PI = TMath::Pi();
using namespace ROOT::Math;
using namespace std;


    TClonesArray *data;
    TMinosClust *minosdata;
    TClonesArray data_result;
    TMinosResult *minosdata_result;

// char* filename = "/data/rootfiles/after/R2013_03_08-16_47_07-000.root";
void fflushstdin (){
    int c;
    while ((c = fgetc (stdin)) != EOF && c != '\n');
}

double FitFunction(double *x, double *p)
{
    double val=p[0]+p[1]*x[0];
    return(val);
}
void FindStart(double pStart[4], double chi[2],  int fitStatus[2],TGraph *grxz, TGraph *gryz)
{
    double par1D[2],ndof[2];
    TF1 *myfit1 = new TF1("myfit1",FitFunction, -100,500,2);
    myfit1->SetParameter(0,0);
    myfit1->SetParameter(1,10);
    fitStatus[0] =0;
    grxz->Fit(myfit1,"RQM");
    chi[0]=myfit1->GetChisquare();
    par1D[0]=myfit1->GetParameter(0);
    par1D[1]=myfit1->GetParameter(1);
    pStart[0]=par1D[0];
    pStart[1]=par1D[1];
    fitStatus[1] =0; 
    gryz->Fit(myfit1,"RQM");
    chi[1]=myfit1->GetChisquare();
    par1D[0]=myfit1->GetParameter(0);
    par1D[1]=myfit1->GetParameter(1);
    pStart[2]=par1D[0];
    pStart[3]=par1D[1];
}
// Calculation of the distance line-point

// Calculation of the distance line-point
double distance2(double x,double y,double z, double *p)
{
    // distance line point is D= | (xp-x0) cross  ux |
    // where ux is direction of line and x0 is a point in the line (like t = 0)
    XYZVector xp(x,y,z); //point of the track
    XYZVector x0(p[0], p[2], 0. );
    XYZVector x1(p[0] + p[1], p[2] + p[3], 1. ); //line 
    XYZVector u = (x1-x0).Unit(); 
    double d2 = ((xp-x0).Cross(u)) .Mag2();
    return d2;
}
void SumDistance(int &, double *, double & sum, double * par,  int)
{
    int nused=0;
    double qtot=0;
    sum = 0;
    //cout<<"sum "<<sum<<" over "<<npoints<<endl;
    //double factor;
    //cout<<"*************after fit "<<endl;
    for(int i=0; i<data_result.GetEntriesFast(); i++) 
    {
    	minosdata_result = (TMinosResult*)data_result.At(i);
    	if(minosdata_result->n_Cluster==1)
    	{
        	float x=minosdata_result->x_mm;
        	float y=minosdata_result->y_mm;
        	float z=minosdata_result->z_mm;
        	float q=minosdata_result->Chargemax;
        	//if(nused<2)
        	//cout<<minosdata_result->n_Cluster<<" "<<x<<" "<<y<<" "<<z<<" "<<q<<endl;
        	double d = distance2(x, y, z, par);
      		sum += d*q;
      		nused++;
		qtot+=q;
    	}
    }
    //sum/=nused;
    sum/=qtot;
}



//void Hough_filter(TCanvas *c1, vector<double> *x,vector<double> *y,vector<double> *z,vector<double> *q, vector<double> *x_out,vector<double> *y_out,vector<double> *z_out,vector<double> *q_out)
void Hough_filter(vector<double> *x,vector<double> *y,vector<double> *z,vector<double> *q, vector<double> *x_out,vector<double> *y_out,vector<double> *z_out,vector<double> *q_out)
{
	int nt_xy=180;
	int nt_xz=180;
	int nt_yz=180;
	int nr_xy=52;
	int nr_xz=300;
	int nr_yz=300;
	double bint_xy=2.;
	double bint_xz=2.;
	double bint_yz=2.;
	double binr_xy=3.;
	double binr_xz=3.;
	double binr_yz=3.;
	int nt,nr;
	
	double rho_xy,rho_xz,rho_yz;
	double theta_xy,theta_xz,theta_yz;
	double tmin, tmax;
	
	TH2F *hp_xy = new TH2F("hp_xy","hp_xy",nt_xy,0,180,nr_xy,-52,52);
	TH2F *hp_xz = new TH2F("hp_xz","hp_xz",nt_xz,0,180,nr_xz,-300,300);
	TH2F *hp_yz = new TH2F("hp_yz","hp_yz",nt_yz,0,180,nr_yz,-300,300);

	TH2F *hc_xy = new TH2F("hc_xy","Track in xy plane",100,-150,150,100,-150,150);
	TH2F *hc_xz = new TH2F("hc_xz","Track in xz plane",250,-250,250,100,-150,150);
	TH2F *hc_yz = new TH2F("hc_yz","Track in yz plane",250,-250,250,100,-150,150);

	TH2F *hcnew_xy = new TH2F("hcnew_xy","Track in xy plane AFTER Hough transform",100,-150,150,100,-150,150);
	TH2F *hcnew_xz = new TH2F("hcnew_xz","Track in xz plane AFTER Hough transform",250,-250,250,100,-150,150);
	TH2F *hcnew_yz = new TH2F("hcnew_yz","Track in yz plane AFTER Hough transform",250,-250,250,100,-150,150);

	int npeaks_xy, npeaks_xz, npeaks_yz;
        vector<double> thetapeaks_xy, rpeaks_xy, thetapeaks_xz, rpeaks_xz, thetapeaks_yz, rpeaks_yz;
	double max_xy, max_xz, max_yz;
	double rmean_xy=0, thetamean_xy=0, rmean_xz=0, thetamean_xz=0, rmean_yz=0, thetamean_yz=0;
//	TF1* line_xy = new TF1("line_xy","[0] + [1]*x",-85,85);
//	TF1* line_xz = new TF1("line_xz","[0] + [1]*x",-50,450);
//	TF1* line_yz = new TF1("line_yz","[0] + [1]*x",-50,450);	

	double r0_xy=0., r0_xz=0., r0_yz=0., rmin_xy=0., rmin_xz=0., rmin_yz=0., rmax_xy=0., rmax_xz=0., rmax_yz=0.;
	double rinf=0., rsup=0.;

	nt=nt_xy;
	nr=nr_xy;
	if(nt<nt_xz)nt=nt_xz;
	if(nr<nr_xz)nr=nr_xz;
	if(nt<nt_yz)nt=nt_yz;
	if(nr<nr_yz)nr=nr_yz;
	
	cout<<"In filter !!!"<<endl;
	cout<<"size: "<<x->size()<<endl;
	
	for(int i=0;i<x->size();i++)
	{
		//Loop of indices and fill Histograms
		for(int j=0;j<nt;j++)
		{
			//Fill coordinate space histograms for plots
			hc_xy->Fill(x->at(i),y->at(i),q->at(i)); 
			hc_xz->Fill(z->at(i),x->at(i),q->at(i)); 
			hc_yz->Fill(z->at(i),y->at(i),q->at(i)); 

			//xy
			theta_xy = j*180./nt_xy;
			rho_xy = x->at(i)*TMath::Cos(theta_xy*PI/180.)+y->at(i)*TMath::Sin(theta_xy*PI/180.);
			if(abs(theta_xy)<180.)// && abs(rho_xy)<52)
			{
				//if(i%40==0) cout<<"i="<<i<<" xy "<<rho_xy<<" "<<theta_xy<<endl;
				hp_xy->Fill(theta_xy,rho_xy);
			}

			//xz
			theta_xz = j*180./nt_xz;
			rho_xz = z->at(i)*TMath::Cos(theta_xz*PI/180.)+x->at(i)*TMath::Sin(theta_xz*PI/180.);
			if(abs(theta_xz)<180.)// && abs(rho_xz)<300)
			{
				//if(i%40==0) cout<<"i="<<i<<" xz "<<rho_xz<<" "<<theta_xz<<endl;
				hp_xz->Fill(theta_xz,rho_xz);
			}

			//yz
			theta_yz = j*180./nt_yz;
			rho_yz = z->at(i)*TMath::Cos(theta_yz*PI/180.)+y->at(i)*TMath::Sin(theta_yz*PI/180.);
			if(abs(theta_yz)<180.)// && abs(rho_yz)<300)
			{
				//if(i%40==0) cout<<"i="<<i<<" yz "<<rho_yz<<" "<<theta_yz<<endl;
				hp_yz->Fill(theta_yz,rho_yz);
			}
		}
	}


	max_xy = hp_xy->GetMaximum();
	max_xz = hp_xz->GetMaximum();
	max_yz = hp_yz->GetMaximum();
	
	for(int ii=0; ii<nt; ii++)
	{
		for(int jj=0; jj<nr; jj++)
		{
			if(hp_xy->GetBinContent(ii+1, jj+1) == max_xy)
			{
				thetapeaks_xy.push_back((ii+0.5)*180./nt);
				rpeaks_xy.push_back((jj+0.5)*2. - 52.);
				rmean_xy += rpeaks_xy.back();
				thetamean_xy += thetapeaks_xy.back();
//				cout << "xy: " << thetapeaks_xy.back() << " , " << rpeaks_xy.back() << endl;
			}
			if(hp_xz->GetBinContent(ii+1, jj+1) == max_xz)
			{
				thetapeaks_xz.push_back((ii+0.5)*180./nt);
				rpeaks_xz.push_back((jj+0.5)*2. - 150.);
				rmean_xz += rpeaks_xz.back();
				thetamean_xz += thetapeaks_xz.back();
				//cout << "xz: " << thetapeaks_xz.back() << " , " << rpeaks_xz.back() << endl;
			}
			if(hp_yz->GetBinContent(ii+1, jj+1) == max_yz)
			{
				thetapeaks_yz.push_back((ii+0.5)*180./nt);
				rpeaks_yz.push_back((jj+0.5)*2. - 150.);
				rmean_yz += rpeaks_yz.back();
				thetamean_yz += thetapeaks_yz.back();
				//cout << "yz: " << thetapeaks_yz.back() << " , " << rpeaks_yz.back() << endl;
			}
		}
	}

//	cout << "Number of max found :::     IN xy = " << rpeaks_xy.size() << " ,     IN xz = " << rpeaks_xz.size() << " ,     IN yz = " << rpeaks_yz.size() << endl;

	// xy PEAK
	rmean_xy = rmean_xy / rpeaks_xy.size();
	thetamean_xy = thetamean_xy / thetapeaks_xy.size();
//	line_xy->SetParameter(0,rmean_xy/(TMath::Sin(thetamean_xy*PI/180)));
//	line_xy->SetParameter(1,( -(TMath::Cos(thetamean_xy*PI/180))/(TMath::Sin(thetamean_xy*PI/180)) ));
//	hc_xy->GetListOfFunctions()->Add(line_xy);

	// xz PEAK
	rmean_xz = rmean_xz / rpeaks_xz.size();
	thetamean_xz = thetamean_xz / thetapeaks_xz.size();
//	line_xz->SetParameter(0,rmean_xz/(TMath::Sin(thetamean_xz*PI/180)));
//	line_xz->SetParameter(1,( -(TMath::Cos(thetamean_xz*PI/180))/(TMath::Sin(thetamean_xz*PI/180)) ));
//	hc_xz->GetListOfFunctions()->Add(line_xz);

	// yz PEAK
	rmean_yz = rmean_yz / rpeaks_yz.size();
	thetamean_yz = thetamean_yz / thetapeaks_yz.size();
//	line_yz->SetParameter(0,rmean_yz/(TMath::Sin(thetamean_yz*PI/180)));
//	line_yz->SetParameter(1,( -(TMath::Cos(thetamean_yz*PI/180))/(TMath::Sin(thetamean_yz*PI/180)) ));
//	hc_yz->GetListOfFunctions()->Add(line_yz);

/*
	line_xy->SetLineWidth(1);
	line_xz->SetLineWidth(1);
	line_yz->SetLineWidth(1);
*/
	//Selection of x,y,z points COMMON to the 3 maxmean+/-1 found in Hough spaces for xy, xz and yz spaces
	for(unsigned int i=0;i<x->size();i++)
	{
		r0_xy = x->at(i)*TMath::Cos(thetamean_xy*PI/180.)+y->at(i)*TMath::Sin(thetamean_xy*PI/180.);
		tmin = thetamean_xy-bint_xy;
		tmax = thetamean_xy+bint_xy;
		if((tmin)<0) tmin = tmin + 180.;
		if((tmax)>180) tmax = tmax - 180.;
		rmin_xy = x->at(i)*TMath::Cos(tmin*PI/180.)+y->at(i)*TMath::Sin(tmin*PI/180.);
		rmax_xy = x->at(i)*TMath::Cos(tmax*PI/180.)+y->at(i)*TMath::Sin(tmax*PI/180.);

		rinf = min( rmean_xy - binr_xy, rmean_xy + binr_xy);
		rsup = max( rmean_xy - binr_xy, rmean_xy + binr_xy);
		cout<<"before 1st if"<<endl;
		if((r0_xy>=rinf || rmin_xy>=rinf || rmax_xy>=rinf) && (r0_xy<=rsup || rmin_xy<=rsup || rmax_xy<=rsup))
		{
			r0_xz = z->at(i)*TMath::Cos(thetamean_xz*PI/180.)+x->at(i)*TMath::Sin(thetamean_xz*PI/180.);
			tmin = thetamean_xz-bint_xz;
			tmax = thetamean_xz+bint_xz;
			if((tmin)<0) tmin = tmin + 180.;
			if((tmax)>180) tmax = tmax - 180.;
			rmin_xz = z->at(i)*TMath::Cos(tmin*PI/180.)+x->at(i)*TMath::Sin(tmin*PI/180.);
			rmax_xz = z->at(i)*TMath::Cos(tmax*PI/180.)+x->at(i)*TMath::Sin(tmax*PI/180.);

			rinf = min( rmean_xz - binr_xz, rmean_xz + binr_xz);
			rsup = max( rmean_xz - binr_xz, rmean_xz + binr_xz);
			cout<<"before 2nd if"<<endl;
			if((r0_xz>=rinf || rmin_xz>=rinf || rmax_xz>=rinf) && (r0_xz<=rsup || rmin_xz<=rsup || rmax_xz<=rsup))
			{
				r0_yz = z->at(i)*TMath::Cos(thetamean_yz*PI/180.)+y->at(i)*TMath::Sin(thetamean_yz*PI/180.);
				tmin = thetamean_yz-bint_yz;
				tmax = thetamean_yz+bint_yz;
				if((tmin)<0) tmin = tmin + 180.;
				if((tmax)>180) tmax = tmax - 180.;
				rmin_yz = z->at(i)*TMath::Cos(tmin*PI/180.)+y->at(i)*TMath::Sin(tmin*PI/180.);
				rmax_yz = z->at(i)*TMath::Cos(tmax*PI/180.)+y->at(i)*TMath::Sin(tmax*PI/180.);

				rinf = min( rmean_yz - binr_yz, rmean_yz + binr_yz);
				rsup = max( rmean_yz - binr_yz, rmean_yz + binr_yz);
				cout<<"before 3rd if"<<endl;
				if((r0_yz>=rinf || rmin_yz>=rinf || rmax_yz>=rinf) && (r0_yz<=rsup || rmin_yz<=rsup || rmax_yz<=rsup))
				{
					cout << "Taken points= " << x->at(i) << " , " << y->at(i) << " , " << z->at(i) << endl;
					hcnew_xy->Fill(x->at(i),y->at(i),q->at(i));
					hcnew_xz->Fill(z->at(i),x->at(i),q->at(i)); 
					hcnew_yz->Fill(z->at(i),y->at(i),q->at(i)); 
					x_out->push_back(x->at(i));
					y_out->push_back(y->at(i));
					z_out->push_back(z->at(i));
					q_out->push_back(q->at(i));
				}

			}
		}



	}

/*
	c1->Divide(3,3);
	// Coordinate space
	c1->cd(1);
	hc_xy->Draw("colz");
	c1->cd(2);
	hc_xz->Draw("colz");
	c1->cd(3);
	hc_yz->Draw("colz");

	// Hough space
	c1->cd(4);
	hp_xy->Draw("colz");
	c1->cd(5);
	hp_xz->Draw("colz");
	c1->cd(6);
	hp_yz->Draw("colz"); 
  		
	// Coordinate space : New plots
	c1->cd(7);
	hcnew_xy->Draw("colz");
	c1->cd(8);
	hcnew_xz->Draw("colz");
	c1->cd(9);
	hcnew_yz->Draw("colz");

	c1->Update();*/

	delete hp_xy;
	delete hp_xz;
	delete hp_yz;

}




// Calculation of the minimal distance between 2 lines in 3D space & calculation of mid-point=>vertex of interaction
void vertex(double *p, double *pp, double &xv,double &yv,double &zv)
{ 
    double a1 = p[0];
    double a2 = p[2];
    double b1 = p[1];
    double b2 = p[3];
    double ap1 = pp[0];
    double ap2 = pp[2];
    double bp1 = pp[1];
    double bp2 = pp[3];
    
    double alpha, beta, A, B, C;
    
    alpha = (bp1*(a1-ap1)+bp2*(a2-ap2))/(bp1*bp1 + bp2*bp2 + 1);
    beta = (bp1*b1+bp2*b2+1)/(bp1*bp1 + bp2*bp2 + 1);
    
    A = beta*(bp1*bp1 + bp2*bp2 + 1) - (bp1*b1 + bp2*b2 + 1);
    B = (b1*b1 + b2*b2 + 1) - beta*(bp1*b1+bp2*b2+1);
    C = beta*(bp1*(ap1-a1) + bp2*(ap2-a2)) - (b1*(ap1-a1) + b2*(ap2-a2));
    
    double sol1, solf1;
    double x,y,z,xp,yp,zp;


    sol1 = -(A*alpha + C)/(A*beta + B);
    solf1 = alpha + beta* sol1;

    x = a1 + b1*sol1;
    y = a2 + b2*sol1;
    z = sol1;
    xp = ap1 + bp1*solf1;
    yp = ap2 + bp2*solf1;
    zp = solf1;
    
    xv = (x+xp)/2.;
    yv = (y+yp)/2.;
    zv = (z+zp)/2.;
    
    //cout << "Vertex 1st :" << x << "," << y << "," << z << endl;
    //cout << "Vertex 2nd :" << xp << "," << yp << "," << zp << endl;
    //cout << "Vertex middle :" << xv << "," << yv << "," << zv << endl;
    
    //cout << "min dist " << sqrt(pow((x-xp),2) + pow((y-yp),2) + pow((z-zp),2)) << endl;
    
}

void plotvertexHoughNew(char* filename)
{
    
    //load file    
    TFile *f=new TFile(filename);
    char path[80];
    data = new TClonesArray("TMinosClust");
    TTree *tree = (TTree*)f->Get("tree_out");
    tree->GetBranch("data_out");
    tree->SetBranchAddress("data_out",&data);
    int trackNbr;
    tree->GetBranch("trackNbr");
    tree->SetBranchAddress("trackNbr",&trackNbr);
    int evtOrig;
    tree->GetBranch("evtOrig");
    tree->SetBranchAddress("evtOrig",&evtOrig);
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

    //freopen("vertexinfo.txt","w",stdout);

    // ROOT output file (locate number of run, time & date)
    string pchroot= filename;
    string pch2 = ".root";
    pchroot.replace(pchroot.find(pch2),pch2.length(),"_result.root");
    const char *OutFile = pchroot.c_str();

    cout<<"***"<<OutFile<<endl;
    
    TFile *fout = new TFile(OutFile,"RECREATE");
    TTree *tree_result = new TTree("tree_result","Result tree");
    data_result.SetClass("TMinosResult");
    tree_result->Branch("data_result",&data_result);
    tree_result->Branch("evtOrig",&evtOrig,"evtOrig/I");
    int trackNbr_FINAL;
    tree_result->Branch("trackNbr_FINAL",&trackNbr_FINAL,"trackNbr_FINAL/I");

    double z_vertex=0., x_vertex=0., y_vertex=0., r_vertex=0.;

    int all_1pclusters=0, all_2pclusters=0, all_1pfiltered=0, all_2pfiltered=0, all_1pvertex=0, all_2pvertex=0;
    double zmax=0.;
    int array_final=0;

    //vector<TCanvas*> Hough_canvas;

    //make a vector with the event number corresponding to good events ie above threshold and with z!=-1 (pads not fitted)   
    vector<int> goodev;//if(event==NULL) cout<<"error"<<endl;
    int goodevcnt1, goodevcnt2;
    for(int ev=0; ev<tree->GetEntries(); ev++)
    {
	tree->GetEntry(ev);
	//if(trackNbr==1 || trackNbr==2 || trackNbr==3 || trackNbr==4 /*|| trackNbr==5*/) goodev.push_back(ev);
	if(trackNbr>0) goodev.push_back(ev);
    }
    if(goodev.size()==0) cout<<"no good events"<<endl;

    goodevcnt1=0;
    goodevcnt2=0;
    double nbrct=0;
    int npoint_temp=0, cluster_temp=0;
    int cluster1=0, cluster2=0;
    int ringsum=0;
    int ringtouch[24]={0};
    double pStart[4];

    //start loop on goodev
    TGraph * gryz_1;
    TGraph * grxz_1;
    TGraph * gryz_2;
    TGraph * grxz_2;
    TF1 *grxz_fit_1;
    TF1 *gryz_fit_1;
    TF1 *grxz_fit_2;
    TF1 *gryz_fit_2;
    vector<double> xin, yin, zin, qin, xout, yout, zout, qout; 
    TGraph* grxz_clone;
    TGraph* gryz_clone;
    vector<Double_t> parFit1,parFit2,parFit3,parFit4;
	double x_vertex_average;
	double y_vertex_average;
	double z_vertex_average;
	vector<double> x_vertex_vect;
	vector<double> y_vertex_vect;
	vector<double> z_vertex_vect;
	vector<double> r_vertex_vect;

    tree_result->Branch("x_vertex_vect",&x_vertex_vect);
    tree_result->Branch("y_vertex_vect",&y_vertex_vect);
    tree_result->Branch("z_vertex_vect",&z_vertex_vect);
    tree_result->Branch("r_vertex_vect",&r_vertex_vect);
    tree_result->Branch("x_vertex_average",&x_vertex_average,"x_vertex_average/D");
    tree_result->Branch("y_vertex_average",&y_vertex_average,"y_vertex_average/D");
    tree_result->Branch("z_vertex_average",&z_vertex_average,"z_vertex_average/D");
    TGraph* gryz_tmp; //= new TGraph();
   	TGraph* grxz_tmp; // = new TGraph();
   	vector<TGraph*> grxz;
   	vector<TGraph*> gryz;
    int npoint=0;
    //for(int ii=0; ii<10;ii++)
    for(int ii=0; ii<goodev.size();ii++)
    {
    	int npoint1,npoint2;
    	int entry=goodev.at(ii);
    	tree->GetEntry(entry);
	cerr << "EVENT " << ii << " * (original evt " << evtOrig << ")" << endl;
//	if(ii%100==0)cout << "EVENT " << ii << " * (original evt " << evtOrig << ")" << endl;

	data_result.Clear();
	trackNbr_FINAL=0;
	array_final=0;
	ringsum=0;
    	npoint1=0;npoint2=0;
	cluster1=0, cluster2=0;
	cluster_temp=0;
	zmax=0.;
	x_vertex=0., y_vertex=0., z_vertex=0.;

    grxz.clear();
    gryz.clear();

   	float rmin1=100; float rmin2=100;
   	float rmax1=00; float rmax2=00;
   	float deltar1=0,deltar2=0,radius=0;
	cluster1 = 0;



    TMinuit *min ;
    Double_t pStart[4];
    double parFit_temp[4], err_temp[4];
    Double_t chi[2];
    Int_t fitStatus[2];
    Double_t arglist[10];
    Int_t iflag;
    int nvpar,nparx;
    double amin,edm, errdef;
	cout<<"entr"<<data->GetEntriesFast()<<endl;
   	for(int i=0;i<(data->GetEntriesFast());i++) 
    {
		minosdata = (TMinosClust*)data->At(i);
		cout<<"nev "<<i<<" "<<cluster_temp<<" "<<minosdata->n_Cluster<<endl;
		if( xin.size()>0 && ((cluster_temp!=int(minosdata->n_Cluster) && i!=0) || i==(data->GetEntriesFast() - 1)))
		{
			cerr << i << " *** cluster_temp:"<<cluster_temp << " "<<minosdata->n_Cluster << " "<<data->GetEntriesFast()-1<<endl;
   			//Hough_canvas.push_back(new TCanvas(Form("Event%d_cluster%d", evtOrig, cluster_temp), Form("Event%d_cluster%d", evtOrig, cluster_temp)));
			//Hough_filter(Hough_canvas.back(), &xin, &yin, &zin, &qin, &xout, &yout, &zout, &qout);
			Hough_filter(&xin, &yin, &zin, &qin, &xout, &yout, &zout, &qout);
			gryz_tmp = new TGraph();
   			grxz_tmp = new TGraph();
   			//Hough_canvas.back()->SaveAs(Form("Event%d_cluster%d.pdf", evtOrig, cluster_temp));

			//Hough_canvas.back()->Write();
			cout<<"size "<<yout.size()<<endl;npoint=0;
			for(unsigned int ij=0; ij<xout.size();ij++)
			{
				if(zout[ij]>zmax) zmax = zout[ij];
				ringtouch[int((sqrt(xout[ij]*xout[ij]+yout[ij]*yout[ij])-52)/3)]++;
			}
			cerr << "!!!evtOrig " << evtOrig << ", cluster " << cluster_temp << " :: Npads=" << xout.size() << ", zmax=" << zmax << endl; 
			for(int ko=0; ko<24; ko++)
			{
				if(ringtouch[ko]>0) ringsum++;
				cerr << "     ------   ringtouch(" << ko << ")=" << ringtouch[ko] << endl;
			}
			cerr << "Ringsum=" << ringsum;
			if(zmax>150||zmax<-150) ringsum=25;
			cerr << ", Newringsum=" << ringsum << endl;
			if(xout.size()>10 && ringsum<=24)
			{
					trackNbr_FINAL++;
					for(unsigned int ij=0; ij<xout.size(); ij++)
					{	
						grxz_tmp->SetPoint(npoint,zout[ij],xout[ij]);
            			gryz_tmp->SetPoint(npoint,zout[ij],yout[ij]);
						minosdata_result = (TMinosResult*)data_result.ConstructedAt(array_final);
                		minosdata_result->Set(xout[ij], yout[ij], zout[ij], qout[ij], 1, xout.size(), zmax);
						array_final++;
						npoint++;
						//cout<<yout[ij]<<" "<<xout[ij]<<" "<<zout[ij]<<endl;
					}
				
				grxz.push_back(grxz_tmp);
				gryz.push_back(gryz_tmp);

			}
			
			//cout<<"Evt "<<entry<<" data point "<<i<<" track Nr = "<<trackNbr_FINAL<<" xout.size ="<<xout.size()<<" data cluster "<<cluster_temp<<endl;

			xin.clear();
			yin.clear();
			zin.clear();
			qin.clear();
			xout.clear();
			yout.clear();
			zout.clear();
			qout.clear();
			npoint_temp=0;
			ringsum=0;
			zmax=0.;
			for(int ko=0; ko<18; ko++) ringtouch[ko] = 0;

		}
		
		cluster_temp = minosdata->n_Cluster;
		
		if(! (minosdata->n_Pads>=10 && minosdata->RingBool==1 && minosdata->z_mm>-10000 && minosdata->z_mm<=320) ) {cout<<"continue"<<minosdata->n_Pads<<" "<<minosdata->RingBool<<endl;continue;}
		else
		{
				//cout<<minosdata->x_mm<<" "<<minosdata->z_mm<<" "<<endl;
				cout<<"trackNbr "<<trackNbr<<endl;
				xin.push_back(minosdata->x_mm);
				yin.push_back(minosdata->y_mm);
				zin.push_back(minosdata->z_mm);
				qin.push_back(minosdata->Chargemax);
				npoint_temp++;
		}

    }// end for GetEntriesFast

//	cerr << "EvtOrig=" << evtOrig << " ; trackNbrFINAL=" << trackNbr_FINAL << endl;
    int NclusterFit=0; 
    cout<<"nbfinal "<<trackNbr_FINAL<<" "<<grxz.size()<<endl;
	if(trackNbr_FINAL>0 )
	{

                    //////////Minimization in 3D to reconstruct track lines
                    
                    
                    for(int itr=0; itr<trackNbr_FINAL; itr++) 
                    {
			            pStart[0]=0; pStart[2]=0; pStart[1]=1; pStart[3]=3;
	                        //cout<<"start Minuit "<<endl;
	                    min = new TMinuit(4);
	                    min->SetPrintLevel(-1);
	                    arglist[0] = 3;

                        cout<<"debug1"<<endl;
                        FindStart(pStart,chi,fitStatus, grxz.at(itr),gryz.at(itr));
                        grxz_clone=(TGraph*)grxz.at(itr)->Clone();
						grxz_clone->SetName(Form("grxz_%d_%d",ii,itr));grxz_clone->Draw("A*");grxz_clone->Write();
						gryz_clone=(TGraph*)gryz.at(itr)->Clone();
						gryz_clone->SetName(Form("gryz_%d_%d",ii,itr));gryz_clone->Draw("A*");gryz_clone->Write();
						

                        NclusterFit = itr+1;
                        min->SetFCN(SumDistance);
                        // Set starting values and step sizes for parameters
                        min->mnparm(0,"x0",pStart[0],0.1,-500,500,iflag);
                        min->mnparm(1,"Ax",pStart[1],0.1,0,0,iflag);
                        min->mnparm(2,"y0",pStart[2],0.1,-500,500,iflag);
                        min->mnparm(3,"Ay",pStart[3],0.1,0,0,iflag);
                        /*min->mnparm(0,"x0",pStart[0],0.1,0,0,iflag);
                        min->mnparm(1,"Ax",pStart[1],0.1,-10,10,iflag);
                        min->mnparm(2,"y0",pStart[2],0.1,0,0,iflag);
                        min->mnparm(3,"Ay",pStart[3],0.1,-10,10,iflag);*/
                        arglist[0] = 100; // number of function calls
                        arglist[1] = 0.000001; // tolerance
                        min->mnexcm("MIGRAD",arglist,2,iflag); // minimization with MIGRAD
                        cout<<"debug1"<<endl;
                        min->mnstat(amin,edm,errdef,nvpar,nparx,iflag);  //returns current status of the minimization
                        // get fit parameters
                        for(int i = 0; i <4; i++) {min->GetParameter(i,parFit_temp[i],err_temp[i]);cout<<parFit_temp[i]<<endl;}
                        
                        parFit1.push_back(parFit_temp[0]);
                        parFit2.push_back(parFit_temp[1]);
                        parFit3.push_back(parFit_temp[2]);
                        parFit4.push_back(parFit_temp[3]);

                        //offset SAMURAI-mid target+mid target-target entrance+target entrance+TPC MM plane (z=0 in local coord)

						delete min;
        			}
    	cout<<"nbfinal "<<trackNbr_FINAL<<" "<<grxz.size()<<endl;
        
        Double_t parFitA[4], parFitB[4];
        if(trackNbr_FINAL>1)
        {
        	        for(int i=0; i<trackNbr_FINAL; i++)
			        {
			        	for(int j=i+1; j<trackNbr_FINAL; j++)
			        	{
			        		parFitA[0]=parFit1[i];parFitA[1]=parFit2[i];parFitA[2]=parFit3[i];parFitA[3]=parFit4[i];
			        		parFitB[0]=parFit1[j];parFitB[1]=parFit2[j];parFitB[2]=parFit3[j];parFitB[3]=parFit4[j];
			        		vertex(parFitA, parFitB, x_vertex, y_vertex, z_vertex);cout<<"vertex "<<z_vertex<<endl;
			        		x_vertex_average+=x_vertex;x_vertex_vect.push_back(x_vertex);
			         		y_vertex_average+=y_vertex;y_vertex_vect.push_back(y_vertex);
			       			z_vertex_average+=z_vertex;z_vertex_vect.push_back(z_vertex);
			      		}
			      	}

      		parFit1.clear();parFit2.clear();parFit3.clear();parFit4.clear();
      		x_vertex_average=x_vertex_average/trackNbr_FINAL;
      		y_vertex_average=y_vertex_average/trackNbr_FINAL;
      		z_vertex_average=z_vertex_average/trackNbr_FINAL;
        }
        else 
        {
      		parFit1.clear();parFit2.clear();parFit3.clear();parFit4.clear();
      		x_vertex_average=-999;
      		y_vertex_average=-999;
      		z_vertex_average=-999;


        }


      		grxz.clear();
      		gryz.clear();



      		
	}// end if trackNbr_FINAL==1 || 2

	tree_result->Fill();
	x_vertex_vect.clear(); y_vertex_vect.clear(); z_vertex_vect.clear();
	x_vertex_average=y_vertex_average=z_vertex_average=0;
	
    }
    tree_result->Write();
    float x1, w1, x2, w2;
    int nentries;
    tree->GetEntry(tree->GetEntries()-1); 
    nentries = evtOrig+1;




}




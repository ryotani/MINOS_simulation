
generatexsec()
{ 

  	TH2F*h=new TH2F("h","h",300,0,3,90,0,90);
  	TRandom *rand=new TRandom();
  	float x,y;
  	TFile*f=new TFile("xsec_K.root","recreate");
  	for(int i=0; i<100000; i++)
  	{
  		x=rand->Gaus()*0.6;
  		y=rand->Gaus()*0.6;
  		if(x*y<.05&&(x+y)<0.7&&x*y>.003) h->Fill(x*3,y*90);
      //if(x*y>.005&&x*y<0.1) h->Fill(x*3,y*90);

  	}
  	h->Write();
  	h->Draw("colz");
  	f->Close();
}
#include "ExN03DataoRings.hh"
ClassImp(ExN03DataoRings)

ExN03DataoRings::ExN03DataoRings(){
ReadConfigurationFile("ConfigFileRings.txt");
}

ExN03DataoRings::~ExN03DataoRings(){
}

void ExN03DataoRings::ClearEvent(){
datai.x_InRoh.clear();
datai.y_InRoh.clear();
datai.z_InRoh.clear();
datai.e_InRoh.clear();

datai.x_OutRoh.clear();
datai.y_OutRoh.clear();
datai.z_OutRoh.clear();
datai.e_OutRoh.clear();

datai.x_Kap.clear();
datai.y_Kap.clear();
datai.z_Kap.clear();
datai.e_Kap.clear();

datai.x_tpc.clear();
datai.y_tpc.clear();
datai.z_tpc.clear();
datai.e_tpc.clear();

datai.x_tar.clear();
datai.y_tar.clear();
datai.z_tar.clear();
datai.e_tar.clear();

datai.x_ch.clear();
datai.y_ch.clear();
datai.z_ch.clear();
datai.e_ch.clear();

x_pad.clear();
y_pad.clear();
t_pad.clear();
q_pad.clear();
dt.clear();

datai.x_win.clear();
datai.y_win.clear();
datai.z_win.clear();
datai.e_win.clear();

datai.A.clear();
datai.Z.clear();
datai.Et_tpc_tot=0;
datai.Et_tar.clear();
datai.Et_ch.clear();
datai.Et_tpc.clear();
datai.Et_win.clear();
Et_pad=0;
datai.Et_InnerRohacell.clear();
datai.Et_OuterRohacell.clear();
datai.Et_Kapton.clear();
nb_pads=0;

datai.x0.clear();
datai.y0.clear();
datai.z0.clear();
datai.theta0.clear();
datai.phi0.clear();
datai.energy0.clear();
datai.detection.clear();
}

void ExN03DataoRings::ReadConfigurationFile(string Config)
{
   string LineBuffer;
   string DataBuffer;

   bool cNRings    = false;
   bool cNSegments = false;
   bool cTimeBin   = false;	
   bool csigL      = false;
   bool csigT      = false;
   bool cdriftV    = false;
   bool cIonis     = false;
   bool cThreshold = false;
   bool cGain      = false;
   bool cNoiseRMS  = false;
   bool cTheta      = false;
   bool cFirstRing = false;
   bool cLastRing = false;
   bool cShapingTime = false;

   string Path = "/Users/acorsi/codes/MINOS_simulation/Inputs/"+Config;
   ifstream ConfigFile;
   ConfigFile.open(Path.c_str());

   if (ConfigFile.is_open())
      cout << " Configuration file " << Config << " loading " << endl;
   else {
      cout << " Error, no configuration file" << Config << " found" << endl;
      return;
   }


   while (!ConfigFile.eof()) {
      //Pick-up next line
      getline(ConfigFile, LineBuffer);
      //Search for comment Symbol: %
      if (LineBuffer.compare(0, 1, "%") == 0) {   /*Do  Nothing*/;}

      else if (LineBuffer.compare(0,6, "NRings") == 0 && cNRings == false) {
         cNRings = true ;
	 string test = LineBuffer.substr(6,LineBuffer.length()-6);
         NRings=atoi(test.c_str());
	 cout << "NRings = " << NRings<<endl;
      }
      else if (LineBuffer.compare(0,9, "NSegments") == 0 && cNSegments == false) {
         cNSegments = true ;
	 string test = LineBuffer.substr(9,LineBuffer.length()-9);
         NSegments=atoi(test.c_str());
	 cout << "NSegments = " << NSegments<<endl;
     }
      else if (LineBuffer.compare(0,9, "FirstRing") == 0 && cFirstRing == false) {
         cFirstRing = true ;
	 string test = LineBuffer.substr(9,LineBuffer.length()-9);
         FirstRing=atof(test.c_str());
         cout << "First Ring Distance = "<<FirstRing<<" mm" << endl   ;
     }
      else if (LineBuffer.compare(0,8, "LastRing") == 0 && cLastRing == false) {
         cLastRing = true ;
	 string test = LineBuffer.substr(8,LineBuffer.length()-8);
         LastRing=atof(test.c_str());
         cout << "Last Ring Distance = "<<LastRing<<" mm" << endl   ;
     }
      else if (LineBuffer.compare(0,7, "TimeBin") == 0 && cTimeBin == false) {
         cTimeBin = true ;
	 string test = LineBuffer.substr(7,LineBuffer.length()-7);
         TimeBinSize=atof(test.c_str());
         cout << "Time bin = "<<TimeBinSize<<" ns" << endl   ;
     }
      else if (LineBuffer.compare(0,4, "sigL") == 0 && csigL == false) {
         csigL = true ;
	 string test = LineBuffer.substr(4,LineBuffer.length()-4);
         sigL=atof(test.c_str());
         cout << "Longitudinal dispersion = "<< sigL<<" mm per sqrt(cm)" << endl   ;
     }
      else if (LineBuffer.compare(0,4, "sigT") == 0 && csigT == false) {
         csigT = true ;
	 string test = LineBuffer.substr(4,LineBuffer.length()-4);
         sigT=atof(test.c_str());
         cout << "Transverse dispersion = "<< sigT<<" mm per sqrt(cm)" << endl   ;
     }
      else if (LineBuffer.compare(0,6, "driftV") == 0 && cdriftV == false) {
         cdriftV = true ;
	 string test = LineBuffer.substr(6,LineBuffer.length()-6);
         driftV=atof(test.c_str());
         cout << "Drift velocity = "<< driftV<<" mm per ns" << endl   ;
     }
      else if (LineBuffer.compare(0,5, "Ionis") == 0 && cIonis == false) {
         cIonis = true ;
	 string test = LineBuffer.substr(5,LineBuffer.length()-5);
         Ionis=atof(test.c_str());
         cout << "Ionisation Energy = "<<Ionis<<" eV" << endl   ;
     }
      else if (LineBuffer.compare(0,9, "Threshold") == 0 && cThreshold == false) {
         cThreshold = true ;
	 string test = LineBuffer.substr(9,LineBuffer.length()-9);
         Threshold=atof(test.c_str());
         cout << "Detection Threshold = "<<Threshold<<" times noise rms" << endl   ;
     }
      else if (LineBuffer.compare(0,4, "Gain") == 0 && cGain == false) {
         cGain = true ;
	 string test = LineBuffer.substr(4,LineBuffer.length()-4);
         Gain=atof(test.c_str());
         cout << "Detection Gain = "<<Gain << endl   ;
     }
      else if (LineBuffer.compare(0,5, "Theta") == 0 && cTheta == false) {
         cTheta = true ;
	 string test = LineBuffer.substr(5,LineBuffer.length()-5);
         Theta=atof(test.c_str());
         cout << "Parameter theta for Polya = "<<Theta << endl   ;
     }
      else if (LineBuffer.compare(0,8, "NoiseRMS") == 0 && cNoiseRMS == false) {
         cNoiseRMS = true ;
	 string test = LineBuffer.substr(8,LineBuffer.length()-8);
         NoiseRMS=atof(test.c_str());
         cout << "Noise background = "<<NoiseRMS<<" e- RMS" << endl   ;
     }
      else if (LineBuffer.compare(0,11, "ShapingTime") == 0 && cShapingTime == false) {
         cShapingTime = true ;
	 string test = LineBuffer.substr(11,LineBuffer.length()-11);
         ShapingTime=atoi(test.c_str());
         cout << "ShapingTime = "<<ShapingTime<<" ns" << endl   ;
     }
   }
   ConfigFile.close();
   return   ;
}



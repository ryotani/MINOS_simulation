#include "ExN03BeamIn.hh"
ClassImp(ExN03BeamIn)

ExN03BeamIn::ExN03BeamIn(){
ReadConfigurationFile("ConfigFileBeam.txt");
}

ExN03BeamIn::~ExN03BeamIn(){
}

void ExN03BeamIn::ReadConfigurationFile(string Config)
{
   string LineBuffer;
   string DataBuffer;

   bool cMeanX       = false;
   bool cMeanY       = false;
   bool cMeanMomentumX       = false;
   bool cMeanMomentumY       = false;
   bool cSigmaX       = false;
   bool cSigmaY       = false;
   bool cSigmaMomentumX       = false;
   bool cSigmaMomentumY       = false;   
   bool cMeanEnergy       = false;
   bool cSigmaEnergy      = false;
   bool cMomentumZ0       = false;   
   bool cZ0       = false;
   bool cZ       = false;  
   bool cA       = false;   
   
   string Path = "/Users/acorsi/codes/MINOS_simulation/Inputs/"+Config;
   ifstream ConfigFile;
   ConfigFile.open(Path.c_str());

   if (ConfigFile.is_open())
      cout << " Configuration file for beam " << Config << " loading " << endl;
   else {
      cout << " Error, no configuration file for beam" << Config << " found" << endl;
      return;
   }


   while (!ConfigFile.eof()) {
      //Pick-up next line
      getline(ConfigFile, LineBuffer);
      //Search for comment Symbol: %
      if (LineBuffer.compare(0, 1, "%") == 0) {   /*Do  Nothing*/;}

      else if (LineBuffer.compare(0,5, "MeanX") == 0 && cMeanX == false) {
         cMeanX = true ;
	 string test = LineBuffer.substr(5,LineBuffer.length()-5);
         MeanX=atof(test.c_str());
	 cout << "MeanX = " << MeanX<<" mm"<<endl;
      }
      else if (LineBuffer.compare(0,5, "MeanY") == 0 && cMeanY == false) {
         cMeanY = true ;
	 string test = LineBuffer.substr(5,LineBuffer.length()-5);
         MeanY=atof(test.c_str());
	 cout << "MeanY = " << MeanY<<" mm"<<endl;
      }
      else if (LineBuffer.compare(0,13, "MeanMomentumX") == 0 && cMeanMomentumX == false) {
         cMeanMomentumX = true ;
	 string test = LineBuffer.substr(13,LineBuffer.length()-13);
         MeanMomentumX=atof(test.c_str());
	 cout << "MeanMomentumX = " << MeanMomentumX<<endl;
      }
      else if (LineBuffer.compare(0,13, "MeanMomentumY") == 0 && cMeanMomentumY == false) {
         cMeanMomentumY = true ;
	 string test = LineBuffer.substr(13,LineBuffer.length()-13);
         MeanMomentumY=atof(test.c_str());
	 cout << "MeanMomentumY = " << MeanMomentumY<<endl;
      }
      else if (LineBuffer.compare(0,6, "SigmaX") == 0 && cSigmaX == false) {
         cSigmaX = true ;
	 string test = LineBuffer.substr(6,LineBuffer.length()-6);
         SigmaX=atof(test.c_str());
	 cout << "SigmaX = " << SigmaX<<" mm"<<endl;
      }
      else if (LineBuffer.compare(0,6, "SigmaY") == 0 && cSigmaY == false) {
         cSigmaY = true ;
	 string test = LineBuffer.substr(6,LineBuffer.length()-6);
         SigmaY=atof(test.c_str());
	 cout << "SigmaY = " << SigmaY<<" mm"<<endl;
      }
      else if (LineBuffer.compare(0,14, "SigmaMomentumX") == 0 && cSigmaMomentumX == false) {
         cSigmaMomentumX = true ;
	 string test = LineBuffer.substr(14,LineBuffer.length()-14);
         SigmaMomentumX=atof(test.c_str());
	 cout << "SigmaMomentumX = " << SigmaMomentumX<<endl;
      }
      else if (LineBuffer.compare(0,14, "SigmaMomentumY") == 0 && cSigmaMomentumY == false) {
         cSigmaMomentumY = true ;
	 string test = LineBuffer.substr(14,LineBuffer.length()-14);
         SigmaMomentumY=atof(test.c_str());
	 cout << "SigmaMomentumY = " << SigmaMomentumY<<endl;
      }
      else if (LineBuffer.compare(0,10, "MeanEnergy") == 0 && cMeanEnergy == false) {
         cMeanEnergy = true ;
	 string test = LineBuffer.substr(10,LineBuffer.length()-10);
         MeanEnergy=atof(test.c_str());
	 cout << "MeanEnergy = " << MeanEnergy<<" MeV"<<endl;
      }
      else if (LineBuffer.compare(0,11, "SigmaEnergy") == 0 && cSigmaEnergy == false) {
         cSigmaEnergy = true ;
	 string test = LineBuffer.substr(11,LineBuffer.length()-11);
         SigmaEnergy=atof(test.c_str());
	 cout << "SigmaEnergy = " << SigmaEnergy<<" MeV"<<endl;
      }
      else if (LineBuffer.compare(0,2, "Z0") == 0 && cZ0 == false) {
         cZ0 = true ;
	 string test = LineBuffer.substr(2,LineBuffer.length()-2);
         Z0=atof(test.c_str());
	 cout << "Z0 = " << Z0<<" mm"<<endl;
      }
      else if (LineBuffer.compare(0,10, "MomentumZ0") == 0 && cMomentumZ0 == false) {
         cMomentumZ0 = true ;
	 string test = LineBuffer.substr(10,LineBuffer.length()-10);
         MomentumZ0=atof(test.c_str());
	 cout << "MomentumZ0 = " << MomentumZ0<<endl;
      }
      else if (LineBuffer.compare(0,1, "Z") == 0 && cZ == false) {
         cZ = true ;
	 string test = LineBuffer.substr(1,LineBuffer.length()-1);
         Z=atoi(test.c_str());
	 cout << "Z = " << Z<<endl;
      }
      else if (LineBuffer.compare(0,1, "A") == 0 && cA == false) {
         cA = true ;
	 string test = LineBuffer.substr(1,LineBuffer.length()-1);
         A=atoi(test.c_str());
	 cout << "A = " << A<<endl;
      }
}
   ConfigFile.close();
   return   ;
}


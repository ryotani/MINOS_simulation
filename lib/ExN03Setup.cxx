#include "ExN03Setup.hh"
ClassImp(ExN03Setup)

ExN03Setup::ExN03Setup(){
  //ReadConfigurationFile("ConfigFileSetup_ACTAR.txt");
  ReadConfigurationFile("ConfigFileSetup.txt");
}

ExN03Setup::~ExN03Setup(){
}

void ExN03Setup::ReadConfigurationFile(string Config)
{
   string LineBuffer;
   string DataBuffer;

   bool cTargetRadius       = false;
   bool cTargetLength       = false;
   bool cChamberInnerRadius       = false;
   bool cChamberThickness       = false;
   bool cChamberLength       = false;
   bool cInnerRohacellThickness       = false;
   bool cKaptonThickness       = false;
   bool cOuterRohacellThickness       = false;   
   bool cTPCRadiusExt       = false;
   bool cWindowThickness       = false;
   
   //string Path = "/Users/acorsi/codes/MINOS_simulation/Inputs/"+Config;
   string Path = "/home/local1/workspace/MINOS_simulation/Inputs/"+Config;
   ifstream ConfigFile;
   ConfigFile.open(Path.c_str());

   if (ConfigFile.is_open())
      cout << " Configuration file for setup " << Config << " loading " << endl;
   else {
      cout << " Error, no configuration file for setup " << Config << " found" << endl;
      return;
   }


   while (!ConfigFile.eof()) {
      //Pick-up next line
      getline(ConfigFile, LineBuffer);
      //Search for comment Symbol: %
      if (LineBuffer.compare(0, 1, "%") == 0) {   /*Do  Nothing*/;}

      else if (LineBuffer.compare(0,12, "TargetRadius") == 0 && cTargetRadius == false) {
         cTargetRadius = true ;
	 string test = LineBuffer.substr(12,LineBuffer.length()-12);
         TargetRadius=atof(test.c_str());
	 cout << "TargetRadius = " << TargetRadius<<" mm"<<endl;
      }
      else if (LineBuffer.compare(0,12, "TargetLength") == 0 && cTargetLength == false) {
         cTargetLength = true ;
	 string test = LineBuffer.substr(12,LineBuffer.length()-12);
         TargetLength=atof(test.c_str());
	 cout << "TargetLength = " << TargetLength<<" mm"<<endl;
      }
      else if (LineBuffer.compare(0,18, "ChamberInnerRadius") == 0 && cChamberInnerRadius == false) {
         cChamberInnerRadius = true ;
	 string test = LineBuffer.substr(18,LineBuffer.length()-18);
         ChamberInnerRadius=atof(test.c_str());
	 cout << "ChamberInnerRadius = " << ChamberInnerRadius<<" mm"<<endl;
      }
      else if (LineBuffer.compare(0,16, "ChamberThickness") == 0 && cChamberThickness == false) {
         cChamberThickness = true ;
	 string test = LineBuffer.substr(16,LineBuffer.length()-16);
         ChamberThickness=atof(test.c_str());
	 cout << "ChamberThickness = " << ChamberThickness<<" mm"<<endl;
      }
      else if (LineBuffer.compare(0,13, "ChamberLength") == 0 && cChamberLength == false) {
         cChamberLength = true ;
	 string test = LineBuffer.substr(13,LineBuffer.length()-13);
         ChamberLength=atof(test.c_str());
	 cout << "ChamberLength = " << ChamberLength<<" mm"<<endl;
      }
      else if (LineBuffer.compare(0,22, "InnerRohacellThickness") == 0 && cInnerRohacellThickness == false) {
         cInnerRohacellThickness = true ;
	 string test = LineBuffer.substr(22,LineBuffer.length()-22);
         InnerRohacellThickness=atof(test.c_str());
	 cout << "InnerRohacellThickness = " << InnerRohacellThickness<<" mm"<<endl;
      }
      else if (LineBuffer.compare(0,15, "KaptonThickness") == 0 && cKaptonThickness == false) {
         cKaptonThickness = true ;
	 string test = LineBuffer.substr(15,LineBuffer.length()-15);
         KaptonThickness=atof(test.c_str());
	 cout << "KaptonThickness = " << KaptonThickness<<" mm"<<endl;
      }
      else if (LineBuffer.compare(0,22, "OuterRohacellThickness") == 0 && cOuterRohacellThickness == false) {
         cOuterRohacellThickness = true ;
	 string test = LineBuffer.substr(22,LineBuffer.length()-22);
         OuterRohacellThickness=atof(test.c_str());
	 cout << "OuterRohacellThickness = " << OuterRohacellThickness<<" mm"<<endl;
      }
      else if (LineBuffer.compare(0,12, "TPCRadiusExt") == 0 && cTPCRadiusExt == false) {
         cTPCRadiusExt = true ;
	 string test = LineBuffer.substr(12,LineBuffer.length()-12);
         TPCRadiusExt=atof(test.c_str());
	 cout << "TPCRadiusExt = " << TPCRadiusExt<<" mm"<<endl;
      }
      else if (LineBuffer.compare(0,15, "WindowThickness") == 0 && cWindowThickness == false) {
         cWindowThickness = true ;
	 string test = LineBuffer.substr(15,LineBuffer.length()-15);
         WindowThickness=atof(test.c_str());
	 cout << "WindowThickness = " << WindowThickness<<" mm"<<endl;
      }
}
   ConfigFile.close();
   return   ;
}


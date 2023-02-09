#ifndef RECON_HCANA_H
#define RECON_HCANA_H

#include<iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

class recon_hcana
{
  
 public:

  //Consructor / Destructor
  recon_hcana();
  ~recon_hcana();
  
  void grabHistData(TString InSIMCHistname);

  //Auxiliary Function Prototypes (obtained from hcana) to calculate Pmx, Pmy, Pmz in the Lab/q-frame correctly
  void GeoToSph( Double_t  th_geo, Double_t  ph_geo, Double_t& th_sph, Double_t& ph_sph);
  void SetCentralAngles(Double_t th_cent, Double_t ph_cent);
  void TransportToLab( Double_t p, Double_t xptar, Double_t yptar, TVector3& pvect ); 
  
  //Utilities Functions for String Parsing
  string getString(char x);
  vector <string> FindString(TString keyword, TString fname);
  vector <string> split(string str, char del=':');
  
  void buildFileName(){

    string kinematics = "Q5p5W3p02_highe";
    string phi_setting = "Right";
    vector<string> kinematics_split;
    stringstream kinematics_stream(kinematics);
    string kinematics_part;

    while (getline(kinematics_stream, kinematics_part, '_')) {
      kinematics_split.push_back(kinematics_part);
    }

    transform(phi_setting.begin(), phi_setting.end(), phi_setting.begin(), [](unsigned char c) { return std::tolower(c); });

    InSIMCFilename = "../OUTPUTS/Prod_Coin_" + kinematics_split[0] + phi_setting + "_" + kinematics_split[1];

  }

  void ReadTree();
  void EventLoop();
  void WriteHist();

  TFile *f;
  TTree *tree;
  TTree *newTree;

  Int_t nentries;
  
  Float_t x;
  
  TString InSIMCFilename;
  TString InSIMCHistname;
  TString InSIMCRootname;

  Int_t simc_nevents;
  Double_t simc_normfactor;

  // Progress bar
  double progress=0.0;

  TString reaction = "heep";

  //Primary Kinematics (electron kinematics) (USED BY DATA AND SIMC)
  Double_t theta_e;              //Central electron arm angle relative to +z (hall coord. system)
  Double_t W;                    //Invariant Mass W (should be proton mass in H(e,e'p))
  Double_t W2;                    //Invariant mass squared
  Double_t Q2;                   //Four-momentum trasfer
  Double_t X;                    //B-jorken X  scaling variable
  Double_t nu;                   //Energy Transfer
  Double_t q;                  //magnitude of the 3-vector q
  Double_t th_q;                 //angle between q and +z (hall coord. system)

  //Secondary Kinematics (USED BY DATA AND SIMC)
  Double_t Ein;                  //single beam energy value (SIMC Uses this energy. If not corr. for energy loss, it should be same as in input file)
  Double_t Ep;                     //proton energy
  Double_t Em;                    //Standard Missing Energy for H(e,e'p)
  Double_t Em_nuc;                //Nuclear definition of Missing Energy (Used for D(e,e'p): B.E. of deuteron)
  Double_t Pm;                    //Missing Momentum (should be zero for H(e,e'p). Should be neutron momentum for D(e,e'p))
  Double_t Pmx_lab;               //X-Component of Missing Momentum (in Lab(or Hall) frame. +X: beam left, +Y: up, +Z: downstream beam) 
  Double_t Pmy_lab;
  Double_t Pmz_lab;
  Double_t Pmx_q;                 //X-Component of Missing Momentum (in frame where +z_lab is rotated to +z_q. Pmz_q is along +z(parallel to q))
  Double_t Pmy_q;
  Double_t Pmz_q;
  Double_t Kp;                    //Kinetic Energy of detected particle (proton)
  Double_t Kn;                    //Kinetic Energy of recoil system (neutron)
  Double_t M_recoil;              //Missing Mass (neutron Mass)
  Double_t MM2;                   //Missing Mass Squared
  Double_t E_recoil;              //Recoil Energy of the system (neutron total energy)
  Double_t En;                    //Same as above
  Double_t th_pq;                  //Polar angle of detected particle with q   ----> th_pq
  Double_t th_nq;                  //Polar angle of recoil system with q (rad)  ---> th_nq (neutreon-q angle. IMPORTANT in D(e,e'p))
  Double_t ph_pq;                  //Azimuth angle of detected particle with q    ----> phi_pq angle between proton and q-vector
  Double_t ph_nq;                  //Azimuth of recoil system with scattering plane (rad) ----> phi_nq angle between neutron and q-vector
  Double_t xangle;                //Angle of detected particle with scattered electron (Used to determine hadron angle)
  Double_t theta_p;               //to be calculated separately (in data)

  //Electron Arm Focal Plane / Reconstructed Quantities (USED BY DATA AND SIMC)
  Double_t e_xfp;
  Double_t e_xpfp;
  Double_t e_yfp;
  Double_t e_ypfp;
  
  Double_t e_ytar;
  Double_t e_yptar;
  Double_t e_xptar;
  Double_t e_delta;
  Double_t kf;                        //final electron momentum
  Double_t ki;                        //initial electron momentum

  //Hadron Arm Focal Plane / Reconstructed Quantities (USED BY DATA AND SIMC)
  Double_t h_xfp;
  Double_t h_xpfp;
  Double_t h_yfp;
  Double_t h_ypfp;
  
  Double_t h_ytar;
  Double_t h_yptar;
  Double_t h_xptar;
  Double_t h_delta;
  Double_t Pf;                 //final proton momentum
  
  //Target Quantities (tarx, tary, tarz) in Hall Coord. System (USED BY DATA AND SIMC)
  Double_t tar_x; //For SIMC ONLY (It is the same for HMS/SHMS)

  Double_t  htar_x;
  Double_t  htar_y;
  Double_t  htar_z;
  
  Double_t  etar_x;
  Double_t  etar_y;
  Double_t  etar_z;

  Double_t ztar_diff;

  //Light-Cone Momentum Variables
  Double_t PmPar;     //parallel component of recoil momentum relative to q-vector
  Double_t PmPerp;    //transverse component of recoil momentum relative to q-vector
  Double_t alpha_n;   //light-cone momentum fraction of the recoil neutron
  Double_t alpha;     //momentum fraction of struck nucleon (normalized such that: alpha + alpha_n = 2)

  //----------Variables Used in Auxiliary Functions--------------------------------------

  TRotation       fToLabRot;              //Rotation matrix from TRANSPORT to lab
  Double_t        fThetaGeo;              //In-plane geographic central angle (rad)
  Double_t        fPhiGeo;                //Out-of-plane geographic central angle (rad)
  Double_t        fThetaSph, fPhiSph;     //Central angles in spherical coords. (rad)
  Double_t        fSinThGeo, fCosThGeo;   //Sine and cosine of central angles
  Double_t        fSinPhGeo, fCosPhGeo;   // in geographical coordinates
  Double_t        fSinThSph, fCosThSph;   //Sine and cosine of central angles in 
  Double_t        fSinPhSph, fCosPhSph;   // spherical coordinates  

  //Declare Neccessary Variables to Determine the 4-Momentum of Recoil System
  TLorentzVector fP0;           // Beam 4-momentum
  TLorentzVector fP1;           // Scattered electron 4-momentum
  TLorentzVector fA;            // Target 4-momentum
  TLorentzVector fA1;           // Final system 4-momentum
  TLorentzVector fQ;            // Momentum transfer 4-vector
  TLorentzVector fX;            // Detected secondary particle 4-momentum (GeV)
  TLorentzVector fB;            // Recoil system 4-momentum (GeV)

  TVector3 Pf_vec;
  TVector3 kf_vec;

  //Declare necessary variables for rotaion from +z to +q
  TVector3 qvec;
  TVector3 kfvec;
  TRotation rot_to_q;
  TVector3 bq;   //recoil system in lab frame (Pmx, Pmy, Pmz)
  TVector3 xq;   //detected system in lab frame
  TVector3 p_miss_q;   //recoil system in q-frame

  //Leaf Variables
  Double_t fTheta_xq;
  Double_t fPhi_xq;
  Double_t fTheta_bq;
  Double_t fPhi_bq;  

  //Set Constants
  Double_t pi; 
  Double_t dtr;
  Double_t MP = 0.938272;     //proton mass
  Double_t MD = 1.87561;      //deuteron mass
  Double_t MN = 0.939566;     //neutron mass
  Double_t me = 0.00051099;   //electron mass
  Double_t MAL = 25.131710;   //aluminum mass
  Double_t tgt_mass;

  //Set Varibales to be read from REPORT_FILE
  Double_t e_th=0.;    //electron arm central angle
  Double_t e_ph=0.;    
  Double_t h_th=0.;    //hadron arm central angle
  Double_t h_ph=0.;

  Double_t xBPM;  //in cm
  Double_t yBPM;  //in cm
  
  Double_t e_xMisPoint;
  Double_t e_yMisPoint;
  Double_t h_xMisPoint;
  Double_t h_yMisPoint;

  //Central Spec. Momenta
  Double_t e_Pcen;
  Double_t h_Pcen;
  
};

#endif  //RECON_HCANA_H

// Testing signal events for process pp -> HZ (-> ejs, ll), from cms paper 2110.13218

#include <Pythia8/Pythia.h> 
#include <math.h>
#include <iostream>
#include <random>
#include <fstream>
#include <experimental/filesystem>
#include <numeric>
#include <algorithm>

using namespace Pythia8;

double cmedia (vector<double> &v){
  int size=v.size();
  sort(v.begin(), v.end());
  if(size % 2 != 0){
    return v[size/2];
  } 
  else{
    return (v[(size-1)/2] + v[size/2])/2.0;
  }
}
///compute transverse impact parameter of track
double d0t ( Particle* part){
  return abs(part->xProd()*part->py() - part->yProd()*part->px())/part->pT();
}
///compute lomgitudinal impact parameter
double d0z ( Particle* part){
  Vec4 a = part->vProd();
  Vec4 b = part->p();
  double dt = abs(part->xProd()*part->py() - part->yProd()*part->px())/part->pT();
  double d0 = sqrt(dot3(cross3(a,b),cross3(a,b))) / part->pAbs();
  double dz = sqrt(pow(d0,2)-pow(dt,2));
  return abs(dz);
}

///compute transverse angle between track direction vector from PV to innermost hit silicon tracker
double thetpx (Particle* part){
  // radial locations of silicon pixels at cms
  vector<double> pxd{44,73,102};
  double theta, rpx;
  double rprod = abs(sqrt(pow(part->xProd(),2) + pow(part->yProd(),2)));
  if (rprod<pxd[0]) rpx = pxd[0]; 
  else if (rprod<pxd[1]) rpx = pxd[1];
  else if (rprod<pxd[2]) rpx = pxd[2];
  else rpx = rprod;
  //cout << imother <<" "<< part->id() <<" -ev"<<iev << endl;
  double calph = ((part->xProd()*part->px()) + (part->yProd()*part->py())) / ((part->pT())*rprod);
  double alph = acos(calph);
  double gamm = M_PI - alph;
  double sthet = sin(gamm) * (rprod /rpx);
  theta = asin(sthet);
  return theta;
}
  
double thetpx3 (Particle* part){
  double rprod = abs(sqrt(pow(part->xProd(),2) + pow(part->yProd(),2)));
  double ctheta=((part->px()*part->xProd()) + (part->py()*part->yProd())) / ((part->pT()) * rprod);
  double theta = acos(ctheta);
  if (theta>M_PI) theta = (2*M_PI)-theta;
  return theta;
}

/// Generate random number from normal distribution
double randN (double const mu, double const sigma, int const seed){  
  //std::random_device rd;
  //std::mt19937 urn(rd());
  std::mt19937 urn(seed);
  std::normal_distribution<double> norm{mu,sigma};
  double value = norm(urn);
  double posval = abs(value);
  return posval;
}
/// Generate random number from binomial distribution
double randBool (double const pv, int const seed){
  //std::random_device rd;  
  //std::mt19937 urn(rd());
  std::mt19937 urn(seed);  
  std::bernoulli_distribution flip{pv};
  bool pass = flip(urn);
  return pass;
}

/// Read track reconstruction effiency as function of radial distance from IP, from iteration 5 (taken from Fig.12 in 1405.6569)
double readEff (double const rad){
  int bins=30;
  double wbin=2;
  int nb=-1;
  for (int i=0; i<bins; i++){
    double low = i*wbin;
    double upp = low + wbin;
    if (rad >= low && rad < upp){
      nb = i;
      break;
    }
  }
  ifstream fin;
  fin.open("Iter_5.csv", ios::in);
  //fin.open("Iter_4.csv", ios::in);
  vector<string> row;
  string line, val;
  double row2=0;
  int nrow=0;
  while (getline (fin, line)) {
    row.clear();
    stringstream s(line);
    while (getline(s, val, ',')) {
      row.push_back(val);
    }
    row2 = stod(row[1]);
    //if (rad >= (round(row1)-1) && rad <(round(row1)+1)) break;
    if (nrow == nb){
      goto end_fun;
    }
    nrow++;
  }
  if (fin.eof( ) && nb==-1) row2=0; 
  end_fun:
  fin.close();
  return row2;
}

/// Read resolution in d0_T of track as function of momentum and pseudorapidity (taken from Fig.15 in 1405.6569)
double resolIPpt (double const xp, double const eta){
  ifstream fin;
  fin.open("fig14_d0_resolution_pT.csv", ios::in);
  vector<string> rowL, rowU;
  string line, val;
  double x1=0,x2=1,y1=0.5,y2=0.5;
  double yp=0;
  int nrow=0;
  int ptl=0, resl=0;
  if (eta>=0 && eta<0.9){
    ptl=0;
    resl=1;
  }else if(eta>=0.9 && eta<1.4){
    ptl=4;
    resl=5;
  }else if(eta>=1.4 && eta<2.5){
    ptl=8;
    resl=9;
  }
  getline(fin, line);
  while (getline (fin, line)) {
    if (nrow>=2){
      rowL.clear();
      rowL = rowU;
      rowU.clear();
    }
    stringstream ss(line);
    while (getline(ss, val, ',')) {
      if (nrow==0){
        rowL.push_back(val);
      }else if (nrow==1){
        rowU.push_back(val);
      }else if (nrow>=2){
        rowU.push_back(val);
      }
    }
    if(nrow>=1){
      x1 = stod(rowL[ptl]);
      y1 = stod(rowL[resl]);
      x2 = stod(rowU[ptl]);
      y2 = stod(rowU[resl]);
      if (xp>=x1 && xp<x2){
        yp = y1 + ((y2-y1)/(x2-x1)) * (xp - x1);
        goto end_fun;
      }
    }   
    nrow++;
  }
  if (fin.eof( ) && xp>=x2 ) yp = y1 + ((y2-y1)/(x2-x1)) * (xp - x1);
  end_fun:
  fin.close();
  return yp * 0.001; //in [mm]
}

double etalayer(double const rl, double const zl){
  double theta = atan(rl / zl);
  double eta = - log(tan(theta/2.));
  return eta;
}

  // This function to determine how many tracker layers are hit

vector<bool> getBhitsV(Particle* part){
  /////###From andre code
  // Pixel
  // https://arxiv.org/pdf/0911.5434.pdf for 
  // Pixel barrel layers at e = 3, 6.8, 10.2, 16 cm
  // pixel z extends to |z| = 53.3/2 cm
  // This means |eta| < 1.3 is entirely inside the pixel barrel
  // Pre-upgrade: The endcap disks, extending from 6 to 15 cm in radius, are placed at z = ±35.5 cm and z = ±48.5 cm.
  // Post-upgrade: three endcap disks, see figure 2.1 in https://lss.fnal.gov/archive/design/fermilab-design-2012-02.pdf
  // let us guess that these are at 30, 40, 50 cm

  //1710.03842: Theradii of the four barrel layers are 2.9 cm, 6.8 cm, 10.9 cm, and 16.0 cm. 
  //The three forward disks arepositioned along z at 32 cm, 39 cm, and 48 cm. 
  // Each disk is composed of two rings of moduleswith average radii of 12.8 cm and 7.8 cm.

  vector<double> pixelbarrelR= {44.0,73.0,102.0}; //1
  //vector<double> pixelbarrelR= {30.0,68.0,102.0,160.0}; //2
  //vector<double> pixelbarrelR= {29.0,68.0,109.0,160.0}; //3 
  vector<double> pixelbarrelE= {2.5,2.0,1.6};
  //vector<double> pixelbarrelE= {2.5,2.1,1.6,1.3};

  vector<double> pixelForwardZ = {345.0,465.0}; //1
  //vector<double> pixelForwardZ = {300.0,400.0,500.0}; //2
  //vector<double> pixelForwardZ = {320.0,390.0,480.0}; //3
  vector<double> pixelForwardE= {1.5,1.8};
  //vector<double> pixelForwardE= {1.5,2.0,1.6};
  //vector<double> pixelForwardE= {60.0,150.0};

  // Now for tracker barrel, see CMS 2008 report for some of these
  /*
  The Tracker Inner Barrel (TIB) and Disks (TID) cover r < 55 cm and | z | < 118 cm, and are composed of four barrel layers,
  supplemented by three disks at each end. The Tracker Outer Barrel (TOB) covers r > 55 cm and | z | < 118 cm and 
  consists of six barrel layers. The Tracker EndCaps (TEC) cover the region 124 < | z | < 282 cm.
  Each TEC is composed of nine disks, each containing up to seven concentric
  rings of silicon strip modules, yielding a range of resolutions similar to that of the TOB.
  */
  //static vector<double> TIBR= {230.0, 300.0, 400.0,500.0};
  // See CMS 2008 report page 64, these exetend -700 to 700 mm
  vector<double> TIBR= {255.0, 339.0, 418.5,498.0};
  //vector<double> TIBZ= {0.0,700.0};
  
  // The TID± are assemblies of three disks placed in z between ±800 mm and ±900 mm. The
  // disks are identical and each one consists of three rings which span the radius from roughly 200 mm
  // to 500 mm.
  // Here I use values inferred from Figure 3.34 in the CMS 2008 report
  //static vector<double> TIDZ= {775.0,900.0,1025.0};
  vector<double> TIDZ= {775.0,900.0,1100.0};
  //vector<double> TIDR = {200.0,500.0};

  //static vector<double> TOBR = {600.0,700.0,800.0,890.0,980.0,1060.0};
  // See CMS 2008 report page 68
  static vector<double> TOBR = {608.0,692.0,780.0,868.0,965.0,1080.0};
  //vector<double> TOBZ = {0.0,1180.0}

  //// Inferred from Figure 3.34 in the CMS 2008 report
  vector<double> TECZ= {1250.0,1400.0,1550.0,1700.0,1950.0,2000.0,2225.0,2450.0,2700.0};
  //vector<double> TECZ = {200.0,1100.0}
  /////###

  //double ptp = part->pT(); 
  //double pzp = part->pz();
  double etap = abs(part->eta());
  double rprod = abs(sqrt(pow(part->xProd(),2) + pow(part->yProd(),2)));
  double zprod = abs(part->zProd());
  //double dprod = part->vProd();
  bool etamax = false;
  if (etap < 2.5) etamax = true;

  vector<bool> Bhits;
  //vector<bool> Dhits;
  int tothit=0;
  if (etamax){
    //// N_hits in pixel
    int pxhits=0;
    for (int i=0; i<int(pixelbarrelR.size()); i++){
      if ((rprod <= pixelbarrelR[i] && zprod <=270.0) && (etap <= pixelbarrelE[i])){
        pxhits++;
        Bhits.push_back(true);
        tothit++;
      }else{
        Bhits.push_back(false);
      }
    }
    //// N_hits in TIB, TOB
    int tibhits=0;
    for (int i=0; i<int(TIBR.size()); i++){
      double tibE = etalayer(TIBR[i], 700.0);
      if ((rprod <= TIBR[i] && zprod <700.0) && (etap <= tibE)){
        tibhits++;
        Bhits.push_back(true);
        tothit++;
      }else{
        Bhits.push_back(false);
      }
    }
    int tobhits=0;
    for (int i=0; i<int(TOBR.size()); i++){
      double tobE = etalayer(TOBR[i], 1100.0);
      if ((rprod <= TOBR[i] && zprod <1100.0) && (etap <= tobE)){
        tobhits++;
        Bhits.push_back(true);
        tothit++;
      }else{
        Bhits.push_back(false);
      }
    }
  }
  return Bhits;
}

vector<bool> getDhitsV(Particle* part){

  vector<double> pixelForwardZ = {345.0,465.0}; //1
  //vector<double> pixelForwardZ = {300.0,400.0,500.0}; //2
  //vector<double> pixelForwardZ = {320.0,390.0,480.0}; //3
  vector<double> pixelForwardE= {1.5,1.8};
  //vector<double> pixelForwardE= {1.5,2.0,1.6};
  //vector<double> pixelForwardE= {60.0,150.0};
  vector<double> TIDZ= {775.0,900.0,1100.0};
  //vector<double> TIDR = {200.0,500.0};
  vector<double> TECZ= {1250.0,1400.0,1550.0,1700.0,1950.0,2000.0,2225.0,2450.0,2700.0};
  //vector<double> TECZ = {200.0,1100.0}

  double etap = abs(part->eta());
  double rprod = abs(sqrt(pow(part->xProd(),2) + pow(part->yProd(),2)));
  double zprod = abs(part->zProd());
  //double dprod = part->vProd();
  
  vector<bool> Dhits;
  vector<bool> tDhits;
  int tothit=0;
  if (etap < 2.5){
  //// N_hits in pixel
    for (int i=0; i<int(pixelForwardZ.size()); i++){
      if ((zprod < pixelForwardZ[i] && rprod <150.0) && (etap > pixelForwardE[i])){
        Dhits.push_back(true);
        tothit++;
      }else{
        Dhits.push_back(false);
      }
    }
    //// N_hits in TID, TEC
    for (int i=0; i<int(TIDZ.size()); i++){  
      double tidE = etalayer(500.0,TIDZ[i]);   
      if ((zprod < TIDZ[i] && rprod <500.0) && (etap > tidE)){
        Dhits.push_back(true);
        tothit++;
      }else{
        Dhits.push_back(false);
      }
      //// check 3D hits
      double tid3d = etalayer(400.0,TIDZ[i]);
      if ((zprod < TIDZ[i] && rprod <400.0) && (etap > tid3d)){
        tDhits.push_back(true);
      }else{
        tDhits.push_back(false);
      }
    }
    for (int i=0; i<int(TECZ.size()); i++){  
      double tecE = etalayer(1100.0,TECZ[i]);   
      if ((zprod < TECZ[i] && rprod <1100.0) && (etap > tecE)){
        Dhits.push_back(true);
        tothit++;
      }else{
        Dhits.push_back(false);
      }
      //// check 3D hits
      double tecE1 = etalayer(400.0,TECZ[i]);
      double tecE2a = etalayer(600.0,TECZ[i]);
      double tecE2b = etalayer(760.0,TECZ[i]); //(estimated)
      if (i<6){
        if ((zprod < TECZ[i] && rprod <400.0) && (etap > tecE1)){
          tDhits.push_back(true);
        }else{
          tDhits.push_back(false);
        }
      }
      if ((zprod < TECZ[i] && rprod <760.0) && (etap < tecE2a && etap > tecE2b)){
        tDhits.push_back(true);
      }else{
        tDhits.push_back(false);
      }
    }
    for (int k=0; k<int(tDhits.size()); k++){
      Dhits.push_back(tDhits[k]);
    }
  }
  return Dhits;
}
////Nhits barrel (tot 13)
int getBhits(vector<bool> vbhit){
  int nBhits=0;
  if (vbhit.size()>0){
    for (int k=0; k<int(vbhit.size()); k++){
      if (vbhit[k]) nBhits++;
    }
  }
  return nBhits;
}
////3Dhits Barrel (tot 7, stage1)
int get3dB(vector<bool> vbhitx){
  vector<int> tdl{0,1,2,3,4,7,8};
  int tdhitB=0;
  if (vbhitx.size()>0){
    for (auto j : tdl){
      if (vbhitx[j]) tdhitB++;
    }
  }
  return tdhitB;
}
////Nhits disk (tot 14)
int getDhits(vector<bool> vdhit){
  int nDhits=0;
  if (vdhit.size()>0){
    for (int k=0; k<14; k++){
      if (vdhit[k]) nDhits++;
    }
  }
  return nDhits;
}
////3Dhits Disk (tot 20)
int get3dD(vector<bool> vdhitx){
  vector<int> tdd{0,1};
  for (int i=14; i<32; i++){
    tdd.push_back(i);
  }
  int tdhitD=0;
  if (vdhitx.size()>0){
    for (auto j : tdd){
      if (vdhitx[j]) tdhitD++; 
    }
  }
  return tdhitD;
}

bool HPtrack(const int nlay, const int n3dlay){
  bool hpPass=false;
  if (nlay >= 6 && n3dlay >=2) hpPass=true;
  return hpPass;
}


int main(int argc, char* argv[]) {

  Pythia pythia;
  Event& event = pythia.event;
  
  if (argc != 2) {
    cout << "Give input card. Usage ./main1002.x inputcard" << endl;
    return 0;
  }

  pythia.readFile(argv[1]);
  pythia.readString("Beams:idA = 2212");
  pythia.readString("Beams:idB = 2212");
  pythia.readString("Beams:eCM = 13000");
  //pythia.readString("Random:setSeed = on");
  //pythia.readString("Random:seed = 1");
  pythia.readString(" PDF:pSet = 20");
  
  //SlowJet slowJet( -1, 0.4, 40, 2.0, 1, 1);
  SlowJet slowJet( -1, 0.4, 35, 2.4, 1, 1);

  pythia.init();

  //Event trimmedevent;
  // Extract settings to be used in the main program.
  int nEvent   = pythia.mode("Main:numberOfEvents");
  //int nAbort   = pythia.mode("Main:timesAllowErrors");
  //int nEvent = 5000;
  int nAbort   = 30;
  
  ofstream myfile1;
  myfile1.open("datHP/jetsR_55ct10a_ipzb.dat", std::ios_base::app);
  ofstream myfile2;
  myfile2.open("datHP/jetsI5_55ct10a_ipzb.dat", std::ios_base::app);
  ofstream myfile3;
  myfile3.open("datHP/evinf_55ct10a_ipzb.dat", std::ios_base::app);

  //ofstream myfile7;
  //myfile7.open("dat3/dilep_s55ct10.dat", std::ios_base::app);

  //ofstream myfile10;
  //myfile10.open("dat4x/effs_HZdj_a.dat", std::ios_base::app);
     
  // set cuts on variables for DJs tagging
  double IPscut = 1.25;
  double alphcut = 0.45;
  double thetcut = -1.5;
  
  bool smear=true;

  int jsize;
  double npf[3]={0};
  double npf2[3]={0};
  double nEvs=0;
  double nEvs2=0;
  double nEvspt=0;

  int mS=-1;
  int ctS=-1;


  /// Begin event loop.
  int iAbort = 0;

  //int iEvent = 0;
  //while (iEvent < nEvent) {
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    // Generate events. Quit if failure.
    if (!pythia.next()) {
      if (++iAbort < nAbort) continue;
      cout << " Event generation aborted prematurely, owing to error!\n";
      break;
    }
    
    /*if (iEvent < 1) {
      pythia.info.list();
      pythia.event.list();
    }*/
    
    slowJet. analyze( pythia.event );	

    int ndpi=0;
    vector<int> idpi;

    vector<int> dtin;
    vector<double> dtr;
    vector<Vec4> dtp;
    vector<int> dtsj;
    vector<double> dtl3;

    vector<int> zlist;
    vector<int> inL;
    double invmL;
    double ptf;
    bool trgll=false;
    bool leta=false;
    bool lpt=false;
    bool trglpt=false;

    vector<int> inLi;
    int nll=0;
    bool lpair=false;

    int ndq=0;
    vector<int> idq;


    // loop over event constituents
    for (int i = 0; i < pythia.event.size(); ++i){

      if (event[i].isFinal() && (abs(event[i].id())==11 || abs(event[i].id())==13)){
        inLi.push_back(i);
        nll++;
      }

      /*if((abs(event[i].id()) ==  4900101 && abs(event[event[i].daughter1()].id())!=4900101)){// || (abs(event[i].id()) ==  4900021 && abs(event[event[i].mother1()].id()) == 4900001)){
        ndq++;
        idq.push_back(i);
      }*/
      if (abs(event[i].id()) ==  5 && abs(event[event[i].mother1()].id())==4900111){
        ndq++;
        idq.push_back(i);
      }
      
/////////// Find final daughters of LLP decaying into qq
      vector< int > dplist;
      if(abs(pythia.event[i].id()) ==  4900111 && event[event[i].daughter1()].id() != 4900111){
      //if(abs(pythia.event[i].id()) ==  54 && event[event[i].daughter1()].id() != 54){
        Particle part = pythia.event[i];
        ndpi++;   
        idpi.push_back(i);       
        double rpi = sqrt(pow(part.xDec(),2) + pow(part.yDec(),2));
        double lxyz = part.vDec().pAbs();
        mS = part.m0();
        ctS = part.tau0(); 

          //////////////
        dplist = part.daughterListRecursive();
        sort( dplist.begin(), dplist.end() );
        dplist.erase( unique( dplist.begin(), dplist.end() ), dplist.end() );          
        for (int j=0; j < dplist.size(); j++) {
          if(event[dplist[j]].isCharged() && event[dplist[j]].isFinal()){
            double rtp = sqrt(pow(event[dplist[j]].xProd(),2) + pow(event[dplist[j]].yProd(),2));
            if (event[dplist[j]].pT()>1.0 && event[dplist[j]].eta()<2.4){
              dtin.push_back(dplist[j]);                              
              Vec4 tp4 = event[dplist[j]].p();
              dtp.push_back(tp4);
              dtr.push_back(rtp);
              dtl3.push_back(event[dplist[j]].vProd().pAbs()); 
              dtsj.push_back(slowJet.jetAssignment(dplist[j]));
            }
          }
        }           
      }       
    /////end event-size          
    }

    vector<int> lls;
    vector<int> v(nll); 
    std::iota(v.begin(),v.end(),0); 
    sort( v.begin(),v.end(), [&](int i,int j){return event[inLi[i]].pT() > event[inLi[j]].pT();});
    for (int x=0; x<nll; x++){
      lls.push_back(inLi[v[x]]);
    }
    if (nll>2){
      if(event[lls[2]].pT()<15){
        for (int s=0; s<2; s++){
          inL.push_back(lls[s]);
        }
        lpair=true;
      }
    }else if (nll==2){
      inL = lls;
      lpair=true;
    }

    if(lpair && int(inL.size())==2 && (event[inL[0]].id() == -(event[inL[1]].id()))){
      Vec4 pl1 = event[inL[0]].p();
      Vec4 pl2 = event[inL[1]].p();
      invmL = (pl1 + pl2).mCalc();
      ptf = pl1.pT() + pl2.pT();

      //myfile7 << iEvent <<" "<< invmL<<" "<< event[inL[0]].id()<<" "<< pl1.pT()<<" "<< pl1.eta()<<" "<< pl1.phi()<<" "<< pl1.pAbs()/pl1.mCalc()<<" " \
          << event[inL[1]].id() <<" "<< pl2.pT()<<" "<< pl2.eta()<<" "<< pl2.phi()<<" "<< pl2.pAbs()/pl2.mCalc() << endl;

      /// look if leptons pass requirements on pt, eta, minv, pt_tot
      if (abs(event[inL[0]].id())==11){
        if (pl1.pT()>25 && pl2.pT()>15) lpt=true;
        if (abs(pl1.eta())<2.5 && abs(pl2.eta())<2.5) leta=true;
      }else if (abs(event[inL[0]].id())==13){
        if (pl1.pT()>25 && pl2.pT()>12) lpt=true;
        if (abs(pl1.eta())<2.4 && abs(pl2.eta())<2.4) leta=true;
      }
      if (lpt && leta && (invmL >70. && invmL <110.)) trgll=true;
      if (trgll &&  ptf>100) trglpt=true; 
    }

    //### Analyse jets    
    jsize = int(slowJet.sizeJet());
    double ptjja=0;
    int ejn=0, ejn2=0;
    vector<bool> ejp(jsize,false);
    vector<bool> ejp2(jsize,false);
    int ndj=0;     
    vector<double> mipj1(jsize,0);
    vector<double> mipj2(jsize,0);
    vector<double> mipjs1(jsize,0);
    vector<double> mipjs2(jsize,0);
    vector<double> alpha1(jsize,0);
    vector<double> alpha2(jsize,0);
    vector<double> mthet1(jsize,0);
    vector<double> mthet2(jsize,0);

    vector<int> jntr1(jsize,0);
    vector<int> jntr2(jsize,0);
    int jpass=0;
    int jpass2=0;

    int epass=0;
    vector<bool> dqinj(jsize,false);
    vector<bool> sqinj(jsize,false);
    vector<int> njdpi(jsize,0);
    vector<int> ndtk1(jsize,0);
    vector<int> ndtk2(jsize,0);
    vector<int> iej;
    vector<int> iej2;

    double alp1=1.2;
    double alp2=1.1;
    double beta=3;

    //### loop over jets
    for (int j=0; j < jsize; j++){
      vector<int> jpart;
      double ptall1=0;
      double ptpv1=0;
      double ptall2=0;
      double ptpv2=0;

      vector<double> ipj;
      vector<double> ipj2;
      vector<double> ipjs;
      vector<double> ipjs2;
      vector<double> theta;
      vector<double> theta2;
      
      double esum=0;
      bool emp=false;
      bool jin=false;
      bool jin2=false;
      bool isol=false;

      for (int m=0; m < int(idq.size()); m++){
        double dR = RRapPhi(slowJet.p(j),event[idq[m]].p());
        if (dR < 0.4) dqinj[j]=true;
      }
      for (int l=0; l < int(idpi.size()); l++){
        double dR = RRapPhi(slowJet.p(j),event[idpi[l]].p());
        if (dR < 0.4) njdpi[j]++;
      }

      //### loop oves jet constituents      
      jpart = slowJet.constituents(j);
      for (int t=0; t< int(jpart.size()); t++){
        int seed = iEvent+j+t;  
        double rprod = abs(sqrt(pow(event[jpart[t]].xProd(),2) + pow(event[jpart[t]].yProd(),2)));
        double dprod = abs(sqrt(pow(event[jpart[t]].xProd(),2) + pow(event[jpart[t]].yProd(),2) + pow(event[jpart[t]].zProd(),2)));
        double rpcm = rprod*0.1;
        double tpt = event[jpart[t]].pT();
        double teta = abs(event[jpart[t]].eta());
        bool tpass=false;
        double eftrk, ipa, z0a, thet;
        bool hptrk=false, passIPZ=false;
        /// select tracks following two approaches (for tracking) as function of transverse distance rprod
        if (event[jpart[t]].isCharged() && event[jpart[t]].isFinal() && event[jpart[t]].pT()>1.0){
          double z0i=d0z(&event[jpart[t]]);
          double ipi=d0t(&event[jpart[t]]);
          double resd0 = resolIPpt(tpt, teta); // in [mm]!
          double theti = thetpx(&event[jpart[t]]); // theta angle
          double sigip = abs(sqrt(pow(0.03,2) + pow((0.01/event[jpart[t]].pT()),2)));
          double sigd0 = resd0; // sigip or resd0
          /// count nHits layers (BorD??) ans check if track is HP
          vector<bool> bhits = getBhitsV(&event[jpart[t]]);
          vector<bool> dhits = getDhitsV(&event[jpart[t]]);
          int nhitB = getBhits(bhits);
          int nhitD = getDhits(dhits);
          int n3dhitB = get3dB(bhits);
          int n3dhitD = get3dD(dhits);
          int nlayer = nhitB + nhitD;
          int n3dlayer = n3dhitB + n3dhitD;////////Nhits
          hptrk = HPtrack(nlayer, n3dlayer);
          /// smearing
          if (smear){            
            ipa= randN(ipi,sigd0,seed); 
            z0a= randN(z0i,sigd0,seed); //or sigma=0.1
            thet= randN(theti,0.03,seed);
          }else{
            ipa = ipi;
            z0a = z0i;             
          }

          double ipasig=ipa/sigd0;
          double z0asig=z0a/sigd0;
          double DNa = sqrt(pow((z0a/0.1),2) + pow((ipasig),2));
          if (ipasig < pow((alp1*nlayer),beta) && z0asig < pow((alp2*nlayer),beta)) passIPZ=true; ///req. for HP track


          ///1. considering rprod<102mm (i.e. with hit in pixel tracker)
          if (rprod<102.0 && hptrk && passIPZ){ 
            ipj.push_back(log10(ipa));
            ipjs.push_back(log10(ipasig)); //tracks IP_significance  (sigip or sigd0)              
            theta.push_back(log10(thet));
            ptall1 += event[jpart[t]].pT();
            if (ipi < 0.1 ) ptpv1 += event[jpart[t]].pT(); ///req for track from PV??
            //if (z0a<25 && DNa<4) ptpv1 += event[jpart[t]].pT();
            jntr1[j]++;
            if (find(dtin.begin(), dtin.end(), jpart[t]) != dtin.end()) ndtk1[j]++; //ipjD.push_back(ipa);
          }
          if (rprod< 600){          
            eftrk = readEff(rpcm); //as rpcm in cm
            tpass = randBool(eftrk,seed);
          }
          ///2. considering tracking efficiency using Iteration5 in [1]  
          if (tpass && hptrk && passIPZ){
            ipj2.push_back(log10(ipa));
            ipjs2.push_back(log10(ipasig));             
            theta2.push_back(log10(thet));
            ptall2 += event[jpart[t]].pT();
            if (ipi < 0.1 ) ptpv2 += event[jpart[t]].pT();
            //if (z0a<25 && DNa<4) ptpv2 += event[jpart[t]].pT();
            jntr2[j]++;
            if (find(dtin.begin(), dtin.end(), jpart[t]) != dtin.end()) ndtk2[j]++;
          }  
        }
        //if (jpart[t]!=inL[0] || jpart[t]!=inL[1]) isol=true;
        //if(abs(event[jpart[t]].id())==11 || event[jpart[t]].id()==22) esum += event[jpart[t]].e();
        if((abs(event[jpart[t]].id())==11 || abs(event[jpart[t]].id()==13)) && event[jpart[t]].pT()>15) esum ++;

      }


      ///# claculate media(IP_trk) and alpha variables
      if (ptall1>1.){
        mipj1[j] = cmedia(ipj);
        mipjs1[j] = cmedia(ipjs);
        alpha1[j] = ptpv1 / ptall1;
        mthet1[j] = cmedia(theta);
        jin = true;        
      }
      if (ptall2>1.){
        mipj2[j] = cmedia(ipj2);
        mipjs2[j] = cmedia(ipjs2);
        alpha2[j] = ptpv2 / ptall2;
        mthet2[j] = cmedia(theta2);
        jin2 = true;       
      }
      
      //# check energy of jets from electrons and photons is less 90%
      //if (esum < (0.9*slowJet.p(j).e())){
      if (esum <1){
        epass++;
        emp=true;
      }

      // Apply cuts for DJet tagging, on log10(med(IP_sig)), alpha and log1(med(theta)) variables
      if (jin && emp){
        jpass++;
        ptjja += slowJet.pT(j);
        // #look if jet is EJ, passing IP, alpha & theta requirements
        if ((mipjs1[j])>IPscut && alpha1[j]<alphcut && (mthet1[j])>thetcut){
          ejp[j]=true;
          ejn++;
          iej.push_back(j);
        }
        myfile1 << iEvent <<" "<< j <<" "<< slowJet.pT(j) <<" "<< slowJet.phi(j) <<" "<< slowJet.p(j).eta() <<" "<< slowJet.m(j) <<" " \
          << alpha1[j] <<" "<< mipj1[j] <<" "<< mipjs1[j] <<" "<< mthet1[j]<<" "<< ejp[j] <<" "<< jntr1[j] <<" "<< ndtk1[j] <<" " \
          << njdpi[j]<<" "<< dqinj[j] << endl;
      }
      if (jin2 && emp){
        jpass2++;
        ptjja += slowJet.pT(j);
        // #look if jet is EJ, passing IP, alpha & theta requirements
        if ((mipjs2[j])>IPscut && alpha2[j]<alphcut && (mthet2[j])>thetcut){
          ejp2[j]=true;
          ejn2++;
          iej2.push_back(j);
        }
        myfile2 << iEvent <<" "<< j <<" "<< slowJet.pT(j) <<" "<< slowJet.phi(j) <<" "<< slowJet.p(j).eta() <<" "<< slowJet.m(j) <<" " \
          << alpha2[j] <<" "<< mipj2[j] <<" "<< mipjs2[j] <<" "<< mthet2[j]<<" "<< ejp2[j] <<" "<< jntr2[j] <<" "<< ndtk2[j]<<" " \
          << njdpi[j]<<" "<< dqinj[j] << endl;
      }
      
    //end loop over jets
    }
    
    double dRej=-1;
    double dRdq=-1;
    if (ejn >= 2){
      dRej = RRapPhi(slowJet.p(iej[0]),slowJet.p(iej[1]));
    }
    if (int(idq.size())>=2){
      dRdq = RRapPhi(event[idq[0]].p(),event[idq[1]].p());
    }

    vector<bool> ejf(3,false);
    vector<bool> evpf(3,false);
    bool trgZ = false;
    vector<bool> ejf2(3,false);
    vector<bool> evpf2(3,false);
    bool trgZ2 = false;


    if (ejn==0) ejf[0]=true;
    if (ejn==1) ejf[1]=true;
    if (ejn>=2) ejf[2]=true;

    if (ejn2==0) ejf2[0]=true;
    if (ejn2==1) ejf2[1]=true;
    if (ejn2>=2) ejf2[2]=true;
    
    if (trgll && jpass>0 && trglpt) trgZ=true;
    if (trgZ) nEvs++;
    if (trgll && jpass>0) nEvspt++;

    if (trgll && jpass2>0 && trglpt) trgZ2=true;
    if (trgZ2) nEvs2++;

    // look if event passes dileppton requirements and count #Djs
    for (int k=0; k<3; ++k){
      if (trgZ && ejf[k] ) evpf[k]=true;
      if (evpf[k]) npf[k]++;
      if (trgZ2 && ejf2[k]) evpf2[k]=true;
      if (evpf2[k]) npf2[k]++;
    }
    
    myfile3 << iEvent <<" "<< evpf[2] <<" "<< evpf2[2] <<" "<< trgll <<" "<< trglpt<<" "<< ptjja <<" "<< jpass <<" "<< jpass2 <<" "<< jsize <<" " \
      << ejn <<" "<< ejn2 <<" "<< dRej <<" "<< dRdq << endl;
    ///##$$

  // end event loop
  }


  double eff1t = npf[2]/nEvent;
  double eff2t = npf2[2]/nEvent;

//#################################
  cout.precision(4);
  cout << "###################"<<endl;
  cout <<" mS="<< mS<<" ct="<<ctS<< endl;

  cout << "TrgZ evs: " << int(nEvs) <<" out of "<<nEvent << endl;
  cout << "R:"<< endl;
  cout << "S eff (+2dj) wrt trgZ: "<< npf[2]/nEvs <<" ,wrt totev "<< eff1t << endl;
  cout << "I5:"<< endl;
  cout << "S eff (+2dj) wrt trgZ: "<< npf2[2]/nEvs2 <<" ,wrt totev "<< eff2t << endl;

  double xs = 880; //fb-1 xs(pp->HZ)
  double lumi = 117; //fb-1
  double BRt = 0.07 * 1; //BR(Z>em)*BR(h>inv) ///CMS takes BR(h>SS)~10%
  double ny1 = xs * lumi * BRt * eff1t;
  double ny2 = xs * lumi * BRt * eff2t;
  cout << "CMS EJ search with L=117fb-1" << endl;
  //cout << "For mX="<< mX<<"GeV, m_dpi=" << mN << "GeV, ctau_dpi="<< int(ctN) <<"mm" << endl; 
  cout << "S yield (BRh=1):: R: " << ny1 << " ,I5: "<< ny2<< endl;

  pythia.stat();
  //myfile1 << pythia.info.sigmaGen() <<" "<< pythia.info.sigmaErr() <<" "<< pythia.info.nAccepted() <<" $$#MXD#$$" << endl; 

  //myfile10 << mS <<" "<< int(ctS) <<" "<< int(nEvent) <<" "<< int(nEvs) <<" "<< int(nEvspt) <<" "<< eff1t  <<" "<< eff2t <<" 11"<< endl;  

  myfile1.close();
	myfile2.close();
  myfile3.close();
  //myfile7.close();
  //myfile10.close();
  return 0;
}

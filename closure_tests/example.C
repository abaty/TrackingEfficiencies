#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include "TCanvas.h"
#include "TError.h"
#include "TPad.h"
#include "TString.h"
#include "TRandom1.h"
#include "TH1F.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TString.h"
#include "TCut.h"
#include "TNtuple.h"
#include "TLine.h"
#include "../ntupler/trackTree.C"

void example(){
 TH1D::SetDefaultSumw2();

 //input file
 //TString directory="root://eoscms//eos/cms/store/group/phys_heavyions/yjlee/PbPbMC2014/";
 //TString infname="HiForest3_HydjetDrum_Pyquen_Dijet80_Embedded_d20140122_Track7_v1";
 //TString directory="/mnt/hadoop/cms/store/user/dgulhan/HIMC/Jet80/Track8_Jet19_STARTHI53_LV1/merged_50k/";
 //TString infname="HiForest_Pythia_Hydjet_Jet80_Track8_Jet19_STARTHI53_LV1_merged_forest_0"; 
 TString directory="/mnt/hadoop/cms/store/user/dgulhan/HIMC/Jet80/Track8_Jet19_STARTHI53_LV1/merged_300k_2/";
 TString infname="HiForest_Pythia_Hydjet_Jet80_Track8_Jet19_STARTHI53_LV1_merged_forest_0";


 trackTree * ftrk = new trackTree(Form("%s/%s.root",directory.Data(),infname.Data()));
 HiTree * fhi = new HiTree(Form("%s/%s.root",directory.Data(),infname.Data()));
 t * fjet = new t(Form("%s/%s.root",directory.Data(),infname.Data()));
 
 //pt bins for track efficiency correction
 const int nhist=14;
 //const int maxcentbins = 5;
 //int ncent[]={5,5,4,3};
 //double centmax[npt][maxcentbins]={{10,20,30,50,100},{10,20,30,50,100},{10,20,50,100,-1},{10,20,100,-1,-1}};
 //double centmin[npt][maxcentbins]={{0,10,20,30,50},{0,10,20,30,50},{0,10,20,50,-1},{0,10,20,-1,-1}};
 double ptmin[]={0.5,0.5,0.5,0.5,0.5,1,1,1,1,1,3,3,3,8};
 double ptmax[]={1,1,1,1,1,3,3,3,3,3,8,8,8,300};
 double centmin[]={0,10,20,30,50,0,10,20,30,50,0,10,20,0};
 double centmax[]={10,20,30,50,100,10,20,30,50,100,10,20,100,100};
 double cent_step[]={4,4,4,4,4,3,3,3,3,3,3,3,3,3};
 double accept_step[]={4,4,4,4,4,3,3,3,3,3,3,3,3,3};
 double pt_step[]={4,4,4,4,4,3,3,3,3,3,3,3,3,3};
 double rmin_step[]={3,3,3,3,3,2,2,2,2,2,2,2,2,2};

 //getting histograms for track efficiency correction 
 TFile *f_eff[nhist];
 TProfile *p_eff_cent[nhist]; 
 TProfile2D *p_eff_accept[nhist]; 
 TProfile *p_eff_pt[nhist]; 
 TProfile *p_eff_rmin[nhist]; 
 for(int ipt=0; ipt<nhist;ipt++){
   f_eff[ipt]= new TFile(Form("../final_hists_vsCalo/eff_pt%d_%d_cent%d_%d_step_cent%daccept%dpt%drmin%d.root",(int)ptmin[ipt],(int)ptmax[ipt],(int)centmin[ipt],(int)centmax[ipt],(int)cent_step[ipt],(int)accept_step[ipt],(int)pt_step[ipt],(int)rmin_step[ipt]));
   // if(ipt==2)f_eff[ipt]= new TFile(Form("../final_hists/eff_pt%d_%d_step_cent2accept2pt2rmin1.root",(int)ptmin[ipt],(int)ptmax[ipt]));
   // else f_eff[ipt]= new TFile(Form("../final_hists/eff_pt%d_%d_step_cent3accept3pt3rmin2.root",(int)ptmin[ipt],(int)ptmax[ipt]));
   p_eff_cent[ipt]=(TProfile*)f_eff[ipt]->Get("p_eff_cent");
   p_eff_pt[ipt]=(TProfile*)f_eff[ipt]->Get("p_eff_pt");
   p_eff_accept[ipt]=(TProfile2D*)f_eff[ipt]->Get("p_eff_acceptance");
   p_eff_rmin[ipt]=(TProfile*)f_eff[ipt]->Get("p_eff_rmin");
 }

 //output file and tree
 TFile *outf= new TFile(Form("/export/d00/scratch/abaty/trackingEff/closure_ntuples/track_ntuple_%s_vsCalo_300k_1.root",infname.Data()),"recreate");
 std::string trackVars="pt:mpt:eta:phi:rmin:trackselect:cent:eff";
 TNtuple *nt_track = new TNtuple("nt_track","",trackVars.data());

 //loop over events
 int nentries = ftrk->GetEntriesFast();
 //for(int jentry=0;jentry<nentries;jentry++){
 for(int jentry=0; jentry<10000;jentry++){	  
  if((jentry%1000)==0) std::cout<<jentry<<"/"<<nentries<<std::endl;

  ftrk->GetEntry(jentry);
  fhi->GetEntry(jentry);
  fjet->GetEntry(jentry);

  float cent=fhi->hiBin;
  //loop over tracks
  for(int itrk=0;itrk<ftrk->nParticle;itrk++){

   float trackselect=(ftrk->mtrkQual[itrk] && fabs((ftrk->mtrkDxy1[itrk]/ftrk->mtrkDxyError1[itrk]))<3.0 && fabs((ftrk->mtrkDz1[itrk]/ftrk->mtrkDzError1[itrk]))<3 && (ftrk->mtrkPtError[itrk]/ftrk->mtrkPt[itrk])<0.1);
   float eta=ftrk->pEta[itrk];

   if(fabs(eta)>2.4) continue; //acceptance of the tracker
   if(ftrk->pPt[itrk]<ptmin[0]) continue; //cuts out particles w/ pt< the smallest pt we have eff corrections for
   float pt=ftrk->pPt[itrk];
   float mpt=ftrk->mtrkPt[itrk];
   float phi=ftrk->pPhi[itrk];
   float rmin=99;
 
   //find rmin; 
   for(int ijet=0;ijet<fjet->nref;ijet++){
    if(fabs(fjet->jteta[ijet])>2 || fjet->jtpt[ijet]<30) continue;
    float r_reco=sqrt(pow(eta-fjet->jteta[ijet],2)+pow(acos(cos(phi-fjet->jtphi[ijet])),2));
    if(r_reco<rmin)rmin=r_reco;
   }
   
   //get efficiency correction for the track
   float eff_accept=1;
   float eff_pt=1;
   float eff_cent=1;
   float eff_rmin=1;
   
   for(int ipt=0;ipt<nhist;ipt++){
    if(pt>=ptmin[ipt] && pt<ptmax[ipt] && 0.5*cent>=centmin[ipt] && 0.5*cent<centmax[ipt]){
      eff_pt=p_eff_pt[ipt]->GetBinContent(p_eff_pt[ipt]->FindBin(pt));
      eff_cent=p_eff_cent[ipt]->GetBinContent(p_eff_cent[ipt]->FindBin(cent));
      eff_accept=p_eff_accept[ipt]->GetBinContent(p_eff_accept[ipt]->GetXaxis()->FindBin(phi),p_eff_accept[ipt]->GetYaxis()->FindBin(eta));
      if(rmin<5)eff_rmin=p_eff_rmin[ipt]->GetBinContent(p_eff_rmin[ipt]->FindBin(rmin));//efficiency for rmin>5 is 1.
     }     
   }
  
   float eff=eff_accept*eff_cent*eff_pt*eff_rmin;
   
   //fill in the output tree
   float entry[]={pt,mpt,eta,phi,rmin,trackselect,cent,eff};
   nt_track->Fill(entry);
  }
 }
 
  nt_track->Write();
  outf->Close();
}

#include <iostream>
#include "TCanvas.h"
#include "TError.h"
#include "TPad.h"
#include "TString.h"
#include "TRandom.h"
#include "TH1F.h"
#include "TF1.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TProfile.h"
#include "TLatex.h"
#include "TString.h"  
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCut.h"
#include "TNtuple.h"
#include "TProfile2D.h"
#include "TMath.h"
#include "TLine.h"

void plot_efficiency(int nstep_cent=2,int nstep_accept=2,int nstep_pt=2,int nstep_rmin=2,double  bin_pt_min=0.4,double bin_pt_max=1,double bin_cent_min=0, double bin_cent_max=10,int nevents=45000,bool is_final=false){
TH1D::SetDefaultSumw2();
TH2D::SetDefaultSumw2(true);
 // int nstep_cent=1;
 // int nstep_accept=1;
 // int nstep_pt=1;

// TFile * f = new TFile(Form("../workDir/tracking_ntuples/new/track_ntuple_Dijet100_HydjetDrum_v27_mergedV1_%devts_cent_%d_accept_%d_pt_%d_ptmin%d_ptmax%d.root",nevents,nstep_cent,nstep_accept,nstep_pt,bin_pt_min,bin_pt_max));
TFile * f = new TFile(Form("/afs/cern.ch/work/a/abaty/private/TrackingEfficiencies/ntuples/track_ntuple_HiForest_Pythia_Hydjet_Jet80_Track8_Jet19_STARTHI53_LV1_merged_forest_0_%devts_cent_%d_accept_%d_pt_%d_rmin_%d_ptmin%d_ptmax%d_centmin%d_centmax%d.root",nevents,nstep_cent,nstep_accept,nstep_pt,nstep_rmin,(int)bin_pt_min,(int)bin_pt_max,(int)bin_cent_min,(int) bin_cent_max));
TTree * nt = (TTree*)f->Get("nt_track");

//##########################################Initial Efficiencies without Correction###########################################################################

/////////eta dependent efficiency////////////////
//I believe this part of the code is now unused as of 2/24/2014
TProfile * p_eta= new TProfile("p_eta",";#eta;efficiency",20,-2.4,2.4);
TH1D * eta_gen = new TH1D("eta_gen",";#eta;N_{part}",20,-2.4,2.4);
TH1D * eta_reco = new TH1D("eta_reco",";#eta;N_{part}",20,-2.4,2.4);
nt->Draw("eta>>eta_gen","");
nt->Draw("eta>>eta_reco","mpt>0");
nt->Draw("(mpt>0 && trackselect):eta>>p_eta","pt>1");
eta_gen->SetMarkerStyle(25);

TLegend *t5=new TLegend(0.4,0.85,0.9,0.95); 
t5->SetFillColor(0);
t5->SetBorderSize(0);
t5->SetFillStyle(0);
t5->SetTextFont(43);
t5->SetTextSize(22);
t5->AddEntry(eta_gen,"gen","p");
t5->AddEntry(eta_reco,"matched gen","p");

TCanvas * c1 = new TCanvas("c1","",600,600);
p_eta->Draw();
c1->SaveAs("../plots/eff_eta.png");

TCanvas *c1_2 = new TCanvas("c1_2","",600,600);
eta_gen->Draw();
eta_gen->SetMinimum(eta_reco->GetMinimum());
eta_gen->SetMaximum(1.2*eta_gen->GetMaximum());
eta_reco->Draw("same");
t5->Draw("same");
c1_2->SaveAs("../plots/eta_dists.png");


///////////acceptance dependent efficiency/////////
TProfile2D * p_eta_phi = new TProfile2D(Form("p_eta_phi"),";#phi;#eta;",20,-TMath::Pi(),TMath::Pi(),20,-2.4,2.4);
nt->Draw(Form("(mpt>0 && trackselect):eta:phi>>p_eta_phi"),Form("((abs(eta)<2.4 && pt>%.1f && pt<%.1f)",bin_pt_min,bin_pt_max));

TCanvas * c2 = new TCanvas("c2","",600,600);
c2->SetRightMargin(0.2);
p_eta_phi->Draw("colz");
c2->SaveAs(Form("../plots/p_eta_phi.png"));


//////////pt dependent efficiency////////////

const int ny=20;
double x[ny+1];
double inix=log(bin_pt_min)/log(10);
double delta=(log(bin_pt_max)-log(bin_pt_min))/(20*log(10));
int maxbin=ny;
for(int ix=0; ix<ny+1;ix++){
 x[ix]=pow(10,inix);
 if(x[ix]>100){
	x[ix]=bin_pt_max;
	maxbin=ix;
	break;
 }
 inix+=delta;
} 

TProfile * p_pt= new TProfile("p_pt",";p_{T}(GeV/c);efficiency",maxbin,x);
TH1D * pt_gen = new TH1D("pt_gen",";p_{T}(GeV/c);N_{part}",maxbin,x);
TH1D * pt_reco = new TH1D("pt_reco",";p_{T}(GeV/c);N_{part}",maxbin,x);
nt->Draw("pt>>pt_gen","");
nt->Draw("pt>>pt_reco","mpt>0");
nt->Draw("(mpt>0 && trackselect):pt>>p_pt","");

TCanvas * c3 = new TCanvas("c3","",600,600);
c3->SetLogx();
p_pt->Draw();
c3->SaveAs("../plots/eff_pt.png");

TCanvas *c3_2 = new TCanvas("c3_2","",600,600);
c3_2->SetLogy();
pt_gen->SetMarkerStyle(25);
pt_gen->SetMinimum(pt_reco->GetMinimum());
pt_gen->SetMaximum(1.2*pt_gen->GetMaximum());
pt_gen->Draw();
pt_reco->Draw("same");
t5->Draw("same");
c3_2->SaveAs("../plots/pt_dists.png");


////////centrality dependent efficiency//////////

TProfile * p_cent= new TProfile("p_cent",";centrality bin;efficiency",40,0,40);
TH1D * cent_gen = new TH1D("cent_gen",";centrality bin;N_{part}",40,0,40);
TH1D * cent_reco = new TH1D("cent_reco",";centrality bin;N_{part}",40,0,40);
nt->Draw("cent>>cent_gen","");
nt->Draw("cent>>cent_reco","mpt>0");
nt->Draw("(mpt>0 && trackselect):cent>>p_cent",Form("(abs(eta)<2.4&& pt>%.1f && pt<%.1f)",bin_pt_min,bin_pt_max));

TH1D * eff_cent=(TH1D*)cent_reco->Clone("eff_cent");
eff_cent->Divide(cent_gen);

TCanvas * c4 = new TCanvas("c4","",600,600);
p_cent->Draw();
c4->SaveAs("../plots/eff_cent.png");


TCanvas *c4_2 = new TCanvas("c4_2","",600,600);
cent_gen->SetMarkerStyle(25);
cent_gen->SetMinimum(cent_gen->GetMinimum());
cent_gen->SetMaximum(1.2*cent_gen->GetMaximum());
cent_gen->Draw();
cent_reco->Draw("same");
t5->Draw("same");
c4_2->SaveAs("../plots/cent_dists.png");

//################################################rmin efficiency (calculated elsewhere, for bookkeeping here#################################################

const int n_rmin_bins=45;
double rmin_bins[n_rmin_bins+1];
for(int i=0; i<=n_rmin_bins; i++){
	if(i<=20) rmin_bins[i]=0.05*i;
	else if(i<=30) rmin_bins[i]=1+0.1*(i-20);
	else if(i<=45) rmin_bins[i]=2+0.2*(i-30);
}

TProfile * p_rmin = new TProfile("p_rmin",";#phi;#eta;efficiency",n_rmin_bins,rmin_bins);
nt->Draw("eff_rmin:rmin_reco>>p_rmin",Form("(abs(eta)<2.4&& pt>%.1f && pt<%.1f)",bin_pt_min,bin_pt_max),"prof");

//###############################################Efficiencies after Correction###########################################################################

///////////////cent dependent///////////////////
TProfile * p_cent_corr= new TProfile("p_cent_corr",";centrality bin;efficiency",200,0,200);
nt->Draw("(1/eff)*(mpt>0 && trackselect):cent>>p_cent_corr",Form("(abs(eta)<2.4&& pt>%.1f && pt<%.1f)",bin_pt_min,bin_pt_max));


TCanvas * c5 = new TCanvas("c5","",600,600);
p_cent_corr->SetMaximum(1.1);
p_cent_corr->SetMinimum(0.9);
p_cent_corr->Draw();
c5->SaveAs(Form("../plots/eff_cent_after_corr_nstep_cent%d_accept%d_pt%d.png",nstep_cent,nstep_accept,nstep_pt));


////////////pt dependent////////////////////////
TProfile * p_pt_corr= new TProfile("p_pt_corr",";p_{T}(GeV/c);efficiency",maxbin,x);
nt->Draw("(1/eff)*(mpt>0 && trackselect):pt>>p_pt_corr",Form("(abs(eta)<2.4&& pt>%.1f && pt<%.1f)",bin_pt_min,bin_pt_max));

TCanvas * c6 = new TCanvas("c6","",600,600);
c6->SetLogx();
p_pt_corr->SetMaximum(1.1);
p_pt_corr->SetMinimum(0.9);
p_pt_corr->Draw();
c6->SaveAs(Form("../plots/eff_pt_after_corr_nstep_cent%d_accept%d_pt%d_rmin%d_pt%d_%d.png",nstep_cent,nstep_accept,nstep_pt,nstep_rmin,(int)bin_pt_min,(int)bin_pt_max));
c6->SaveAs(Form("../plots/eff_pt_after_corr_nstep_cent%d_accept%d_pt%d_rmin%d_pt%d_%d.pdf",nstep_cent,nstep_accept,nstep_pt,nstep_rmin,(int)bin_pt_min,(int)bin_pt_max));

//////////acceptance dependent///////////////////

TProfile2D * p_eta_phi_corr = new TProfile2D("p_eta_phi_corr",";#phi;#eta;",20,-TMath::Pi(),TMath::Pi(),20,-2.4,2.4);
nt->Draw("(1/eff)*(mpt>0 && trackselect):eta:phi>>p_eta_phi_corr",Form("(abs(eta)<2.4&& pt>%.1f && pt<%.1f)",bin_pt_min,bin_pt_max));

TCanvas * c7 = new TCanvas("c7","",600,600);
c7->SetRightMargin(0.2);
p_eta_phi_corr->Draw("colz");
c7->SaveAs(Form("../plots/eff_accept_after_corr_nstep_cent%d_accept%d_pt%d_rmin%d_pt%d_%d.png",nstep_cent,nstep_accept,nstep_pt,nstep_rmin,(int)bin_pt_min,(int)bin_pt_max));
c7->SaveAs(Form("../plots/eff_accept_after_corr_nstep_cent%d_accept%d_pt%d_rmin%d_pt%d_%d.pdf",nstep_cent,nstep_accept,nstep_pt,nstep_rmin,(int)bin_pt_min,(int)bin_pt_max));

///////////eta dependent/////////////////////////
TProfile * p_eta_corr = new TProfile("p_eta_corr",";#eta;efficiency",20,-2.4,2.4);

nt->Draw("(1/eff)*(mpt>0 && trackselect):eta>>p_eta_corr",Form("(abs(eta)<2.4&& pt>%.1f && pt<%.1f)",bin_pt_min,bin_pt_max));

TCanvas * c8 = new TCanvas("c8","",600,600);
c8->SetRightMargin(0.2);
p_eta_corr->SetMaximum(1.1);
p_eta_corr->SetMinimum(0.9);
p_eta_corr->Draw();
c8->SaveAs(Form("../plots/eff_eta_after_corr_nstep_cent%d_accept%d_pt%d_rmin%d_pt%d_%d.png",nstep_cent,nstep_accept,nstep_pt,nstep_rmin,(int)bin_pt_min,(int)bin_pt_max));
c8->SaveAs(Form("../plots/eff_eta_after_corr_nstep_cent%d_accept%d_pt%d_rmin%d_pt%d_%d.pdf",nstep_cent,nstep_accept,nstep_pt,nstep_rmin,(int)bin_pt_min,(int)bin_pt_max));


////////////phi dependent/////////////////////////
TProfile * p_phi_corr = new TProfile("p_phi_corr",";#phi;efficiency",20,-TMath::Pi(),TMath::Pi());
nt->Draw("(1/eff)*(mpt>0 && trackselect):phi>>p_phi_corr",Form("(abs(eta)<2.4&& pt>%.1f && pt<%.1f)",bin_pt_min,bin_pt_max));

TCanvas * c9 = new TCanvas("c9","",600,600);
c9->SetRightMargin(0.2);

p_phi_corr->SetMaximum(1.1);
p_phi_corr->SetMinimum(0.9);
p_phi_corr->Draw();
c9->SaveAs(Form("../plots/eff_phi_after_corr_nstep_cent%d_accept%d_pt%d_rmin%d_pt%d_%d.png",nstep_cent,nstep_accept,nstep_pt,nstep_rmin,(int)bin_pt_min,(int)bin_pt_max));
c9->SaveAs(Form("../plots/eff_phi_after_corr_nstep_cent%d_accept%d_pt%d_rmin%d_pt%d_%d.pdf",nstep_cent,nstep_accept,nstep_pt,nstep_rmin,(int)bin_pt_min,(int)bin_pt_max));

////////////rmin dependent/////////////////////////
TProfile * p_rmin_corr = new TProfile("p_rmin_corr",";r_{min};efficiency",n_rmin_bins,rmin_bins);
nt->Draw("(1/eff)*(mpt>0 && trackselect):rmin_reco>>p_rmin_corr",Form("(abs(eta)<2.4&& pt>%.1f && pt<%.1f)",bin_pt_min,bin_pt_max));

TCanvas * c10 = new TCanvas("c10","",600,600);
c10->SetRightMargin(0.2);

p_rmin_corr->SetMaximum(1.1);
p_rmin_corr->SetMinimum(0.9);
p_rmin_corr->Draw();
c10->SaveAs(Form("../plots/eff_rmin_after_corr_nstep_cent%d_accept%d_pt%d_rmin%d_pt%d_%d.png",nstep_cent,nstep_accept,nstep_pt,nstep_rmin,(int)bin_pt_min,(int)bin_pt_max));
c10->SaveAs(Form("../plots/eff_rmin_after_corr_nstep_cent%d_accept%d_pt%d_rmin%d_pt%d_%d.pdf",nstep_cent,nstep_accept,nstep_pt,nstep_rmin,(int)bin_pt_min,(int)bin_pt_max));

/////output file to be used in the ntupler in the next step/////
TFile *outf = new TFile(Form("../stepwise_hists/eff_corr_nstep_cent%d_accept%d_pt%d_rmin%d_pt%d_%d_cent%d_%d.root",nstep_cent,nstep_accept,nstep_pt,nstep_rmin,(int)bin_pt_min,(int)bin_pt_max,(int)bin_cent_min,(int)bin_cent_max),"recreate");
p_cent->Write();
p_pt->Write();
p_eta_phi->Write();
p_cent_corr->Write();
p_pt_corr->Write();
p_eta_phi_corr->Write();
p_rmin_corr->Write();
p_rmin->Write();

outf->Close();

////overall efficiency histograms///////////////////////////////
TProfile * p_eff_cent = new TProfile("p_eff_cent",";centrality bin;efficiency",200,0,200);
nt->Draw("eff_cent:cent>>p_eff_cent",Form("((mpt<=0 || (trackselect && mpt>0))&& abs(eta)<2.4&& pt>%.1f && pt<%.1f)",bin_pt_min,bin_pt_max),"prof");

TProfile * p_eff_pt = new TProfile("p_eff_pt",";p_{T} bin;efficiency",maxbin,x);
nt->Draw("eff_pt:pt>>p_eff_pt",Form("((mpt<=0 || (trackselect && mpt>0))&& abs(eta)<2.4&& pt>%.1f && pt<%.1f)",bin_pt_min,bin_pt_max),"prof");

TProfile2D * p_eff_acceptance = new TProfile2D("p_eff_acceptance",";#phi;#eta;efficiency",20,-TMath::Pi(),TMath::Pi(),20,-2.4,2.4);
nt->Draw("eff_accept:eta:phi>>p_eff_acceptance",Form("((mpt<=0 || (trackselect && mpt>0))&& abs(eta)<2.4&& pt>%.1f && pt<%.1f)",bin_pt_min,bin_pt_max),"prof");

TProfile * p_eff_rmin = new TProfile("p_eff_rmin",";#phi;#eta;efficiency",n_rmin_bins,rmin_bins);
nt->Draw("eff_rmin:rmin_reco>>p_eff_rmin",Form("((mpt<=0 || (trackselect && mpt>0))&& abs(eta)<2.4&& pt>%.1f && pt<%.1f)",bin_pt_min,bin_pt_max),"prof");


TFile *f_efficiency;
 if(is_final){
 f_efficiency= new TFile(Form("../final_hists_puCalo/eff_pt%d_%d_cent%d_%d_step_cent%daccept%dpt%drmin%d.root",(int)bin_pt_min,(int)bin_pt_max,(int)bin_cent_min,(int)bin_cent_max,nstep_cent,nstep_accept,nstep_pt,nstep_rmin),"recreate");

 p_eff_cent->Write();
 p_eff_pt->Write();
 p_eff_acceptance->Write();
 p_eff_rmin->Write();

 f_efficiency->Close();
 f->Close();

//memory cleanup
/*
 delete p_eff_cent;
 delete p_eff_pt;
 delete p_eff_acceptance;
 delete p_eff_rmin;
 
 delete p_eta;
 delete eta_gen;
 delete eta_reco;
 delete p_eta_phi;
 delete pt_gen;
 delete pt_reco;
 delete p_pt;
 delete p_cent;
 delete cent_gen;
 delete cent_reco;
 delete eff_cent;
 delete p_rmin;

 delete p_cent_corr;
 delete p_pt_corr;
 delete p_eta_phi_corr;
 delete p_eta_corr;
 delete p_phi_corr;
 delete p_rmin_corr;
*/
}
}

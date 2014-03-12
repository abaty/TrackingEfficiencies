#include "../derivation/plot_efficiency.C"
#include "../ntupler/track_ntupler.C"

void run(){
 int nevents=50000;
 int npt=4;
 int ncent=5;
 double ptmin[]={0.5,1,3,8};
 double ptmax[]={1,3,8,300};
 double centmin[]={0,10,20,30,50};
 double centmax[]={10,20,30,50,100};

//changed these to arrays of size npt
//change these if npt changes!
 int ncent_step[]={5,5,5,4};
 int naccept_step[]={5,5,5,4};
 int npt_step[]={5,5,5,4};
 int nrmin_step[]={4,4,4,4};

//change to 0 normally
 cout<<"before the pt loop"<<endl;
 for(int ipt=0; ipt<npt;ipt++){
 	for(int icent=0; icent<ncent; icent++){
  		int icent_step=0;
  		int iaccept_step=0;
  		int ipt_step=0;
  		int irmin_step=0;
  		int istep=0;
		
		//changes the width of the high pt centrality bin (merges all bins together)
		int maxtemp = 100;
		if((ipt == npt-2) && (icent==2)){
			maxtemp = centmax[icent];
			centmax[icent] = centmax[ncent-1];			
			}
		if((ipt == npt-1) && (icent==0)){
			maxtemp = centmax[icent];
			centmax[icent] = centmax[ncent-1];
			}

  		while(icent_step<ncent_step[ipt] && iaccept_step<naccept_step[ipt] && ipt_step<npt_step[ipt] && irmin_step<nrmin_step[ipt]){
  			cout<<"\n\n icent_step= "<<icent_step<<" iaccept_step= "<<iaccept_step<<" ipt_step= "<<ipt_step<<" irmin_step= "<<irmin_step<<" pt_range= "<<ptmin[ipt]<<"-"<<ptmax[ipt]<<" cent_range= "<<centmin[icent]<<"-"<<centmax[icent]<<"\n"<<endl;
  		 	track_ntupler(icent_step,iaccept_step,ipt_step,irmin_step,ptmin[ipt],ptmax[ipt],centmin[icent],centmax[icent],nevents);
  		 	plot_efficiency(icent_step,iaccept_step,ipt_step,irmin_step,ptmin[ipt],ptmax[ipt],centmin[icent],centmax[icent],nevents);
  		 	if(istep%4==0) icent_step++;
  		 	if(istep%4==1) iaccept_step++;
  		 	if(istep%4==2) ipt_step++;
  		 	if(istep%4==3) irmin_step++;
  		 	istep++;
  		}
  		if((istep-1)%4==0) icent_step--;
  		if((istep-1)%4==1) iaccept_step--;
  		if((istep-1)%4==2) ipt_step--;
  		if((istep-1)%4==3) irmin_step--;
  		plot_efficiency(icent_step,iaccept_step,ipt_step,irmin_step,ptmin[ipt],ptmax[ipt],centmin[icent],centmax[icent],nevents,1);
		
		//terminates loop after 1 iteration only in the high pt bin
		if((ipt == npt-2) && (icent==2)){
			centmax[icent]=maxtemp;
			break;
			}
		if(ipt == npt-1){
 			centmax[icent]=maxtemp;
			break;
	 		}	
 }
 cout<<"after the pt loop"<<endl;
}
}


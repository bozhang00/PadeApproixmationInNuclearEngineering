#include <iostream>
#include <fstream>
#include <stdio.h>    
#include <vector>  
#include <stdlib.h>    
#include <time.h>
#include <string>
#include <algorithm>
#include <math.h>
#include <cmath>
#include "IC.h"
#include "INPUT.h"
#include "DIVIDELABLE.h"
#include "BIN.h"
#include "PRINT.h"
#include "RECOVERDATA.h"
using namespace std;
#define MATHLIB_STANDALONE
																											
int main() {
	
	int i, j, l, h, v, k, st,sv, record, note, LABLE_1 = 12, binh1, binh2, binh3, peak, bincheck=0, region_c;
	int	N_IC, STEP, BIN_I, BINSIZEA, BINSIZEB, BINSTEP;
	int PRE_PRE, POINT_N, N_H_ROW, N_V_COL, CAL_ROW, N_LABLE, BLOCK_COL;
	double ERROR_T, BIAS, REF_T, REF_DENSITY, dydt, dt, high_den, Densit_Threshold = 0.00000001;
	double sinke, sinkv1, sinkv2, sinkv3, sinkv4, sinkh1, sinkh2, sinkh3;
	clock_t t1, t2; 	  t1 = clock();
	
	double GE, GV1, GV2, GV3, GV4, GV5, GV6, GV7, GV8, GV9, GH, DE ,DV1, DV2,DV3,DV4,DH1, DH2, DH3;
	double e, c, cb, v1,v2, v3,v4,h1,h2,h3,vfc, vf1, vf2, vf3, vf4, vb1, hfc,hf1,hf2,hf3,hb1,te,tv1,tv2,tv3,tv4,th1,th2,th3;
	double kvf1c, kvf2c, kvf3c, kvf4c, kvf1, kvf2, kvf3, kvf4, khf1c, khf2c, khf3c, khf1, khf2, khf3, ke, keb; 
	double per;
	
	vector<vector<vector<double> > > a, av1, av2, av3, av4, ah1, ah2, ah3, ae, pa, pb, KV, KH;
	vector<vector<double> > KE, KBV, KBH, KBE, FULL;
	vector<double> high_density, se, sv1, sv2, sv3, sv4, sh1, sh2, sh3;
	vector<vector<int> > s, indicator;
	vector<int> L;
	vector<double> vtemp;
    ifstream infile;
		
	INPUTIC(N_IC, BIAS, BIN_I, BINSIZEA, BINSIZEB, BINSTEP, POINT_N, STEP, ERROR_T, REF_T, REF_DENSITY, dt, record);
	
	REACTION(N_IC, REF_DENSITY, REF_T, KV, KH, KE, KBV, KBH, KBE);
	
	DENSITY(N_IC, a, ae, av1, av2, av3, av4, ah1, ah2, ah3, pa, pb, L, s, se, sv1, sv2, sv3, sv4, sh1, sh2, sh3);
	
	GENERATION(GE, GV1, GV2, GV3, GV4, GV5, GV6, GV7, GV8, GV9, GH, REF_T, REF_DENSITY); 
	
	DIFFUSION(DE, DV1, DV2, DV3, DV4, DH1, DH2, DH3); 
	
	st = 0.0;		N_H_ROW = a.size();  	N_LABLE = 0;

	a[0][0][0] = a[0][0][0] + GE * dt; 		 			 
	a[0][0][1] = a[0][0][1] + GV1 * dt;		 			
	a[0][0][2] = a[0][0][2] + GV2 * dt;		 			
	a[0][0][3] = a[0][0][3] + GV3 * dt; 			
	a[0][0][4] = a[0][0][4] + GV4 * dt;			
	a[0][0][5] = a[0][0][5] + GV5 * dt;
	a[0][0][6] = a[0][0][6] + GV6 * dt;
	a[0][0][7] = a[0][0][7] + GV7 * dt;
	a[0][0][8] = a[0][0][8] + GV8 * dt;
	a[0][0][9] = a[0][0][9] + GV9 * dt;
	a[1][0][0] = a[1][0][0] + GH * dt;
		
	CAL_ROW = 2;
	indicator.clear(); 		
	indicator.resize(N_H_ROW); 		
	
	for(h =0; h<N_H_ROW; h++ )	indicator[h].resize(1);	

	if(record == 1){

		infile.open("data.txt");
        double tmp;
        while(!infile.eof())
        {
            infile>>tmp;
            vtemp.push_back(tmp);
        } 
        infile.close();
		//---------------------------------------------
		note = 0;
     	for(h =0; h<N_H_ROW; h++){
			for( l = 0; l<N_LABLE; l++ ) {	
				N_V_COL = a[h][l].size(); 
				for(v =0; v< N_V_COL ; v++) {
					a[h][l][v] =  vtemp[note];
					note ++;
				}						
			}
		}
		//-----------------------------------------------
	}
	
	
	do{ 
		st++; 			
		
		te = 0.0; 		
		tv1 = 0.0; 		 tv2 = 0.0;  	tv3 = 0.0;  	tv4 = 0.0;		th1 = 0.0; 			 th2 = 0.0; 	th3 = 0.0;
		
		e  = a[0][0][0];  	v1 = a[0][0][1];  	v2 = a[0][0][2];  	 v3 = a[0][0][3];  		v4 = a[0][0][4];   
		h1 = a[1][0][0];	h2 = a[2][0][0];	h3 = a[3][0][0];
		
		N_LABLE = a[0].size();
		
		if(st % 1000 ==0 ){	       // check the Helium size and decide the number of the raw......
			high_density.clear(); 			
			high_density.resize(N_H_ROW);
			for( h = 0; h < N_H_ROW; h++){
				high_density[h] = 0.0;
				for( l = 0; l<N_LABLE; l++ ) {	
					if(indicator[h][l] == 0) {  	
						N_V_COL = a[h][l].size();
						for(v =0; v< N_V_COL ; v++){
							 if( high_density[h] < a[h][l][v] )	high_density[h] = a[h][l][v];
						}
					}
				}
			}
			CAL_ROW = 1;
			for( h = 1; h<N_H_ROW; h++){ if( high_density[h] > Densit_Threshold ) CAL_ROW = h+1; }
			CAL_ROW ++ ;
			if( CAL_ROW > N_H_ROW ) CAL_ROW = N_H_ROW;
			
			high_den = 0.0;
			l=N_LABLE-1;
			for( h = 0; h < N_H_ROW; h++){
				if ( high_den < a[h][l][a[h][l].size()-3] ) high_den = a[h][l][a[h][l].size()-3];	
			}
			if(high_den > Densit_Threshold){
				for( h =0; h<N_H_ROW; h++) {
				a[h][N_LABLE-1].push_back(0.0);
				ae[h][N_LABLE-1].push_back(0.0);
				av1[h][N_LABLE-1].push_back(0.0); 
				av2[h][N_LABLE-1].push_back(0.0); 
				av3[h][N_LABLE-1].push_back(0.0); 
				av4[h][N_LABLE-1].push_back(0.0); 	
				
				ah1[h][N_LABLE-1].push_back(0.0);
				ah2[h][N_LABLE-1].push_back(0.0);
				ah3[h][N_LABLE-1].push_back(0.0);
				}
				s[N_LABLE-1].push_back( s[N_LABLE-1][s[N_LABLE-1].size()-1]+1 );
				L[N_LABLE - 1] = s[N_LABLE-1][s[N_LABLE-1].size()-1];
			}
		}		
	
		for(h =0; h<CAL_ROW; h++){	
			N_V_COL = a[h][0].size();
			#pragma omp parallel for			
			for(v =0; v< N_V_COL ; v++){
			
				c = a[h][0][v];	
				ke = KE[h][v]; 	
				keb = KE[h][v+1];
				
				if(v == N_V_COL-1){ 
					if(a[h].size() ==1 ) cb= 0.0;
					else cb = a[h][1][0];
				}
				else{
					cb = a[h][0][v+1];
				}
				
				kvf1c = KV[0][h][v]; 	if(v - 1 < 0) { kvf1 = 0.0; } else { kvf1 = KV[0][h][v-1]; }
				kvf2c = KV[1][h][v];	if(v - 2 < 0) { kvf2 = 0.0; } else { kvf2 = KV[1][h][v-2]; }
				kvf3c = KV[2][h][v];	if(v - 3 < 0) { kvf3 = 0.0; } else { kvf3 = KV[2][h][v-3]; }	
				kvf4c = KV[3][h][v];	if(v - 4 < 0) { kvf4 = 0.0; } else { kvf4 = KV[3][h][v-4]; }
				
				khf1c = KH[0][h][v];	if(h - 1 < 0) { khf1 = 0.0; } else { khf1 = KH[0][h-1][v]; }
				khf2c = KH[1][h][v]; 	if(h - 2 < 0) { khf2 = 0.0; } else { khf2 = KH[1][h-2][v]; }	
				khf3c = KH[2][h][v]; 	if(h - 3 < 0) { khf3 = 0.0; } else { khf3 = KH[2][h-3][v]; }	
				
				if(v < 4 ){
					switch (v){
					case 0: vf1 = 0.0;				vf2 = 0.0; 				vf3 = 0.0;				vf4 = 0.0;	break;
					case 1:	vf1 = a[h][0][v-1];		vf2 = 0.0; 				vf3 = 0.0;				vf4 = 0.0;	break;
					case 2: vf1 = a[h][0][v-1];		vf2 = a[h][0][v-2];		vf3 = 0.0;				vf4 = 0.0;	break;	
					case 3: vf1 = a[h][0][v-1];		vf2 = a[h][0][v-2];	 	vf3 = a[h][0][v-3];		vf4 = 0.0;	break;
					}
				}
				else{ 	vf1 = a[h][0][v-1];		vf2 = a[h][0][v-2];	    	vf3 = a[h][0][v-3]; 	vf4 = a[h][0][v-4];	}				
		
				if(h < 3 ){
								hf1 = 0.0;  	  		 hf2 = 0.0; 		 	hf3 = 0.0;
					switch (h){
					case 1: 	hf1 = a[h-1][0][v]; 							break;
					case 2:  	hf1 = a[h-1][0][v];		hf2 = a[h-2][0][v];		break;
					}
				}
				else{ 	hf1 = a[h-1][0][v]; 	hf2 = a[h-2][0][v]; 		hf3 = a[h-3][0][v]; 	}
				
				av1[h][0][v]= kvf1*vf1 - kvf1c*c;		
				av2[h][0][v]= kvf2*vf2 - kvf2c*c;			
				av3[h][0][v]= kvf3*vf3 - kvf3c*c;				
				av4[h][0][v]= kvf4*vf4 - kvf4c*c;
				
				ah1[h][0][v]= khf1*hf1 - khf1c*c; 						
				ah2[h][0][v]= khf2*hf2 - khf2c*c; 				
				ah3[h][0][v]= khf3*hf3 - khf3c*c;		
				
				ae[h][0][v] = keb*cb - ke*c;
				
				te = te - c*ke;
				tv1= tv1 - c*kvf1c; 						
				tv2= tv2 - c*kvf2c;
				tv3= tv3 - c*kvf3c;						
				tv4= tv4 - c*kvf4c;
				th1= th1 - c*khf1c;						
				th2= th2 - c*khf2c;
				th3= th3 - c*khf3c;
			}
			te = te + se[h];
			tv1= tv1 + sv1[h];	
			tv2= tv2 + sv2[h];	 
			tv3= tv3 + sv3[h];		
			tv4= tv4 + sv4[h];		
			th1= th1 + sh1[h];	
			th2= th2 + sh2[h];    
			th3= th3 + sh3[h];
		}
		
		{        // ******************** update the moving terms and source term ************************
		
		a[0][0][0] = a[0][0][0] + GE * dt; 		 			 
		a[0][0][1] = a[0][0][1] + GV1 * dt;		 			
		a[0][0][2] = a[0][0][2] + GV2 * dt;		 			
		a[0][0][3] = a[0][0][3] + GV3 * dt; 			
		a[0][0][4] = a[0][0][4] + GV4 * dt;			
		a[0][0][5] = a[0][0][5] + GV5 * dt;
		a[0][0][6] = a[0][0][6] + GV6 * dt;
		a[0][0][7] = a[0][0][7] + GV7 * dt;
		a[0][0][8] = a[0][0][8] + GV8 * dt;
		a[0][0][9] = a[0][0][9] + GV9 * dt;
		a[1][0][0] = a[1][0][0] + GH * dt;	
		
		per = 1.0;  if(te != 0 & e !=0 )   {	per = a[0][0][0]/( - te * e * dt); 		if( per >1 ){ per = 1.0; } 	}   e = per *e;
		per = 1.0;  if(tv1 != 0 & v1 !=0 ) {	per = a[0][0][1]/( - tv1 * v1 * dt); 	if( per >1 ){ per = 1.0; } 	}   v1 = per *v1;
		per = 1.0; 	if(tv2 != 0 & v2 !=0 ) {	per = a[0][0][2]/( - tv2 * v2 * dt); 	if( per >1 ){ per = 1.0; }	}   v2 = per *v2;
		per = 1.0; 	if(tv3 != 0 & v3 !=0 ) {	per = a[0][0][3]/( - tv3 * v3 * dt); 	if( per >1 ){ per = 1.0; } 	}   v3 = per *v3;
		per = 1.0; 	if(tv4 != 0 & v4 !=0 ) {	per = a[0][0][4]/( - tv4 * v4 * dt); 	if( per >1 ){ per = 1.0; } 	}   v4 = per *v4;
		per = 1.0;	if(th1 != 0 & h1 !=0 ) {	per = a[1][0][0]/( - th1 * h1 * dt);	if( per >1 ){ per = 1.0; } 	}	h1 = per *h1; 
		per = 1.0; 	if(th2 != 0 & h2 !=0 ) {	per = a[2][0][0]/( - th2 * h2 * dt); 	if( per >1 ){ per = 1.0; } 	}	h2 = per *h2;
		per = 1.0; 	if(th3 != 0 & h3 !=0 ) {	per = a[3][0][0]/( - th3 * h3 * dt); 	if( per >1 ){ per = 1.0; } 	}	h3 = per *h3; 
		
		
		ae[0][0][0] = ae[0][0][0] + te;
		av1[0][0][1] = av1[0][0][1] + tv1;
		av2[0][0][2] = av2[0][0][2] + tv2;
		av3[0][0][3] = av3[0][0][3] + tv3;
		av4[0][0][4] = av4[0][0][4] + tv4;
		ah1[1][0][0] = ah1[1][0][0] + th1;
		ah2[2][0][0] = ah2[2][0][0] + th2;
		ah3[3][0][0] = ah3[3][0][0] + th3;
		
		}

	
		for( h =0; h < CAL_ROW; h++){ 
			N_V_COL = a[h][0].size(); 
			#pragma omp parallel for
			for(v =0; v< N_V_COL ; v++)	 { 
				a[h][0][v] = a[h][0][v]+( v1*av1[h][0][v]+v2*av2[h][0][v]+v3*av3[h][0][v]+v4*av4[h][0][v]+
										  e*ae[h][0][v]+
										  h1*ah1[h][0][v]+h2*ah2[h][0][v]+h3*ah3[h][0][v])*dt; 
										  
				if( a[h][0][v] < 0.0 ) a[h][0][v] = 0.0;
			}
		}
		
	// Binning Region
		if(st%BINSTEP == 0 && L[N_LABLE - 1]>80 && BIN_I == 1 ){
			
			FULL.clear(); FULL.resize(N_H_ROW); for(h=0; h<N_H_ROW; h++){ RECOVER_KEEP_DATA(FULL[h], a[h], s, pa[h], pb[h], indicator[h] );	}
		
			NEW_LABLE(L, LABLE_1, s, BINSIZEB, FULL[0].size());
			
			PADESTEP(FULL, KE, KV, KH, a, se, sv1, sv2, sv3, sv4, sh1, sh2, sh3, dt, BINSTEP, CAL_ROW, L,
					 GE, GV1, GV2, GV3, GV4, GV5, GV6, GV7, GV8, GV9, GH);
			
			ae.clear();     
			av1.clear();	 		 av2.clear(); 	 		av3.clear(); 			 av4.clear();	 
			ah1.clear();	 		 ah2.clear();			 ah3.clear();
			pa.clear(); 			 pb.clear();
			
			ae.resize(N_H_ROW);
			av1.resize(N_H_ROW);	 av2.resize(N_H_ROW); 	 av3.resize(N_H_ROW); 	 av4.resize(N_H_ROW);	 
			ah1.resize(N_H_ROW);	 ah2.resize(N_H_ROW);	 ah3.resize(N_H_ROW);
			pa.resize(N_H_ROW); 	 pb.resize(N_H_ROW);
			
			indicator.clear();  
			indicator.resize(N_H_ROW);
			
			N_LABLE = L.size();
			
			for(h =0; h<N_H_ROW; h++ ){
				indicator[h].resize(N_LABLE);
				ae[h].resize(N_LABLE);
				av1[h].resize(N_LABLE);	 av2[h].resize(N_LABLE); 		 av3[h].resize(N_LABLE); 	 	 av4[h].resize(N_LABLE);	 
				ah1[h].resize(N_LABLE);	 ah2[h].resize(N_LABLE);	 	 ah3[h].resize(N_LABLE);
				pa[h].resize(N_LABLE);	 pb[h].resize(N_LABLE);
				#pragma omp parallel for
				for(l=0; l< N_LABLE; l++){
					N_V_COL = a[h][l].size(); 
					ae[h][0].resize(N_V_COL);
					av1[h][0].resize(N_V_COL);		 av2[h][0].resize(N_V_COL); 	 av3[h][0].resize(N_V_COL); 	 av4[h][0].resize(N_V_COL);	 
					ah1[h][0].resize(N_V_COL);	 	 ah2[h][0].resize(N_V_COL);		 ah3[h][0].resize(N_V_COL);					
				}	
				
			}
			
			for(h =0; h<CAL_ROW; h++ ){
				for(l=1; l< N_LABLE-1; l++) GET_PADE(a[h][l], s[l], pa[h][l], pb[h][l], POINT_N,REF_T, indicator[h][l] );
			}			
		}
		
		if (st % outputstep(st) == 0) 	print(a, s, pa, pb, st, indicator); 
	}while(st < STEP);	
	
	t2= clock();
	float diff_time((float)t2 - (float)t1);   
	float spend_time = diff_time/ CLOCKS_PER_SEC;
	
	finalprint(N_IC, BIAS, BIN_I, BINSIZEA, BINSIZEB, BINSTEP, POINT_N, STEP, 
				ERROR_T, REF_T, REF_DENSITY, dt, spend_time, a, s, pa, pb, st, indicator);
	
	cout<<"the total time is:  "<<spend_time<<endl;
	return 0;
}


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
#include <omp.h>
using namespace std;
#define MATHLIB_STANDALONE


void GENERATION(double& GE, double& GV1, double& GV2, double& GV3, double& GV4, double& GV5, double& GV6, double& GV7,
				double& GV8, double& GV9, double& GH)
{
	GE  = 1.49/pow(10, 5);
	GV1 = 9.91/pow(10, 6);
	GV2 = 1.51/pow(10, 6);
	GV3 = 2.60/pow(10, 7);
	GV4 = 1.58/pow(10, 7);
	GV5 = 6.29/pow(10, 8);
	GV6 = 0.0;
	GV7 = 0.0;
	GV8 = 0.0;
	GV9 = 3.16/pow(10, 8);
	GH  = 2.11/pow(10, 11);
}	

void DIFFUSION(double& DE, double& DV1, double& DV2, double& DV3, double& DV4, double& DH1, double& DH2, double& DH3)
{
	double Ev_KB_T = (16202.1/(1.38*698.15));
	
	DE   = 4.0 * 3.1415926 * 1.00 * pow(10, 11) * exp(-0.34*Ev_KB_T);
	DV1  = 4.0 * 3.1415926 * 1.00 * pow(10, 11) * exp(-0.67*Ev_KB_T);
	DV2  = 4.0 * 3.1415926 * 0.50 * pow(10, 11) * exp(-0.62*Ev_KB_T);	
	DV3  = 4.0 * 3.1415926 * 0.33 * pow(10, 11) * exp(-0.37*Ev_KB_T);	
	DV4  = 4.0 * 3.1415926 * 0.25 * pow(10, 11) * exp(-0.48*Ev_KB_T);	
	DH1  = 4.0 * 3.1415926 * 1.00 * pow(10, 11) * exp(-0.06*Ev_KB_T);	
	DH2  = 4.0 * 3.1415926 * 0.50 * pow(10, 11) * exp(-0.06*Ev_KB_T);	
	DH3  = 4.0 * 3.1415926 * 0.33 * pow(10, 11) * exp(-0.06*Ev_KB_T);
	
}	

																											
int main() {
	
	int i, j, h, v, st = 0, N_IC, N_ICh, N_ICv,  first10=10, STEP, writecheck = 0;
	clock_t t1, t2; 	  t1 = clock();
	
	double GE, GV1, GV2, GV3, GV4, GV5, GV6, GV7, GV8, GV9, GH;
	double DE ,DV1, DV2, DV3, DV4, DH1, DH2, DH3, Ev_KB_T, dydt, dt, Ev, Ei, reference_density;
	double e, c, v1,v2, v3,v4, h1, h2, h3, totaltime = 0.0, begintime=0.0;
	double vf1, vf2, vf3, vf4, hf1, hf2, hf3, cvb, cfb;
	double te, tv1, tv2, tv3, tv4, th1, th2, th3;
	double kvf1c, kvf2c, kvf3c, kvf4c, kvf1, kvf2, kvf3, kvf4, khf1c, khf2c, khf3c, khf1, khf2, khf3, ke, keb; 
	double per, sinkr1, sinkr2, sinkr3, sinkr4, a0, Omiga, R0, R1, R2, R3, R4, R;
	char filename[100] = "00000000000000000000000000000";
	
	ofstream note;
	ifstream myfile;
	myfile.open("input.txt");  
	myfile >>N_IC;
	myfile >>N_ICh;
	myfile >>N_ICv;	
	myfile >> dt;   	
	myfile >> STEP;	
	myfile >> reference_density;	
	myfile.close();
		
	reference_density = pow(10, -reference_density);

	vector<vector<double> > a, k1, k2, k3, k4, KV1, KV2, KV3, KV4, KH1, KH2, KH3, KE, KEI, KEV, local_lost, gain_e;
	
	GENERATION(GE, GV1, GV2, GV3, GV4, GV5, GV6, GV7, GV8, GV9, GH);  
	DIFFUSION(DE, DV1, DV2, DV3, DV4, DH1, DH2, DH3);
	
	DE = DE * reference_density;
	DV1 = DV1 * reference_density;
	DV2 = DV2 * reference_density;
	DV3 = DV3 * reference_density;
	DV4 = DV4 * reference_density;
	DH1 = DH1 * reference_density;
	DH2 = DH2 * reference_density;
	DH3 = DH3 * reference_density;
	
	GE = GE/reference_density;
	GV1 = GV1/reference_density;
	GV2 = GV2/reference_density;
	GV3 = GV3/reference_density;
	GV4 = GV4/reference_density;
	GV5 = GV5/reference_density;
	GV6 = GV6/reference_density;
	GV7 = GV7/reference_density;
	GV8 = GV8/reference_density;
	GV9 = GV9/reference_density;
	GH = GH/reference_density;
	
	a.resize(N_IC);		

	KV1.resize(N_IC); 		KV2.resize(N_IC); 		KV3.resize(N_IC); 		KV4.resize(N_IC); 
	KH1.resize(N_IC); 		KH2.resize(N_IC); 		KH3.resize(N_IC); 
	KE.resize(N_IC);
	KEI.resize(N_IC);
	KEV.resize(N_IC);
	local_lost.resize(N_IC);
	gain_e.resize(N_IC);

	for(h=0; h<N_IC; h++){
		a[h].resize(N_IC);			
	
		KV1[h].resize(N_IC); 		KV2[h].resize(N_IC); 		KV3[h].resize(N_IC);		KV4[h].resize(N_IC); 
		KH1[h].resize(N_IC); 		KH2[h].resize(N_IC); 		KH3[h].resize(N_IC); 
		KE[h].resize(N_IC);	
		KEI[h].resize(N_IC);
		KEV[h].resize(N_IC);
		local_lost[h].resize(N_IC);
		gain_e[h].resize(N_IC);	
		
	}
	
	a0 = 0.287;  			
	Omiga=a0*a0*a0/2.0;				
	R0 = 0.62;
	
	Ev_KB_T = (16202.1/(1.38*698.15));  // transfer unit to eV
	
	R1 = 0.287*cbrt(3.0/(8*3.1415926));
	R2 = 0.287*cbrt(6.0/(8*3.1415926));
	R3 = 0.287*cbrt(9.0/(8*3.1415926));
	R4 = 0.287*cbrt(12.0/(8*3.1415926));	
	
	for(h=0; h<N_IC; h++){
		for(v=0; v<N_IC; v++){
			R = 0.287*cbrt(3.0*(h+v)/(8*3.1415926));
			KV1[h][v] = DV1*(R1 + R + R0);
			KV2[h][v] = DV2*(R2 + R + R0);
			KV3[h][v] = DV3*(R3 + R + R0);
			KV4[h][v] = DV4*(R4 + R + R0);
			KH1[h][v] = DH1*(R1 + R + R0);
			KH2[h][v] = DH2*(R2 + R + R0);
			KH3[h][v] = DH3*(R3 + R + R0);
			
			KE[h][v] = DE*(R1 + R + R0);
			if(v==0 ) KE[h][v] =0.0;
		}
	}

// reactio with E
	KE[0][0] = 0.0;
	KV1[0][0] = (DV1+ DE)*(R1 + R1 + R0);
	KV2[0][0] = (DV2+ DE)*(R2 + R1 + R0);
	KV3[0][0] = (DV3+ DE)*(R3 + R1 + R0);
	KV4[0][0] = (DV4+ DE)*(R4 + R1 + R0);
	KH1[0][0] = 0.0;
	KH2[0][0] = 0.0;
	KH3[0][0] = 0.0;

// react with V1
	KE[0][1]  = (DV1+ DE)*(R1 + R1 + R0);
	KV1[0][1] = (DV1+ DV1)*(R1 + R1 + R0);
	KV2[0][1] = (DV1+ DV2)*(R2 + R1 + R0);
	KV3[0][1] = (DV1+ DV3)*(R3 + R1 + R0);
	KV4[0][1] = (DV1+ DV4)*(R4 + R1 + R0);
	KH1[0][1] = (DV1+ DH1)*(R1 + R1 + R0);
	KH2[0][1] = (DV1+ DH2)*(R2 + R1 + R0);
	KH3[0][1] = (DV1+ DH3)*(R3 + R1 + R0);
	
// react with V2
	KE[0][2]  = (DV2+ DE)*(R2 + R1 + R0);	
	KV1[0][2] = (DV2+ DV1)*(R1 + R2 + R0);
	KV2[0][2] = (DV2+ DV2)*(R2 + R2 + R0);
	KV3[0][2] = (DV2+ DV3)*(R3 + R2 + R0);
	KV4[0][2] = (DV2+ DV4)*(R4 + R2 + R0);
	KH1[0][2] = (DV2+ DH1)*(R1 + R2 + R0);
	KH2[0][2] = (DV2+ DH2)*(R2 + R2 + R0);
	KH3[0][2] = (DV2+ DH3)*(R3 + R2 + R0);

// react with V3	
	KE[0][3]  = (DV3+ DE)*(R3 + R1 + R0);
	KV1[0][3] = (DV3 + DV1)*(R1 + R3 + R0);
	KV2[0][3] = (DV3 + DV2)*(R2 + R3 + R0);
	KV3[0][3] = (DV3 + DV3)*(R3 + R3 + R0);
	KV4[0][3] = (DV3 + DV4)*(R4 + R3 + R0);
	KH1[0][3] = (DV3 + DH1)*(R1 + R3 + R0);
	KH2[0][3] = (DV3 + DH2)*(R2 + R3 + R0);
	KH3[0][3] = (DV3 + DH3)*(R3 + R3 + R0);
	
// react with V4
	KE[0][4]  = (DV4+ DE)*(R4 + R1 + R0);	
	KV1[0][4] = (DV4 + DV1)*(R1 + R4 + R0);
	KV2[0][4] = (DV4 + DV2)*(R2 + R4 + R0);
	KV3[0][4] = (DV4 + DV3)*(R3 + R4 + R0);
	KV4[0][4] = (DV4 + DV4)*(R4 + R4 + R0);
	KH1[0][4] = (DV4 + DH1)*(R1 + R4 + R0);
	KH2[0][4] = (DV4 + DH2)*(R2 + R4 + R0);
	KH3[0][4] = (DV4 + DH3)*(R3 + R4 + R0);

// react with H1	
	KV1[1][0] = (DH1+ DV1)*(R1 + R1 + R0);
	KV2[1][0] = (DH1+ DV2)*(R2 + R1 + R0);
	KV3[1][0] = (DH1+ DV3)*(R3 + R1 + R0);
	KV4[1][0] = (DH1+ DV4)*(R4 + R1 + R0);
	KH1[1][0] = (DH1+ DH1)*(R1 + R1 + R0);
	KH2[1][0] = (DH1+ DH2)*(R2 + R1 + R0);
	KH3[1][0] = (DH1+ DH3)*(R3 + R1 + R0);
	
// react with H2	
	KV1[2][0] = (DH2+ DV1)*(R1 + R2 + R0);
	KV2[2][0] = (DH2+ DV2)*(R2 + R2 + R0);
	KV3[2][0] = (DH2+ DV3)*(R3 + R2 + R0);
	KV4[2][0] = (DH2+ DV4)*(R4 + R2 + R0);
	KH1[2][0] = (DH2+ DH1)*(R1 + R2 + R0);
	KH2[2][0] = (DH2+ DH2)*(R2 + R2 + R0);
	KH3[2][0] = (DH2+ DH3)*(R3 + R2 + R0);
		
// react with H3	
	KV1[3][0] = (DH3 + DV1)*(R1 + R3 + R0);
	KV2[3][0] = (DH3 + DV2)*(R2 + R3 + R0);
	KV3[3][0] = (DH3 + DV3)*(R3 + R3 + R0);
	KV4[3][0] = (DH3 + DV4)*(R4 + R3 + R0);
	KH1[3][0] = (DH3 + DH1)*(R1 + R3 + R0);
	KH2[3][0] = (DH3 + DH2)*(R2 + R3 + R0);
	KH3[3][0] = (DH3 + DH3)*(R3 + R3 + R0);	

	a[0][0] = a[0][0] + GE *dt; 		 			 
	a[0][1] = a[0][1] + GV1*dt;		 			
	a[0][2] = a[0][2] + GV2*dt ;		 			
	a[0][3] = a[0][3] + GV3*dt; 			
	a[0][4] = a[0][4] + GV4*dt;			
	a[0][5] = a[0][5] + GV5*dt;
	a[0][6] = a[0][6] + GV6*dt;
	a[0][7] = a[0][7] + GV7*dt;
	a[0][8] = a[0][8] + GV8*dt;
	a[0][9] = a[0][9] + GV9*dt;
	a[1][0] = a[1][0] + GH*dt;
	
	sinkr1 = 1.01*(-4*3.1415926*2.5*0.0001/log(3.1415926*2.5*0.0001*(0.62+R1)*(0.62+R1)))/reference_density;
	
	for(h=0; h<N_IC; h++){
		for(v=1; v<N_IC; v++){	
			Ev = 1.73 - 2.59*(pow(v, 0.6666) - pow(v-1, 0.6666)) + 2.50*log(1.0 + 1.0*h/v);
			Ev = exp(-1.0 * Ev * Ev_KB_T)/Omiga;
			if(v == 1 && h ==0 ) Ev = 0.0;	
			KEV[h][v] =  KV1[h][v-1] * Ev;
		}
	}

	for(h=0; h<N_IC; h++){
		for(v=0; v<N_IC-1; v++){
			Ei = 4.88 + 2.59*(pow(v, 0.6666) - pow(v-1, 0.6666)) - 2.50*log(1.0 + 1.0*h/v);
			Ei = exp(-1.0 * Ei * Ev_KB_T)/Omiga;	
			if(v == 0 ) Ei = 0.0;			
			KEI[h][v] =  KE[h][v+1] * Ei;
		}
	}	
	
	/*
	ofstream kvalue;	
	kvalue.open("reaction.txt");
	
	kvalue<<"step "<<dt<<endl;
	kvalue<<"reference density  "<<reference_density<<endl;
	kvalue<<"generation "<<endl;
	kvalue<< GE<<" "<<GV1<<" "<<GV2<<" "<<GV3<<" "<<GV4<<" "<<GV5<<" "<<GV6<<" "<<GV7<<" "<<GV8<<" "<<GV9<<"  "<<GH<<endl;
	
	kvalue<<"sink strength for v1 and h1"<<endl;
	kvalue<<"for E  "<< DE *sinkr1<<endl;	
	kvalue<<"for v1  "<<DV1*sinkr1<<endl;
	kvalue<<"for h1  "<<DH1*sinkr1<<endl;
	
	kvalue<<" e  "<<endl;
	for( h =0; h < first10; h++){ 	for(v =0; v< first10; v++) kvalue<<KE[h][v]<<"  "; kvalue<<endl; } kvalue<<endl;
	
	kvalue<<" kv1  "<<endl;
	for( h =0; h < first10; h++){ 	for(v =0; v< first10; v++) kvalue<<KV1[h][v]<<"  "; kvalue<<endl;	}kvalue<<endl;
	
	kvalue<<" kv2  "<<endl;
	for( h =0; h < first10; h++){ 	for(v =0; v< first10; v++) kvalue<<KV2[h][v]<<"  "; kvalue<<endl;  }kvalue<<endl;
	
	kvalue<<" kv3  "<<endl;
	for( h =0; h < first10; h++){	for(v =0; v< first10; v++) kvalue<<KV3[h][v]<<"  "; kvalue<<endl;	}kvalue<<endl;
	
	kvalue<<" kv4  "<<endl;
	for( h =0; h < first10; h++){ 	for(v =0; v< first10; v++) kvalue<<KV4[h][v]<<"  "; kvalue<<endl;	}kvalue<<endl;
	
	kvalue<<" kh1  "<<endl;
	for( h =0; h < first10; h++){ 	for(v =0; v< first10; v++) kvalue<<KH1[h][v]<<"  "; kvalue<<endl;	}kvalue<<endl;
	
	kvalue<<" kh2  "<<endl;
	for( h =0; h < first10; h++){ 	for(v =0; v< first10; v++) kvalue<<KH2[h][v]<<"  "; kvalue<<endl;	}kvalue<<endl;
	
	kvalue<<" kh3  "<<endl;
	for( h =0; h < first10; h++){ 	for(v =0; v< first10; v++) kvalue<<KH3[h][v]<<"  "; kvalue<<endl;	}kvalue<<endl;
	
	kvalue.close();
	
	ofstream pkvf1c, pkvf2c, pkvf3c, pkvf4c, pkvf1, pkvf2, pkvf3, pkvf4, pkhf1c, pkhf2c, pkhf3c, pkhf1, pkhf2, pkhf3, pke, pkeb;
	pkvf1c.open("pkvf1c.txt");
	pkvf2c.open("pkvf2c.txt");
	pkvf3c.open("pkvf3c.txt");
	pkvf4c.open("pkvf4c.txt");	
	pkvf1.open("pkvf1.txt");
	pkvf2.open("pkvf2.txt");
	pkvf3.open("pkvf3.txt");
	pkvf4.open("pkvf4.txt");		
	pkhf1c.open("pkhf1c.txt");
	pkhf2c.open("pkhf2c.txt");
	pkhf3c.open("pkhf3c.txt");
	pkhf1.open("pkhf1.txt");
	pkhf2.open("pkhf2.txt");
	pkhf3.open("pkhf3.txt");
	pke.open("pke.txt");
	pkeb.open("pkeb.txt");	
	*/
	
	ifstream inputdata;
	inputdata.open("olddata.txt"); 
	inputdata>>begintime;
	for(h =0; h<N_ICh; h++){	
		for(v =0; v< N_ICv ; v++)	inputdata>>a[h][v];
	}
	inputdata.close();

	for(h =0; h<N_ICh; h++){	
		for(v =0; v< N_ICv ; v++)	cout<< a[h][v]<<"  ";  	cout<<endl;
	}	
	

	
	sprintf(filename,"date begin time %.3f target time %.3f dt %f.txt", begintime, 1.0*begintime + 1.0*STEP*dt, dt);		
	ofstream indata; 		
	indata.open(filename);
	
	do{ 
		st++; 	
		k1.clear();				k2.clear();				k3.clear();				k4.clear();	
		k1.resize(N_IC); 		k2.resize(N_IC);		k3.resize(N_IC); 		k4.resize(N_IC);
		for(h=0; h<N_IC; h++){
			k1[h].resize(N_IC); 	k2[h].resize(N_IC);		k3[h].resize(N_IC); 	k4[h].resize(N_IC); 
		}
		
// the k1 term 
{
		te = 0.0; 		
		tv1 = 0.0; 		 tv2 = 0.0;  	tv3 = 0.0;  	tv4 = 0.0;		
		th1 = 0.0; 		 th2 = 0.0; 	th3 = 0.0;
	
		e  = a[0][0];
		v1 = a[0][1];  	v2 = a[0][2];  	v3 = a[0][3];  		v4 = a[0][4];   
		h1 = a[1][0];	h2 = a[2][0];	h3 = a[3][0];
		
		omp_set_num_threads(2);
		#pragma omp parallel for private(c) reduction(+: te, tv1,tv2, tv3, tv4, th1, th2, th3)	
		for(int oh =0; oh<N_ICh; oh++){	
			for(int ov =0; ov< N_ICv ; ov++){
				c = a[oh][ov];
				k1[oh][ov] = - KV1[oh][ov]*c*v1 - KV2[oh][ov]*c*v2 -KV3[oh][ov]*c*v3 - KV4[oh][ov]*c*v4
							 - KH1[oh][ov]*c*h1 - KH2[oh][ov]*c*h2 -KH3[oh][ov]*c*h3 - KE[oh][ov]*c*e;
						 
						 
				if(ov!=N_IC-1){ k1[oh][ov] = k1[oh][ov] + KE[oh][ov+1] * a[oh][ov+1] * e; }
				
				if(oh > 0 ){
					if(ov>=1){ k1[oh][ov] = k1[oh][ov] + KV1[oh][ov-1] * a[oh][ov-1] * v1; }
					if(ov>=2){ k1[oh][ov] = k1[oh][ov] + KV2[oh][ov-2] * a[oh][ov-2] * v2; }
					if(ov>=3){ k1[oh][ov] = k1[oh][ov] + KV3[oh][ov-3] * a[oh][ov-3] * v3; }	
					if(ov>=4){ k1[oh][ov] = k1[oh][ov] + KV4[oh][ov-4] * a[oh][ov-4] * v4; }
				}
				else{
					if(ov>1){ k1[oh][ov] = k1[oh][ov] + KV1[oh][ov-1] * a[oh][ov-1] * v1; }
					if(ov>2){ k1[oh][ov] = k1[oh][ov] + KV2[oh][ov-2] * a[oh][ov-2] * v2; }
					if(ov>3){ k1[oh][ov] = k1[oh][ov] + KV3[oh][ov-3] * a[oh][ov-3] * v3; }	
					if(ov>4){ k1[oh][ov] = k1[oh][ov] + KV4[oh][ov-4] * a[oh][ov-4] * v4; }
				}
				
				if(oh>=1){ k1[oh][ov] = k1[oh][ov] + KH1[oh-1][ov] * a[oh-1][ov] * h1; }
				if(oh>=2){ k1[oh][ov] = k1[oh][ov] + KH2[oh-2][ov] * a[oh-2][ov] * h2; }	
				if(oh>=3){ k1[oh][ov] = k1[oh][ov] + KH3[oh-3][ov] * a[oh-3][ov] * h3; }
				
				te += -1.0*KE[oh][ov]*c*e;
				tv1 += -1.0*KV1[oh][ov]*c*v1; 						
				tv2 += -1.0*KV2[oh][ov]*c*v2;
				tv3 += -1.0*KV3[oh][ov]*c*v3;						
				tv4 += -1.0*KV4[oh][ov]*c*v4;
				th1 += -1.0*KH1[oh][ov]*c*h1;						
				th2 += -1.0*KH2[oh][ov]*c*h2;
				th3 += -1.0*KH3[oh][ov]*c*h3;
			}
		}
		
		tv1= tv1 + KV2[0][1]*v1*v2 + KV3[0][1]*v1*v3 + KV4[0][1]*v1*v4
				 + KH1[0][1]*v1*h1 + KH2[0][1]*v1*h2 + KH3[0][1]*v1*h3 + KE[0][1]*v1*e; 
				 
		tv2= tv2 + KV1[0][2]*v2*v1 + KV3[0][2]*v2*v3 + KV4[0][2]*v2*v4
				 + KH1[0][2]*v2*h1 + KH2[0][2]*v2*h2 + KH3[0][2]*v2*h3 + KE[0][2]*v2*e;
				 
		tv3= tv3 + KV1[0][3]*v3*v1 + KV2[0][3]*v3*v2 + KV4[0][3]*v3*v4
				 + KH1[0][3]*v3*h1 + KH2[0][3]*v3*h2 + KH3[0][3]*v3*h3 + KE[0][3]*v3*e;	
				 
		tv4= tv4 + KV1[0][4]*v4*v1 + KV2[0][4]*v4*v2 + KV3[0][4]*v4*v3 
				 + KH1[0][4]*v4*h1 + KH2[0][4]*v4*h2 + KH3[0][4]*v4*h3 + KE[0][4]*v4*e;
				
		th1= th1 + KV1[1][0]*h1*v1 + KV2[1][0]*h1*v2 + KV3[1][0]*h1*v3 + KV4[1][0]*h1*v4
				 + KH2[1][0]*h1*h2 + KH3[1][0]*h1*h3 + KE[1][0]*h1*e;	
				 
		th2= th2 + KV1[2][0]*h2*v1 + KV2[2][0]*h2*v2 + KV3[2][0]*h2*v3 + KV4[2][0]*h2*v4
				 + KH1[2][0]*h2*h1 + KH3[2][0]*h2*h3 + KE[2][0]*h2*e;
				 
		th3= th3 + KV1[3][0]*h3*v1 + KV2[3][0]*h3*v2 + KV3[3][0]*h3*v3 + KV4[3][0]*h3*v4
				 + KH1[3][0]*h3*h1 + KH2[3][0]*h3*h2  + KE[3][0]*h3*e;
		
		k1[0][0] = 			  GE - DE*sinkr1 * e  + te; 
		
		k1[0][1] = k1[0][1] + GV1 + tv1;		 			
		k1[0][2] = k1[0][2] + GV2 + tv2;			 			
		k1[0][3] = k1[0][3] + GV3 + tv3;	 			
		k1[0][4] = k1[0][4] + GV4 + tv4;				
		k1[0][5] = k1[0][5] + GV5;
		k1[0][6] = k1[0][6] + GV6;
		k1[0][7] = k1[0][7] + GV7;
		k1[0][8] = k1[0][8] + GV8;
		k1[0][9] = k1[0][9] + GV9;
		k1[1][0] = k1[1][0] + GH + th1;   
		k1[2][0] = k1[2][0] + th2;  
		k1[3][0] = k1[3][0] + th3; 		
}	
// the k2 term 
{	
		te = 0.0; 		tv1 = 0.0; 		 tv2 = 0.0;  	tv3 = 0.0;  	tv4 = 0.0;		th1 = 0.0; 		 th2 = 0.0; 	th3 = 0.0;
		
		e  = a[0][0]+ 0.5 * k1[0][0]*dt;
		v1 = a[0][1]+ 0.5 * k1[0][1]*dt;  	
		v2 = a[0][2]+ 0.5 * k1[0][2]*dt;  	
		v3 = a[0][3]+ 0.5 * k1[0][3]*dt;  		
		v4 = a[0][4]+ 0.5 * k1[0][4]*dt;   
		h1 = a[1][0]+ 0.5 * k1[1][0]*dt;	
		h2 = a[2][0]+ 0.5 * k1[2][0]*dt;	
		h3 = a[3][0]+ 0.5 * k1[3][0]*dt;
		
		omp_set_num_threads(2);		
		#pragma omp parallel for private(c) reduction(+: te, tv1,tv2, tv3, tv4, th1, th2, th3)
		for(int oh =0; oh<N_ICh; oh++){	
			for(int ov =0; ov< N_ICv ; ov++){
				c = a[oh][ov] + 0.5 * k1[oh][ov]*dt ;
				
				k2[oh][ov] = - KV1[oh][ov]*c*v1 - KV2[oh][ov]*c*v2 -KV3[oh][ov]*c*v3 - KV4[oh][ov]*c*v4
							 - KH1[oh][ov]*c*h1 - KH2[oh][ov]*c*h2 -KH3[oh][ov]*c*h3 - KE[oh][ov]*c*e;
						 
				
				if(ov!=N_IC-1){ k2[oh][ov] = k2[oh][ov] + KE[oh][ov+1]*(a[oh][ov+1]+0.5*k1[oh][ov+1]*dt)*e; }
				
				if(oh > 0 ){
					if(ov>=1){ k2[oh][ov] = k2[oh][ov] + KV1[oh][ov-1]*(a[oh][ov-1]+0.5*k1[oh][ov-1]*dt)*v1; }
					if(ov>=2){ k2[oh][ov] = k2[oh][ov] + KV2[oh][ov-2]*(a[oh][ov-2]+0.5*k1[oh][ov-2]*dt)*v2; }
					if(ov>=3){ k2[oh][ov] = k2[oh][ov] + KV3[oh][ov-3]*(a[oh][ov-3]+0.5*k1[oh][ov-3]*dt)*v3; }	
					if(ov>=4){ k2[oh][ov] = k2[oh][ov] + KV4[oh][ov-4]*(a[oh][ov-4]+0.5*k1[oh][ov-4]*dt)*v4; }
				}
				else{
					if(ov>1){ k2[oh][ov] = k2[oh][ov] + KV1[oh][ov-1]*(a[oh][ov-1]+0.5*k1[oh][ov-1]*dt)*v1; }
					if(ov>2){ k2[oh][ov] = k2[oh][ov] + KV2[oh][ov-2]*(a[oh][ov-2]+0.5*k1[oh][ov-2]*dt)*v2; }
					if(ov>3){ k2[oh][ov] = k2[oh][ov] + KV3[oh][ov-3]*(a[oh][ov-3]+0.5*k1[oh][ov-3]*dt)*v3; }	
					if(ov>4){ k2[oh][ov] = k2[oh][ov] + KV4[oh][ov-4]*(a[oh][ov-4]+0.5*k1[oh][ov-4]*dt)*v4; }
				}
				
				if(oh>=1){ k2[oh][ov] = k2[oh][ov] + KH1[oh-1][ov]*(a[oh-1][ov]+0.5*k1[oh-1][ov]*dt)*h1; }
				if(oh>=2){ k2[oh][ov] = k2[oh][ov] + KH2[oh-2][ov]*(a[oh-2][ov]+0.5*k1[oh-2][ov]*dt)*h2; }	
				if(oh>=3){ k2[oh][ov] = k2[oh][ov] + KH3[oh-3][ov]*(a[oh-3][ov]+0.5*k1[oh-3][ov]*dt)*h3; }
				
				te += -1.0*KE[oh][ov]*c*e;
				tv1 += -1.0*KV1[oh][ov]*c*v1; 						
				tv2 += -1.0*KV2[oh][ov]*c*v2;
				tv3 += -1.0*KV3[oh][ov]*c*v3;						
				tv4 += -1.0*KV4[oh][ov]*c*v4;
				th1 += -1.0*KH1[oh][ov]*c*h1;						
				th2 += -1.0*KH2[oh][ov]*c*h2;
				th3 += -1.0*KH3[oh][ov]*c*h3;
			}
		}
		
		tv1= tv1 + KV2[0][1]*v1*v2 + KV3[0][1]*v1*v3 + KV4[0][1]*v1*v4
				 + KH1[0][1]*v1*h1 + KH2[0][1]*v1*h2 + KH3[0][1]*v1*h3 + KE[0][1]*v1*e; 
				 
		tv2= tv2 + KV1[0][2]*v2*v1 + KV3[0][2]*v2*v3 + KV4[0][2]*v2*v4
				 + KH1[0][2]*v2*h1 + KH2[0][2]*v2*h2 + KH3[0][2]*v2*h3 + KE[0][2]*v2*e;
				 
		tv3= tv3 + KV1[0][3]*v3*v1 + KV2[0][3]*v3*v2 + KV4[0][3]*v3*v4
				 + KH1[0][3]*v3*h1 + KH2[0][3]*v3*h2 + KH3[0][3]*v3*h3 + KE[0][3]*v3*e;	
				 
		tv4= tv4 + KV1[0][4]*v4*v1 + KV2[0][4]*v4*v2 + KV3[0][4]*v4*v3 
				 + KH1[0][4]*v4*h1 + KH2[0][4]*v4*h2 + KH3[0][4]*v4*h3 + KE[0][4]*v4*e;
				
		th1= th1 + KV1[1][0]*h1*v1 + KV2[1][0]*h1*v2 + KV3[1][0]*h1*v3 + KV4[1][0]*h1*v4
				 + KH2[1][0]*h1*h2 + KH3[1][0]*h1*h3 + KE[1][0]*h1*e;	
				 
		th2= th2 + KV1[2][0]*h2*v1 + KV2[2][0]*h2*v2 + KV3[2][0]*h2*v3 + KV4[2][0]*h2*v4
				 + KH1[2][0]*h2*h1 + KH3[2][0]*h2*h3 + KE[2][0]*h2*e;
				 
		th3= th3 + KV1[3][0]*h3*v1 + KV2[3][0]*h3*v2 + KV3[3][0]*h3*v3 + KV4[3][0]*h3*v4
				 + KH1[3][0]*h3*h1 + KH2[3][0]*h3*h2  + KE[3][0]*h3*e;
		
		k2[0][0] = 			  GE - DE*sinkr1 *  e  + te; 		 			 
		k2[0][1] = k2[0][1] + GV1 + tv1;		 			
		k2[0][2] = k2[0][2] + GV2 + tv2;			 			
		k2[0][3] = k2[0][3] + GV3 + tv3;	 			
		k2[0][4] = k2[0][4] + GV4 + tv4;				
		k2[0][5] = k2[0][5] + GV5;
		k2[0][6] = k2[0][6] + GV6;
		k2[0][7] = k2[0][7] + GV7;
		k2[0][8] = k2[0][8] + GV8;
		k2[0][9] = k2[0][9] + GV9;
		k2[1][0] = k2[1][0] + GH + th1;   
		k2[2][0] = k2[2][0] + th2;  
		k2[3][0] = k2[3][0] + th3; 		
}	
// the k3 term 		
{	
		te = 0.0; 		tv1 = 0.0; 		 tv2 = 0.0;  	tv3 = 0.0;  	tv4 = 0.0;		th1 = 0.0; 		 th2 = 0.0; 	th3 = 0.0;
		
		e  = a[0][0]+ 0.5 * k2[0][0]*dt;
		v1 = a[0][1]+ 0.5 * k2[0][1]*dt;  	
		v2 = a[0][2]+ 0.5 * k2[0][2]*dt;  	
		v3 = a[0][3]+ 0.5 * k2[0][3]*dt;  		
		v4 = a[0][4]+ 0.5 * k2[0][4]*dt;   
		h1 = a[1][0]+ 0.5 * k2[1][0]*dt;	
		h2 = a[2][0]+ 0.5 * k2[2][0]*dt;	
		h3 = a[3][0]+ 0.5 * k2[3][0]*dt;
		
		omp_set_num_threads(2);		
		#pragma omp parallel for private(c) reduction(+: te, tv1,tv2, tv3, tv4, th1, th2, th3)
		for(int oh =0; oh<N_ICh; oh++){	
			for(int ov =0; ov< N_ICv ; ov++){
				c = a[oh][ov] + 0.5 * k2[oh][ov]*dt ;
				
				k3[oh][ov] = - KV1[oh][ov]*c*v1 - KV2[oh][ov]*c*v2 -KV3[oh][ov]*c*v3 - KV4[oh][ov]*c*v4
							 - KH1[oh][ov]*c*h1 - KH2[oh][ov]*c*h2 -KH3[oh][ov]*c*h3 - KE[oh][ov]*c*e;
						 
				
				if(ov!=N_IC-1){ k3[oh][ov] = k3[oh][ov] + KE[oh][ov+1]*(a[oh][ov+1]+0.5*k2[oh][ov+1]*dt)*e; }
				
				if(oh > 0 ){
					if(ov>=1){ k3[oh][ov] = k3[oh][ov] + KV1[oh][ov-1]*(a[oh][ov-1]+0.5*k2[oh][ov-1]*dt)*v1; }
					if(ov>=2){ k3[oh][ov] = k3[oh][ov] + KV2[oh][ov-2]*(a[oh][ov-2]+0.5*k2[oh][ov-2]*dt)*v2; }
					if(ov>=3){ k3[oh][ov] = k3[oh][ov] + KV3[oh][ov-3]*(a[oh][ov-3]+0.5*k2[oh][ov-3]*dt)*v3; }	
					if(ov>=4){ k3[oh][ov] = k3[oh][ov] + KV4[oh][ov-4]*(a[oh][ov-4]+0.5*k2[oh][ov-4]*dt)*v4; }
				}
				else{
					if(ov>1){ k3[oh][ov] = k3[oh][ov] + KV1[oh][ov-1]*(a[oh][ov-1]+0.5*k2[oh][ov-1]*dt)*v1; }
					if(ov>2){ k3[oh][ov] = k3[oh][ov] + KV2[oh][ov-2]*(a[oh][ov-2]+0.5*k2[oh][ov-2]*dt)*v2; }
					if(ov>3){ k3[oh][ov] = k3[oh][ov] + KV3[oh][ov-3]*(a[oh][ov-3]+0.5*k2[oh][ov-3]*dt)*v3; }	
					if(ov>4){ k3[oh][ov] = k3[oh][ov] + KV4[oh][ov-4]*(a[oh][ov-4]+0.5*k2[oh][ov-4]*dt)*v4; }
				}
				if(oh>=1){ k3[oh][ov] = k3[oh][ov] + KH1[oh-1][ov]*(a[oh-1][ov]+0.5*k2[oh-1][ov]*dt)*h1; }
				if(oh>=2){ k3[oh][ov] = k3[oh][ov] + KH2[oh-2][ov]*(a[oh-2][ov]+0.5*k2[oh-2][ov]*dt)*h2; }	
				if(oh>=3){ k3[oh][ov] = k3[oh][ov] + KH3[oh-3][ov]*(a[oh-3][ov]+0.5*k2[oh-3][ov]*dt)*h3; }
				
				te += -1.0*KE[oh][ov]*c*e;
				tv1 += -1.0*KV1[oh][ov]*c*v1; 						
				tv2 += -1.0*KV2[oh][ov]*c*v2;
				tv3 += -1.0*KV3[oh][ov]*c*v3;						
				tv4 += -1.0*KV4[oh][ov]*c*v4;
				th1 += -1.0*KH1[oh][ov]*c*h1;						
				th2 += -1.0*KH2[oh][ov]*c*h2;
				th3 += -1.0*KH3[oh][ov]*c*h3;
			}
		}
		
		tv1= tv1 + KV2[0][1]*v1*v2 + KV3[0][1]*v1*v3 + KV4[0][1]*v1*v4
				 + KH1[0][1]*v1*h1 + KH2[0][1]*v1*h2 + KH3[0][1]*v1*h3 + KE[0][1]*v1*e; 
				 
		tv2= tv2 + KV1[0][2]*v2*v1 + KV3[0][2]*v2*v3 + KV4[0][2]*v2*v4
				 + KH1[0][2]*v2*h1 + KH2[0][2]*v2*h2 + KH3[0][2]*v2*h3 + KE[0][2]*v2*e;
				 
		tv3= tv3 + KV1[0][3]*v3*v1 + KV2[0][3]*v3*v2 + KV4[0][3]*v3*v4
				 + KH1[0][3]*v3*h1 + KH2[0][3]*v3*h2 + KH3[0][3]*v3*h3 + KE[0][3]*v3*e;	
				 
		tv4= tv4 + KV1[0][4]*v4*v1 + KV2[0][4]*v4*v2 + KV3[0][4]*v4*v3 
				 + KH1[0][4]*v4*h1 + KH2[0][4]*v4*h2 + KH3[0][4]*v4*h3 + KE[0][4]*v4*e;
				
		th1= th1 + KV1[1][0]*h1*v1 + KV2[1][0]*h1*v2 + KV3[1][0]*h1*v3 + KV4[1][0]*h1*v4
				 + KH2[1][0]*h1*h2 + KH3[1][0]*h1*h3 + KE[1][0]*h1*e;	
				 
		th2= th2 + KV1[2][0]*h2*v1 + KV2[2][0]*h2*v2 + KV3[2][0]*h2*v3 + KV4[2][0]*h2*v4
				 + KH1[2][0]*h2*h1 + KH3[2][0]*h2*h3 + KE[2][0]*h2*e;
				 
		th3= th3 + KV1[3][0]*h3*v1 + KV2[3][0]*h3*v2 + KV3[3][0]*h3*v3 + KV4[3][0]*h3*v4
				 + KH1[3][0]*h3*h1 + KH2[3][0]*h3*h2  + KE[3][0]*h3*e;
		
		k3[0][0] = 			  GE - DE*sinkr1 *  e  + te; 		 			 
		k3[0][1] = k3[0][1] + GV1 + tv1;		 			
		k3[0][2] = k3[0][2] + GV2 + tv2;			 			
		k3[0][3] = k3[0][3] + GV3 + tv3;	 			
		k3[0][4] = k3[0][4] + GV4 + tv4;				
		k3[0][5] = k3[0][5] + GV5;
		k3[0][6] = k3[0][6] + GV6;
		k3[0][7] = k3[0][7] + GV7;
		k3[0][8] = k3[0][8] + GV8;
		k3[0][9] = k3[0][9] + GV9;
		k3[1][0] = k3[1][0] + GH + th1;   
		k3[2][0] = k3[2][0] + th2;  
		k3[3][0] = k3[3][0] + th3; 		
}			
// the k4 term 
{	
		te = 0.0; 		
		tv1 = 0.0; 		 tv2 = 0.0;  	tv3 = 0.0;  	tv4 = 0.0;		
		th1 = 0.0; 		 th2 = 0.0; 	th3 = 0.0;
		
		e  = a[0][0]+ k3[0][0]*dt;
		v1 = a[0][1]+ k3[0][1]*dt;  	
		v2 = a[0][2]+ k3[0][2]*dt;  	
		v3 = a[0][3]+ k3[0][3]*dt;  		
		v4 = a[0][4]+ k3[0][4]*dt;   
		h1 = a[1][0]+ k3[1][0]*dt;	
		h2 = a[2][0]+ k3[2][0]*dt;	
		h3 = a[3][0]+ k3[3][0]*dt;
		omp_set_num_threads(2);		
		#pragma omp parallel for private(c) reduction(+: te, tv1, tv2, tv3, tv4, th1, th2, th3)
		for(int oh =0; oh<N_ICh; oh++){	
			for(int ov =0; ov< N_ICv ; ov++){
				c = a[oh][ov] + k3[oh][ov]*dt ;
				
				k4[oh][ov] = - KV1[oh][ov]*c*v1 - KV2[oh][ov]*c*v2 -KV3[oh][ov]*c*v3 - KV4[oh][ov]*c*v4
							 - KH1[oh][ov]*c*h1 - KH2[oh][ov]*c*h2 -KH3[oh][ov]*c*h3 - KE[oh][ov]*c*e;
						 
				
				if(ov!=N_IC-1){ k4[oh][ov] = k4[oh][ov] + KE[oh][ov+1]*(a[oh][ov+1]+k3[oh][ov+1]*dt)*e; }
				
				if(oh > 0 ){					
					if(ov>=1){ k4[oh][ov] = k4[oh][ov] + KV1[oh][ov-1]*(a[oh][ov-1]+k3[oh][ov-1]*dt)*v1; }
					if(ov>=2){ k4[oh][ov] = k4[oh][ov] + KV2[oh][ov-2]*(a[oh][ov-2]+k3[oh][ov-2]*dt)*v2; }
					if(ov>=3){ k4[oh][ov] = k4[oh][ov] + KV3[oh][ov-3]*(a[oh][ov-3]+k3[oh][ov-3]*dt)*v3; }	
					if(ov>=4){ k4[oh][ov] = k4[oh][ov] + KV4[oh][ov-4]*(a[oh][ov-4]+k3[oh][ov-4]*dt)*v4; }
				}
				else{
					if(ov>1){ k4[oh][ov] = k4[oh][ov] + KV1[oh][ov-1]*(a[oh][ov-1]+k3[oh][ov-1]*dt)*v1; }
					if(ov>2){ k4[oh][ov] = k4[oh][ov] + KV2[oh][ov-2]*(a[oh][ov-2]+k3[oh][ov-2]*dt)*v2; }
					if(ov>3){ k4[oh][ov] = k4[oh][ov] + KV3[oh][ov-3]*(a[oh][ov-3]+k3[oh][ov-3]*dt)*v3; }	
					if(ov>4){ k4[oh][ov] = k4[oh][ov] + KV4[oh][ov-4]*(a[oh][ov-4]+k3[oh][ov-4]*dt)*v4; }	
				}
				
				if(oh>=1){ k4[oh][ov] = k4[oh][ov] + KH1[oh-1][ov]*(a[oh-1][ov]+k3[oh-1][ov]*dt)*h1; }
				if(oh>=2){ k4[oh][ov] = k4[oh][ov] + KH2[oh-2][ov]*(a[oh-2][ov]+k3[oh-2][ov]*dt)*h2; }	
				if(oh>=3){ k4[oh][ov] = k4[oh][ov] + KH3[oh-3][ov]*(a[oh-3][ov]+k3[oh-3][ov]*dt)*h3; }
				
				te += -1.0*KE[oh][ov]*c*e;
				tv1 += -1.0*KV1[oh][ov]*c*v1; 						
				tv2 += -1.0*KV2[oh][ov]*c*v2;
				tv3 += -1.0*KV3[oh][ov]*c*v3;						
				tv4 += -1.0*KV4[oh][ov]*c*v4;
				th1 += -1.0*KH1[oh][ov]*c*h1;						
				th2 += -1.0*KH2[oh][ov]*c*h2;
				th3 += -1.0*KH3[oh][ov]*c*h3;
			}
		}
		
		tv1= tv1 + KV2[0][1]*v1*v2 + KV3[0][1]*v1*v3 + KV4[0][1]*v1*v4
				 + KH1[0][1]*v1*h1 + KH2[0][1]*v1*h2 + KH3[0][1]*v1*h3 + KE[0][1]*v1*e; 
				 
		tv2= tv2 + KV1[0][2]*v2*v1 + KV3[0][2]*v2*v3 + KV4[0][2]*v2*v4
				 + KH1[0][2]*v2*h1 + KH2[0][2]*v2*h2 + KH3[0][2]*v2*h3 + KE[0][2]*v2*e;
				 
		tv3= tv3 + KV1[0][3]*v3*v1 + KV2[0][3]*v3*v2 + KV4[0][3]*v3*v4
				 + KH1[0][3]*v3*h1 + KH2[0][3]*v3*h2 + KH3[0][3]*v3*h3 + KE[0][3]*v3*e;	
				 
		tv4= tv4 + KV1[0][4]*v4*v1 + KV2[0][4]*v4*v2 + KV3[0][4]*v4*v3 
				 + KH1[0][4]*v4*h1 + KH2[0][4]*v4*h2 + KH3[0][4]*v4*h3 + KE[0][4]*v4*e;
				
		th1= th1 + KV1[1][0]*h1*v1 + KV2[1][0]*h1*v2 + KV3[1][0]*h1*v3 + KV4[1][0]*h1*v4
				 + KH2[1][0]*h1*h2 + KH3[1][0]*h1*h3 + KE[1][0]*h1*e;	
				 
		th2= th2 + KV1[2][0]*h2*v1 + KV2[2][0]*h2*v2 + KV3[2][0]*h2*v3 + KV4[2][0]*h2*v4
				 + KH1[2][0]*h2*h1 + KH3[2][0]*h2*h3 + KE[2][0]*h2*e;
				 
		th3= th3 + KV1[3][0]*h3*v1 + KV2[3][0]*h3*v2 + KV3[3][0]*h3*v3 + KV4[3][0]*h3*v4
				 + KH1[3][0]*h3*h1 + KH2[3][0]*h3*h2  + KE[3][0]*h3*e;
		
		k4[0][0] = 			  GE - DE*sinkr1 *  e  + te; 		 			 
		k4[0][1] = k4[0][1] + GV1 + tv1;		 			
		k4[0][2] = k4[0][2] + GV2 + tv2;			 			
		k4[0][3] = k4[0][3] + GV3 + tv3;	 			
		k4[0][4] = k4[0][4] + GV4 + tv4;				
		k4[0][5] = k4[0][5] + GV5;
		k4[0][6] = k4[0][6] + GV6;
		k4[0][7] = k4[0][7] + GV7;
		k4[0][8] = k4[0][8] + GV8;
		k4[0][9] = k4[0][9] + GV9;
		k4[1][0] = k4[1][0] + GH + th1;   
		k4[2][0] = k4[2][0] + th2;  
		k4[3][0] = k4[3][0] + th3; 			
}
// update the value 
		
		omp_set_num_threads(2);
		#pragma omp parallel for
		for(int mh =0; mh<N_ICh; mh++){	
			for(int mv =0; mv< N_ICv ; mv++){		
				a[mh][mv] = a[mh][mv]+ (k1[mh][mv] + 2.0*k2[mh][mv] +2.0*k3[mh][mv] +k4[mh][mv])*dt/6.0;
				if( a[mh][mv] < 0.0 ) writecheck = 1;
			}
		}
		
		omp_set_num_threads(2);
		#pragma omp parallel for
		for(int mh =0; mh<N_ICh; mh++){	
			for(int mv =0; mv< N_ICv ; mv++){
				local_lost[mh][mv] = a[mh][mv] -   a[mh][mv] * exp(-1.0 * KEI[mh][mv]*dt); 
				a[mh][mv] = a[mh][mv] * exp(-1.0 * KEI[mh][mv]*dt);
			}
		}
		
		omp_set_num_threads(2);
		#pragma omp parallel for
		for(int mh =0; mh<N_ICh; mh++){	
			for(int mv =1; mv< N_ICv ; mv++){		
				a[mh][mv] = a[mh][mv] + local_lost[mh][mv-1];
				a[0][0] = a[0][0] + local_lost[mh][mv-1];				
			}
		}
		
		if(writecheck == 1){
				totaltime = 1.0* st * dt;			
				st = STEP + 1.0;
		}
		else{
			totaltime = 1.0* st * dt;
			if(st%(STEP/100) == 0){
				
					note.open("olddata.txt");
					note<<1.0*begintime+totaltime<<endl;
					for(h =0; h<N_ICh; h++){	
						for(v =0; v< N_ICv ; v++)		note<<a[h][v]<<endl;
					}
					note.close();
				
				indata<<"runing step  " <<st<< "  the dt is  "<< dt <<" ==> total time is  " << 1.0*st*dt + 1.0*begintime
					  <<" included with the begin time ( "<< 1.0*begintime <<" ) "<<endl;
				for( h = 0; h < N_IC; h++){ 	for(v =0; v< N_IC; v++) indata<<a[h][v]<<"  "; indata<<endl; } indata<<endl; 
			}
		}

	}while(st < STEP);
	
	indata.close();

	t2= clock();
	float diff_time((float)t2 - (float)t1);   
	float spend_time = diff_time/ CLOCKS_PER_SEC;

	sprintf(filename,"RK running time from %.3f to %.3f, using cpu time %.3f.txt", begintime, 1.0*begintime + totaltime, spend_time);	
	ofstream summary;
	summary.open(filename);
	summary<<"runign step  "<<st<<" dt: " <<dt<< " ==> total time: "<< 1.0*begintime + totaltime
			<<" included with begin time ( " << 1.0*begintime<< " ) "<<endl;
	for( h =0; h < N_IC; h++){ 	for(v =0; v< N_IC; v++) summary<<a[h][v]<<"  "; summary<<endl; 		}summary<<endl;
	summary.close();
	
		
	summary.open("olddata.txt");
	summary<<1.0*begintime+totaltime<<endl;
	for(h =0; h<N_ICh; h++){	
		for(v =0; v< N_ICv ; v++)		summary<<a[h][v]<<endl;
	}
	summary.close();

	
	return 0;
}


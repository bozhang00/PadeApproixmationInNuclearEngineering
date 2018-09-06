#include <iostream>
#include <fstream>
#include <stdio.h>    
#include <vector>  
#include <stdlib.h>    
#include <time.h>
#include <string>
#include <math.h>
#include <cmath>

using namespace std;
#define MATHLIB_STANDALONE

void GENERATION(double& GE, double& GV1, double& GV2, double& GV3, double& GV4, double& GV5, double& GV6, double& GV7,
				double& GV8, double& GV9, double& GH, double REF_T, double REF_DENSITY)
{

	
	GE  = REF_T * (1.49/pow(10, 5))/REF_DENSITY;
	GV1 = REF_T * (9.91/pow(10, 6))/REF_DENSITY;
	GV2 = REF_T * (1.51/pow(10, 6))/REF_DENSITY;
	GV3 = REF_T * (2.60/pow(10, 7))/REF_DENSITY;
	GV4 = REF_T * (1.58/pow(10, 7))/REF_DENSITY;
	GV5 = REF_T * (6.29/pow(10, 8))/REF_DENSITY;
	GV6 = REF_T * 0.0;
	GV7 = REF_T * 0.0;
	GV8 = REF_T * 0.0;
	GV9 = REF_T * (3.16/pow(10, 8))/REF_DENSITY;
	GH  = REF_T * (2.11/pow(10, 11))/REF_DENSITY;

}	

void DIFFUSION(double& DE, double& DV1, double& DV2, double& DV3, double& DV4, double& DH1, double& DH2, double& DH3, 
			  double REF_T, double REF_DENSITY)
{
	double Ev_KB_T = (16202.1/(1.38*425));
	
	DE   = 4.0 * 3.1415926 * 1.00 * pow(10, 11) * exp(-0.34*Ev_KB_T)* REF_DENSITY * REF_T;
	DV1  = 4.0 * 3.1415926 * 1.00 * pow(10, 11) * exp(-0.67*Ev_KB_T)* REF_DENSITY * REF_T;
	DV2  = 4.0 * 3.1415926 * 0.50 * pow(10, 11) * exp(-0.62*Ev_KB_T)* REF_DENSITY * REF_T;
	DV3  = 4.0 * 3.1415926 * 0.33 * pow(10, 11) * exp(-0.37*Ev_KB_T)* REF_DENSITY * REF_T;
	DV4  = 4.0 * 3.1415926 * 0.25 * pow(10, 11) * exp(-0.48*Ev_KB_T)* REF_DENSITY * REF_T;
	DH1  = 4.0 * 3.1415926 * 1.00 * pow(10, 11) * exp(-0.06*Ev_KB_T)* REF_DENSITY * REF_T;
	DH2  = 4.0 * 3.1415926 * 0.50 * pow(10, 11) * exp(-0.06*Ev_KB_T)* REF_DENSITY * REF_T;
	DH3  = 4.0 * 3.1415926 * 0.33 * pow(10, 11) * exp(-0.06*Ev_KB_T)* REF_DENSITY * REF_T;
	
}	
		
void REACTION(int N_ic, double REF_DENSITY, double REF_T, 
			vector<vector<double> >& KV1, vector<vector<double> >& KV2, vector<vector<double> >& KV3, vector<vector<double> >& KV4, 
			vector<vector<double> >& KH1, vector<vector<double> >& KH2, vector<vector<double> >& KH3,
			vector<vector<double> >& KE,  vector<vector<double> >& KEI, vector<vector<double> >& KEV)
{
	int i, j, n, k, N_IC;
	double 	EV, EE, Omiga, a0, Ev_KB_T, DE, DV1, DV2, DV3, DV4, DH1, DH2, DH3, a0, Omiga, R0, R, R1, R2, R3, R4;
	
	N_IC = N_ic + 170;
	
	a0 = 0.287;
	Omiga=a0*a0*a0/2.0;
	R0 = 0.62;
	
	Ev_KB_T = (16202.1/(1.38*425));  // transfer unit to eV

	DE   = 4.0*3.1415926* 1.00 * pow(10, 11) * exp(-0.34*Ev_KB_T) * REF_DENSITY * REF_T;
	DV1  = 4.0*3.1415926* 1.00 * pow(10, 11) * exp(-0.67*Ev_KB_T) * REF_DENSITY * REF_T;
	DV2  = 4.0*3.1415926* 0.50 * pow(10, 11) * exp(-0.62*Ev_KB_T) * REF_DENSITY * REF_T;
	DV3  = 4.0*3.1415926* 0.33 * pow(10, 11) * exp(-0.37*Ev_KB_T) * REF_DENSITY * REF_T;	
	DV4  = 4.0*3.1415926* 0.25 * pow(10, 11) * exp(-0.48*Ev_KB_T) * REF_DENSITY * REF_T;
	DH1  = 4.0*3.1415926* 1.00 * pow(10, 11) * exp(-0.06*Ev_KB_T) * REF_DENSITY * REF_T;
	DH2  = 4.0*3.1415926* 0.50 * pow(10, 11) * exp(-0.06*Ev_KB_T) * REF_DENSITY * REF_T;
	DH3  = 4.0*3.1415926* 0.33 * pow(10, 11) * exp(-0.06*Ev_KB_T) * REF_DENSITY * REF_T;

	R1 = 0.287*cbrt(3.0/(8*3.1415926));
	R2 = 0.287*cbrt(6.0/(8*3.1415926));
	R3 = 0.287*cbrt(9.0/(8*3.1415926));
	R4 = 0.287*cbrt(12.0/(8*3.1415926));	

	KV1.resize(N_IC); 		KV2.resize(N_IC); 		KV3.resize(N_IC); 		KV4.resize(N_IC); 
	KH1.resize(N_IC); 		KH2.resize(N_IC); 		KH3.resize(N_IC); 
	KE.resize(N_IC);
	KEI.resize(N_IC);
	KEV.resize(N_IC);

	for(h=0; h<N_IC; h++){			
		KV1[h].resize(N_IC); 		KV2[h].resize(N_IC); 		KV3[h].resize(N_IC);		KV4[h].resize(N_IC); 
		KH1[h].resize(N_IC); 		KH2[h].resize(N_IC); 		KH3[h].resize(N_IC); 
		KE[h].resize(N_IC);	
		KEI[h].resize(N_IC);
		KEV[h].resize(N_IC);	
	}
	
	a0 = 0.287;  			
	Omiga=a0*a0*a0/2.0;				
	R0 = 0.62;
	
	Ev_KB_T = (16202.1/(1.38*425));  // transfer unit to eV
	
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
	
}
		

void DENSITY(int N_ic, vector<vector<vector<double> > >& a, vector<vector<vector<double> > >& ae, 
				vector<vector<vector<double> > >& av1, vector<vector<vector<double> > >& av2,
				vector<vector<vector<double> > >& av3, vector<vector<vector<double> > >& av4, vector<vector<vector<double> > >& ah1,
				vector<vector<vector<double> > >& ah2, vector<vector<vector<double> > >& ah3,
				vector<vector<vector<double> > >& gain_e, vector<vector<vector<double> > >& local_lost,
				vector<vector<vector<double> > >& pa, vector<vector<vector<double> > >& pb,	vector<int>& L, vector<vector<int> >& s,
				vector<double>& se, vector<double>& sv1, vector<double>& sv2, vector<double>& sv3, 
				vector<double>& sv4,vector<double>& sh1, vector<double>& sh2, vector<double>& sh3)
	{		
	int i, j, h,l, v, k, N_H_ROW, N_V_COL;
	
// **************** the summation of the rest ****************************************************	
	se.clear(); 			
	sv1.clear();			sv2.clear();			sv3.clear();			sv4.clear();	
	sh1.clear();			sh2.clear();			sh3.clear();
	
	se.resize(N_ic);		
	sv1.resize(N_ic);		sv2.resize(N_ic);		sv3.resize(N_ic);		sv4.resize(N_ic);	
	sh1.resize(N_ic);		sh2.resize(N_ic);		sh3.resize(N_ic);
	
// ******************   Label and s  ** *********************************************************	
	L.clear();				s.clear(); 		s.resize(1);
	L.push_back(N_ic-1);	for(v =0; v<N_ic; v++) s[0].push_back(v);
	
// *****************	the density term ********************************************************
	a.clear();	   	  a.resize(N_ic);	ae.clear();	   	  ae.resize(N_ic);
	av1.clear();	  av2.clear();		av3.clear(); 	  av4.clear();		ah1.clear(); 	  ah2.clear(); 		 ah3.clear(); 	
	av1.resize(N_ic); av2.resize(N_ic); av3.resize(N_ic); av4.resize(N_ic); ah1.resize(N_ic); ah2.resize(N_ic);  ah3.resize(N_ic); 	
	
	pa.clear();				pb.clear();
	pa.resize(N_ic);    	pb.resize(N_ic);	

	for(h=0;h<N_ic; h++){
	
		a[h].resize(1);		ae[h].resize(1);
		
		av1[h].resize(1); av2[h].resize(1);	av3[h].resize(1);  av4[h].resize(1); ah1[h].resize(1); ah2[h].resize(1); ah3[h].resize(1); 	
				
		pa[h].resize(1);		pb[h].resize(1);
		
		s.resize(1);   // the s only need to respons to 1 h
		
		//**********************************************************
		a[h][0].resize(N_ic);			ae[h][0].resize(N_ic);

		av1[h][0].resize(N_ic);			av2[h][0].resize(N_ic);			 av3[h][0].resize(N_ic); 		av4[h][0].resize(N_ic);	
		ah1[h][0].resize(N_ic); 		ah2[h][0].resize(N_ic); 		 ah3[h][0].resize(N_ic); 	
		
			
		pa[h][0].push_back(0.0);		pb[h][0].push_back(0.0);
	}
	
}





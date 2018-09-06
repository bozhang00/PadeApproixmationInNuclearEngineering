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

void DIFFUSION(double& DE, double& DV1, double& DV2, double& DV3, double& DV4, double& DH1, double& DH2, double& DH3)
{
	double Ev_KB_T = (16202.1/(1.38*425));
	DE   = 4.0 * 3.1415926 * 1.00 * pow(10, 11) * exp(-0.34*Ev_KB_T);
	DV1  = 4.0 * 3.1415926 * 1.00 * pow(10, 11) * exp(-0.67*Ev_KB_T);
	DV2  = 4.0 * 3.1415926 * 0.50 * pow(10, 11) * exp(-0.62*Ev_KB_T);	
	DV3  = 4.0 * 3.1415926 * 0.33 * pow(10, 11) * exp(-0.37*Ev_KB_T);	
	DV4  = 4.0 * 3.1415926 * 0.25 * pow(10, 11) * exp(-0.48*Ev_KB_T);	
	DH1  = 4.0 * 3.1415926 * 1.00 * pow(10, 11) * exp(-0.06*Ev_KB_T);	
	DH2  = 4.0 * 3.1415926 * 0.50 * pow(10, 11) * exp(-0.06*Ev_KB_T);	
	DH3  = 4.0 * 3.1415926 * 0.33 * pow(10, 11) * exp(-0.06*Ev_KB_T);
	
}	
		
void REACTION(int N_ic, double REF_DENSITY, double REF_T, vector<vector<vector<double> > >& KV, vector<vector<vector<double> > >& KH,
	vector<vector<double> >& KE, vector<vector<double> >& KBV, vector<vector<double> >& KBH, vector<vector<double> >& KBE)
{
	int i, j, n, k;
	double 	EV, EE, Omiga, a0, Ev_KB_T, DE, DV1, DV2, DV3, DV4, DH1, DH2, DH3, R0, R, RE, R1, R2, R3, R4;
	
	N_ic = N_ic + 170;
	
	KV.resize(4);
	KH.resize(3);
	
	KV[0].resize(N_ic);		// the reaction with v1
	KV[1].resize(N_ic);		// the reaction with v2		
	KV[2].resize(N_ic);		// the reaction with v3
	KV[3].resize(N_ic);		// the reaction with v4
	
	KH[0].resize(N_ic);
	KH[1].resize(N_ic);
	KH[2].resize(N_ic);
	
	KE.resize(N_ic);
	KBV.resize(N_ic);
	KBE.resize(N_ic);
	
	for(i=0;i<N_ic;i++){
		KV[0][i].resize(N_ic);
		KV[1][i].resize(N_ic);
		KV[2][i].resize(N_ic);
		KV[3][i].resize(N_ic);
		KH[0][i].resize(N_ic);
		KH[1][i].resize(N_ic);
		KH[2][i].resize(N_ic);
		KE[i].resize(N_ic);
		KBV[i].resize(N_ic);
		KBE[i].resize(N_ic);
	}
	
	a0 = 0.287;
	Omiga=a0*a0*a0/2.0;
	R0 = 0.62;
	
	Ev_KB_T = (16202.1/(1.38*425));  // transfer unit to eV

	DE   = 4.0*3.1415926* 1.0*pow(10, 11)* exp(-0.34*Ev_KB_T);
	DV1  = 4.0*3.1415926* 1.0*pow(10, 11)* exp(-0.67*Ev_KB_T);
	DV2  = 4.0*3.1415926* 0.5*pow(10, 11)* exp(-0.62*Ev_KB_T);	
	DV3  = 4.0*3.1415926* 0.33*pow(10, 11)*exp(-0.37*Ev_KB_T);	
	DV4  = 4.0*3.1415926* 0.25*pow(10, 11)*exp(-0.48*Ev_KB_T);	
	DH1  = 4.0*3.1415926* 1.0*pow(10, 11)* exp(-0.06*Ev_KB_T);	
	DH2  = 4.0*3.1415926* 0.5*pow(10, 11)* exp(-0.06*Ev_KB_T);	
	DH3  = 4.0*3.1415926* 0.33*pow(10, 11)*exp(-0.06*Ev_KB_T);
	
	RE = 0.287*cbrt(3.0/(8*3.1415926));
	R1 = 0.287*cbrt(3.0/(8*3.1415926));
	R2 = 0.287*cbrt(6.0/(8*3.1415926));
	R3 = 0.287*cbrt(9.0/(8*3.1415926));
	R4 = 0.287*cbrt(12.0/(8*3.1415926));
	
	
	for(i=0;i<N_ic;i++){
		for(j=0;j<N_ic;j++){
			R = 0.287*cbrt(3.0*(i+j)/(8*3.1415926));
			KV[0][i][j] = DV1*(R1 + R + R0)*REF_DENSITY*REF_T;
			KV[1][i][j] = DV2*(R2 + R + R0)*REF_DENSITY*REF_T;
			KV[2][i][j] = DV3*(R3 + R + R0)*REF_DENSITY*REF_T;
			KV[3][i][j] = DV4*(R4 + R + R0)*REF_DENSITY*REF_T;
			KH[0][i][j] = DH1*(R1 + R + R0)*REF_DENSITY*REF_T;
			KH[1][i][j] = DH2*(R2 + R + R0)*REF_DENSITY*REF_T;
			KH[2][i][j] = DH3*(R3 + R + R0)*REF_DENSITY*REF_T;
			KE[i][j] = DE*(R1 + R + R0)*REF_DENSITY*REF_T;

		}
	}
/*	
// reactio with E
	KV[0][0][0] = (DE+DV1 )*(R1+R0)*REF_DENSITY*REF_T;
	KV[1][0][0] = (DE+DV2 )*(R2+R0)*REF_DENSITY*REF_T;
	KV[2][0][0] = (DE+DV3 )*(R3+R0)*REF_DENSITY*REF_T;
	KV[3][0][0] = (DE+DV4 )*(R4+R0)*REF_DENSITY*REF_T;
	KH[0][0][0] = (DE+DH1 )*(R1+R0)*REF_DENSITY*REF_T;
	KH[1][0][0] = (DE+DH2 )*(R2+R0)*REF_DENSITY*REF_T;
	KH[2][0][0] = (DE+DH3 )*(R3+R0)*REF_DENSITY*REF_T;
	
// react with V1	
	KV[0][0][1] = (DV1+ DV1)*(R1 + R1 + R0)*REF_DENSITY*REF_T;
	KV[1][0][1] = (DV1+ DV2)*(R2 + R1 + R0)*REF_DENSITY*REF_T;
	KV[2][0][1] = (DV1+ DV3)*(R3 + R1 + R0)*REF_DENSITY*REF_T;
	KV[3][0][1] = (DV1+ DV4)*(R4 + R1 + R0)*REF_DENSITY*REF_T;
	KH[0][0][1] = (DV1+ DH1)*(R1 + R1 + R0)*REF_DENSITY*REF_T;
	KH[1][0][1] = (DV1+ DH2)*(R2 + R1 + R0)*REF_DENSITY*REF_T;
	KH[2][0][1] = (DV1+ DH3)*(R3 + R1 + R0)*REF_DENSITY*REF_T;
	
// react with V2	
	KV[0][0][2] = (DV2+ DV1)*(R1 + R2 + R0)*REF_DENSITY*REF_T;
	KV[1][0][2] = (DV2+ DV2)*(R2 + R2 + R0)*REF_DENSITY*REF_T;
	KV[2][0][2] = (DV2+ DV3)*(R3 + R2 + R0)*REF_DENSITY*REF_T;
	KV[3][0][2] = (DV2+ DV4)*(R4 + R2 + R0)*REF_DENSITY*REF_T;
	KH[0][0][2] = (DV2+ DH1)*(R1 + R2 + R0)*REF_DENSITY*REF_T;
	KH[1][0][2] = (DV2+ DH2)*(R2 + R2 + R0)*REF_DENSITY*REF_T;
	KH[2][0][2] = (DV2+ DH3)*(R3 + R2 + R0)*REF_DENSITY*REF_T;

// react with V3	
	KV[0][0][3] = (DV3 + DV1)*(R1 + R3 + R0)*REF_DENSITY*REF_T;
	KV[1][0][3] = (DV3 + DV2)*(R2 + R3 + R0)*REF_DENSITY*REF_T;
	KV[2][0][3] = (DV3 + DV3)*(R3 + R3 + R0)*REF_DENSITY*REF_T;
	KV[3][0][3] = (DV3 + DV4)*(R4 + R3 + R0)*REF_DENSITY*REF_T;
	KH[0][0][3] = (DV3 + DH1)*(R1 + R3 + R0)*REF_DENSITY*REF_T;
	KH[1][0][3] = (DV3 + DH2)*(R2 + R3 + R0)*REF_DENSITY*REF_T;
	KH[2][0][3] = (DV3 + DH3)*(R3 + R3 + R0)*REF_DENSITY*REF_T;
	
// react with V4	
	KV[0][0][4] = (DV4 + DV1)*(R1 + R4 + R0)*REF_DENSITY*REF_T;
	KV[1][0][4] = (DV4 + DV2)*(R2 + R4 + R0)*REF_DENSITY*REF_T;
	KV[2][0][4] = (DV4 + DV3)*(R3 + R4 + R0)*REF_DENSITY*REF_T;
	KV[3][0][4] = (DV4 + DV4)*(R4 + R4 + R0)*REF_DENSITY*REF_T;
	KH[0][0][4] = (DV4 + DH1)*(R1 + R4 + R0)*REF_DENSITY*REF_T;
	KH[1][0][4] = (DV4 + DH2)*(R2 + R4 + R0)*REF_DENSITY*REF_T;
	KH[2][0][4] = (DV4 + DH3)*(R3 + R4 + R0)*REF_DENSITY*REF_T;


// react with H1	
	KV[0][1][0] = (DH1+ DV1)*(R1 + R1 + R0)*REF_DENSITY*REF_T;
	KV[1][1][0] = (DH1+ DV2)*(R2 + R1 + R0)*REF_DENSITY*REF_T;
	KV[2][1][0] = (DH1+ DV3)*(R3 + R1 + R0)*REF_DENSITY*REF_T;
	KV[3][1][0] = (DH1+ DV4)*(R4 + R1 + R0)*REF_DENSITY*REF_T;
	KH[0][1][0] = (DH1+ DH1)*(R1 + R1 + R0)*REF_DENSITY*REF_T;
	KH[1][1][0] = (DH1+ DH2)*(R2 + R1 + R0)*REF_DENSITY*REF_T;
	KH[2][1][0] = (DH1+ DH3)*(R3 + R1 + R0)*REF_DENSITY*REF_T;
	
// react with V2	
	KV[0][2][0] = (DH2+ DV1)*(R1 + R2 + R0)*REF_DENSITY*REF_T;
	KV[1][2][0] = (DH2+ DV2)*(R2 + R2 + R0)*REF_DENSITY*REF_T;
	KV[2][2][0] = (DH2+ DV3)*(R3 + R2 + R0)*REF_DENSITY*REF_T;
	KV[3][2][0] = (DH2+ DV4)*(R4 + R2 + R0)*REF_DENSITY*REF_T;
	KH[0][2][0] = (DH2+ DH1)*(R1 + R2 + R0)*REF_DENSITY*REF_T;
	KH[1][2][0] = (DH2+ DH2)*(R2 + R2 + R0)*REF_DENSITY*REF_T;
	KH[2][2][0] = (DH2+ DH3)*(R3 + R2 + R0)*REF_DENSITY*REF_T;
		
// react with V3	
	KV[0][3][0] = (DH3 + DV1)*(R1 + R3 + R0)*REF_DENSITY*REF_T;
	KV[1][3][0] = (DH3 + DV2)*(R2 + R3 + R0)*REF_DENSITY*REF_T;
	KV[2][3][0] = (DH3 + DV3)*(R3 + R3 + R0)*REF_DENSITY*REF_T;
	KV[3][3][0] = (DH3 + DV4)*(R4 + R3 + R0)*REF_DENSITY*REF_T;
	KH[0][3][0] = (DH3 + DH1)*(R1 + R3 + R0)*REF_DENSITY*REF_T;
	KH[1][3][0] = (DH3 + DH2)*(R2 + R3 + R0)*REF_DENSITY*REF_T;
	KH[2][3][0] = (DH3 + DH3)*(R3 + R3 + R0)*REF_DENSITY*REF_T;
	
	
	for(i=0; i<N_ic;i++){
		for(j=1; j<N_ic;j++){
			
			EE = 4.33-5.76*(pow(j, 0.66) - pow(j-1, 0.666));
			EV = 4.88+2.59*(pow(j, 0.66) - pow(j-1, 0.666))-2.5*log(1+i/j);
			EV = 1.73-2.59*(pow(j, 0.66) - pow(j-1, 0.666))+2.5*log(1+i/j);
			
			KE[i][j] = DE*(R1 + R + R0)*exp(EE *Ev_KB_T)/Omiga;
		//	KBV[i][j] = KV[0][i][j]*exp(EV *Ev_KB_T)/Omiga;
		//	KBE[i][j] = DV1*(R1 + R + R0)*exp(0.06*Ev_KB_T)/Omiga;
		}
	}
	
	*/
}
		

void DENSITY(int N_ic, vector<vector<vector<double> > >& a, vector<vector<vector<double> > >& ae, 
				vector<vector<vector<double> > >& av1, vector<vector<vector<double> > >& av2,
				vector<vector<vector<double> > >& av3, vector<vector<vector<double> > >& av4, vector<vector<vector<double> > >& ah1,
				vector<vector<vector<double> > >& ah2, vector<vector<vector<double> > >& ah3,
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





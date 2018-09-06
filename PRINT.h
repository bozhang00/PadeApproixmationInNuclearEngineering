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

using namespace std;
#define MATHLIB_STANDALONE

int outputstep(int a)
{
	int b = 10;
	if (a <= 1000) { b = 1000; }
	else if (a > 1000 			&& a <= 100000) 			{ b = 10000; }
	else if (a > 100000 		&& a <= 1000000)			{ b = 100000; }
	else if (a > 1000000 		&& a <= 10000000) 			{ b = 1000000; }
	else if (a > 10000000 		&& a <= 1000000000) 		{ b = 10000000; }
	else if (a > 100000000 		&& a <= 10000000000)		{ b = 100000000; }
	else if (a > 1000000000		&& a <= 100000000000) 		{ b = 1000000000; }
	
	else { b = 1000000000; }
	return b;
} 


void print(vector<vector<vector<double> > >& a, vector<vector<int> >& s, vector<vector<vector<double> > >& pa, 
			vector<vector<vector<double> > >& pb, int st, vector<vector<int> >& indicator)
{
	int h, l, v, N_H_ROW, N_LABLE, N_V_COL; 
	double pv;
	
	N_H_ROW = a.size();
	N_LABLE = a[0].size();

	ofstream summary, inforecord;
	
	inforecord.open("record.txt", std::ios_base::app);
	inforecord << "step is,   " << st << ",  the total is "<< s[N_LABLE-1][s[N_LABLE-1].size()-1] <<"   divided into "<<N_LABLE<<" regions "<< endl;
	for (l = 0; l < N_LABLE; l++) { 
		for (v = 0; v < s[l].size(); v++) 		inforecord << s[l][v] << ", ";
	}
	inforecord <<"  \n "<< endl;
	
	inforecord.close();
	
	
	//summary.open("data.txt", std::ios_base::app);
	summary.open("data.txt");
	N_H_ROW = a.size();
	N_LABLE = a[0].size();
			
	for(h =0; h<N_H_ROW; h++){
		for( l = 0; l<N_LABLE; l++ ) {
			if( indicator[h][l] == 0 ) {	
				N_V_COL = a[h][l].size(); 
				for(v =0; v< N_V_COL ; v++)			summary << a[h][l][v] << " ";					
			}
			else{
				for (v = s[l][0]; v <= s[l][s[l].size() - 1]; v++) summary <<getvalue(pa[h][l], pb[h][l], v) << " ";
			}
		}
		//summary<<endl;
	}

	//summary << endl;
	summary.close();
	
}


void finalprint(int N_IC, double BIAS, int BIN_I, int BINSIZEA, int BINSIZEB, int BINSTEP, int POINT_N, int STEP, 
			double ERROR_T, double REF_T, double REF_DENSITY, double dt, float time,  
			vector<vector<vector<double> > >& a, vector<vector<int> >& s, vector<vector<vector<double> > >& pa, 
			vector<vector<vector<double> > >& pb, int st, vector<vector<int> >& indicator)
{
	int h, l, v, N_H_ROW, result, N_LABLE, N_V_COL; 
	double pv;
	char oldname[] ="data.txt";
	char filename[100] = "00000000000000000000000000000";
	char newname[100] = "00000000000000000000000000000";	
		
	ofstream summary;
	
	N_H_ROW = a.size();
	N_LABLE = a[0].size();
	
	if(BIN_I == 1) sprintf(filename,"data/BIN STEP %d dt %f REF_t %f  Total %d Label %d Spend_time %f.txt",STEP, dt, REF_T, s[N_LABLE-1][s[N_LABLE-1].size()-1], N_LABLE, time);
	else sprintf(filename,"data/UN-BIN Step %d Dt %f Ref_t %f  Total %d Label %d Spend_time %f.txt", STEP, dt, REF_T, s[N_LABLE-1][s[N_LABLE-1].size()-1], N_LABLE, time);

	if(BIN_I == 1) sprintf(newname,"data/BIN step %d dt %f ref_t %f  total %d label %d spend_time %f data.txt",  STEP, dt, REF_T, s[N_LABLE-1][s[N_LABLE-1].size()-1], N_LABLE, time);
	else sprintf(newname,"data/UN-BIN step %d dt %f ref_t %f  total %d label %d spend_time %f data.txt",  STEP, dt, REF_T, s[N_LABLE-1][s[N_LABLE-1].size()-1], N_LABLE, time);
	
	
	result= rename( oldname , newname );

	
	summary.open(filename);
	
	summary<<"check bin or not  ";
	if(BIN_I == 1) summary<<" -- BIN "<<endl;
	else summary<<" -- NOT BIN "<<endl;

	summary<<"bin size   "<<BINSIZEA<<endl;
	summary<<"update every    "<<BINSTEP<<"  step  "<<endl;
	summary<<" There are  "<<POINT_N<<"  points in the pade equation "<<endl;
	summary<<" total step is     "<<STEP<<endl;
	summary<<" the reference time is (in second)  "<<REF_T<<endl;
	summary<<" the error threshold is  "<<ERROR_T<<endl;
	summary<<" the reference density is  "<<REF_DENSITY<<endl;
	summary<<" IC # of h and v  is  "<<N_IC<<endl;
	summary<<" Dislocation bias   "<<BIAS <<endl;	
	summary<<" step dt is   "<<dt <<endl;		
	summary<<"**********************************************"<<endl;	
	
	summary << "step is,   " << st << ",  the total is "<< s[N_LABLE-1][s[N_LABLE-1].size()-1] <<"   divided into "<<N_LABLE<<" regions "<< endl;
	
	
	
	for (l = 0; l < N_LABLE; l++) { 
		for (v = 0; v < s[l].size(); v++) 		summary << s[l][v] << ", ";
	}
	summary <<"  \n "<< endl;
		
	for(h =0; h<N_H_ROW; h++){
		for( l = 0; l<N_LABLE; l++ ) {
			if(indicator[h][l] == 0 ) {	
				N_V_COL = a[h][l].size(); 
				for(v =0; v< N_V_COL ; v++){ 	summary << a[h][l][v] << ", ";			}
			} 
			else{
				for (v = s[l][0]; v <= s[l][s[l].size() - 1]; v++) summary <<getvalue(pa[h][l], pb[h][l], v) << ", ";
			}
		}
		summary<<endl;
	}
	
	
	for(h =0; h<N_H_ROW; h++){
		for( l = 0; l<N_LABLE; l++ ) {
			if(indicator[h][l] == 0 ) {		summary<<" U ";	}
			else{		 					summary<<" B ";	}
		}
		summary<<endl;
	}	

	summary << endl;
	summary.close();
	
}

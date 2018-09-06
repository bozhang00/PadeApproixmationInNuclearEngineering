
double getvalue(vector<double>& a,vector<double>& b, double unknown){
    int i,  aN = a.size(),  bM = b.size();
    double up=0.0, down=1.0, y=0.0;
    for(i=0;i<aN;i++) up=up+a[i]*pow(unknown,i);
    for(i=0;i<bM;i++) down=down+b[i]*pow(unknown,i+1);
    return y=up/down;
}

vector<double> gauss(vector< vector<double> > A) {
    int n = A.size(), i,j,maxRow,k;
	double maxEl,tmp,c;
    vector<double> x(n);
    for (int i = 0; i<n; i++) {
		// Search for maximum in this column
        maxEl = fabs(A[i][i]);  		
        maxRow = i;
        for (k = i + 1; k<n; k++) {
            if (abs(A[k][i]) > maxEl) {
                maxEl = abs(A[k][i]);
                maxRow = k;
            }
        }   // Swap maximum row with current row (column by column)
        for (k = i; k<n + 1; k++) {  	 
            tmp = A[maxRow][k];
            A[maxRow][k] = A[i][k];
            A[i][k] = tmp;
        } // Make all rows below this one 0 in current column
        for (k = i + 1; k<n; k++) {   
            c = -A[k][i] / A[i][i];
            for (j = i; j<n + 1; j++) {
                if (i == j) {  A[k][j] = 0; }
                else {	A[k][j] += c * A[i][j]; }
            }
        }
    }	// Solve equation Ax=b for an upper triangular matrix A
    for (i = n - 1; i >= 0; i--) { 
        x[i] = A[i][n] / A[i][i];
        for (k = i - 1; k >= 0; k--) { A[k][n] -= A[k][i] * x[i]; }
    }
    return x;
}

void getfactor(vector<double>& x, vector<double>& f,vector<double>& a, vector<double>& b){
    int i,j;
	int N=a.size(),M=b.size(),k;
	k=N+M;
    vector<vector<double> > A(k);	for (i = 0; i < k; i++) { A[i].resize(k+1); }
    for(i=0;i<k;i++){
        for(j=0;j<N;j++) A[i][j]=pow(x[i],j);
        for(j=N;j<k;j++) A[i][j]=-1.0*f[i]*pow(x[i],(j-N+1));
        A[i][k]= f[i];
    }
    x = gauss(A);
    for(i=0;i<N;i++)   a[i]=x[i];
    for(i=0;i<M;i++)   b[i]=x[i+N]; 
}
					
void PADESTEP(vector<vector<double> >& FULL, vector<vector<double> >& KE,
			vector<vector<vector<double> > >& KV, vector<vector<vector<double> > >& KH, 
			vector<vector<vector<double> > >& a, vector<double>& se, vector<double>& sv1, vector<double>& sv2, 
			vector<double>& sv3, vector<double>& sv4, vector<double>& sh1, vector<double>& sh2, vector<double>& sh3, 
			double dt, int step,  int CAL_ROW, vector<int>& L,   
			double GE, double GV1,double GV2, double GV3, double GV4, double GV5, double GV6, double GV7, double GV8, double GV9, double GH)
{
	
	int h, v, i, j,l, k, N_COL, N_LABLE, LABLE_1, N_H_ROW;
	double te, tv1, tv2, tv3, tv4, th1, th2, th3, c, cb, e, v1, v2, v3, v4, h1, h2, h3, vf1, vf2, vf3, vf4, hf1, hf2, hf3, per;
	double kvf1c, kvf2c, kvf3c, kvf4c, kvf1, kvf2, kvf3, kvf4, khf1c, khf2c, khf3c, khf1, khf2, khf3, ke, keb;
	double t;
	vector<vector<double> > pe, pv1, pv2, pv3, pv4, ph1, ph2, ph3;
	
	t=1.0*dt*step*1.0;
	
	N_H_ROW = FULL.size();
	N_COL = FULL[0].size();
	LABLE_1	= L[0];
	N_LABLE = L.size();
	
	te = 0.0; 	tv1 = 0.0; 	 tv2 = 0.0;  	tv3 = 0.0;  	tv4 = 0.0;	th1 = 0.0; 		 th2 = 0.0; 	th3 = 0.0;
		

	pe.resize(N_H_ROW);
	pv1.resize(N_H_ROW); 	pv2.resize(N_H_ROW); 	pv3.resize(N_H_ROW);	pv4.resize(N_H_ROW);
	ph1.resize(N_H_ROW); 	ph2.resize(N_H_ROW); 	ph3.resize(N_H_ROW); 
	
	for(h=0; h<N_H_ROW; h++){
		pe[h].resize(N_COL); 
		pv1[h].resize(N_COL); 		pv2[h].resize(N_COL); 		pv3[h].resize(N_COL);			pv4[h].resize(N_COL);
		ph1[h].resize(N_COL); 		ph2[h].resize(N_COL); 		ph3[h].resize(N_COL); 		
	}
	
	FULL[0][0] = FULL[0][0] + GE * dt;		 			 
	FULL[0][1] = FULL[0][1] + GV1 * dt;		 			
	FULL[0][2] = FULL[0][2] + GV2 * dt;		 			
	FULL[0][3] = FULL[0][3] + GV3 * dt; 			
	FULL[0][4] = FULL[0][4] + GV4 * dt;			
	FULL[0][5] = FULL[0][5] + GV5 * dt;
	FULL[0][6] = FULL[0][6] + GV6 * dt;
	FULL[0][7] = FULL[0][7] + GV7 * dt;
	FULL[0][8] = FULL[0][8] + GV8 * dt;
	FULL[0][9] = FULL[0][9] + GV9 * dt;
	FULL[1][0] = FULL[1][0] + GH * dt;
	
	e  = FULL[0][0];  	v1 = FULL[0][1];  	v2 = FULL[0][2];   v3 = FULL[0][3];  v4 = FULL[0][4];   
	h1 = FULL[1][0];	h2 = FULL[2][0];	h3 = FULL[3][0];	
	
	for(h=0; h<CAL_ROW; h++){
		for(v=0; v<N_COL; v++){
			c = FULL[h][v];	
			
			ke = KE[h][v]; 	keb = KE[h][v+1];
			
			if(v == FULL[h].size() ==1 ) cb= 0.0;
			else cb = FULL[h][v+1];
				
			kvf1c = KV[0][h][v]; 		if(v - 1 < 0) { kvf1 = 0.0; } else { kvf1 = KV[0][h][v-1]; }
			kvf2c = KV[1][h][v];		if(v - 2 < 0) { kvf2 = 0.0; } else { kvf2 = KV[1][h][v-2]; }
			kvf3c = KV[2][h][v];		if(v - 3 < 0) { kvf3 = 0.0; } else { kvf3 = KV[2][h][v-3]; }	
			kvf4c = KV[3][h][v];		if(v - 4 < 0) { kvf4 = 0.0; } else { kvf4 = KV[3][h][v-4]; }
				
			khf1c = KH[0][h][v];		if( h -1 < 0) { khf1 = 0.0; } else { khf1 = KH[0][h-1][v]; }
			khf2c = KH[1][h][v]; 		if( h -2 < 0) { khf2 = 0.0; } else { khf2 = KH[1][h-2][v]; }	
			khf3c = KH[2][h][v]; 		if( h -3 < 0) { khf3 = 0.0; } else { khf3 = KH[2][h-3][v]; }	
				
			if(v < 4 ){
				switch (v){
					case 0: vf1 = 0.0;				vf2 = 0.0; 				vf3 = 0.0;				vf4 = 0.0;	break;
					case 1:	vf1 = FULL[h][v-1];		vf2 = 0.0; 				vf3 = 0.0;				vf4 = 0.0;	break;
					case 2: vf1 = FULL[h][v-1];		vf2 = FULL[h][v-2];		vf3 = 0.0;				vf4 = 0.0;	break;	
					case 3: vf1 = FULL[h][v-1];		vf2 = FULL[h][v-2];	 	vf3 = FULL[h][v-3];		vf4 = 0.0;	break;
				}
			}
			else{ 	vf1 = FULL[h][v-1];		vf2 = FULL[h][v-2];	    	vf3 = FULL[h][v-3]; 	vf4 = FULL[h][v-4];	}				
						
			if(h<3 ){
				hf1 = 0.0;  	   hf2 = 0.0; 		 hf3 = 0.0;
				switch (h){
					case 1: 	hf1 = FULL[h-1][v]; 	break;
					case 2:  	hf1 = FULL[h-1][v];		hf2 = FULL[h-2][v];		break;
				}
			}
			else{ 	hf1 = FULL[h-1][v]; 	hf2 = FULL[h-2][v]; 		hf3 = FULL[h-3][v]; }
				
			pv1[h][v]= kvf1*vf1 - kvf1c*c; 		
			pv2[h][v]= kvf2*vf2 - kvf2c*c;
			pv3[h][v]= kvf3*vf3 - kvf3c*c;
			pv4[h][v]= kvf4*vf4 - kvf4c*c;						
				
			ph1[h][v]= khf1*hf1 - khf1c*c;				
			ph2[h][v]= khf2*hf2 - khf2c*c;
			ph3[h][v]= khf3*hf3 - khf3c*c;		
			pe[h][v] = cb*keb-c*ke;
			
			te = te - c*ke;
			tv1= tv1 - c*kvf1c; 						
			tv2= tv2 - c*kvf2c;
			tv3= tv3 - c*kvf3c;						
			tv4= tv4 - c*kvf4c;
				
			th1= th1 - c*khf1c;						
			th2= th2 - c*khf2c;
			th3= th3 - c*khf3c;
		}
		
	}
	per = 1.0; 	if(te != 0 & e !=0 )   { per = FULL[0][0]/( - te * e * dt ); if( per >1 ){ per = 1.0; }   }  e = per *e;	
	per = 1.0; 	if(tv1 != 0 & v1 !=0 ) { per = FULL[0][1]/( - tv1 * v1 * dt ); if( per >1 ){ per = 1.0; } } v1 = per *v1;
	per = 1.0; 	if(tv2 != 0 & v2 !=0 ) { per = FULL[0][2]/( - tv2 * v2 * dt ); if( per >1 ){ per = 1.0; } } v2 = per *v2;
	per = 1.0; 	if(tv3 != 0 & v3 !=0 ) { per = FULL[0][3]/( - tv3 * v3 * dt ); if( per >1 ){ per = 1.0; } } v3 = per *v3;
	per = 1.0; 	if(tv4 != 0 & v4 !=0 ) { per = FULL[0][4]/( - tv4 * v4 * dt ); if( per >1 ){ per = 1.0; } } v4 = per *v4;
		
	per = 1.0;	if(th1 != 0 & h1 !=0 ) { per = FULL[1][0]/( - th1 * h1 * dt ); if( per >1 ){ per = 1.0; } }	h1 = per *h1; 
	per = 1.0; 	if(th2 != 0 & h2 !=0 ) { per = FULL[2][0]/( - th2 * h2 * dt ); if( per >1 ){ per = 1.0; } }	h2 = per *h2;
	per = 1.0; 	if(th3 != 0 & h3 !=0 ) { per = FULL[3][0]/( - th3 * h3 * dt ); if( per >1 ){ per = 1.0; } }	h3 = per *h3; 
	
	pe[0][0] = pe[0][0] + te;	
	pv1[0][1] = pv1[0][1] + tv1;
	pv2[0][2] = pv2[0][2] + tv2;
	pv3[0][3] = pv3[0][3] + tv3;
	pv4[0][4] = pv4[0][4] + tv4;

	ph1[1][0] = ph1[1][0] + th1;
	ph2[2][0] = ph2[2][0] + th2;
	ph3[3][0] = ph3[3][0] + th3;

// sink strenght  at the foot part of this file blocked 
	for( h =0; h<CAL_ROW; h++){ 
		for(v =LABLE_1+1; v<N_COL; v++){	  	
			FULL[h][v] = FULL[h][v]+(v1*pv1[h][v] + v2*pv2[h][v] + v3*pv3[h][v] + v4*pv4[h][v]+
									 e*pe[h][v]+
									 h1*ph1[h][v] + h2*ph2[h][v] + h3*ph3[h][v])*t; 
			if( FULL[h][v] < 0.0 ) FULL[h][v] = 0.0;
		}
	}

	for( h =0; h<N_H_ROW; h++){ 	
		se[h] = 0.0;    sv1[h]=0.0;		sv2[h]=0.0;		sv3[h]=0.0;		sv4[h]=0.0;		sh1[h]=0.0;		sh2[h]=0.0;		sh3[h]=0.0;
		for(v =LABLE_1; v<N_COL; v++){	
			c = FULL[h][v];
			se[h] = se[h]-e*KE[h][v];
			sv1[h]=sv1[h]-c*KV[0][h][v];	sv2[h]=sv2[h]-c*KV[1][h][v];			
			sv3[h]=sv3[h]-c*KV[2][h][v];	sv4[h]=sv4[h]-c*KV[3][h][v];			
			sh1[h]=sh1[h]-c*KH[0][h][v];	sh2[h]=sh2[h]-c*KH[1][h][v];	sh3[h]=sh3[h]-c*KH[2][h][v];			
		}
	}	
	
	a.clear();	 	
	a.resize(N_H_ROW);		
	for(h=0; h<N_H_ROW; h++){	a[h].resize(N_LABLE);  }
	
	for(h=0; h<N_H_ROW; h++){
		for(v=0; v<=L[0]; v++)  a[h][0].push_back( FULL[h][v] );
		for(l=1; l<N_LABLE; l++){
			for(v=L[l-1]+1; v<=L[l]; v++) a[h][l].push_back( FULL[h][v] );
		}
	}
	
}


void GET_PADE(vector<double>& a, vector<int>& s, vector<double>& pa, vector<double>& pb, int POINT_N, double REF_T, int& indic )
{
	int i, j, k, v, index, bin_space, num, begin, end;
	double error_check, error_value;
	vector<double> x, fx, y, s_note;
	
	for(i =0 ; i<a.size(); i++) y.push_back( a[i] );

	begin = s[0];   	end = s[s.size()-1];		num = 1+ end - begin;
	
	x.resize(POINT_N);					fx.resize(POINT_N); 				s_note.resize(POINT_N);
	x[0]=begin;							fx[0]=y[0];							s_note[0] = begin;
		
	for(i = 1; i< POINT_N-1; i++){
		index = int(bin_space*(i)); 	x[i] = begin + index; 		fx[i] = y[index];		s_note[i] = begin+index;
	}
	
	x[POINT_N - 1] = end;				fx[POINT_N - 1]= y[y.size()-1];		s_note[POINT_N - 1] = end;	

	pa.resize(4);   	pb.resize(2);
	
	getfactor(x, fx, pa, pb);

	/*
	error_check = 0.0;	 error_value = 0.0;
	
	for(j=0; j<y.size(); j++){
		error_value = fabs( y[j] - getvalue(pa, pb, s_note[j])) / y[j];
		if( error_value > error_check) error_check = error_value;
	}
	*/
}



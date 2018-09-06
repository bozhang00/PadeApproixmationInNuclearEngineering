
void RECOVER_KEEP_DATA(vector<double>& full, vector<vector<double> >& a, vector<vector<int> > s, 
				 vector<vector<double> >& pa, vector<vector<double> >& pb, vector<int>& indicator)
{
	int i, l; 
	
	for (l = 0; l < a.size(); l++) {
		if(indicator[l] == 0 ){
			for (i = 0; i < a[l].size(); i++) full.push_back( a[l][i] );	
		}
		if(indicator[l] == 1){
			for (i = s[l][0]; i < s[l][s[l].size()-1]; i++) full.push_back( getvalue(pa[l], pb[l], i) ); 
		}
	}
}

void NEW_LABLE(vector<int>& label, int first_label, vector<vector<int> >& s, int size, int total_number)
{
	int i, j, l;
	
	label.clear();	
	label.push_back(first_label-1);
	
	for(i = first_label; i<total_number; ){
		i=i+size;
		
		if(i<( total_number-5 ) ) {	label.push_back(i-1); }
		else {
			label.push_back(total_number-1); 
			i = total_number+12; 
		}
	}

	s.clear();
	s.resize(label.size());
	for(i = 0; i<= label[0]; i++ ) s[0].push_back(i);
	
	for(l = 1 ; l< label.size(); l++){
		for(i = label[l-1]+1; i<=label[l]; i++ )  s[l].push_back(i);
	}

}

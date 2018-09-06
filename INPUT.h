void INPUTIC(int& N_IC, double& BIAS, int& BIN_I, int& BINSIZEA, int& BINSIZEB, int& BINSTEP, int& POINT_N, int& STEP, 
			double& ERROR_T, double& REF_T, double& REF_DENSITY, double& dt, int& record)
{
    char note_of_number[20];  // char should have enough space to read the letter in the note book
	
	ifstream myfile;
	
	myfile.open("input.txt");  
	myfile >> note_of_number;
	myfile >> BIN_I;
	myfile >> note_of_number;	
	myfile >> BINSIZEA;  
	myfile >> note_of_number;	
	myfile >> BINSIZEB; 	
	myfile >> note_of_number;	
	myfile >> BINSTEP;
	myfile >> note_of_number;	
	myfile >> POINT_N;
	myfile >> note_of_number;	
	myfile >> STEP;
	myfile >> note_of_number;	
	myfile >> REF_T;
	myfile >> note_of_number;	
	myfile >> ERROR_T;
	myfile >> note_of_number;	
	myfile >> REF_DENSITY;	
	myfile >> note_of_number;	
	myfile >> N_IC;	
	myfile >> note_of_number;	
	myfile >> BIAS;
	myfile >> note_of_number;	
	myfile >> dt;
	myfile >> note_of_number;	
	myfile >> record;
	myfile.close();
	
	ofstream ICfile;
 	ICfile.open("data.txt");  	ICfile.close();
 	ICfile.open("data.txt");  	ICfile.close();	
	

	
}

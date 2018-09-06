#include <iostream>
#include <fstream>
#include <vector>
 
using namespace std;
 
int main()
{
    vector<double> v;
     
    ifstream infile;
     
    infile.open("data.txt");
     
    double tmp;
    while(!infile.eof())
    {
        infile>>tmp;
        v.push_back(tmp);
    }
     
    infile.close();
     
	for( int i = 0; i < 6; i++){
       cout<< v[i];
	}
	 
    return 0;
}
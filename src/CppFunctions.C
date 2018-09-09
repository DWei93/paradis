#include <cstring>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <sstream>

#include <stdlib.h>
#include <algorithm>
#include <numeric>

using namespace std;

vector<string> split(const string& str, const string& delim) {
    vector<string> res;
    if("" == str) return res;
    char * strs = new char[str.length() + 1] ; 
    strcpy(strs, str.c_str()); 
    char * d = new char[delim.length() + 1];
    strcpy(d, delim.c_str());
    char *p = strtok(strs, d);
    while(p) {
        string s = p; 
        res.push_back(s); 
        p = strtok(NULL, d);

    }
    return res;
}

int ReadRandomForceTable(double **randomForceTable, double para[5])
{
    
    string          str;
    vector<string>  line;
    int             lineNo = 0, i = 0, numData;
    ifstream        infile;

    
    infile.open("Tec.plt");
    if(!infile){
        cout << "ReadRandomForceTable: can not open thefile" << endl;
        exit(-1);
    }

    while(getline(infile,str)){   
        line = split(str, " ");

        lineNo++;

        if(lineNo == 1){
            numData = atoi(line[0].c_str());
            *randomForceTable = (double *)malloc(numData*sizeof(double));
            para[0] = atof(line[0].c_str());
            para[1] = atof(line[1].c_str());
            para[2] = atof(line[2].c_str());
            para[3] = atof(line[3].c_str());
            para[4] = atof(line[4].c_str());
        }else{
            (*randomForceTable)[i++] = atof(line[0].c_str());
            (*randomForceTable)[i++] = atof(line[1].c_str());
            (*randomForceTable)[i++] = atof(line[2].c_str());
            (*randomForceTable)[i++] = atof(line[3].c_str());
        }
    }
    infile.close();

    return(numData);
}

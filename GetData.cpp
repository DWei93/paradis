#include <cstring>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <sstream>

using namespace std;

typedef struct {
   double   x,y;
   int      burgID;
} Data_t;


typedef struct {
   double   x,y;
} Point_t;


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

int main(int argc,char** argv)
{
    ifstream infile;
    stringstream strStream;
    string str, Avef, P1f, P2f, RDf, arg2;
    char filename[100], c1[100], c2[100], c3[100];
    int lineNo=0, i, j, k;
    double dis, rd;
    double aveY1=0.0, aveY2=0.0;
    double    n1=0, n2=0;
    Point_t point;
    
    double length = 20020;
    double seg = 20;

    if(argc < 2)return(1);

    int    nSegs = (int)(length/seg+1);
    vector<Point_t> p1, p2;

    vector<vector<double> > Rd1(nSegs), Rd2(nSegs), Rd(nSegs);
    for(i=0; i<nSegs; i++){
        Rd1[i].resize(1);
        Rd1[i][0] = ((double)i)*seg + seg*0.5;
        Rd2[i].resize(1);
        Rd2[i][0] = ((double)i)*seg + seg*0.5;
        Rd[i].resize(3);
        Rd[i][0] = ((double)i)*seg + seg*0.5;
    }

    Data_t          d;
    vector<Data_t>  data;
    vector<string> line;

    strcpy(filename,argv[1]);
    arg2 = argv[2];

    infile.open(filename);
    if(!infile) cout<<"error"<<endl;

    point.x = -100000;
    point.y = -100000;
    while(getline(infile,str))
    {   
        lineNo++;
        if(lineNo < 3)continue;
        
        line = split(str, " ");
        if(line.size() < 11)continue;

        d.x = atof(line[0].c_str());
        d.y = atof(line[1].c_str());
        d.burgID = atoi(line[10].c_str());
        data.push_back(d);

        if(point.x == d.x && point.y ==d.y)continue;

        point.x = d.x;
        point.y = d.y;

        if(d.burgID == 12){
            aveY1 += d.y;
            n1 += 1.0;
            p1.push_back(point);
        }else{
            aveY2 += d.y;
            n2 += 1.0;
            p2.push_back(point);
        }
    }

    Avef = arg2 + "Ave.plt";
    ofstream Ave(Avef);
    gcvt(aveY1/n1,10,c1);
    gcvt(aveY2/n2,10,c2);
    Ave <<  c1 << " " << c2 << endl; 
    Ave.close();

    
    P1f = arg2 + "P1.plt";
    ofstream P1(P1f);
    for(i=0; i<p1.size(); i++){
        gcvt(p1[i].x,15,c1);
        gcvt(p1[i].y,15,c2);
        P1 << c1 << " " << c2 << endl;
    }
    P1.close();

    P2f = arg2 + "P2.plt";
    ofstream P2(P2f);
    for(i=0; i<p2.size(); i++){
        gcvt(p2[i].x,15,c1);
        gcvt(p2[i].y,15,c2);
        P2 << c1 << " " << c2 << endl;
    }
    P2.close();

    for(i=0; i<data.size(); i++){
        for(j=0; j<data.size(); j++){
            dis = data[i].x - data[j].x;
            if(dis <= 0.0 || data[i].burgID != data[j].burgID)continue;
            
            k = (int)(dis/seg);
            if(k>Rd2.size())k=Rd2.size()-1;

            rd = fabs(data[i].y - data[j].y);
            if(data[i].burgID == 12){
                Rd1[k].push_back(rd); 
            }else{
                Rd2[k].push_back(rd); 
            }
        }
    } 

    RDf = arg2 + "Rd.plt";
    ofstream RD(RDf);
    for(i=0; i<Rd1.size(); i++){
        Rd[i][1] = 0.0;
        Rd[i][2] = 0.0;
        for(j=1; j<Rd1[i].size(); j++){
            Rd[i][1] += (Rd1[i][j]/(double)(Rd1[i].size()));
        }
        for(j=1; j<Rd2[i].size(); j++){
            Rd[i][2] += (Rd2[i][j]/(double)(Rd2[i].size()));
        }
        gcvt(Rd[i][0],15,c1);
        gcvt(Rd[i][1],15,c2);
        gcvt(Rd[i][2],15,c3);
        RD << c1 << " " << c2 << " "<< c3 << endl;
    }

    infile.close();
    RD.close();
    return 0;
}

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

typedef struct {
   double   x,y;
   int      burgID;
} Data_t;


typedef struct {
   double   x,y;
} Point_t;


bool cmp(const Point_t &p, const Point_t &q)
{
    return(p.x < q.x);
}

void FindIntersect(Point_t &p0, Point_t &pout, vector<Point_t> &list, int &ID, double length)
{
    int         i,j;
    Point_t     p1, p2;
    if(list.size()<1){
        pout.x = p0.x;
        pout.y = 0.0;
    }

    i = (ID < list.size()) ? ID : 0;
    for( ; i<list.size(); i++){
        if(i<list.size()-1){
            j = i+1;
            p2.x = list[j].x;
            p2.y = list[j].y;
        }else{
            j=0;
            p2.x = list[j].x+length;
            p2.y = list[j].y;
        }

        p1.x = list[i].x;
        p1.y = list[i].y;

        if(ID==0){
            if(p0.x < p1.x && p0.x < p2.x){
                p1.x = list[list.size()-1].x-length;
                p1.y = list[list.size()-1].y;

                p2.x = list[0].x;
                p2.y = list[0].y;
            }            
        }

        if(p0.x >= p1.x && p0.x <= p2.x){
            pout.x = p0.x;
            if(p2.x>p1.x){
                pout.y = p1.y + (p2.y-p1.y)*(p0.x-p1.x)/(p2.x-p1.x);
            }else{
                pout.y = 0.5*(p1.y+p2.y);
            }
        }else{
            ID++;
        }
       
    }
    
}
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
    char filename[100], c1[100], c2[100], c3[100], c4[100];
    int lineNo=0, i, j, k;
    double dis, rd;
    double aveY1=0.0, aveY2=0.0;
    double    n1=0, n2=0, area = 0.0;
    Point_t point, pout;
    
    double length = 20000;
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

        if(d.x < -0.5*length)d.x += length;
        if(d.x > 0.5*length)d.x -= length;

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

    sort(p1.begin(), p1.end(), cmp); 
    sort(p2.begin(), p2.end(), cmp); 

    vector<Point_t> p3(p2.size());
    int ID = 0;
    for(i=0; i<p2.size(); i++){
        point.x = p2[i].x;
        point.y = p2[i].y;
        FindIntersect(point, pout, p1, ID, length); 
        p3[i].x = pout.x;
        p3[i].y = pout.y;
    }

    
#if 0
    P1f = arg2 + "P1.plt";
    ofstream P1(P1f.c_str());
    for(i=0; i<p1.size(); i++){
        gcvt(p1[i].x,15,c1);
        gcvt(p1[i].y,15,c2);
        P1 << c1 << " " << c2 << endl;
    }
    P1.close();
#endif

    P2f = arg2 + "P.plt";
    ofstream P2(P2f.c_str());
    cout << "size of P: " << p2.size() << endl;
    for(i=0; i<p2.size(); i++){
        gcvt(p2[i].x,15,c1);
        gcvt(p2[i].y,15,c2);
        gcvt(p3[i].y,15,c3);
        gcvt(fabs(p3[i].y-p2[i].y),15,c4);
        P2 << c1 << " " << c2 << " " << c3 << " " << c4 <<endl;
        
        if(i == p2.size()-1)continue;
        area += ((p2[i+1].x-p2[i].x)*(p3[i+1].y-p2[i+1].y+p3[i].y-p2[i].y)*0.5);
    }
    P2.close();

    Avef = arg2 + "Ave.plt";
    ofstream Ave(Avef.c_str());
    gcvt(aveY1/n1,15,c1);
    gcvt(aveY2/n2,15,c2);
    gcvt(fabs(area)/length,15,c3);
    Ave << "average y of line 1: " << c1 << endl;
    Ave << "average y of line 2: " << c2 << endl;
    Ave << "average distance: " << c3 <<endl;
    Ave.close();

    for(i=0; i<data.size(); i++){
        for(j=0; j<data.size(); j++){
            dis = data[i].x - data[j].x;
            if(dis <= 0.0 || data[i].burgID != data[j].burgID)continue;
            
            k = (int)(dis/seg);
            if(k>=Rd2.size())k=0;
            if(k<0)k=Rd2.size()-1;

            rd = fabs(data[i].y - data[j].y);
            if(data[i].burgID == 12){
                Rd1[k].push_back(rd); 
            }else{
                Rd2[k].push_back(rd); 
            }
        }
    } 

    RDf = arg2 + "Rd.plt";
    ofstream RD(RDf.c_str());
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

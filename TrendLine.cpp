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
} Point_t;

bool cmp(const Point_t &p, const Point_t &q)
{
    return(p.x < q.x);
}

void RefineCurve(vector<Point_t>& curve, int nPoints){
    int         i, n, j;
    Point_t     point;
    double      x0, meshSize;

    n = curve.size()-1;
    for(i=0; i<n; i++){
        if(fabs(curve[i+1].x - curve[i].x) < 1.0E-3)continue;
        meshSize = (curve[i+1].x - curve[i].x)/(double)(nPoints+1);
        for(j=0; j<nPoints; j++){
            point.x = curve[i].x + (double)(j+1)*meshSize;
            point.y = (curve[i+1].y-curve[i].y)*(point.x-curve[i].x)/(curve[i+1].x - curve[i].x) + curve[i].y;
            curve.push_back(point);    
        }
    }
    sort(curve.begin(), curve.end(), cmp); 
    return;
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
    if(argc < 3)return(1);

    int     i, j, ID;
    int     xLine, yLine;
    string  str;
    double  length = 4000.0, aveY, nY;
    vector<string> line;
    char    c[100];
    
    xLine = atoi(argv[1]);
    yLine = atoi(argv[2]);

    ifstream                infile;

    Point_t                 point, pout;
    vector<Point_t>         curve;
    vector<vector<Point_t> > curveList;


    for(i=0; i<argc-3; i++){
        infile.open(argv[i+3]);
        if(!infile)return(1);

        vector<Point_t>().swap(curve);
        while(getline(infile,str)){
            line = split(str, " ");
            if(line.size()< yLine+1)continue;

            point.x = atof(line[xLine].c_str());
            point.y = atof(line[yLine].c_str());
            curve.push_back(point);
        }
        if(i<argc-4){
            curveList.push_back(curve);
        }
        infile.close();
    }
    
    RefineCurve(curve, 5);
    vector<vector<double> >      outArray(curve.size());
    for(i=0; i<outArray.size(); i++){
        outArray[i].push_back(curve[i].x);
        outArray[i].push_back(curve[i].y);
    }

    for(i=0; i<curveList.size(); i++){
        ID = 0;
        for(j=0; j<curve.size(); j++){
            FindIntersect(curve[j], pout, curveList[i], ID, length);
            outArray[j].push_back(pout.y);
        }
    }

    for(i=0; i<outArray.size(); i++){
        aveY = 0.0;
        nY = 0.0;
        for(j=1; j<outArray[i].size(); j++){
            aveY += outArray[i][j];
            nY += 1.0;
        }
        outArray[i].push_back(aveY/nY);
    }

    ofstream outfile("trendLine.plt");
    for(i=0; i<outArray.size(); i++){
        for(j=0; j<outArray[i].size(); j++){
            gcvt(outArray[i][j],15,c);
            outfile << c << " "; 
        }
        outfile <<endl;
    }
    outfile.close();
    return 0;
}

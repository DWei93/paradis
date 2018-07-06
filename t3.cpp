#include <stdio.h>
#include <stdlib.h>
#include <vector>

using namespace std;

void SumOfValuesIsConstant(int n, int sum, vector<int>& valueList, int tot, 
                            vector<vector<int> >& valueLists, int maxVar)
{
    int i, j, nextSum;

    if(n == 0)return;
    for(i=0; i<= maxVar && i<= sum; i++){
        valueList[n-1] = i;
        nextSum = sum - i;
        if(nextSum > (n-1)*maxVar)continue;
        SumOfValuesIsConstant(n-1, nextSum, valueList, tot
                                 ,valueLists, maxVar);
    }

    if(n-1 == 0 && nextSum <= maxVar){
        valueLists.push_back(valueList);
    }
    return;
}

int main(){
    int i,j, n, sum;
    

    n = 9;
    vector<int>         valueList(n);
    vector<vector<int> > valueLists;


    for(sum=0; sum<3*n+1; sum++){
        printf("sum %d\n", sum);
        SumOfValuesIsConstant(n, sum, valueList, n, valueLists, 3);
        printf("\n");
    }

    for(i=0; i<valueLists.size(); i++){
        printf("\ni = %d:", i);
        for(j=0; j<valueLists[i].size(); j++){
            printf("%-2d ", valueLists[i][j]);
        }
    }

   vector<vector<int> > vec;

    vector<int> a;
    a.push_back(1);
    a.push_back(2);
    a.push_back(3);
    a.push_back(3);

    vector<int> b;
    b.push_back(4);
    b.push_back(5);
    b.push_back(6);

    vec.push_back(a);
    vec.push_back(b); 
    for(i=0; i<vec.size(); i++){
        printf("i = %d: ", i);
        for(j=0; j<vec[i].size(); j++){
            printf("%d ",vec[i][j]);
        }
        printf("\n");
    }

    return(0);
}





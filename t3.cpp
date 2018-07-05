#include <stdio.h>
#include <stdlib.h>
#include <vector>

using namespace std;

void SumOfValuesIsConstant(int n, int sum, vector<int> valueList, int tot, 
                            vector<vector<int>> valueLists, int maxV = 3)
{
    int i, j, nextSum;

    if(n == 0)return;
    for(i=0; i<= maxV && i<= sum; i++){
        valueList[n-1] = i;
        nextSum = sum - i;
        if(nextSum > (n-1)*maxV)continue;
        SumOfValuesIsConstant(n-1, nextSum, valueList, tot
                                 ,valueLists);
    }

    if(n-1 == 0 && nextSum <= 3){
        for(j=0; j<tot; j++){
            printf("%2d ", valueList[j]);
        }
        printf("\n");
        valueLists.push_back(valueList);
    }
    return;
}

int main(){
    int i, n, sum;
    

    n = 3;
    vector<int>         valueList(n);
    vector<vector<int>> valueLists;


    for(sum=0; sum<3*n+1; sum++){
        printf("sum %d\n", sum);
        SumOfValuesIsConstant(n, sum, valueList, n, valueLists);
        printf("\n");
    }

    
    return(0);
}

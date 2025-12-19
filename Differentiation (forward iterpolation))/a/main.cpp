#include <bits/stdc++.h>

using namespace std;

vector<vector<double>>getminor(int r, int c,vector<vector<double>>&A ){
int n =A.size();
vector<vector<double>>minor;
for(int i=0;i<n;i++){
    if(i==r)continue;
    vector<double>row;
    for(int j =0;j<n;j++){
    if(j==c)continue;
        row.push_back(A[i][j]);
    }
    minor.push_back(row);
}

return minor;
}
double det(vector<vector<double>>&A){
int n =A.size();
if(n==1){
    return A[0][0];
}
if(n==2){
    return A[0][0]*A[1][1]-A[0][1]*A[1][0];
}
double deter=0;
for(int col=0;col<n;col++){
    vector<vector<double>>m=getminor(0,col,A);
    if(col %2==0)deter+=A[0][col]*det(m);
    else deter-=A[0][col]*det(m);
}

return deter;
}
vector<vector<double>>cofactor(vector<vector<double>>&A){
int n =A.size();vector<vector<double>>c(n,vector<double>(n));

for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
        vector<vector<double>>m=getminor(i,j,A);
        if((i+j)%2==0)c[i][j]=det(m);
        else c[i][j]=-det(m);
    }
}
return c;

}
vector<vector<double>>Transpose(vector<vector<double>>&A){
int n =A.size();
vector<vector<double>>T(n,vector<double>(n));
for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
        T[j][i]=A[i][j];
    }
}
return T;
}
vector<vector<double>>Inverse(vector<vector<double>>&A){
int n=A.size();
double detA=det(A);
if(detA==0)cout<<"inverse dont exist"<<endl;

vector<vector<double>>C=cofactor(A);
vector<vector<double>>adj=Transpose(C);
vector<vector<double>>inv(n,vector<double>(n));
for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
        inv[i][j]=adj[i][j]/detA;
    }
}
return inv;
}

vector<double>mul(vector<vector<double>>&A,vector<double>&B){
int n=A.size();
vector<double>X(n);
for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
     X[i]+=A[i][j]*B[j];
    }

}
return X;

}
void rankofM(vector<vector<double>>A){
int n =A.size();

for(int i=0;i<n-1;i++){
    int p=i;
    for(int j=i+1;j<n;j++){
        if(fabs(A[j][i]>A[p][i]))p=j;
        if(p!=i)swap(A[i],A[p]);

        if(A[i][i]==0)continue;
    }
    for(int j=i+1;j<n;j++){
        double f=A[j][i]/A[i][i];
        for(int k=i;k<=n;k++){
            A[j][k]-=f*A[i][k];
        }

    }

}
int rankA=0;
int rankAug=0;
for(int i=0;i<n;i++){
        bool nz=false;
    for(int j=0;j<n;j++){
        if(A[i][j]!=0){
            nz=true;
        }

    }
if(nz)rankA++;

}
for(int i=0;i<n;i++){
        bool nz=false;
    for(int j=0;j<n+1;j++){
        if(A[i][j]!=0){
            nz=true;
        }

    }
if(nz)rankAug++;

}
cout<<"rank of A matrix:  "<<rankA<<endl;
cout<<"rank of Augmented matrix:  "<<rankAug<<endl;
if(rankA<rankAug)cout<<" no solution exists "<<endl;
if(rankA==rankAug && rankA <n )cout<<" infinite solution exists "<<endl;
else if(rankA==rankAug && rankA ==n)cout<<" unique solution exists "<<endl;

}

int main()
{
   int n;cout<<"enter no eqn" <<endl;cin>>n;
    vector<vector<double>>A(n,vector<double>(n));
     vector<double>B(n);

    for(int i =0;i<n;i++){
         for(int j =0;j<n;j++){
            cin>>A[i][j];
         }
    }
      for(int j =0;j<n;j++){
            cin>>B[j];
         }
        vector<vector<double>>Aug(n,vector<double>(n+1));
         for(int i =0;i<n;i++){
         for(int j =0;j<n;j++){
            Aug[i][j]=A[i][j];
         }
         Aug[i][n]=B[i];
    }
    rankofM(Aug);

    vector<vector<double>>inverseM=Inverse(A);
    vector<double>Multiply=mul(inverseM,B);
    cout<<"Solution vector printed --->>>>> "<<endl;
     for(int j =0;j<n;j++){
            cout<<Multiply[j]<<" ";
         }
         cout<<endl;

}

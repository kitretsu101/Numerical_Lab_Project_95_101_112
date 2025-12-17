#include <bits/stdc++.h>
using namespace std;

int type(vector<vector<double>>&a,int n)
{
    int r=0,c=0;
    for(int i=0;i<n;i++){
        bool nz=false;
        for(int j=0;j<n;j++)
            if(fabs(a[i][j])>1e-9) nz=true;
        if(nz) r++;
        else if(fabs(a[i][n])>1e-9) return 0;
    }
    for(int j=0;j<n;j++){
        for(int i=0;i<n;i++)
            if(fabs(a[i][j])>1e-9){
                c++;
                break;
            }
    }
    if(r==c) return 1;
    return 2;
}

int main()
{
    ifstream fin("input.txt");
    ofstream fout("output.txt");
    if(!fin){
        cout<<"input.txt not found"<<endl;
        return 0;
    }

    while(true)
    {
        int n;
        fin>>n;
        vector<vector<double>>a(n,vector<double>(n+1));
        for(int i=0;i<n;i++)
            for(int j=0;j<=n;j++)
                fin>>a[i][j];

        for(int k=0;k<n;k++){
            int mx=k;
            for(int i=k+1;i<n;i++)
                if(fabs(a[i][k])>fabs(a[mx][k])) mx=i;
            swap(a[k],a[mx]);

            if(fabs(a[k][k])<1e-9) continue;

            for(int i=k+1;i<n;i++){
                double f=a[i][k]/a[k][k];
                for(int j=k;j<=n;j++)
                    a[i][j]-=f*a[k][j];
            }
        }

        int t=type(a,n);

        fout<<"System size: "<<n<<endl;

        if(t==0){
            fout<<"Solution type: No solution"<<endl;
        }
        else if(t==2){
            fout<<"Solution type: Infinite solutions"<<endl;
        }
        else{
            fout<<"Solution type: Unique solution"<<endl;
            vector<double>x(n);
            for(int i=n-1;i>=0;i--){
                x[i]=a[i][n];
                for(int j=i+1;j<n;j++)
                    x[i]-=a[i][j]*x[j];
                x[i]/=a[i][i];
            }
            for(int i=0;i<n;i++)
                fout<<"x"<<i+1<<" = "<<x[i]<<endl;
        }

        char ch;
        fin>>ch;
        if(ch!='y' && ch!='Y') break;
        fout<<endl;
    }

    fin.close();
    fout.close();
    return 0;
}

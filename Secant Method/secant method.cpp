#include <bits/stdc++.h>
using namespace std;

double poly(const vector<double>&c,double x){
    double r=0;
    int n=c.size();
    for(int i=0;i<n;i++) r+=c[i]*pow(x,n-1-i);
    return r;
}

void printPoly(const vector<double>&c,ostream &out){
    int n=c.size();
    for(int i=0;i<n;i++){
        if(c[i]==0) continue;
        if(i>0&&c[i]>0) out<<" + ";
        if(c[i]<0) out<<" - ";
        if(abs(c[i])!=1||i==n-1) out<<abs(c[i]);
        if(i<n-1){
            out<<"x";
            if(n-1-i>1) out<<"^"<<n-1-i;
        }
    }
    out<<endl;
}

int main(){
    ifstream fin("input.txt");
    ofstream fout("output.txt");
    if(!fin) return 0;

    int degree;
    fin>>degree;
    vector<double>coef(degree+1);
    for(int i=0;i<=degree;i++) fin>>coef[i];

    double minX,maxX,step,tol;
    fin>>minX>>maxX>>step>>tol;

    fout<<"Degree of Polynomial: "<<degree<<endl;
    fout<<"Polynomial: ";
    printPoly(coef,fout);

    vector<double>roots;
    double x=minX;

    while(x<maxX){
        double x0=x,x1=x+step;
        if(x1>maxX) x1=maxX;

        double f0=poly(coef,x0),f1=poly(coef,x1);

        if(f0*f1<=0){
            double curr=x0,prev=x1;
            double fcurr=poly(coef,curr),fprev=poly(coef,prev);
            double next;

            do{
                if(fcurr==fprev) break;
                next=curr-fcurr*(curr-prev)/(fcurr-fprev);
                prev=curr;
                fprev=fcurr;
                curr=next;
                fcurr=poly(coef,curr);
            }while(fabs(fcurr)>tol);

            bool exists=false;
            for(double r:roots) if(fabs(r-curr)<tol) exists=true;
            if(!exists) roots.push_back(curr);
        }

        x+=step;
    }

    if(roots.empty()) fout<<"No roots found in given range"<<endl;
    else{
        fout<<"Roots found:"<<endl;
        for(int i=0;i<roots.size();i++) fout<<"Root "<<i+1<<": "<<roots[i]<<endl;
    }

    fin.close();
    fout.close();
    return 0;
}

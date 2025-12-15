#include <bits/stdc++.h>
using namespace std;

double poly(const vector<double>&c,double x)
{
    double s=0;
    int n=c.size();
    for(int i=0;i<n;i++)s+=c[i]*pow(x,n-1-i);
    return s;
}

vector<double> derivative(const vector<double>&c)
{
    int n=c.size();
    vector<double>d;
    for(int i=0;i<n-1;i++)
        d.push_back(c[i]*(n-1-i));
    return d;
}

void printPoly(const vector<double>& c,ostream& out)
{
    int n=c.size();
    for(int i=0;i<n;i++)
    {
        if(c[i]==0)continue;
        if(i>0 && c[i]>0)out<<" + ";
        if(c[i]<0)out<<" - ";
        if(abs(c[i])!=1 || i==n-1)out<<abs(c[i]);
        if(i<n-1)
        {
            out<<"x";
            if(n-1-i>1)out<<"^"<<n-1-i;
        }
    }
    out<<endl;
}

int main()
{
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    if(!fin)
    {
        cout<<"Input file not found"<<endl;
        return 0;
    }

    int degree;
    fin>>degree;

    vector<double>coef(degree+1);
    for(int i=0;i<=degree;i++)fin>>coef[i];

    double a,b,step,tol;
    fin>>a>>b>>step>>tol;

    vector<double>dcoef=derivative(coef);
    vector<double>roots;

    for(double x=a;x<=b;x+=step)
    {
        double x0=x;

        for(int i=0;i<100;i++)
        {
            double f=poly(coef,x0);
            double df=poly(dcoef,x0);

            if(fabs(df)<1e-8)break;

            double x1=x0-f/df;

            if(fabs(x1-x0)<tol)
            {
                if(x1>=a && x1<=b)
                {
                    bool exists=false;
                    for(double r:roots)
                        if(fabs(r-x1)<tol) exists=true;

                    if(!exists) roots.push_back(x1);
                }
                break;
            }
            x0=x1;
        }
    }

    cout<<"Degree of Polynomial: "<<degree<<endl;
    fout<<"Degree of Polynomial: "<<degree<<endl;

    cout<<"Polynomial: ";
    fout<<"Polynomial: ";
    printPoly(coef,cout);
    printPoly(coef,fout);

    cout<<"Derivative: ";
    fout<<"Derivative: ";
    printPoly(dcoef,cout);
    printPoly(dcoef,fout);

    if(roots.empty())
    {
        cout<<"No roots found in given range"<<endl;
        fout<<"No roots found in given range"<<endl;
    }
    else
    {
        cout<<"Roots found:"<<endl;
        fout<<"Roots found:"<<endl;
        for(int i=0;i<roots.size();i++)
        {
            cout<<"Root "<<i+1<<": "<<roots[i]<<endl;
            fout<<"Root "<<i+1<<": "<<roots[i]<<endl;
        }
    }

    fin.close();
    fout.close();
    return 0;
}

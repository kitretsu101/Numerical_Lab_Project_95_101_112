#include <bits/stdc++.h>
using namespace std;

float dydx(float x,float y)
{
    return (x-y)/2;
}

int main()
{
    ifstream fin("input.txt");
    ofstream fout("output.txt");
    if(!fin)
    {
        cout<<"input.txt not found"<<endl;
        return 0;
    }

    float x0,y0,x,h;
    fin>>x0>>y0>>x>>h;
    float X0=x0;

    cout<<"Input Read From File:"<<endl;
    cout<<"x0 = "<<x0<<endl;
    cout<<"y0 = "<<y0<<endl;
    cout<<"x  = "<<x<<endl;
    cout<<"h  = "<<h<<endl;

    int n=(int)((x-x0)/h);
    float y=y0;
    for(int i=1;i<=n;i++)
    {
        float k1=h*dydx(x0,y);
        float k2=h*dydx(x0+0.5*h,y+0.5*k1);
        float k3=h*dydx(x0+0.5*h,y+0.5*k2);
        float k4=h*dydx(x0+h,y+k3);
        y=y+(k1+2*k2+2*k3+k4)/6.0;
        x0+=h;
    }

    cout<<"The value of y at x is: "<<y<<endl;
    fout<<"Initial x0: "<<X0<<endl;
    fout<<"Initial y0: "<<y0<<endl;
    fout<<"Final x: "<<x<<endl;
    fout<<"Step h: "<<h<<endl;
    fout<<"The value of y at x is: "<<y<<endl;

    fin.close();
    fout.close();
    return 0;
}

#include <iostream>
#include <bits/stdc++.h>
#include <fstream>
#include <math.h>

using namespace std;

int main()
{
    int n;
    //cin >> n;

    ifstream infile("input.txt");
    ofstream outfile("output.txt");
    double sumx = 0;
    double sumy = 0;
    double sumxy = 0;
    double sumx2 = 0;

    while(infile >> n)
    {
        for(int i = 0; i < n; i++)
        {
            double x, y;
            //cin >> x >> y;
            infile >> x >> y;
            sumx += log(x);
            sumy += log(y);
            sumxy += log(x) * log(y);
            sumx2 += log(x) * log(x);
        }

        double b = (n*sumxy - sumx*sumy) / (n*sumx2 - sumx*sumx);
        double a = (sumy - b*sumx) / n;
        double A = exp(a);

        ostringstream os;
        os << "y = " << A << "x^" << b << endl;
        cout << os.str();
        outfile << os.str();
    }

    outfile.close();
    infile.close();
}

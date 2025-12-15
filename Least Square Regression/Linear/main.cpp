#include <iostream>
#include <bits/stdc++.h>
#include <fstream>

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
            sumx += x;
            sumy += y;
            sumxy += x * y;
            sumx2 += x * x;
        }

        double b = (n*sumxy - sumx*sumy) / (n*sumx2 - sumx*sumx);
        double a = (sumy - b*sumx) / n;

        ostringstream os;
        os << "y = " << a << " + " << b << "x" << endl;
        cout << os.str();
        outfile << os.str();
    }

    outfile.close();
    infile.close();
}

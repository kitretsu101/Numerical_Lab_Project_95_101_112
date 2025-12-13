#include <iostream>
#include <bits/stdc++.h>
#include <fstream>

using namespace std;

int main()
{
    ifstream infile("input.txt");
    ofstream outfile("output.txt");


    int n;
    //cin >> n;
    while(infile >> n) {
    vector<vector<double>> y(n, vector<double>(n, 0));
    vector<double> x(n, 0);

    for(int i = 0; i < n; i++) {
        infile >> x[i] >> y[i][0];
    }

    for(int j = 1; j < n; j++) {
        for(int i = 0; i < n-j; i++) {
            y[i][j] = y[i+1][j-1] - y[i][j-1];
        }
    }

    double interpolate_val;
    infile >> interpolate_val;

    double u = (interpolate_val - x[0]) / (x[1] - x[0]);
    double sum = y[0][0];
    double u_sum = 1;
    for(int i = 0; i < n-1; i++) {
        u_sum *= (u - i) / (i+1);
        sum += u_sum * y[0][i+1];
    }

    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            cout << y[i][j] << " ";
        }
        cout << endl;
    }

    stringstream os;
    os << "value of y at " << interpolate_val << " is : " << sum << endl;
    cout << os.str();
    outfile << os.str();
    }
}

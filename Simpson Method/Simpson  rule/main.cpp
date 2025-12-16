#include <bits/stdc++.h>
using namespace std;

double f(double x)
{
    double sin_x = sin(x);
    return pow(sin_x, 5) + 4.0 * pow(sin_x, 4) + 1.0;
}

double simp_one_third(double a, double b, int N)
{
    if (N % 2 != 0)
    {
        return NAN;
    }

    double h = (b - a) / N;
    double integral = f(a) + f(b);

    for (int i = 1; i < N; ++i)
    {
        double x = a + i * h;
        if (i % 2 == 0)
        {
            integral += 2.0 * f(x);
        }
        else
        {
            integral += 4.0 * f(x);
        }
    }

    return (h / 3.0) * integral;
}

double simpsons_three_eighths(double a, double b, int N)
{
    if (N % 3 != 0)
    {
        return NAN;
    }

    double h = (b - a) / N;
    double integral = f(a) + f(b);

    for (int i = 1; i < N; ++i)
    {
        double x = a + i * h;
        if (i % 3 == 0)
        {
            integral += 2.0 * f(x);
        }
        else
        {
            integral += 3.0 * f(x);
        }
    }

    return (3.0 * h / 8.0) * integral;
}

void run_simpsons_rules()
{
    double a, b;
    int N;

    ifstream inputFile("input.txt");
    if (!inputFile.is_open())
    {
        cerr << "Error: Could not open input.txt" << endl;
        return;
    }

    if (!(inputFile >> a >> b >> N))
    {
        cerr << "Error: Failed to read a, b, and N from input.txt. Check the file format." << endl;
        inputFile.close();
        return;
    }
    inputFile.close();

    double integral_1_3 = simp_one_third(a, b, N);
    double integral_3_8 = simpsons_three_eighths(a, b, N);

    ofstream outputFile("output.txt");
    if (!outputFile.is_open())
    {
        cerr << "Error: Could not open output.txt for writing." << endl;
        return;
    }

    outputFile << fixed << setprecision(10);

    outputFile << "--- Numerical Integration Results (Simpson's Rules) ---" << endl;
    outputFile << "Function: f(x) = sin^5(x) + 4*sin^4(x) + 1" << endl;
    outputFile << "Limits: a = " << a << ", b = " << b << endl;
    outputFile << "Number of Intervals (N) = " << N << endl;
    outputFile << "Interval width (h) = " << (b - a) / N << endl;
    outputFile << "--------------------------------------------------------" << endl;

    if (N % 2 == 0)
    {
        outputFile << "Simpson's 1/3 Rule (N is even):" << endl;
        outputFile << "Integral Estimate: " << integral_1_3 << endl;
    }
    else
    {
        outputFile << "Simpson's 1/3 Rule: N=" << N << " is not even. Calculation skipped." << endl;
    }

    outputFile << "--------------------------------------------------------" << endl;

    if (N % 3 == 0)
    {
        outputFile << "Simpson's 3/8 Rule (N is multiple of 3):" << endl;
        outputFile << "Integral Estimate: " << integral_3_8 << endl;
    }
    else
    {
        outputFile << "Simpson's 3/8 Rule: N=" << N << " is not a multiple of 3. Calculation skipped." << endl;
    }

    outputFile.close();

    cout << "Calculation complete. Results written to output.txt." << endl;

    cout << "\n--- Generated output.txt Content ---" << endl;
    ifstream mockOutput("output.txt");
    string line;
    while (getline(mockOutput, line))
    {
        cout << line << endl;
    }
    mockOutput.close();
}

int main()
{
    cout << fixed << setprecision(10);
    run_simpsons_rules();
    return 0;
}
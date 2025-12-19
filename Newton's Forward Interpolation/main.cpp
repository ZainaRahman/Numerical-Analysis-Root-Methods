#include <bits/stdc++.h>
using namespace std;


int fact(int n)
{
    int f = 1;
    for(int i = 2; i <= n; i++)
        f *= i;
    return f;
}


int detectOrder( vector<vector<double>>& diff, double eps = 1e-9)
{
    int n = diff.size();

    for(int j = 1; j < n; j++)
    {
        bool allZero = true;
        for(int i = 0; i < n - j; i++)
        {
            if (fabs(diff[i][j]) > eps)
            {
                allZero = false;
                break;
            }
        }
        if (allZero)
            return j - 1;
    }
    return n - 1;
}

int main()
{
    ifstream inputFile("Input.txt");
    ofstream outputFile("Output.txt");

    if (!inputFile.is_open())
    {
        cout << "Error opening Input.txt\n";
        return 1;
    }

    int testCases;
    inputFile >> testCases;

    for (int tc = 1; tc <= testCases; tc++)
    {
        int n;
        inputFile >> n;

        vector<double> x(n), y(n);
        for (int i = 0; i < n; i++) inputFile >> x[i];
        for (int i = 0; i < n; i++) inputFile >> y[i];

        double X;
        inputFile >> X;

        cout << "Test Case " << tc << "\n";
        outputFile << "Test Case " << tc << "\n";

        double h = x[1] - x[0];
        double u = (X - x[0]) / h;


        vector<vector<double>> diff(n, vector<double>(n, 0.0));
        for (int i = 0; i < n; i++)
            diff[i][0] = y[i];

        for (int j = 1; j < n; j++)
            for (int i = 0; i < n - j; i++)
                diff[i][j] = diff[i + 1][j - 1] - diff[i][j - 1];

        cout << "Forward Difference Table:\n";
        outputFile << "Forward Difference Table:\n";

       cout << fixed << setprecision(2);
outputFile << fixed << setprecision(2);

for (int i = 0; i < n; i++)
{

    cout << setw(10) << diff[i][0];
    outputFile << setw(10) << diff[i][0];


    for (int j = 1; j < n - i; j++)
    {
        cout << setw(10) << diff[i][j];
        outputFile << setw(10) << diff[i][j];
    }

    cout << "\n";
    outputFile << "\n";
}


        double fx = diff[0][0];
        double u_term = 1.0;

        for (int k = 1; k < n; k++)
        {
            u_term *= (u - (k - 1));
            fx += (u_term * diff[0][k]) / fact(k);
        }

        int order = detectOrder(diff);


        double error = 0.0;
        if (order + 1 < n)
        {
            double u_err = 1.0;
            for (int i = 0; i <= order; i++)
                u_err *= (u - i);

            error = fabs((u_err * diff[0][order + 1]) / fact(order + 1));
        }

        cout << "Step size (h) = " << h << "\n";
        cout << "u = " << u << "\n";
        cout << "Detected Polynomial Order = " << order << "\n";
        cout << "Interpolated value f(" << X << ") = "
             << fixed << setprecision(2) << fx << "\n";
        cout << "Estimated Forward Interpolation Error = "
             << fixed << setprecision(2) << error << "\n\n";

        outputFile << "Step size (h) = " << h << "\n";
        outputFile << "u = " << u << "\n";
        outputFile << "Detected Polynomial Order = " << order << "\n";
        outputFile << "Interpolated value f(" << X << ") = "
                   << fixed << setprecision(2) << fx << "\n";
        outputFile << "Estimated Forward Interpolation Error = "
                   << fixed << setprecision(2) << error << "\n\n";
    }

    inputFile.close();
    outputFile.close();
    return 0;
}

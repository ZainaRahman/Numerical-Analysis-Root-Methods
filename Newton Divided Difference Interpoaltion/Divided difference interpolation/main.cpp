#include <bits/stdc++.h>
using namespace std;

int main()
{
    ifstream inputFile("Input.txt");
    ofstream outputFile("Output.txt");

    if(!inputFile.is_open())
        return 1;

    int T;
    inputFile >> T;

    cout << fixed << setprecision(2);
    outputFile << fixed << setprecision(2);

    for(int tc = 1; tc <= T; tc++)
    {
        int n;
        inputFile >> n;

        vector<double> x(n), y(n);
        for(int i = 0; i < n; i++) inputFile >> x[i];
        for(int i = 0; i < n; i++) inputFile >> y[i];

        double x_extra, y_extra;
        inputFile >> x_extra >> y_extra;

        double X;
        inputFile >> X;

        vector<double> x_all = x;
        vector<double> y_all = y;
        x_all.push_back(x_extra);
        y_all.push_back(y_extra);

        vector<vector<double>> dd(n+1, vector<double>(n+1, 0.0));

        for(int i = 0; i <= n; i++)
            dd[i][0] = y_all[i];

        for(int j = 1; j <= n; j++)
            for(int i = 0; i <= n - j; i++)
                dd[i][j] = (dd[i+1][j-1] - dd[i][j-1]) / (x_all[i+j] - x_all[i]);

        double fx = dd[0][0];
        double term = 1.0;

        for(int i = 1; i < n; i++)
        {
            term *= (X - x[i-1]);
            fx += term * dd[0][i];
        }

        double prod = 1.0;
        for(int i = 0; i < n; i++)
            prod *= (X - x[i]);

        double error = fabs(dd[0][n] * prod);

        cout << "Test Case " << tc << "\n";
        cout << "Divided Difference Table:\n";

        for(int i = 0; i <= n; i++)
        {
            for(int j = 0; j <= n - i; j++)
                cout << setw(12) << dd[i][j];
            cout << "\n";
        }

        cout << "Interpolated value f(" << X << ") = " << fx << "\n";
        cout << "Truncation Error = " << error << "\n\n";

        outputFile << "Test Case " << tc << "\n";
        outputFile << "Interpolated value f(" << X << ") = " << fx << "\n";
        outputFile << "Truncation Error = " << error << "\n\n";
    }

    inputFile.close();
    outputFile.close();
    return 0;
}

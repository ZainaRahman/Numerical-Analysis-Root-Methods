#include <bits/stdc++.h>
using namespace std;

int main() {

    // ---- File Handling ----
    ifstream in("input.txt");
    ofstream out("output.txt");

    if (!in) {
        cerr << "Error: input.txt not found\n";
        return 1;
    }

    int n;
    in >> n;

    vector<double> x(n), y(n);
    for (int i = 0; i < n; i++) {
        in >> x[i] >> y[i];
    }

    // Take log of y
    vector<double> Y(n);
    for (int i = 0; i < n; i++) {
        Y[i] = log(y[i]);   // ln(y)
    }

    double sumx = 0, sumY = 0, sumxY = 0, sumx2 = 0;

    for (int i = 0; i < n; i++) {
        sumx  += x[i];
        sumY  += Y[i];
        sumxY += x[i] * Y[i];
        sumx2 += x[i] * x[i];
    }

    // Least squares for Y = A + b x
    double b = (n * sumxY - sumx * sumY) / (n * sumx2 - sumx * sumx);
    double A = (sumY - b * sumx) / n;
    double a = exp(A);

    // ---- Output ----
    out << "Transcendental (Exponential) Fit:\n";
    out << "y = " << a << " * e^(" << b << "x)\n";

    return 0;
}

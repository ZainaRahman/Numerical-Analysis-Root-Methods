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

    double sumx = 0, sumy = 0, sumxy = 0, sumx2 = 0;

    for (int i = 0; i < n; i++) {
        sumx  += x[i];
        sumy  += y[i];
        sumxy += x[i] * y[i];
        sumx2 += x[i] * x[i];
    }

    // Least squares coefficients
    double b = (n * sumxy - sumx * sumy) / (n * sumx2 - sumx * sumx);
    double a = (sumy - b * sumx) / n;

    // ---- Output ----
    out << "Linear Fit Equation:\n";
    out << "y = " << a << " + " << b << "x\n";

    return 0;
}

#include <bits/stdc++.h>
using namespace std;

// Solve nXn linear system using Gauss–Jordan elimination
vector<double> solveN(vector<vector<double>> A, vector<double> B) {

    int n = A.size();

    for (int i = 0; i < n; i++) {

        // Pivot selection
        int pivot = i;
        for (int j = i + 1; j < n; j++) {
            if (fabs(A[j][i]) > fabs(A[pivot][i])) {
                pivot = j;
            }
        }

        swap(A[i], A[pivot]);
        swap(B[i], B[pivot]);

        double div = A[i][i];
        if (fabs(div) < 1e-12) {
            cout << "Singular system detected. No unique solution." << endl;
            return {};   // return empty vector
        }

        // Normalize pivot row
        for (int j = 0; j < n; j++)
            A[i][j] /= div;
        B[i] /= div;

        // Eliminate other rows
        for (int j = 0; j < n; j++) {
            if (j != i) {
                double factor = A[j][i];
                for (int k = 0; k < n; k++)
                    A[j][k] -= factor * A[i][k];
                B[j] -= factor * B[i];
            }
        }
    }

    return B;
}


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
    for (int i = 0; i < n; i++)
        in >> x[i] >> y[i];

    // Required summations
    double sx = 0, sx2 = 0, sx3 = 0, sx4 = 0;
    double sy = 0, sxy = 0, sx2y = 0;

    for (int i = 0; i < n; i++) {
        sx   += x[i];
        sx2  += x[i] * x[i];
        sx3  += x[i] * x[i] * x[i];
        sx4  += x[i] * x[i] * x[i] * x[i];
        sy   += y[i];
        sxy  += x[i] * y[i];
        sx2y += x[i] * x[i] * y[i];
    }

    // Normal equations
    vector<vector<double>> A = {
        { double(n), sx,  sx2 },
        { sx,  sx2, sx3 },
        { sx2, sx3, sx4 }
    };

    vector<double> B = { sy, sxy, sx2y };

    vector<double> sol = solveN(A, B);

    if (sol.empty()) {
    cout << "Solution could not be computed.\n";
    return 0;
    }

    // ---- Output ----
    out << "Quadratic Polynomial Fit:\n";
    out << "y = " << sol[0]
        << " + " << sol[1] << "x"
        << " + " << sol[2] << "x^2\n";

    return 0;
}

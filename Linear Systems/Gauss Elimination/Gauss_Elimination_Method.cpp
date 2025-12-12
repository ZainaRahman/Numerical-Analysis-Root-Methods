#include <bits/stdc++.h>
using namespace std;

int main() {

    // File handling
    ifstream in("input.txt");
    ofstream out("output.txt");

    if (!in) {
        cerr << "Error: input.txt not found\n";
        return 1;
    }

    int n;
    in >> n;

    vector<vector<float>> a(n, vector<float>(n + 1));

    // Read augmented matrix
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= n; j++) {
            in >> a[i][j];
        }
    }

    // -------- Forward Elimination --------
    for (int i = 0; i < n - 1; i++) {
        for (int k = i + 1; k < n; k++) {

            if (a[i][i] == 0) {
                out << "Error: Zero pivot encountered.\n";
                return 1;
            }

            float factor = a[k][i] / a[i][i];

            for (int j = i; j <= n; j++) {
                a[k][j] -= factor * a[i][j];
            }
        }
    }

    // -------- Back Substitution --------
    vector<float> x(n);

    for (int i = n - 1; i >= 0; i--) {
        x[i] = a[i][n];

        for (int j = i + 1; j < n; j++) {
            x[i] -= a[i][j] * x[j];
        }

        if (a[i][i] == 0) {
            out << "Error: Division by zero (infinite or no solution).\n";
            return 1;
        }

        x[i] /= a[i][i];
    }

    // -------- Output --------
    out << "Solution:\n";
    for (int i = 0; i < n; i++) {
        out << "x" << i + 1 << " = " << x[i] << "\n";
    }

    return 0;
}

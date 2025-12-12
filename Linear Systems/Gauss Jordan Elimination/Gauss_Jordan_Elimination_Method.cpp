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

    vector<vector<float>> a(n, vector<float>(n + 1));

    // Read augmented matrix
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= n; j++) {
            in >> a[i][j];
        }
    }

    // Copy original matrix for printing echelon form
    vector<vector<float>> echelon = a;

    // ---- Forward Elimination for Echelon Form ----
    for (int i = 0; i < n; i++) {

        // Pivot row selection
        int maxRow = i;
        for (int k = i + 1; k < n; k++) {
            if (fabs(echelon[k][i]) > fabs(echelon[maxRow][i]))
                maxRow = k;
        }

        swap(echelon[i], echelon[maxRow]);

        if (fabs(echelon[i][i]) < 1e-9)
            continue; // Singular or zero pivot; skip elimination

        // Eliminate below pivot
        for (int k = i + 1; k < n; k++) {
            float factor = echelon[k][i] / echelon[i][i];
            for (int j = i; j <= n; j++) {
                echelon[k][j] -= factor * echelon[i][j];
            }
        }
    }

    // ---- Print Echelon Form ----
    out << "Echelon Form (Row-Reduced to Upper Triangular):\n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= n; j++)
            out << fixed << setprecision(4) << echelon[i][j] << " ";
        out << "\n";
    }
    out << "\n";

    // ---- Full Gauss–Jordan Elimination ----
    for (int i = 0; i < n; i++) {

        // Find pivot row again
        int maxRow = i;
        for (int k = i + 1; k < n; k++) {
            if (fabs(a[k][i]) > fabs(a[maxRow][i]))
                maxRow = k;
        }

        swap(a[i], a[maxRow]);

        if (fabs(a[i][i]) < 1e-9) {
            out << "Error: Zero pivot encountered. System may be dependent or inconsistent.\n";
            return 1;
        }

        // Normalize pivot row
        float pivot = a[i][i];
        for (int j = 0; j <= n; j++)
            a[i][j] /= pivot;

        // Eliminate all other rows
        for (int k = 0; k < n; k++) {
            if (k != i) {
                float factor = a[k][i];
                for (int j = 0; j <= n; j++)
                    a[k][j] -= factor * a[i][j];
            }
        }
    }

    // ---- Output Solution ----
    out << "Solution:\n";
    for (int i = 0; i < n; i++) {
        out << "x" << i + 1 << " = " << a[i][n] << "\n";
    }

    return 0;
}

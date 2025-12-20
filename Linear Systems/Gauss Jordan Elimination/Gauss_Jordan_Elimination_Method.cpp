#include <bits/stdc++.h>
using namespace std;

int main() {

    // File Handling
    ifstream in("input.txt");
    ofstream out("output.txt");

    if (!in) {
        cerr << "Error: input.txt not found\n";
        return 1;
    }

    int casing;
    in >> casing;
    int r=1;

    while (casing--){
    out << "----- Case " << r++ << " -----\n\n";
    int n;
    in >> n;

    vector<vector<float>> a(n, vector<float>(n + 1));
    cout<<endl;

    // Reading augmented matrix
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= n; j++) {
            in >> a[i][j];
        }
    }

    //Copy of matrix for echelon form
    vector<vector<float>> echelon = a;

    // Forward Elimination (Echelon Form)
    for (int i = 0; i < n; i++) {

        int maxRow = i;
        for (int k = i + 1; k < n; k++) {
            if (fabs(echelon[k][i]) > fabs(echelon[maxRow][i]))
                maxRow = k;
        }

        swap(echelon[i], echelon[maxRow]);

        if (fabs(echelon[i][i]) < 1e-9)
            continue;

        for (int k = i + 1; k < n; k++) {
            float factor = echelon[k][i] / echelon[i][i];
            for (int j = i; j <= n; j++) {
                echelon[k][j] -= factor * echelon[i][j];
            }
        }
    }

    // Printing Echelon Form
    out << "Echelon Form (Upper Triangular):\n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= n; j++)
            out << fixed << setprecision(4) << echelon[i][j] << " ";
        out << "\n";
    }
    out << "\n";

   
    int rankA = 0, rankAug = 0;
    const float EPS = 1e-9;

    for (int i = 0; i < n; i++) {
        bool nonZeroCoeff = false;
        bool nonZeroAug = false;

        for (int j = 0; j < n; j++) {
            if (fabs(echelon[i][j]) > EPS)
                nonZeroCoeff = true;
        }

        if (fabs(echelon[i][n]) > EPS)
            nonZeroAug = true;

        if (nonZeroCoeff)
            rankA++;

        if (nonZeroCoeff || nonZeroAug)
            rankAug++;
    }

    out << "System Classification:\n";

    if (rankA < rankAug) {
        out << "→ No Solution (Inconsistent System)\n";
        continue;
    }
    else if (rankA < n) {
        out << "→ Infinite Solutions (Dependent System)\n";
        continue;
    }
    else {
        out << "→ Unique Solution Exists\n\n";
    }

    
    for (int i = 0; i < n; i++) {

        int maxRow = i;
        for (int k = i + 1; k < n; k++) {
            if (fabs(a[k][i]) > fabs(a[maxRow][i]))
                maxRow = k;
        }

        swap(a[i], a[maxRow]);

        float pivot = a[i][i];
        if (fabs(pivot) < EPS) {
            out << "Numerical instability detected.\n";
            return 1;
        }

        // Normalize pivot row
        for (int j = 0; j <= n; j++)
            a[i][j] /= pivot;

        // Eliminate other rows
        for (int k = 0; k < n; k++) {
            if (k != i) {
                float factor = a[k][i];
                for (int j = 0; j <= n; j++)
                    a[k][j] -= factor * a[i][j];
            }
        }
    }

     // Printing Row Reduced Echelon Form
    out << "The Row Reduced Echelon Form:\n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= n; j++)
            out << fixed << setprecision(4) << a[i][j] << " ";
        out << "\n";
    }
    out << "\n";


    //Output Solution
    out << "Solution:\n";
    for (int i = 0; i < n; i++) {
        out << "x" << i + 1 << " = " << a[i][n] << "\n";
    }
    }

    return 0;
}

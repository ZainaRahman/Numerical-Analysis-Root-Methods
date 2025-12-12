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

    int k = 1;  // B is always n x 1

    // Augmented matrix (n x (n+1))
    vector<vector<double>> aug(n+1, vector<double>(n+2, 0));

    // A matrix (n x n)
    vector<vector<double>> a(n+1, vector<double>(n+1, 0));

    // B matrix (n x 1)
    vector<vector<double>> B(n+1, vector<double>(k+1, 0));

    vector<vector<double>> C(n+1, vector<double>(n+1, 0));
    vector<vector<double>> C1(n+1, vector<double>(n+1, 0));
    vector<vector<double>> res(n+1, vector<double>(k+1, 0));


    // ----- Read Augmented Matrix -----
    for (int i = 1; i <= n; i++)
        for (int j = 1; j <= n + 1; j++)
            in >> aug[i][j];

    // Split into A and B
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++)
            a[i][j] = aug[i][j];

        B[i][1] = aug[i][n+1];
    }

    // ------------------ Determinant ------------------
    int det1, ans = 0;
    int i = 1;

    for (int j = 1; j <= n; j++) {
        vector<int> R1, Cc1;

        for (int r = 1; r <= n; r++)
            if (r != i) R1.push_back(r);

        for (int c = 1; c <= n; c++)
            if (c != j) Cc1.push_back(c);

        det1 = a[R1[0]][Cc1[0]] * a[R1[1]][Cc1[1]]
             - a[R1[0]][Cc1[1]] * a[R1[1]][Cc1[0]];

        if ((i + j) % 2 != 0)
            det1 = -det1;

        det1 *= a[i][j];
        ans += det1;
    }

    // ------------------ DET = 0 → check solution type ------------------
    if (ans == 0) {

        bool noSolution = false, infinite = true;

        for (int r = 1; r <= n; r++) {
            bool allZeroA = true;

            for (int c = 1; c <= n; c++) {
                if (aug[r][c] != 0)
                    allZeroA = false;
            }

            if (allZeroA && aug[r][n+1] != 0) {
                noSolution = true;
                infinite = false;
                break;
            }
            if (!allZeroA)
                infinite = false;
        }

        if (noSolution) {
            out << "Determinant = 0 → No Solution (Inconsistent System)\n";
        }
        else {
            out << "Determinant = 0 → Infinite Solutions (Dependent System)\n";
        }

        return 0;
    }

    out << "Inverse exists. Determinant = " << ans << "\n\n";

    // ------------------ Cofactor Matrix ------------------
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            vector<int> R, Cc;

            for (int r = 1; r <= n; r++)
                if (r != i) R.push_back(r);

            for (int c = 1; c <= n; c++)
                if (c != j) Cc.push_back(c);

            int det = a[R[0]][Cc[0]] * a[R[1]][Cc[1]]
                    - a[R[0]][Cc[1]] * a[R[1]][Cc[0]];

            if ((i + j) % 2 != 0)
                det = -det;

            C[i][j] = det;
        }
    }

    out << "Cofactor Matrix:\n";
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++)
            out << C[i][j] << " ";
        out << "\n";
    }

    out << "\nAdjoint / Inverse Matrix:\n";
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            C1[i][j] = C[j][i] / ans;
            out << C1[i][j] << " ";
        }
        out << "\n";
    }

    // ------------------ Multiply Inverse * B ------------------
    for (int i = 1; i <= n; i++) {
        double sum = 0;

        for (int t = 1; t <= n; t++)
            sum += C1[i][t] * B[t][1];

        res[i][1] = sum;
    }

    out << "\nSolution Vector X:\n";
    for (int i = 1; i <= n; i++)
        out << "x" << i << " = " << res[i][1] << "\n";

    return 0;
}

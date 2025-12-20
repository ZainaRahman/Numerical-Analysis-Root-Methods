#include <bits/stdc++.h>
using namespace std;

// Cofactor
void getCofactor(const vector<vector<double>>& A,
                 vector<vector<double>>& temp,
                 int p, int q, int n)
{
    int i = 1, j = 1;
    for (int row = 1; row <= n; row++) {
        for (int col = 1; col <= n; col++) {
            if (row != p && col != q) {
                temp[i][j++] = A[row][col];
                if (j == n) {
                    j = 1;
                    i++;
                }
            }
        }
    }
}

//Recursive Determinant
double determinant(const vector<vector<double>>& A, int n)
{
    if (n == 1)
        return A[1][1];

    double det = 0;
    int sign = 1;
    vector<vector<double>> temp(n, vector<double>(n, 0));

    for (int i = 1; i <= n; i++) {
        getCofactor(A, temp, 1, i, n);
        det += sign * A[1][i] * determinant(temp, n - 1);
        sign = -sign;
    }
    return det;
}

int main()
{
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    if (!fin) {
        cout << "Error opening input file.\n";
        return 1;
    }

    int casing;
    fin>> casing;
    int r=1;
    while(casing--){
    fout<<endl;
    fout<<"Case "<<r++<<":\n";
    int n;
    fin >> n;

    vector<vector<double>> aug(n+1, vector<double>(n+2, 0));
    vector<vector<double>> a(n+1, vector<double>(n+1, 0));
    vector<vector<double>> B(n+1, vector<double>(2, 0));
    vector<vector<double>> C(n+1, vector<double>(n+1, 0));
    vector<vector<double>> C1(n+1, vector<double>(n+1, 0));
    vector<vector<double>> res(n+1, vector<double>(2, 0));

    for (int i = 1; i <= n; i++)
        for (int j = 1; j <= n + 1; j++)
            fin >> aug[i][j];

    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++)
            a[i][j] = aug[i][j];
        B[i][1] = aug[i][n+1];
    }

    double detA = determinant(a, n);

    if (fabs(detA) < 1e-9) {
        // Row reduce to check rank by gauss jordan elimination
        vector<vector<double>> tempAug = aug;
        const double EPS = 1e-9;
        
        for (int i = 1; i <= n; i++) {
            // Find pivot
            int maxRow = i;
            for (int k = i + 1; k <= n; k++) {
                if (fabs(tempAug[k][i]) > fabs(tempAug[maxRow][i]))
                    maxRow = k;
            }
            swap(tempAug[i], tempAug[maxRow]);
            
            if (fabs(tempAug[i][i]) < EPS) continue;
            
            // Eliminate below
            for (int k = i + 1; k <= n; k++) {
                double factor = tempAug[k][i] / tempAug[i][i];
                for (int j = i; j <= n + 1; j++) {
                    tempAug[k][j] -= factor * tempAug[i][j];
                }
            }
        }
        
        // Check for inconsistency:  row with all zeros in A but non-zero in b
        bool noSol = false;
        for (int i = 1; i <= n; i++) {
            bool allZero = true;
            for (int j = 1; j <= n; j++) {
                if (fabs(tempAug[i][j]) > EPS) {
                    allZero = false;
                    break;
                }
            }
            if (allZero && fabs(tempAug[i][n+1]) > EPS) {
                noSol = true;
                break;
            }
        }

        if (noSol)
            fout << "Determinant = 0 → No Solution (Inconsistent System)\n";
        else
            fout << "Determinant = 0 → Infinite Solutions (Dependent System)\n";

        continue;
    }

    fout << "Determinant = " << detA << "\n\n";

    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            vector<vector<double>> temp(n, vector<double>(n, 0));
            getCofactor(a, temp, i, j, n);
            C[i][j] = pow(-1, i + j) * determinant(temp, n - 1);
        }
    }

    fout << "Inverse Matrix:\n";
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            C1[i][j] = C[j][i] / detA;
            fout << C1[i][j] << " ";
        }
        fout << "\n";
    }

    for (int i = 1; i <= n; i++) {
        for (int t = 1; t <= n; t++)
            res[i][1] += C1[i][t] * B[t][1];
    }

    fout << "\nSolution Vector:\n";
    for (int i = 1; i <= n; i++)
        fout << "x" << i << " = " << res[i][1] << "\n";
    fout << "\n";
}

    fin.close();
    fout.close();


    return 0;
}
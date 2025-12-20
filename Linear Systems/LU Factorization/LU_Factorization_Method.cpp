#include <bits/stdc++.h>
using namespace std;

int main()
{

    // ---- File Handling ----
    ifstream in("input.txt");
    ofstream out("output.txt");

    if (!in)
    {
        cerr << "Error: input.txt not found\n";
        return 1;
    }

    int casing;
    in >> casing;
    int r = 1;

    while (casing--)
    {
        cout << endl;
        out << "----- Case " << r++ << " -----\n\n";
        int n;
        in >> n;

        vector<vector<double>> a(n + 1, vector<double>(n + 2));

        // Read augmented matrix A|b
        for (int i = 1; i <= n; i++)
        {
            for (int j = 1; j <= n + 1; j++)
            {
                in >> a[i][j];
            }
        }

        vector<vector<double>> u(n + 1, vector<double>(n + 1, 0));
        vector<vector<double>> l(n + 1, vector<double>(n + 1, 0));

        for (int i = 1; i <= n; i++)
        {
            l[i][i] = 1;
        }

        // ---- LU Decomposition ----
        for (int i = 1; i <= n; i++)
        {
            for (int j = 1; j <= n; j++)
            {

                if (i <= j)
                {
                    u[i][j] = a[i][j];
                    for (int k = 1; k < i; k++)
                        u[i][j] -= l[i][k] * u[k][j];
                }
                else
                {
                    l[i][j] = a[i][j];
                    for (int k = 1; k < j; k++)
                        l[i][j] -= l[i][k] * u[k][j];

                    if (u[j][j] == 0)
                    {
                        out << "Matrix is singular. Cannot compute LU decomposition.\n";
                        return 0;
                    }

                    l[i][j] /= u[j][j];
                }
            }
        }

        // ---- Print U Matrix ----
        out << "U Matrix:\n";
        for (int i = 1; i <= n; i++)
        {
            for (int j = 1; j <= n; j++)
            {
                out << u[i][j] << " ";
            }
            out << "\n";
        }

        // ---- Print L Matrix ----
        out << "\nL Matrix:\n";
        for (int i = 1; i <= n; i++)
        {
            for (int j = 1; j <= n; j++)
            {
                out << l[i][j] << " ";
            }
            out << "\n";
        }

        // ---- Forward Substitution: Ly = b ----
        vector<double> y(n + 1, 0);

        for (int i = 1; i <= n; i++)
        {
            y[i] = a[i][n + 1];
            for (int k = 1; k < i; k++)
            {
                y[i] -= l[i][k] * y[k];
            }
        }

        // ---- Check Solution Type ----
        bool noSolution = false;
        bool infiniteSolution = false;

        for (int i = 1; i <= n; i++)
        {
            bool allZero = true;

            for (int j = 1; j <= n; j++)
            {
                if (fabs(u[i][j]) > 1e-9)
                {
                    allZero = false;
                    break;
                }
            }

            if (allZero)
            {
                if (fabs(y[i]) > 1e-9)
                {
                    noSolution = true;
                }
                else
                {
                    infiniteSolution = true;
                }
            }
        }

        if (noSolution)
        {
            out << "\nThe system has NO SOLUTION (Inconsistent equations).\n";
            continue;
        }

        if (infiniteSolution)
        {
            out << "\nThe system has INFINITE SOLUTIONS (Dependent equations).\n";
            continue;
        }

        out << "\nThe system has a UNIQUE SOLUTION.\n";

        // ---- Backward Substitution: Ux = y ----
        vector<double> ans(n + 1, 0);

        for (int i = n; i >= 1; i--)
        {
            ans[i] = y[i];
            for (int k = i + 1; k <= n; k++)
            {
                ans[i] -= u[i][k] * ans[k];
            }
            ans[i] /= u[i][i];
        }

        // ---- Print Solution ----
        out << "\nFinal Solution (x values):\n";
        for (int i = 1; i <= n; i++)
        {
            out << "x" << i << " = " << ans[i] << "\n";
        }
    }

    return 0;
}

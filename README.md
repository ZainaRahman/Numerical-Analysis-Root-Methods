# Numerical-Analysis-Root-Methods
# Table of Contents

- [Solution of Linear Equations](#solution-of-linear-equations)
  - [Gauss Elimination Method](#gauss-elimination-method)
    - [Theory](#gauss-elimination-theory)
    - [Code](#gauss-elimination-code)
    - [Input](#gauss-elimination-input)
    - [Output](#gauss-elimination-output)
  - [Gauss Jordan Elimination Method](#gauss-jordan-elimination-method)
    - [Theory](#gauss-jordan-theory)
    - [Code](#gauss-jordan-code)
    - [Input](#gauss-jordan-input)
    - [Output](#gauss-jordan-output)
  - [LU Decomposition Method](#lu-decomposition-method)
    - [Theory](#lu-decomposition-theory)
    - [Code](#lu-decomposition-code)
    - [Input](#lu-decomposition-input)
    - [Output](#lu-decomposition-output)
  - [Matrix Inversion](#matrix-inversion)
    - [Theory](#matrix-inversion-theory)
    - [Code](#matrix-inversion-code)
    - [Input](#matrix-inversion-input)
    - [Output](#matrix-inversion-output)

- [Solution of Non-Linear Equations](#solution-of-non-linear-equations)
  - [Bisection Method](#bisection-method)
    - [Theory](#bisection-theory)
    - [Code](#bisection-code)
    - [Input](#bisection-input)
    - [Output](#bisection-output)
  - [False Position Method](#false-position-method)
    - [Theory](#false-position-theory)
    - [Code](#false-position-code)
    - [Input](#false-position-input)
    - [Output](#false-position-output)
  - [Secant Method](#secant-method)
    - [Theory](#secant-theory)
    - [Code](#secant-code)
    - [Input](#secant-input)
    - [Output](#secant-output)
  - [Newton Raphson Method](#newton-raphson-method)
    - [Theory](#newton-raphson-theory)
    - [Code](#newton-raphson-code)
    - [Input](#newton-raphson-input)
    - [Output](#newton-raphson-output)

- [Solution of Interpolation](#solution-of-interpolation)
  - [Newton's Forward Interpolation Method](#newtons-forward-interpolation-method)
    - [Theory](#newtons-forward-interpolation-theory)
    - [Code](#newtons-forward-interpolation-code)
    - [Input](#newtons-forward-interpolation-input)
    - [Output](#newtons-forward-interpolation-output)
  - [Newton's Backward Interpolation Method](#newtons-backward-interpolation-method)
    - [Theory](#newtons-backward-interpolation-theory)
    - [Code](#newtons-backward-interpolation-code)
    - [Input](#newtons-backward-interpolation-input)
    - [Output](#newtons-backward-interpolation-output)
  - [Divided Difference Method](#divided-difference-method)
    - [Theory](#divided-difference-theory)
    - [Code](#divided-difference-code)
    - [Input](#divided-difference-input)
    - [Output](#divided-difference-output)

- [Solution of Curve Fitting Model](#solution-of-curve-fitting-model)
  - [Least Square Regression Method For Linear Equations](#least-square-regression-method-for-linear-equations-method)
    - [Theory](#least-square-regression-method-for-linear-equations-theory)
    - [Code](#least-square-regression-method-for-linear-equations-code)
    - [Input](#least-square-regression-method-for-linear-equations-input)
    - [Output](#least-square-regression-method-for-linear-equations-output)
  - [Least Square Regression Method For Transcendental Equations](#least-square-regression-method-for-transcendental-equations)
    - [Theory](#least-square-regression-method-for-transcendental-equations-theory)
    - [Code](#least-square-regression-method-for-transcendental-equations-code)
    - [Input](#least-square-regression-method-for-transcendental-equations-input)
    - [Output](#least-square-regression-method-for-transcendental-equations-output)
  - [Least Square Regression Method For Polynomial Equations](#least-square-regression-method-for-polynomial-equations)
    - [Theory](#least-square-regression-method-for-polynomial-equations-theory)
    - [Code](#least-square-regression-method-for-polynomial-equations-code)
    - [Input](#least-square-regression-method-for-polynomial-equations-input)
    - [Output](#least-square-regression-method-for-polynomial-equations-output)

- [Solution of Differential Equations](#solution-of-differential-equations)
  - [Equal Interval Interpolation Method](#equal-interval-interpolation-method)
    - [Theory](#equal-interval-interpolation-theory)
    - [Code](#equal-interval-interpolation-code)
    - [Input](#equal-interval-interpolation-input)
    - [Output](#equal-interval-interpolation-output)
  - [Second Order Derivative Method](#second-order-derivative-method)
    - [Theory](#second-order-derivative-theory)
    - [Code](#second-order-derivative-code)
    - [Input](#second-order-derivative-input)
    - [Output](#second-order-derivative-output)
  - [Runge Kutta Method](#runge-kutta-method)
    - [Theory](#runge-kutta-theory)
    - [Code](#runge-kutta-code)
    - [Input](#runge-kutta-input)
    - [Output](#runge-kutta-output)
  - [Numerical Differentiation Method](#numerical-differentiation-method)
    - [Theory](#numerical-differentiation-theory)
    - [Code](#numerical-differentiation-code)
    - [Input](#numerical-differentiation-input)
    - [Output](#numerical-differentiation-output)

- [Solution of Numerical Integrations](#solution-of-numerical-integrations)
  - [Simpson's One-Third Rule](#simpsons-one-third-rule)
    - [Theory](#simpsons-one-third-rule-theory)
    - [Code](#simpsons-one-third-rule-code)
    - [Input](#simpsons-one-third-rule-input)
    - [Output](#simpsons-one-third-rule-output)
  - [Simpson's Three-Eighths Rule](#simpsons-three-eighths-rule)
    - [Theory](#simpsons-three-eighths-rule-theory)
    - [Code](#simpsons-three-eighths-rule-code)
    - [Input](#simpsons-three-eighths-rule-input)
    - [Output](#simpsons-three-eighths-rule-output)
  


---

### Solution of Linear Equations

### Gauss Elimination Method

#### Gauss Elimination Theory
[Add your theory content here]

#### Gauss Elimination Code
```cpp
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

    int n;
    in >> n;

    vector<vector<float>> a(n, vector<float>(n + 1));

    // augmented matrix
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= n; j++) {
            in >> a[i][j];
        }
    }

    // Copying original matrix for echelon form
    vector<vector<float>> echelon = a;

    //Forward Elimination (Echelon Form)
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
        return 0;
    }
    else if (rankA < n) {
        out << "→ Infinite Solutions (Dependent System)\n";
        return 0;
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

        for (int j = 0; j <= n; j++)
            a[i][j] /= pivot;

        for (int k = 0; k < n; k++) {
            if (k != i) {
                float factor = a[k][i];
                for (int j = 0; j <= n; j++)
                    a[k][j] -= factor * a[i][j];
            }
        }
    }

    //  Output Solution 
    out << "Solution:\n";
    for (int i = 0; i < n; i++) {
        out << "x" << i + 1 << " = " << a[i][n] << "\n";
    }

    return 0;
}

```

#### Gauss Elimination Input
```
3
2 1 -1 8
-3 -1 2 -11
-2 1 2 -3
```

#### Gauss Elimination Output
```
Echelon Form (Upper Triangular):
-3.0000 -1.0000 2.0000 -11.0000 
0.0000 1.6667 0.6667 4.3333 
0.0000 0.0000 0.2000 -0.2000 

System Classification:
→ Unique Solution Exists

Solution:
x1 = 2.0000
x2 = 3.0000
x3 = -1.0000

```

---

### Gauss Jordan Elimination Method

#### Gauss Jordan Theory
[Add your theory content here]

#### Gauss Jordan Code
```cpp
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

    int n;
    in >> n;

    vector<vector<float>> a(n, vector<float>(n + 1));

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
        return 0;
    }
    else if (rankA < n) {
        out << "→ Infinite Solutions (Dependent System)\n";
        return 0;
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

    //Output Solution
    out << "Solution:\n";
    for (int i = 0; i < n; i++) {
        out << "x" << i + 1 << " = " << a[i][n] << "\n";
    }

    return 0;
}

```

#### Gauss Jordan Input
```
5
2 1 -1 3 2 9
1 3 2 -1 1 8
3 2 4 1 -2 20
2 1 3 2 1 17
1 -1 2 3 4 15
```

#### Gauss Jordan Output
```
Echelon Form (Upper Triangular):
3.0000 2.0000 4.0000 1.0000 -2.0000 20.0000 
0.0000 2.3333 0.6667 -1.3333 1.6667 1.3333 
0.0000 0.0000 -3.5714 2.1429 3.5714 -4.1429 
0.0000 0.0000 0.0000 2.4000 7.0000 7.9600 
0.0000 0.0000 0.0000 0.0000 -1.0833 -1.2833 

System Classification:
→ Unique Solution Exists

Solution:
x1 = 5.1538
x2 = -1.0000
x3 = 2.2615
x4 = -0.1385
x5 = 1.1846

```

---

### LU Decomposition Method

#### LU Decomposition Theory
[Add your theory content here]

#### LU Decomposition Code
```cpp
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

    int n;
    in >> n;

    vector<vector<double>> a(n + 1, vector<double>(n + 2));

    // Reading augmented matrix A|b
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n + 1; j++) {
            in >> a[i][j];
        }
    }

    vector<vector<double>> u(n + 1, vector<double>(n + 1, 0));
    vector<vector<double>> l(n + 1, vector<double>(n + 1, 0));

    for (int i = 1; i <= n; i++) {
        l[i][i] = 1;
    }

    //  LU Decomposition
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {

            if (i <= j) {
                u[i][j] = a[i][j];
                for (int k = 1; k < i; k++)
                    u[i][j] -= l[i][k] * u[k][j];
            }
            else {
                l[i][j] = a[i][j];
                for (int k = 1; k < j; k++)
                    l[i][j] -= l[i][k] * u[k][j];

                if (u[j][j] == 0) {
                    out << "Matrix is singular. Cannot compute LU decomposition.\n";
                    return 0;
                }

                l[i][j] /= u[j][j];
            }
        }
    }

    // Printing U Matrix
    out << "U Matrix:\n";
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            out << u[i][j] << " ";
        }
        out << "\n";
    }

    // Printing L Matrix
    out << "\nL Matrix:\n";
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            out << l[i][j] << " ";
        }
        out << "\n";
    }

    // Checking Solution Type
    bool noSolution = false;
    bool infiniteSolution = false;

    for (int i = 1; i <= n; i++) {
        bool allZero = true;

        for (int j = 1; j <= n; j++) {
            if (fabs(u[i][j]) > 1e-9) {
                allZero = false;
                break;
            }
        }

        if (allZero) {
            double rhs = a[i][n + 1];

            if (fabs(rhs) > 1e-9) {
                noSolution = true;
            }
            else {
                infiniteSolution = true;
            }
        }
    }

    if (noSolution) {
        out << "\nThe system has NO SOLUTION (Inconsistent equations).\n";
        return 0;
    }

    if (infiniteSolution) {
        out << "\nThe system has INFINITE SOLUTIONS (Dependent equations).\n";
        return 0;
    }

    out << "\nThe system has a UNIQUE SOLUTION.\n";

    //  Forward Substitution: Ly = b
    vector<double> y(n + 1, 0);

    for (int i = 1; i <= n; i++) {
        y[i] = a[i][n + 1];
        for (int k = 1; k < i; k++) {
            y[i] -= l[i][k] * y[k];
        }
    }

    //  Backward Substitution: Ux = y 
    vector<double> ans(n + 1, 0);

    for (int i = n; i >= 1; i--) {
        ans[i] = y[i];
        for (int k = i + 1; k <= n; k++) {
            ans[i] -= u[i][k] * ans[k];
        }
        ans[i] /= u[i][i];
    }

    // Printing Solution
    out << "\nFinal Solution (x values):\n";
    for (int i = 1; i <= n; i++) {
        out << "x" << i << " = " << ans[i] << "\n";
    }

    return 0;
}

```

#### LU Decomposition Input
```
5
2 1 -1 3 2 9
1 3 2 -1 1 8
3 2 4 1 -2 20
2 1 3 2 1 17
1 -1 2 3 4 15
```

#### LU Decomposition Output
```
U Matrix:
2 1 -1 3 2 
0 2.5 2.5 -2.5 0 
0 0 5 -3 -5 
0 0 0 1.4 3 
0 0 0 0 1.85714 

L Matrix:
1 0 0 0 0 
0.5 1 0 0 0 
1.5 0.2 1 0 0 
1 0 0.8 1 0 
0.5 -0.6 0.8 1.71429 1 

The system has a UNIQUE SOLUTION.

Final Solution (x values):
x1 = 5.15385
x2 = -1
x3 = 2.26154
x4 = -0.138462
x5 = 1.18462

```

---

### Matrix Inversion

#### Matrix Inversion Theory
[Add your theory content here]

#### Matrix Inversion Code
```cpp
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

    if (detA == 0) {
        bool noSol = false, infiniteSol = true;

        for (int r = 1; r <= n; r++) {
            bool allZero = true;
            for (int c = 1; c <= n; c++)
                if (aug[r][c] != 0)
                    allZero = false;

            if (allZero && aug[r][n+1] != 0) {
                noSol = true;
                infiniteSol = false;
                break;
            }
            if (!allZero)
                infiniteSol = false;
        }

        if (noSol)
            fout << "Determinant = 0 → No Solution (Inconsistent System)\n";
        else
            fout << "Determinant = 0 → Infinite Solutions (Dependent System)\n";

        return 0;
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

    fin.close();
    fout.close();

    return 0;
}

```

#### Matrix Inversion Input
```
5
2 1 -1 3 2 9
1 3 2 -1 1 8
3 2 4 1 -2 20
2 1 3 2 1 17
1 -1 2 3 4 15
```

#### Matrix Inversion Output
```
Determinant = 65

Inverse Matrix:
0.384615 0.384615 1.92308 -3.76923 1.61538 
-0 0 -1 2 -1 
-0.246154 -0.0461538 -0.230769 0.692308 -0.153846 
-0.0461538 -0.446154 -1.23077 2.69231 -1.15385 
0.0615385 0.261538 0.307692 -0.923077 0.538462 

Solution Vector:
x1 = 5.15385
x2 = -1
x3 = 2.26154
x4 = -0.138462
x5 = 1.18462

```

---

### Solution of Non-Linear Equations

### Bisection Method

#### Bisection Theory
[Add your theory content here]

#### Bisection Code
```python
# Add your code here
```

#### Bisection Input
```
[Add your input format here]
```

#### Bisection Output
```
[Add your output format here]
```

---

### False Position Method

#### False Position Theory
[Add your theory content here]

#### False Position Code
```python
# Add your code here
```

#### False Position Input
```
[Add your input format here]
```

#### False Position Output
```
[Add your output format here]
```

---

### Secant Method

#### Secant Theory
[Add your theory content here]

#### Secant Code
```python
# Add your code here
```

#### Secant Input
```
[Add your input format here]
```

#### Secant Output
```
[Add your output format here]
```

---

### Newton Raphson Method

#### Newton Raphson Theory
[Add your theory content here]

#### Newton Raphson Code
```python
# Add your code here
```

#### Newton Raphson Input
```
[Add your input format here]
```

#### Newton Raphson Output
```
[Add your output format here]
```

---

### Solution of Interpolation

### Newton's Forward Interpolation Method

#### Newton's Forward Interpolation Theory
[Add your theory content here]

#### Newton's Forward Interpolation Code
```python
# Add your code here
```

#### Newton's Forward Interpolation Input
```
[Add your input format here]
```

#### Newton's Forward Interpolation Output
```
[Add your output format here]
```

---

### Newton's Backward Interpolation Method

#### Newton's Backward Interpolation Theory
[Add your theory content here]

#### Newton's Backward Interpolation Code
```python
# Add your code here
```

#### Newton's Backward Interpolation Input
```
[Add your input format here]
```

#### Newton's Backward Interpolation Output
```
[Add your output format here]
```

---

### Divided Difference Method

#### Divided Difference Theory
[Add your theory content here]

#### Divided Difference Code
```python
# Add your code here
```

#### Divided Difference Input
```
[Add your input format here]
```

#### Divided Difference Output
```
[Add your output format here]
```

---

### Solution of Curve Fitting Model

### Least Square Regression Method For Linear Equations Method

#### Least Square Regression Method For Linear Equations Theory
[Add your theory content here]

#### Least Square Regression Method For Linear Equations Code
```python
# Add your code here
```

#### Least Square Regression Method For Linear Equations Input
```
[Add your input format here]
```

#### Least Square Regression Method For Linear Equations Output
```
[Add your output format here]
```

---

### Least Square Regression Method For Transcendental Equations 

#### Least Square Regression Method For Transcendental Equations Theory
[Add your theory content here]

#### Least Square Regression Method For Transcendental Equations Code
```python
# Add your code here
```

#### Least Square Regression Method For Transcendental Equations Input
```
[Add your input format here]
```

#### Least Square Regression Method For Transcendental Equations Output
```
[Add your output format here]
```

---

### Least Square Regression Method For Polynomial Equations 

#### Least Square Regression Method For Polynomial Equations Theory
[Add your theory content here]

#### Least Square Regression Method For Polynomial Equations Code
```python
# Add your code here
```

#### Least Square Regression Method For Polynomial Equations Input
```
[Add your input format here]
```

#### Least Square Regression Method For Polynomial Equations Output
```
[Add your output format here]
```

---

### Solution of Differential Equations

### Equal Interval Interpolation Method

#### Equal Interval Interpolation Theory
[Add your theory content here]

#### Equal Interval Interpolation Code
```python
# Add your code here
```

#### Equal Interval Interpolation Input
```
[Add your input format here]
```

#### Equal Interval Interpolation Output
```
[Add your output format here]
```

---

### Second Order Derivative Method 

#### Second Order Derivative Theory
[Add your theory content here]

#### Second Order Derivative Code
```python
# Add your code here
```

#### Second Order Derivative Input
```
[Add your input format here]
```

#### Second Order Derivative Output
```
[Add your output format here]
```

---

### Runge Kutta Method 

#### Runge Kutta Theory
[Add your theory content here]

#### Runge Kutta Code
```python
# Add your code here
```

#### Runge Kutta Input
```
[Add your input format here]
```

#### Runge Kutta Output
```
[Add your output format here]
```

---

### Numerical Differentiation Method

#### Numerical Differentiation Theory
[Add your theory content here]

#### Numerical Differentiation Code
```python
# Add your code here
```

#### Numerical Differentiation Input
```
[Add your input format here]
```

#### Numerical Differentiation Output
```
[Add your output format here]
```

---

### Solution of Numerical Integrations

### Simpson's One-Third Rule

#### Simpson's One-Third Rule Theory
[Add your theory content here]

#### Simpson's One-Third Rule Code
```python
# Add your code here
```

#### Simpson's One-Third Rule Input
```
[Add your input format here]
```

#### Simpson's One-Third Rule Output
```
[Add your output format here]
```

---

### Simpson's Three-Eighths Rule 

#### Simpson's Three-Eighths Rule Theory
[Add your theory content here]

#### Simpson's Three-Eighths Rule Code
```python
# Add your code here
```

#### Simpson's Three-Eighths Rule Input
```
[Add your input format here]
```

#### Simpson's Three-Eighths Rule Output
```
[Add your output format here]
```

---


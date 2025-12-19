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

## Solution of Linear Equations

### Gauss Elimination Method

### Gauss Elimination Theory

#### Method used
**Gauss Elimination Method**

#### Objective
To solve a system of linear algebraic equations by transforming it into an **upper triangular system**, followed by **back substitution**.

#### Data Requirement
A system of `n` linear equations:
```
a11x1 + a12x2 + ... + a1nxn = b1
a21x1 + a22x2 + ... + a2nxn = b2
...
an1x1 + an2x2 + ... + annxn = bn
```

Matrix form:
```
AX = B
```

#### Notation
- `A = [aij]` : coefficient matrix of order `n × n`
- `X = [x1, x2, ..., xn]ᵀ` : vector of unknowns
- `B = [b1, b2, ..., bn]ᵀ` : constant vector

#### Core Idea
The system is simplified by eliminating variables using **elementary row operations** to obtain an upper triangular matrix.

#### Elimination Approach (Formula)
To eliminate element `aij` (where `j < i`):
```
Ri ← Ri − (aij / ajj) Rj
```


#### Phases Involved

##### Forward Elimination
Transforms the augmented matrix `[A | B]` into an **upper triangular form**.

##### Back Substitution
Solutions are obtained using:
```
xn = bn / ann
xi = (1 / aii) [ bi − Σ (aij xj) ], j = i+1 to n
```

#### Accuracy Considerations
- Exact in theory
- Rounding errors may occur in floating-point arithmetic
- Pivoting improves numerical stability

#### Applicability
- Suitable for small to medium-sized systems
- Widely used due to conceptual simplicity



### Gauss Elimination Code
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

### Gauss Elimination Input
```
3
2 1 -1 8
-3 -1 2 -11
-2 1 2 -3
```

### Gauss Elimination Output
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

### Gauss Jordan Theory
#### Method used
**Gauss–Jordan Elimination Method**

#### Objective
To solve a system of linear equations by reducing the augmented matrix directly to **Reduced Row Echelon Form (RREF)**.

#### Data Requirement
Augmented matrix form:
```
[A | B]
```

#### Notation
- `aij` : element of coefficient matrix
- `bi`  : element of constant vector

#### Core Idea
Each pivot element is made **unity**, and all other elements in its column are eliminated, producing an identity matrix.

#### Elimination Approach (Formulae)

**Normalization of pivot row:**
```
Ri ← Ri / aii
```

**Elimination of other rows:**
```
Rj ← Rj − aji Ri (j ≠ i)
```
#### Matrix Form Obtained
```
[I | X]
```

where `I` is the identity matrix and `X` contains the solution.

#### Evaluation Process
The solution is obtained directly as:
```
xi = bi
```
#### Accuracy Considerations
- More computationally expensive than Gauss Elimination
- Sensitive to rounding errors for large systems

#### Applicability
- Used when a direct solution or matrix inverse is required


### Gauss Jordan Code
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

### Gauss Jordan Input
```
5
2 1 -1 3 2 9
1 3 2 -1 1 8
3 2 4 1 -2 20
2 1 3 2 1 17
1 -1 2 3 4 15
```

### Gauss Jordan Output
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

### LU Decomposition Theory
#### Method used
**LU Decomposition Method**

#### Objective
To solve a system of linear equations by factorizing the coefficient matrix into lower and upper triangular matrices.

#### Data Requirement
A square matrix with non-zero pivots.

#### Core Idea (Formula)
```
A = LU
```

#### Notation
- `L = [lij]` : lower triangular matrix
- `U = [uij]` : upper triangular matrix
- `A` : coefficient matrix
- `X` : solution vector
- `B` : constant vector

#### Solution Process

**Step 1:** Solve
```
LY = B
```
using forward substitution:
```
yi = bi − Σ (lij yj), j = 1 to i−1
```

**Step 2:** Solve
```
UX = Y
```
using back substitution:
```
xi = (1 / uii) [ yi − Σ (uij xj) ], j = i+1 to n
```

#### Accuracy Considerations
- More efficient than repeated Gauss Elimination
- Numerical stability improves with pivoting

#### Applicability
- Ideal for solving multiple systems with the same coefficient matrix


### LU Decomposition Code
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

### LU Decomposition Input
```
5
2 1 -1 3 2 9
1 3 2 -1 1 8
3 2 4 1 -2 20
2 1 3 2 1 17
1 -1 2 3 4 15
```

### LU Decomposition Output
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

### Matrix Inversion Theory

#### Method used
**Matrix Inversion Method**

#### Objective
To solve a system of linear equations using the inverse of the coefficient matrix.

#### Data Requirement
A square, non-singular matrix (`det(A) ≠ 0`).

#### Core Idea (Formula)
Given:
```
AX = B
```
The solution is:
```
X = A⁻¹ B
```
#### Notation
- `A⁻¹` : inverse of matrix `A`
- `I` : identity matrix

#### Inversion Approach
The inverse is computed using Gauss–Jordan elimination:
```
[A | I] → [I | A⁻¹]
```

#### Evaluation Process
Once `A⁻¹` is obtained, the solution vector is computed using matrix multiplication.

#### Accuracy Considerations
- Computationally expensive
- Sensitive to rounding errors
- Not recommended for large systems

#### Applicability
- Useful for theoretical analysis
- Suitable for small systems only


### Matrix Inversion Code
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

### Matrix Inversion Input
```
5
2 1 -1 3 2 9
1 3 2 -1 1 8
3 2 4 1 -2 20
2 1 3 2 1 17
1 -1 2 3 4 15
```

### Matrix Inversion Output
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

## Solution of Non-Linear Equations

### Bisection Method

### Bisection Theory
[Add your theory content here]

### Bisection Code
```cpp
# Add your code here
```

### Bisection Input
```
[Add your input format here]
```

### Bisection Output
```
[Add your output format here]
```

---

### False Position Method

### False Position Theory
[Add your theory content here]

### False Position Code
```cpp
# Add your code here
```

### False Position Input
```
[Add your input format here]
```

### False Position Output
```
[Add your output format here]
```

---

### Secant Method

#### Secant Theory
[Add your theory content here]

#### Secant Code
```cpp
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
```cpp
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

### Newton's Forward Interpolation Theory
[Add your theory content here]

### Newton's Forward Interpolation Code
```cpp
# Add your code here
```

### Newton's Forward Interpolation Input
```
[Add your input format here]
```

### Newton's Forward Interpolation Output
```
[Add your output format here]
```

---

### Newton's Backward Interpolation Method

### Newton's Backward Interpolation Theory
[Add your theory content here]

### Newton's Backward Interpolation Code
```cpp
#include <bits/stdc++.h>
using namespace std;

int fact(int n)
{
    int f = 1;
    for(int i = 2; i <= n; i++)
        f *= i;
    return f;
}

int detectOrderBackward(const vector<vector<double>>& diff, double eps = 1e-9)
{
    int n = diff.size();

    for(int j = 1; j < n; j++)
    {
        bool allZero = true;
        for(int i = j; i < n; i++)
        {
            if (fabs(diff[i][j]) > eps)
            {
                allZero = false;
                break;
            }
        }
        if (allZero)
            return j - 1;
    }
    return n - 1;
}

int main()
{
    ifstream inputFile("Input.txt");
    ofstream outputFile("Output.txt");

    if (!inputFile.is_open())
        return 1;

    int testCases;
    inputFile >> testCases;

    for (int tc = 1; tc <= testCases; tc++)
    {
        int n;
        inputFile >> n;

        vector<double> x(n), y(n);
        for (int i = 0; i < n; i++) inputFile >> x[i];
        for (int i = 0; i < n; i++) inputFile >> y[i];

        double X;
        inputFile >> X;

        double h = x[1] - x[0];
        double u = (X - x[n-1]) / h;

        vector<vector<double>> diff(n, vector<double>(n, 0.0));
        for (int i = 0; i < n; i++)
            diff[i][0] = y[i];

        for (int j = 1; j < n; j++)
            for (int i = n-1; i >= j; i--)
                diff[i][j] = diff[i][j-1] - diff[i-1][j-1];

        cout << fixed << setprecision(2);
        outputFile << fixed << setprecision(2);

        cout << "Test Case " << tc << "\n";
        outputFile << "Test Case " << tc << "\n";

        cout << "Backward Difference Table:\n";
        outputFile << "Backward Difference Table:\n";

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j <= i; j++)
            {
                cout << setw(10) << diff[i][j];
                outputFile << setw(10) << diff[i][j];
            }
            cout << "\n";
            outputFile << "\n";
        }

        double fx = diff[n-1][0];
        double u_term = 1.0;

        for (int k = 1; k < n; k++)
        {
            u_term *= (u + (k - 1));
            fx += (u_term * diff[n-1][k]) / fact(k);
        }

        int order = detectOrderBackward(diff);

        double error = 0.0;
        if (order + 1 < n)
        {
            double u_err = 1.0;
            for (int i = 0; i <= order; i++)
                u_err *= (u + i);

            error = fabs((u_err * diff[n-1][order + 1]) / fact(order + 1));
        }

        cout << "h = " << h << "\n";
        cout << "u = " << u << "\n";
        cout << "Detected Polynomial Order = " << order << "\n";
        cout << "Interpolated value f(" << X << ") = " << fx << "\n";
        cout << "Estimated Backward Interpolation Error = " << error << "\n\n";

        outputFile << "h = " << h << "\n";
        outputFile << "u = " << u << "\n";
        outputFile << "Detected Polynomial Order = " << order << "\n";
        outputFile << "Interpolated value f(" << X << ") = " << fx << "\n";
        outputFile << "Estimated Backward Interpolation Error = " << error << "\n\n";
    }

    inputFile.close();
    outputFile.close();
    return 0;
}

```

### Newton's Backward Interpolation Input
```
[Add your input format here]
```

### Newton's Backward Interpolation Output
```
[Add your output format here]
```

---

### Divided Difference Method

### Divided Difference Theory
[Add your theory content here]

### Divided Difference Code
```cpp
# Add your code here
```

### Divided Difference Input
```
[Add your input format here]
```

### Divided Difference Output
```
[Add your output format here]
```

---

## Solution of Curve Fitting Model

### Least Square Regression Method For Linear Equations Method

### Least Square Regression Method For Linear Equations Theory
#### Method used
**Least Squares Regression – Linear Equation**

#### Objective
To fit a straight line of the form 
```
y = a + bx
```
that best represents the given experimental data.

#### Data Requirement
A set of `n` observed data points:
```
(x1, y1), (x2, y2), ..., (xn, yn)
```

#### Core Idea
The best-fitting straight line is obtained by minimizing the **sum of squares of vertical deviations** (errors) between the observed data points and the assumed line.

#### Assumed Model
```
y = a + bx
```

#### Least Squares Conditions
To minimize error:
```
∑(y − a − bx)² → minimum
```
This leads to the **normal equations**:
```
∑y = na + b∑x
∑xy = a∑x + b∑x²
```
#### Evaluation Process
The two normal equations are solved simultaneously to determine constants `a` and `b`.

#### Accuracy Considerations
- Simple and effective for linear trends  
- Not suitable for nonlinear data  
- Accuracy depends on data distribution

#### Applicability
- Widely used for trend analysis
- Useful when data follows an approximately linear pattern

### Least Square Regression Method For Linear Equations Code
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

    // Output 
    out << "Linear Fit Equation:\n";
    out << "y = " << a << " + " << b << "x\n";

    return 0;
}

```

### Least Square Regression Method For Linear Equations Input
```
5
1 50
2 80
3 96
4 120
5 45
```

### Least Square Regression Method For Linear Equations Output
```
Linear Fit Equation:
y = 69.2 + 3x

```

---

### Least Square Regression Method For Transcendental Equations 

### Least Square Regression Method For Transcendental Equations Theory
#### Method used
**Least Squares Regression – Transcendental Equation**

#### Objective
To fit nonlinear relationships by transforming them into linear form so that least squares method can be applied.

#### Common Transcendental Forms
1. **Exponential model**
```
y = ae^(bx)
```

2. **Power model**
```
y = ax^b
```

#### Linearization Technique

##### Exponential Equation
Taking natural logarithm:
```
ln y = ln a + bx
```
Let:
```
Y = ln y , A = ln a
```
Then:
```
Y = A + bx
```

##### Power Equation
Taking logarithm on both sides:
```
log y = log a + b log x
```

Let:
```
Y = log y , X = log x , A = log a
```
Then:
```
Y = A + bX
```
#### Evaluation Process
- Transform the given data
- Apply linear least squares regression
- Compute constants
- Convert back to original form

#### Accuracy Considerations
- Transformation may amplify errors
- Requires positive data values
- Fit quality depends on correct model assumption

#### Applicability
- Used in population growth, decay processes, and empirical laws
- Suitable for nonlinear experimental data


### Least Square Regression Method For Transcendental Equations Code
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

    // Output
    out << "Transcendental (Exponential) Fit:\n";
    out << "y = " << a << " * e^(" << b << "x)\n";

    return 0;
}

```

### Least Square Regression Method For Transcendental Equations Input
```
5
1 50
2 80
3 96
4 120
5 45
```

### Least Square Regression Method For Transcendental Equations Output
```
Transcendental (Exponential) Fit:
y = 68.8608 * e^(0.0194744x)

```

---

### Least Square Regression Method For Polynomial Equations 

### Least Square Regression Method For Polynomial Equations Theory
#### Method used
**Least Squares Regression – Polynomial Equation**

#### Objective
To fit a polynomial curve when linear regression is insufficient to represent data trends.

#### Assumed Model (Second Order Polynomial)
```
y = a + bx + cx²
```

#### Data Requirement
A set of experimental observations:
```
(x1, y1), (x2, y2), ..., (xn, yn)
```

#### Core Idea
The coefficients of the polynomial are determined by minimizing the sum of squared deviations between observed and computed values.

#### Normal Equations
For a second-degree polynomial:
```
∑y = na + b∑x + c∑x²
∑xy = a∑x + b∑x² + c∑x³
∑x²y = a∑x² + b∑x³ + c∑x⁴
```
#### Evaluation Process
- Compute required summations from data table
- Solve the system of equations
- Substitute coefficients into polynomial

#### Accuracy Considerations
- Higher degree improves fit but may cause overfitting
- Computational complexity increases with degree
- Sensitive to data errors

#### Applicability
- Used when data shows curvature
- Suitable for engineering and experimental modeling

### Least Square Regression Method For Polynomial Equations Code
```cpp
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
            return {};  
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

    // File Handling
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


```

### Least Square Regression Method For Polynomial Equations Input
```
5
1 6 
2 11
3 18
4 27
5 38
```

### Least Square Regression Method For Polynomial Equations Output
```
Quadratic Polynomial Fit:
y = 3 + 2x + 1x^2

```

---

### Solution of Differential Equations

### Equal Interval Interpolation Method

### Equal Interval Interpolation Theory
[Add your theory content here]

### Equal Interval Interpolation Code
```cpp
# Add your code here
```

### Equal Interval Interpolation Input
```
[Add your input format here]
```

### Equal Interval Interpolation Output
```
[Add your output format here]
```

---

### Second Order Derivative Method 

### Second Order Derivative Theory
[Add your theory content here]

### Second Order Derivative Code
```cpp
# Add your code here
```

### Second Order Derivative Input
```
[Add your input format here]
```

### Second Order Derivative Output
```
[Add your output format here]
```

---

### Runge Kutta Method 

### Runge Kutta Theory
[Add your theory content here]

### Runge Kutta Code
```cpp
# Add your code here
```

### Runge Kutta Input
```
[Add your input format here]
```

### Runge Kutta Output
```
[Add your output format here]
```

---

### Numerical Differentiation Method

### Numerical Differentiation Theory
[Add your theory content here]

### Numerical Differentiation Code
```cpp
# Add your code here
```

### Numerical Differentiation Input
```
[Add your input format here]
```

### Numerical Differentiation Output
```
[Add your output format here]
```

---

## Solution of Numerical Integrations

### Simpson's One-Third Rule

### Simpson's One-Third Rule Theory
[Add your theory content here]

### Simpson's One-Third Rule Code
```cpp
# Add your code here
```

### Simpson's One-Third Rule Input
```
[Add your input format here]
```

### Simpson's One-Third Rule Output
```
[Add your output format here]
```

---

### Simpson's Three-Eighths Rule 

### Simpson's Three-Eighths Rule Theory
[Add your theory content here]

### Simpson's Three-Eighths Rule Code
```cpp
# Add your code here
```

### Simpson's Three-Eighths Rule Input
```
[Add your input format here]
```

### Simpson's Three-Eighths Rule Output
```
[Add your output format here]
```

---


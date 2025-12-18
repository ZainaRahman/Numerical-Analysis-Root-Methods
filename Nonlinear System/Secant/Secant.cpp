#include <bits/stdc++.h>
using namespace std;

double f(double x, vector<double>& coef) {
    int n = coef.size() - 1;
    double s = 0;
    for (int i = 0; i <= n; i++)
        s += coef[i] * pow(x, n - i);
    return s;
}
void printEq(ofstream &out, vector<double>& coef) {
    int n = coef.size() - 1;
    out << "Equation: ";
    for (int i = 0; i <= n; i++) {
        out << coef[i] << "*x^" << (n - i);
        if (i != n) out << " + ";
    }
    out << endl;
}
int main() {
    ifstream in("input.txt");
    ofstream out("output.txt");

    int deg, maxit, it;
    double x0, x1, root, err;

    in >> deg;
    vector<double> coef(deg + 1);
    for (int i = 0; i <= deg; i++) in >> coef[i];
    in >> x0 >> x1 >> err >> maxit;
    printEq(out,coef);
    for (it = 1; it <= maxit; it++) {
        root = x1 - f(x1, coef) * (x1 - x0) / (f(x1, coef) - f(x0, coef));
        if (fabs(root - x1) < err || fabs(f(root, coef)) < err) break;
        x0 = x1;
        x1 = root;
    }
    out << "Secant Root = " << root << "\n";
    out << "Iterations = " << it << "\n";
    return 0;
}

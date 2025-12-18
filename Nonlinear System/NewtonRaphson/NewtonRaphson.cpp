#include <bits/stdc++.h>
using namespace std;

double f(double x, vector<double>& c) {
    int n = c.size() - 1;
    double s = 0;
    for (int i = 0; i <= n; i++)
        s += c[i] * pow(x, n - i);
    return s;
}

double df(double x, vector<double>& c) {
    int n = c.size() - 1;
    double s = 0;
    for (int i = 0; i < n; i++)
        s += c[i] * (n - i) * pow(x, n - i - 1);
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
    double x0, root, err;

    in >> deg;
    vector<double> coef(deg + 1);
    for (int i = 0; i <= deg; i++) in >> coef[i];
    in >> x0 >> err >> maxit;
    printEq(out,coef);
    for (it = 1; it <= maxit; it++) {
        if (fabs(df(x0, coef)) < 1e-12) {
            out << "Derivative too small\n";
            return 0;
        }
        root = x0 - f(x0, coef) / df(x0, coef);
        if (fabs(root - x0) < err || fabs(f(root, coef)) < err) break;
        x0 = root;
    }
    out << "Newton-Raphson Root = " << root << "\n";
    out << "Iterations = " << it << "\n";
    return 0;
}

#include <bits/stdc++.h>
using namespace std;

double f(double x, vector<double>& c) {
    int n = c.size() - 1;
    double s = 0;
    for (int i = 0; i <= n; i++)
        s += c[i] * pow(x, n - i);
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
    double lowl, upl, err, root, prev;

    in >> deg;
    vector<double> coef(deg + 1);
    for (int i = 0; i <= deg; i++) in >> coef[i];
    in >> lowl >> upl >> err >> maxit;
    printEq(out,coef);
    if (f(lowl, coef) * f(upl, coef) >= 0) {
        out << "Invalid interval\n";
        return 0;
    }

    for (it = 1; it <= maxit; it++) {
        prev = root;
        root = (lowl + upl) / 2.0;

        if (fabs(root - prev) < err || fabs(f(root, coef)) < err)
            break;

        f(lowl, coef) * f(root, coef) < 0 ? upl = root : lowl = root;
    }

    out << "Bisection Root = " << root << "\n";
    out << "Iterations = " << it << "\n";
    return 0;
}

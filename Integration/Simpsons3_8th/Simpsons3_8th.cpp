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

    int deg, n;
    double a, b, h, result = 0;

    in >> deg;
    vector<double> coef(deg + 1);
    for (int i = 0; i <= deg; i++) in >> coef[i];

    in >> a >> b >> n;
    printEq(out,coef);
    if (n % 3 != 0) {
        out << "Number of subintervals must be multiple of 3\n";
        return 0;
    }

    h = (b - a) / n;

    result = f(a, coef) + f(b, coef);

    for (int i = 1; i < n; i++) {
        double x = a + i * h;
        if (i % 3 == 0)
            result += 2 * f(x, coef);
        else
            result += 3 * f(x, coef);
    }

    result *= 3 * h / 8.0;

    out << "Simpson 3/8 Rule Result = " << result << endl;
    return 0;
}

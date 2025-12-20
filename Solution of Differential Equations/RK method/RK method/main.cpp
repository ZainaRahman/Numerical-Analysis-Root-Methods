#include <iostream>
#include <fstream>
using namespace std;

double f(double x, double y)
{
    return x + y;
}

double rungeKutta(double x0, double y0, double h, double xn)
{
    double x = x0;
    double y = y0;

    while (x  < xn)
    {
        double k1 = h * f(x, y);
        double k2 = h * f(x + h/2.0, y + k1/2.0);
        double k3 = h * f(x + h/2.0, y + k2/2.0);
        double k4 = h * f(x + h, y + k3);

        y = y + (k1 + 2*k2 + 2*k3 + k4) / 6.0;
        x = x + h;
    }

    return y;
}

int main()
{
    ifstream fin("Input.txt");
    ofstream fout("Output.txt");

    if (!fin.is_open())
    {
        cout << "Error: Input.txt not found\n";
        return 0;
    }

    double x0, y0, h, xn;
    fin >> x0 >> y0 >> h >> xn;

    double result = rungeKutta(x0, y0, h, xn);

    cout << "Runge-Kutta 4th Order Method\n";
    cout << "----------------------------\n";
    cout << "Differential equation: dy/dx = x + y\n";
    cout << "Initial condition: y(" << x0 << ") = " << y0 << "\n";
    cout << "Step size (h): " << h << "\n";
    cout << "Required value at x = " << xn << "\n";
    cout << "Computed solution: y(" << xn << ") = " << result << "\n";

    fout << "Runge-Kutta 4th Order Method\n";
    fout << "----------------------------\n";
    fout << "Differential equation: dy/dx = x + y\n";
    fout << "Initial condition: y(" << x0 << ") = " << y0 << "\n";
    fout << "Step size (h): " << h << "\n";
    fout << "Required value at x = " << xn << "\n";
    fout << "Computed solution: y(" << xn << ") = " << result << "\n";

    fin.close();
    fout.close();

    return 0;
}

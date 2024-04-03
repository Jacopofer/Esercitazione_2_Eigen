#include <iostream>
#include "Eigen/Eigen"
#include <iomanip>

using namespace std;
using namespace Eigen;

int main()
{
    MatrixXd A1(2,2);
    VectorXd b1(2);
    A1 << 5.547001962252291e-01, -3.770900990025203e-02,
        8.320502943378437e-01, -9.992887623566787e-01;
    b1 << -5.169911863249772e-01, 1.672384680188350e-01;

    MatrixXd A2(2,2);
    VectorXd b2(2);
    A2 << 5.547001962252291e-01, -5.540607316466765e-01,
        8.320502943378437e-01, -8.324762492991313e-01;
    b2 << -6.394645785530173e-04, 4.25954962877223e-04;

    MatrixXd A3(2,2);
    VectorXd b3(2);
    A3 << 5.547001962252291e-01, -5.547001955851905e-01,
        8.320502943378437e-01, -8.320502947645361e-01;
    b3 << -6.400391328043042e-10, 4.266924591433963e-10;

    VectorXd exactSolution = -VectorXd::Ones(2);

    //METODO QR

    VectorXd x1(2);
    VectorXd x2(2);
    VectorXd x3(2);

    x1 = A1.householderQr().solve(b1);
    x2 = A2.householderQr().solve(b2);
    x3 = A3.colPivHouseholderQr().solve(b3);

    double rel_err1 = (exactSolution - x1).norm() / exactSolution.norm();
    double rel_err2 = (exactSolution - x2).norm() / exactSolution.norm();
    double rel_err3 = (exactSolution - x3).norm() / exactSolution.norm();

    cout << "QR:"<< endl;

    cout << "Soluzione del primo sistema: " << scientific << setprecision(1) << x1.transpose() << endl;
    cout << "Errrore relativo del primo sistema: " << rel_err1 << endl;

    cout << "Soluzione del secondo sistema: " << scientific << setprecision(1) << x2.transpose() << endl;
    cout << "Errrore relativo del secondo sistema: " << rel_err2 << endl;

    cout << "Soluzione del terzo sistema: " << scientific << setprecision(1) << x3.transpose() << endl;
    cout << "Errrore relativo del terzo sistema: " << rel_err3 << endl;

    //METODO PALU

    VectorXd y1(2);
    VectorXd y2(2);
    VectorXd y3(2);

    y1 = A1.fullPivLu().solve(b1);
    y2 = A2.fullPivLu().solve(b2);
    y3 = A3.fullPivLu().solve(b3);

    double err_rel1 = (exactSolution - y1).norm() / exactSolution.norm();
    double err_rel2 = (exactSolution - y2).norm() / exactSolution.norm();
    double err_rel3 = (exactSolution - y3).norm() / exactSolution.norm();

    cout << "PALU:" << endl;

    cout << "Soluzione del primo sistema: " << scientific << setprecision(1) << y1.transpose() << endl;
    cout << "Errrore relativo del primo sistema: " << err_rel1 << endl;

    cout << "Soluzione del secondo sistema: " << scientific << setprecision(1) << y2.transpose() << endl;
    cout << "Errrore relativo del secondo sistema: " << err_rel2 << endl;

    cout << "Soluzione del terzo sistema: " << scientific << setprecision(1) << y3.transpose() << endl;
    cout << "Errrore relativo del terzo sistema: " << err_rel3 << endl;

    return 0;
}

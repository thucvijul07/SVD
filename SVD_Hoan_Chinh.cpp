#include <iostream>
#include <Eigen/Dense>
#include <iomanip>
#include <math.h>

using namespace std;
using namespace Eigen;

void Swap(double &a, double &b)
{
    double temp = a;
    a = b;
    b = temp;
}
void Show(double a[][10], int row, int col)
{
    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
            cout << setw(9) << setprecision(4) << a[i][j];
        cout << endl;
    }
}

void ChuyenVi(double a[][10], double at[][10], int row, int col)
{
    for (int i = 0; i < row; i++)
        for (int j = 0; j < col; j++)
            at[j][i] = a[i][j];
}

void NhanMaTran(MatrixXd &x, double a[][10], double b[][10], int row1, int col1, int col2)
{
    for (int i = 0; i < row1; i++)
        for (int j = 0; j < col2; j++)
        {
            x(i, j) = 0;
            for (int k = 0; k < col1; k++)
                x(i, j) = x(i, j) + a[i][k] * b[k][j];
        }
}

void GetEigenValuesAndVector(MatrixXd x, MatrixXd &lambda, MatrixXd &vector)
{
    SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(x);
    vector = eigensolver.eigenvectors();
    lambda = eigensolver.eigenvalues();
    int k = lambda.rows(), l = vector.rows();
    for (int i = 0; i < k; i++)
        if (lambda(i, 0) < 0.000001)
            lambda(i, 0) = 0;
    // dua ve lambda tu lon den be
    for (int i = 0; i < k; i++)
        for (int j = i + 1; j < k; j++)
        {
            if (lambda(j, 0) > lambda(i, 0))
            {
                Swap(lambda(j, 0), lambda(i, 0));
                for (int h = 0; h < l; h++)
                    Swap(vector(h, i), vector(h, j));
            }
        }
}

void MaTranS(MatrixXd lambda, double sigma[][10], int rows, int cols)
{
    int k = 0;
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
        {
            sigma[i][j] = (i != j) ? 0 : sqrt(lambda(k, 0));
            if (i == j)
                k++;
        }
}

void MaTranU(MatrixXd lambda, MatrixXd vector, double U[][10], double a[][10], int rows, int cols)
{
    MatrixXd x(rows, 1);
    double Vi[cols][10];

    for (int i = 0; i < cols; i++)
    { // i la so vector
        for (int j = 0; j < cols; j++)
            Vi[j][0] = vector(j, i); // lay ra tung cot cua ma tran vector rieng
        NhanMaTran(x, a, Vi, rows, cols, 1);
        for (int k = 0; k < rows; k++)
            if (lambda(i, 0) != 0)
                U[k][i] = (1 / sqrt(lambda(i, 0))) * x(k, 0);
            else
                k = rows;
    }

    for (int i = 0; i < rows; i++)
    {
        if (lambda(i, 0) == 0)
        {
            MatrixXd A(rows - 1, rows);
            VectorXd B(rows - 1, 1), x1(rows, 1);
            for (int k = 0; k < rows; k++)
                B(k, 0) = 0;
            for (int j = 0; j < rows - 1; j++)
                for (int h = 0; h < rows; h++)
                    if (j != i)
                        A(j, h) = U[h][j];
            x1 = A.colPivHouseholderQr().solve(B);
            for (int k = 0; k < rows; k++)
                U[k][i] = x1(k, 0);
        }
    }
}
void MaTranV(MatrixXd vector, double V[][10])
{
    int row = vector.rows(), col = vector.cols();
    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        {
            V[i][j] = vector(i, j);
        }
    }
}

int main()
{
    int rows, cols;
    cout << "Nhap so hang va so cot " << endl;
    double a[10][10], at[10][10];
    cin >> rows >> cols;
    cout << "Nhap ma tran A" << endl;
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
        {
            cout << "a[" << i + 1 << "][" << j + 1 << "] = ";
            cin >> a[i][j];
        }
    ChuyenVi(a, at, rows, cols);
    MatrixXd x(cols, cols), lambda(cols, 1), vector(cols, cols);
    NhanMaTran(x, at, a, cols, rows, cols);
    GetEigenValuesAndVector(x, lambda, vector);
    double sigma[10][10], U[10][10], V[10][10];
    MaTranU(lambda, vector, U, a, rows, cols);
    cout << "Ma tran U" << endl;
    Show(U, rows, rows);
    MaTranS(lambda, sigma, rows, cols);
    cout << "Ma tran sigma" << endl;
    Show(sigma, rows, cols);
    MaTranV(vector, V);
    cout << "Ma tran VT" << endl;
    double VT[10][10];
    ChuyenVi(V, VT, cols, cols);
    Show(VT, cols, cols);
    
}

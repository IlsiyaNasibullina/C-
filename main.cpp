//Ilsiia Nasibullina CS-05
#include<iostream>
#include <utility>
#include <vector>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <cmath>

#define GNUPLOT_NAME "gnuplot -persist"

using namespace std;
typedef long long ll;
typedef vector<ll> vll;
typedef long double ld;

class ColumnVector {
public:
    int a;
    vector<double> myColumnVector;

    ColumnVector() {};

    ColumnVector(int n) {
        a = n;
        for (int i = 0; i < a; i++) {
            myColumnVector.push_back(0);
        }
    }

    void generateColumnVector(int n) {
        for (int i = 0; i < n; i++) {
            myColumnVector.push_back(0);
        }
    }

    void generateColumnVector(vector<pair<double, double>> input) {
        a = input.size();
        for (int i = 0; i < a; i++) {
            myColumnVector.push_back(0);
        }
        for (int i = 0; i < a; i++) {
            myColumnVector[i] = input[i].second;
        }
    }

    void operator=(ColumnVector A) {
        for (int i = 0; i < a; i++) {
            myColumnVector[i] = A.myColumnVector[i];
        }
    }

    // counting the norm
    double norm() {
        double countingNorm = 0;
        for (int i = 0; i < a; i++) {
            countingNorm += myColumnVector[i] * myColumnVector[i];
        }
        countingNorm = sqrt(countingNorm);
        return countingNorm;
    }

    ColumnVector operator+(ColumnVector A) {
        ColumnVector B(a);
        for (int i = 0; i < a; i++) {
            B.myColumnVector[i] = myColumnVector[i] + A.myColumnVector[i];
        }
        return B;
    }

    ColumnVector operator-(ColumnVector A) {
        ColumnVector B(a);
        for (int i = 0; i < a; i++) {
            B.myColumnVector[i] = myColumnVector[i] - A.myColumnVector[i];
        }
        return B;
    }

    ColumnVector operator*(ColumnVector A) {
        ColumnVector B(a);
        for (int i = 0; i < a; i++) {
            B.myColumnVector[i] = myColumnVector[i] * A.myColumnVector[i];
        }
        return B;
    }

    friend ostream &operator<<(ostream &cout, ColumnVector &A) {
        for (int i = 0; i < A.a; i++) {
            if (abs(A.myColumnVector[i]) < 1e-10) {
                double val = 0;
                cout << fixed << setprecision(4) << val << "\n";
            } else cout << fixed << setprecision(4) << A.myColumnVector[i] << "\n";
        }
        return cout;
    }

    friend istream &operator>>(istream &cin, ColumnVector &A) {
        cin >> A.a;
        A.generateColumnVector(A.a);
        int value;
        for (int i = 0; i < A.a; i++) {
            cin >> value;
            A.myColumnVector[i] = value;
        }
        return cin;
    }

};

class Matrix {
public:
    int rows;
    int columns;
    vector<vector<double>> myMatrix;

    Matrix() {}

    Matrix(int n, int m) {
        rows = n;
        columns = m;
        vector<double> temp(columns, 0);
        for (int i = 0; i < rows; i++) {
            myMatrix.push_back(temp);
        }
    }

    void generateMatrix(int r, int c) {
        vector<double> temp(c, 0);
        for (int i = 0; i < r; i++) {
            myMatrix.push_back(temp);
        }
    }

    void generateMatrix(Matrix A) {
        rows = A.rows;
        columns = A.columns / 2;
        vector<double> temp(columns, 0);
        for (int i = 0; i < rows; i++) {
            myMatrix.push_back(temp);
        }
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                myMatrix[i][j] = A.myMatrix[i][j + A.rows];
            }
        }
    }


    void generateMatrix(Matrix A, Matrix B) {
        rows = A.rows;
        columns = A.columns + B.columns;
        vector<double> temp(columns, 0);
        for (int i = 0; i < rows; i++) {
            myMatrix.push_back(temp);
        }
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < A.columns; j++) {
                myMatrix[i][j] = A.myMatrix[i][j];
            }
            for (int j = A.columns; j < columns; j++) {
                myMatrix[i][j] = B.myMatrix[i][j - A.columns];
            }
        }
    }

    void generateMatrix(ColumnVector b) {
        rows = b.a;
        columns = 1;
        vector<double> temp(columns, 0);
        for (int i = 0; i < rows; i++) {
            myMatrix.push_back(temp);
        }
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                myMatrix[i][j] = b.myColumnVector[i];
            }
        }
    }

    void generateMatrix(vector<pair<double, double>> input, int n) {
        columns = n + 1;
        rows = input.size();
        vector<double> temp(columns, 0);
        for (int i = 0; i < rows; i++) {
            myMatrix.push_back(temp);
        }
        for (int i = 0; i <= n; i++) {
            for (int j = 0; j < input.size(); j++) {
                myMatrix[j][i] = pow(input[j].first, i);
            }
        }
    }

    void operator=(Matrix A) {
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                myMatrix[i][j] = A.myMatrix[i][j];
            }
        }
    }

    Matrix operator+(Matrix A) {
        Matrix B(rows, columns);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                B.myMatrix[i][j] = myMatrix[i][j] + A.myMatrix[i][j];
            }
        }
        return B;
    }

    Matrix operator-(Matrix A) {
        Matrix B(rows, columns);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                B.myMatrix[i][j] = myMatrix[i][j] - A.myMatrix[i][j];
            }
        }
        return B;
    }

    Matrix operator*(Matrix A) {
        Matrix B(rows, A.columns);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < A.columns; j++) {
                for (int k = 0; k < columns; k++) {
                    B.myMatrix[i][j] += myMatrix[i][k] * A.myMatrix[k][j];
                }
            }
        }
        return B;
    }

    int indexOfMaxPivot(int i) {
        double pivot = myMatrix[i][i];
        int indexOfMaxPivot = i;
        for (int j = i; j < rows; j++) {
            if (abs(myMatrix[j][i]) > abs(pivot)) {
                pivot = myMatrix[j][i];
                indexOfMaxPivot = j;
            }
        }
        return indexOfMaxPivot;
    }

    bool allZeros(int i, bool upper) {
        if (upper) {
            for (int j = i + 1; j < rows; j++) {
                if (abs(myMatrix[j][i]) < 1e-3) continue;
                else return false;
            }
        } else {
            for (int j = i - 1; j >= 0; j--) {
                if (abs(myMatrix[j][i]) < 1e-3) continue;
                else return false;
            }
        }
        return true;
    }

    void inverse() {
        for (int i = 0; i < rows - 1; i++) { //pivot
            if (i != indexOfMaxPivot(i)) {
                swap(myMatrix[i], myMatrix[indexOfMaxPivot(i)]);
            }
            if (allZeros(i, true)) {
                continue;
            } else {
                for (int j = i + 1; j < rows; j++) { // rows for reduction
                    double ratio = myMatrix[j][i] / myMatrix[i][i];
                    for (int k = 0; k < columns; k++) { // reducing a row
                        myMatrix[j][k] = myMatrix[j][k] - ratio * myMatrix[i][k];
                    }
                }
            }
        }
        for (int i = rows - 1; i >= 0; i--) { //pivot
            if (allZeros(i, false)) {
                continue;
            } else {
                for (int j = i - 1; j >= 0; j--) { // rows for reduction
                    double ratio = myMatrix[j][i] / myMatrix[i][i];
                    for (int k = 0; k < columns; k++) { // reducing a row
                        myMatrix[j][k] = myMatrix[j][k] - ratio * myMatrix[i][k];
                    }
                }
            }
        }
        for (int i = 0; i < rows; i++) {
            double ratio = myMatrix[i][i];
            for (int j = 0; j < columns; j++) {
                if (ratio != 0) myMatrix[i][j] = myMatrix[i][j] / ratio;
            }
        }
    }

    Matrix transpose() {
        Matrix B(columns, rows);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                B.myMatrix[j][i] = myMatrix[i][j];
            }
        }
        return B;
    }

    friend ostream &operator<<(ostream &cout, Matrix &A) {
        for (int i = 0; i < A.rows; i++) {
            for (int j = 0; j < A.columns; j++) {
                if (abs(A.myMatrix[i][j]) < 1e-10) {
                    double val = 0;
                    cout << fixed << setprecision(4) << val << " ";
                } else cout << fixed << setprecision(4) << A.myMatrix[i][j] << " ";
            }
            cout << "\n";
        }
        return cout;
    }

    friend istream &operator>>(istream &cin, Matrix &A) {
        A.columns = 2;
        cin >> A.rows;
        A.generateMatrix(A.rows, A.columns);
        int value;
        for (int i = 0; i < A.rows; i++) {
            for (int j = 0; j < A.columns; j++) {
                cin >> value;
                A.myMatrix[i][j] = value;
            }
        }
        return cin;
    }
};

class SquareMatrix : public Matrix {
public:
    SquareMatrix() = default;

    SquareMatrix(int n) : Matrix(n, n) {
    }

    void generateSquareMatrix(int n) {
        vector<double> temp(n, 0);
        for (int i = 0; i < n; i++) {
            myMatrix.push_back(temp);
        }
    }


    void operator=(SquareMatrix A) {
        for (int i = 0; i < A.rows; i++) {
            for (int j = 0; j < A.columns; j++) {
                this->myMatrix[i][j] = A.myMatrix[i][j];
            }
        }
    }

    SquareMatrix operator+(SquareMatrix &A) {
        Matrix *first = this;
        Matrix *second = &A;
        Matrix result = *first + *second;
        SquareMatrix *answer = (SquareMatrix *) &result;
        return *answer;
    }

    SquareMatrix operator-(SquareMatrix &A) {
        Matrix *first = this;
        Matrix *second = &A;
        Matrix result = *first - *second;
        SquareMatrix *answer = (SquareMatrix *) &result;
        return *answer;
    }

    SquareMatrix operator*(SquareMatrix &A) {
        Matrix *first = this; // upcasting
        Matrix *second = &A;
        Matrix result = *first * *second;
        SquareMatrix *answer = (SquareMatrix *) &result; // downcasting
        return *answer;
    }

    SquareMatrix transpose() {
        Matrix *first = this;
        Matrix second = (*first).transpose();
        SquareMatrix *result = (SquareMatrix *) &second;
        return *result;
    }

    friend ostream &operator<<(ostream &cout, SquareMatrix &A) {
        for (int i = 0; i < A.rows; i++) {
            for (int j = 0; j < A.columns; j++) {
                if (abs(A.myMatrix[i][j]) < 1e-10) {
                    double val = 0;
                    cout << fixed << setprecision(4) << val << " ";
                } else cout << fixed << setprecision(4) << A.myMatrix[i][j] << " ";
            }
            cout << "\n";
        }
        return cout;
    }

    friend istream &operator>>(istream &cin, SquareMatrix &A) {
        int n;
        cin >> n;
        A.rows = n;
        A.columns = n;
        A.generateSquareMatrix(n);
        int value;
        for (int i = 0; i < A.rows; i++) {
            for (int j = 0; j < A.columns; j++) {
                cin >> value;
                A.myMatrix[i][j] = value;
            }
        }
        return cin;
    }
};

class IdentityMatrix : public SquareMatrix {
public:
    IdentityMatrix() = default;

    IdentityMatrix(int n) {
        rows = n;
        columns = n;
        vector<double> temp(columns, 0);
        for (int i = 0; i < n; i++) {
            myMatrix.push_back(temp);
        }
        for (int i = 0; i < n; i++) {
            myMatrix[i][i] = 1;
        }
    }

};

Matrix LeastSquareApproximation(vector<pair<double, double>> input, int n) {
    ColumnVector b;
    b.generateColumnVector(input);

    Matrix A;
    A.generateMatrix(input, n);
    cout << "A:\n";
    cout << A;

    cout << "A_T*A:\n";
    Matrix C = A.transpose() * A;
    cout << C;

    Matrix B;
    B.generateMatrix(C, IdentityMatrix(C.rows));
    B.inverse();

    cout << "(A_T*A)^-1:\n";
    Matrix Inverse;
    Inverse.generateMatrix(B);
    cout << Inverse;

    cout << "A_T*b:\n";
    Matrix G = A.transpose();
    Matrix E;
    E.generateMatrix(b);
    Matrix F = G * E;
    cout << F;


    cout << "x~:\n";
    Matrix answer = Inverse * F;
    cout << answer;
    return answer;
}

int main() {
    int m, n;
    cin >> m;
    vector<pair<double, double>> input;
    for (int i = 0; i < m; i++) {
        double value1, value2;
        cin >> value1 >> value2;
        input.push_back({value1, value2});
    }
    cin >> n;

    Matrix answer = LeastSquareApproximation(input, n);

    string formula = "f(x) = ";
    for (int i = 0; i < answer.rows - 1; i++) {
        formula += to_string(round(answer.myMatrix[i][0] * 10000) / 10000)
                + " * x**" + to_string(i) + " + ";
    }
    formula += to_string(round(answer.myMatrix[answer.rows - 1][0] * 10000) / 10000)
            + " * x**" + to_string(answer.rows - 1);

    double minX = INFINITY;
    double maxX = -INFINITY;
    double minY = INFINITY;
    double maxY = -INFINITY;
    for (pair<double, double> p: input) {
        if (p.first < minX) {
            minX = p.first;
        }
        if (p.first > maxX) {
            maxX = p.first;
        }
        if (p.second < minY) {
            minY = p.second;
        }
        if (p.second > maxY) {
            maxY = p.second;
        }
    }
    FILE *pipe = popen(GNUPLOT_NAME, "w");
    if (pipe != nullptr) {
        ::fprintf(pipe, "set xrange [%f:%f]\n", minX - 10.0, maxX + 10.0);
        ::fprintf(pipe, "set yrange [%f:%f]\n", minY - 10.0, maxY + 10.0);
        ::fprintf(pipe, "%s\n",
                  formula.c_str());

        ::fprintf(pipe, "%s\n",
                  "plot '-' using 1:2 title 'Points' with points, "
                  "f(x) title 'y = 2.6525 + 0.9351*x + "
                  "-0.1165*x^2 + 0.0083*x^3 + -0.0002*x^4' with lines");
        for (int i = 0; i < input.size(); i++) {
            double x = input[i].first;
            double y1 = input[i].second;
            fprintf(pipe, "%f\t%f\n", x, y1);
        }
        ::fprintf(pipe, "%s\n", "e");
        ::fflush(pipe);
        pclose(pipe);
    }
    return 0;

}



#include "octave/oct.h"
#include "octave/parse.h"
#include <ctime>
#include <cstdio>

DEFUN_DLD(fast_sor, args, nargout, "C++ implementation of the SOR method")
{
    octave_value_list retval;

    Matrix A = args(0).matrix_value();
    ColumnVector b = args(1).column_vector_value();
    double tol = args(2).double_value();
    int nMaxIter = args(3).int_value();
    double w = args(4).double_value();
    int n = A.rows();

    ColumnVector x0(n);
    x0.fill(0.0f);

    ColumnVector x(n);
    x = x0;

    ColumnVector er(n);
    er.fill(0.0f);
    er(0) = 1.0f;

    octave_idx_type i = 0;

    float timeStart = (float)clock() / (float)CLOCKS_PER_SEC;

    while (er(i) > tol && i < nMaxIter)
    {
        for (octave_idx_type j = 0; j < n; ++j)
        {
            float soma = 0.0f;

            for (octave_idx_type k = 0; k < j - 1; ++k)
            {
                soma += A(j, k) * x(k);
            }

            for (octave_idx_type k = j + 1; k < n; ++k)
            {
                soma += A(j, k) * x0(k);
            }

            if (A(i, j) != 0)
                x(j) = w * (b(j) - soma) / A(j, j) + (1 - w) * x0(j);
            else
                x(j) = (1 - w) * x0(j);
        }

        ++i;

        octave_value_list a_in;
        a_in(0) = (x - x0);
        a_in(1) = "inf";
        double a = (octave::feval("norm", a_in))(0).double_value();

        octave_value_list b_in;
        b_in(0) = x;
        b_in(1) = "inf";
        double b = (octave::feval("norm", b_in))(0).double_value();

        er(i) = b != 0.0f ? a / b : 0.0f;
        x0 = x;
    }

    float timeEnd = (float)clock() / (float)CLOCKS_PER_SEC;

    retval(0) = x;
    retval(1) = er;
    retval(2) = i;
    retval(3) = timeEnd - timeStart;

    return retval;
}

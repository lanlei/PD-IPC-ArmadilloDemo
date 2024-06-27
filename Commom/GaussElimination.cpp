#include "GaussElimination.h"
#include <iostream>

void gaussElimination(MatrixX& A, VectorX& b, VectorX& x)
{
	int dim = x.size();
	MatrixX extend_A(dim, dim + 1);
	extend_A.setZero();
	extend_A.block(0, 0, A.rows(), A.cols()) = A;
	int nnz_rows = A.rows();
	VectorX extend_b(dim);
	extend_b.setZero();
	extend_b.head(b.size()) = b;
	extend_A.col(dim) = extend_b;
	MatrixX Q(dim, dim);
	Q.setIdentity();
	std::cout << extend_A.block(0, 0, dim, dim + 1) << std::endl << std::endl;
	
	for (int i = 0; i < nnz_rows; i++)
	{	
		int j = i;

		if (extend_A.row(i).isZero())
		{			
			bool has_nnz = false;
			for (int j = i + 1; j < nnz_rows; j++)
			{
				if (extend_A.row(j).isZero())
					continue;
				else
				{
					has_nnz = true;
					VectorX ri = extend_A.row(i);
					VectorX rj = extend_A.row(j);
					extend_A.row(i) = rj;
					extend_A.row(j) = ri;
					break;
				}
			}
			if (!has_nnz)
				break;
		}

		for (int c = 0; c < i; c++)
		{
			qeal v_ic = extend_A(i, c);
			if (v_ic != 0)
			{
				for (int k = 0; k < extend_A.cols(); k++)
					extend_A(i, k) = extend_A(i, k) - v_ic * extend_A(c, k);
			}
		}
		qeal vii = extend_A(i, i);
		if (vii != 0)
		{
			for (int k = 0; k < extend_A.cols(); k++)
				extend_A(i, k) /= vii;
		}
		else
		{
			qeal v_ik;
			int nk = -1;
			for (int k = i + 1; k < extend_A.cols() - 1; k++)
			{
				v_ik = extend_A(i, k);
				if (v_ik != 0)
				{
					nk = k;
					break;
				}
			}
			if (nk != -1)
			{
				VectorX ci = extend_A.col(i);
				VectorX ck = extend_A.col(nk);

				extend_A.col(i) = ck;
				extend_A.col(nk) = ci;

				VectorX ri = Q.row(i);
				VectorX rk =Q.row(nk);
				Q.row(i) = rk;
				Q.row(nk) = ri;

				for (int k = 0; k < extend_A.cols(); k++)
					extend_A(i, k) /= v_ik;
			}
		}

	}

	std::cout << Q << std::endl << std::endl;
	MatrixX R = extend_A.block(0, 0, dim, dim);
	b = extend_A.col(dim);
	std::cout << R << std::endl << std::endl;



	nnz_rows = 0;
	for (int i = 0; i < R.rows(); i++)
	{
		if (R.row(i).isZero())
			continue;
		nnz_rows++;
	}

	MatrixX sub_R = R.block(0, 0, nnz_rows, R.cols());
	std::cout << sub_R << std::endl << std::endl;

	MatrixX R1 = R.block(0, 0, nnz_rows, nnz_rows);
	MatrixX R2 = R.block(0, nnz_rows, nnz_rows, R.cols() - nnz_rows);
	std::cout << R1 << std::endl << std::endl;
	std::cout << R2 << std::endl << std::endl;

	MatrixX U(R.rows(), R2.cols());
	U.block(0, 0, R2.rows(), R2.cols()) = -1 * R2;

	MatrixX I2(R2.cols(), R2.cols());
	I2.setIdentity();
	U.block(R1.rows(), 0, R2.cols(), R2.cols()).setIdentity();
	std::cout << U << std::endl << std::endl;
	
	int d = 0;

}

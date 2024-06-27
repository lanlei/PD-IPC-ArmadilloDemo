#pragma once
#include "MatrixCore.h"

namespace PD_IPC
{
	class TetConstraints
	{
	public:
		Matrix3 rest;
		qeal w;
		Vector4i _idI;
		int _idO;

		qeal svclamp(qeal v, qeal vMin, qeal vMax) {
			qeal result = v > vMin ? v : vMin;
			return result > vMax ? vMax : result;
		}

		TetConstraints(const Vector4i& idI, qeal stiffness, const MatrixX& points)
		{
			_idI = idI;
			Matrix3 edges;
			for (int i = 0; i < 3; ++i) 
				edges.col(i) = points.col(idI[i]) - points.col(idI[3]);
			rest = edges.inverse();

			qeal V = (edges).determinant() / 6.0f;
			w = stiffness;
			w *= std::sqrt(std::abs(V));
			
		}

		void project(const MatrixX& points, MatrixX &projections, qeal& e)
		{
			Matrix3 edges;
			for (int i = 0; i < 3; ++i) 
				edges.col(i) = points.col(_idI[i]) - points.col(_idI[3]);
			Matrix3 F = edges * rest;
			Eigen::JacobiSVD<Matrix3> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
			//Vector3 S = svd.singularValues();
			//S(0) = svclamp(S(0), 0.95, 1.05);
			//S(1) = svclamp(S(1), 0.95, 1.05);
			//S(2) = svclamp(S(2), 0.95, 1.05);
			//if (svd.matrixU().determinant()*svd.matrixV().determinant() < 0.0f) S(2) = -S(2);
			//F = svd.matrixU() * S.asDiagonal() * svd.matrixV().transpose();
			//projections.block<3, 3>(0, _idO) = w * F;

			Vector3 S = svd.singularValues();
			qeal U[9], V[9];
			for (int k = 0; k < 9; k++)
			{
				U[k] = svd.matrixU().data()[k];
				V[k] = svd.matrixV().data()[k];
			}
			qeal det_F = S[0] * S[1] * S[2];
			if (det_F < 0)
			{
				S[2] *= -1;
				qeal high_val = S[0];
				qeal mid_val = S[1];
				qeal low_val = S[2];

				if (mid_val < low_val) {
					qeal temp = low_val;
					low_val = mid_val;
					mid_val = temp;
				}

				if (high_val < low_val) {
					qeal temp = low_val;
					low_val = high_val;
					high_val = temp;
				}
				if (high_val < mid_val) {
					qeal temp = mid_val;
					mid_val = high_val;
					high_val = temp;
				}

				S[0] = high_val;
				S[1] = mid_val;
				S[2] = low_val;
			}

			qeal sigma_new[3];
			sigma_new[0] = svclamp(S(0), 0.95, 1.05);
			sigma_new[1] = svclamp(S(1), 0.95, 1.05);
			sigma_new[2] = svclamp(S(2), 0.95, 1.05);

			e = (sigma_new[0] - S(0)) *  (sigma_new[0] - S(0)) +
				(sigma_new[1] - S(1)) *  (sigma_new[1] - S(1)) +
				(sigma_new[2] - S(2)) *  (sigma_new[2] - S(2));
			e *= 0.5 * w;

			qeal SVt[9]; // = new_sigma * vt
			SVt[0] = sigma_new[0] * V[0];
			SVt[1] = sigma_new[1] * V[3];
			SVt[2] = sigma_new[2] * V[6];

			SVt[3] = sigma_new[0] * V[1];
			SVt[4] = sigma_new[1] * V[4];
			SVt[5] = sigma_new[2] * V[7];

			SVt[6] = sigma_new[0] * V[2];
			SVt[7] = sigma_new[1] * V[5];
			SVt[8] = sigma_new[2] * V[8];

			qeal F_new[9]; // u * s * vt
			F_new[0] = (U[0] * SVt[0] + U[3] * SVt[1] + U[6] * SVt[2]);
			F_new[1] = (U[1] * SVt[0] + U[4] * SVt[1] + U[7] * SVt[2]);
			F_new[2] = (U[2] * SVt[0] + U[5] * SVt[1] + U[8] * SVt[2]);

			F_new[3] = (U[0] * SVt[3] + U[3] * SVt[4] + U[6] * SVt[5]);
			F_new[4] = (U[1] * SVt[3] + U[4] * SVt[4] + U[7] * SVt[5]);
			F_new[5] = (U[2] * SVt[3] + U[5] * SVt[4] + U[8] * SVt[5]);

			F_new[6] = (U[0] * SVt[6] + U[3] * SVt[7] + U[6] * SVt[8]);
			F_new[7] = (U[1] * SVt[6] + U[4] * SVt[7] + U[7] * SVt[8]);
			F_new[8] = (U[2] * SVt[6] + U[5] * SVt[7] + U[8] * SVt[8]);

			Matrix3 FF;
			for (int k = 0; k < 9; k++)
				FF.data()[k] = F_new[k];
			
			projections.block<3, 3>(0, _idO) = w * FF;
		}

		void addConstraint(std::vector<TripletX> &triplets, int& idO)
		{
	
			_idO = idO;
			int n = 3;
			for (int i = 0; i < n; ++i) {
				triplets.push_back(TripletX(_idO + i, _idI[0], w * rest(0, i)));
				triplets.push_back(TripletX(_idO + i, _idI[1], w * rest(1, i)));
				triplets.push_back(TripletX(_idO + i, _idI[2], w * rest(2, i)));
				triplets.push_back(TripletX(_idO + i, _idI[3], -w * (rest(0, i) + rest(1, i) + rest(2, i))));
			}
			idO += n;
		}



	};



}
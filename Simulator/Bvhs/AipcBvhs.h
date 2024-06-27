#pragma once
#ifndef AIPC_BVHS_H
#define AIPC_BVHS_H
#include <BVH.hpp>
#include <Morton.hpp>
#include "Model\GeometryElement.h"
#include "MatrixCore.h"
namespace AIPC
{
	bool box_box_intersection(const Vector3 &min1, const Vector3 &max1,
		const Vector3 &min2, const Vector3 &max2);


	bool box_box_intersection(const Vector3 &min1, const Vector3 &max1,
		const Vector3 &min2, const Vector3 &max2, Vector3 &overlap_min, Vector3 &overlap_max);
	
	class AipcBvhs
	{
	//private:
	public:
		std::vector<std::array<Vector3, 2>> boxlist;
		std::vector<int> new2old;
		std::vector<int> oldboxId;
		size_t n_corners = -1;
		std::vector<int> leafNodeFlag;
		std::vector<int> nodeIdToElementId;

		int maxLevelId = 0;
		std::vector<int> nodeLevelId;
		std::vector<std::vector<int>> nodeLevelGroup;

		void init_boxes_recursive(const std::vector<std::array<Vector3, 2>> &cornerlist, int node_index, int b, int e);

		void box_search_recursive(
			const Vector3 &bbd0, const Vector3 &bbd1,
			std::vector<unsigned int> &list,
			int n, int b, int e) const;

		static int max_node_index(int node_index, int b, int e);

		bool box_intersects_box(const Vector3 &bbd0, const Vector3 &bbd1, int index) const;

		bool is_initialized = false;

	public:
		void init(const std::vector<std::array<Vector3, 2>> &cornerlist);
		void init(const MatrixX &V, const MatrixXi &F, const qeal tol);

		int signBvhsLevel(int nodeId, int maxLevelId)
		{
			if (nodeIdToElementId[nodeId] == -2)
			{
				nodeLevelId[nodeId] = -1;
				return maxLevelId;
			}
			nodeLevelId[nodeId] = maxLevelId;
			int childl = 2 * nodeId;
			int childr = 2 * nodeId + 1;

			int nl = maxLevelId, nr = maxLevelId;
			if (childl < boxlist.size())
			{
				nl = signBvhsLevel(childl, maxLevelId + 1);
			}
			if (childr < boxlist.size())
			{
				nr = signBvhsLevel(childr, maxLevelId + 1);
			}
			return nl > nr ? nl : nr;
		}

		inline void intersect_box(const Vector3 &bbd0, const Vector3 &bbd1, std::vector<unsigned int> &list) const
		{
			std::vector<unsigned int> tmp;
			assert(n_corners >= 0);
			box_search_recursive(bbd0, bbd1, tmp, 1, 0, n_corners);

			list.resize(tmp.size());
			for (int i = 0; i < tmp.size(); ++i)
				list[i] = new2old[tmp[i]];
		}

		bool box_intersects_leaf_box(const Vector3 &bbd0, const Vector3 &bbd1, int index) const;

		void refine_box(const std::vector < std::array<Vector3, 2>>& cornerlist);

		inline std::array<Vector3, 2> refine_box_recursive(const std::vector < std::array<Vector3, 2>>& cornerlist, int n, int b, int e);

		void show(int n)
		{
			if (n >= boxlist.size() || n <= 0) n = 1;
			Box box;
			QVector3D bmin(0, 0, 0);
			QVector3D bmax(0, 0, 0);
			for (int j = 0; j < 3; j++)
			{
				bmin[j] = boxlist[n][0][j];
				bmax[j] = boxlist[n][1][j];
			}
			box.setFromBounding(bmin, bmax);
			box.show();


			if (n > 1 && n % 2 == 1 && n + 1 < boxlist.size())
			{
				for (int j = 0; j < 3; j++)
				{
					bmin[j] = boxlist[n + 1][0][j];
					bmax[j] = boxlist[n + 1][1][j];
				}
				box.setFromBounding(bmin, bmax);
				box.show();
			}
		}
	};



}


#endif
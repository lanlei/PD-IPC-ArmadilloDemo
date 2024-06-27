#include "AipcBvhs.h"
#include <iostream>

namespace AIPC
{
	bool box_box_intersection(const Vector3 &min1, const Vector3 &max1,
		const Vector3 &min2, const Vector3 &max2)
	{
		if (max1[0] < min2[0] || max1[1] < min2[1] || max1[2] < min2[2])
			return 0;
		if (max2[0] < min1[0] || max2[1] < min1[1] || max2[2] < min1[2])
			return 0;
		return 1;
	}

	bool box_box_intersection(const Vector3 &min1, const Vector3 &max1,
		const Vector3 &min2, const Vector3 &max2, Vector3 &overlap_min, Vector3 &overlap_max)
	{
		if (max1[0] < min2[0] || max1[1] < min2[1] || max1[2] < min2[2])
			return 0;
		if (max2[0] < min1[0] || max2[1] < min1[1] || max2[2] < min1[2])
			return 0;

		qeal cx0, cx1, cx2, cx3;
		qeal cy0, cy1, cy2, cy3;
		qeal cz0, cz1, cz2, cz3;
		std::vector<qeal> overlap_x = { min1[0], max1[0], min2[0], max2[0] };
		std::vector<qeal> overlap_y = { min1[1], max1[1], min2[1], max2[1] };
		std::vector<qeal> overlap_z = { min1[2], max1[2], min2[2], max2[2] };
		std::sort(overlap_x.begin(), overlap_x.end());
		std::sort(overlap_y.begin(), overlap_y.end());
		std::sort(overlap_z.begin(), overlap_z.end());

		overlap_min[0] = overlap_x[1];
		overlap_max[0] = overlap_x[2];

		overlap_min[1] = overlap_y[1];
		overlap_max[1] = overlap_y[2];

		overlap_min[2] = overlap_z[1];
		overlap_max[2] = overlap_z[2];

		return 1;
	}


	void AipcBvhs::init_boxes_recursive(const std::vector<std::array<Vector3, 2>> &cornerlist,
		int node_index, int b, int e)
	{
		assert(b != e);
		assert(node_index < boxlist.size());

		if (b + 1 == e)
		{
			boxlist[node_index] = cornerlist[b];
			oldboxId[new2old[b]] = node_index;
			leafNodeFlag[node_index] = 1;
			nodeIdToElementId[node_index] = new2old[b];
			return;
		}
		int m = b + (e - b) / 2;
		int childl = 2 * node_index;
		int childr = 2 * node_index + 1;

		assert(childl < boxlist.size());
		assert(childr < boxlist.size());

		init_boxes_recursive(cornerlist, childl, b, m);
		init_boxes_recursive(cornerlist, childr, m, e);

		assert(childl < boxlist.size());
		assert(childr < boxlist.size());
		for (int c = 0; c < 3; ++c)
		{
			boxlist[node_index][0][c] = std::min(boxlist[childl][0][c], boxlist[childr][0][c]);
			boxlist[node_index][1][c] = std::max(boxlist[childl][1][c], boxlist[childr][1][c]);
		}

		if (nodeIdToElementId[childl] != -2 || nodeIdToElementId[childr] != 2)
			nodeIdToElementId[node_index] = -1;

	}

	void AipcBvhs::box_search_recursive(const Vector3 &bbd0, const Vector3 &bbd1, std::vector<unsigned int> &list,
		int n, int b, int e) const
	{
		assert(e != b);

		assert(n < boxlist.size());
		bool cut = box_intersects_box(bbd0, bbd1, n);

		if (cut == false)
			return;

		// Leaf case
		if (e == b + 1)
		{
			list.emplace_back(b);
			return;
		}

		int m = b + (e - b) / 2;
		int childl = 2 * n;
		int childr = 2 * n + 1;

		//assert(childl < boxlist.size());
		//assert(childr < boxlist.size());

		// Traverse the "nearest" child first, so that it has more chances
		// to prune the traversal of the other child.
		box_search_recursive(
			bbd0, bbd1, list,
			childl, b, m);
		box_search_recursive(
			bbd0, bbd1, list,
			childr, m, e);
	}

	int AipcBvhs::max_node_index(int node_index, int b, int e) // 1, 0, n_corners
	{
		assert(e > b);
		if (b + 1 == e)
		{
			return node_index;
		}
		int m = b + (e - b) / 2;
		int childl = 2 * node_index;
		int childr = 2 * node_index + 1;
		return std::max(
			max_node_index(childl, b, m),
			max_node_index(childr, m, e));
	}

	void AipcBvhs::init(const std::vector<std::array<Vector3, 2>> &cornerlist)
	{
		n_corners = cornerlist.size();

		MatrixX box_centers(n_corners, 3);
		for (int i = 0; i < n_corners; ++i)
		{
			box_centers.row(i) = (cornerlist[i][0] + cornerlist[i][1]) / 2;
		}

		const VectorR3 vmin = box_centers.colwise().minCoeff();
		const VectorR3 vmax = box_centers.colwise().maxCoeff();
		const VectorR3 center = (vmin + vmax) / 2;
		for (int i = 0; i < n_corners; i++)
		{
			// make box centered at origin
			box_centers.row(i) -= center;
		}

		// after placing box at origin, vmax and vmin are symetric.
		const Vector3 scale_point = vmax - center;
		const qeal scale = scale_point.lpNorm<Eigen::Infinity>();
		// if the box is too big, resize it
		if (scale > 100)
		{
			box_centers /= scale;
		}

		struct sortstruct
		{
			int order; // face id
			Resorting::MortonCode64 morton;
		};
		std::vector<sortstruct> list;
		const int multi = 1000;
		list.resize(n_corners);

		for (int i = 0; i < n_corners; i++)
		{
			const MatrixX tmp = box_centers.row(i) * multi;

			list[i].morton = Resorting::MortonCode64(int(tmp(0)), int(tmp(1)), int(tmp(2)));
			list[i].order = i; // face id
		}

		const auto morton_compare = [](const sortstruct &a, const sortstruct &b) {
			return (a.morton < b.morton);
		};
		std::sort(list.begin(), list.end(), morton_compare);

		new2old.resize(n_corners);
		for (int i = 0; i < n_corners; i++)
		{
			new2old[i] = list[i].order;
		}
		oldboxId.resize(n_corners);
		std::vector<std::array<Vector3, 2>> sorted_cornerlist(n_corners);

		for (int i = 0; i < n_corners; i++)
		{
			sorted_cornerlist[i] = cornerlist[list[i].order];
		}

		boxlist.resize(max_node_index(1, 0, n_corners) + 1); // <-- this is because size == max_index + 1 !!!
		leafNodeFlag.resize(boxlist.size(), 0);
		nodeIdToElementId.resize(boxlist.size(), -2);
		init_boxes_recursive(sorted_cornerlist, 1, 0, n_corners);

		//

		nodeLevelId.resize(boxlist.size(), -1);
		maxLevelId = signBvhsLevel(1, 0);
		nodeLevelGroup.resize(maxLevelId + 1);

		for (int i = 0; i < nodeLevelId.size(); i++)
		{
			if (nodeLevelId[i] < 0)
				continue;
			nodeLevelGroup[maxLevelId - nodeLevelId[i]].push_back(i);
		}
	}

	void AipcBvhs::init(const MatrixX &V, const MatrixXi &F, const qeal tol)
	{
		assert(F.cols() == 3);
		assert(V.cols() == 3);

		std::vector<std::array<Vector3, 2>> cornerlist(F.rows());

		for (int i = 0; i < F.rows(); i++)
		{
			const Vector3i face = F.row(i);
			const VectorR3 v0 = V.row(face(0));
			const VectorR3 v1 = V.row(face(1));
			const VectorR3 v2 = V.row(face(2));

			Matrix3 tmp;
			tmp.row(0) = v0;
			tmp.row(1) = v1;
			tmp.row(2) = v2;

			const VectorR3 min = tmp.colwise().minCoeff().array() - tol;
			const VectorR3 max = tmp.colwise().maxCoeff().array() + tol;

			cornerlist[i][0] = min.transpose();
			cornerlist[i][1] = max.transpose();
		}

		init(cornerlist);
	}

	bool AipcBvhs::box_intersects_box(const Vector3 &bbd0, const Vector3 &bbd1, int index) const
	{
		const auto &bmin = boxlist[index][0];
		const auto &bmax = boxlist[index][1];

		return box_box_intersection(bbd0, bbd1, bmin, bmax);
	}

	bool AipcBvhs::box_intersects_leaf_box(const Vector3 & bbd0, const Vector3 & bbd1, int index) const
	{
		assert(index < oldboxId.size());
		int bbox_id = oldboxId[index];
		return box_box_intersection(bbd0, bbd1, boxlist[bbox_id][0], boxlist[bbox_id][1]);
	}

	void AipcBvhs::refine_box(const std::vector < std::array<Vector3, 2>>& cornerlist)
	{
		refine_box_recursive(cornerlist, 1, 0, n_corners);
	}

	std::array<Vector3, 2> AipcBvhs::refine_box_recursive(const std::vector < std::array<Vector3, 2>>& cornerlist, int n, int b, int e)
	{
		assert(e != b);

		assert(n < boxlist.size());
		// Leaf case
		if (e == b + 1)
		{
			// update face bbox
			const auto cmin = cornerlist[new2old[b]][0];
			const auto cmax = cornerlist[new2old[b]][1];

			for (int c = 0; c < 3; ++c)
			{
				boxlist[n][0][c] = cmin[c];
				boxlist[n][1][c] = cmax[c];
			}
		}
		else
		{
			int m = b + (e - b) / 2;
			int childl = 2 * n;
			int childr = 2 * n + 1;

			//assert(childl < boxlist.size());
			//assert(childr < boxlist.size());

			// Traverse the "nearest" child first, so that it has more chances
			// to prune the traversal of the other child.
			std::array<Vector3, 2> lb = refine_box_recursive(
				cornerlist,
				childl, b, m);
			std::array<Vector3, 2> rb = refine_box_recursive(
				cornerlist,
				childr, m, e);

			for (int c = 0; c < 3; ++c)
			{
				boxlist[n][0][c] = std::min(lb[0][c], rb[0][c]);
				boxlist[n][1][c] = std::max(lb[1][c], rb[1][c]);
			}
		}
		std::array<Vector3, 2> result;
		result[0] = boxlist[n][0];
		result[1] = boxlist[n][1];
		return result;
	}
}
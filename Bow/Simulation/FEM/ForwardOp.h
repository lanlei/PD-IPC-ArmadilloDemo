#ifndef FEM_DRIVER_H
#define FEM_DRIVER_H
#include <Bow/Macros.h>
#include <Bow/Types.h>
#include <Bow/Physics/ConstitutiveModel.h>
#include <Bow/Physics/FixedCorotated.h>
#include <oneapi/tbb.h>
#include <Eigen/Sparse>
#include <Bow/Utils/Timer.h>
#include <Bow/Optimization/Newton.h>
#include <Bow/Simulation/FEM/InertialEnergy.h>
#include <Bow/Simulation/FEM/ElasticityEnergy.h>
#include <Bow/Simulation/FEM/ExternalForceEnergy.h>
#include <Bow/Geometry/AnalyticalLevelSet.h>

#include "Utils.h"
#include "../Prototypes.h"
namespace Bow {
namespace FEM {

class AbstractOp {
};

template <class T, int dim>
class InitializeOp : public AbstractOp {
public:
    // Inputs:
    const Field<Vector<T, dim>>& m_X; // material coordinate;
    const Field<Vector<int, dim + 1>>& m_elem; // element vertex indices
    const Field<T>& m_density; // element density

    // Outputs:
    Field<T>& m_mass; // node mass
    Field<T>& m_vol; // element volume
    Field<Matrix<T, dim, dim>>& m_IB; // inverse basis

    inline void operator()();
};

template <class T, int dim>
class BoundaryConditionUpdateOp : public AbstractOp {
public:
    Geometry::BoundaryConditionManager<T, dim>& BC;
    const Field<Vector<T, dim>>& m_x;
    Field<Matrix<T, dim, dim>>& BC_basis;
    std::vector<int>& BC_order;
    Field<Vector<T, dim>>& BC_target;
    T dt;

    void operator()()
    {
        tbb::parallel_for(size_t(0), m_x.size(), [&](size_t i) {
            if (!BC.level_set_based_update(m_x[i], dt, BC_basis[i], BC_order[i], BC_target[i]))
                BC.index_based_update(i, m_x[i], dt, BC_basis[i], BC_order[i], BC_target[i]);
        });
    }
};

template <class T, int dim, class _StorageIndex = int>
class BackwardEulerUpdateOp : public Optimization::Newton<T, dim, _StorageIndex> {
public:
    using StorageIndex = _StorageIndex;
    using Base = Optimization::Newton<T, dim, StorageIndex>;
    using Vec = Bow::Vector<T, Eigen::Dynamic>;
    using Mat = Eigen::SparseMatrix<T, Eigen::ColMajor, StorageIndex>;

    bool project_dirichlet = true;

    const Field<Matrix<T, dim, dim>>& BC_basis;
    const std::vector<int>& BC_order;
    const Field<Vector<T, dim>>& BC_target;

    Field<Vector<T, dim>>& m_x;
    Field<Vector<T, dim>>& m_v;
    Field<Vector<T, dim>>& m_a;
    Field<Vector<T, dim>>& m_x_tilde;
    Field<T>& m_mass;

    Mat m_transform_matrix;

    enum TSMethod {
        BE,
        NM
    };
    TSMethod tsMethod = BE;
    T tsParam[2][3] = {
        { 1, 0.5, 1 },
        { 0.5, 0.25, 0.5 }
    };

    T dt = 0.02;

    BackwardEulerUpdateOp(const Field<Matrix<T, dim, dim>>& BC_basis, const std::vector<int>& BC_order, const Field<Vector<T, dim>>& BC_target, Field<T>& mass, Field<Vector<T, dim>>& x, Field<Vector<T, dim>>& v, Field<Vector<T, dim>>& a, Field<Vector<T, dim>>& x_tilde);
    virtual void gradient(const Vec& x_vec, Vec& grad_vec) override;
    virtual void hessian(const Vec& x_vec, Mat& hess, const bool project_pd) override;
    virtual T residual(const Vec& x, const Vec& grad, const Vec& direction) override;
    void postprocess(Vec& direction) override;
    void set_ts_weights();
    void initialize_acceleration();
    void update_predictive_pos();
    void update_transformation_matrix();
    void advance(const Vec& new_x);
    void operator()();
};
}
} // namespace Bow::FEM

#include "ForwardOp.tpp"

#endif
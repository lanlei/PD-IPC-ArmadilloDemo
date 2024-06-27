#ifndef SIMULATION_OP_PROTOTYPES_H
#define SIMULATION_OP_PROTOTYPES_H
#include <Bow/Types.h>
#include <Bow/Utils/Logging.h>
#include <Eigen/SparseCore>

namespace Bow {

template <class T, int dim, class StorageIndex = int>
class EnergyOp {
public:
    T energy_scale = 1;
    /* called whenever x is changed */
    virtual void precompute(const Field<Vector<T, dim>>& x){};
    /* called before each timestep */
    virtual void callback(const Field<Vector<T, dim>>& xn){};
    virtual T stepsize_upperbound(const Field<Vector<T, dim>>& x, const Field<Vector<T, dim>>& dx) { return T(1); }
    virtual T energy(const Field<Vector<T, dim>>& x) { return 0; }
    virtual void gradient(const Field<Vector<T, dim>>& x, Field<Vector<T, dim>>& grad){};
    virtual void multiply(const Field<Vector<T, dim>>& x, Field<Vector<T, dim>>& Ax, bool project_pd = true) { BOW_NOT_IMPLEMENTED }
    virtual void project(Field<Vector<T, dim>>& b) { BOW_NOT_IMPLEMENTED }
    virtual void precondition(Field<Vector<T, dim>>& diagonal) { BOW_NOT_IMPLEMENTED }
    virtual void hessian(const Field<Vector<T, dim>>& x, Eigen::SparseMatrix<T, Eigen::ColMajor, StorageIndex>& hess, bool project_pd = true) { BOW_NOT_IMPLEMENTED }
    virtual void postprocess(Field<Vector<T, dim>>& direction){};
    virtual void internal_force(const Field<Vector<T, dim>>& xn, const Field<Vector<T, dim>>& vn, Field<Vector<T, dim>>& force)
    {
        force.resize(xn.size(), Vector<T, dim>::Zero());
        std::fill(force.begin(), force.end(), Vector<T, dim>::Zero());
    }
};
} // namespace Bow

#endif
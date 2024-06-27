#pragma once

#include <Bow/Math/SVD.h>
#include <Bow/Simulation/MPM/ElasticityOp.h>

namespace Bow::MPM {

template <class T, int dim>
class PlasticityOp {
public:
    virtual void project_strain() = 0;
};

template <class T, int dim>
class VonMisesStvkHencky : public PlasticityOp<T, dim> {
    using TM = Matrix<T, dim, dim>;
    using TV = Vector<T, dim>;

    std::shared_ptr<StvkWithHenckyOp<T, dim>> stvk;
    T yield_stress, fail_stress, xi;

public:
    VonMisesStvkHencky(std::shared_ptr<ElasticityOp<T, dim>> stvk, T yield_stress, T fail_stress, T xi)
        : stvk(std::dynamic_pointer_cast<StvkWithHenckyOp<T, dim>>(stvk)), yield_stress(yield_stress), fail_stress(fail_stress), xi(xi) {}

    void project_strain() override
    {
        tbb::parallel_for(size_t(0), stvk->m_F.size(), [&](size_t i) {
            TM& F = stvk->m_F[i];
            TM U, V;
            TV sigma;

            // TODO: this is inefficient because next time step updateState will do the svd again!
            Math::svd(F, U, sigma, V);

            //TV epsilon = sigma.array().log();
            TV epsilon = sigma.array().max(1e-4).log(); //TODO: need the max part?
            T trace_epsilon = epsilon.sum();
            TV epsilon_hat = epsilon - (trace_epsilon / (T)dim) * TV::Ones();
            T epsilon_hat_squared_norm = epsilon_hat.squaredNorm();
            T epsilon_hat_norm = std::sqrt(epsilon_hat_squared_norm);
            T delta_gamma = epsilon_hat_norm - yield_stress / (2 * stvk->mu);
            if (delta_gamma <= 0) // case I
            {
                return;
            }
            //hardening
            yield_stress -= xi * delta_gamma; //supposed to only increase yield_stress
            //yield_stress = std::max((T)0, yield_stress);

            TV H = epsilon - (delta_gamma / epsilon_hat_norm) * epsilon_hat; // case II
            TV exp_H = H.array().exp();
            F = U * exp_H.asDiagonal() * V.transpose();
        });
    }
};

} // namespace Bow::MPM
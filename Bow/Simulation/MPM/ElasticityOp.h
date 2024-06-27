#ifndef ELASTICITY_OP_H
#define ELASTICITY_OP_H

#include <Bow/IO/partio.h>
#include <Bow/Types.h>
#include <Bow/Physics/FixedCorotated.h>
#include <Bow/Physics/NeoHookean.h>
#include <Bow/Physics/StvkWithHenckyIsotropic.h>
#include <Bow/Physics/EquationOfState.h>
#include <Bow/Simulation/MPM/MPMTransfer.h>
#include <Bow/Math/Utils.h>
#include <Bow/Utils/Serialization.h>

namespace Bow::MPM {

template <class T, int dim>
class ElasticityOp {
public:
    virtual void append(int start, int end, T vol) = 0;
    virtual void compute_stress(Field<Matrix<T, dim, dim>>& stress) = 0;
    virtual void compute_cauchy(Field<Matrix<T, dim, dim>>& stress) { BOW_NOT_IMPLEMENTED }
    virtual void compute_criticalStress(T percent, Field<Matrix<T, dim, dim>>& stretchedCauchy) { BOW_NOT_IMPLEMENTED }
    virtual void evolve_strain(const Field<Matrix<T, dim, dim>>& m_gradXp) = 0;
    virtual void trial_strain(const Field<Matrix<T, dim, dim>>& m_gradXp) { BOW_NOT_IMPLEMENTED }
    virtual void trial_energy(Field<T>& t_energy) { BOW_NOT_IMPLEMENTED }
    virtual void trial_gradient(Field<Matrix<T, dim, dim>>& t_gradient) { BOW_NOT_IMPLEMENTED }
    virtual void trial_differential(const Field<Matrix<T, dim, dim>>& d_F, Field<Matrix<T, dim, dim>>& t_differential, bool project_pd) { BOW_NOT_IMPLEMENTED }
    virtual void trial_hessian(Field<Matrix<T, dim * dim, dim * dim>>& t_hessian, bool project_pd) { BOW_NOT_IMPLEMENTED }
    virtual T stepsize_upperbound(const Field<Matrix<T, dim, dim>>& m_gradDVp, T dt) { return 1.0; }
    virtual void collect_strain(Field<Matrix<T, dim, dim>>& m_Fs) { BOW_NOT_IMPLEMENTED }
};

template <class T, int dim>
class FixedCorotatedOp : public ElasticityOp<T, dim> {
public:
    using TV = Vector<T, dim>;
    using TM = Matrix<T, dim, dim>;
    Field<Matrix<T, dim, dim>> m_F;
    Field<Matrix<T, dim, dim>> t_F; // only used in implicit
    std::vector<T> m_vol;
    T mu, lambda;
    int index_offset;

    SERIALIZATION_REGISTER(m_F)
    SERIALIZATION_REGISTER(m_vol)

    FixedCorotatedOp(T E, T nu)
    {
        std::tie(mu, lambda) = Bow::ConstitutiveModel::lame_paramters(E, nu);
    }

    void append(int start, int end, T vol) override
    {
        index_offset = start;
        for (int i = start; i < end; ++i) {
            m_F.push_back(Matrix<T, dim, dim>::Identity());
            m_vol.push_back(vol);
        }
    }

    void evolve_strain(const Field<Matrix<T, dim, dim>>& m_gradXp) override
    {
        tbb::parallel_for(size_t(0), m_F.size(), [&](size_t i) {
            m_F[i] = (m_gradXp[index_offset + i]) * m_F[i];
        });
    }

    void trial_strain(const Field<Matrix<T, dim, dim>>& m_gradXp) override
    {
        t_F = m_F;
        tbb::parallel_for(size_t(0), m_F.size(), [&](size_t i) {
            t_F[i] = (m_gradXp[index_offset + i]) * m_F[i];
        });
    }

    void trial_energy(Field<T>& t_energy) override
    {
        tbb::parallel_for(size_t(0), m_F.size(), [&](size_t i) {
            T energy = m_vol[i] * Bow::ConstitutiveModel::FixedCorotated::psi(t_F[i], mu, lambda);
            t_energy[index_offset + i] += energy;
        });
    }

    void trial_gradient(Field<TM>& t_gradient) override
    {
        tbb::parallel_for(size_t(0), m_F.size(), [&](size_t i) {
            Matrix<T, dim, dim> first_piola;
            Bow::ConstitutiveModel::FixedCorotated::first_piola(t_F[i], mu, lambda, first_piola);
            // Eqn 194. https://www.seas.upenn.edu/~cffjiang/research/mpmcourse/mpmcourse.pdf
            Matrix<T, dim, dim> stress = m_vol[i] * first_piola * m_F[i].transpose();
            t_gradient[index_offset + i] += stress;
        });
    }

    void trial_differential(const Field<Matrix<T, dim, dim>>& d_F, Field<Matrix<T, dim, dim>>& t_differential, bool project_pd) override
    {
        tbb::parallel_for(size_t(0), m_F.size(), [&](size_t i) {
            TM A = TM::Zero();
            TM D = d_F[index_offset + i] * m_F[i];
            if (project_pd)
                Bow::ConstitutiveModel::FixedCorotated::first_piola_differential<true>(t_F[i], D, mu, lambda, A);
            else
                Bow::ConstitutiveModel::FixedCorotated::first_piola_differential<false>(t_F[i], D, mu, lambda, A);
            t_differential[index_offset + i] = m_vol[i] * A * m_F[i].transpose();
        });
    }

    void trial_hessian(Field<Matrix<T, dim * dim, dim * dim>>& t_hessian, bool project_pd) override
    {
        tbb::parallel_for(size_t(0), m_F.size(), [&](size_t i) {
            Matrix<T, dim * dim, dim * dim> dPdF;
            if (project_pd)
                ConstitutiveModel::FixedCorotated::first_piola_derivative<true>(t_F[i], mu, lambda, dPdF);
            else
                ConstitutiveModel::FixedCorotated::first_piola_derivative<false>(t_F[i], mu, lambda, dPdF);
            TM FT = m_F[i].transpose();
            Matrix<T, dim * dim, dim* dim> deformed_dPdF = Matrix<T, dim * dim, dim * dim>::Zero();
            for (int u = 0; u < dim; ++u)
                for (int v = 0; v < dim; ++v)
                    for (int x = 0; x < dim; ++x)
                        for (int p = 0; p < dim; ++p)
                            for (int q = 0; q < dim; ++q)
                                for (int y = 0; y < dim; ++y)
                                    deformed_dPdF(u + x * dim, p + y * dim) += dPdF(u + v * dim, p + q * dim) * FT(v, x) * FT(q, y);
            t_hessian[index_offset + i] += m_vol[i] * deformed_dPdF;
        });
    }

    void compute_cauchy(Field<Matrix<T, dim, dim>>& cauchy)
    {
        BOW_TIMER_FLAG("compute cauchy");
        tbb::parallel_for(size_t(0), m_F.size(), [&](size_t i) {
            Matrix<T, dim, dim> F = m_F[i];
            Matrix<T, dim, dim> first_piola;
            T J = F.determinant();
            Bow::ConstitutiveModel::FixedCorotated::first_piola(F, mu, lambda, first_piola);
            cauchy[index_offset + i] = (1.0 / J) * first_piola * F.transpose();
        });
    }

    void compute_stress(Field<Matrix<T, dim, dim>>& stress) override
    {
        BOW_TIMER_FLAG("compute elasticity");
        tbb::parallel_for(size_t(0), m_F.size(), [&](size_t i) {
            Matrix<T, dim, dim> F = m_F[i];
            Matrix<T, dim, dim> first_piola;
            Bow::ConstitutiveModel::FixedCorotated::first_piola(F, mu, lambda, first_piola);
            stress[index_offset + i] = m_vol[i] * first_piola * F.transpose();
        });
    }

    //this will have to be written for each elasticity model we want to use with MPM damage
    void compute_criticalStress(T percent, Field<Matrix<T, dim, dim>>& stretchedCauchy)
    {
        BOW_TIMER_FLAG("compute sigmaC");
        tbb::parallel_for(size_t(0), stretchedCauchy.size(), [&](size_t i) {
            TM F = TM::Identity() * (1.0 + percent); //stretched F
            T J = F.determinant();
            Matrix<T, dim, dim> first_piola;
            Bow::ConstitutiveModel::FixedCorotated::first_piola(F, mu, lambda, first_piola);
            stretchedCauchy[i] = (1.0 / J) * first_piola * F.transpose();
        });
    }

    void collect_strain(Field<Matrix<T, dim, dim>>& m_Fs) override
    {
        tbb::parallel_for(size_t(0), m_F.size(), [&](size_t i) {
            m_Fs[index_offset + i] = m_F[i];
        });
    }
};

template <class T, int dim>
class NeoHookeanOp : public ElasticityOp<T, dim> {
public:
    using TV = Vector<T, dim>;
    using TM = Matrix<T, dim, dim>;
    Field<Matrix<T, dim, dim>> m_F;
    Field<Matrix<T, dim, dim>> t_F; // only used in implicit
    std::vector<T> m_vol;
    T mu, lambda;
    int index_offset;

    SERIALIZATION_REGISTER(m_F)
    SERIALIZATION_REGISTER(m_vol)

    NeoHookeanOp(T E, T nu)
    {
        std::tie(mu, lambda) = Bow::ConstitutiveModel::lame_paramters(E, nu);
    }

    void append(int start, int end, T vol) override
    {
        index_offset = start;
        for (int i = start; i < end; ++i) {
            m_F.push_back(Matrix<T, dim, dim>::Identity());
            m_vol.push_back(vol);
        }
    }

    void evolve_strain(const Field<Matrix<T, dim, dim>>& m_gradXp) override
    {
        tbb::parallel_for(size_t(0), m_F.size(), [&](size_t i) {
            m_F[i] = (m_gradXp[index_offset + i]) * m_F[i];
        });
    }

    void trial_strain(const Field<Matrix<T, dim, dim>>& m_gradXp) override
    {
        t_F = m_F;
        tbb::parallel_for(size_t(0), m_F.size(), [&](size_t i) {
            t_F[i] = (m_gradXp[index_offset + i]) * m_F[i];
        });
    }

    void trial_energy(Field<T>& t_energy) override
    {
        tbb::parallel_for(size_t(0), m_F.size(), [&](size_t i) {
            T energy = m_vol[i] * Bow::ConstitutiveModel::NeoHookean::psi(t_F[i], mu, lambda);
            t_energy[index_offset + i] += energy;
        });
    }

    void trial_gradient(Field<TM>& t_gradient) override
    {
        tbb::parallel_for(size_t(0), m_F.size(), [&](size_t i) {
            Matrix<T, dim, dim> first_piola;
            Bow::ConstitutiveModel::NeoHookean::first_piola(t_F[i], mu, lambda, first_piola);
            // Eqn 194. https://www.seas.upenn.edu/~cffjiang/research/mpmcourse/mpmcourse.pdf
            Matrix<T, dim, dim> stress = m_vol[i] * first_piola * m_F[i].transpose();
            t_gradient[index_offset + i] += stress;
        });
    }

    void trial_differential(const Field<Matrix<T, dim, dim>>& d_F, Field<Matrix<T, dim, dim>>& t_differential, bool project_pd) override
    {
        tbb::parallel_for(size_t(0), m_F.size(), [&](size_t i) {
            TM A = TM::Zero();
            TM D = d_F[index_offset + i] * m_F[i];
            if (project_pd)
                Bow::ConstitutiveModel::NeoHookean::first_piola_differential<true>(t_F[i], D, mu, lambda, A);
            else
                Bow::ConstitutiveModel::NeoHookean::first_piola_differential<false>(t_F[i], D, mu, lambda, A);
            t_differential[index_offset + i] = m_vol[i] * A * m_F[i].transpose();
        });
    }

    void trial_hessian(Field<Matrix<T, dim * dim, dim * dim>>& t_hessian, bool project_pd) override
    {
        tbb::parallel_for(size_t(0), m_F.size(), [&](size_t i) {
            Matrix<T, dim * dim, dim * dim> dPdF;
            if (project_pd)
                ConstitutiveModel::NeoHookean::first_piola_derivative<true>(t_F[i], mu, lambda, dPdF);
            else
                ConstitutiveModel::NeoHookean::first_piola_derivative<false>(t_F[i], mu, lambda, dPdF);
            TM FT = m_F[i].transpose();
            Matrix<T, dim * dim, dim* dim> deformed_dPdF = Matrix<T, dim * dim, dim * dim>::Zero();
            for (int u = 0; u < dim; ++u)
                for (int v = 0; v < dim; ++v)
                    for (int x = 0; x < dim; ++x)
                        for (int p = 0; p < dim; ++p)
                            for (int q = 0; q < dim; ++q)
                                for (int y = 0; y < dim; ++y)
                                    deformed_dPdF(u + x * dim, p + y * dim) += dPdF(u + v * dim, p + q * dim) * FT(v, x) * FT(q, y);
            t_hessian[index_offset + i] += m_vol[i] * deformed_dPdF;
        });
    }

    void compute_cauchy(Field<Matrix<T, dim, dim>>& cauchy)
    {
        BOW_TIMER_FLAG("compute cauchy");
        tbb::parallel_for(size_t(0), m_F.size(), [&](size_t i) {
            Matrix<T, dim, dim> F = m_F[i];
            Matrix<T, dim, dim> first_piola;
            T J = F.determinant();
            Bow::ConstitutiveModel::NeoHookean::first_piola(F, mu, lambda, first_piola);
            cauchy[index_offset + i] = (1.0 / J) * first_piola * F.transpose();
        });
    }

    void compute_stress(Field<Matrix<T, dim, dim>>& stress) override
    {
        BOW_TIMER_FLAG("compute elasticity");
        tbb::parallel_for(size_t(0), m_F.size(), [&](size_t i) {
            Matrix<T, dim, dim> F = m_F[i];
            Matrix<T, dim, dim> first_piola;
            Bow::ConstitutiveModel::NeoHookean::first_piola(F, mu, lambda, first_piola);
            stress[index_offset + i] = m_vol[i] * first_piola * F.transpose();
        });
    }

    //this will have to be written for each elasticity model we want to use with MPM damage
    void compute_criticalStress(T percent, Field<Matrix<T, dim, dim>>& stretchedCauchy)
    {
        BOW_TIMER_FLAG("compute sigmaC");
        tbb::parallel_for(size_t(0), stretchedCauchy.size(), [&](size_t i) {
            TM F = TM::Identity() * (1.0 + percent); //stretched F
            T J = F.determinant();
            Matrix<T, dim, dim> first_piola;
            Bow::ConstitutiveModel::NeoHookean::first_piola(F, mu, lambda, first_piola);
            stretchedCauchy[i] = (1.0 / J) * first_piola * F.transpose();
        });
    }

    T stepsize_upperbound(const Field<Matrix<T, dim, dim>>& m_gradDVp, T dt)
    {
        Vector<T, Eigen::Dynamic> alphas(m_F.size());
        alphas.setOnes();
        tbb::parallel_for(size_t(0), m_F.size(), [&](size_t i) {
            Matrix<T, dim, dim> A = t_F[i].transpose().partialPivLu().solve((dt * m_gradDVp[index_offset + i] * m_F[i]).transpose());
            T a, b, c, d;
            if constexpr (dim == 2) {
                a = 0;
                b = A.determinant();
            }
            else {
                a = A.determinant();
                b = A(0, 0) * A(1, 1) - A(0, 1) * A(1, 0) + A(0, 0) * A(2, 2) - A(0, 2) * A(2, 0) + A(1, 1) * A(2, 2) - A(1, 2) * A(2, 1);
            }
            c = A.diagonal().sum();
            d = 0.9;

            T t = Math::get_smallest_positive_real_cubic_root(a, b, c, d);
            if (t < 0 || t > 1) t = 1;
            alphas(i) = t;
        });
        return alphas.minCoeff();
    }
};

template <class T, int dim>
class StvkWithHenckyOp : public ElasticityOp<T, dim> {
public:
    using TV = Vector<T, dim>;
    using TM = Matrix<T, dim, dim>;
    Field<Matrix<T, dim, dim>> m_F;
    Field<Matrix<T, dim, dim>> t_F; // only used in implicit
    std::vector<T> m_vol;
    T mu, lambda;
    int index_offset;

    SERIALIZATION_REGISTER(m_F)
    SERIALIZATION_REGISTER(m_vol)

    StvkWithHenckyOp(T E, T nu)
    {
        std::tie(mu, lambda) = Bow::ConstitutiveModel::lame_paramters(E, nu);
    }

    void append(int start, int end, T vol) override
    {
        index_offset = start;
        for (int i = start; i < end; ++i) {
            m_F.push_back(Matrix<T, dim, dim>::Identity());
            m_vol.push_back(vol);
        }
    }

    void evolve_strain(const Field<Matrix<T, dim, dim>>& m_gradXp) override
    {
        tbb::parallel_for(size_t(0), m_F.size(), [&](size_t i) {
            m_F[i] = (m_gradXp[index_offset + i]) * m_F[i];
        });
    }

    void trial_strain(const Field<Matrix<T, dim, dim>>& m_gradXp) override
    {
        t_F = m_F;
        tbb::parallel_for(size_t(0), m_F.size(), [&](size_t i) {
            t_F[i] = (m_gradXp[index_offset + i]) * m_F[i];
        });
    }

    void trial_energy(Field<T>& t_energy) override
    {
        tbb::parallel_for(size_t(0), m_F.size(), [&](size_t i) {
            T energy = m_vol[i] * Bow::ConstitutiveModel::StvkWithHenckyIsotropic::psi(t_F[i], mu, lambda);
            t_energy[index_offset + i] += energy;
        });
    }

    void trial_gradient(Field<TM>& t_gradient) override
    {
        tbb::parallel_for(size_t(0), m_F.size(), [&](size_t i) {
            Matrix<T, dim, dim> first_piola;
            Bow::ConstitutiveModel::StvkWithHenckyIsotropic::first_piola(t_F[i], mu, lambda, first_piola);
            // Eqn 194. https://www.seas.upenn.edu/~cffjiang/research/mpmcourse/mpmcourse.pdf
            Matrix<T, dim, dim> stress = m_vol[i] * first_piola * m_F[i].transpose();
            t_gradient[index_offset + i] += stress;
        });
    }

    void trial_differential(const Field<Matrix<T, dim, dim>>& d_F, Field<Matrix<T, dim, dim>>& t_differential, bool project_pd) override
    {
        tbb::parallel_for(size_t(0), m_F.size(), [&](size_t i) {
            TM A = TM::Zero();
            TM D = d_F[index_offset + i] * m_F[i];
            if (project_pd)
                Bow::ConstitutiveModel::StvkWithHenckyIsotropic::first_piola_differential<true>(t_F[i], D, mu, lambda, A);
            else
                Bow::ConstitutiveModel::StvkWithHenckyIsotropic::first_piola_differential<false>(t_F[i], D, mu, lambda, A);
            t_differential[index_offset + i] = m_vol[i] * A * m_F[i].transpose();
        });
    }

    void trial_hessian(Field<Matrix<T, dim * dim, dim * dim>>& t_hessian, bool project_pd) override
    {
        tbb::parallel_for(size_t(0), m_F.size(), [&](size_t i) {
            Matrix<T, dim * dim, dim * dim> dPdF;
            if (project_pd)
                ConstitutiveModel::StvkWithHenckyIsotropic::first_piola_derivative<true>(t_F[i], mu, lambda, dPdF);
            else
                ConstitutiveModel::StvkWithHenckyIsotropic::first_piola_derivative<false>(t_F[i], mu, lambda, dPdF);
            TM FT = m_F[i].transpose();
            Matrix<T, dim * dim, dim* dim> deformed_dPdF = Matrix<T, dim * dim, dim * dim>::Zero();
            for (int u = 0; u < dim; ++u)
                for (int v = 0; v < dim; ++v)
                    for (int x = 0; x < dim; ++x)
                        for (int p = 0; p < dim; ++p)
                            for (int q = 0; q < dim; ++q)
                                for (int y = 0; y < dim; ++y)
                                    deformed_dPdF(u + x * dim, p + y * dim) += dPdF(u + v * dim, p + q * dim) * FT(v, x) * FT(q, y);
            t_hessian[index_offset + i] += m_vol[i] * deformed_dPdF;
        });
    }

    void compute_cauchy(Field<Matrix<T, dim, dim>>& cauchy)
    {
        BOW_TIMER_FLAG("compute cauchy");
        tbb::parallel_for(size_t(0), m_F.size(), [&](size_t i) {
            Matrix<T, dim, dim> F = m_F[i];
            Matrix<T, dim, dim> first_piola;
            T J = F.determinant();
            Bow::ConstitutiveModel::StvkWithHenckyIsotropic::first_piola(F, mu, lambda, first_piola);
            cauchy[index_offset + i] = (1.0 / J) * first_piola * F.transpose();
        });
    }

    void compute_stress(Field<Matrix<T, dim, dim>>& stress) override
    {
        BOW_TIMER_FLAG("compute elasticity");
        tbb::parallel_for(size_t(0), m_F.size(), [&](size_t i) {
            Matrix<T, dim, dim> F = m_F[i];
            Matrix<T, dim, dim> first_piola;
            Bow::ConstitutiveModel::StvkWithHenckyIsotropic::first_piola(F, mu, lambda, first_piola);
            stress[index_offset + i] = m_vol[i] * first_piola * F.transpose();
        });
    }

    //this will have to be written for each elasticity model we want to use with MPM damage
    void compute_criticalStress(T percent, Field<Matrix<T, dim, dim>>& stretchedCauchy)
    {
        BOW_TIMER_FLAG("compute sigmaC");
        tbb::parallel_for(size_t(0), stretchedCauchy.size(), [&](size_t i) {
            TM F = TM::Identity() * (1.0 + percent); //stretched F
            T J = F.determinant();
            Matrix<T, dim, dim> first_piola;
            Bow::ConstitutiveModel::StvkWithHenckyIsotropic::first_piola(F, mu, lambda, first_piola);
            stretchedCauchy[i] = (1.0 / J) * first_piola * F.transpose();
        });
    }

    T stepsize_upperbound(const Field<Matrix<T, dim, dim>>& m_gradDVp, T dt)
    {
        Vector<T, Eigen::Dynamic> alphas(m_F.size());
        alphas.setOnes();
        tbb::parallel_for(size_t(0), m_F.size(), [&](size_t i) {
            Matrix<T, dim, dim> A = t_F[i].transpose().partialPivLu().solve((dt * m_gradDVp[index_offset + i] * m_F[i]).transpose());
            T a, b, c, d;
            if constexpr (dim == 2) {
                a = 0;
                b = A.determinant();
            }
            else {
                a = A.determinant();
                b = A(0, 0) * A(1, 1) - A(0, 1) * A(1, 0) + A(0, 0) * A(2, 2) - A(0, 2) * A(2, 0) + A(1, 1) * A(2, 2) - A(1, 2) * A(2, 1);
            }
            c = A.diagonal().sum();
            d = 0.9;

            T t = Math::get_smallest_positive_real_cubic_root(a, b, c, d);
            if (t < 0 || t > 1) t = 1;
            alphas(i) = t;
        });
        return alphas.minCoeff();
    }
};

template <class T, int dim>
class EquationOfStateOp : public ElasticityOp<T, dim> {
public:
    std::vector<T> m_J;
    std::vector<T> t_J; // only used in implicit
    std::vector<T> m_vol;
    T bulk, gamma;
    int index_offset;

    SERIALIZATION_REGISTER(m_J)
    SERIALIZATION_REGISTER(m_vol)

    EquationOfStateOp(T bulk, T gamma)
        : bulk(bulk), gamma(gamma) {}

    void append(int start, int end, T vol) override
    {
        index_offset = start;
        for (int i = start; i < end; ++i) {
            m_J.push_back(1);
            m_vol.push_back(vol);
        }
    }

    void evolve_strain(const Field<Matrix<T, dim, dim>>& m_gradXp) override
    {
        tbb::parallel_for(size_t(0), m_J.size(), [&](size_t i) {
            m_J[i] = (1 + (m_gradXp[index_offset + i].trace() - dim)) * m_J[i];
        });
    }

    void trial_strain(const Field<Matrix<T, dim, dim>>& m_gradXp) override
    {
        t_J = m_J;
        tbb::parallel_for(size_t(0), m_J.size(), [&](size_t i) {
            t_J[i] = (1 + (m_gradXp[index_offset + i].trace() - dim)) * m_J[i];
        });
    }

    void trial_energy(Field<T>& t_energy) override
    {
        tbb::parallel_for(size_t(0), m_J.size(), [&](size_t i) {
            T energy = m_vol[i] * ConstitutiveModel::EquationOfState::psi(t_J[i], bulk, gamma);
            t_energy[index_offset + i] += energy;
        });
    }

    void trial_gradient(Field<Matrix<T, dim, dim>>& t_gradient) override
    {
        tbb::parallel_for(size_t(0), m_J.size(), [&](size_t i) {
            T first_piola;
            ConstitutiveModel::EquationOfState::first_piola(t_J[i], bulk, gamma, first_piola);
            // Eqn 194. https://www.seas.upenn.edu/~cffjiang/research/mpmcourse/mpmcourse.pdf
            Matrix<T, dim, dim> stress = m_vol[i] * first_piola * m_J[i] * Matrix<T, dim, dim>::Identity();
            t_gradient[index_offset + i] += stress;
        });
    }

    void trial_hessian(Field<Matrix<T, dim * dim, dim * dim>>& t_hessian, bool project_pd) override
    {
        tbb::parallel_for(size_t(0), m_J.size(), [&](size_t i) {
            T P, dPdJ;
            ConstitutiveModel::EquationOfState::first_piola(t_J[i], bulk, gamma, P);
            ConstitutiveModel::EquationOfState::first_piola_derivative(t_J[i], bulk, gamma, dPdJ);
            Matrix<T, dim * dim, dim* dim> deformed_dPdF = Matrix<T, dim * dim, dim * dim>::Zero();
            for (int u = 0; u < dim; ++u)
                for (int x = 0; x < dim; ++x)
                    for (int p = 0; p < dim; ++p)
                        for (int y = 0; y < dim; ++y)
                            if (u == x && p == y)
                                deformed_dPdF(u + x * dim, p + y * dim) += dPdJ * m_J[i] * m_J[i];
            t_hessian[index_offset + i] += m_vol[i] * deformed_dPdF;
        });
    }

    void compute_stress(Field<Matrix<T, dim, dim>>& stress) override
    {
        BOW_TIMER_FLAG("compute elasticity");
        tbb::parallel_for(size_t(0), m_J.size(), [&](size_t i) {
            T J = m_J[i];
            T first_piola;
            Bow::ConstitutiveModel::EquationOfState::first_piola(J, bulk, gamma, first_piola);
            stress[index_offset + i] = m_vol[i] * first_piola * J * Matrix<T, dim, dim>::Identity();
        });
    }

    T stepsize_upperbound(const Field<Matrix<T, dim, dim>>& m_gradDVp, T dt)
    {
        Vector<T, Eigen::Dynamic> alphas(m_J.size());
        alphas.setOnes();
        tbb::parallel_for(size_t(0), m_J.size(), [&](size_t i) {
            T A = (dt * m_gradDVp[index_offset + i].trace() * m_J[i] / t_J[i]);
            T d = 0.9;
            T alpha = d / A;
            if (alpha <= 0 || alpha > 1) alpha = 1;
            alphas(i) = alpha;
        });
        return alphas.minCoeff();
    }
};

} // namespace Bow::MPM

#endif

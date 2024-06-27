#pragma once

#include <Bow/Types.h>
#include <Partio.h>

namespace Bow {

template <class T, int dim>
void writePositionVectorToPartio(const std::string& particleFile,
    const Field<Vector<T, dim>>& particles)
{
    Partio::ParticlesDataMutable* parts = Partio::create();
    Partio::ParticleAttribute posH;
    posH = parts->addAttribute("position", Partio::VECTOR, 3);
    for (size_t iter = 0; iter < particles.size(); ++iter) {
        int idx = parts->addParticle();
        float* p = parts->dataWrite<float>(posH, idx);
        for (int k = 0; k < 3; k++)
            p[k] = (T)0;
        for (int k = 0; k < dim; k++)
            p[k] = particles[iter](k);
    }
    Partio::write(particleFile.c_str(), *parts);
    parts->release();
}

template <class T, int dim, typename OP>
void writePositionVectorToPartio(const std::string& particleFile,
    const Field<Vector<T, dim>>& particles, const OP& transform)
{
    Partio::ParticlesDataMutable* parts = Partio::create();
    Partio::ParticleAttribute posH;
    posH = parts->addAttribute("position", Partio::VECTOR, 3);
    for (size_t iter = 0; iter < particles.size(); ++iter) {
        int idx = parts->addParticle();
        float* p = parts->dataWrite<float>(posH, idx);
        for (int k = 0; k < 3; k++)
            p[k] = (T)0;
        auto pos_after_trans = transform(particles[iter]);
        for (int k = 0; k < dim; k++)
            p[k] = pos_after_trans(k);
    }
    Partio::write(particleFile.c_str(), *parts);
    parts->release();
}

template <class T, int dim>
void writePositionVectorToPartio(const std::string& particleFile,
    const Field<Vector<T, dim>>& particles,
    const std::vector<T>& info_list)
{
    Partio::ParticlesDataMutable* parts = Partio::create();
    Partio::ParticleAttribute posH, infoH;
    posH = parts->addAttribute("position", Partio::VECTOR, 3);
    infoH = parts->addAttribute("info", Partio::VECTOR, 1);
    for (size_t iter = 0; iter < particles.size(); ++iter) {
        int idx = parts->addParticle();
        float* p = parts->dataWrite<float>(posH, idx);
        float* info = parts->dataWrite<float>(infoH, idx);
        for (int k = 0; k < 3; k++)
            p[k] = (T)0;
        for (int k = 0; k < dim; k++)
            p[k] = particles[iter](k);
        info[0] = info_list[iter];
    }
    Partio::write(particleFile.c_str(), *parts);
    parts->release();
}

template <class T, int dim>
void writePositionVectorToPartio(const std::string& particleFile,
    const Field<Vector<T, dim>>& particles,
    const Field<Vector<T, dim>>& info_list)
{
    Partio::ParticlesDataMutable* parts = Partio::create();
    Partio::ParticleAttribute posH, infoH;
    posH = parts->addAttribute("position", Partio::VECTOR, 3);
    infoH = parts->addAttribute("info", Partio::VECTOR, 3);
    for (size_t iter = 0; iter < particles.size(); ++iter) {
        int idx = parts->addParticle();
        float* p = parts->dataWrite<float>(posH, idx);
        float* info = parts->dataWrite<float>(infoH, idx);
        for (int k = 0; k < 3; k++)
            p[k] = (T)0;
        for (int k = 0; k < dim; k++)
            p[k] = particles[iter](k);
        for (int k = 0; k < 3; k++)
            info[k] = (T)0;
        for (int k = 0; k < dim; ++k)
            info[k] = info_list[iter](k);
    }
    Partio::write(particleFile.c_str(), *parts);
    parts->release();
}

template <class T, int dim>
void writeDFGMPMToPartio(const std::string& particleFile,
    const Field<Vector<T, dim>>& particles,
    const Field<Vector<T, dim>>& velocities,
    const Field<Vector<T, dim>>& damageGradients,
    const std::vector<T>& masses,
    const std::vector<T>& damage,
    const std::vector<int>& sp)
{
    Partio::ParticlesDataMutable* parts = Partio::create();
    Partio::ParticleAttribute posH, velH, dgH, massH, damageH, spH;
    posH = parts->addAttribute("position", Partio::VECTOR, 3);
    velH = parts->addAttribute("v", Partio::VECTOR, 3);
    dgH = parts->addAttribute("DG", Partio::VECTOR, 3);
    massH = parts->addAttribute("mp", Partio::VECTOR, 1);
    damageH = parts->addAttribute("Dp", Partio::VECTOR, 1);
    spH = parts->addAttribute("sp", Partio::VECTOR, 1);
    for (size_t iter = 0; iter < particles.size(); ++iter) {
        int idx = parts->addParticle();
        float* p = parts->dataWrite<float>(posH, idx);
        float* vel = parts->dataWrite<float>(velH, idx);
        float* DG = parts->dataWrite<float>(dgH, idx);
        float* mp = parts->dataWrite<float>(massH, idx);
        float* Dp = parts->dataWrite<float>(damageH, idx);
        float* Sp = parts->dataWrite<float>(spH, idx);

        for (int k = 0; k < 3; k++)
            p[k] = (T)0;
        for (int k = 0; k < dim; k++)
            p[k] = particles[iter](k);

        for (int k = 0; k < 3; k++)
            vel[k] = (T)0;
        for (int k = 0; k < dim; ++k)
            vel[k] = velocities[iter](k);

        for (int k = 0; k < 3; k++)
            DG[k] = (T)0;
        for (int k = 0; k < dim; ++k)
            DG[k] = damageGradients[iter](k);

        mp[0] = masses[iter];
        Dp[0] = damage[iter];
        Sp[0] = sp[iter];
    }
    Partio::write(particleFile.c_str(), *parts);
    parts->release();
}

template <class T, int dim>
void writeDFGMPMNodesToPartio(const std::string& particleFile,
    const Field<Vector<T, dim>>& particles,
    const Field<Vector<T, dim>>& damageGradients,
    const Field<Vector<T, dim>>& v1,
    const Field<Vector<T, dim>>& v2,
    const Field<Vector<T, dim>>& fct1,
    const Field<Vector<T, dim>>& fct2,
    const std::vector<T>& m1,
    const std::vector<T>& m2,
    const std::vector<T>& sep1,
    const std::vector<T>& sep2,
    const std::vector<int>& separable)
{
    Partio::ParticlesDataMutable* parts = Partio::create();
    Partio::ParticleAttribute posH, dgH, v1H, v2H, fct1H, fct2H, m1H, m2H, sep1H, sep2H, sepH;
    posH = parts->addAttribute("position", Partio::VECTOR, 3);
    dgH = parts->addAttribute("DG", Partio::VECTOR, 3);
    v1H = parts->addAttribute("v1", Partio::VECTOR, 3);
    v2H = parts->addAttribute("v2", Partio::VECTOR, 3);
    fct1H = parts->addAttribute("fct1", Partio::VECTOR, 3);
    fct2H = parts->addAttribute("fct2", Partio::VECTOR, 3);
    m1H = parts->addAttribute("m1", Partio::VECTOR, 1);
    m2H = parts->addAttribute("m2", Partio::VECTOR, 1);
    sep1H = parts->addAttribute("sep1", Partio::VECTOR, 1);
    sep2H = parts->addAttribute("sep2", Partio::VECTOR, 1);
    sepH = parts->addAttribute("separable", Partio::VECTOR, 1);
    for (size_t iter = 0; iter < particles.size(); ++iter) {
        int idx = parts->addParticle();
        float* p = parts->dataWrite<float>(posH, idx);
        float* DG = parts->dataWrite<float>(dgH, idx);
        float* _v1 = parts->dataWrite<float>(v1H, idx);
        float* _v2 = parts->dataWrite<float>(v2H, idx);
        float* _fct1 = parts->dataWrite<float>(fct1H, idx);
        float* _fct2 = parts->dataWrite<float>(fct2H, idx);
        float* _m1 = parts->dataWrite<float>(m1H, idx);
        float* _m2 = parts->dataWrite<float>(m2H, idx);
        float* _sep1 = parts->dataWrite<float>(sep1H, idx);
        float* _sep2 = parts->dataWrite<float>(sep2H, idx);
        float* sep = parts->dataWrite<float>(sepH, idx);

        for (int k = 0; k < 3; k++)
            p[k] = (T)0;
        for (int k = 0; k < dim; k++)
            p[k] = particles[iter](k);

        for (int k = 0; k < 3; k++)
            DG[k] = (T)0;
        for (int k = 0; k < dim; ++k)
            DG[k] = damageGradients[iter](k);

        for (int k = 0; k < 3; k++)
            _v1[k] = (T)0;
        for (int k = 0; k < dim; ++k)
            _v1[k] = v1[iter](k);

        for (int k = 0; k < 3; k++)
            _v2[k] = (T)0;
        for (int k = 0; k < dim; ++k)
            _v2[k] = v2[iter](k);

        for (int k = 0; k < 3; k++)
            _fct1[k] = (T)0;
        for (int k = 0; k < dim; ++k)
            _fct1[k] = fct1[iter](k);

        for (int k = 0; k < 3; k++)
            _fct2[k] = (T)0;
        for (int k = 0; k < dim; ++k)
            _fct2[k] = fct2[iter](k);

        _m1[0] = m1[iter];
        _m2[0] = m2[iter];
        _sep1[0] = sep1[iter];
        _sep2[0] = sep2[iter];
        sep[0] = separable[iter];
    }
    Partio::write(particleFile.c_str(), *parts);
    parts->release();
}
} // namespace Bow
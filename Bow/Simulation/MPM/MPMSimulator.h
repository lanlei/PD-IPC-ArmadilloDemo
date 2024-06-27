#pragma once

#include <Bow/Simulation/MPM/MPMTransfer.h>
#include <Bow/Simulation/MPM/ElasticityOp.h>
#include <Bow/Simulation/MPM/PlasticityOp.h>
#include <Bow/Simulation/MPM/MPMGrid.h>
#include <Bow/Simulation/MPM/BackwardEuler.h>
#include <Bow/Simulation/MPM/HodgepodgeEnergy.h>
#include <Bow/Simulation/PhysicallyBasedSimulator.h>
#include <Bow/IO/partio.h>
#include <Bow/IO/tetwild.h>
#include <Bow/IO/vtk.h>
#include <Bow/Geometry/AnalyticalLevelSet.h>
#include <Bow/Utils/Serialization.h>

namespace Bow::MPM {

template <class T, int dim>
class MPMSimulator : virtual public PhysicallyBasedSimulator<T, dim> {
public:
    using TV = Vector<T, dim>;
    using TM = Matrix<T, dim, dim>;
    Field<TV> m_X;
    Field<TV> m_V;
    Field<TM> m_C;
    std::vector<T> m_mass;

    SERIALIZATION_REGISTER(m_X)
    SERIALIZATION_REGISTER(m_V)
    SERIALIZATION_REGISTER(m_C)
    SERIALIZATION_REGISTER(m_mass)

    Field<TM> stress;
    MPMGrid<T, dim> grid;
    T dx = 0.02;
    T ppc = (T)(1 << dim);
    TV gravity = TV::Zero();
    bool symplectic = true;
    bool backward_euler = true;
    bool apic = true;
    bool dump_F_for_meshing = false;
    std::string output_directory = "mpm_output/";

    std::vector<std::shared_ptr<ElasticityOp<T, dim>>> elasticity_models;
    std::vector<std::shared_ptr<PlasticityOp<T, dim>>> plasticity_models;
    Geometry::BoundaryConditionManager<T, dim> BC;

    void add_particles(std::shared_ptr<ElasticityOp<T, dim>> model, const TV& min_corner, const TV& max_corner, const TV& velocity = TV::Zero(), T density = 1000.)
    {
        // sample particles
        T vol = dim == 2 ? dx * dx / 4 : dx * dx * dx / 8;
        T interval = dx / std::pow((T)ppc, (T)1 / dim);
        Vector<int, dim> region = ((max_corner - min_corner) / interval).template cast<int>();
        int start = m_X.size();
        iterateRegion(region, [&](const Vector<int, dim>& offset) {
            TV position = min_corner + offset.template cast<T>() * interval;
            position += TV::Ones() * 0.5 * interval + TV::Random() * 0.5 * interval;
            m_X.push_back(position);
            m_V.push_back(velocity);
            m_C.push_back(TM::Zero());
            m_mass.push_back(density * vol);
            stress.push_back(TM::Zero());
        });
        int end = m_X.size();
        model->append(start, end, vol);
    }

    void add_particles_from_tetwild(std::shared_ptr<ElasticityOp<T, dim>> model, const std::string mesh_path, const std::string vtk_path = "tet.vtk", const TV& center = TV::Zero(), const TV& velocity = TV::Zero(), T density = 1000.)
    {
        Field<Vector<T, 3>> X;
        Field<Vector<int, 4>> cells;
        IO::read_mesh(mesh_path, X, cells);
        IO::write_vtk(vtk_path, X, cells, false);
        T vol = dim == 2 ? dx * dx / 4 : dx * dx * dx / 8;

        T total_volume = 0;
        for (size_t i = 0; i < cells.size(); i++) {
            TV p0 = X[cells[i](0)], p1 = X[cells[i](1)],
               p2 = X[cells[i](2)], p3 = X[cells[i](3)];
            Matrix<T, 4, 4> A;
            A << 1, p0(0), p0(1), p0(2), 1, p1(0), p1(1), p1(2), 1, p2(0), p2(1), p2(2), 1, p3(0), p3(1), p3(2);
            T temp = A.determinant() / (T)6;
            total_volume += (temp > 0 ? temp : (-temp));
        }
        T total_mass = total_volume * density;
        T mass_per_particle = total_mass / (T)X.size();

        int start = m_X.size();
        for (size_t i = 0; i < X.size(); i++) {
            m_X.push_back(X[i] + center);
            m_V.push_back(velocity);
            m_C.push_back(TM::Zero());
            m_mass.push_back(mass_per_particle);
            stress.push_back(TM::Zero());
        }
        int end = m_X.size();
        model->append(start, end, vol);
    }

    template <class ELASTICITY>
    std::shared_ptr<ElasticityOp<T, dim>> create_elasticity(ELASTICITY* e)
    {
        elasticity_models.push_back(std::shared_ptr<ELASTICITY>(e));
        return elasticity_models.back();
    }

    template <class PLASTICITY>
    void create_plasticity(PLASTICITY* e)
    {
        plasticity_models.push_back(std::shared_ptr<PLASTICITY>(e));
    }

    template <class BOUNDARY_CONDITION>
    void add_boundary_condition(BOUNDARY_CONDITION* bc)
    {
        BC.add(std::shared_ptr<BOUNDARY_CONDITION>(bc));
    }

    // https://stackoverflow.com/questions/47333843/using-initializer-list-for-a-struct-with-inheritance
    void advance(T dt) override
    {
        grid.sortParticles(m_X, dx);
        if (symplectic) {
            for (auto& model : elasticity_models)
                model->compute_stress(stress);
            ParticlesToGridOp<T, dim, true> p2g{ {}, m_X, m_V, m_mass, m_C, stress, grid, dx, dt };
            p2g();
            BoundaryConditionUpdateOp<T, dim> bc_update{ {}, grid, gravity, BC, dx, dt };
            bc_update();
        }
        else {
            ParticlesToGridOp<T, dim, false> p2g{ {}, m_X, m_V, m_mass, m_C, stress, grid, dx, dt };
            p2g();
            Bow::MPM::BackwardEulerUpdateOp<T, dim, int> implicit_mpm(grid, BC, dx, dt);
            Bow::MPM::HodgepodgeEnergy<T, dim> energy(grid, gravity, dx, dt,
                implicit_mpm.BC_basis, implicit_mpm.BC_order, m_X, elasticity_models);
            implicit_mpm.m_energy_terms.push_back(&energy);
            if (backward_euler) {
                energy.tsMethod = HodgepodgeEnergy<T, dim>::BE;
                implicit_mpm.tsMethod = BackwardEulerUpdateOp<T, dim>::BE;
            }
            else {
                energy.tsMethod = HodgepodgeEnergy<T, dim>::NM;
                implicit_mpm.tsMethod = BackwardEulerUpdateOp<T, dim>::NM;
            }
            implicit_mpm();
        }
        GridToParticlesOp<T, dim> g2p{ {}, m_X, m_V, m_C, grid, dx, dt };
        g2p(apic);
        for (auto& model : elasticity_models)
            model->evolve_strain(g2p.m_gradXp);
        for (auto& model : plasticity_models)
            model->project_strain();
    }

    void dump_output(int frame_num) override
    {
        Bow::writePositionVectorToPartio(output_directory + std::to_string(frame_num) + ".bgeo", m_X);

        if (dump_F_for_meshing) {
            Field<TM> m_Fs;
            m_Fs.resize(m_X.size(), TM::Identity());
            for (auto& model : elasticity_models)
                model->collect_strain(m_Fs);
            Bow::IO::write_meshing_data("mpm_output/" + std::to_string(frame_num) + ".dat", m_X, m_Fs);
        }
    }
};

} // namespace Bow::MPM
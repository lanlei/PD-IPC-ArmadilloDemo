#pragma once

#include <Bow/Simulation/PhysicallyBasedSimulator.h>
#include <Bow/Simulation/FEM/FEMEnergies.h>
#include <Bow/Simulation/FEM/ForwardOp.h>
#include <Bow/IO/vtk.h>
#include <Bow/IO/ply.h>
#include <Bow/Utils/Serialization.h>
#include <Bow/Geometry/AnalyticalLevelSet.h>
#include <Bow/Geometry/Query.h>

namespace Bow::FEM {

template <class T, int dim>
class FEMSimulator : public PhysicallyBasedSimulator<T, dim> {
public:
    using TV = Vector<T, dim>;
    using TM = Matrix<T, dim, dim>;

    // constant
    TV gravity = TV::Zero();
    Field<TV> m_X; // material coordinate;
    Field<T> m_mass; // node mass
    Field<T> m_density; // element density
    Field<T> m_vol; // element volume
    Field<Vector<int, dim + 1>> m_elem; // element vertex indices
    Field<TM> m_IB;
    Field<T> m_mu;
    Field<T> m_lam;
    std::unordered_map<int, std::vector<std::pair<int, int>>> m_obj_divider;

    // state
    Field<TV> m_x;
    Field<TV> m_v;
    Field<TV> m_a;

    // control
    Field<TM> BC_basis;
    std::vector<int> BC_order;
    Field<TV> BC_target;
    Field<TV> m_f;

    // intermediate
    Field<TV> m_x_tilde;

    SERIALIZATION_REGISTER(m_x)
    SERIALIZATION_REGISTER(m_v)
    SERIALIZATION_REGISTER(m_a)
    SERIALIZATION_REGISTER(BC_basis)
    SERIALIZATION_REGISTER(BC_order)
    SERIALIZATION_REGISTER(BC_target)
    SERIALIZATION_REGISTER(m_f)
    SERIALIZATION_REGISTER(m_x_tilde)

    std::vector<std::shared_ptr<EnergyOp<T, dim>>> m_energy_terms;
    std::shared_ptr<Bow::FEM::BackwardEulerUpdateOp<T, dim>> m_advance_op;

    bool backward_euler = true;
    bool static_sim = false;
    T tol = 1e-3;

    // level-set based boundary condition
    Geometry::BoundaryConditionManager<T, dim> BC;

    // index based boundary condition
    std::map<int, Field<TV>> node_velocities;

    virtual void append(const Field<TV>& x, const Field<Vector<int, dim + 1>>& elem, const int type, const T E, const T nu, const T density)
    {
        Field<Vector<int, dim + 1>> shifted_elem = elem;
        Eigen::Map<Matrix<int, dim + 1, Eigen::Dynamic>>(&(shifted_elem[0][0]), dim + 1, shifted_elem.size()).array() += m_X.size();
        T mu, lam;
        std::tie(mu, lam) = ConstitutiveModel::lame_paramters(E, nu);

        m_density.resize(m_density.size() + elem.size(), density);
        m_X.insert(m_X.end(), x.begin(), x.end());
        m_elem.insert(m_elem.end(), shifted_elem.begin(), shifted_elem.end());
        m_mu.resize(m_mu.size() + elem.size(), mu);
        m_lam.resize(m_lam.size() + elem.size(), lam);
        m_obj_divider[type].push_back(std::make_pair(m_elem.size() - elem.size(), m_elem.size()));

        m_x.insert(m_x.end(), x.begin(), x.end());
        m_v.resize(m_v.size() + x.size(), Vector<T, dim>::Zero());
        m_a.resize(m_a.size() + x.size(), Vector<T, dim>::Zero());
        m_x_tilde.insert(m_x_tilde.end(), x.begin(), x.end());
        BC_basis.resize(m_x.size(), TM::Identity());
        BC_order.resize(m_x.size(), 0);
        BC_target.resize(m_x.size(), TV::Zero());
        m_f.resize(m_x.size(), TV::Zero());
    }

    /** initialize energy terms and operators */
    virtual void initialize() override
    {
        Bow::FEM::InitializeOp<T, dim> fem_initialize{ {}, m_X, m_elem, m_density, m_mass, m_vol, m_IB };
        fem_initialize();
        m_energy_terms.clear();
        if (!static_sim)
            m_energy_terms.push_back(std::make_shared<InertialEnergyOp<T, dim>>(m_mass, m_x_tilde));
        m_energy_terms.push_back(std::make_shared<StaticForceEnergy<T, dim>>(m_f, m_mass, gravity));
        if (m_obj_divider.find(Bow::ConstitutiveModel::FIXED_COROTATED) != m_obj_divider.end()) {
            auto energy_ptr = std::make_shared<FixedCorotatedEnergyOp<T, dim>>(m_elem, m_vol, m_mu, m_lam, m_IB, m_obj_divider[Bow::ConstitutiveModel::FIXED_COROTATED]);
            m_energy_terms.push_back(energy_ptr);
        }
        if (m_obj_divider.find(Bow::ConstitutiveModel::NEO_HOOKEAN) != m_obj_divider.end()) {
            auto energy_ptr = std::make_shared<NeoHookeanEnergyOp<T, dim>>(m_elem, m_vol, m_mu, m_lam, m_IB, m_obj_divider[Bow::ConstitutiveModel::NEO_HOOKEAN]);
            m_energy_terms.push_back(energy_ptr);
        }
        if (m_obj_divider.find(Bow::ConstitutiveModel::LINEAR_ELASTICITY) != m_obj_divider.end()) {
            auto energy_ptr = std::make_shared<LinearElasticityEnergyOp<T, dim>>(m_elem, m_vol, m_mu, m_lam, m_IB, m_obj_divider[Bow::ConstitutiveModel::LINEAR_ELASTICITY]);
            m_energy_terms.push_back(energy_ptr);
        }

        m_advance_op = std::make_shared<Bow::FEM::BackwardEulerUpdateOp<T, dim>>(BC_basis, BC_order, BC_target, m_mass, m_x, m_v, m_a, m_x_tilde);
        for (auto energy_ptr : m_energy_terms) {
            m_advance_op->m_energy_terms.push_back(energy_ptr.get());
        }
        if (backward_euler)
            m_advance_op->tsMethod = BackwardEulerUpdateOp<T, dim>::BE;
        else
            m_advance_op->tsMethod = BackwardEulerUpdateOp<T, dim>::NM;
        m_advance_op->set_ts_weights();
        m_advance_op->tol = tol;
        m_advance_op->line_search = this->line_search;
    }

    template <class BOUNDARY_CONDITION>
    void add_boundary_condition(BOUNDARY_CONDITION* bc)
    {
        BC.add(std::shared_ptr<BOUNDARY_CONDITION>(bc));
    }
    virtual void advance(T dt) override
    {
        this->m_advance_op->dt = dt;
        this->m_advance_op->set_ts_weights();
        BoundaryConditionUpdateOp<T, dim> bc_update{ {}, BC, m_x, BC_basis, BC_order, BC_target, dt };
        bc_update();
        (*m_advance_op)();
    }

    virtual void dump_output(int frame_num) override
    {
        if constexpr (dim == 2)
            IO::write_ply(this->output_directory + "/" + std::to_string(frame_num) + ".ply", m_x, m_elem);
        else {
            IO::write_vtk(this->output_directory + "/" + std::to_string(frame_num) + ".vtk", m_x, m_elem);
            IO::write_ply(this->output_directory + "/" + std::to_string(frame_num) + ".ply", m_x);
        }
    }
};

// TODO: 3D
template <class T, int dim>
class IPCSimulator : public FEMSimulator<T, dim> {
public:
    using Base = FEMSimulator<T, dim>;
    using TV = Vector<T, dim>;
    using TM = Matrix<T, dim, dim>;
    std::shared_ptr<IPC::IpcEnergyOp<T, dim>> ipc_energy;
    Field<int> m_boundary_points;
    Field<Vector<int, 2>> m_boundary_edges;
    Field<Vector<int, 3>> m_boundary_faces;
    std::map<int, T> m_snode_area;
    T dHat = 1e-3;
    T kappa = 1e4;

    // friction related
    T mu = 0;
    T epsv = 1e-7;
    bool lag = false;

    virtual void append(const Field<TV>& x, const Field<Vector<int, dim + 1>>& elem, const int type, const T E, const T nu, const T density) override
    {
        Field<int> new_boundary_points;
        Field<Vector<int, 2>> new_boundary_edges;
        Field<Vector<int, 3>> new_boundary_faces;
        if constexpr (dim == 2)
            Geometry::find_boundary(elem, new_boundary_edges, new_boundary_points);
        else
            Geometry::find_boundary(elem, new_boundary_faces, new_boundary_edges, new_boundary_points);
        new_boundary_points += int(this->m_x.size());
        new_boundary_edges += Vector<int, 2>(this->m_x.size(), this->m_x.size());
        new_boundary_faces += Vector<int, 3>(this->m_x.size(), this->m_x.size(), this->m_x.size());

        Base::append(x, elem, type, E, nu, density);

        m_boundary_edges.insert(m_boundary_edges.end(), new_boundary_edges.begin(), new_boundary_edges.end());
        m_boundary_points.insert(m_boundary_points.end(), new_boundary_points.begin(), new_boundary_points.end());
        m_boundary_faces.insert(m_boundary_faces.end(), new_boundary_faces.begin(), new_boundary_faces.end());
        if constexpr (dim == 2) {
            for (const auto& eI : new_boundary_edges) {
                T eLen = (this->m_x[eI[0]] - this->m_x[eI[1]]).norm();
                m_snode_area[eI[0]] += eLen / 2;
                m_snode_area[eI[1]] += eLen / 2;
            }
        }
        else {
            for (const auto& fI : new_boundary_faces) {
                T face_area = (this->m_x[fI[1]] - this->m_x[fI[0]]).cross(this->m_x[fI[2]] - this->m_x[fI[0]]).norm() / 2;
                m_snode_area[fI[0]] += face_area / 3;
                m_snode_area[fI[1]] += face_area / 3;
                m_snode_area[fI[2]] += face_area / 3;
            }
        }
        // TODO edge pair weight
    }

    virtual void initialize() override
    {
        Base::initialize();
        if constexpr (dim == 2)
            ipc_energy = std::make_shared<IPC::IpcEnergyOp<T, 2>>(m_boundary_points, m_boundary_edges, this->m_mass, m_snode_area);
        else
            ipc_energy = std::make_shared<IPC::IpcEnergyOp<T, 3>>(m_boundary_points, m_boundary_edges, m_boundary_faces, this->m_X, this->m_mass, m_snode_area);
        this->m_energy_terms.push_back(ipc_energy);
        ipc_energy->dHat = dHat;
        ipc_energy->kappa = kappa;
        this->m_advance_op->m_energy_terms.push_back(ipc_energy.get());
        this->m_advance_op->set_ts_weights();
    }

    virtual void advance(T dt) override
    {
        this->m_advance_op->dt = dt;
        this->m_advance_op->set_ts_weights();
        if constexpr (dim == 2) {
            ipc_energy->update_weight_and_xhat(this->m_x, this->m_v, this->m_a, dt, this->m_advance_op->tsParam[this->m_advance_op->tsMethod]);
            ipc_energy->initialize_friction(mu, epsv, dt);
            if (!lag)
                Base::advance(dt);
            else {
                this->m_advance_op->set_ts_weights();
                this->m_advance_op->initialize_acceleration();
                BoundaryConditionUpdateOp<T, dim> bc_update{ {}, this->BC, this->m_x, this->BC_basis, this->BC_order, this->BC_target, dt };
                bc_update();
                this->m_advance_op->update_transformation_matrix();
                this->m_advance_op->update_predictive_pos();
                Bow::Vector<T, Eigen::Dynamic> new_x = Bow::to_vec(this->m_x);
                this->m_advance_op->optimize(new_x);
                while (true && !this->static_sim) {
                    ipc_energy->initialize_friction(mu, epsv, dt);
                    int iter_count = this->m_advance_op->optimize(new_x);
                    if (iter_count == 1 || iter_count == this->m_advance_op->max_iter)
                        break;
                }
                this->m_advance_op->advance(new_x);
            }
        }
        else {
            Base::advance(dt);
        }
    }

    virtual void dump_output(int frame_num) override
    {
        if constexpr (dim == 2)
            IO::write_ply(this->output_directory + "/" + std::to_string(frame_num) + ".ply", this->m_x, this->m_elem);
        else {
            IO::write_ply(this->output_directory + "/" + std::to_string(frame_num) + ".ply", this->m_x, m_boundary_faces);
        }
    }
};

} // namespace Bow::FEM
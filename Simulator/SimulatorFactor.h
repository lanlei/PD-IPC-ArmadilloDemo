#pragma once
//#include "ProjectiveDynamics\ProjectiveDynamicsSimulator.h"
//#include "FiniteElementMethod\FemSimulator.h"
//#include "FiniteElementMethod\Reduced\FemReducedSimulator.h"
//#include "ReducedIpc\ReducedIPCSimulator.h"
//#include "ReducedIpc\BaseReducedIPCSimulator.h"
//#include "ReducedIpc\BarbarianShipDemo\BarbarianShipDemoSimulator.h"
//#include "ReducedIpc\CactusDemo\CactusDemoSimulator.h"
//#include "ReducedIpc\ArmadilloDemo\ArmadilloDemoSimulator.h"
#include "ObjSequenceRender.h"
#include "PD-IPC\PdIpcSimulator.h"

static BaseSimulator* createSimulator(std::string& simName, RunPlatform platform)
{
	/*if (simName == "fem_simulator")
	{
		return new FiniteElementMethod::FemSimulator(simName, platform);
	}
	else if (simName == "fem_reduced_simulator")
	{
		return new FiniteElementMethod::FemReducedSimulator(simName, platform);
	}
	else if (simName == "reduced_ipc_simulator")
	{
		return new ReducedIPC::ReducedIPCSimulator(simName, platform);
	}
	else if (simName == "base_reduced_ipc_simulator")
	{
		return new ReducedIPC::BaseReducedIPCSimulator(simName, platform);
	}
	else if (simName == "barbarian_ship_demo_simulator")
	{
		return new ReducedIPC::BarbarianShipDemoSimulator(simName, platform);
	}
	else if (simName == "cactus_demo_simulator")
	{
		return new ReducedIPC::CactusDemoSimulator(simName, platform);
	}
	else if (simName == "armadillo_demo_simulator")
	{
		return new ReducedIPC::ArmadilloDemoSimulator(simName, platform);
	}*/
	

	if (simName == "render_obj_sequence")
	{
		return new ObjSequenceRender();	
	}
	else if (simName == "pd_ipc_simulator")
	{
		return new PD_IPC::PdIpcSimulator(simName, platform);
	}
	else
	{
		return new BaseSimulator("base_simulator", platform);
	}
}

static bool readSimulatorFromConfigFile(BaseSimulator* sim, const std::string filename, TiXmlElement* item)
{
	std::string simName = sim->getSimulatorName();
	/*if (simName == "fem_simulator")
	{
		return dynamic_cast<FiniteElementMethod::FemSimulator*>(sim)->readSimulatorFromConfigFile(filename, item);
	}
	else if (simName == "fem_reduced_simulator")
	{
		return dynamic_cast<FiniteElementMethod::FemReducedSimulator*>(sim)->readSimulatorFromConfigFile(filename, item);
	}
	else if (simName == "reduced_ipc_simulator")
	{
		return dynamic_cast<ReducedIPC::ReducedIPCSimulator*>(sim)->readSimulatorFromConfigFile(filename, item);
	}
	else if (simName == "base_reduced_ipc_simulator")
	{
		return dynamic_cast<ReducedIPC::BaseReducedIPCSimulator*>(sim)->readSimulatorFromConfigFile(filename, item);
	}
	else if (simName == "barbarian_ship_demo_simulator")
	{
		return dynamic_cast<ReducedIPC::BarbarianShipDemoSimulator*>(sim)->readSimulatorFromConfigFile(filename, item);
	}
	else if (simName == "cactus_demo_simulator")
	{
		return dynamic_cast<ReducedIPC::CactusDemoSimulator*>(sim)->readSimulatorFromConfigFile(filename, item);
	}
	else if (simName == "armadillo_demo_simulator")
	{
		return dynamic_cast<ReducedIPC::ArmadilloDemoSimulator*>(sim)->readSimulatorFromConfigFile(filename, item);
	}
	*/	
	if (simName == "render_obj_sequence")
	{
		return dynamic_cast<ObjSequenceRender*>(sim)->readSimulatorFromConfigFile(filename, item);
	}
	else if (simName == "pd_ipc_simulator")
	{
		return dynamic_cast<PD_IPC::PdIpcSimulator*>(sim)->readSimulatorFromConfigFile(filename, item);
	}
	else if (simName == "base_simulator")
	{
		return sim->readSimulatorFromConfigFile(filename, item);
	}

	return false;
}

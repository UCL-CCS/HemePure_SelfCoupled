
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "net/mpi.h"
#include "net/IOCommunicator.h"
#include "configuration/CommandLine.h"
#include "SimulationMaster.h"

extern "C" void split_comm_world_(
		MPI_Fint *FCOMM_UNIVERSE, MPI_Fint *FCOMM_WORLD, int *MY_UNIVERSE_RANK, int *MY_WORLD_RANK, int *INSTANCE);

int main(int argc, char *argv[])
{
	// main function needed to perform the entire simulation. Some
	// simulation parameters and performance statistics are output on
	// standard output.

	// Bring up MPI.
	hemelb::net::MpiEnvironment mpi(argc, argv);
	try
	{
		MPI_Fint FCOMM_UNIVERSE, FCOMM_WORLD;
		MPI_Comm CCOMM_UNIVERSE, CCOMM_WORLD;
		MPI_Status istat;

		int MY_UNIVERSE_RANK, MY_WORLD_RANK;;

		std::string str = argv[0];
		int LASTINDX = str.find_last_not_of("0123456789");
		int INSTANCE = std::stoi(str.substr(LASTINDX + 1));

		split_comm_world_(&FCOMM_UNIVERSE, &FCOMM_WORLD, &MY_UNIVERSE_RANK, &MY_WORLD_RANK, &INSTANCE);

		// Convert from a Fortran handle (which is an integer) to a C handle.
		CCOMM_WORLD    = MPI_Comm_f2c(FCOMM_WORLD);
		CCOMM_UNIVERSE = MPI_Comm_f2c(FCOMM_UNIVERSE);

		hemelb::net::MpiCommunicator commWorld    = hemelb::net::MpiCommunicator::Get(CCOMM_WORLD);
//		hemelb::net::MpiCommunicator commUniverse = hemelb::net::MpiCommunicator::Get(CCOMM_UNIVERSE);

		hemelb::log::Logger::Init(MY_WORLD_RANK);

		hemelb::net::IOCommunicator hemelbCommunicatorWorld(
				commWorld);
//		hemelb::net::IOCommunicator hemelbCommunicatorUniverse(
//				commUniverse);

		try
		{
			// Parse command line.
			hemelb::configuration::CommandLine options = hemelb::configuration::CommandLine(argc, argv);

			// Prepare main simulation object...
			SimulationMaster master = SimulationMaster(options, hemelbCommunicatorWorld, CCOMM_UNIVERSE);

			// ..and run it.
			master.RunSimulation();
		}

		// Interpose this catch to print usage before propagating the error.
		catch (hemelb::configuration::CommandLine::OptionError& e)
		{
			hemelb::log::Logger::Log<hemelb::log::Critical, hemelb::log::Singleton>(hemelb::configuration::CommandLine::GetUsage());
			throw;
		}
	}
	catch (std::exception& e)
	{
		hemelb::log::Logger::Log<hemelb::log::Critical, hemelb::log::OnePerCore>(e.what());
		mpi.Abort(-1);
	}
	// MPI gets finalised by MpiEnv's d'tor.
}

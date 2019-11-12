
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "SimulationMaster.h"
#include "configuration/SimConfig.h"
#include "extraction/PropertyActor.h"
#include "extraction/LbDataSourceIterator.h"
#include "io/writers/xdr/XdrFileWriter.h"
#include "util/utilityFunctions.h"
#include "geometry/GeometryReader.h"
#include "geometry/LatticeData.h"
#include "util/fileutils.h"
#include "log/Logger.h"
#include "lb/HFunction.h"
#include "io/xml/XmlAbstractionLayer.h"
#include "colloids/ColloidController.h"
#include "net/BuildInfo.h"
#include "net/IOCommunicator.h"
#include "colloids/BodyForces.h"
#include "colloids/BoundaryConditions.h"

#include <boost/algorithm/string.hpp>
#include <cstdlib>
#include <limits>
#include <map>

/**
 * Constructor for the SimulationMaster class.
 *
 * Initialises member variables including the network topology
 * object.
 */

double hemelb::tau_;

SimulationMaster::SimulationMaster(hemelb::configuration::CommandLine & options,
		const hemelb::net::IOCommunicator& ioComm,
		const MPI_Comm& unComm) :
	ioComms(ioComm), unComms(unComm), timings(ioComm), build_info(), communicationNet(ioComm) {
	timings[hemelb::reporting::Timers::total].Start();

	latticeData = NULL;

	colloidController = NULL;
	latticeBoltzmannModel = NULL;
	propertyDataSource = NULL;
	propertyExtractor = NULL;
	simulationState = NULL;
	stepManager = NULL;
	netConcern = NULL;
	neighbouringDataManager = NULL;
	imagesPerSimulation = options.NumberOfImages();

	fileManager = new hemelb::io::PathManager(options, IsCurrentProcTheIOProc(), GetProcessorCount());
	simConfig = hemelb::configuration::SimConfig::New(fileManager->GetInputFile());
	unitConverter = &simConfig->GetUnitConverter();
	monitoringConfig = simConfig->GetMonitoringConfiguration();

	fileManager->SaveConfiguration(simConfig);
	Initialise();

	if (IsCurrentProcTheIOProc()) {
		reporter = new hemelb::reporting::Reporter(
			fileManager->GetReportPath(),
			fileManager->GetInputFile());
		reporter->AddReportable(&build_info);

		if (monitoringConfig->doIncompressibilityCheck) {
			reporter->AddReportable(incompressibilityChecker);
		}

		reporter->AddReportable(&timings);
		reporter->AddReportable(latticeData);
		reporter->AddReportable(simulationState);
	}
}

/**
 * Destructor for the SimulationMaster class.
 *
 * Deallocates dynamically allocated memory to contained classes.
 */
SimulationMaster::~SimulationMaster() {

	delete latticeData;
	delete colloidController;
	delete latticeBoltzmannModel;
	delete inletValues;
	delete outletValues;
	delete propertyExtractor;
	delete propertyDataSource;
	delete stabilityTester;
	delete entropyTester;
	delete simulationState;
	delete incompressibilityChecker;
	delete neighbouringDataManager;

	delete simConfig;
	delete fileManager;

	if (IsCurrentProcTheIOProc()) {
		delete reporter;
	}

	delete stepManager;
	delete netConcern;
}

/**
 * Returns true if the current processor is the dedicated I/O
 * processor.
 */
bool SimulationMaster::IsCurrentProcTheIOProc() {
	return ioComms.OnIORank();
}

/**
 * Returns the number of processors involved in the simulation.
 */
int SimulationMaster::GetProcessorCount() {
	return ioComms.Size();
}

/**
 * Initialises various elements of the simulation.
 */
void SimulationMaster::Initialise() {

	simulationState = new hemelb::lb::SimulationState(
		simConfig->GetTimeStepLength(),
		simConfig->GetTotalTimeSteps());

	timings[hemelb::reporting::Timers::latDatInitialise].Start();

	int sizeUni, rankUni;
	MPI_Comm_size(unComms, &sizeUni);
	MPI_Comm_rank(unComms, &rankUni);

	// Use a reader to read in the file.
	hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("(%06i) size of this universe: %05i", getpid(), sizeUni);
	hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("(%06i) size of this world:    %05i", getpid(), ioComms.Size());
	hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("INITIALISE");
	hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("----------");
	hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("--> loading input and decomposing geometry");
	hemelb::geometry::GeometryReader reader(
		latticeType::GetLatticeInfo(),
		timings, ioComms);
	hemelb::geometry::Geometry readGeometryData =
		reader.LoadAndDecompose(simConfig->GetDataFilePath());

	// Create a new lattice based on that info and return it.
	hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("--> lattice data");
	latticeData = new hemelb::geometry::LatticeData(latticeType::GetLatticeInfo(),
			readGeometryData,
			ioComms);

	timings[hemelb::reporting::Timers::latDatInitialise].Stop();

	hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("--> neighbouring data manager");
	neighbouringDataManager = new hemelb::geometry::neighbouring::NeighbouringDataManager(*latticeData,
				latticeData->GetNeighbouringData(),
				communicationNet);

	hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("--> lattice-Boltzmann model");
	latticeBoltzmannModel = new hemelb::lb::LBM<latticeType>(simConfig,
			&communicationNet,
			latticeData,
			simulationState,
			timings,
			neighbouringDataManager);

	hemelb::lb::MacroscopicPropertyCache& propertyCache = latticeBoltzmannModel->GetPropertyCache();

	hemelb::tau_ = latticeBoltzmannModel->GetLbmParams()->GetTau();
	
	if (simConfig->HasColloidSection()) {
		hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("--> colloid section present");

		timings[hemelb::reporting::Timers::colloidInitialisation].Start();
		hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("----> loading colloid configuration");
		std::string colloidConfigPath = simConfig->GetColloidConfigPath();
		hemelb::io::xml::Document xml(colloidConfigPath);

		hemelb::colloids::BodyForces::InitBodyForces(xml);
		hemelb::colloids::BoundaryConditions::InitBoundaryConditions(latticeData, xml);

		hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("----> initialising colloid controller");
		colloidController =
			new hemelb::colloids::ColloidController(*latticeData,
					*simulationState,
					simConfig,
					readGeometryData,
					xml,
					propertyCache,
					//latticeBoltzmannModel->GetLbmParams(),
					fileManager->GetColloidPath(),
					ioComms,
					timings);
		timings[hemelb::reporting::Timers::colloidInitialisation].Stop();
	}

	stabilityTester = new hemelb::lb::StabilityTester<latticeType>(latticeData,
			&communicationNet,
			simulationState,
			timings,
			monitoringConfig);
	entropyTester = NULL;

	if (monitoringConfig->doIncompressibilityCheck) {
		incompressibilityChecker = new hemelb::lb::IncompressibilityChecker <
		hemelb::net::PhasedBroadcastRegular<> > (latticeData,
				&communicationNet,
				simulationState,
				latticeBoltzmannModel->GetPropertyCache(),
				timings);
	} else {
		incompressibilityChecker = NULL;
	}

	inletValues = new hemelb::lb::iolets::BoundaryValues(hemelb::geometry::INLET_TYPE,
			latticeData,
			simConfig->GetInlets(),
			simulationState,
			ioComms,
			*unitConverter);

	outletValues = new hemelb::lb::iolets::BoundaryValues(hemelb::geometry::OUTLET_TYPE,
			latticeData,
			simConfig->GetOutlets(),
			simulationState,
			ioComms,
			*unitConverter);

	latticeBoltzmannModel->Initialise(inletValues, outletValues, unitConverter);
	latticeBoltzmannModel->SetInitialConditions(ioComms); //JM Checkpoint addition
	neighbouringDataManager->ShareNeeds();
	neighbouringDataManager->TransferNonFieldDependentInformation();

	propertyDataSource =
		new hemelb::extraction::LbDataSourceIterator(latticeBoltzmannModel->GetPropertyCache(),
				*latticeData,
				ioComms.Rank(),
				*unitConverter);

	if (simConfig->PropertyOutputCount() > 0) {

		for (unsigned outputNumber = 0; outputNumber < simConfig->PropertyOutputCount(); ++outputNumber) {
			simConfig->GetPropertyOutput(outputNumber)->filename = fileManager->GetDataExtractionPath()
					+ simConfig->GetPropertyOutput(outputNumber)->filename;
		}

		propertyExtractor = new hemelb::extraction::PropertyActor(*simulationState,
				simConfig->GetPropertyOutputs(),
				*propertyDataSource,
				timings, ioComms, unComms);

		rank2IoletMap = propertyExtractor->GetIoletMap();

		// Invert the map.
		for (int rank = 0; rank < ioComms.Size(); rank++)
		{
			if (rank2IoletMap.find(rank) != rank2IoletMap.end())
			{
				for (std::vector<int>::iterator it = rank2IoletMap.at(rank).begin();
						it != rank2IoletMap.at(rank).end(); ++it)
				{
					iolet2RankMap[(*it)].push_back(rank);
				}
			}
		}

		if (IsCurrentProcTheIOProc() && !rank2IoletMap.empty())
		{
			std::ifstream file(simConfig->GetMapFilePath());
			std::string line = "";

			// Iterate through each line and split the content using delimiter.
			while (getline(file, line))
			{
				std::vector<std::string> vec;
				boost::algorithm::split(vec, line, boost::is_any_of(","));
				couplingMap.push_back(vec);
			}
			file.close();
		}
	}

	imagesPeriod = OutputPeriod(imagesPerSimulation);

	stepManager = new hemelb::net::phased::StepManager(2,
			&timings,
			hemelb::net::separate_communications);
	netConcern = new hemelb::net::phased::NetConcern(communicationNet);
	stepManager->RegisterIteratedActorSteps(*neighbouringDataManager, 0);

	if (colloidController != NULL) {
		stepManager->RegisterIteratedActorSteps(*colloidController, 1);
	}

	stepManager->RegisterIteratedActorSteps(*latticeBoltzmannModel, 1);

	stepManager->RegisterIteratedActorSteps(*inletValues, 1);
	stepManager->RegisterIteratedActorSteps(*outletValues, 1);
	stepManager->RegisterIteratedActorSteps(*stabilityTester, 1);

	if (entropyTester != NULL) {
		stepManager->RegisterIteratedActorSteps(*entropyTester, 1);
	}

	if (monitoringConfig->doIncompressibilityCheck) {
		stepManager->RegisterIteratedActorSteps(*incompressibilityChecker, 1);
	}

	if (propertyExtractor != NULL) {
		stepManager->RegisterIteratedActorSteps(*propertyExtractor, 1);
	}

	stepManager->RegisterCommsForAllPhases(*netConcern);

	MPI_Barrier(unComms);

	hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("-------------------");
	hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("INITIALISE FINISHED");
}

unsigned int SimulationMaster::OutputPeriod(unsigned int frequency) {
	if (frequency == 0) {
		return 1000000000;
	}

	unsigned long roundedPeriod = simulationState->GetTotalTimeSteps() / frequency;
	return hemelb::util::NumericalFunctions::max(1U, (unsigned int) roundedPeriod);
}

void SimulationMaster::HandleActors() {
	stepManager->CallActions();
}

void SimulationMaster::OnUnstableSimulation() {
	LogStabilityReport();
	hemelb::log::Logger::Log<hemelb::log::Warning, hemelb::log::Singleton>("ABORTING :: time step length: %f",
			simulationState->GetTimeStepLength());
	Finalise();
	Abort();
}

/**
 * Begin the simulation.
 */
void SimulationMaster::RunSimulation() {
	hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("SIMULATION STARTING");
	hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("-------------------");
	timings[hemelb::reporting::Timers::simulation].Start();

	while (simulationState->GetTimeStep() <= simulationState->GetTotalTimeSteps()) {
		DoTimeStep();

		if (simulationState->IsTerminating()) {
			break;
		}
	}

	timings[hemelb::reporting::Timers::simulation].Stop();
	Finalise();
}

void SimulationMaster::Finalise() {
	timings[hemelb::reporting::Timers::total].Stop();
	timings.Reduce();

	if (IsCurrentProcTheIOProc()) {
		reporter->FillDictionary();
		reporter->Write();
	}

	// DTMP: Logging output on communication as debug output for now.
	hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::OnePerCore>("sync points: %lld, bytes sent: %lld",
			communicationNet.SyncPointsCounted,
			communicationNet.BytesSent);

	hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("-------------------");
	hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("SIMULATION FINISHED");
}

void SimulationMaster::DoTimeStep() {
	bool writeImage = ((simulationState->GetTimeStep() % imagesPeriod) == 0) ?
					true :
					false;

	if (simulationState->GetTimeStep() % 200 == 0) {
		hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("time step %07i :: write_image_to_disk %i",
				simulationState->GetTimeStep(),
				writeImage);
		LogStabilityReport();
	}

	RecalculatePropertyRequirements();

	HandleActors();
	
	if (simulationState->GetStability() == hemelb::lb::Unstable) {
		OnUnstableSimulation();
	}

	// If the user requested to terminate converged steady flow simulations, mark
	// simulation to be finished at the end of the current timestep.
	if ((simulationState->GetStability() == hemelb::lb::StableAndConverged)
			&& monitoringConfig->convergenceTerminate) {
		LogStabilityReport();
		simulationState->SetIsTerminating(true);
	}

	// Colloid output.
	//if ((simulationState->GetTimeStep() % 100 == 0) && colloidController != NULL)
	//	colloidController->OutputInformation(simulationState->GetTimeStep());

	if (simulationState->GetTimeStep() % FORCE_FLUSH_PERIOD == 0 && IsCurrentProcTheIOProc()) {
		fflush(NULL);
	}
	
	// Coupling HemePure instances.
	if (simulationState->GetTimeStep() > 1 && propertyExtractor != NULL) {

		// Get coupledFieldsThere from other HemePure instance.
		std::unordered_map<int, std::vector<double> > coupledFieldsThere;
		bool coupling = propertyExtractor->GetCoupledFields(coupledFieldsThere);

		if (coupling)
		{
			// Get number of outlets	on lattice with ID = 0.
			int nOutlets = couplingMap.size();

			// Get number of inlets		on lattice with ID = 1.
			std::vector<int> nInlets(nOutlets);
			int cumInlets = 0;
			for (int i = 0; i < nOutlets; ++i)
			{	
				nInlets[i] = (couplingMap[i].size() - 3)/3;
				cumInlets += (couplingMap[i].size() - 3)/3;
			}

			std::vector<int> inletOrder(cumInlets);
			std::vector<int> outletOrder(cumInlets);

			// outletLatticeID outletID inletLatticeID (inletID velocityScaleFactor pressureScaleFactor)
		
			std::vector<std::vector<double> > inletVelocity(nOutlets);
			std::vector<std::vector<double> > inletPressure(nOutlets);
			std::vector<std::vector<double> > inletLBMPressure(nOutlets);
			std::vector<std::vector<double> > outletVelocity(nOutlets);
			std::vector<std::vector<double> > outletPressure(nOutlets);
			std::vector<std::vector<double> > outletLBMPressure(nOutlets);
			
			if ((simulationState->GetTimeStep()/simConfig->GetPropertyOutput(0)->frequency)%2 == 1) //JM to oscillate velocity/pressure transfer - ASSUMES coupling is first property
			{
				// Pressure transfer from mesh 1 to mesh 0
				if (IsCurrentProcTheIOProc() && simConfig->GetLatticeId() == 0)
				{

					// ONLY FOR GOOD FOR 2 MESHES!
					int k=0;
					for (int i = 0; i < nOutlets; ++i)
					{
						double tmpVelocity = 0;
						double tmpPressure = 0;
						
						outletOrder[k] = i;
						k +=1;
						
						for (int j = 0; j < nInlets[i]; ++j)
						{
							// Global IOlet ID.
							int iolet = stoi(couplingMap[i][3*j+3]);
							
							// Voutlet = Vinlet/scaleFactor
							tmpVelocity += coupledFieldsThere.at(iolet)[0]/stod(couplingMap[i][3*j+4]);
																					
							// Poutlet = Pinlet/ (1.0 - scaleFactor)
							//tmpPressure += coupledFieldsThere.at(iolet)[1]/(1.0 - stod(couplingMap[i][3*j+5]));
						}
							
						// Voutlet = Vinlet/scaleFactor
						//outletVelocity[i].push_back(tmpVelocity/nInlets[i]);
							
						// Poutlet = Pinlet/ (1.0 - scaleFactor)
						//outletLBMPressure[i].push_back(tmpPressure/nInlets[i]); //JM original option that uses local LBM pressure values
						outletPressure[i].push_back(0.5*hemelb::BLOOD_DENSITY_Kg_per_m3*tmpVelocity*tmpVelocity/(nInlets[i]*nInlets[i])); //JM option that uses dynamic pressure based on velocity
						
					}
				}

				if (simConfig->GetLatticeId() == 0)
				{
					if (IsCurrentProcTheIOProc())
					{
						// ONLY FOR GOOD FOR 2 MESHES!
						for (int i = 0; i < nOutlets; ++i)
						{
							for (std::unordered_map<int, std::vector<int> >::iterator it1 = iolet2RankMap.begin();
									it1 != iolet2RankMap.end(); ++it1)
							{
								for (std::vector<int>::iterator it2 = it1->second.begin();
										it2 != it1->second.end(); ++it2)
								{
									// Send IOlet index first.
									ioComms.Send(it1->first, *it2);
									// Send pressure for this IOlet.
									ioComms.Send(outletPressure[outletOrder[it1->first]][0], *it2);
									// Send LBM pressure for this IOlet.
									//ioComms.Send(outletLBMPressure[outletOrder[it1->first]][0], *it2);
									// Send Velocity for this IOlet.
									//ioComms.Send(outletVelocity[outletOrder[it1->first]][0], *it2);
								}
							}
						}
					}
					else
					{
						double tmp_pressure, tmp_LBMpressure, tmp_Velocity; int tmp_iolet;
						for (std::vector<int>::iterator it = rank2IoletMap.at(-1).begin();
								it != rank2IoletMap.at(-1).end(); ++it)
						{
							// Receive IOlet index first.
							ioComms.Receive(tmp_iolet,ioComms.GetIORank());
							// Receive pressure for this IOlet.
							ioComms.Receive(tmp_pressure,ioComms.GetIORank());
							// Receive LBM pressure for this IOlet.
							//ioComms.Receive(tmp_LBMpressure,ioComms.GetIORank());
							// Receive Velocity for this IOlet.
							//ioComms.Receive(tmp_Velocity,ioComms.GetIORank());
							//tmp_iolet = 0; tmp_velocity = sin(simulationState->GetTimeStep());
							
							((hemelb::lb::iolets::InOutLetCosine*)outletValues->GetLocalIolet(tmp_iolet))->SetForce(outletValues->GetLocalIolet(tmp_iolet)->GetNormal()*unitConverter->ConvertPressureDifferenceToLatticeUnits(tmp_pressure)*unitConverter->GetVoxelSize()); //Dynamic pressure; tmp_pressure*area_dot_normal*dx==forceDensityvector);
						}
					}
					//if (IsCurrentProcTheIOProc())
					//{
					//	int i = 0;
					//	char hostname[256];
					//	gethostname(hostname, sizeof(hostname));
					//	printf("PID %d (%d) on %s ready for attach\n", getpid(), 1, hostname);
					//	fflush(stdout);
					//	while (0 == i)
					//		sleep(5);
					//}
			
				}
			}
		
			if ((simulationState->GetTimeStep()/simConfig->GetPropertyOutput(0)->frequency)%2 == 0) //JM to oscillate pressure/velocity transfer - ASSUMES coupling is first property
			{
				// Velocity transfer from mesh 0 to mesh 1
				if (IsCurrentProcTheIOProc() && simConfig->GetLatticeId() == 1)
				{

					// ONLY FOR GOOD FOR 2 MESHES!
					int k=0;
					for (int i = 0; i < nOutlets; ++i)
					{
						for (int j = 0; j < nInlets[i]; ++j)
						{
							// Global IOlet ID.
							int iolet = stoi(couplingMap[i][1]);
							inletOrder[k] = j;
							outletOrder[k] = i;
							k +=1;

							// Vinlet = scaleFactor*Voutlet
							inletVelocity[i].push_back(stod(couplingMap[i][3*j+4])*coupledFieldsThere.at(iolet)[0]);

							// Pinlet = Poutlet - scaleFactor*Poutlet = (1.0 - scaleFactor)*Poutlet
							//inletLBMPressure[i].push_back((1.0 - stod(couplingMap[i][3*j+5]))*coupledFieldsThere.at(iolet)[1]);
							//inletPressure[i].push_back(0.5*hemelb::BLOOD_DENSITY_Kg_per_m3*stod(couplingMap[i][3*j+4])*coupledFieldsThere.at(iolet)[0]*stod(couplingMap[i][3*j+4])*coupledFieldsThere.at(iolet)[0]); //JM Dynamic pressure at inlet
					
						}
					}
				}

				if (simConfig->GetLatticeId() == 1)
				{
					if (IsCurrentProcTheIOProc())
					{
						// ONLY FOR GOOD FOR 2 MESHES!
						
						for (std::unordered_map<int, std::vector<int> >::iterator it1 = iolet2RankMap.begin();
								it1 != iolet2RankMap.end(); ++it1)
						{
							for (std::vector<int>::iterator it2 = it1->second.begin();
									it2 != it1->second.end(); ++it2)
							{
							
								// Send IOlet index first.
								ioComms.Send(it1->first, *it2);
								// Send velocity for this IOlet.
								ioComms.Send(inletVelocity[outletOrder[it1->first]][inletOrder[it1->first]], *it2);
								// Send pressure for this IOlet.
								//ioComms.Send(inletLBMPressure[outletOrder[it1->first]][inletOrder[it1->first]], *it2);
							}
						}
						
					}
					else
					{
						double tmp_velocity, tmp_pressure; int tmp_iolet;
						for (std::vector<int>::iterator it = rank2IoletMap.at(-1).begin();
								it != rank2IoletMap.at(-1).end(); ++it)
						{
							// Receive IOlet index first.
							ioComms.Receive(tmp_iolet,ioComms.GetIORank());
							// Receive velocity for this IOlet.
							ioComms.Receive(tmp_velocity,ioComms.GetIORank());
							// Receive velocity for this IOlet.
							//ioComms.Receive(tmp_pressure,ioComms.GetIORank());
							
						
#ifdef HEMELB_USE_VELOCITY_WEIGHTS_FILE
							//JM File Inlets
							((hemelb::lb::iolets::InOutLetFileVelocity*)inletValues->GetLocalIolet(tmp_iolet))->SetMaxSpeed(unitConverter->ConvertSpeedToLatticeUnits(tmp_velocity/(((hemelb::lb::iolets::InOutLetFileVelocity*)inletValues->GetLocalIolet(tmp_iolet))->weightAve)));
#endif
#ifndef HEMELB_USE_VELOCITY_WEIGHTS_FILE
							//JM Parabolic (circular) inlets
							((hemelb::lb::iolets::InOutLetParabolicVelocity*)inletValues->GetLocalIolet(tmp_iolet))->SetMaxSpeed(unitConverter->ConvertSpeedToLatticeUnits(2.0*tmp_velocity));
#endif
						}
					}
					//if (IsCurrentProcTheIOProc())
					//{
					//	int i = 0;
					//	char hostname[256];
					//	gethostname(hostname, sizeof(hostname));
					//	printf("PID %d (%d) on %s ready for attach\n", getpid(), 1, hostname);
					//	fflush(stdout);
					//	while (0 == i)
					//		sleep(5);
					//}
				}
			}
		}
	}
	simulationState->Increment();

	if (simulationState->GetTimeStep()%100000 == 0) {
		MPI_Barrier(ioComms);
	}

}

void SimulationMaster::RecalculatePropertyRequirements() {
	// Get the property cache & reset its list of properties to get.
	hemelb::lb::MacroscopicPropertyCache& propertyCache = latticeBoltzmannModel->GetPropertyCache();

	propertyCache.ResetRequirements();

	if (monitoringConfig->doIncompressibilityCheck) {
		propertyCache.densityCache.SetRefreshFlag();
		propertyCache.velocityCache.SetRefreshFlag();
	}

	// If extracting property results, check what's required by them.
	if (propertyExtractor != NULL) {
		propertyExtractor->SetRequiredProperties(propertyCache);
	}
}

/**
 * Called on error to abort the simulation and pull-down the MPI environment.
 */
void SimulationMaster::Abort() {
	// This gives us something to work from when we have an error - we get the rank
	// that calls abort, and we get a stack-trace from the exception having been thrown.
	hemelb::log::Logger::Log<hemelb::log::Critical, hemelb::log::Singleton>("ABORTING");
	hemelb::net::MpiEnvironment::Abort(1);

	exit(1);
}

void SimulationMaster::LogStabilityReport() {
	if (monitoringConfig->doIncompressibilityCheck
			&& incompressibilityChecker->AreDensitiesAvailable()) {
		hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("time step %07i :: tau: %.6f, max_relative_press_diff: %.3f, Ma: %.3f, max_vel_phys: %e",
				simulationState->GetTimeStep(),
				latticeBoltzmannModel->GetLbmParams()->GetTau(),
				incompressibilityChecker->GetMaxRelativeDensityDifference(),
				incompressibilityChecker->GetGlobalLargestVelocityMagnitude()
				/ hemelb::Cs,
				unitConverter->ConvertVelocityToPhysicalUnits(incompressibilityChecker->GetGlobalLargestVelocityMagnitude()));
	}

	if (simulationState->GetStability() == hemelb::lb::StableAndConverged) {
		hemelb::log::Logger::Log<hemelb::log::Info, hemelb::log::Singleton>("time step %07i :: steady flow simulation converged",
				simulationState->GetTimeStep());
	}
}


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

			//JM Preliminary houskeeping and datastructures

			// Get number of beds
			nBeds = couplingMap.size()/3;

			// Get number of inlets/outlest for each Bed
			nOutlets.resize(nBeds);  //Lattice = 0
			nInlets.resize(nBeds); //Lattice = 1
			
			inletVelocity.resize(nBeds);
			outletPressure.resize(nBeds);
		
			Qbed.resize(nBeds);

			UdA.resize(nBeds);
			UdV.resize(nBeds);

			UdAold.resize(nBeds);
			UdVold.resize(nBeds);
			
			int cumInlets = 0;
			int cumOutlets = 0;

			for (int i = 0; i < nBeds; ++i)
			{	
				nOutlets[i] = (couplingMap[3*i].size() - 2)/3;
				cumOutlets += (couplingMap[3*i].size() - 2)/3;
				nInlets[i] = (couplingMap[3*i+1].size() - 2)/3;
				cumInlets += (couplingMap[3*i+1].size() - 2)/3;

				inletVelocity[i].resize(nInlets[i]);
				outletPressure[i].resize(nOutlets[i]);
			}

			inletBedOrder.resize(cumInlets);
			inletOrder.resize(cumInlets);
			inletTotalOrder.resize(cumInlets);
			outletBedOrder.resize(cumOutlets);
			outletOrder.resize(cumOutlets);
			outletTotalOrder.resize(cumOutlets);

			int k=0, kk=0;

			for (int i = 0; i < nBeds; ++i)
			{
				for (int j=0; j<nOutlets[i]; ++j)
				{
					//outletBedOrder[k] = i;
					//outletOrder[k] = j;
					outletBedOrder[stoi(couplingMap[3*i][3*j+2])] = stoi(couplingMap[3*i][0]);
					outletOrder[stoi(couplingMap[3*i][3*j+2])] = j;
					//outletTotalOrder[k] = stoi(couplingMap[3*i][3*j+2]);
					k+=1;
					
				}
				
				for (int j=0; j<nInlets[i]; ++j)
				{
					//inletBedOrder[kk] = i;
					//inletOrder[kk] = j;
					inletBedOrder[stoi(couplingMap[3*i+1][3*j+2])] = stoi(couplingMap[3*i+1][0]);
					inletOrder[stoi(couplingMap[3*i+1][3*j+2])] = j;
					//inletTotalOrder[kk] = stoi(couplingMap[3*i+1][3*j+2]);
					kk+=1;
				}
			}
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

		// Get coupledFieldsHere from current HemePure instance.
		std::unordered_map<int, std::vector<double> > coupledFieldsHere;
		// Get coupledFieldsThere from other HemePure instance.
		std::unordered_map<int, std::vector<double> > coupledFieldsThere;
		bool coupling = propertyExtractor->GetCoupledFields(coupledFieldsHere,coupledFieldsThere);

		if (coupling)
		{
			double tmpVel, tmpVel1, tmpVel2;
			int iolet;
			int k=0;
			int kk=0;
			
			if (IsCurrentProcTheIOProc() && (simulationState->GetTimeStep()%simConfig->GetPropertyOutput(0)->frequency) == 0)
			{
				std::cout <<std::endl;
				for (int i=0; i<nBeds; ++i)
				{
					//std::cout << "Step " << simulationState->GetTimeStep() << " Mesh " << simConfig->GetLatticeId() << " Bed " << i << ": UdA " << UdA[i] << ", UdV " << UdV[i] << std::endl;
					std::cout << "Step " << simulationState->GetTimeStep() << " Mesh " << simConfig->GetLatticeId() << " Bed " << i << ": UdA " << UdA[i] << ", UdAold " << UdAold[i] << ", UdV " << UdV[i] << ", UdVold " << UdVold[i] << std::endl;
					
				/*	if (simConfig->GetLatticeId()==0)
					{
						for (int j = 0; j < nInlets[i]; ++j)
						{	
							// Global IOlet ID.
							iolet = stoi(couplingMap[3*i+1][3*j+2]);
						
							std::cout << "Mesh 0 Step " << simulationState->GetTimeStep() << ", inlet " << iolet << " " << coupledFieldsThere.at(iolet)[0] << std::endl;		
						}
					}
					else
					{
						for (int j = 0; j < nInlets[i]; ++j)
						{	
							// Global IOlet ID.
							iolet = stoi(couplingMap[3*i+1][3*j+2]);
						
							std::cout << "Mesh 1 Step " << simulationState->GetTimeStep() << ", inlet " << iolet << " " << coupledFieldsHere.at(iolet)[0] << std::endl;		
						}
					}*/ 

				}
			}
		
			if (simConfig->GetLatticeId()==0 && (simulationState->GetTimeStep()/simConfig->GetPropertyOutput(0)->frequency)%2 == 1) //JM to oscillate velocity/pressure transfer - ASSUMES coupling is first property
			{
				// Pressure transfer from mesh 1 to mesh 0
				if (IsCurrentProcTheIOProc())
				{	
					//std::cout << "coupling housekeeping at Step " << simulationState->GetTimeStep() << " on mesh 0" << std::endl;
					// ONLY FOR GOOD FOR 2 MESHES!
					for (int i = 0; i < nBeds; ++i)
					{
						//Calc Udv and Uda
						//Calc d(Uda)/dt
						//get qA_i for each inlet
						tmpVel1 = 0;
						tmpVel2 = 0;
				
						for (int j = 0; j < nOutlets[i]; ++j)
						{
							// Global IOlet ID.
							iolet = stoi(couplingMap[3*i][3*j+2]);

							// Voutlet = Vinlet/scaleFactor
							tmpVel1 += stod(couplingMap[3*i][3*j+3])*coupledFieldsHere.at(iolet)[0]*coupledFieldsHere.at(iolet)[0];
							tmpVel2 += stod(couplingMap[3*i][3*j+3])*coupledFieldsHere.at(iolet)[0];
							//std::cout << "Mesh 0 Step " << simulationState->GetTimeStep() << ", outlet " << iolet << " " << coupledFieldsHere.at(iolet)[0] << std::endl;		
						}
						UdAold[i] = UdA[i];
					//	UdA[i] = tmpVel1/tmpVel2;
						if (tmpVel2 == 0)
						{
							UdA[i] = 0;
						}
					       	else
						{
							UdA[i] = tmpVel1/tmpVel2;
						}
					
						tmpVel1 = 0;
						tmpVel2 = 0;

						for (int j = 0; j < nInlets[i]; ++j)
						{	
							// Global IOlet ID.
							iolet = stoi(couplingMap[3*i+1][3*j+2]);

							// Voutlet = Vinlet/scaleFactor
							tmpVel1 += stod(couplingMap[3*i+1][3*j+3])*coupledFieldsThere.at(iolet)[0]*coupledFieldsThere.at(iolet)[0];
							tmpVel2 += stod(couplingMap[3*i+1][3*j+3])*coupledFieldsThere.at(iolet)[0];
							
							//std::cout << "Mesh 0 Step " << simulationState->GetTimeStep() << ", inlet " << iolet << " " << coupledFieldsThere.at(iolet)[0] << std::endl;		
						}
						UdVold[i] = UdV[i];
						//UdV[i] = tmpVel1/tmpVel2;
						if (tmpVel2 == 0)
						{
							UdV[i] = 0;
						}
					       	else
						{
							UdV[i] = tmpVel1/tmpVel2;
						}
						
				
					}
			
					//std::cout << "Mesh housekeeping at Step " << simulationState->GetTimeStep() << " on " << simConfig->GetLatticeId() << std::endl;
					for (int i = 0; i < nBeds; ++i)
					{
						if (UdAold[i] < 1e-10 or UdVold[i] < 1e-10)
						{
							std::cout << "0: adjust for gradients" << std::endl;
							UdAold[i] = UdA[i];
							UdVold[i] = UdV[i];
						//	Qbed[i] = UdV[i]*(stod(couplingMap[3*i+1][3])/stod(couplingMap[3*i+1][4]));	
						}
				
						//PhysUnits
						
						//0 tmpVel1 = UdV[i]*(stod(couplingMap[3*i+1][3])/stod(couplingMap[3*i+1][4])) + (stod(couplingMap[3*i+2][2]) + stod(couplingMap[3*i+2][1]))/(4.0*stod(couplingMap[3*i+2][3])) *((UdV[i] - UdVold[i])*(stod(couplingMap[3*i+1][3])/stod(couplingMap[3*i+1][4])) + (UdA[i] - UdAold[i])*(stod(couplingMap[3*i][3])/stod(couplingMap[3*i][4])))/(2.0*unitConverter->ConvertTimeToPhysicalUnits(simConfig->GetPropertyOutput(0)->frequency))  ; //Qbed
						
						//1 tmpVel1 = UdV[i]*(stod(couplingMap[3*i+1][3])/stod(couplingMap[3*i+1][4]));
						//1 tmpVel2 = (stod(couplingMap[3*i+2][2]) + stod(couplingMap[3*i+2][1]))/(8.0*stod(couplingMap[3*i+2][3])) *((UdV[i] - UdVold[i])*(stod(couplingMap[3*i+1][3])/stod(couplingMap[3*i+1][4])) + (UdA[i] - UdAold[i])*(stod(couplingMap[3*i][3])/stod(couplingMap[3*i][4])))/(2.0*unitConverter->ConvertTimeToPhysicalUnits(simConfig->GetPropertyOutput(0)->frequency))  ; //Qbed;

						//tmpVel1 = std::min(UdV[i]*(stod(couplingMap[3*i+1][3])/stod(couplingMap[3*i+1][4])), UdA[i]*(stod(couplingMap[3*i][3])/stod(couplingMap[3*i][4])));
						tmpVel1 = 0.5*(UdV[i]*(stod(couplingMap[3*i+1][3])/stod(couplingMap[3*i+1][4])) + UdA[i]*(stod(couplingMap[3*i][3])/stod(couplingMap[3*i][4])));
						tmpVel2 = (stod(couplingMap[3*i+2][2]) + stod(couplingMap[3*i+2][1]))/(4.0*stod(couplingMap[3*i+2][3])) *(tmpVel1 - Qbed[i])/(2.0*unitConverter->ConvertTimeToPhysicalUnits(simConfig->GetPropertyOutput(0)->frequency))  ; //Qbed;



/* Original versions - pressure based, less LBM applicable.
						//tmpVel1 = stod(couplingMap[3*i+2][3])*0.5*hemelb::BLOOD_DENSITY_Kg_per_m3*(UdA[i]*UdA[i] - UdV[i]*UdV[i]);
					
						tmpVel1 = UdV[i]*(stod(couplingMap[3*i+1][3])/stod(couplingMap[3*i+1][4])) + stod(couplingMap[3*i+2][2]) * 0.5*hemelb::BLOOD_DENSITY_Kg_per_m3*(UdV[i]*UdV[i] - UdVold[i]*UdVold[i])/(2.0*unitConverter->ConvertTimeToPhysicalUnits(simConfig->GetPropertyOutput(0)->frequency))  ; //Qbed
						

						if (tmpVel1<0)
						{
							tmpVel1=0;
							std::cout << "0: Zeroing bed flow"<<std::endl;
						}
							
						tmpVel2 =  stod(couplingMap[3*i+2][1]) * (0.5*hemelb::BLOOD_DENSITY_Kg_per_m3*(UdA[i]*UdA[i] - UdAold[i]*UdAold[i])/(2.0*unitConverter->ConvertTimeToPhysicalUnits(simConfig->GetPropertyOutput(0)->frequency)));//QV_pt2
						//tmpVel2 =  stod(couplingMap[3*i+2][1]) * (0.5*hemelb::BLOOD_DENSITY_Kg_per_m3*(UdV[i]*UdV[i] - UdVold[i]*UdVold[i])/(2.0*unitConverter->ConvertTimeToPhysicalUnits(simConfig->GetPropertyOutput(0)->frequency)) + (tmpVel1 - Qbed[i])/(stod(couplingMap[3*i+2][3]) * 2.0*unitConverter->ConvertTimeToPhysicalUnits(simConfig->GetPropertyOutput(0)->frequency)) );//QV_pt2
						
						//tmpVel2=0;
						if (tmpVel2<0)
						{
							//tmpVel2=0; //To zero the below
							std::cout << "0: Zeroing bed flow pt 2"<<std::endl;
						}
						//std::cout << "Mesh 0 Bed = " << i << ", Qbed = " << tmpVel1 << ", Qbed old = " << Qbed[i] << std::endl;		*/
						Qbed[i] = tmpVel1;

						for (int j = 0; j < nOutlets[i]; ++j)
						{ 
							// Global IOlet ID.
							iolet = stoi(couplingMap[3*i][3*j+2]);
					   		//std::cout << "Outlet: " << iolet << std::endl;		
					   	
							// Original - error: tmpVel = (tmpVel1 + tmpVel2)/(stod(couplingMap[3*i][3*j+3])*stod(couplingMap[3*i][3*j+4])); 
							tmpVel = (tmpVel1 - tmpVel2)*stod(couplingMap[3*i][3*j+4])/stod(couplingMap[3*i][3*j+3]); 	

							if ((tmpVel1-tmpVel2)<0)
							{
								tmpVel=0; //To zero the below
								std::cout << "0: Zeroing bed flow pt 3"<<std::endl;
							}
							outletPressure[i][j] = 0.5*hemelb::BLOOD_DENSITY_Kg_per_m3*tmpVel*tmpVel;
							//outletPressure[i][j] = 0.5*tmpVel1/stod(couplingMap[3*i+2][3]);
							
							std::cout << "0: Storing: Bed " << i << ", let " << iolet << ", pressure " << outletPressure[i][j] << "(vel = " << tmpVel << ") from Qbed " << tmpVel1 << " and Qa " << tmpVel1-tmpVel2 << std::endl;
						}
					
						
					}
					
					for (std::unordered_map<int, std::vector<int> >::iterator it1 = iolet2RankMap.begin();
							it1 != iolet2RankMap.end(); ++it1)
					{
						for (std::vector<int>::iterator it2 = it1->second.begin();
								it2 != it1->second.end(); ++it2)
						{
							// Send IOlet index first.
							ioComms.Send(it1->first, *it2);					
							//ioComms.Send(outletTotalOrder[it1->first], *it2);					
							// Send pressure for this IOlet.
							ioComms.Send(outletPressure[outletBedOrder[it1->first]][outletOrder[it1->first]], *it2);
							
							//std::cout << "Sending to itfirst" << it1->first << ": Bed " << outletBedOrder[it1->first] << ", let " << outletOrder[it1->first] << ", pressure " << outletPressure[outletBedOrder[it1->first]][outletOrder[it1->first]] << std::endl;
						}
					}
				}
				else
				{
				
					double tmp_pressure; int tmp_iolet;
					for (std::vector<int>::iterator it = rank2IoletMap.at(-1).begin();
							it != rank2IoletMap.at(-1).end(); ++it)
					{
						// Receive IOlet index first.
						ioComms.Receive(tmp_iolet,ioComms.GetIORank());
						// Receive pressure for this IOlet.
						ioComms.Receive(tmp_pressure,ioComms.GetIORank());
						//std::cout << "Received iolet " << tmp_iolet << ", pressure " << tmp_pressure << std::endl;

						((hemelb::lb::iolets::InOutLetCosine*)outletValues->GetLocalIolet(tmp_iolet))->SetForce(outletValues->GetLocalIolet(tmp_iolet)->GetNormal()*unitConverter->ConvertPressureDifferenceToLatticeUnits(tmp_pressure)*unitConverter->GetVoxelSize()); //Dynamic pressure; tmp_pressure*area_dot_normal*dx==forceDensityvector);
						
					}
				}
			}
		
			if (simConfig->GetLatticeId()==1 && (simulationState->GetTimeStep()/simConfig->GetPropertyOutput(0)->frequency)%2 == 0) //JM to oscillate velocity/pressure transfer - ASSUMES coupling is first property
			{
				// Velocity transfer from mesh 0 to mesh 1
				if (IsCurrentProcTheIOProc())
				{	
					//std::cout << "coupling housekeeping at Step " << simulationState->GetTimeStep() << " on mesh 1" << std::endl;
					// ONLY FOR GOOD FOR 2 MESHES!
					for (int i = 0; i < nBeds; ++i)
					{
						//Calc Udv and Uda
						//Calc d(Uda)/dt
						//get qA_i for each inlet
						tmpVel1 = 0;
						tmpVel2 = 0;
				
						for (int j = 0; j < nOutlets[i]; ++j)
						{
							// Global IOlet ID.
							iolet = stoi(couplingMap[3*i][3*j+2]);

							// Calculation for datum velocity
							tmpVel1 += stod(couplingMap[3*i][3*j+3])*coupledFieldsThere.at(iolet)[0]*coupledFieldsThere.at(iolet)[0];
							tmpVel2 += stod(couplingMap[3*i][3*j+3])*coupledFieldsThere.at(iolet)[0];
							
							//std::cout << "Mesh 1 Step " << simulationState->GetTimeStep() << ", outlet " << iolet << " " << coupledFieldsThere.at(iolet)[0] << std::endl;		
						}
						UdAold[i] = UdA[i];
						//UdA[i] = tmpVel1/tmpVel2;
						if (tmpVel2 == 0)
						{
							UdA[i] = 0;
						}
					       	else
						{
							UdA[i] = tmpVel1/tmpVel2;
						}
					
						
						tmpVel1 = 0;
						tmpVel2 = 0;

						for (int j = 0; j < nInlets[i]; ++j)
						{
						//	std::cout << "Mesh 1 Bed " <<  i << ", nInlet " << nInlets[i] << " inlet " << j << std::endl;		
							// Global IOlet ID.
							iolet = stoi(couplingMap[3*i+1][3*j+2]);

							// Calculation for datum velocity
							tmpVel1 += stod(couplingMap[3*i+1][3*j+3])*coupledFieldsHere.at(iolet)[0]*coupledFieldsHere.at(iolet)[0];
							tmpVel2 += stod(couplingMap[3*i+1][3*j+3])*coupledFieldsHere.at(iolet)[0];
							
							//std::cout << "Mesh 1 Step " << simulationState->GetTimeStep() << ", inlet " << iolet << " " << coupledFieldsHere.at(iolet)[0] << std::endl;		

						}
						UdVold[i] = UdV[i];
						//UdV[i] = tmpVel1/tmpVel2;
						if (tmpVel2 == 0)
						{
							UdV[i] = 0;
						}
					       	else
						{
							UdV[i] = tmpVel1/tmpVel2;
						}
					}
			
					//std::cout << "Mesh housekeeping at Step " << simulationState->GetTimeStep() << " on " << simConfig->GetLatticeId() << std::endl;
					
					for (int i = 0; i < nBeds; ++i)
					{
						if (UdAold[i] < 1e-10 or UdVold[i] < 1e-10)
						{
							std::cout << "1: adjust for gradients" << std::endl;
							UdAold[i] = UdA[i];
							UdVold[i] = UdV[i];
						//	Qbed[i] = UdV[i]*(stod(couplingMap[3*i+1][3])/stod(couplingMap[3*i+1][4]));	
						}
						
						
						//0 tmpVel1 = UdA[i]*(stod(couplingMap[3*i][3])/stod(couplingMap[3*i][4])) - (stod(couplingMap[3*i+2][1]) + stod(couplingMap[3*i+2][2]))/(4.0*stod(couplingMap[3*i+2][3])) *((UdV[i] - UdVold[i])*(stod(couplingMap[3*i+1][3])/stod(couplingMap[3*i+1][4])) + (UdA[i] - UdAold[i])*(stod(couplingMap[3*i][3])/stod(couplingMap[3*i][4])))/(2.0*unitConverter->ConvertTimeToPhysicalUnits(simConfig->GetPropertyOutput(0)->frequency))  ; //Qbed
						//1 tmpVel1 = UdA[i]*(stod(couplingMap[3*i][3])/stod(couplingMap[3*i][4]));
						//1 tmpVel2 = (stod(couplingMap[3*i+2][1]) + stod(couplingMap[3*i+2][2]))/(8.0*stod(couplingMap[3*i+2][3])) *((UdV[i] - UdVold[i])*(stod(couplingMap[3*i+1][3])/stod(couplingMap[3*i+1][4])) + (UdA[i] - UdAold[i])*(stod(couplingMap[3*i][3])/stod(couplingMap[3*i][4])))/(2.0*unitConverter->ConvertTimeToPhysicalUnits(simConfig->GetPropertyOutput(0)->frequency))  ; //Qbed

						//tmpVel1 = std::min(UdV[i]*(stod(couplingMap[3*i+1][3])/stod(couplingMap[3*i+1][4])), UdA[i]*(stod(couplingMap[3*i][3])/stod(couplingMap[3*i][4])));
						tmpVel1 = 0.5*(UdV[i]*(stod(couplingMap[3*i+1][3])/stod(couplingMap[3*i+1][4])) + UdA[i]*(stod(couplingMap[3*i][3])/stod(couplingMap[3*i][4])));
						tmpVel2 = (stod(couplingMap[3*i+2][2]) + stod(couplingMap[3*i+2][1]))/(4.0*stod(couplingMap[3*i+2][3])) *(tmpVel1 - Qbed[i])/(2.0*unitConverter->ConvertTimeToPhysicalUnits(simConfig->GetPropertyOutput(0)->frequency))  ; //Qbed;



/* Original approaches on pressure basis - less LBM applicable
 *
						//tmpVel1 = stod(couplingMap[3*i+2][3])*0.5*hemelb::BLOOD_DENSITY_Kg_per_m3*(UdA[i]*UdA[i] - UdV[i]*UdV[i]);
						
						
						tmpVel1 = UdA[i]*(stod(couplingMap[3*i][3])/stod(couplingMap[3*i][4])) - stod(couplingMap[3*i+2][1]) * 0.5*hemelb::BLOOD_DENSITY_Kg_per_m3*(UdA[i]*UdA[i] - UdAold[i]*UdAold[i])/(2.0*unitConverter->ConvertTimeToPhysicalUnits(simConfig->GetPropertyOutput(0)->frequency))  ; //Qbed
						

						if (tmpVel1<0)
						{
							tmpVel1=0;
							std::cout << "1: Zeroing bed flow"<<std::endl;
						}
						
						tmpVel2 =  stod(couplingMap[3*i+2][2]) * (0.5*hemelb::BLOOD_DENSITY_Kg_per_m3*(UdV[i]*UdV[i] - UdVold[i]*UdVold[i])/(2.0*unitConverter->ConvertTimeToPhysicalUnits(simConfig->GetPropertyOutput(0)->frequency)));//QV_pt2
						//tmpVel2 =  stod(couplingMap[3*i+2][2]) * (0.5*hemelb::BLOOD_DENSITY_Kg_per_m3*(UdA[i]*UdA[i] - UdAold[i]*UdAold[i])/(2.0*unitConverter->ConvertTimeToPhysicalUnits(simConfig->GetPropertyOutput(0)->frequency)) - (tmpVel1 - Qbed[i])/(stod(couplingMap[3*i+2][3]) * 2.0*unitConverter->ConvertTimeToPhysicalUnits(simConfig->GetPropertyOutput(0)->frequency)) );//QV_pt2
					

//						std::cout << "Mesh 1 Bed = " << i << ", Qbed = " << tmpVel1 << ", Qbed old = " << Qbed[i] << std::endl;		
						//tmpVel2=0;
						if (tmpVel2<0)
						{
							//tmpVel2=0; //To zero the below.
							std::cout << "1: Zeroing bed flowi pt2"<<std::endl;
						}
*/						
						Qbed[i] = tmpVel1;
						
						for (int j = 0; j < nInlets[i]; ++j)
						{
							// Global IOlet ID.
							iolet = stoi(couplingMap[3*i+1][3*j+2]);
							
					   		//std::cout << "Inlet: " << iolet << " ij " << i << " " << j << "  at line " << 3*i+1 << " element " << 3*j+2 << std::endl;		
					   		
							// original error: inletVelocity[i][j] = (tmpVel1 - tmpVel2)/(stod(couplingMap[3*i+1][3*j+3])*stod(couplingMap[3*i+1][3*j+4]));						
							inletVelocity[i][j] = (tmpVel1 + tmpVel2)*stod(couplingMap[3*i+1][3*j+4])/stod(couplingMap[3*i+1][3*j+3]);							
							if ((tmpVel1+tmpVel2)<0)
							{
								inletVelocity[i][j] = 0; //To zero the below.
								std::cout << "1: Zeroing bed flowi pt3"<<std::endl;
							}									
							std::cout << "1: Storing: Bed " << i << ", iolet " << iolet << ", velocity " << inletVelocity[i][j] << " from Qbed " << tmpVel1 << " and Qv " << tmpVel1 + tmpVel2 << std::endl;
						}

						
					}
				
					for (std::unordered_map<int, std::vector<int> >::iterator it1 = iolet2RankMap.begin();
							it1 != iolet2RankMap.end(); ++it1)
					{
						for (std::vector<int>::iterator it2 = it1->second.begin();
								it2 != it1->second.end(); ++it2)
						{				
							//std::cout << "Sending to itfirst " << it1->first << ": Bed " << inletBedOrder[it1->first] << ", let " << inletOrder[it1->first] << ", velocity " << inletVelocity[inletBedOrder[it1->first]][inletOrder[it1->first]] << std::endl;
							
							// Send IOlet index first.
							ioComms.Send(it1->first, *it2);	
							//ioComms.Send(inletTotalOrder[it1->first], *it2);	
							// Send velocity for this IOlet.
							ioComms.Send(inletVelocity[inletBedOrder[it1->first]][inletOrder[it1->first]], *it2);

						//	std::cout << "Done sending to itfirst " << it1->first << std::endl;
						}
					}
					
				}
				else
				{
					double tmp_velocity; int tmp_iolet;
					for (std::vector<int>::iterator it = rank2IoletMap.at(-1).begin();
							it != rank2IoletMap.at(-1).end(); ++it)
					{
						// Receive IOlet index first.
						ioComms.Receive(tmp_iolet,ioComms.GetIORank());
						// Receive velocity for this IOlet.
						ioComms.Receive(tmp_velocity,ioComms.GetIORank());

						//std::cout << "Received iolet " << tmp_iolet << ", velocity " << tmp_velocity << std::endl;

#ifdef HEMELB_USE_VELOCITY_WEIGHTS_FILE
						//JM File Inlets
						//std::cout << "max speed is " << unitConverter->ConvertSpeedToPhysicalUnits(((hemelb::lb::iolets::InOutLetFileVelocity*)inletValues->GetLocalIolet(tmp_iolet))->GetMaxSpeed()) << std::endl;
						//((hemelb::lb::iolets::InOutLetFileVelocity*)inletValues->GetLocalIolet(tmp_iolet))->SetMaxSpeed(unitConverter->ConvertSpeedToLatticeUnits(2.0*tmp_velocity));
						((hemelb::lb::iolets::InOutLetFileVelocity*)inletValues->GetLocalIolet(tmp_iolet))->SetMaxSpeed(unitConverter->ConvertSpeedToLatticeUnits(tmp_velocity/(((hemelb::lb::iolets::InOutLetFileVelocity*)inletValues->GetLocalIolet(tmp_iolet))->weightAve)));
						//std::cout << "max speed is " << unitConverter->ConvertSpeedToPhysicalUnits(((hemelb::lb::iolets::InOutLetFileVelocity*)inletValues->GetLocalIolet(tmp_iolet))->GetMaxSpeed()) << std::endl;
#endif
#ifndef HEMELB_USE_VELOCITY_WEIGHTS_FILE
//						std::cout << "max speed is " << unitConverter->ConvertSpeedToPhysicalUnits(((hemelb::lb::iolets::InOutLetParabolicVelocity*)inletValues->GetLocalIolet(tmp_iolet))->GetMaxSpeed()) << std::endl;
						//JM Parabolic (circular) inlets
						((hemelb::lb::iolets::InOutLetParabolicVelocity*)inletValues->GetLocalIolet(tmp_iolet))->SetMaxSpeed(unitConverter->ConvertSpeedToLatticeUnits(2.0*tmp_velocity));
						//((hemelb::lb::iolets::InOutLetParabolicVelocity*)inletValues->GetLocalIolet(tmp_iolet))->SetMaxSpeed(unitConverter->ConvertSpeedToLatticeUnits(tmp_velocity));
						//std::cout << "max speed is " << unitConverter->ConvertSpeedToPhysicalUnits(((hemelb::lb::iolets::InOutLetParabolicVelocity*)inletValues->GetLocalIolet(tmp_iolet))->GetMaxSpeed()) << std::endl;
						
#endif 

					}
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

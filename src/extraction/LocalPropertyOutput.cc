
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <cassert>
#include <numeric>

#include "extraction/LocalPropertyOutput.h"
#include "io/formats/formats.h"
#include "io/formats/extraction.h"
#include "io/formats/offset.h"
#include "io/writers/xdr/XdrMemWriter.h"
#include "io/writers/xdr/XdrVectorWriter.h"
#include "net/IOCommunicator.h"
#include "constants.h"
#include "units.h"

namespace hemelb
{
	namespace extraction
	{
	    // Declare recursive helper
	    template <typename... Ts>
	    io::writers::Writer& encode(io::writers::Writer& enc, Ts... args);
	    // Terminating case - one arg
	    template <typename T>
	    io::writers::Writer& encode(io::writers::Writer& enc, T arg) {
	      return enc << arg;
	    }
	    // Recursive case - N + 1 args
	    template <typename T, typename... Ts>
	    io::writers::Writer& encode(io::writers::Writer& enc, T arg, Ts... args) {
	      return encode(enc << arg, args...);
	    }


	    // XDR encode some values and return the result buffer
	    template <typename... Ts>
	    std::vector<char> quick_encode(Ts... args) {
	      io::writers::xdr::XdrVectorWriter encoder;
	      encode(encoder, args...);
	      auto ans = encoder.GetBuf();
	      return ans;
	    }
	
    	    LocalPropertyOutput::LocalPropertyOutput(IterableDataSource& dataSource,
				const PropertyOutputFile* outputSpec,
				const net::IOCommunicator& ioComms,
				const MPI_Comm& unComms) :
			comms(ioComms), commsUni(unComms), dataSource(dataSource), outputSpec(outputSpec)
		{
			// Open the file as write-only, create it if it doesn't exist, don't create if the file
			// already exists.
			outputFile = net::MpiFile::Open(comms, outputSpec->filename,
					MPI_MODE_WRONLY | MPI_MODE_CREATE | MPI_MODE_EXCL);

			MPI_Comm_size(commsUni, &sizeUni);
			MPI_Comm_rank(commsUni, &rankUni);

			//// Extract filename to check if we're coupling.
			//std::string filePath = outputSpec->filename;
			//int dotPos = filePath.rfind('.');
			//int sepPos = filePath.rfind('/');
			//std::string filename = filePath.substr(sepPos + 1, filePath.size() - dotPos + 1);

			
			// Count sites on this task.
			uint64_t siteCount = 0;
			dataSource.Reset();
			while (dataSource.ReadNext())
			{
				if (outputSpec->geometry->Include(dataSource, dataSource.GetPosition()))
				{
					++siteCount;
				}
			}

			// Calculate how long local writes need to be.

			// First get the length per-site.
			// Always have 3 uint32's for the position of a site.
			writeLength = 3 * 4;

			// Then add each field's length.
			for (unsigned outputNumber = 0; outputNumber < outputSpec->fields.size(); ++outputNumber)
			{
				writeLength += sizeof(WrittenDataType)
					* GetFieldLength(outputSpec->fields[outputNumber].type);
			}

			// Check if this rank is to do any coupling work.
			for (unsigned outputNumber = 0; outputNumber < outputSpec->fields.size(); ++outputNumber)
			{
				if (outputSpec->fields[outputNumber].type == OutputField::Coupled)
				{
					dataSource.Reset();
					while (dataSource.ReadNext())
					{
						if (outputSpec->geometry->Include(dataSource, dataSource.GetPosition()))
						{
							// Ensure that we're coupling at IOlets!
							if (dataSource.GetIoletId() == -1)
							{
								hemelb::log::Logger::Log<hemelb::log::Critical, hemelb::log::OnePerCore>(
										"Attempt to derive coupling field from site data that is not at IOlet!");
								hemelb::net::MpiEnvironment::Abort(1);
								exit(1);
							}

							// Prepare for coupling.
							if (coupledFields[dataSource.GetIoletId()].empty())
							{
								// Flag to indicate that this rank will take part in coupling.
								coupling = true;
								// Create a list of local IOlets.
								myLets.insert(dataSource.GetIoletId());
								// Prepare for summation.
								// Three values per IOlet.
								// Averaged 1) velocity, 2) pressure, 3) local site count at IOlet.
								coupledFields[dataSource.GetIoletId()] = {0.0, 0.0, 0.0};
							}
							coupledFields[dataSource.GetIoletId()][2]++;
						}
					}
				}
			}

			// Let IO rank know how many coupled IOlets I own.
			if (comms.OnIORank())
			{
				nCoupledLetsWorld.resize(comms.Size());
				for (proc_t rank = 0; rank < comms.Size(); rank++)
				{
					if (comms.GetIORank() != rank)
					{
						comms.Receive(nCoupledLets, rank);
						nCoupledLetsWorld[rank] = nCoupledLets;
					}
				}
			}
			else
			{
				nCoupledLets = myLets.size();
				comms.Send(nCoupledLets, comms.GetIORank());
			}

			if (comms.OnIORank())
			{
				int sum = std::accumulate(
						nCoupledLetsWorld.begin(),
						nCoupledLetsWorld.end(), 0);
				if (sum > 0) coupling = true;
			}

			// We use standard MPI calls.
			if (coupling)
			{
				// Let IO rank know the IDs of coupled IOlets I own.
				if (comms.OnIORank())
				{
					for (proc_t rank = 0; rank < comms.Size(); rank++)
					{
						if (comms.GetIORank() != rank && nCoupledLetsWorld[rank] > 0)
						{
							// Receive coupled IOlet IDs.
							std::vector<int> ioletIds_tmp;
							ioletIds_tmp.resize(nCoupledLetsWorld[rank]);
							MPI_Recv(&ioletIds_tmp[0],
									ioletIds_tmp.size(),
									MPI_INT, rank, 0, comms, MPI_STATUS_IGNORE);

							// Create IOlet map.
							letMap[rank] = ioletIds_tmp;
						}
					}
				}
				else
				{
					// If I own sites at coupled IOlets.
					if (nCoupledLets > 0)
					{
						// Convert unordered set to vector for transmission.
						std::vector<int> ioletIds_tmp;
						for (std::unordered_set<int>::iterator it = myLets.begin();
								it != myLets.end(); ++it)
						{
							ioletIds_tmp.push_back(*it);
						}

						// Send coupled IOlet IDs as a vector to IO rank.
						MPI_Send(&ioletIds_tmp[0],
								ioletIds_tmp.size(),
								MPI_INT, comms.GetIORank(), 0, comms);

						// Put ioletIds_tmp vector into letMap.
						// Needed for retrieval in SimulationMaster.cc.
						// Use key = -1 to indicate incomplete map.
						letMap[-1] = ioletIds_tmp;
					}
				}
			}

			if (comms.OnIORank())
			{
				double partialField;
				for (proc_t rank = 0; rank < comms.Size(); rank++)
				{
					if (comms.GetIORank() != rank && nCoupledLetsWorld[rank] > 0)
					{
						for (std::vector<int>::iterator it = letMap[rank].begin();
								it != letMap[rank].end(); ++it)
						{
							// If this entry doesn't exist, create it.
							if (coupledFields.find(*it) == coupledFields.end())
							{
								coupledFields[*it] = {0.0, 0.0, 0.0};
							}

							comms.Receive(partialField, rank);
							coupledFields.at(*it)[2] += partialField;
						}
					}
				}
			}
			else
			{
				if (nCoupledLets > 0)
				{
					for (std::unordered_set<int>::iterator it = myLets.begin();
							it != myLets.end(); ++it)
					{
						comms.Send(coupledFields[*it][2], comms.GetIORank());
					}
				}
			}

			if (comms.OnIORank())
			{
				coupledField.resize(coupledFields.size());
			}

			//if (comms.OnIORank())
			//{
			//	int i = 0;
			//	char hostname[256];
			//	gethostname(hostname, sizeof(hostname));
			//	printf("PID %d (%d) on %s ready for attach\n", getpid(), 1, hostname);
			//	fflush(stdout);
			//	while (0 == i)
			//		sleep(5);
			//}

			//  Now multiply by local site count.
			writeLength *= siteCount;

			// The IO rank also writes the iteration number.
			if (comms.OnIORank())
			{
				writeLength += 8;
			}

			//! @TODO: These two MPI calls can be replaced with one

			// Everyone needs to know the total length written during one iteration.
			allCoresWriteLength = comms.AllReduce(writeLength, MPI_SUM);

			// Only the root process needs to know the total number of sites written.
			// Note this has a garbage value on other procs.
			allSiteCount = comms.Reduce(siteCount, MPI_SUM,
					comms.GetIORank());

			unsigned totalHeaderLength = 0;

			// Write the header information on the IO rank.
			if (comms.OnIORank())
			{
				// Compute the length of the field header.
				unsigned fieldHeaderLength = 0;
				for (unsigned outputNumber = 0; outputNumber < outputSpec->fields.size(); ++outputNumber)
				{
					// Name.
					fieldHeaderLength
						+= io::formats::extraction::GetStoredLengthOfString(outputSpec->fields[outputNumber].name);
					// uint32 for number of fields.
					fieldHeaderLength += 4;
					// Double for the offset in each field.
					fieldHeaderLength += 8;
				}

				// Create a header buffer
				totalHeaderLength = io::formats::extraction::MainHeaderLength + fieldHeaderLength;
				std::vector<char> headerBuffer(totalHeaderLength);

				{
					// Encoder for ONLY the main header (note shorter length).
					io::writers::xdr::XdrMemWriter
						mainHeaderWriter(&headerBuffer[0], io::formats::extraction::MainHeaderLength);

					// Fill it.
					mainHeaderWriter << uint32_t(io::formats::HemeLbMagicNumber)
						<< uint32_t(io::formats::extraction::MagicNumber)
						<< uint32_t(io::formats::extraction::VersionNumber);
					mainHeaderWriter << double(dataSource.GetVoxelSize());
					const util::Vector3D<distribn_t> &origin = dataSource.GetOrigin();
					mainHeaderWriter << double(origin[0]) << double(origin[1]) << double(origin[2]);

					// Write the total site count and number of fields.
					mainHeaderWriter << uint64_t(allSiteCount) << uint32_t(outputSpec->fields.size())
						<< uint32_t(fieldHeaderLength);
					// Main header now finished.
					// Exiting the block kills the mainHeaderWriter.
				}

				{
					// Create the field header writer.
					io::writers::xdr::XdrMemWriter
						fieldHeaderWriter(&headerBuffer[io::formats::extraction::MainHeaderLength],
								fieldHeaderLength);

					// Write it.
					for (unsigned outputNumber = 0; outputNumber < outputSpec->fields.size(); ++outputNumber)
					{
						fieldHeaderWriter << outputSpec->fields[outputNumber].name
							<< uint32_t(GetFieldLength(outputSpec->fields[outputNumber].type))
							<< GetOffset(outputSpec->fields[outputNumber].type);
					}
					// Exiting the block cleans up the writer.
				}

				// Write from the buffer.
				outputFile.WriteAt(0, headerBuffer);
			}

			// Calculate where each rank should start writing.
			if (comms.OnIORank())
			{
				// For rank 0 this is easy: it passes the value for rank 1 to the rank.
				localDataOffsetIntoFile = totalHeaderLength;

				if (comms.Size() > 1)
				{
					comms.Send(localDataOffsetIntoFile+writeLength, 1, 1);
				}
			}
			else
			{
				// Receive the writing start position from the previous rank.
				comms.Receive(localDataOffsetIntoFile, comms.Rank()-1, 1);

				// Send the next rank its start position.
				if (comms.Rank() != (comms.Size() - 1))
				{
					comms.Send(localDataOffsetIntoFile+writeLength, comms.Rank() + 1, 1);
				}
			}

			// Create the buffer that we'll write each iteration's data into.
			buffer.resize(writeLength);
		

			for (unsigned outputNumber = 0; outputNumber < outputSpec->fields.size(); ++outputNumber)
			{
				if (outputSpec->fields[outputNumber].type == OutputField::Distributions){
					// Create a new file name by changing the suffix - outputSpec
					// must have a .-separated extension!
				
					const auto offsetFileName = io::formats::offset::ExtractionToOffset(outputSpec->filename);

					// Now create the file.
					offsetFile = net::MpiFile::Open(comms, offsetFileName,
							MPI_MODE_WRONLY | MPI_MODE_CREATE | MPI_MODE_EXCL);
					WriteOffsetFile();
				}
			}
		}

		LocalPropertyOutput::~LocalPropertyOutput()
		{
		}

		bool LocalPropertyOutput::ShouldWrite(unsigned long timestepNumber) const
		{
			return ((timestepNumber % outputSpec->frequency) == 0);
		}

		bool LocalPropertyOutput::Coupling() const
		{
			return coupling;
		}

		const PropertyOutputFile* LocalPropertyOutput::GetOutputSpec() const
		{
			return outputSpec;
		}

		void LocalPropertyOutput::Write(unsigned long timestepNumber)
		{
			// Don't write if we shouldn't this iteration.
			if (!ShouldWrite(timestepNumber))
			{
				return;
			}

			// Don't write if this rank doesn't do anything.
			if (writeLength <= 0)
			{
				return;
			}

			// Create the buffer.
			io::writers::xdr::XdrMemWriter xdrWriter(&buffer[0], buffer.size());

			// Firstly, the IO rank must write the iteration number.
			if (comms.OnIORank())
			{
				xdrWriter << (uint64_t) timestepNumber;
			}

			dataSource.Reset();

			// Temporary velocity vector.
			PhysicalVelocity Velocity;

			if (coupling)
			{
				for (std::unordered_map<int, std::vector<double> >::iterator it = coupledFields.begin();
						it != coupledFields.end(); ++it)
				{
					coupledFields.at(it->first)[0] = 0.0;
					coupledFields.at(it->first)[1] = 0.0;
				}
			}

			// Check if we're coupling this step.
			bool couplingNow = false;
			for (unsigned outputNumber = 0; outputNumber < outputSpec->fields.size(); ++outputNumber)
			{
				if (outputSpec->fields[outputNumber].type == OutputField::Coupled)
				{
					couplingNow = true;
				}
			}

			if (writeLength > 0)
			{
				while (dataSource.ReadNext())
				{
					const util::Vector3D<site_t>& position = dataSource.GetPosition();
					if (outputSpec->geometry->Include(dataSource, position))
					{
						// Write the position.
						xdrWriter << (uint32_t) position.x << (uint32_t) position.y << (uint32_t) position.z;

						// Write for each field.
						for (unsigned outputNumber = 0; outputNumber < outputSpec->fields.size(); ++outputNumber)
						{
							switch (outputSpec->fields[outputNumber].type)
							{
								case OutputField::Coupled:
									Velocity[0] = dataSource.GetVelocity().x;
									Velocity[1] = dataSource.GetVelocity().y;
									Velocity[2] = dataSource.GetVelocity().z;
									coupledFields.at(dataSource.GetIoletId())[0] += Velocity.GetMagnitude();
									coupledFields.at(dataSource.GetIoletId())[1] += (
											dataSource.GetPressure() - REFERENCE_PRESSURE_mmHg);
									break;
								case OutputField::Pressure:
									xdrWriter << static_cast<WrittenDataType> (
											dataSource.GetPressure() - REFERENCE_PRESSURE_mmHg);
									
									//printf("(x,y,z) = %e, %e, %e; pressure = %e \n", position.x, position.y, position.z, dataSource.GetPressure());
									
									break;
								case OutputField::Velocity:
									xdrWriter
										<< static_cast<WrittenDataType> (dataSource.GetVelocity().x)
										<< static_cast<WrittenDataType> (dataSource.GetVelocity().y)
										<< static_cast<WrittenDataType> (dataSource.GetVelocity().z);
									break;
									//! @TODO: Work out how to handle the different stresses.
								case OutputField::VonMisesStress:
									xdrWriter << static_cast<WrittenDataType> (dataSource.GetVonMisesStress());
									break;
								case OutputField::ShearStress:
									xdrWriter << static_cast<WrittenDataType> (dataSource.GetShearStress());
									break;
								case OutputField::ShearRate:
									xdrWriter << static_cast<WrittenDataType> (dataSource.GetShearRate());
									break;
								case OutputField::StressTensor:
									{
										util::Matrix3D tensor = dataSource.GetStressTensor();
										// Only the upper triangular part of the symmetric tensor is stored. Storage is row-wise.
										xdrWriter << static_cast<WrittenDataType> (tensor[0][0])
											<< static_cast<WrittenDataType> (tensor[0][1])
											<< static_cast<WrittenDataType> (tensor[0][2])
											<< static_cast<WrittenDataType> (tensor[1][1])
											<< static_cast<WrittenDataType> (tensor[1][2])
											<< static_cast<WrittenDataType> (tensor[2][2]);
										break;
									}
								case OutputField::Traction:
									xdrWriter
										<< static_cast<WrittenDataType> (dataSource.GetTraction().x)
										<< static_cast<WrittenDataType> (dataSource.GetTraction().y)
										<< static_cast<WrittenDataType> (dataSource.GetTraction().z);
									break;
								case OutputField::TangentialProjectionTraction:
									xdrWriter
										<< static_cast<WrittenDataType> (dataSource.GetTangentialProjectionTraction().x)
										<< static_cast<WrittenDataType> (dataSource.GetTangentialProjectionTraction().y)
										<< static_cast<WrittenDataType> (dataSource.GetTangentialProjectionTraction().z);
									break;
								case OutputField::Distributions:
   
							      		unsigned numComponents;
									const distribn_t *d_ptr;
									numComponents = dataSource.GetNumVectors();
									d_ptr = dataSource.GetDistribution();
									for (int i = 0; i < numComponents; i++)
									{
							      			xdrWriter << static_cast<WrittenDataType> (*d_ptr);
							      			d_ptr++;
									}
									break;
								
								case OutputField::MpiRank:
									xdrWriter
										<< static_cast<WrittenDataType> (comms.Rank());
									break;
								default:
									// This should never trip. It only occurs when a new OutputField field is added and no
									// implementation is provided for its serialisation.
									assert(false);
							}
						}
					}
				}
			}

			// If I'm to do coupling work this step.
			if (coupling && couplingNow)
			{
				// Clear coupledFieldsThere to prepare for receipt of new fields from other HemePure instance.
				coupledFieldsThere.clear();

				for (unsigned outputNumberCoupled = 0; outputNumberCoupled < 2; ++outputNumberCoupled)
				{
					if (comms.OnIORank())
					{
						if (coupledFields.at(0)[outputNumberCoupled] > 0.0)
						{
							hemelb::log::Logger::Log<hemelb::log::Critical, hemelb::log::OnePerCore>(
									"coupledFields no reset!");
							hemelb::net::MpiEnvironment::Abort(1);
							exit(1);
						}

						double partialField;
						for (proc_t rank = 0; rank < comms.Size(); rank++)
						{
							if (comms.GetIORank() != rank && nCoupledLetsWorld[rank] > 0)
							{
								for (std::vector<int>::iterator it = letMap[rank].begin();
										it != letMap[rank].end(); ++it)
								{
									comms.Receive(partialField, rank);
									coupledFields.at(*it)[outputNumberCoupled] += partialField;
								}
							}
						}
					}
					else
					{
						if (nCoupledLets > 0)
						{
							for (std::unordered_set<int>::iterator it = myLets.begin();
									it != myLets.end(); ++it)
							{
								comms.Send(coupledFields.at(*it)[outputNumberCoupled], comms.GetIORank());
							}
						}
					}

					// Will send to other HemePure instance.
					// We use standard MPI calls.
					// Split coupledFields to make it easier to send field.
					if (comms.OnIORank())
					{
						for (proc_t rank = 0; rank < comms.Size(); rank++)
						{
							if (nCoupledLetsWorld[rank] > 0)
							{
								// Calculate average.
								for (std::vector<int>::iterator it = letMap[rank].begin();
										it != letMap[rank].end(); ++it)
								{
									// Note, indices of coupledField = IOlet indices.
									coupledField[*it] =
										coupledFields.at(*it)[outputNumberCoupled]/
										coupledFields.at(*it)[2];
								}
							}
						}

						// Number of coupled IOlets in my world.
						int nCoupledLetsHere  = coupledField.size();
						// Number of coupled IOlets in other world.
						int nCoupledLetsThere = 0;

						// Exchange counts.
						if (rankUni == 0)
						{
							MPI_Send(&nCoupledLetsHere,
									1,
									MPI_INT, comms.Size(),
									0, commsUni);
							MPI_Recv(&nCoupledLetsThere,
									1,
									MPI_DOUBLE, comms.Size(),
									1, commsUni, MPI_STATUS_IGNORE);
						}
						else
						{
							MPI_Send(&nCoupledLetsHere,
									1,
									MPI_INT, 0,
									1, commsUni);
							MPI_Recv(&nCoupledLetsThere,
									1,
									MPI_INT, 0,
									0, commsUni, MPI_STATUS_IGNORE);
						}

						// Temporary vector to receive a coupled field.
						std::vector<double> coupledFieldThere(nCoupledLetsThere);

						if (rankUni == 0)
						{
							MPI_Send(&coupledField[0],
									coupledField.size(),
									MPI_DOUBLE, comms.Size(), 0, commsUni);
							MPI_Recv(&coupledFieldThere[0],
									coupledFieldThere.size(),
									MPI_DOUBLE, comms.Size(),
									1, commsUni, MPI_STATUS_IGNORE);
						}
						else
						{
							MPI_Send(&coupledField[0],
									coupledField.size(),
									MPI_DOUBLE, 0, 1, commsUni);
							MPI_Recv(&coupledFieldThere[0],
									coupledFieldThere.size(),
									MPI_DOUBLE, 0,
									0, commsUni, MPI_STATUS_IGNORE);
						}

						for (int i = 0; i < nCoupledLetsThere; i++)
						{
							coupledFieldsThere[i].push_back(coupledFieldThere[i]);
						}
					}
					//if (comms.OnIORank())
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

			// Actually do the MPI writing.
			outputFile.WriteAt(localDataOffsetIntoFile, buffer);

			// Set the offset to the right place for writing on the next iteration.
			localDataOffsetIntoFile += allCoresWriteLength;
		}
		
		// Write the offset file.
		    void LocalPropertyOutput::WriteOffsetFile() {
		      namespace fmt = io::formats;

		      // On process 0 only, write the header
		      if (comms.OnIORank()) {
			auto buf = quick_encode(
						uint32_t(fmt::HemeLbMagicNumber),
						uint32_t(fmt::offset::MagicNumber),
						uint32_t(fmt::offset::VersionNumber),
						uint32_t(comms.Size())
						);
			assert(buf.size() == fmt::offset::HeaderLength);
			offsetFile.WriteAt(0, buf);
		      }
		      // Every rank writes its offset
		      uint64_t offsetForOffset = comms.Rank() * sizeof(localDataOffsetIntoFile)
			+ fmt::offset::HeaderLength;
		      offsetFile.WriteAt(offsetForOffset, quick_encode(localDataOffsetIntoFile));

		      // Last process writes total
		      if (comms.Rank() == (comms.Size()-1)) {
			offsetFile.WriteAt(offsetForOffset + sizeof(localDataOffsetIntoFile),
					   quick_encode(localDataOffsetIntoFile + writeLength));
		      }
		    }
		
		unsigned LocalPropertyOutput::GetFieldLength(OutputField::FieldType field)
		{
			switch (field)
			{
				case OutputField::Coupled:
				case OutputField::Pressure:
				case OutputField::VonMisesStress:
				case OutputField::ShearStress:
				case OutputField::ShearRate:
				case OutputField::MpiRank:
					return 1;
				case OutputField::Velocity:
				case OutputField::Traction:
				case OutputField::TangentialProjectionTraction:
					return 3;
				case OutputField::StressTensor:
					return 6; // We only store the upper triangular part of the symmetric tensor.
				case OutputField::Distributions:
					return latticeType::NUMVECTORS;
				default:
					// This should never trip. Only occurs if someone adds a new field and forgets
					// to add to this method.
					assert(false);
					return 0;
			}
		}

		double LocalPropertyOutput::GetOffset(OutputField::FieldType field)
		{
			switch (field)
			{
				case OutputField::Coupled:
					return REFERENCE_PRESSURE_mmHg;
				case OutputField::Pressure:
					return REFERENCE_PRESSURE_mmHg;
				default:
					return 0.;
			}
		}
	}
}

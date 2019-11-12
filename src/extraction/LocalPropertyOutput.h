
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_EXTRACTION_LOCALPROPERTYOUTPUT_H
#define HEMELB_EXTRACTION_LOCALPROPERTYOUTPUT_H

#include <unordered_map>
#include <unordered_set>

#include "extraction/IterableDataSource.h"
#include "extraction/PropertyOutputFile.h"
#include "lb/lattices/Lattices.h"
#include "net/mpi.h"
#include "net/MpiFile.h"

namespace hemelb
{
	namespace net
	{
		class IOCommunicator;
	}
	namespace extraction
	{
		/**
		 * Stores sufficient information to output property information from this core.
		 */
		class LocalPropertyOutput
		{
			public:
				/**
				 * Initialises a LocalPropertyOutput. Required so we can use const reference types.
				 * @param file
				 * @param offset
				 * @return
				 */
				LocalPropertyOutput(IterableDataSource& dataSource, const PropertyOutputFile* outputSpec, const net::IOCommunicator& ioComms, const MPI_Comm& unComms);

				/**
				 * Tidies up the LocalPropertyOutput (close files etc).
				 * @return
				 */
				~LocalPropertyOutput();

				/**
				 * True if this property output should be written on the current iteration.
				 * @return
				 */
				bool ShouldWrite(unsigned long timestepNumber) const;

				/**
				 * True if this rank is involved in coupling.
				 * @return
				 */
				bool Coupling() const;

				/**
				 * Returns the property output file object to be written.
				 * @return
				 */
				const PropertyOutputFile* GetOutputSpec() const;

				/**
				 * Write this core's section of the data file. Only writes if appropriate for the current
				 * iteration number.
				 */
				void Write(unsigned long timestepNumber);

				/**
				 * Write the offset file
				 */
				void WriteOffsetFile();

				/**
				 * Return coupled fields.
				 */
				const std::unordered_map<int, std::vector<double> > GetCoupledFields() const
				{
					return coupledFieldsThere;
				}

				/**
				 * Return IOlet map.
				 */
				const std::unordered_map<int, std::vector<int> > GetIoletMap() const
				{
					return letMap;
				}

			private:
				/**
				 * Returns the number of floats written for the field.
				 * @param field
				 */
				static unsigned GetFieldLength(OutputField::FieldType field);

				/**
				 * Returns the offset to the field, as it should be written to file.
				 * @param field
				 * @return
				 */
				static double GetOffset(OutputField::FieldType field);

				const net::IOCommunicator& comms;
				const MPI_Comm&            commsUni;
				/**
				 * The MPI file to write into.
				 */
				net::MpiFile outputFile;

				/**
				 * The data source to use for file output.
				 */
				IterableDataSource& dataSource;

				/**
				 * PropertyOutputFile spec.
				 */
				const PropertyOutputFile* outputSpec;

				/**
				 * Total number of sites written.
				 */
				uint64_t allSiteCount;

				/**
				 * Where to begin writing into the file.
				 */
				uint64_t localDataOffsetIntoFile;

				/**
				 * The length, in bytes, of the local write.
				 */
				uint64_t writeLength;

				/**
				 * The length, in bytes, of the total write length.
				 */
				uint64_t allCoresWriteLength;

				/**
				 * Buffer to write into before writing to disk.
				 */
				std::vector<char> buffer;

				/**
				 * The MPI file to write the offsets into.
				 */
				net::MpiFile offsetFile;
				
				/**
				 * Type of written values.
				 */
				typedef float WrittenDataType;

				/**
				 * My place in the universe.
				 */
				proc_t rankUni;

				/**
				 * Size of the universe.
				 */
				proc_t sizeUni;

				/**
				 * Flag to indicate that this rank will take part in coupling.
				 */
				bool coupling = false;

				/**
				 * How many (possibly partial) coupled IOlets on this rank?
				 */
				int nCoupledLets;

				/**
				 * How many (possibly partial) coupled IOlets on each rank?
				 * Only meaningful on IO rank.
				 */
				std::vector<int> nCoupledLetsWorld;

				/**
				 * IDs of my coupled IOlets.
				 */
				std::unordered_set<int> myLets;

				/**
				 * Unordered map of IOlet placement.
				 * Will be a complete map on IO rank:
				 * letMap.first  is rank;
				 * letMap.second is vector of IOlet IDs.
				 * Only meaningful on IO rank.
				 */
				std::unordered_map<int, std::vector<int> > letMap;

				/**
				 * Vector containing averaged
				 * 1) velocities, or
				 * 2) pressures.
				 * Only meaningful on IO rank.
				 */
				std::vector<double> coupledField;

				/**
				 * Unordered map containing (ultimately) averaged
				 * 1) velocities and
				 * 2) pressures
				 * on all IOlets. Also contains count of
				 * 3) sites on any coupled IOlets.
				 * Averaging only occurs on IO rank.
				 */
				std::unordered_map<int, std::vector<double> > coupledFields;
				std::unordered_map<int, std::vector<double> > coupledFieldsThere;

				typedef hemelb::lb::lattices:: HEMELB_LATTICE latticeType;
		};
	}
}

#endif /* HEMELB_EXTRACTION_LOCALPROPERTYOUTPUT_H */

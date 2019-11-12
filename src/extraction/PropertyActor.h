
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_EXTRACTION_PROPERTYACTOR_H
#define HEMELB_EXTRACTION_PROPERTYACTOR_H

#include "extraction/PropertyWriter.h"
#include "io/PathManager.h"
#include "lb/MacroscopicPropertyCache.h"
#include "lb/SimulationState.h"
#include "net/IteratedAction.h"

namespace hemelb
{
	namespace extraction
	{
		class PropertyActor : public net::IteratedAction
		{
			public:
				/**
				 * Constructor, gets the class ready for reading.
				 * @param simulationState
				 * @param propertyOutputs
				 * @param dataSource
				 * @return
				 */
				PropertyActor(const lb::SimulationState& simulationState,
						const std::vector<PropertyOutputFile*>& propertyOutputs,
						IterableDataSource& dataSource,
						reporting::Timers& timers,
						const net::IOCommunicator& ioComms,
						const MPI_Comm& unComms);

				~PropertyActor();

				/**
				 * Set which properties will be required this iteration.
				 * @param propertyCache
				 */
				void SetRequiredProperties(lb::MacroscopicPropertyCache& propertyCache);

				/**
				 * Override the iterated actor end of iteration method to perform writing.
				 */
				void EndIteration();

				/**
				 * Get coupled fields.
				 */
				bool GetCoupledFields(std::unordered_map<int, std::vector<double> >& coupledFields) const;

				/**
				 * Get unordered IOlet map.
				 */
				const std::unordered_map<int, std::vector<int> > GetIoletMap() const;

			private:
				const lb::SimulationState& simulationState;
				const net::IOCommunicator& ioComms;
				PropertyWriter* propertyWriter;
				reporting::Timers& timers;
		};
	}
}

#endif /* HEMELB_EXTRACTION_PROPERTYACTOR_H */

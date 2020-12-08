
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "extraction/PropertyWriter.h"

namespace hemelb
{
	namespace extraction
	{
		PropertyWriter::PropertyWriter(IterableDataSource& dataSource,
				const std::vector<PropertyOutputFile*>& propertyOutputs,
				const net::IOCommunicator& ioComms,
				const MPI_Comm& unComms)
		{
			for (unsigned outputNumber = 0; outputNumber < propertyOutputs.size(); ++outputNumber)
			{
				localPropertyOutputs.push_back(new LocalPropertyOutput(dataSource, propertyOutputs[outputNumber], ioComms, unComms));
			}
		}

		PropertyWriter::~PropertyWriter()
		{
			for (unsigned outputNumber = 0; outputNumber < localPropertyOutputs.size(); ++outputNumber)
			{
				delete localPropertyOutputs[outputNumber];
			}
		}

		const std::vector<LocalPropertyOutput*>& PropertyWriter::GetPropertyOutputs() const
		{
			return localPropertyOutputs;
		}

		void PropertyWriter::Write(unsigned long iterationNumber) const
		{
			for (unsigned outputNumber = 0; outputNumber < localPropertyOutputs.size(); ++outputNumber)
			{
				localPropertyOutputs[outputNumber]->Write((uint64_t) iterationNumber);
			}
		}
		
		bool PropertyWriter::GetCoupledFields(
				unsigned long iterationNumber,
				bool masterRank,
				std::unordered_map<int, std::vector<double> >& coupledFields) const
		{
			for (unsigned outputNumber = 0; outputNumber < localPropertyOutputs.size(); ++outputNumber)
			{
				if (
						localPropertyOutputs[outputNumber]->Coupling() &&
						localPropertyOutputs[outputNumber]->ShouldWrite(iterationNumber))
				{
					if (masterRank)
					{
						coupledFields = localPropertyOutputs[outputNumber]->GetCoupledFieldsThere();
					}
					return true;
				}
				else
				{
					return false;
				}
			}
		}
		
		bool PropertyWriter::GetCoupledFields(
				unsigned long iterationNumber,
				bool masterRank,
				std::unordered_map<int, std::vector<double> >& coupledFieldsHere,
				std::unordered_map<int, std::vector<double> >& coupledFieldsThere) const
		{
			for (unsigned outputNumber = 0; outputNumber < localPropertyOutputs.size(); ++outputNumber)
			{
				if (
						localPropertyOutputs[outputNumber]->Coupling() &&
						localPropertyOutputs[outputNumber]->ShouldWrite(iterationNumber))
				{
					if (masterRank)
					{
						coupledFieldsHere = localPropertyOutputs[outputNumber]->GetCoupledFieldsHere();
						coupledFieldsThere = localPropertyOutputs[outputNumber]->GetCoupledFieldsThere();
					}
					return true;
				}
				else
				{
					return false;
				}
			}
		}

		std::unordered_map<int, std::vector<int> > PropertyWriter::GetIoletMap() const
		{
			std::unordered_map<int, std::vector<int> > ioletMap;
			for (unsigned outputNumber = 0; outputNumber < localPropertyOutputs.size(); ++outputNumber)
			{
				if (localPropertyOutputs[outputNumber]->Coupling())
				{
					ioletMap = localPropertyOutputs[outputNumber]->GetIoletMap();
				}
			}
			return ioletMap;
		}
	}
}

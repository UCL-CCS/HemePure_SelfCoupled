
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#ifndef HEMELB_LB_IOLETS_INOUTLETFILEVELOCITY_H
#define HEMELB_LB_IOLETS_INOUTLETFILEVELOCITY_H

#include "lb/iolets/InOutLetVelocity.h"

namespace hemelb
{
  namespace lb
  {
    namespace iolets
    {

      class InOutLetFileVelocity : public InOutLetVelocity
      {
        public:
          InOutLetFileVelocity();

          InOutLet* Clone() const;
          void Reset(SimulationState &state)
          {
            CalculateTable(state.GetTotalTimeSteps(), state.GetTimeStepLength());
          }

          const std::string& GetFilePath()
          {
            return velocityFilePath;
          }
          void SetFilePath(const std::string& path)
          {
            velocityFilePath = path;
          }

          LatticeVelocity GetVelocity(const LatticePosition& x, const LatticeTimeStep t) const;
          
          const LatticeSpeed& GetMaxSpeed() const //JM
          {
		return maxSpeed;
          }
	  void SetMaxSpeed(const LatticeSpeed& v)  //JM
	  {
		maxSpeed = v;
	  }	

          void Initialise(const util::UnitConverter* unitConverter);

          bool useWeightsFromFile;
          
          bool coupledMaxVelocity; //JM

	  //JM from parabolic
       	  void SetForce(const LatticeForceVector& v) //JM
	  {
	      latticeForce = v;
	  }

	  LatticeForceVector& GetForce() //JM
	  {
	      return latticeForce;
	  }

	  LatticeForceVector GetForceOnBoundary() const;
	
	  LatticePressure GetPressureMean() const //JM1
	  {
	      return densityMean * Cs2;
	  }

	  void SetPressureMean(const LatticePressure& pressure) //JM1
	  {
    	      densityMean = pressure / Cs2;
   	  }
       

	  double weightAve;

	private:
          std::string velocityFilePath;
          std::string velocityWeightsFilePath;
          void CalculateTable(LatticeTimeStep totalTimeSteps, PhysicalTime timeStepLength);
          std::vector<LatticeSpeed> velocityTable;
          const util::UnitConverter* units;
          
	  LatticeSpeed maxSpeed; //JM
          
	  std::map<std::vector<int>, double> weights_table;

   	  LatticeForceVector latticeForce; //JM
	  LatticeDensity densityMean;
       

	  // double weightAve;
	  
	  //double calcVTot(std::vector<double> v);

          //std::vector<double> updateV(std::vector<double> v, std::vector<int> xyz, std::map<std::vector<int>, double> weights_table);

      };

    }
  }
}
#endif /* HEMELB_LB_IOLETS_INOUTLETFILEVELOCITY_H */

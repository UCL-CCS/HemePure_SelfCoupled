
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_IOLETS_INOUTLETPARABOLICVELOCITY_H
#define HEMELB_LB_IOLETS_INOUTLETPARABOLICVELOCITY_H
#include "lb/iolets/InOutLetVelocity.h"
//JM1 #include "lb/iolets/InOutLet.h"

namespace hemelb
{
	namespace lb
	{
		namespace iolets
		{
			class InOutLetParabolicVelocity : public InOutLetVelocity
			{
				public:
					InOutLetParabolicVelocity();
					virtual ~InOutLetParabolicVelocity();
					InOutLet* Clone() const;
					LatticeVelocity GetVelocity(const LatticePosition& x, const LatticeTimeStep t) const;

					const LatticeSpeed& GetMaxSpeed() const
					{
						return maxSpeed;
					}
					void SetMaxSpeed(const LatticeSpeed& v)
					{
						maxSpeed = v;
					}

					void SetWarmup(unsigned int warmup)
					{
						warmUpLength = warmup;
					}

					void SetForce(const LatticeForceVector& v) //JM
                                        {

					    //printf("Set force in PV of %e, %e, %e \n", v.x, v.y, v.z);
					    latticeForce = v;
				        }

				        //const LatticeForceVector& GetForce() const //JM
				        LatticeForceVector& GetForce() //JM
				        {
					    // printf("Got force in PV of %e, %e, %e \n", latticeForce.x, latticeForce.y, latticeForce.z);
					    return latticeForce;
                                        }

					LatticeForceVector GetForceOnBoundary() const;
					
					LatticePressure GetPressureMean() const //JM1
					{
					    return densityMean * Cs2;
					}

					void SetPressureMean(const LatticePressure& pressure) //JM1
					{
					            //printf("setting pressure of %e \n", pressure);
						    densityMean = pressure / Cs2;
					}

				 
				protected:
					LatticeSpeed maxSpeed;
					unsigned int warmUpLength;
					LatticeForceVector latticeForce; //JM

					LatticeDensity densityMean;
			};
		}
	}
}
#endif // HEMELB_LB_IOLETS_INOUTLETPARABOLICVELOCITY_H

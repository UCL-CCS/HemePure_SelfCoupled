
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_IOLETS_INOUTLETVELOCITY_H
#define HEMELB_LB_IOLETS_INOUTLETVELOCITY_H
#include "lb/iolets/InOutLet.h"

namespace hemelb
{
	namespace lb
	{
		namespace iolets
		{
			class InOutLetVelocity : public InOutLet
			{
				public:
					InOutLetVelocity();
					virtual ~InOutLetVelocity();
					LatticeDensity GetDensityMin() const;
					LatticeDensity GetDensityMax() const;
					LatticeDensity GetDensity(LatticeTimeStep time_step) const;

					void Reset(SimulationState &state)
					{
						//pass;
					}
					/**
					 * Note that the radius and max speed for these are specified in LATTICE UNITS in the XML file.
					 * This is indeed a horrible hack.
					 * @return
					 */
					const LatticeDistance& GetRadius() const
					{
						return radius;
					}
					void SetRadius(const LatticeDistance& r)
					{
						radius = r;
					}

					virtual LatticeVelocity GetVelocity(
							const LatticePosition& x, const LatticeTimeStep t) const = 0;
					
					//JM1 LatticeDensity GetBoundaryDensity(const int index);

					void SetForce(const LatticeForceVector& v) //JM
				        {
				         	latticeForce = v;
				        }

				        //const LatticeForceVector& GetForce() const //JM
				        LatticeForceVector GetForce() const//JM
				        {
				 	    return latticeForce;
				        }
 
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
					LatticeDistance radius;
				        LatticeForceVector latticeForce; //JM

					LatticeDensity densityMean; //JM1
			};
		}
	}
}
#endif // HEMELB_LB_IOLETS_INOUTLETVELOCITY_H

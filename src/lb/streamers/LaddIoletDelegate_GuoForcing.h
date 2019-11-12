
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_STREAMERS_LADDIOLETDELEGATE_GUOFORCING_H
#define HEMELB_LB_STREAMERS_LADDIOLETDELEGATE_GUOFORCING_H

#include "lb/streamers/SimpleBounceBackDelegate.h"

namespace hemelb
{
  namespace lb
  {
    namespace streamers
    {

      template<typename CollisionImpl>
      class LaddIoletDelegate_GuoForcing : public SimpleBounceBackDelegate<CollisionImpl>
      {
        public:
          typedef CollisionImpl CollisionType;
          typedef typename CollisionType::CKernel::LatticeType LatticeType;

          LaddIoletDelegate_GuoForcing(CollisionType& delegatorCollider, kernels::InitParams& initParams) :
              SimpleBounceBackDelegate<CollisionType>(delegatorCollider, initParams),
                  bValues(initParams.boundaryObject)
          {
          }

          inline void StreamLink(const LbmParameters* lbmParams,
                                 geometry::LatticeData* const latticeData,
                                 const geometry::Site<geometry::LatticeData>& site,
                                 kernels::HydroVars<typename CollisionType::CKernel>& hydroVars,
                                 const Direction& ii)
          {
            // Translating from Ladd, J. Fluid Mech. "Numerical simulations
            // of particulate suspensions via a discretized Boltzmann
            // equation. Part 1. Theoretical foundation", 1994
            // Eq (3.2) -- simple bounce-back -- becomes:
            //   f_i'(r, t+1) = f_i(r, t*)
            // Eq (3.3) --- modified BB -- becomes:
            //   f_i'(r, t+1) = f_i(r, t*) - 2 a1_i \rho u . c_i
            // where u is the velocity of the boundary half way along the
            // link and a1_i = w_1 / cs2

            int boundaryId = site.GetIoletId();
            iolets::InOutLetVelocity* iolet =
                dynamic_cast<iolets::InOutLetVelocity*>(bValues->GetLocalIolet(boundaryId));
            LatticePosition sitePos(site.GetGlobalSiteCoords());

            LatticePosition halfWay(sitePos);
            halfWay.x += 0.5 * LatticeType::CX[ii];
            halfWay.y += 0.5 * LatticeType::CY[ii];
            halfWay.z += 0.5 * LatticeType::CZ[ii];

            LatticeVelocity wallMom(iolet->GetVelocity(halfWay, bValues->GetTimeStep()));
            //TODO: Add site.GetGlobalSiteCoords() as a first argument?

            //if (CollisionType::CKernel::LatticeType::IsLatticeCompressible())
            //{
            //  wallMom *= hydroVars.density;
            //}

#ifdef HEMELB_USE_VELOCITY_WEIGHTS_FILE
            //JM File inlets
	    wallMom *= ((hemelb::lb::iolets::InOutLetFileVelocity*)iolet)->GetPressureMean()/Cs2; //JM1
            hydroVars.density = ((hemelb::lb::iolets::InOutLetFileVelocity*)iolet)->GetPressureMean()/Cs2; //JM1
	    hydroVars.force = ((hemelb::lb::iolets::InOutLetFileVelocity*)(bValues->GetLocalIolet(boundaryId)))->GetForce(); //JM some Guo Forcing terms as a check

#endif
#ifndef HEMELB_USE_VELOCITY_WEIGHTS_FILE
	    //JM Parabolic (circular) inlets
	    wallMom *= ((hemelb::lb::iolets::InOutLetParabolicVelocity*)iolet)->GetPressureMean()/Cs2; //JM1
            hydroVars.density = ((hemelb::lb::iolets::InOutLetParabolicVelocity*)iolet)->GetPressureMean()/Cs2; //JM1
	    hydroVars.force = ((hemelb::lb::iolets::InOutLetParabolicVelocity*)(bValues->GetLocalIolet(boundaryId)))->GetForce(); //JM some Guo Forcing terms as a check
#endif
	     LatticeType::CalculateForceDistribution(lbmParams->GetTau(),
			hydroVars.velocity.x,
			hydroVars.velocity.y,
			hydroVars.velocity.z,
			hydroVars.force.x,
			hydroVars.force.y,
			hydroVars.force.z,
			hydroVars.forceDist.f);



            distribn_t correction = 2. * LatticeType::EQMWEIGHTS[ii]
                * (wallMom.x * LatticeType::CX[ii] + wallMom.y * LatticeType::CY[ii]
                    + wallMom.z * LatticeType::CZ[ii]) / Cs2;

            * (latticeData->GetFNew(SimpleBounceBackDelegate<CollisionImpl>::GetBBIndex(site.GetIndex(),
                                                                                        ii))) =
                //JM hydroVars.GetFPostCollision()[ii] - correction;
                hydroVars.GetFPostCollision()[ii]*hydroVars.density - correction + hydroVars.forceDist.f[ii];
            
          }

        private:
          iolets::BoundaryValues* bValues;
      };

    }
  }
}

#endif /* HEMELB_LB_STREAMERS_LADDIOLETDELEGATE_H */

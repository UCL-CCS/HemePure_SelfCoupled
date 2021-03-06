
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_STREAMERS_NASHZEROTHORDERPRESSUREDELEGATE_GUOFORCING_H
#define HEMELB_LB_STREAMERS_NASHZEROTHORDERPRESSUREDELEGATE_GUOFORCING_H

#include "util/utilityFunctions.h"
#include "lb/streamers/BaseStreamerDelegate.h"

namespace hemelb
{
	namespace lb
	{
		namespace streamers
		{
			template<typename CollisionImpl>
				class NashZerothOrderPressureDelegate_GuoForcing : public BaseStreamerDelegate<CollisionImpl>
			{
				public:
					typedef CollisionImpl CollisionType;
					typedef typename CollisionType::CKernel::LatticeType LatticeType;

					NashZerothOrderPressureDelegate_GuoForcing(CollisionType& delegatorCollider, kernels::InitParams& initParams) :
						collider(delegatorCollider), iolet(*initParams.boundaryObject)
				{
				}

					inline void StreamLink(const LbmParameters* lbmParams,
							geometry::LatticeData* const latticeData,
							const geometry::Site<geometry::LatticeData>& site,
							kernels::HydroVars<typename CollisionType::CKernel>& hydroVars,
							const Direction& direction)
					{
						int boundaryId = site.GetIoletId();

						// Set the density at the "ghost" site to be the density of the iolet.
						distribn_t ghostDensity = iolet.GetBoundaryDensity(boundaryId);

						// Calculate the velocity at the ghost site, as the component normal to the iolet.
						util::Vector3D<float> ioletNormal = iolet.GetLocalIolet(boundaryId)->GetNormal();

						// Note that the division by density compensates for the fact that v_x etc have momentum
						// not velocity.
						distribn_t component = (hydroVars.momentum / hydroVars.density).Dot(ioletNormal);

						// TODO it's ugly that we have to do this.
						// TODO having to give 0 as an argument is also ugly.
						// TODO it's ugly that we have to give hydroVars a nonsense distribution vector
						// that doesn't get used.
						kernels::HydroVars<typename CollisionType::CKernel> ghostHydrovars(site);

						ghostHydrovars.density = ghostDensity;
						ghostHydrovars.momentum = ioletNormal * component * ghostDensity;
						
					//	ghostHydrovars.force =  iolet.GetBoundaryForce(boundaryId); //JM
						ghostHydrovars.force =  ((hemelb::lb::iolets::InOutLetCosine*)iolet.GetLocalIolet(boundaryId))->GetForce();
						//ghostHydrovars.force =  iolet.GetBoundaryForce(boundaryId);
									
						collider.kernel.CalculateFeq(ghostHydrovars, 0);

						Direction unstreamed = LatticeType::INVERSEDIRECTIONS[direction];

						LatticeType::CalculateForceDistribution(lbmParams->GetTau(),
							ghostHydrovars.velocity.x,
							ghostHydrovars.velocity.y,
							ghostHydrovars.velocity.z,
							ghostHydrovars.force.x,
							ghostHydrovars.force.y,
							ghostHydrovars.force.z,
							ghostHydrovars.forceDist.f);

						*latticeData->GetFNew(site.GetIndex() * LatticeType::NUMVECTORS + unstreamed)
							= ghostHydrovars.GetFEq()[unstreamed]
							+ ghostHydrovars.forceDist.f[unstreamed];
					}
				protected:
					CollisionType& collider;
					iolets::BoundaryValues& iolet;
			};
		}
	}
}

#endif // HEMELB_LB_STREAMERS_NASHZEROTHORDERPRESSUREDELEGATE_GUOFORCING_H

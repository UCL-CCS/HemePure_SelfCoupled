# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

add_library(hemelb_geometry
	Block.cc
	BlockTraverser.cc
	BlockTraverserWithVisitedBlockTracker.cc
	GeometryReader.cc
	LatticeData.cc
	SiteData.cc
	SiteDataBare.cc
	SiteTraverser.cc
	VolumeTraverser.cc
	decomposition/BasicDecomposition.cc
	decomposition/OptimisedDecomposition.cc
	needs/Needs.cc
	neighbouring/NeighbouringDataManager.cc
	neighbouring/NeighbouringLatticeData.cc
	neighbouring/RequiredSiteInformation.cc)

configure_file(
	"${PROJECT_SOURCE_DIR}/geometry/decomposition/DecompositionWeights.h.in"
	"${PROJECT_BINARY_DIR}/geometry/decomposition/DecompositionWeights.h")

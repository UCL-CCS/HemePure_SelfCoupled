
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <iostream>
#include <sstream>
#include <iomanip>
#include <cstdio>
#include <cstdarg>
#include <sys/time.h>
#include <sys/resource.h>
#include <string.h>

#include "util/utilityFunctions.h"
#include "net/mpi.h"
#include "log/Logger.h"

namespace hemelb
{
	namespace log
	{
		const LogLevel Logger::currentLogLevel = HEMELB_LOG_LEVEL;
		// Use negative value to indicate uninitialised.
		int Logger::thisRank = -1;
		double Logger::startTime = -1.0;

		void Logger::Init(int MY_WORLD_RANK)
		{	// If Logger uninitialised
			if (thisRank < 0)
			{
				// Check that MPI is ready
				if (net::MpiEnvironment::Initialized())
				{
					thisRank = MY_WORLD_RANK;
				}
				startTime = util::myClock();
			}
		}

		template<>
			void Logger::LogInternal<OnePerCore>(std::string format, std::va_list args)
			{
				std::stringstream output;
				// Set the fill digit to be 0, so the integer 1 renders as 0000001
				output.fill('0');
				output << "[Rank " << std::setw(7) << thisRank << ", " << std::setiosflags(std::ios::scientific)
					<< std::setw(8) << (util::myClock() - startTime) << " s";

#ifdef HAVE_RUSAGE
				rusage usage;
				getrusage(RUSAGE_SELF, &usage);
				output << ", " << std::setw(14) << usage.ru_maxrss;
#endif
				output << " kB] :: " << format << '\n';

				std::string overFormat(output.str());
				std::vprintf(overFormat.c_str(), args);
			}

		template<>
			void Logger::LogInternal<Singleton>(std::string format, std::va_list args)
			{
				if (thisRank == 0)
				{
					std::stringstream output;
					// Set the fill digit to be 0, so the integer 1 renders as 0000001
					output.fill('0');
					output << "[Rank " << std::setw(7) << thisRank << ", " << std::setiosflags(std::ios::scientific)
						<< std::setw(8) << (util::myClock() - startTime) << " s";

#ifdef HAVE_RUSAGE
					rusage usage;
					getrusage(RUSAGE_SELF, &usage);
					output << ", " << std::setw(14) << usage.ru_maxrss;
#endif
					output << " kB] :: " << format << '\n';

					std::string overFormat(output.str());
					std::vprintf(overFormat.c_str(), args);
				}
			}
	}
}

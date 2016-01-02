/*  Copyright 2014 IST Austria

Contributed by: Jan Reininghaus

This file is part of DIPHA.

DIPHA is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

DIPHA is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with DIPHA.  If not, see <http://www.gnu.org/licenses/>. */

#pragma once

#include <dipha/includes.h>

#define DIPHA_MACROS_BENCHMARK(Code) \
if (dipha::globals::benchmark) \
{ \
  MPI_Barrier( MPI_COMM_WORLD ); \
  int64_t cur_mem_before = getCurrentRSS() >> 20; \
  double start = MPI_Wtime(); \
  int64_t bytes_sent_before = dipha::globals::bytes_received; \
  Code \
  MPI_Barrier( MPI_COMM_WORLD ); \
  dipha::mpi_utils::cout_if_root() << std::setw(10) << std::setiosflags(std::ios::fixed) << std::setiosflags(std::ios::showpoint) \
                                    << std::setprecision(1) << MPI_Wtime() - start << "s" << std::setw(10) << cur_mem_before << " MB" \
                                    << std::setw(10) << (getPeakRSS() >> 20) << " MB" \
                                    << std::setw(10) << ((dipha::globals::bytes_received - bytes_sent_before) >> 20) << " MB   " \
                                    << #Code << std::endl; \
} \
else \
{ \
  Code \
}
  

// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

#include "parallel.h"

#ifndef USE_MPI

auto MPI_Barrier(MPI_Comm) -> void {}
auto MPI_Irecv(const double*, int, int, int, int, MPI_Comm, MPI_Request*)
  -> void {}
auto MPI_Isend(const double*, int, int, int, int, MPI_Comm, MPI_Request*)
  -> void {}
auto MPI_Irecv(const float*, int, int, int, int, MPI_Comm, MPI_Request*)
  -> void {}
auto MPI_Isend(const float*, int, int, int, int, MPI_Comm, MPI_Request*)
  -> void {}
auto MPI_Test(MPI_Request*, int*, MPI_Status*) -> void {}

#endif

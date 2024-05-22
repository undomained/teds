// This source code is licensed under the 3-clause BSD license found
// in the LICENSE file in the root directory of this project.

// Provides stub functions and constants for MPI

#pragma once

using MPI_Comm = int;

constexpr MPI_Comm MPI_COMM_NULL {};
constexpr MPI_Comm MPI_COMM_WORLD { MPI_COMM_NULL };
constexpr int MPI_DOUBLE {};
constexpr int MPI_IN_PLACE {};
constexpr int MPI_SUM {};
constexpr int MPI_STATUS_IGNORE {};

auto MPI_Comm_rank(MPI_Comm, int*) -> void;
auto MPI_Comm_size(MPI_Comm, int*) -> void;
auto MPI_Initialized(int*) -> void;
auto MPI_Finalized(int*) -> void;
auto MPI_Allreduce(int, double*, int, int, int, MPI_Comm) -> void;
auto MPI_Comm_split(MPI_Comm, int, int, MPI_Comm*) -> void;
auto MPI_Recv(void*, int, int, int, int, MPI_Comm, int) -> void;
auto MPI_Send(const void*, int, int, int, int, MPI_Comm) -> void;

using MPI_Request = int;
using MPI_Status = int;

constexpr int MPI_REQUEST_NULL {};
constexpr int MPI_INT {};
constexpr float MPI_FLOAT {};

auto MPI_Barrier(MPI_Comm) -> void;
auto MPI_Irecv(const double*, int, int, int, int, MPI_Comm, MPI_Request*)
  -> void;
auto MPI_Isend(const double*, int, int, int, int, MPI_Comm, MPI_Request*)
  -> void;
auto MPI_Irecv(const float*, int, int, int, int, MPI_Comm, MPI_Request*)
  -> void;
auto MPI_Isend(const float*, int, int, int, int, MPI_Comm, MPI_Request*)
  -> void;
auto MPI_Test(MPI_Request*, int*, MPI_Status*) -> void;

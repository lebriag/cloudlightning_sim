#ifndef COMMUNICATOR_H
#define COMMUNICATOR_H
#include <mpi.h>
#include <list>

class cell;
class gs;
class task;

using std::list;

class communicator
{
 public:
  /// Sends configuration parameters from the gateway to the cells
  void simulationParameters(class siminputs& si, const int& rank, const int& clusterSize, const MPI_Comm& Comm);

  /// Retrieves statistics from the cells
  void cellStatistics(const gs* gates, const cell* clCell, const int& rank, const int& clusterSize,
                      const MPI_Comm& Comm);

  /// Sends tasks to the appropriate cells for further processing
  void taskParameters(list<task>& jobs, const int& rank, const int& clusterSize, const int* commCell,
                      const MPI_Comm& Comm);
};
#endif

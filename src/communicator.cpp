#include <cell.h>
#include <communicator.h>
#include <gs.h>
#include <stat.h>
#include <task.h>

void communicator::simulationParameters(struct siminputs& si, const int& rank, const int& clusterSize,
                                        const MPI_Comm& Comm)
{
  int i, j;
  int one = 1;
  MPI_Status status;

  if (rank == 0) {
    for (i = 0; i < clusterSize - 1; i++) {
      // Send simulator-specific info
      MPI_Send(&si.maxTime, 1, MPI_DOUBLE, i + 1, i + 1, Comm);
      MPI_Send(&si.sosmIntegration, 1, MPI_INT, i + 1, i + 1, Comm);
      MPI_Send(&si.updateInterval, 1, MPI_DOUBLE, i + 1, i + 1, Comm);
      MPI_Send(&one, 1, MPI_INT, i + 1, i + 1, Comm);

      // Send cell-specific info
      MPI_Send(&si.cinp[i].ID, 1, MPI_INT, i + 1, i + 1, Comm);
      MPI_Send(&si.cinp[i].numberOfTypes, 1, MPI_INT, i + 1, i + 1, Comm);
      MPI_Send(si.cinp[i].numberOfResourcesPerType, si.cinp[i].numberOfTypes, MPI_INT, i + 1, i + 1, Comm);
      MPI_Send(si.cinp[i].types, si.cinp[i].numberOfTypes, MPI_INT, i + 1, i + 1, Comm);

      MPI_Send(&si.cinp[i].binp[0].numberOfFunctions, 1, MPI_INT, i + 1, i + 1, Comm);
      MPI_Send(si.cinp[i].binp[0].Ws, si.cinp[i].binp[0].numberOfFunctions, MPI_DOUBLE, i + 1, i + 1, Comm);
      MPI_Send(&si.cinp[i].binp[0].initResPervRM, 1, MPI_INT, i + 1, i + 1, Comm);
      MPI_Send(&si.cinp[i].binp[0].initvRMPerpSwitch, 1, MPI_INT, i + 1, i + 1, Comm);
      MPI_Send(&si.cinp[i].binp[0].initpSwitchPerpRouter, 1, MPI_INT, i + 1, i + 1, Comm);
      MPI_Send(&si.cinp[i].binp[0].pollIntervalCellM, 1, MPI_DOUBLE, i + 1, i + 1, Comm);
      MPI_Send(&si.cinp[i].binp[0].pollIntervalpRouter, 1, MPI_DOUBLE, i + 1, i + 1, Comm);
      MPI_Send(&si.cinp[i].binp[0].pollIntervalpSwitch, 1, MPI_DOUBLE, i + 1, i + 1, Comm);
      MPI_Send(&si.cinp[i].binp[0].pollIntervalvRM, 1, MPI_DOUBLE, i + 1, i + 1, Comm);
      MPI_Send(&si.cinp[i].binp[0].vRMdeploystrategy, 1, MPI_INT, i + 1, i + 1, Comm);

      MPI_Send(&si.cinp[i].ninp[0].netBW, 1, MPI_DOUBLE, i + 1, i + 1, Comm);
      MPI_Send(&si.cinp[i].ninp[0].overCommitmentNetwork, 1, MPI_DOUBLE, i + 1, i + 1, Comm);

      // Send resource-specific info
      for (j = 0; j < si.cinp[i].numberOfTypes; j++) {
        MPI_Send(&si.cinp[i].rinp[j].type, 1, MPI_INT, i + 1, i + 1, Comm);
        MPI_Send(&si.cinp[i].rinp[j].numOfProcUnits, 1, MPI_DOUBLE, i + 1, i + 1, Comm);
        MPI_Send(&si.cinp[i].rinp[j].totalMemory, 1, MPI_DOUBLE, i + 1, i + 1, Comm);
        MPI_Send(&si.cinp[i].rinp[j].totalStorage, 1, MPI_DOUBLE, i + 1, i + 1, Comm);
        MPI_Send(&si.cinp[i].rinp[j].overcommitmentProcessors, 1, MPI_DOUBLE, i + 1, i + 1, Comm);
        MPI_Send(&si.cinp[i].rinp[j].overcommitmentMemory, 1, MPI_DOUBLE, i + 1, i + 1, Comm);
        MPI_Send(&si.cinp[i].rinp[j].computeCapability, 1, MPI_DOUBLE, i + 1, i + 1, Comm);
        MPI_Send(&si.cinp[i].rinp[j].accelerator, 1, MPI_INT, i + 1, i + 1, Comm);
        MPI_Send(&si.cinp[i].rinp[j].totalAccelerators, 1, MPI_INT, i + 1, i + 1, Comm);
        MPI_Send(&si.cinp[i].rinp[j].acceleratorComputeCapability, 1, MPI_DOUBLE, i + 1, i + 1, Comm);

        MPI_Send(&si.cinp[i].pinp[j].typeCpu, 1, MPI_INT, i + 1, i + 1, Comm);
        MPI_Send(&si.cinp[i].pinp[j].cpuPmin, 1, MPI_DOUBLE, i + 1, i + 1, Comm);
        MPI_Send(&si.cinp[i].pinp[j].cpuPmax, 1, MPI_DOUBLE, i + 1, i + 1, Comm);
        MPI_Send(&si.cinp[i].pinp[j].numOfPoints, 1, MPI_INT, i + 1, i + 1, Comm);
        if (si.cinp[i].pinp[j].numOfPoints > 0) {
          MPI_Send(si.cinp[i].pinp[j].cpubins, si.cinp[i].pinp[j].numOfPoints, MPI_DOUBLE, i + 1, i + 1, Comm);

          MPI_Send(si.cinp[i].pinp[j].cpuP, si.cinp[i].pinp[j].numOfPoints, MPI_DOUBLE, i + 1, i + 1, Comm);
        }
        MPI_Send(&si.cinp[i].pinp[j].cpuC, 1, MPI_DOUBLE, i + 1, i + 1, Comm);
        MPI_Send(&si.cinp[i].pinp[j].typeAcc, 1, MPI_INT, i + 1, i + 1, Comm);
        MPI_Send(&si.cinp[i].pinp[j].accPmin, 1, MPI_DOUBLE, i + 1, i + 1, Comm);
        MPI_Send(&si.cinp[i].pinp[j].accPmax, 1, MPI_DOUBLE, i + 1, i + 1, Comm);
        MPI_Send(&si.cinp[i].pinp[j].accC, 1, MPI_DOUBLE, i + 1, i + 1, Comm);
      }
    }
  } else {
    // Receive Simulator specific info
    si.alloc = 1;
    MPI_Recv(&si.maxTime, 1, MPI_DOUBLE, 0, rank, Comm, &status);
    MPI_Recv(&si.sosmIntegration, 1, MPI_INT, 0, rank, Comm, &status);
    MPI_Recv(&si.updateInterval, 1, MPI_DOUBLE, 0, rank, Comm, &status);
    MPI_Recv(&si.numOfCells, 1, MPI_INT, 0, rank, Comm, &status);

    // Receive Cell specific info
    si.cinp = new cellinputs[1];
    si.cinp->alloc = 1;
    MPI_Recv(&si.cinp->ID, 1, MPI_INT, 0, rank, Comm, &status);
    MPI_Recv(&si.cinp->numberOfTypes, 1, MPI_INT, 0, rank, Comm, &status);
    si.cinp->numberOfResourcesPerType = new int[si.cinp->numberOfTypes];
    si.cinp->types = new int[si.cinp->numberOfTypes];
    MPI_Recv(si.cinp->numberOfResourcesPerType, si.cinp->numberOfTypes, MPI_INT, 0, rank, Comm, &status);
    MPI_Recv(si.cinp->types, si.cinp->numberOfTypes, MPI_INT, 0, rank, Comm, &status);

    // Receive Broker info
    si.cinp->binp = new brinputs[1];
    si.cinp->binp[0].alloc = 1;
    MPI_Recv(&si.cinp->binp[0].numberOfFunctions, 1, MPI_INT, 0, rank, Comm, &status);
    si.cinp->binp[0].Ws = new double[si.cinp->binp[0].numberOfFunctions];
    MPI_Recv(si.cinp->binp[0].Ws, si.cinp->binp[0].numberOfFunctions, MPI_DOUBLE, 0, rank, Comm, &status);
    MPI_Recv(&si.cinp->binp[0].initResPervRM, 1, MPI_INT, 0, rank, Comm, &status);
    MPI_Recv(&si.cinp->binp[0].initvRMPerpSwitch, 1, MPI_INT, 0, rank, Comm, &status);
    MPI_Recv(&si.cinp->binp[0].initpSwitchPerpRouter, 1, MPI_INT, 0, rank, Comm, &status);
    MPI_Recv(&si.cinp->binp[0].pollIntervalCellM, 1, MPI_DOUBLE, 0, rank, Comm, &status);
    MPI_Recv(&si.cinp->binp[0].pollIntervalpRouter, 1, MPI_DOUBLE, 0, rank, Comm, &status);
    MPI_Recv(&si.cinp->binp[0].pollIntervalpSwitch, 1, MPI_DOUBLE, 0, rank, Comm, &status);
    MPI_Recv(&si.cinp->binp[0].pollIntervalvRM, 1, MPI_DOUBLE, 0, rank, Comm, &status);
    MPI_Recv(&si.cinp->binp[0].vRMdeploystrategy, 1, MPI_INT, 0, rank, Comm, &status);

    // Receive Interconnection info
    si.cinp->ninp = new netinputs[1];
    si.cinp->ninp[0].alloc = 1;
    MPI_Recv(&si.cinp->ninp[0].netBW, 1, MPI_DOUBLE, 0, rank, Comm, &status);
    MPI_Recv(&si.cinp->ninp[0].overCommitmentNetwork, 1, MPI_DOUBLE, 0, rank, Comm, &status);

    // Receive Resource specific info
    si.cinp->rinp = new resinputs[si.cinp->numberOfTypes];
    si.cinp->pinp = new powinputs[si.cinp->numberOfTypes];
    for (j = 0; j < si.cinp->numberOfTypes; j++) {
      si.cinp->rinp[j].alloc = 1;
      MPI_Recv(&si.cinp->rinp[j].type, 1, MPI_INT, 0, rank, Comm, &status);
      MPI_Recv(&si.cinp->rinp[j].numOfProcUnits, 1, MPI_DOUBLE, 0, rank, Comm, &status);
      MPI_Recv(&si.cinp->rinp[j].totalMemory, 1, MPI_DOUBLE, 0, rank, Comm, &status);
      MPI_Recv(&si.cinp->rinp[j].totalStorage, 1, MPI_DOUBLE, 0, rank, Comm, &status);
      MPI_Recv(&si.cinp->rinp[j].overcommitmentProcessors, 1, MPI_DOUBLE, 0, rank, Comm, &status);
      MPI_Recv(&si.cinp->rinp[j].overcommitmentMemory, 1, MPI_DOUBLE, 0, rank, Comm, &status);
      MPI_Recv(&si.cinp->rinp[j].computeCapability, 1, MPI_DOUBLE, 0, rank, Comm, &status);
      MPI_Recv(&si.cinp->rinp[j].accelerator, 1, MPI_INT, 0, rank, Comm, &status);
      MPI_Recv(&si.cinp->rinp[j].totalAccelerators, 1, MPI_INT, 0, rank, Comm, &status);
      MPI_Recv(&si.cinp->rinp[j].acceleratorComputeCapability, 1, MPI_DOUBLE, 0, rank, Comm, &status);

      si.cinp->pinp[j].alloc = 1;
      si.cinp->pinp[j].accelerator = si.cinp->rinp[j].accelerator;
      MPI_Recv(&si.cinp->pinp[j].typeCpu, 1, MPI_INT, 0, rank, Comm, &status);
      MPI_Recv(&si.cinp->pinp[j].cpuPmin, 1, MPI_DOUBLE, 0, rank, Comm, &status);
      MPI_Recv(&si.cinp->pinp[j].cpuPmax, 1, MPI_DOUBLE, 0, rank, Comm, &status);
      MPI_Recv(&si.cinp->pinp[j].numOfPoints, 1, MPI_INT, 0, rank, Comm, &status);
      if (si.cinp->pinp[j].numOfPoints > 0) {
        si.cinp->pinp[j].cpubins = new double[si.cinp->pinp[j].numOfPoints];
        MPI_Recv(si.cinp->pinp[j].cpubins, si.cinp->pinp[j].numOfPoints, MPI_DOUBLE, 0, rank, Comm, &status);
        si.cinp->pinp[j].cpuP = new double[si.cinp->pinp[j].numOfPoints];
        MPI_Recv(si.cinp->pinp[j].cpuP, si.cinp->pinp[j].numOfPoints, MPI_DOUBLE, 0, rank, Comm, &status);
      }
      MPI_Recv(&si.cinp->pinp[j].cpuC, 1, MPI_DOUBLE, 0, rank, Comm, &status);
      MPI_Recv(&si.cinp->pinp[j].typeAcc, 1, MPI_INT, 0, rank, Comm, &status);
      MPI_Recv(&si.cinp->pinp[j].accPmin, 1, MPI_DOUBLE, 0, rank, Comm, &status);
      MPI_Recv(&si.cinp->pinp[j].accPmax, 1, MPI_DOUBLE, 0, rank, Comm, &status);
      MPI_Recv(&si.cinp->pinp[j].accC, 1, MPI_DOUBLE, 0, rank, Comm, &status);
    }
  }
}

void communicator::taskParameters(list<task>& jobs, const int& rank, const int& clusterSize, const int* commCell,
                                  const MPI_Comm& Comm)
{
  int i, j, k, who;
  int L_type, L_numberOfAvailableImplementations;
  int* L_availableImplementations;
  double L_requestedInstructions;
  int L_numberOfVMs;
  double* L_reqPMNS;
  int* L_typeactPMN;
  double** L_minmaxactPMN;
  int* L_avAcc;
  double* L_rhoAcc;
  MPI_Status status;
  int establish;
  if (rank == 0) {
    establish = jobs.size();
    for (i = 0; i < clusterSize - 1; i++) {
      MPI_Send(&establish, 1, MPI_INT, i + 1, i + 1, Comm);
    }
  } else {
    MPI_Recv(&establish, 1, MPI_INT, 0, rank, Comm, &status);
  }
  if (!establish)
    return;
  if (rank == 0) {
    list<task>::iterator it = jobs.begin();
    for (i = 0; i < establish; i++) {
      who = commCell[i];
      for (j = 0; j < clusterSize - 1; j++) {
        MPI_Send(&who, 1, MPI_INT, j + 1, j + 1, Comm);
      }
      if (who > 0) {
        L_type = (*it).getType();
        MPI_Send(&L_type, 1, MPI_INT, who, who, Comm);

        L_numberOfAvailableImplementations = (*it).getNumberOfAvailableImplementations();
        MPI_Send(&L_numberOfAvailableImplementations, 1, MPI_INT, who, who, Comm);

        L_availableImplementations = new int[L_numberOfAvailableImplementations];
        for (j = 0; j < L_numberOfAvailableImplementations; j++)
          L_availableImplementations[j] = (*it).getAvailableImplementations()[j];
        MPI_Send(L_availableImplementations, L_numberOfAvailableImplementations, MPI_INT, who, who, Comm);

        L_requestedInstructions = (*it).grequestedInstructions();
        MPI_Send(&L_requestedInstructions, 1, MPI_DOUBLE, who, who, Comm);

        L_numberOfVMs = (*it).getNumberOfVMs();
        MPI_Send(&L_numberOfVMs, 1, MPI_INT, who, who, Comm);

        L_reqPMNS = new double[4];
        for (j = 0; j < 4; j++)
          L_reqPMNS[j] = (*it).greqPMNS()[j];
        MPI_Send(L_reqPMNS, 4, MPI_DOUBLE, who, who, Comm);

        L_typeactPMN = new int[3];
        for (j = 0; j < 3; j++)
          L_typeactPMN[j] = (*it).getTypeactPMN()[j];
        MPI_Send(L_typeactPMN, 3, MPI_INT, who, who, Comm);

        L_minmaxactPMN = new double*[3];
        for (j = 0; j < 3; j++)
          L_minmaxactPMN[j] = new double[2];
        for (j = 0; j < 3; j++) {
          for (k = 0; k < 2; k++) {
            L_minmaxactPMN[j][k] = (*it).gminmaxactPMN()[j][k];
          }
        }
        for (j = 0; j < 3; j++) {
          MPI_Send(&L_minmaxactPMN[j][0], 2, MPI_DOUBLE, who, who, Comm);
        }

        L_avAcc = new int[L_numberOfAvailableImplementations];
        for (j = 0; j < L_numberOfAvailableImplementations; j++) {
          L_avAcc[j] = (*it).gavAcc()[j];
        }
        MPI_Send(L_avAcc, L_numberOfAvailableImplementations, MPI_INT, who, who, Comm);

        L_rhoAcc = new double[L_numberOfAvailableImplementations];
        for (j = 0; j < L_numberOfAvailableImplementations; j++) {
          L_rhoAcc[j] = (*it).grhoAcc()[j];
        }
        MPI_Send(L_rhoAcc, L_numberOfAvailableImplementations, MPI_DOUBLE, who, who, Comm);

        delete[] L_availableImplementations;
        delete[] L_reqPMNS;
        delete[] L_typeactPMN;
        for (j = 0; j < 3; j++) {
          delete[] L_minmaxactPMN[j];
        }
        delete[] L_minmaxactPMN;
        delete[] L_avAcc;
        delete[] L_rhoAcc;
      }
      it++;
    }
  } else {
    for (i = 0; i < establish; i++) {
      MPI_Recv(&who, 1, MPI_INT, 0, rank, Comm, &status);
      if (who == rank) {
        MPI_Recv(&L_type, 1, MPI_INT, 0, rank, Comm, &status);
        MPI_Recv(&L_numberOfAvailableImplementations, 1, MPI_INT, 0, rank, Comm, &status);
        L_availableImplementations = new int[L_numberOfAvailableImplementations];
        MPI_Recv(L_availableImplementations, L_numberOfAvailableImplementations, MPI_INT, 0, rank, Comm, &status);
        MPI_Recv(&L_requestedInstructions, 1, MPI_DOUBLE, 0, rank, Comm, &status);
        MPI_Recv(&L_numberOfVMs, 1, MPI_INT, 0, rank, Comm, &status);
        L_reqPMNS = new double[4];
        MPI_Recv(L_reqPMNS, 4, MPI_DOUBLE, 0, rank, Comm, &status);
        L_typeactPMN = new int[3];
        MPI_Recv(L_typeactPMN, 3, MPI_INT, 0, rank, Comm, &status);
        L_minmaxactPMN = new double*[3];

        for (j = 0; j < 3; j++) {
          L_minmaxactPMN[j] = new double[2];
        }
        for (j = 0; j < 3; j++) {
          MPI_Recv(&L_minmaxactPMN[j][0], 2, MPI_DOUBLE, 0, rank, Comm, &status);
        }
        L_avAcc = new int[L_numberOfAvailableImplementations];
        MPI_Recv(L_avAcc, L_numberOfAvailableImplementations, MPI_INT, 0, rank, Comm, &status);
        L_rhoAcc = new double[L_numberOfAvailableImplementations];
        MPI_Recv(L_rhoAcc, L_numberOfAvailableImplementations, MPI_DOUBLE, 0, rank, Comm, &status);

        jobs.push_back(task(L_type, L_numberOfAvailableImplementations, L_availableImplementations,
                            L_requestedInstructions, L_numberOfVMs, L_reqPMNS[0], L_reqPMNS[1], L_reqPMNS[2],
                            L_reqPMNS[3], L_typeactPMN[0], L_typeactPMN[1], L_typeactPMN[2], &L_minmaxactPMN[0][0],
                            &L_minmaxactPMN[1][0], &L_minmaxactPMN[2][0], L_avAcc, L_rhoAcc));

        delete[] L_availableImplementations;
        delete[] L_reqPMNS;
        delete[] L_typeactPMN;
        for (j = 0; j < 3; j++) {
          delete[] L_minmaxactPMN[j];
        }
        delete[] L_minmaxactPMN;
        delete[] L_avAcc;
        delete[] L_rhoAcc;
      }
    }
  }
}

void communicator::cellStatistics(const gs* gates, const cell* clCell, const int& rank, const int& clusterSize,
                                  const MPI_Comm& Comm)
{
  int i, j;
  MPI_Status status;
  if (rank == 0) {
    for (i = 0; i < clusterSize - 1; i++) {
      for (j = 0; j < gates[0].gsi()[0].cinp[i].numberOfTypes; j++) {
        MPI_Recv(&gates[0].getStats()[i][j].alloc, 1, MPI_INT, i + 1, i + 1, Comm, &status);
        MPI_Recv(&gates[0].getStats()[i][j].currentTimestep, 1, MPI_DOUBLE, i + 1, i + 1, Comm, &status);
        MPI_Recv(&gates[0].getStats()[i][j].processorsOverActiveServers, 1, MPI_INT, i + 1, i + 1, Comm, &status);
        MPI_Recv(&gates[0].getStats()[i][j].memoryOverActiveServers, 1, MPI_DOUBLE, i + 1, i + 1, Comm, &status);
        MPI_Recv(&gates[0].getStats()[i][j].storageOverActiveServers, 1, MPI_DOUBLE, i + 1, i + 1, Comm, &status);
        MPI_Recv(&gates[0].getStats()[i][j].acceleratorsOverActiveServers, 1, MPI_INT, i + 1, i + 1, Comm, &status);
        MPI_Recv(&gates[0].getStats()[i][j].physicalMemory, 1, MPI_DOUBLE, i + 1, i + 1, Comm, &status);
        MPI_Recv(&gates[0].getStats()[i][j].physicalProcessors, 1, MPI_DOUBLE, i + 1, i + 1, Comm, &status);
        MPI_Recv(&gates[0].getStats()[i][j].physicalStorage, 1, MPI_DOUBLE, i + 1, i + 1, Comm, &status);
        MPI_Recv(&gates[0].getStats()[i][j].physicalNetwork, 1, MPI_DOUBLE, i + 1, i + 1, Comm, &status);
        MPI_Recv(&gates[0].getStats()[i][j].totalMemory, 1, MPI_DOUBLE, i + 1, i + 1, Comm, &status);
        MPI_Recv(&gates[0].getStats()[i][j].totalProcessors, 1, MPI_DOUBLE, i + 1, i + 1, Comm, &status);
        MPI_Recv(&gates[0].getStats()[i][j].totalStorage, 1, MPI_DOUBLE, i + 1, i + 1, Comm, &status);
        MPI_Recv(&gates[0].getStats()[i][j].availableMemory, 1, MPI_DOUBLE, i + 1, i + 1, Comm, &status);
        MPI_Recv(&gates[0].getStats()[i][j].availableProcessors, 1, MPI_DOUBLE, i + 1, i + 1, Comm, &status);
        MPI_Recv(&gates[0].getStats()[i][j].availableStorage, 1, MPI_DOUBLE, i + 1, i + 1, Comm, &status);

        MPI_Recv(&gates[0].getStats()[i][j].actualUtilizedMemory, 1, MPI_DOUBLE, i + 1, i + 1, Comm, &status);
        MPI_Recv(&gates[0].getStats()[i][j].actualUtilizedProcessors, 1, MPI_DOUBLE, i + 1, i + 1, Comm, &status);

        gates[0].getStats()[i][j].utilizedProcessors =
          gates[0].getStats()[i][j].totalProcessors - gates[0].getStats()[i][j].availableProcessors;
        gates[0].getStats()[i][j].utilizedMemory =
          gates[0].getStats()[i][j].totalMemory - gates[0].getStats()[i][j].availableMemory;
        gates[0].getStats()[i][j].utilizedStorage =
          gates[0].getStats()[i][j].totalStorage - gates[0].getStats()[i][j].availableStorage;

        //				MPI_Recv(&gates[0].getStats()[i][j].utilizedProcessors,1,MPI_DOUBLE,i+1,i+1,Comm,&status);
        //				MPI_Recv(&gates[0].getStats()[i][j].utilizedMemory,1,MPI_DOUBLE,i+1,i+1,Comm,&status);
        //				MPI_Recv(&gates[0].getStats()[i][j].utilizedStorage,1,MPI_DOUBLE,i+1,i+1,Comm,&status);

        MPI_Recv(&gates[0].getStats()[i][j].totalNetwork, 1, MPI_DOUBLE, i + 1, i + 1, Comm, &status);
        MPI_Recv(&gates[0].getStats()[i][j].availableNetwork, 1, MPI_DOUBLE, i + 1, i + 1, Comm, &status);

        MPI_Recv(&gates[0].getStats()[i][j].actualUtilizedNetwork, 1, MPI_DOUBLE, i + 1, i + 1, Comm, &status);

        //				MPI_Recv(&gates[0].getStats()[i][j].utilizedNetwork,1,MPI_DOUBLE,i+1,i+1,Comm,&status);

        gates[0].getStats()[i][j].utilizedNetwork =
          gates[0].getStats()[i][j].totalNetwork - gates[0].getStats()[i][j].availableNetwork;

        MPI_Recv(&gates[0].getStats()[i][j].totalPowerConsumption, 1, MPI_DOUBLE, i + 1, i + 1, Comm, &status);
        MPI_Recv(&gates[0].getStats()[i][j].totalAccelerators, 1, MPI_INT, i + 1, i + 1, Comm, &status);
        MPI_Recv(&gates[0].getStats()[i][j].availableAccelerators, 1, MPI_INT, i + 1, i + 1, Comm, &status);

        //				MPI_Recv(&gates[0].getStats()[i][j].utilizedAccelerators,1,MPI_INT,i+1,i+1,Comm,&status);

        gates[0].getStats()[i][j].utilizedAccelerators =
          gates[0].getStats()[i][j].totalAccelerators - gates[0].getStats()[i][j].availableAccelerators;

        MPI_Recv(&gates[0].getStats()[i][j].activeServers, 1, MPI_INT, i + 1, i + 1, Comm, &status);
        MPI_Recv(&gates[0].getStats()[i][j].runningVMs, 1, MPI_INT, i + 1, i + 1, Comm, &status);
        MPI_Recv(&gates[0].getStats()[i][j].rejectedTasks, 1, MPI_INT, i + 1, i + 1, Comm, &status);
        MPI_Recv(&gates[0].getStats()[i][j].acceptedTasks, 1, MPI_INT, i + 1, i + 1, Comm, &status);
      }
    }
  } else {
    for (i = 0; i < clCell[0].getNumberOfTypes(); i++) {
      MPI_Send(&clCell[0].getStats()[i].alloc, 1, MPI_INT, 0, rank, Comm);
      MPI_Send(&clCell[0].getStats()[i].currentTimestep, 1, MPI_DOUBLE, 0, rank, Comm);
      MPI_Send(&clCell[0].getStats()[i].processorsOverActiveServers, 1, MPI_INT, 0, rank, Comm);
      MPI_Send(&clCell[0].getStats()[i].memoryOverActiveServers, 1, MPI_DOUBLE, 0, rank, Comm);
      MPI_Send(&clCell[0].getStats()[i].storageOverActiveServers, 1, MPI_DOUBLE, 0, rank, Comm);
      MPI_Send(&clCell[0].getStats()[i].acceleratorsOverActiveServers, 1, MPI_INT, 0, rank, Comm);
      MPI_Send(&clCell[0].getStats()[i].physicalMemory, 1, MPI_DOUBLE, 0, rank, Comm);
      MPI_Send(&clCell[0].getStats()[i].physicalProcessors, 1, MPI_DOUBLE, 0, rank, Comm);
      MPI_Send(&clCell[0].getStats()[i].physicalStorage, 1, MPI_DOUBLE, 0, rank, Comm);
      MPI_Send(&clCell[0].getStats()[i].physicalNetwork, 1, MPI_DOUBLE, 0, rank, Comm);
      MPI_Send(&clCell[0].getStats()[i].totalMemory, 1, MPI_DOUBLE, 0, rank, Comm);
      MPI_Send(&clCell[0].getStats()[i].totalProcessors, 1, MPI_DOUBLE, 0, rank, Comm);
      MPI_Send(&clCell[0].getStats()[i].totalStorage, 1, MPI_DOUBLE, 0, rank, Comm);
      MPI_Send(&clCell[0].getStats()[i].availableMemory, 1, MPI_DOUBLE, 0, rank, Comm);
      MPI_Send(&clCell[0].getStats()[i].availableProcessors, 1, MPI_DOUBLE, 0, rank, Comm);
      MPI_Send(&clCell[0].getStats()[i].availableStorage, 1, MPI_DOUBLE, 0, rank, Comm);

      MPI_Send(&clCell[0].getStats()[i].actualUtilizedMemory, 1, MPI_DOUBLE, 0, rank, Comm);
      MPI_Send(&clCell[0].getStats()[i].actualUtilizedProcessors, 1, MPI_DOUBLE, 0, rank, Comm);

      //			MPI_Send(&clCell[0].getStats()[i].utilizedProcessors,1,MPI_DOUBLE,0,rank,Comm);
      //			MPI_Send(&clCell[0].getStats()[i].utilizedMemory,1,MPI_DOUBLE,0,rank,Comm);
      //			MPI_Send(&clCell[0].getStats()[i].utilizedStorage,1,MPI_DOUBLE,0,rank,Comm);

      MPI_Send(&clCell[0].getStats()[i].totalNetwork, 1, MPI_DOUBLE, 0, rank, Comm);
      MPI_Send(&clCell[0].getStats()[i].availableNetwork, 1, MPI_DOUBLE, 0, rank, Comm);

      MPI_Send(&clCell[0].getStats()[i].actualUtilizedNetwork, 1, MPI_DOUBLE, 0, rank, Comm);

      //			MPI_Send(&clCell[0].getStats()[i].utilizedNetwork,1,MPI_DOUBLE,0,rank,Comm);

      MPI_Send(&clCell[0].getStats()[i].totalPowerConsumption, 1, MPI_DOUBLE, 0, rank, Comm);
      MPI_Send(&clCell[0].getStats()[i].totalAccelerators, 1, MPI_INT, 0, rank, Comm);
      MPI_Send(&clCell[0].getStats()[i].availableAccelerators, 1, MPI_INT, 0, rank, Comm);

      //			MPI_Send(&clCell[0].getStats()[i].utilizedAccelerators,1,MPI_INT,0,rank,Comm);

      MPI_Send(&clCell[0].getStats()[i].activeServers, 1, MPI_INT, 0, rank, Comm);
      MPI_Send(&clCell[0].getStats()[i].runningVMs, 1, MPI_INT, 0, rank, Comm);
      MPI_Send(&clCell[0].getStats()[i].rejectedTasks, 1, MPI_INT, 0, rank, Comm);
      MPI_Send(&clCell[0].getStats()[i].acceptedTasks, 1, MPI_INT, 0, rank, Comm);
    }
  }
}

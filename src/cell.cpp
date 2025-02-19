#include <cell.h>
#include <inputs.h>     // for cellinputs, brinputs, netinputs, powi...
#include <netw.h>       // for netw
#include <power.h>      // for power
#include <resource.h>   // for resource
#include <sosmBroker.h> // for sosmBroker
#include <stat.h>       // for stat
#include <task.h>
#include <traditionalBroker.h> // for traditionalBroker
#include <cstdlib>             // for atoi, getenv
#include <iostream>            // for cout, endl

using std::cout;
using std::endl;

cell::cell()
  : ID(0),
    alloc(0),
    types(nullptr),
    numberOfTypes(0),
    sosmIntegration(0),
    numberOfResourcesPerType(nullptr),
    resources(nullptr),
    powerComp(nullptr),
    network(nullptr),
    broker(nullptr),
    stats(nullptr)
{
}

cell::cell(const cellinputs& setup, int _sosmIntegration) : alloc(1), sosmIntegration(_sosmIntegration)
{
  ID = setup.ID;
  numberOfTypes = setup.numberOfTypes;
  types = new int[numberOfTypes];
  numberOfResourcesPerType = new int[numberOfTypes];

  for (int i = 0; i < numberOfTypes; i++) {
    types[i] = setup.types[i];
    numberOfResourcesPerType[i] = setup.numberOfResourcesPerType[i];
  }

  resources = new resource*[numberOfTypes];
  for (int i = 0; i < numberOfTypes; i++) {
    resources[i] = new resource[numberOfResourcesPerType[i]];
  }
  for (int i = 0; i < numberOfTypes; i++) {
    for (int j = 0; j < numberOfResourcesPerType[i]; j++) {
      // Copy information about computer resources in the resources[hardware type][resource] array
      resources[i][j] = resource(setup.rinp[i], j);
    }
  }

  powerComp = new power[numberOfTypes];
  for (int i = 0; i < numberOfTypes; i++) {
    powerComp[i] = power(setup.pinp[i]);
  }
  network = new netw[1];
  network[0] = netw(setup.ninp[0]);

  // Initialize the *broker classes polymorphism based on the value of sosmIntegration
  if (sosmIntegration) {
    broker = new sosmBroker[1];
  } else {
    broker = new traditionalBroker[1];
  }

  stats = new stat[numberOfTypes];
  for (int i = 0; i < numberOfTypes; i++) {
    stats[i] = stat();
  }

  updateStats(0.0);
}

cell::cell(const cell& t)
{
  if (t.galloc()) {
    ID = t.gID();
    alloc = t.galloc();
    numberOfTypes = t.getNumberOfTypes();
    sosmIntegration = t.getSosmIntegration();
    types = new int[numberOfTypes];
    numberOfResourcesPerType = new int[numberOfTypes];

    for (int i = 0; i < numberOfTypes; i++) {
      types[i] = t.getTypes()[i];
      numberOfResourcesPerType[i] = t.getNumberOfResourcesPerType()[i];
    }

    resources = new resource*[numberOfTypes];
    for (int i = 0; i < numberOfTypes; i++) {
      resources[i] = new resource[numberOfResourcesPerType[i]];
    }
    for (int i = 0; i < numberOfTypes; i++) {
      for (int j = 0; j < numberOfResourcesPerType[i]; j++) {
        resources[i][j] = t.getResources()[i][j];
      }
    }

    powerComp = new power[numberOfTypes];
    for (int i = 0; i < numberOfTypes; i++) {
      powerComp[i] = t.getPowerConsumption()[i];
    }

    network = new netw[1];
    network[0] = t.getNetwork()[0];

    if (sosmIntegration) {
      broker = new sosmBroker[1];
    } else {
      broker = new traditionalBroker[1];
    }
    broker[0] = t.getBroker()[0];

    stats = new stat[numberOfTypes];

    for (int i = 0; i < numberOfTypes; i++) {
      stats[i] = t.getStats()[i];
    }
  }
}

cell& cell::operator=(const cell& t)
{
  if (this != &t) {
    if (alloc) {
      ID = 0;
      alloc = 0;
      delete[] types;
      types = nullptr;
      delete[] numberOfResourcesPerType;
      numberOfResourcesPerType = nullptr;

      for (int i = 0; i < numberOfTypes; i++) {
        delete[] resources[i];
      }

      delete[] resources;
      delete[] powerComp;
      delete[] network;
      delete[] broker;
      delete[] stats;
      resources = nullptr;
      powerComp = nullptr;
      network = nullptr;
      broker = nullptr;
      stats = nullptr;
      numberOfTypes = 0;
      sosmIntegration = 0;
    }
    alloc = t.galloc();
    if (alloc) {
      ID = t.gID();
      numberOfTypes = t.getNumberOfTypes();
      sosmIntegration = t.getSosmIntegration();
      types = new int[numberOfTypes];

      numberOfResourcesPerType = new int[numberOfTypes];
      for (int i = 0; i < numberOfTypes; i++) {
        types[i] = t.getTypes()[i];
        numberOfResourcesPerType[i] = t.getNumberOfResourcesPerType()[i];
      }

      resources = new resource*[numberOfTypes];
      for (int i = 0; i < numberOfTypes; i++) {
        resources[i] = new resource[numberOfResourcesPerType[i]];
      }
      for (int i = 0; i < numberOfTypes; i++) {
        for (int j = 0; j < numberOfResourcesPerType[i]; j++)
          resources[i][j] = t.getResources()[i][j];
      }

      powerComp = new power[numberOfTypes];
      for (int i = 0; i < numberOfTypes; i++) {
        powerComp[i] = t.getPowerConsumption()[i];
      }

      network = new netw[1];
      network[0] = t.getNetwork()[0];

      if (sosmIntegration) {
        broker = new sosmBroker[1];
      } else {
        broker = new traditionalBroker[1];
      }
      broker[0] = t.getBroker()[0];

      stats = new stat[numberOfTypes];
      for (int i = 0; i < numberOfTypes; i++) {
        stats[i] = t.getStats()[i];
      }
    }
  }
  return *this;
}

cell::~cell()
{
  if (alloc) {
    ID = 0;
    alloc = 0;

    delete[] types;
    delete[] numberOfResourcesPerType;
    types = nullptr;
    numberOfResourcesPerType = nullptr;
    for (int i = 0; i < numberOfTypes; i++) {
      delete[] resources[i];
    }
    delete[] resources;
    delete[] powerComp;
    delete[] network;
    delete[] broker;
    delete[] stats;
    resources = nullptr;
    powerComp = nullptr;
    network = nullptr;
    broker = nullptr;
    stats = nullptr;
    numberOfTypes = 0;
    sosmIntegration = 0;
  }
}

int cell::gID() const { return ID; }
int cell::galloc() const { return alloc; }
int cell::getNumberOfTypes() const { return numberOfTypes; }
int cell::getSosmIntegration() const { return sosmIntegration; }
int* cell::getTypes() const { return types; }
int* cell::getNumberOfResourcesPerType() const { return numberOfResourcesPerType; }
resource** cell::getResources() const { return resources; }
power* cell::getPowerConsumption() const { return powerComp; }
baseBroker* cell::getBroker() const { return broker; }
netw* cell::getNetwork() const { return network; }
stat* cell::getStats() const { return stats; }
void cell::deploy(list<task>& jobs)
{
  if (alloc) {
    for (auto& it : jobs) {
      broker[0].deploy(resources, network, stats, it);
    }
  }
}

void cell::updateStats(const double& tstep)
{
  int i, j;
  int omp_thr = atoi(getenv("OMP_NUM_THREADS"));
  for (i = 0; i < numberOfTypes; i++) {
    stats[i].alloc = 1;
    stats[i].currentTimestep = tstep;

    stats[i].totalNetwork = network[0].getTotalNetwork();
    stats[i].availableNetwork = network[0].getAvailableNetwork();
    stats[i].utilizedNetwork = stats[i].totalNetwork - stats[i].availableNetwork;
    stats[i].actualUtilizedNetwork = network[0].getActualUtilizedNetwork();

    int processorsPerServer = resources[i][0].getTotalProcessors();
    int memoryPerServer = resources[i][0].getTotalMemory();
    int storagePerServer = resources[i][0].getTotalStorage();
    int acceleratorsPerServer = resources[i][0].getTotalAccelerators();

    int processorsOverActiveServers = 0;
    int memoryOverActiveServers = 0;
    int storageOverActiveServers = 0;
    int acceleratorsOverActiveServers = 0;

    double availableMemory = 0.0;
    double physicalMemory = 0.0;
    double totalMemory = 0.0;

    double availableProcessors = 0.0;
    double physicalProcessors = 0.0;
    double totalProcessors = 0.0;

    double availableStorage = 0.0;
    double physicalStorage = 0.0;
    double totalStorage = 0.0;

    int availableAccelerators = 0;
    int totalAccelerators = 0;

    int activeServers = 0;
    int runningVMs = 0;

    double actualUtilizedProcessors = 0.0;
    double actualUtilizedMemory = 0.0;

#pragma omp parallel for default(shared) private(j) num_threads(omp_thr) schedule(static)              \
  reduction(+ : physicalProcessors, totalProcessors, availableProcessors, physicalMemory, totalMemory, \
            availableMemory, physicalStorage, totalStorage, availableStorage, totalAccelerators,       \
            availableAccelerators, activeServers, runningVMs, actualUtilizedProcessors, actualUtilizedMemory)
    for (j = 0; j < numberOfResourcesPerType[i]; j++) {
      physicalProcessors += resources[i][j].getPhysicalProcessors();
      totalProcessors += resources[i][j].getTotalProcessors();
      availableProcessors += resources[i][j].getAvailableProcessors();
      physicalMemory += resources[i][j].getPhysicalMemory();
      totalMemory += resources[i][j].getTotalMemory();
      availableMemory += resources[i][j].getAvailableMemory();
      physicalStorage += resources[i][j].getPhysicalStorage();
      totalStorage += resources[i][j].getTotalStorage();
      availableStorage += resources[i][j].getAvailableStorage();
      totalAccelerators += resources[i][j].getTotalAccelerators();
      availableAccelerators += resources[i][j].getAvailableAccelerators();
      activeServers += resources[i][j].getActive();
      runningVMs += resources[i][j].getRunningVMs();

      actualUtilizedProcessors += resources[i][j].getActualUtilizedProcessors();
      actualUtilizedMemory += resources[i][j].getActualUtilizedMemory();
    }
    stats[i].physicalProcessors = physicalProcessors;
    stats[i].totalProcessors = totalProcessors;
    stats[i].availableProcessors = availableProcessors;
    stats[i].utilizedProcessors = totalProcessors - availableProcessors;
    stats[i].physicalMemory = physicalMemory;
    stats[i].totalMemory = totalMemory;
    stats[i].availableMemory = availableMemory;
    stats[i].utilizedMemory = totalMemory - availableMemory;
    stats[i].physicalStorage = physicalStorage;
    stats[i].totalStorage = totalStorage;
    stats[i].availableStorage = availableStorage;
    stats[i].utilizedStorage = totalStorage - availableStorage;
    stats[i].totalAccelerators = totalAccelerators;
    stats[i].availableAccelerators = availableAccelerators;
    stats[i].utilizedAccelerators = totalAccelerators - availableAccelerators;
    stats[i].activeServers = activeServers;
    stats[i].runningVMs = runningVMs;
    stats[i].actualUtilizedProcessors = actualUtilizedProcessors;
    stats[i].actualUtilizedMemory = actualUtilizedMemory;

    processorsOverActiveServers = activeServers * processorsPerServer;
    memoryOverActiveServers = activeServers * memoryPerServer;
    storageOverActiveServers = activeServers * storagePerServer;
    acceleratorsOverActiveServers = activeServers * acceleratorsPerServer;

    stats[i].processorsOverActiveServers = processorsOverActiveServers;
    stats[i].memoryOverActiveServers = memoryOverActiveServers;
    stats[i].storageOverActiveServers = storageOverActiveServers;
    stats[i].acceleratorsOverActiveServers = acceleratorsOverActiveServers;
  }
}

void cell::print()
{
  if (alloc) {
    cout << "Cell ID: " << ID << endl;
    cout << "Number of HW types: " << numberOfTypes << endl;
    cout << "HW types: ";

    for (int i = 0; i < numberOfTypes; i++) {
      cout << types[i] << " ";
    }
    cout << endl;

    cout << "Number of Resources Per Type: ";
    for (int i = 0; i < numberOfTypes; i++) {
      cout << numberOfResourcesPerType[i] << " ";
    }
    cout << endl;

    cout << "---------------------------------------------" << endl;

    network[0].print();
    broker[0].print();

    for (int i = 0; i < numberOfTypes; i++) {
      cout << "     Resource Type: " << types[i] << endl;
      powerComp[i].print();
      stats[i].print();
    }
  }
}

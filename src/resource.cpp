#include <resource.h>
#include <algorithm> // for max
#include <iostream>  // for cout, endl
#include "inputs.h"  // for resinputs
#include "task.h"    // for task

using std::cout;
using std::endl;
using std::max;

resource::resource()
  : alloc(0),
    active(0),
    movable(0),
    type(-1),
    ID(-1),
    totalProcessors(0.0),
    availableProcessors(0.0),
    utilizedProcessors(0.0),
    totalMemory(0.0),
    availableMemory(0.0),
    utilizedMemory(0.0),
    totalStorage(0.0),
    availableStorage(0.0),
    utilizedStorage(0.0),
    physicalProcessors(0.0),
    physicalMemory(0.0),
    physicalStorage(0.0),
    computeCapability(0.0),
    accelerator(0),
    totalAccelerators(0),
    availableAccelerators(0),
    utilizedAccelerators(0),
    acceleratorComputeCapability(0.0),
    overcommitmentProcessors(1.0),
    overcommitmentMemory(1.0),
    actualUtilizedProcessors(0.0),
    actualUtilizedMemory(0.0),
    actualRhoAccelerators(0.0),
    runningVMs(0),
    currentCompCapPerProc(0.0),
    currentCompCapPerAcc(0.0)
{
}

resource::resource(const resinputs& setup, const int& _ID)
  : alloc(1),
    active(0),
    movable(1),
    utilizedProcessors(0.0),
    utilizedMemory(0.0),
    utilizedStorage(0.0),
    utilizedAccelerators(0),
    actualUtilizedProcessors(0.0),
    actualUtilizedMemory(0.0),
    actualRhoAccelerators(0.0),
    runningVMs(0),
    currentCompCapPerProc(0.0),
    currentCompCapPerAcc(0.0)
{
  type = setup.type;
  ID = _ID;

  physicalProcessors = setup.numOfProcUnits;
  physicalMemory = setup.totalMemory;
  physicalStorage = setup.totalStorage;

  overcommitmentProcessors = setup.overcommitmentProcessors;
  totalProcessors = ((double)setup.numOfProcUnits) * overcommitmentProcessors;
  availableProcessors = totalProcessors;

  overcommitmentMemory = setup.overcommitmentMemory;
  totalMemory = ((double)setup.totalMemory) * overcommitmentMemory;
  availableMemory = totalMemory;

  totalStorage = ((double)setup.totalStorage);
  availableStorage = totalStorage;

  computeCapability = setup.computeCapability;

  accelerator = setup.accelerator;
  acceleratorComputeCapability = setup.acceleratorComputeCapability;
  totalAccelerators = setup.totalAccelerators;
  availableAccelerators = totalAccelerators;
}

resource::resource(const resource& t)
{
  if (t.galloc()) {
    alloc = t.galloc();
    active = t.getActive();
    movable = t.getMovable();
    type = t.getType();
    ID = t.gID();

    totalProcessors = t.getTotalProcessors();
    availableProcessors = t.getAvailableProcessors();
    utilizedProcessors = t.getUtilizedProcessors();

    totalMemory = t.getTotalMemory();
    availableMemory = t.getAvailableMemory();
    utilizedMemory = t.getUtilizedMemory();

    physicalProcessors = t.getPhysicalProcessors();
    physicalMemory = t.getPhysicalMemory();

    computeCapability = t.getComputeCapability();

    accelerator = t.getAccelerator();
    totalAccelerators = t.getTotalAccelerators();
    utilizedAccelerators = t.getUtilizedAccelerators();
    availableAccelerators = t.getAvailableAccelerators();
    acceleratorComputeCapability = t.getAcceleratorComputeCapability();

    overcommitmentProcessors = t.getOvercommitmentProcessors();
    overcommitmentMemory = t.getOvercommitmentMemory();

    actualUtilizedProcessors = t.getActualUtilizedProcessors();
    actualUtilizedMemory = t.getActualUtilizedMemory();
    actualRhoAccelerators = t.getActualRhoAccelerators();

    runningVMs = t.getRunningVMs();
    currentCompCapPerProc = t.getCurrentCompCapPerProc();
    currentCompCapPerAcc = t.getCurrentCompCapPerAcc();
  }
}

resource& resource::operator=(const resource& t)
{
  if (this != &t) {
    if (alloc) {
      alloc = 0;
      active = 0;
      movable = 0;
      type = -1;
      ID = -1;

      totalProcessors = 0.0;
      availableProcessors = 0.0;
      utilizedProcessors = 0.0;

      totalMemory = 0.0;
      availableMemory = 0.0;
      utilizedMemory = 0.0;

      totalStorage = 0.0;
      availableStorage = 0.0;
      utilizedStorage = 0.0;

      physicalProcessors = 0.0;
      physicalMemory = 0.0;
      physicalStorage = 0.0;

      computeCapability = 0.0;

      accelerator = 0;
      totalAccelerators = 0;
      utilizedAccelerators = 0;
      availableAccelerators = 0;
      acceleratorComputeCapability = 0.0;

      overcommitmentProcessors = 1.0;
      overcommitmentMemory = 1.0;

      actualUtilizedProcessors = 0.0;
      actualUtilizedMemory = 0.0;
      actualRhoAccelerators = 0.0;

      runningVMs = 0;
      currentCompCapPerProc = 0.0;
      currentCompCapPerAcc = 0.0;
    }
    alloc = t.galloc();
    if (alloc) {
      active = t.getActive();
      movable = t.getMovable();
      type = t.getType();
      ID = t.gID();

      totalProcessors = t.getTotalProcessors();
      availableProcessors = t.getAvailableProcessors();
      utilizedProcessors = t.getUtilizedProcessors();

      totalMemory = t.getTotalMemory();
      availableMemory = t.getAvailableMemory();
      utilizedMemory = t.getUtilizedMemory();

      totalStorage = t.getTotalStorage();
      availableStorage = t.getAvailableStorage();
      utilizedStorage = t.getUtilizedStorage();

      physicalProcessors = t.getPhysicalProcessors();
      physicalMemory = t.getPhysicalMemory();
      physicalStorage = t.getPhysicalStorage();

      computeCapability = t.getComputeCapability();

      accelerator = t.getAccelerator();
      totalAccelerators = t.getTotalAccelerators();
      availableAccelerators = t.getAvailableAccelerators();
      utilizedAccelerators = t.getUtilizedAccelerators();
      acceleratorComputeCapability = t.getAcceleratorComputeCapability();

      overcommitmentProcessors = t.getOvercommitmentProcessors();
      overcommitmentMemory = t.getOvercommitmentMemory();

      actualUtilizedProcessors = t.getActualUtilizedProcessors();
      actualUtilizedMemory = t.getActualUtilizedMemory();
      actualRhoAccelerators = t.getActualRhoAccelerators();

      runningVMs = t.getRunningVMs();
      currentCompCapPerProc = t.getCurrentCompCapPerProc();
      currentCompCapPerAcc = t.getCurrentCompCapPerAcc();
    }
  }
  return *this;
}

resource::~resource()
{
  if (alloc) {
    alloc = 0;
    active = 0;
    movable = 0;
    type = -1;
    ID = -1;

    totalProcessors = 0.0;
    availableProcessors = 0.0;
    utilizedProcessors = 0.0;

    totalMemory = 0.0;
    availableMemory = 0.0;
    utilizedMemory = 0.0;

    totalStorage = 0.0;
    availableStorage = 0.0;
    utilizedStorage = 0.0;

    physicalProcessors = 0.0;
    physicalMemory = 0.0;
    physicalStorage = 0.0;

    computeCapability = 0.0;

    accelerator = 0;
    totalAccelerators = 0;
    utilizedAccelerators = 0;
    availableAccelerators = 0;
    acceleratorComputeCapability = 0.0;

    overcommitmentProcessors = 1.0;
    overcommitmentMemory = 1.0;

    actualUtilizedProcessors = 0.0;
    actualUtilizedMemory = 0.0;
    actualRhoAccelerators = 0.0;

    runningVMs = 0;
    currentCompCapPerProc = 0.0;
    currentCompCapPerAcc = 0.0;
  }
}

void resource::initializeRunningQuantities()
{
  if (alloc) {
    actualUtilizedProcessors = 0.0;
    actualUtilizedMemory = 0.0;
    actualRhoAccelerators = 0.0;
  }
}

void resource::incrementRunningQuantities(const double& uProc, const double& uMem, const double& rAcc)
{
  if (alloc) {
    actualUtilizedProcessors += uProc;
    actualUtilizedMemory += uMem;
    actualRhoAccelerators += rAcc;
  }
}

void resource::compcurrentCompCapPerProc()
{
  double ratio;
  if (alloc) {
    ratio = (max(utilizedProcessors, 1.0) / physicalProcessors);
    currentCompCapPerProc = (computeCapability / physicalProcessors) / ratio;
  }
  // cap = (processor MIPS/physical processors) / (util processors/physical processsors)
}

void resource::compcurrentCompCapPerAcc()
{
  double ratio;
  if (alloc && accelerator) {
    ratio = ceil(max(utilizedProcessors, 1.0) / physicalProcessors);
    currentCompCapPerAcc = (acceleratorComputeCapability / ratio);
  }
}

void resource::deploy(const task& task_)
{
  runningVMs++;
  active = 1;

  availableProcessors -= task_.greqPMNS()[0];
  availableMemory -= task_.greqPMNS()[1];
  availableStorage -= task_.greqPMNS()[3];
  availableAccelerators -= task_.gavAcc()[0];

  utilizedMemory = totalMemory - availableMemory;
  utilizedProcessors = totalProcessors - availableProcessors;
  utilizedStorage = totalStorage - availableStorage;
  utilizedAccelerators = totalAccelerators - availableAccelerators;

  // If the task requests more than 1 VMs, dont allow the resource to be moved for SOSM
  if (task_.getNumberOfVMs() > 1) {
    movable = 0;
  }
}

void resource::unload(const list<task>::iterator& t)
{
  runningVMs--;

  // There are no remaining tasks assigned on the resource
  if (runningVMs == 0) {
    active = 0;

    availableProcessors = totalProcessors;
    availableMemory = totalMemory;
    availableStorage = totalStorage;
    availableAccelerators = totalAccelerators;

    utilizedMemory = 0.0;
    utilizedProcessors = 0.0;
    utilizedStorage = 0.0;
    utilizedAccelerators = 0;

    actualUtilizedProcessors = 0.0;
    actualUtilizedMemory = 0.0;
    actualRhoAccelerators = 0.0;

    currentCompCapPerProc = 0.0;
    currentCompCapPerAcc = 0.0;

    // The resource is allowed to be moved to a different vRM during the Self-Organization phase
    movable = 1;
  } else {
    // There are still tasks assigned on the resource
    availableProcessors += t->greqPMNS()[0];
    availableMemory += t->greqPMNS()[1];
    availableStorage += t->greqPMNS()[3];
    availableAccelerators += t->gavAcc()[0];

    utilizedMemory = totalMemory - availableMemory;
    utilizedProcessors = totalProcessors - availableProcessors;
    utilizedStorage = totalStorage - availableStorage;
    utilizedAccelerators = totalAccelerators - availableAccelerators;

    actualUtilizedProcessors -= t->gcUtilPMNr()[0];
    actualUtilizedMemory -= t->gcUtilPMNr()[1];
    actualRhoAccelerators -= t->gcUtilPMNr()[3];
  }
}

int resource::probe(const double& reqProc, const double& reqMem, const double& reqSto, const int& reqAcc)
{
  int choice = -1;
  if (reqProc <= availableProcessors && reqMem <= availableMemory && reqSto <= availableStorage &&
      reqAcc <= availableAccelerators)
    choice = ID;
  return choice;
}

void resource::print() const
{
  if (alloc) {
    cout << "ID: " << ID << endl;
    cout << "Active: " << active << endl;
    cout << "Movable: " << movable << endl;
    cout << "Type: " << type << endl;

    cout << "Total Processors: " << totalProcessors << endl;
    cout << "Available Processors: " << availableProcessors << endl;
    cout << "Utilized Processors: " << utilizedProcessors << endl;

    cout << "Total Memory: " << totalMemory << endl;
    cout << "Available Memory: " << availableMemory << endl;
    cout << "Utilized Memory: " << utilizedMemory << endl;

    cout << "Total Storage: " << totalStorage << endl;
    cout << "Available Storage: " << availableStorage << endl;
    cout << "Utilized Storage: " << utilizedStorage << endl;

    cout << "Physical Processors: " << physicalProcessors << endl;
    cout << "Physical Memory: " << physicalMemory << endl;
    cout << "Physical Storage: " << physicalStorage << endl;

    cout << "Processors Computational Capability: " << computeCapability << endl;

    cout << "Accelerator Availability: " << accelerator << endl;
    cout << "Total number of accelerators: " << totalAccelerators << endl;
    cout << "Utilized accelerators: " << utilizedAccelerators << endl;
    cout << "Available accelerators: " << availableAccelerators << endl;
    cout << "Accelerator Computational Capability: " << acceleratorComputeCapability << endl;

    cout << "Processors overcommitment ratio: " << overcommitmentProcessors << endl;
    cout << "Memory overcommitment ratio: " << overcommitmentMemory << endl;

    cout << "Actual Processor Utilization: " << actualUtilizedProcessors << endl;
    cout << "Actual Memory Utilization: " << actualUtilizedMemory << endl;
    cout << "Actual rho of Accelerators: " << actualRhoAccelerators << endl;

    cout << "Running VMs: " << runningVMs << endl;
    cout << "Current Processors Computational Capability (per unit): " << currentCompCapPerProc << endl;
    cout << "Current Acc Computational Capability (per unit): " << currentCompCapPerAcc << endl;
  }
}

int resource::galloc() const { return alloc; }
int resource::getActive() const { return active; }
int resource::getMovable() const { return movable; }
int resource::getType() const { return type; }
int resource::gID() const { return ID; }
double resource::getTotalProcessors() const { return totalProcessors; }
double resource::getAvailableProcessors() const { return availableProcessors; }
double resource::getUtilizedProcessors() const { return utilizedProcessors; }
double resource::getTotalMemory() const { return totalMemory; }
double resource::getAvailableMemory() const { return availableMemory; }
double resource::getUtilizedMemory() const { return utilizedMemory; }
double resource::getTotalStorage() const { return totalStorage; }
double resource::getAvailableStorage() const { return availableStorage; }
double resource::getUtilizedStorage() const { return utilizedStorage; }
double resource::getPhysicalProcessors() const { return physicalProcessors; }
double resource::getPhysicalMemory() const { return physicalMemory; }
double resource::getPhysicalStorage() const { return physicalStorage; }
double resource::getComputeCapability() const { return computeCapability; }
int resource::getAccelerator() const { return accelerator; }
int resource::getTotalAccelerators() const { return totalAccelerators; }
int resource::getAvailableAccelerators() const { return availableAccelerators; }
int resource::getUtilizedAccelerators() const { return utilizedAccelerators; }
double resource::getAcceleratorComputeCapability() const { return acceleratorComputeCapability; }
double resource::getOvercommitmentProcessors() const { return overcommitmentProcessors; }
double resource::getOvercommitmentMemory() const { return overcommitmentMemory; }
double resource::getActualUtilizedProcessors() const { return actualUtilizedProcessors; }
double resource::getActualUtilizedMemory() const { return actualUtilizedMemory; }
double resource::getActualRhoAccelerators() const { return actualRhoAccelerators; }
int resource::getRunningVMs() const { return runningVMs; }
double resource::getCurrentCompCapPerProc() const { return currentCompCapPerProc; }
double resource::getCurrentCompCapPerAcc() const { return currentCompCapPerAcc; }

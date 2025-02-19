#include <inputs.h> // for netinputs
#include <netw.h>
#include <task.h> // for task
#include <iostream>

using std::cout;
using std::endl;

netw::netw()
  : alloc(0),
    numberOfTasks(0),
    physicalNetwork(0.0),
    totalNetwork(0.0),
    utilizedNetwork(0.0),
    availableNetwork(0.0),
    actualUtilizedNetwork(0.0),
    overCommitmentNetwork(1.0)
{
}

netw::netw(const netinputs& setup)
{
  alloc = 1;
  numberOfTasks = 0;
  physicalNetwork = setup.netBW;
  totalNetwork = setup.netBW * setup.overCommitmentNetwork;
  utilizedNetwork = 0.0;
  availableNetwork = totalNetwork;
  actualUtilizedNetwork = 0.0;
  overCommitmentNetwork = setup.overCommitmentNetwork;
}

netw::netw(const netw& t)
{
  if (t.galloc()) {
    alloc = t.galloc();
    numberOfTasks = t.getNumberOfTasks();
    physicalNetwork = t.getPhysicalNetwork();
    totalNetwork = t.getTotalNetwork();
    utilizedNetwork = t.getUtilizedNetwork();
    availableNetwork = t.getAvailableNetwork();
    actualUtilizedNetwork = t.getActualUtilizedNetwork();
    overCommitmentNetwork = t.getOverCommitmentNetwork();
  }
}

netw& netw::operator=(const netw& t)
{
  if (this != &t) {
    if (alloc) {
      alloc = 0;
      numberOfTasks = 0;
      physicalNetwork = 0.0;
      totalNetwork = 0.0;
      utilizedNetwork = 0.0;
      availableNetwork = 0.0;
      actualUtilizedNetwork = 0.0;
      overCommitmentNetwork = 1.0;
    }
    alloc = t.galloc();
    if (alloc) {
      numberOfTasks = t.getNumberOfTasks();
      physicalNetwork = t.getPhysicalNetwork();
      totalNetwork = t.getTotalNetwork();
      utilizedNetwork = t.getUtilizedNetwork();
      availableNetwork = t.getAvailableNetwork();
      actualUtilizedNetwork = t.getActualUtilizedNetwork();
      overCommitmentNetwork = t.getOverCommitmentNetwork();
    }
  }
  return *this;
}

netw::~netw()
{
  if (alloc) {
    alloc = 0;
    numberOfTasks = 0;
    physicalNetwork = 0.0;
    totalNetwork = 0.0;
    utilizedNetwork = 0.0;
    availableNetwork = 0.0;
    actualUtilizedNetwork = 0.0;
    overCommitmentNetwork = 1.0;
  }
}

int netw::galloc() const { return alloc; }
double netw::getPhysicalNetwork() const { return physicalNetwork; }
double netw::getTotalNetwork() const { return totalNetwork; }
double netw::getUtilizedNetwork() const { return utilizedNetwork; }
double netw::getAvailableNetwork() const { return availableNetwork; }
double netw::getActualUtilizedNetwork() const { return actualUtilizedNetwork; }
int netw::getNumberOfTasks() const { return numberOfTasks; }
double netw::getOverCommitmentNetwork() const { return overCommitmentNetwork; }
void netw::initializeRunningQuantities()
{
  if (alloc) {
    actualUtilizedNetwork = 0.0;
  }
}

void netw::incrementRunningQuantities(const double& utilizedNetwork)
{
  if (alloc) {
    actualUtilizedNetwork += utilizedNetwork;
  }
}

void netw::sutilizedNetwork(const double& L_utilizedNetwork)
{
  utilizedNetwork = L_utilizedNetwork;
  availableNetwork = (totalNetwork - utilizedNetwork) < 0 ? 0 : (totalNetwork - utilizedNetwork);
}

void netw::print() const
{
  if (alloc) {
    cout << "     Utilized Interconnection Bandwidth: " << utilizedNetwork << " Gbps" << endl;
    cout << "     Available Interconnection Bandwidth: " << availableNetwork << " Gbps" << endl;
    cout << "     Total Interconnection Bandwidth: " << totalNetwork << " Gbps" << endl;
    cout << "     Physical Interconnection Bandwidth: " << physicalNetwork << " Gbps" << endl;
    cout << "     Interconnection Bandwidth Overcommitment ratio: " << overCommitmentNetwork << endl;
  }
}

int netw::probe(const double& requestedNetwork) const
{
  int choice = -1;
  if (requestedNetwork <= availableNetwork) {
    choice = 1;
  }
  return choice;
}

void netw::deploy(const task& task_)
{
  availableNetwork -= task_.greqPMNS()[2];
  utilizedNetwork = totalNetwork - availableNetwork;
  numberOfTasks++;
}

void netw::unload(const double& L_availableNetwork, const double& L_actualUtilizedNetwork, const int& L_numberOfTasks)
{
  availableNetwork += L_availableNetwork;
  utilizedNetwork = totalNetwork - availableNetwork;
  actualUtilizedNetwork -= L_actualUtilizedNetwork;
  numberOfTasks -= L_numberOfTasks;
  if (numberOfTasks == 0) {
    availableNetwork = totalNetwork;
    utilizedNetwork = 0.0;
    actualUtilizedNetwork = 0.0;
  }
}

void netw::unload(list<task>::iterator& t)
{
  availableNetwork += t->greqPMNS()[2];
  utilizedNetwork = totalNetwork - availableNetwork;
  actualUtilizedNetwork -= t->gcUtilPMNr()[2];
  numberOfTasks--;
  if (numberOfTasks == 0) {
    availableNetwork = totalNetwork;
    utilizedNetwork = 0.0;
    actualUtilizedNetwork = 0.0;
  }
}

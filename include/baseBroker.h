#ifndef BASEBROKER_H
#define BASEBROKER_H

class cell;
class netw;
class power;
class resource;
class siminputs;
class stat;
class task;

class baseBroker
{
 public:
  baseBroker();
  // baseBroker(const int& _alloc, const int& _numberOfTypes, int* _types);

  virtual ~baseBroker();

  virtual void deploy(resource** resources, netw* network, stat* stats, task& _task) = 0;

  virtual void init(const cell* clCell);

  virtual void init(const cell* clCell, const siminputs* si) = 0;

  virtual void print() const = 0;

  virtual void timestep(const cell* clCell) = 0;

  virtual void updateStateInfo(const cell* clCell, const double& tstep) = 0;

 protected:
  int alloc;
  int numberOfTypes;
  int* types;
  int* numberOfResourcesPerType;

  double availableNetwork;
  double totalNetwork;

  int galloc() const;
  int* getTypes() const;
  int getNumberOfTypes() const;
  int* getNumberOfResourcesPerType() const;
  double getAvailableNetwork() const;
  double getTotalNetwork() const;
};

#endif

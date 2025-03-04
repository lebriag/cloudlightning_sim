#ifndef VRM_H
#define VRM_H
#include <list>

class netw;
class resource;
class stat;
class task;

using std::list;

class vRM
{
 private:
  int alloc;
  int numberOfResources;
  int numberOfFunctions;
  int optNumOfRes;
  double pollIntervalvRM;
  list<task>* queue;
  list<resource*>* res;
  double* Fs;
  double* Ws;
  double* availableProcessors;
  double* totalProcessors;
  double* availableMemory;
  double* totalMemory;
  double* availableAccelerators;
  double* totalAccelerators;
  double* availableStorage;
  double* totalStorage;
  double* sPMSA;
  double C, P, Pi;
  double SI;
  int dep_strategy;

 public:
  vRM();

  ~vRM();

  vRM(const int& start, const int& end, const int& type, resource** resources, const double& L_pollIntervalvRM,
      const double& L_C, const double& L_P, const double& L_Pi, const int& L_optNumOfRes,
      const int& L_numberOfFunctions, const double* L_Ws, const int& L_dep_strategy);

  vRM(const vRM& t);

  vRM& operator=(const vRM& t);

  /// Computes the assesment functions
  void computeFs();

  /// Computes the suitability index
  void computeSI();

  /// Updates the state information of the cell
  void updateStateInfo(const double& tstep);

  /// Removes the list of resources that satisfy the requested units from a vRM and adds them to the ores list
  /// \param ores The resources list to store the resources removed from the vRM
  /// \param remProc The number of requested processor cores
  /// \param remMem The size of requested memory
  /// \param remSto The size of requested storage
  /// \param remAcc The number of requested accelerators
  void obtainresources(list<resource*>& ores, double& remProc, double& remMem, double& remSto, double& remAcc);

  /// Adds the resources contained in the ores list to the vRM
  /// \param The resources list to be attached to the vRM
  void attachresources(list<resource*>& ores);

  double assessfuncs(const int& choice);

  double deassessmentFunctions(const double& dNu, const double& dNmem, const int& choice);

  /// Returns success if the task's processes, memory, storage and accelerator are less or equal than the vRM's
  int probe(const double& Proc, const double& Mem, const double& Sto, const int& Acc);

  /// Defines the strategies used to position VMs on resources
  int deploy_strategy(list<resource*>::iterator* it, int* IDs, const int& nVMs, const double& Proc, const double& Mem,
                      const double& Sto, const int& Acc);

  /// Deploys the tasks
  void deploy(resource** resources, netw* network, stat* stats, task& task_);

  void enque(const task& task_);

  void print();

  int galloc() const;
  int getNumberOfResources() const;
  int goptNumOfRes() const;
  double gpollIntervalvRM() const;
  list<task>* gqueue() const;
  list<resource*>* gres() const;
  double* getAvailableProcessors() const;
  double* getTotalProcessors() const;
  double* getAvailableMemory() const;
  double* getTotalMemory() const;
  double* getAvailableAccelerators() const;
  double* getTotalAccelerators() const;
  double* getAvailableStorage() const;
  double* getTotalStorage() const;
  double* gsPMSA() const;
  double* gFs() const;
  double* gWs() const;
  double gC() const;
  double gP() const;
  double gPi() const;
  double gSI() const;
  int gdep_strategy() const;
  int getNumberOfFunctions() const;
};

#endif

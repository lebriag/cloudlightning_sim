#ifndef SOSMBROKER_H
#define SOSMBROKER_H

#include <baseBroker.h>
#include <list>

class brinputs;
class cell;
class pRouter;
class pSwitch;
class siminputs;
class vRM;

using std::list;

class sosmBroker : public baseBroker
{
 private:
  int numberOfvRMs;
  int numberOfpSwitches;
  int numberOfpRouters;

  double pollIntervalCellM;
  double pollIntervalpRouter;
  double pollIntervalpSwitch;
  double pollIntervalvRM;

  double** sPMSA;
  double* SIs;
  double *Cs, *Ps, *Pis;
  double* Ws;

  int numberOfFunctions;

  list<vRM>** vRMs;
  list<pSwitch>** pSwitches;
  list<pRouter>** pRouters;

 public:
  sosmBroker();

  /*sosmBroker(const int& L_numberOfTypes, const int* L_types, const int* L_numberOfResourcesPerType, resource**
     resources,
         power* powerComp, netw* network, const brinputs& binp);*/

  sosmBroker(const sosmBroker& t);

  sosmBroker& operator=(const sosmBroker& t);

  ~sosmBroker();

  void init(const cell* clCell, const siminputs* si);

  void print() const;

  /// Calculates the de-assessment functions
  double deassessmentFunctions(const double& dNu, const double& dNmem, const int& choice, const int& type);

  /// Updates the state information of the cell, by calling recursively the pRouter::updateStateInfo method
  void updateStateInfo(const cell* clCell, const double& tstep);

  /// Deploys the tasks to the appropriate vRMs, by calling recursively the pRouter::deploy method
  void deploy(resource** resources, netw* network, stat* stats, task& _task);

  /// Performs the simulation phase
  void timestep(const cell* clCell);

  int getNumberOfvRMs() const;
  int getNumberOfpSwitches() const;
  int getNumberOfpRouters() const;

  double gpollIntervalCellM() const;
  double gpollIntervalpRouter() const;
  double gpollIntervalpSwitch() const;
  double gpollIntervalvRM() const;

  double** gsPMSA() const;
  double* gSIs() const;
  double* gCs() const;
  double* gPs() const;
  double* gPis() const;
  double* gWs() const;

  int getNumberOfFunctions() const;

  list<vRM>** getvRMs() const;
  list<pSwitch>** gpSwitches() const;
  list<pRouter>** gpRouters() const;
};

#endif

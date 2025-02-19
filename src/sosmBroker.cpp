#include <cell.h>
#include <inputs.h>
#include <netw.h> // for netw
#include <omp.h>
#include <pRouter.h>
#include <pSwitch.h>
#include <power.h>    // for power
#include <resource.h> // for resource
#include <sosmBroker.h>
#include <stat.h> // for stat
#include <task.h> // for task
#include <vRM.h>
#include <algorithm> // for min
#include <cstdlib>   // for atoi, getenv
#include <iostream>

using std::cout;
using std::endl;
using std::min;

sosmBroker::sosmBroker()
  : baseBroker(),
    numberOfvRMs(0),
    numberOfpSwitches(0),
    numberOfpRouters(0),
    pollIntervalCellM(0.0),
    pollIntervalpRouter(0.0),
    pollIntervalpSwitch(0.0),
    pollIntervalvRM(0.0),
    sPMSA(nullptr),
    vRMs(nullptr),
    SIs(nullptr),
    pSwitches(nullptr),
    pRouters(nullptr),
    Cs(nullptr),
    Ps(nullptr),
    Pis(nullptr),
    Ws(nullptr),
    numberOfFunctions(0)
{
}

void sosmBroker::init(const cell* clCell, const siminputs* si)
{
  baseBroker::init(clCell);

  // Copy the interval time for polling from the BrokerData configuration file
  pollIntervalCellM = si->cinp->binp[0].pollIntervalCellM;
  pollIntervalpRouter = si->cinp->binp[0].pollIntervalpRouter;
  pollIntervalpSwitch = si->cinp->binp[0].pollIntervalpSwitch;
  pollIntervalvRM = si->cinp->binp[0].pollIntervalvRM;

  double* tempC = new double[numberOfTypes];
  double* tempP = new double[numberOfTypes];
  double* tempPi = new double[numberOfTypes];

  // Calculate the values for C, P and Pi to be used for the assessment functions
  for (int i = 0; i < numberOfTypes; i++) {
    tempC[i] = clCell->getResources()[i][0].getComputeCapability() +
               clCell->getResources()[i][0].getAcceleratorComputeCapability();
    double oz = 1.0;
    int active = clCell->getResources()[i][0].getActive();
    int totalAccelerators = clCell->getResources()[i][0].getTotalAccelerators();
    tempP[i] = clCell->getPowerConsumption()[i].consumption(oz, oz, active, totalAccelerators);
    oz = 0.0;
    tempPi[i] = clCell->getPowerConsumption()[i].consumption(oz, oz, active, totalAccelerators);
  }

  int tminC = 0;
  int tminP = 0;
  for (int i = 1; i < numberOfTypes; i++) {
    if (tempC[tminC] > tempC[i])
      tminC = i;
    if (tempP[tminP] > tempP[i])
      tminP = i;
  }
  double minC = tempC[tminC];
  double minP = tempP[tminP];

  for (int i = 0; i < numberOfTypes; i++) {
    tempPi[i] = tempPi[i] / tempP[i];
  }

  for (int i = 0; i < numberOfTypes; i++) {
    tempC[i] /= minC;
    tempP[i] /= minP;
  }

  // Set the number of assesment functions
  numberOfFunctions = si->cinp->binp[0].numberOfFunctions;

  // Allocate an array of weights for the assesment functions
  Ws = new double[numberOfFunctions];

  // Set the individual weights
  for (int i = 0; i < numberOfFunctions; i++) {
    Ws[i] = si->cinp->binp[0].Ws[i];
  }

  // Set the number of pRouters as equal to the number of hardware types
  numberOfpRouters = numberOfTypes;
  numberOfvRMs = 0;
  numberOfpSwitches = 0;
  vRMs = new list<vRM>*[numberOfTypes];
  pSwitches = new list<pSwitch>*[numberOfTypes];
  pRouters = new list<pRouter>*[numberOfTypes];

  for (int i = 0; i < numberOfTypes; i++) {
    vRMs[i] = new list<vRM>[1];
    pSwitches[i] = new list<pSwitch>[1];
    pRouters[i] = new list<pRouter>[1];
  }

  // Create hierarchical the SOSM vRMs
  for (int i = 0; i < numberOfTypes; i++) {
    int temp = numberOfResourcesPerType[i] / si->cinp->binp[0].initResPervRM;
    for (int j = 0; j < temp; j++) {
      vRMs[i]->push_back(vRM(j * si->cinp->binp[0].initResPervRM, (j + 1) * si->cinp->binp[0].initResPervRM, i,
                             clCell->getResources(), pollIntervalvRM, tempC[i], tempP[i], tempPi[i],
                             si->cinp->binp[0].initResPervRM, si->cinp->binp[0].numberOfFunctions, si->cinp->binp[0].Ws,
                             si->cinp->binp[0].vRMdeploystrategy));
    }
    numberOfvRMs += temp;
    if (temp * si->cinp->binp[0].initResPervRM != numberOfResourcesPerType[i]) {
      vRMs[i]->push_back(vRM(temp * si->cinp->binp[0].initResPervRM, numberOfResourcesPerType[i], i,
                             clCell->getResources(), pollIntervalvRM, tempC[i], tempP[i], tempPi[i],
                             si->cinp->binp[0].initResPervRM, si->cinp->binp[0].numberOfFunctions, si->cinp->binp[0].Ws,
                             si->cinp->binp[0].vRMdeploystrategy));
      numberOfvRMs++;
    }
  }

  // Create hierarchical the SOSM pSwitches
  for (int i = 0; i < numberOfTypes; i++) {
    int temp = ((int)vRMs[i]->size()) / si->cinp->binp[0].initvRMPerpSwitch;
    for (int j = 0; j < temp; j++) {
      pSwitches[i]->push_back(pSwitch(
        j * si->cinp->binp[0].initvRMPerpSwitch, (j + 1) * si->cinp->binp[0].initvRMPerpSwitch, i, vRMs,
        pollIntervalpSwitch, tempC[i], tempP[i], tempPi[i], si->cinp->binp[0].numberOfFunctions, si->cinp->binp[0].Ws));
    }
    numberOfpSwitches += temp;
    if (temp * si->cinp->binp[0].initvRMPerpSwitch != ((int)vRMs[i]->size())) {
      pSwitches[i]->push_back(pSwitch(temp * si->cinp->binp[0].initvRMPerpSwitch, ((int)vRMs[i]->size()), i, vRMs,
                                      pollIntervalpSwitch, tempC[i], tempP[i], tempPi[i],
                                      si->cinp->binp[0].numberOfFunctions, si->cinp->binp[0].Ws));
      numberOfpSwitches++;
    }
  }

  // Create hierarchical the SOSM pRouters
  for (int i = 0; i < numberOfTypes; i++) {
    pRouters[i]->push_back(pRouter(0, ((int)pSwitches[i]->size()), i, pSwitches, pollIntervalpRouter, tempC[i],
                                   tempP[i], tempPi[i], si->cinp->binp[0].numberOfFunctions, si->cinp->binp[0].Ws));
  }

  // Allocate the sPMSA array
  sPMSA = new double*[numberOfTypes];
  for (int i = 0; i < numberOfTypes; i++) {
    sPMSA[i] = new double[8];
  }

  for (int i = 0; i < numberOfTypes; i++) {
    for (list<pRouter>::iterator it = pRouters[i]->begin(); it != pRouters[i]->end(); it++) {
      for (int j = 0; j < 8; j++) {
        sPMSA[i][j] = it->gsPMSA()[j];
      }
    }
  }
  availableNetwork = clCell->getNetwork()->getAvailableNetwork();
  totalNetwork = clCell->getNetwork()->getTotalNetwork();
  SIs = new double[numberOfTypes];
  for (int i = 0; i < numberOfTypes; i++) {
    SIs[i] = 0.0;
  }
  Cs = tempC;
  tempC = nullptr;
  Ps = tempP;
  tempP = nullptr;
  Pis = tempPi;
  tempPi = nullptr;
}
/*
sosmBroker::sosmBroker(const int &L_numberOfTypes, const int *L_types, const int *L_numberOfResourcesPerType, resource
**resources,
power *powerComp, netw *network, const brinputs &binp)
{
  int i=0,j=0;
  double *tempC=nullptr,*tempP=nullptr,*tempPi=nullptr;
  double minP,minC;
  int tminC,tminP;
  alloc=1;
  numberOfTypes=L_numberOfTypes;
  types=new int[numberOfTypes];
  numberOfResourcesPerType=new int[numberOfTypes];
  for(i=0;i<numberOfTypes;i++)
  {
    types[i]=L_types[i];
    numberOfResourcesPerType[i]=L_numberOfResourcesPerType[i];
  }
  pollIntervalCellM=binp.pollIntervalCellM;
  pollIntervalpRouter=binp.pollIntervalpRouter;
  pollIntervalpSwitch=binp.pollIntervalpSwitch;
  pollIntervalvRM=binp.pollIntervalvRM;

  tempC=new double[numberOfTypes];
  tempP=new double[numberOfTypes];
  tempPi=new double[numberOfTypes];
  for(i=0;i<numberOfTypes;i++)
  {
    tempC[i]=resources[i][0].getComputeCapability()+resources[i][0].getAcceleratorComputeCapability();
    double oz=1.0;
    int active=resources[i][0].getActive();
    int totalAccelerators=resources[i][0].getTotalAccelerators();
    tempP[i]=powerComp[i].consumption(oz,oz,active,totalAccelerators);
    oz=0.0;
    tempPi[i]=powerComp[i].consumption(oz,oz,active,totalAccelerators);
  }
  tminC=0;
  tminP=0;
  for(i=1;i<numberOfTypes;i++)
  {
    if(tempC[tminC]>tempC[i])
      tminC=i;
    if(tempP[tminP]>tempP[i])
      tminP=i;
  }
  minC=tempC[tminC];
  minP=tempP[tminP];

  for(i=0;i<numberOfTypes;i++)
  {
    tempPi[i]=tempPi[i]/tempP[i];
  }

  for(i=0;i<numberOfTypes;i++)
  {
    tempC[i]/=minC;
    tempP[i]/=minP;
  }

  numberOfFunctions=binp.numberOfFunctions;
  Ws=new double[numberOfFunctions];
  for(i=0;i<numberOfFunctions;i++)
    Ws[i]=binp.Ws[i];

  numberOfpRouters=numberOfTypes;
  numberOfvRMs=0;
  numberOfpSwitches=0;
  vRMs=new list<vRM>*[numberOfTypes];
  pSwitches=new list<pSwitch>*[numberOfTypes];
  pRouters=new list<pRouter>*[numberOfTypes];
  for(i=0;i<numberOfTypes;i++)
  {
    vRMs[i]=new list<vRM>[1];
    pSwitches[i]=new list<pSwitch>[1];
    pRouters[i]=new list<pRouter>[1];
  }
  for(i=0;i<numberOfTypes;i++)
  {
    int temp=numberOfResourcesPerType[i]/binp.initResPervRM;
    for(j=0;j<temp;j++)
    {
      vRMs[i]->push_back(vRM(j*binp.initResPervRM,(j+1)*binp.initResPervRM,i,resources,pollIntervalvRM,tempC[i],tempP[i],tempPi[i],binp.initResPervRM,binp.numberOfFunctions,binp.Ws,binp.vRMdeploystrategy));
    }
    numberOfvRMs+=temp;
    if (temp*binp.initResPervRM!=numberOfResourcesPerType[i])
    {
      vRMs[i]->push_back(vRM(temp*binp.initResPervRM,numberOfResourcesPerType[i],i,resources,pollIntervalvRM,tempC[i],tempP[i],tempPi[i],binp.initResPervRM,binp.numberOfFunctions,binp.Ws,binp.vRMdeploystrategy));
      numberOfvRMs++;
    }
  }

  for(i=0;i<numberOfTypes;i++)
  {
    int temp=((int)vRMs[i]->size())/binp.initvRMPerpSwitch;
    for(j=0;j<temp;j++)
    {
      pSwitches[i]->push_back(pSwitch(j*binp.initvRMPerpSwitch,(j+1)*binp.initvRMPerpSwitch,i,vRMs,pollIntervalpSwitch,tempC[i],tempP[i],tempPi[i],binp.numberOfFunctions,binp.Ws));
    }
    numberOfpSwitches+=temp;
    if(temp*binp.initvRMPerpSwitch!=((int)vRMs[i]->size()))
    {
      pSwitches[i]->push_back(pSwitch(temp*binp.initvRMPerpSwitch,((int)vRMs[i]->size()),i,vRMs,pollIntervalpSwitch,tempC[i],tempP[i],tempPi[i],binp.numberOfFunctions,binp.Ws));
      numberOfpSwitches++;
    }
  }

  for(i=0;i<numberOfTypes;i++)
  {
    pRouters[i]->push_back(pRouter(0,((int)pSwitches[i]->size()),i,pSwitches,pollIntervalpRouter,tempC[i],tempP[i],tempPi[i],binp.numberOfFunctions,binp.Ws));
  }
  sPMSA=new double*[numberOfTypes];
  for(i=0;i<numberOfTypes;i++)
    sPMSA[i]=new double[8];

  for(i=0;i<numberOfTypes;i++)
  {
    for(list<pRouter>::iterator it = pRouters[i]->begin(); it != pRouters[i]->end(); it++)
    {
      for(j=0;j<8;j++)
        sPMSA[i][j]=it->gsPMSA()[j];
    }
  }
  availableNetwork=network->getAvailableNetwork();
  totalNetwork=network->getTotalNetwork();
  SIs=new double[numberOfTypes];
  for(i=0;i<numberOfTypes;i++)
    SIs[i]=0.0;
  Cs=tempC;
  tempC=nullptr;
  Ps=tempP;
  tempP=nullptr;
  Pis=tempPi;
  tempPi=nullptr;
}
*/
sosmBroker::sosmBroker(const sosmBroker& t)
{
  int i, j;
  if (t.galloc()) {
    alloc = 1;
    numberOfTypes = t.getNumberOfTypes();
    pollIntervalCellM = t.gpollIntervalCellM();
    pollIntervalpRouter = t.gpollIntervalpRouter();
    pollIntervalpSwitch = t.gpollIntervalpSwitch();
    pollIntervalvRM = t.gpollIntervalvRM();
    numberOfpRouters = t.getNumberOfpRouters();
    numberOfpSwitches = t.getNumberOfpSwitches();
    numberOfvRMs = t.getNumberOfvRMs();
    types = new int[numberOfTypes];
    numberOfResourcesPerType = new int[numberOfTypes];

    for (i = 0; i < numberOfTypes; i++) {
      types[i] = t.getTypes()[i];
      numberOfResourcesPerType[i] = t.getNumberOfResourcesPerType()[i];
    }

    vRMs = new list<vRM>*[numberOfTypes];
    pSwitches = new list<pSwitch>*[numberOfTypes];
    pRouters = new list<pRouter>*[numberOfTypes];

    for (i = 0; i < numberOfTypes; i++) {
      vRMs[i] = new list<vRM>[1];
      pSwitches[i] = new list<pSwitch>[1];
      pRouters[i] = new list<pRouter>[1];
    }

    for (i = 0; i < numberOfTypes; i++) {
      for (list<vRM>::iterator it = t.getvRMs()[i]->begin(); it != t.getvRMs()[i]->end(); it++) {
        vRMs[i]->push_back(*it);
      }
    }
    for (i = 0; i < numberOfTypes; i++) {
      pSwitches[i][0] = t.gpSwitches()[i][0];
    }

    for (list<pSwitch>::iterator it = t.gpSwitches()[i]->begin(); it != t.gpSwitches()[i]->end(); it++) {
      pSwitches[i]->push_back(*it);
    }

    for (i = 0; i < numberOfTypes; i++) {
      pRouters[i][0] = t.gpRouters()[i][0];
    }

    for (list<pRouter>::iterator it = t.gpRouters()[i]->begin(); it != t.gpRouters()[i]->end(); it++) {
      pRouters[i]->push_back(*it);
    }

    sPMSA = new double*[numberOfTypes];

    for (i = 0; i < numberOfTypes; i++) {
      sPMSA[i] = new double[8];
    }

    for (i = 0; i < numberOfTypes; i++) {
      for (j = 0; j < 8; j++) {
        sPMSA[i][j] = t.gsPMSA()[i][j];
      }
    }

    availableNetwork = t.getAvailableNetwork();
    totalNetwork = t.getTotalNetwork();
    SIs = new double[numberOfTypes];

    for (i = 0; i < numberOfTypes; i++) {
      SIs[i] = t.gSIs()[i];
    }

    Cs = new double[numberOfTypes];
    Ps = new double[numberOfTypes];
    Pis = new double[numberOfTypes];

    for (i = 0; i < numberOfTypes; i++) {
      Cs[i] = t.gCs()[i];
      Ps[i] = t.gPs()[i];
      Pis[i] = t.gPis()[i];
    }

    numberOfFunctions = t.getNumberOfFunctions();
    Ws = new double[numberOfFunctions];

    for (i = 0; i < numberOfFunctions; i++) {
      Ws[i] = t.gWs()[i];
    }
  }
}

sosmBroker& sosmBroker::operator=(const sosmBroker& t)
{
  int i, j;
  if (this != &t) {
    if (alloc) {
      alloc = 0;
      delete[] types;
      delete[] numberOfResourcesPerType;
      types = nullptr;
      numberOfResourcesPerType = nullptr;
      for (i = 0; i < numberOfTypes; i++) {
        vRMs[i]->clear();
        pSwitches[i]->clear();
        pRouters[i]->clear();
        delete[] vRMs[i];
        delete[] pSwitches[i];
        delete[] pRouters[i];
        delete[] sPMSA[i];
      }
      delete[] pRouters;
      delete[] pSwitches;
      delete[] vRMs;
      delete[] sPMSA;
      delete[] SIs;
      SIs = nullptr;
      sPMSA = nullptr;
      vRMs = nullptr;
      pSwitches = nullptr;
      pRouters = nullptr;
      numberOfvRMs = 0;
      numberOfpSwitches = 0;
      numberOfpRouters = 0;
      pollIntervalCellM = 0.0;
      pollIntervalpRouter = 0.0;
      pollIntervalpSwitch = 0.0;
      pollIntervalvRM = 0.0;
      numberOfTypes = 0;
      availableNetwork = 0.0;
      totalNetwork = 0.0;
      delete[] Cs;
      delete[] Ps;
      delete[] Pis;
      Cs = nullptr;
      Ps = nullptr;
      Pis = nullptr;
      numberOfFunctions = 0;
      delete[] Ws;
      Ws = nullptr;
    }
    alloc = t.galloc();
    if (alloc) {
      numberOfTypes = t.getNumberOfTypes();
      pollIntervalCellM = t.gpollIntervalCellM();
      pollIntervalpRouter = t.gpollIntervalpRouter();
      pollIntervalpSwitch = t.gpollIntervalpSwitch();
      pollIntervalvRM = t.gpollIntervalvRM();
      numberOfpRouters = t.getNumberOfpRouters();
      numberOfpSwitches = t.getNumberOfpSwitches();
      numberOfvRMs = t.getNumberOfvRMs();
      types = new int[numberOfTypes];
      numberOfResourcesPerType = new int[numberOfTypes];

      for (i = 0; i < numberOfTypes; i++) {
        types[i] = t.getTypes()[i];
        numberOfResourcesPerType[i] = t.getNumberOfResourcesPerType()[i];
      }

      vRMs = new list<vRM>*[numberOfTypes];
      pSwitches = new list<pSwitch>*[numberOfTypes];
      pRouters = new list<pRouter>*[numberOfTypes];

      for (i = 0; i < numberOfTypes; i++) {
        vRMs[i] = new list<vRM>[1];
        pSwitches[i] = new list<pSwitch>[1];
        pRouters[i] = new list<pRouter>[1];
      }

      for (i = 0; i < numberOfTypes; i++) {
        for (list<vRM>::iterator it = t.getvRMs()[i]->begin(); it != t.getvRMs()[i]->end(); it++) {
          vRMs[i]->push_back(*it);
        }
      }

      for (i = 0; i < numberOfTypes; i++) {
        for (list<pSwitch>::iterator it = t.gpSwitches()[i]->begin(); it != t.gpSwitches()[i]->end(); it++) {
          pSwitches[i]->push_back(*it);
        }
      }

      for (i = 0; i < numberOfTypes; i++) {
        for (list<pRouter>::iterator it = t.gpRouters()[i]->begin(); it != t.gpRouters()[i]->end(); it++) {
          pRouters[i]->push_back(*it);
        }
      }

      sPMSA = new double*[numberOfTypes];

      for (i = 0; i < numberOfTypes; i++) {
        sPMSA[i] = new double[8];
      }

      for (i = 0; i < numberOfTypes; i++) {
        for (j = 0; j < 8; j++) {
          sPMSA[i][j] = t.gsPMSA()[i][j];
        }
      }

      availableNetwork = t.getAvailableNetwork();
      totalNetwork = t.getTotalNetwork();
      SIs = new double[numberOfTypes];

      for (i = 0; i < numberOfTypes; i++) {
        SIs[i] = t.gSIs()[i];
      }

      Cs = new double[numberOfTypes];
      Ps = new double[numberOfTypes];
      Pis = new double[numberOfTypes];

      for (i = 0; i < numberOfTypes; i++) {
        Cs[i] = t.gCs()[i];
        Ps[i] = t.gPs()[i];
        Pis[i] = t.gPis()[i];
      }

      numberOfFunctions = t.getNumberOfFunctions();
      Ws = new double[numberOfFunctions];

      for (i = 0; i < numberOfFunctions; i++) {
        Ws[i] = t.gWs()[i];
      }
    }
  }
  return *this;
}

sosmBroker::~sosmBroker()
{
  if (alloc) {
    alloc = 0;
    delete[] types;
    delete[] numberOfResourcesPerType;
    types = nullptr;
    numberOfResourcesPerType = nullptr;

    for (int i = 0; i < numberOfTypes; i++) {
      vRMs[i]->clear();
      pSwitches[i]->clear();
      pRouters[i]->clear();
      delete[] vRMs[i];
      delete[] pSwitches[i];
      delete[] pRouters[i];
      delete[] sPMSA[i];
    }

    delete[] vRMs;
    delete[] pSwitches;
    delete[] pRouters;
    delete[] sPMSA;
    delete[] SIs;
    SIs = nullptr;
    vRMs = nullptr;
    pSwitches = nullptr;
    pRouters = nullptr;
    numberOfvRMs = 0;
    numberOfpSwitches = 0;
    numberOfpRouters = 0;
    pollIntervalCellM = 0.0;
    pollIntervalpRouter = 0.0;
    pollIntervalpSwitch = 0.0;
    pollIntervalvRM = 0.0;
    numberOfTypes = 0;
    availableNetwork = 0.0;
    totalNetwork = 0.0;
    delete[] Cs;
    delete[] Ps;
    delete[] Pis;
    Cs = nullptr;
    Ps = nullptr;
    Pis = nullptr;
    delete[] Ws;
    Ws = nullptr;
    numberOfFunctions = 0;
  }
}

void sosmBroker::print() const
{
  if (alloc) {
    for (int i = 0; i < numberOfTypes; i++) {
      cout << "No of Type: " << i << endl;
      cout << "     Number of available processing units: " << sPMSA[i][0] << endl;
      cout << "     Number of total processing units: " << sPMSA[i][1] << endl;
      cout << "     Number of available memory: " << sPMSA[i][2] << endl;
      cout << "     Number of total memory: " << sPMSA[i][3] << endl;
      cout << "     Number of available storage: " << sPMSA[i][4] << endl;
      cout << "     Number of total storage: " << sPMSA[i][5] << endl;
      cout << "     Number of available accelerators: " << sPMSA[i][6] << endl;
      cout << "     Number of total accelerators: " << sPMSA[i][7] << endl;
    }
    cout << "Available network bandwidth: " << availableNetwork << endl;
    cout << "Total network bandwidth: " << totalNetwork << endl;
  }
}

double sosmBroker::deassessmentFunctions(const double& dNu, const double& dNmem, const int& choice, const int& type)
{
  if (sPMSA[type][7] > 0) {
    switch (choice) {
      case 0:
        return (dNu * Cs[type] / sPMSA[type][7]);
        break;
      case 1:
        return (dNmem / sPMSA[type][3]);
        break;
      case 2:
        return (dNu * Pis[type] * Ps[type] * sPMSA[type][7]) /
               ((Ps[type] * (sPMSA[type][7] - sPMSA[type][6]) + Pis[type] * sPMSA[type][6]) *
                (Ps[type] * (sPMSA[type][7] - sPMSA[type][6]) + Pis[type] * sPMSA[type][6]));
        break;
      case 3:
        return (0.2 * dNu / (sPMSA[type][7]));
        break;
      default:
        return 0.0;
    }
  } else {
    switch (choice) {
      case 0:
        return (dNu * Cs[type] / sPMSA[type][1]);
        break;
      case 1:
        return (dNmem / sPMSA[type][3]);
        break;
      case 2:
        return (dNu * Pis[type] * Ps[type] * sPMSA[type][1]) /
               ((Ps[type] * (sPMSA[type][1] - sPMSA[type][0]) + Pis[type] * sPMSA[type][0]) *
                (Ps[type] * (sPMSA[type][1] - sPMSA[type][0]) + Pis[type] * sPMSA[type][0]));
        break;
      case 3:
        return (0.2 * dNu / (sPMSA[type][1]));
        break;
      default:
        return 0.0;
    }
  }
}

void sosmBroker::updateStateInfo(const cell* clCell, const double& tstep)
{
  int i, j, k;
  list<pRouter>::iterator ittt;
  list<pSwitch>::iterator itt;
  list<vRM>::iterator it;
  if (alloc) {
    // On every time interval
    if (((int)tstep % (int)pollIntervalvRM) == 0) {
      int omp_thr, len;
      // For each hardware type
      for (i = 0; i < numberOfTypes; i++) {
        // How many vRMs correspond to each hardware type
        len = (int)vRMs[i]->size();
        omp_thr = atoi(getenv("OMP_NUM_THREADS"));
#pragma omp parallel default(shared) private(j, k, it) num_threads(omp_thr)
        {
          int tid = omp_get_thread_num();
          it = vRMs[i]->begin();
          for (j = 0; j < tid; j++) {
            it++;
          }
          j = tid;
          while (j < len) {
            // Update state info on each vRM in parallel
            it->updateStateInfo(tstep);
            j += omp_thr;
            for (k = 0; k < omp_thr; k++) {
              it++;
            }
          }
#pragma omp barrier
        }
      }
    }
    if (((int)tstep % (int)pollIntervalpSwitch) == 0) {
      int omp_thr, len;
      for (i = 0; i < numberOfTypes; i++) {
        len = (int)pSwitches[i]->size();
        omp_thr = atoi(getenv("OMP_NUM_THREADS"));
#pragma omp parallel default(shared) private(j, k, itt) num_threads(omp_thr)
        {
          int tid = omp_get_thread_num();
          itt = pSwitches[i]->begin();
          for (j = 0; j < tid; j++)
            itt++;
          j = tid;
          while (j < len) {
            itt->updateStateInfo(tstep);
            for (k = 0; k < omp_thr; k++)
              itt++;
            j += omp_thr;
          }
#pragma omp barrier
        }
      }
    }
    if (((int)tstep % (int)pollIntervalpRouter) == 0) {
      int omp_thr = atoi(getenv("OMP_NUM_THREADS"));
#pragma omp parallel for default(shared) private(i, ittt) num_threads(omp_thr) schedule(static, 1)
      for (i = 0; i < numberOfTypes; i++) {
        ittt = pRouters[i]->begin();
        ittt->updateStateInfo(tstep);
      }
    }
    if (((int)tstep % (int)pollIntervalCellM) == 0) {
      for (i = 0; i < numberOfTypes; i++) {
        for (list<pRouter>::iterator itttt = pRouters[i]->begin(); itttt != pRouters[i]->end(); itttt++) {
          for (j = 0; j < 8; j++)
            sPMSA[i][j] = itttt->gsPMSA()[j];
          SIs[i] = itttt->gSI();
        }
      }
    }
    availableNetwork = clCell[0].getNetwork()->getAvailableNetwork();
    totalNetwork = clCell[0].getNetwork()->getTotalNetwork();
  }
}

void sosmBroker::deploy(resource** resources, netw* network, stat* stats, task& _task)
{
  int *rem = new int[numberOfTypes];
  int *rem2 = new int[_task.getNumberOfAvailableImplementations()];
  int count = 0;

  // Determine how many (count) of the pRouter hardware types match the tasks' requested type, and store the types to
  // the rem and rem2 arrays
  for (int j = 0; j < _task.getNumberOfAvailableImplementations(); j++) {
    for (int i = 0; i < numberOfTypes; i++) {
      if (types[i] == _task.getAvailableImplementations()[j]) {
        rem[count] = i;
        rem2[count] = j;
        count++;
        break;
      }
    }
  }

  if (availableNetwork < _task.greqPMNS()[2]) {
    stats[rem[0]].rejectedTasks++;
    delete[] rem;
    delete[] rem2;
    return;
  }
  availableNetwork -= _task.greqPMNS()[2];

  double maxSI = 0.0;
  int type = -1;

  // Lets find the pRouter out of the pRouters with a compatible Hardware Type with the highest Suitability Index
  for (int i = 0; i < count; i++) {
    list<pRouter>::iterator itt = pRouters[rem[i]]->begin();
    // Compare the current found maximum SI to the SI of each pRouter with a matching Hardware Type and double check
    // that the pRouter offers the requested PMNS resources
    if (maxSI < SIs[rem[i]] && _task.getNumberOfVMs() * _task.greqPMNS()[0] <= sPMSA[rem[i]][0] &&
        _task.getNumberOfVMs() * _task.greqPMNS()[1] <= sPMSA[rem[i]][2] &&
        _task.getNumberOfVMs() * _task.greqPMNS()[3] <= sPMSA[rem[i]][4] &&
        _task.getNumberOfVMs() * _task.gavAcc()[rem2[i]] <= sPMSA[rem[i]][6]) {
      if (itt->probe(_task.getNumberOfVMs() * _task.greqPMNS()[0], _task.getNumberOfVMs() * _task.greqPMNS()[1],
                     _task.getNumberOfVMs() * _task.greqPMNS()[3], _task.getNumberOfVMs() * _task.gavAcc()[rem2[i]]) != -1) {
        maxSI = SIs[rem[i]];
        type = i;
      }
    }
  }
  if (type == -1) {
    // Reject tasks if type is still -1
    stats[rem[0]].rejectedTasks++;
    return;
  }
  _task.reduceImpl(&rem2[type]);
  _task.remapType(&rem[type], 1);
  type = rem[type];
  sPMSA[type][0] -= _task.getNumberOfVMs() * _task.greqPMNS()[0];
  sPMSA[type][2] -= _task.getNumberOfVMs() * _task.greqPMNS()[1];
  sPMSA[type][4] -= _task.getNumberOfVMs() * _task.greqPMNS()[3];
  sPMSA[type][6] -= _task.getNumberOfVMs() * _task.gavAcc()[0];

  // Update SI for the chosen pRouter
  for (int i = 0; i < 4; i++) {
    SIs[type] += Ws[i] * deassessmentFunctions(-_task.getNumberOfVMs() * _task.greqPMNS()[0],
                                               -_task.getNumberOfVMs() * _task.greqPMNS()[1], i, type);
  }

  list<pRouter>::iterator it = pRouters[type]->begin();
  it->deploy(resources, network, stats, _task);
  delete[] rem;
  delete[] rem2;
}

void sosmBroker::timestep(const cell* clCell)
{
  // The simulation phase consists of 6 execution steps
  // 1. Initialize the running quantities: Proc and memory utilization, rho, network. Set initial values to 0.
  // 2. Calculate for all tasks per vRM their running quantites (except network) and add them to calulcate the total
  // 3. Calculate the number of instructions per processor can be execute for CPUs and accelerator
  // 4. Calculate the power consumption on each timestep
  // 5. Execute the tasks, and subtract the instructions already calculated
  // 6. Remove the completed tasks

  int i = 0, j, rID, len, k;
  double insR, insRa;
  double procUtil;
  double rhoAcc, L_totalPowerConsumption;
  int active, L_numberOfVMs;
  int totalAcc;
  int omp_thr = atoi(getenv("OMP_NUM_THREADS"));
  int chunk = 50;
  double* L_net = nullptr;
  list<vRM>::iterator itt, **itf;
  list<task>::iterator ittt;
  double ocP;
  int tid, jj;
  if (alloc) {
    // Step 1
    // Create a two-dimensional array of pointers to different vRMs (itf) in order to asssign them later to different
    // OpenMP threads
    itf = new list<vRM>::iterator*[numberOfTypes];
    for (i = 0; i < numberOfTypes; i++) {
      itf[i] = new list<vRM>::iterator[omp_thr];
    }
    for (i = 0; i < numberOfTypes; i++) {
      itf[i][0] = vRMs[i]->begin();
    }
    for (i = 0; i < numberOfTypes; i++) {
      len = (int)vRMs[i]->size();
      for (j = 1; j < omp_thr; j++) {
        itt = itf[i][j - 1];
        for (k = 0; k < (j * len) / omp_thr - ((j - 1) * len) / omp_thr; k++) {
          itt++;
        }
        itf[i][j] = itt;
      }
    }

    for (i = 0; i < numberOfTypes; i++) {
// Assign different resources (j) of the same hardware type (i) to different threads
#pragma omp parallel for default(shared) private(j) num_threads(omp_thr) schedule(static, chunk)
      for (j = 0; j < numberOfResourcesPerType[i]; j++)
        // If a given resource has one or more tasks assigned
        if (clCell->getResources()[i][j].getRunningVMs() > 0) {
          // Reset the overcommitment value in terms of processors and memory
          clCell->getResources()[i][j].initializeRunningQuantities();
        }
    }
    // Reset the overcommitment value in terms of network
    clCell->getNetwork()[0].initializeRunningQuantities();

    // Step 2
    omp_thr = 1;
    L_net = new double[omp_thr];
    for (i = 0; i < omp_thr; i++) {
      L_net[i] = 0.0;
    }
    // Scan all vRMs and for each vRM task retrieve mem and proc utilization and rho (the parallel percentage)
    // Calculate total utilization (by increasing the running quantities) in order to caluclate power consumption and
    // intsructions per seconds to reduce from each task
    for (i = 0; i < numberOfTypes; i++) {
      len = (int)vRMs[i]->size();
#pragma omp parallel default(shared) private(itt, ittt, k, rID, j, tid, L_numberOfVMs) num_threads(omp_thr)
      {
        tid = omp_get_thread_num();
        itt = itf[i][tid];
        // For every vRM..
        for (k = (tid * len) / omp_thr; k < ((tid + 1) * len) / omp_thr; k++) {
          // Scan the task queue of the vRM..
          for (ittt = (itt->gqueue())->begin(); ittt != (itt->gqueue())->end(); ittt++) {
            // Compute the utilization
            ittt->compcUtilPMNr();
            double* gcU = ittt->gcUtilPMNr();
            int* gr = ittt->gresourceIDs();
            L_net[tid] += gcU[2];
            L_numberOfVMs = ittt->getNumberOfVMs();
            for (j = 0; j < L_numberOfVMs; j++) {
              rID = gr[j];
              // For every [hardware type][server per hardware type] increase the mem and proc utilization and rho
              // 1 vRm has 5 tasks; for every task see which servers it uses, and on every server calculate the running
              // quantities it utilizes
              clCell->getResources()[i][rID].incrementRunningQuantities(gcU[0], gcU[1], gcU[3]);
            }
          }
          itt++;
        }
#pragma omp barrier
      }
    }

    for (i = 1; i < omp_thr; i++) {
      L_net[0] += L_net[i];
    }
    clCell->getNetwork()[0].incrementRunningQuantities(L_net[0]);
    delete[] L_net;

    // Step 3
    omp_thr = atoi(getenv("OMP_NUM_THREADS"));

    // Run through all resources and calculate the number of instructions to calculate per processing unit (taking into
    // account the overcommitment of the processor resources)
    for (i = 0; i < numberOfTypes; i++) {
#pragma omp parallel for default(shared) private(j) num_threads(omp_thr) schedule(static, chunk)
      for (j = 0; j < numberOfResourcesPerType[i]; j++) {
        if (clCell->getResources()[i][j].getRunningVMs() > 0) {
          clCell->getResources()[i][j].compcurrentCompCapPerProc();
          clCell->getResources()[i][j].compcurrentCompCapPerAcc();
        }
      }
    }
    // Step 4
    // Calculate power consumption per pc and increase consumption on the stats engine
    for (i = 0; i < numberOfTypes; i++) {
      L_totalPowerConsumption = 0.0;
#pragma omp parallel for default(shared) private(j, procUtil, rhoAcc, active, totalAcc) num_threads(omp_thr) \
                                                   schedule(static) reduction(+ : L_totalPowerConsumption)
      for (j = 0; j < numberOfResourcesPerType[i]; j++) {
        procUtil = clCell->getResources()[i][j].getActualUtilizedProcessors() /
                   clCell->getResources()[i][j].getTotalProcessors();
        rhoAcc = clCell->getResources()[i][j].getActualRhoAccelerators();
        active = clCell->getResources()[i][j].getActive();
        totalAcc = clCell->getResources()[i][j].getTotalAccelerators();
        L_totalPowerConsumption += clCell->getPowerConsumption()[i].consumption(procUtil, rhoAcc, active, totalAcc);
      }
      clCell->getStats()[i].totalPowerConsumption += L_totalPowerConsumption;
    }
    // Step 5
    // Run through all tasks to calculate the minimum of the instructions that can be executed based on the utilization
    // and the overcommitment of each node
    for (i = 0; i < numberOfTypes; i++) {
      ocP = clCell->getResources()[i][0].getOvercommitmentProcessors();
      len = (int)vRMs[i]->size();
// For all vRms..
#pragma omp parallel default(shared) private(itt, ittt, jj, j, tid, rID, insR, insRa, \
                                             L_numberOfVMs) num_threads(omp_thr)
      {
        tid = omp_get_thread_num();
        itt = itf[i][tid];
        // For all tasks of the vRM
        for (jj = (tid * len) / omp_thr; jj < ((tid + 1) * len) / omp_thr; jj++) {
          for (ittt = itt->gqueue()->begin(); ittt != itt->gqueue()->end(); ittt++) {
            rID = (ittt->gresourceIDs())[0];
            insR = clCell->getResources()[i][rID].getCurrentCompCapPerProc();
            insRa = clCell->getResources()[i][rID].getCurrentCompCapPerAcc();
            L_numberOfVMs = ittt->getNumberOfVMs();
            // For each VM of the server (resource)
            for (j = 1; j < L_numberOfVMs; j++) {
              rID = (ittt->gresourceIDs())[j];
              // Get the minimum CPU task instructions that can be calculated on the given server (resource) between the
              // requested instructions and the servers' available instructions based on the servers overcommitment
              // value
              insR = min(insR, clCell->getResources()[i][rID].getCurrentCompCapPerProc());
              // Get the minimum accelerator task instructions that can be calculated on the given server (resource)
              // between the requested instructions and the servers' available instructions based on the servers
              // overcommitment value
              insRa = min(insRa, clCell->getResources()[i][rID].getCurrentCompCapPerAcc());
            }
            // Reduce the instructions left to compute of the given task
            ittt->reduceIns(L_numberOfVMs * insR * min(ittt->gcUtilPMNr()[0] * ocP, 1.0) +
                            L_numberOfVMs * insRa * ((ittt->gcUtilPMNr())[3]));
          }
          itt++;
        }
#pragma omp barrier
      }
    }
    // Step 6
    // Estimate when a task has finished, by calculating when the number of instructions left goes to zero. (on the
    // previous step the resources were reduced gradually). Run through all vRMs and free the task' resources that are
    // done and the network
    int* numberOfTasks;
    double** avau;
    numberOfTasks = new int[omp_thr];
    avau = new double*[omp_thr];

    for (i = 0; i < omp_thr; i++) {
      avau[i] = new double[2];
    }

    for (i = 0; i < numberOfTypes; i++) {
      for (j = 0; j < omp_thr; j++) {
        avau[j][0] = 0.0;
        avau[j][1] = 0.0;
        numberOfTasks[j] = 0;
      }

      len = (int)vRMs[i]->size();
#pragma omp parallel default(shared) private(itt, ittt, jj, j, tid, rID, L_numberOfVMs) num_threads(omp_thr)
      {
        tid = omp_get_thread_num();
        itt = itf[i][tid];
        for (jj = (tid * len) / omp_thr; jj < ((tid + 1) * len) / omp_thr; jj++) {
          ittt = (itt->gqueue())->begin();
          while (ittt != (itt->gqueue())->end()) {
            // If the requested instruction of a task were dropped to zero
            if ((ittt->grequestedInstructions()) <= 0.0) {
              L_numberOfVMs = ittt->getNumberOfVMs();
              // For all VMs
              for (j = 0; j < L_numberOfVMs; j++) {
                rID = (ittt->gresourceIDs())[j];
                // Mark the resources used by the task as not used any more and therefore "free" them to be used by
                // another task in a future timestep
                clCell->getResources()[i][rID].unload(ittt);
              }
              numberOfTasks[tid]++;
              avau[tid][0] += (ittt->greqPMNS())[2];
              avau[tid][1] += (ittt->gcUtilPMNr())[2];
              // Erase the task (ittt) from the task queue of the specific vRM (itt)
              ittt = (itt->gqueue())->erase(ittt);
            } else
              ++ittt;
          }
          itt++;
        }
#pragma omp barrier
      }
      for (j = 1; j < omp_thr; j++) {
        avau[0][0] += avau[j][0];
        avau[0][1] += avau[j][1];
        numberOfTasks[0] += numberOfTasks[j];
      }
      clCell->getNetwork()[0].unload(avau[0][0], avau[0][1], numberOfTasks[0]);
    }
    for (i = 0; i < omp_thr; i++) {
      delete[] avau[i];
    }
    delete[] avau;
    delete[] numberOfTasks;

    for (i = 0; i < numberOfTypes; i++) {
      delete[] itf[i];
    }
    delete[] itf;
  }
}

int sosmBroker::getNumberOfvRMs() const { return numberOfvRMs; }
int sosmBroker::getNumberOfpSwitches() const { return numberOfpSwitches; }
int sosmBroker::getNumberOfpRouters() const { return numberOfpRouters; }
double sosmBroker::gpollIntervalCellM() const { return pollIntervalCellM; }
double sosmBroker::gpollIntervalpRouter() const { return pollIntervalpRouter; }
double sosmBroker::gpollIntervalpSwitch() const { return pollIntervalpSwitch; }
double sosmBroker::gpollIntervalvRM() const { return pollIntervalvRM; }
double** sosmBroker::gsPMSA() const { return sPMSA; }
double* sosmBroker::gSIs() const { return SIs; }
double* sosmBroker::gCs() const { return Cs; }
double* sosmBroker::gPs() const { return Ps; }
double* sosmBroker::gPis() const { return Pis; }
double* sosmBroker::gWs() const { return Ws; }
int sosmBroker::getNumberOfFunctions() const { return numberOfFunctions; }
list<vRM>** sosmBroker::getvRMs() const { return vRMs; }
list<pSwitch>** sosmBroker::gpSwitches() const { return pSwitches; }
list<pRouter>** sosmBroker::gpRouters() const { return pRouters; }

#include <gs.h>
#include <task.h>
#include <tce.h>
#include <cmath>
#include <cstdlib> // for RAND_MAX

void taskCreationEngine(std::list<task>& jobs, const struct appinputs& app)
{
  srand(0);
  double r = ((double)rand()) / RAND_MAX;

  int numberOfTasks = round(app.minmaxJobsPerSec[0] + r * (app.minmaxJobsPerSec[1] - app.minmaxJobsPerSec[0]));

  for (int i = 0; i < numberOfTasks; i++) {
    r = ((double)rand()) / RAND_MAX;
    int applicationType = floor(app.numOfApps * r);   // Choose randomly an application type of the input AppData.json file
    r = ((double)rand()) / RAND_MAX;
    
    double curr_requestedInstructions =
      app.minmaxInsPerApp[applicationType][0] + r * (app.minmaxInsPerApp[applicationType][1] - app.minmaxInsPerApp[applicationType][0]);

    int curr_numberOfVMs = round(app.minmaxVMPerApp[applicationType][0] +
                                 r * (app.minmaxVMPerApp[applicationType][1] - app.minmaxVMPerApp[applicationType][0]));

    // Requirements for resources
    double curr_reqP = round(app.minmaxProcPerVM[applicationType][0] +
                             r * (app.minmaxProcPerVM[applicationType][1] - app.minmaxProcPerVM[applicationType][0]));
    double curr_reqM =
      app.minmaxMemPerVM[applicationType][0] + r * (app.minmaxMemPerVM[applicationType][1] - app.minmaxMemPerVM[applicationType][0]);
    double curr_reqN =
      app.minmaxNetPerApp[applicationType][0] + r * (app.minmaxNetPerApp[applicationType][1] - app.minmaxNetPerApp[applicationType][0]);
    double curr_reqS =
      app.minmaxStoPerVM[applicationType][0] + r * (app.minmaxStoPerVM[applicationType][1] - app.minmaxStoPerVM[applicationType][0]);

    jobs.push_back(task(applicationType, app.numberOfAvailableImplementationsPerApp[applicationType],
                        &app.availableImplementationsPerApp[applicationType][0], curr_requestedInstructions, curr_numberOfVMs,
                        curr_reqP, curr_reqM, curr_reqN, curr_reqS, app.typeOfActP[applicationType],
                        app.typeOfActM[applicationType], app.typeOfActN[applicationType], &app.minmaxActP[applicationType][0],
                        &app.minmaxActM[applicationType][0], &app.minmaxActN[applicationType][0], &app.accelerator[applicationType][0],
                        &app.rhoAcc[applicationType][0]));
  }
}

void taskImplSelect(std::list<task>& jobs)
{
  int numberOfTasks = jobs.size(), implementationType, i, j, k;
  double r;

  int L_type, L_numberOfAvailableImplementations;
  int* L_availableImplementations;
  double L_requestedInstructions;
  int L_numberOfVMs;
  double* L_reqPMNS;
  int* L_typeactPMN;
  double** L_minmaxactPMN;
  int* L_avAcc;
  double* L_rhoAcc;

  std::list<task>::iterator it = jobs.begin();

  for (i = 0; i < numberOfTasks; i++) {
    r = ((double)rand()) / RAND_MAX;
    implementationType = floor(r * ((double)(*it).getNumberOfAvailableImplementations()));
    L_numberOfAvailableImplementations = 1;
    L_type = (*it).getType();
    L_availableImplementations = new int[1];
    L_availableImplementations[0] = (*it).getAvailableImplementations()[implementationType];
    L_requestedInstructions = (*it).grequestedInstructions();
    L_numberOfVMs = (*it).getNumberOfVMs();
    L_reqPMNS = new double[4];
    for (j = 0; j < 4; j++) {
      L_reqPMNS[j] = (*it).greqPMNS()[j];
    }
    L_typeactPMN = new int[3];
    for (j = 0; j < 3; j++) {
      L_typeactPMN[j] = (*it).getTypeactPMN()[j];
    }
    L_minmaxactPMN = new double*[3];
    for (j = 0; j < 3; j++) {
      L_minmaxactPMN[j] = new double[2];
    }
    for (j = 0; j < 3; j++) {
      for (k = 0; k < 2; k++) {
        L_minmaxactPMN[j][k] = (*it).gminmaxactPMN()[j][k];
      }
    }
    L_avAcc = new int[1];
    L_rhoAcc = new double[1];
    L_avAcc[0] = (*it).gavAcc()[implementationType];
    L_rhoAcc[0] = (*it).grhoAcc()[implementationType];

    (*it) =
      (task(L_type, L_numberOfAvailableImplementations, L_availableImplementations, L_requestedInstructions,
            L_numberOfVMs, L_reqPMNS[0], L_reqPMNS[1], L_reqPMNS[2], L_reqPMNS[3], L_typeactPMN[0], L_typeactPMN[1],
            L_typeactPMN[2], &L_minmaxactPMN[0][0], &L_minmaxactPMN[1][0], &L_minmaxactPMN[2][0], L_avAcc, L_rhoAcc));
    delete[] L_availableImplementations;
    delete[] L_reqPMNS;
    delete[] L_typeactPMN;
    for (j = 0; j < 3; j++) {
      delete[] L_minmaxactPMN[j];
    }
    delete[] L_avAcc;
    delete[] L_rhoAcc;
    it++;
  }
}

void taskCellSelect(std::list<task>& jobs, const gs* gates, int** commCells)
{
  if (*commCells != nullptr) {
    delete[](*commCells);
    (*commCells) = nullptr;
  }
  if (jobs.size() != 0 && gates->galloc()) {
    (*commCells) = new int[jobs.size()];
    for (int i = 0; i < (int)jobs.size(); i++) {
      (*commCells)[i] = -1;
    }
    std::list<task>::iterator it = jobs.begin();
    for (int i = 0; i < (int)jobs.size(); i++) {
      // Store an index to commCells to the most appropriate cell for each task,
      // as found by the gs::findCell method
      (*commCells)[i] = gates->findCell(
        (*it).getAvailableImplementations(), (*it).getNumberOfAvailableImplementations(), (*it).getNumberOfVMs(),
        (*it).greqPMNS()[0], (*it).greqPMNS()[1], (*it).greqPMNS()[2], (*it).greqPMNS()[3], (*it).gavAcc());
      it++;
    }
  }
}

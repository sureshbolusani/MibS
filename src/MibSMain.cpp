/*===========================================================================*/
/* This file is part of a Mixed Integer Bilevel Solver                       */
/* developed using the BiCePS Linear Integer Solver (BLIS).                  */
/*                                                                           */
/* Authors: Scott DeNegre, Lehigh University                                 */
/*          Ted Ralphs, Lehigh University                                    */
/*          Sahar Tahernajad, Lehigh University                              */
/*                                                                           */
/* Copyright (C) 2007-2017 Lehigh University, Scott DeNegre, and Ted Ralphs. */
/* All Rights Reserved.                                                      */
/*                                                                           */
/* This software is licensed under the Eclipse Public License. Please see    */
/* accompanying file for terms.                                              */
/*===========================================================================*/

#define COIN_HAS_CPLEX 1

#include <iostream>
#include <vector>

#include "CoinError.hpp"

#include "OsiSolverInterface.hpp"
#include "OsiClpSolverInterface.hpp"
#include "OsiCbcSolverInterface.hpp"

#ifdef COIN_HAS_SYMPHONY
#include "symphony.h"
#include "OsiSymSolverInterface.hpp"
#endif

#ifdef COIN_HAS_CPLEX
#include "cplex.h"
#include "OsiCpxSolverInterface.hpp"
#endif

#include "MibSModel.hpp"

#if  COIN_HAS_MPI
#include "AlpsKnowledgeBrokerMPI.h"
#else
#include "AlpsKnowledgeBrokerSerial.h"
#endif


bool solveMasterProblem(double *masterBestSolution, std::string masterProblemSolver, int masterMaxThreads,
        int masterColNum, int masterRowNum, double *masterObjCoef, double upperObjSense,
        double *masterColLb, double *masterColUb, char *masterColType,
        CoinPackedMatrix *masterMat, double *masterRowLb, double *masterRowUb) {
    //Master problem solver
    OsiSolverInterface *solver;
    if (masterProblemSolver == "Cbc") {
        solver = new OsiCbcSolverInterface();
    } else if (masterProblemSolver == "SYMPHONY") {
#ifdef COIN_HAS_SYMPHONY
        solver = new OsiSymSolverInterface();
#else
        std::cout <<
            "Error: SYMPHONY chosen as solver but it has not been enabled."
            << std::endl;
        return 0;
#endif

    } else if (masterProblemSolver == "CPLEX") {
#ifdef COIN_HAS_CPLEX
        solver = new OsiCpxSolverInterface();
#else
        std::cout <<
            "Error: CPLEX chosen as solver but it has not been enabled."
            << std::endl;
        return 0;
#endif
    } else {
        std::cout << 
            "Error: Unknown solver chosen for solving master problem."
            << std::endl;
        return 0;
    }
    int i;
    // Loading master problem data and related information
    solver->loadProblem(*masterMat, masterColLb, masterColUb,
            masterObjCoef, masterRowLb, masterRowUb);
    for (i = 0; i < masterColNum; i++) {
        if (masterColType[i] == 'C') {
            //do nothing because variables are continuous by default
        } else if (masterColType[i] == 'I' or masterColType[i] == 'B') {
            solver->setInteger(i);
        } else {
            throw CoinError("Unknown column type in the master problem",
                    "solveMasterProblem",
                    "MibSMain");
        }
    }
    solver->setObjSense(upperObjSense);
    solver->setHintParam(OsiDoReducePrint, true, OsiHintDo);

    if (masterProblemSolver == "Cbc") {
        dynamic_cast<OsiCbcSolverInterface *> 
            (solver)->getModelPtr()->messageHandler()->setLogLevel(0);
    } else if (masterProblemSolver == "SYMPHONY") {
#if COIN_HAS_SYMPHONY
        sym_environment *env = dynamic_cast<OsiSymSolverInterface *> 
            (solver)->getSymphonyEnvironment();

        sym_set_int_param(env, "verbosity", -2);
        sym_set_int_param(env, "max_active_nodes", masterMaxThreads);
#else
        throw CoinError("SYMPHONY chosen as solver but it has not been enabled",
                "solveMasterProblem",
                "MibSMain");
#endif
    } else if (masterProblemSolver == "CPLEX") {
#ifdef COIN_HAS_CPLEX
        solver->setHintParam(OsiDoReducePrint);
        solver->messageHandler()->setLogLevel(0);
        CPXENVptr cpxEnv = 
            dynamic_cast<OsiCpxSolverInterface*>(solver)->getEnvironmentPtr();
        assert(cpxEnv);
        CPXsetintparam(cpxEnv, CPX_PARAM_SCRIND, CPX_OFF);
        CPXsetintparam(cpxEnv, CPX_PARAM_THREADS, masterMaxThreads);
#else
        throw CoinError("CPLEX chosen as solver but it has not been enabled",
                "solveMasterProblem",
                "MibSMain");
#endif
    }

    if (1) {
        solver->writeLp("master");
    }

    /** Solving the master problem **/
    solver->branchAndBound();


    /** Getting solution to master problem **/
    bool isInfeasible = solver->isProvenPrimalInfeasible();
    if (isInfeasible) {
        masterBestSolution = NULL;
    } else {
        memcpy(masterBestSolution, solver->getColSolution(), sizeof(double)*masterColNum);
    }

    delete solver;

    return isInfeasible;
}


//#############################################################################
//#############################################################################

int main(int argc, char* argv[])
{

    //FIXME: use paramters such as subproblemRowNum/ColNum, subproblemLowerRowNum/
    //  LowerColNum, masterRowNum/ColNum, etc to make the implementation generic!

    try{
       
      /*** Original problem data setup ***/
      // Create a MibS model to read arguments and problem data
      MibSModel origMibsModel;
      origMibsModel.readParameters(argc, argv);
      std::string dataFile = origMibsModel.AlpsPar()->entry(AlpsParams::instance);
      origMibsModel.readInstance(dataFile.c_str());
      OsiCpxSolverInterface origLpSolver;

      // Gather various matrices and vectors of MibS model
      CoinPackedMatrix rowCoefMatrixByCol = *origMibsModel.getOrigConstCoefMatrix();

      //FIXME: following assertion fails with these formulae!
      //    Need to modify lower dimension somewhere due to one extra variable?
//      int upperColNum = origMibsModel.getUpperDim();
//      int lowerColNum = origMibsModel.getLowerDim();
//      assert((upperColNum + lowerColNum) == origMibsModel.getNumOrigVars());

      int upperColNum = origMibsModel.getUpperDim();
      int lowerColNum = origMibsModel.getNumOrigVars() - upperColNum;
      assert((upperColNum + lowerColNum) == rowCoefMatrixByCol.getMajorDim());

      int upperRowNum = origMibsModel.getOrigUpperRowNum();
      int lowerRowNum = origMibsModel.getLowerRowNum();
      assert((upperRowNum + lowerRowNum) == origMibsModel.getNumOrigCons());
      assert((upperRowNum + lowerRowNum) == rowCoefMatrixByCol.getMinorDim());

      double *origColLb = origMibsModel.getOrigColLb();
      double *origColUb = origMibsModel.getOrigColUb();
      char *colType = origMibsModel.getColType();
      int *upperColInd = origMibsModel.getUpperColInd();
      //FIXME: following lowerColInd is incorrect for interdiction problems
      int *lowerColInd = origMibsModel.getLowerColInd();
      //FIXME: Remove following four lines after debugging above FIXME
      /*
      int *incorrectLowerColInd = origMibsModel.getLowerColInd();
      int *lowerColInd = new int[lowerColNum];
      memcpy(lowerColInd, incorrectLowerColInd, sizeof(int)*(lowerColNum - 1));
      lowerColInd[lowerColNum-1] = upperColNum + lowerColNum - 1;
      */

      //FIXME: since intColIndices_ is being set when broker is initialized,
      //    following line doesn't work. Hence, manual code below.
      //int *intColInd = origMibsModel.getIntColIndices();
      int i, numLowerIntCols = 0;
      int *lowerIntColInd = new int[lowerColNum];
      memset(lowerIntColInd, 0, sizeof(int)*lowerColNum);
      for (i = upperColNum; i < (upperColNum + lowerColNum); i++) {
          if (colType[i] == 'I' || colType[i] == 'B') {
              lowerIntColInd[numLowerIntCols++] = i - upperColNum;
          }
      }

      int *upperRowInd = origMibsModel.getOrigUpperRowInd();
      int *lowerRowInd = origMibsModel.getLowerRowInd();
      double *rowLb = origMibsModel.getOrigRowLb();
      double *rowUb = origMibsModel.getOrigRowUb();

      int structRowNum = origMibsModel.getStructRowNum();
      int *structRowInd = origMibsModel.getStructRowInd();

      double upperObjSense = origMibsModel.BlisPar()->entry(BlisParams::objSense);
      double *upperObjCoef = origMibsModel.getObjCoef();
      double lowerObjSense = origMibsModel.getLowerObjSense();
      //FIXME: following lowerObjCoef is incorrect for interdiction problems
      double *lowerObjCoef = origMibsModel.getLowerObjCoeffs();
      //FIXME: Remove following four lines after debugging above FIXME
      /*
      double *incorrectLowerObjCoef = origMibsModel.getLowerObjCoeffs();
      double *lowerObjCoef = new double[lowerColNum];
      memcpy(lowerObjCoef, incorrectLowerObjCoef, sizeof(double)*(lowerColNum - 1));
      //TODO: Is following coeff correct for extra col?
      lowerObjCoef[lowerColNum-1] = 1;
      */



      /*** do-while loop for decomposition algorithm ***/
      /** Initial setup for MIBLP subproblem **/
      //FIXME: Ideally, infinity should be retrieved from MibSModel itself. 
      double infinity = origLpSolver.getInfinity();

      MibSModel *subproblem;

      AlpsSubTree *ast;

#ifdef  COIN_HAS_MPI
              AlpsKnowledgeBrokerMPI *broker;
#else
              AlpsKnowledgeBrokerSerial *broker;
#endif

      //Various arrays and matrices for the problem setup
      int subproblemColNum = lowerColNum, subproblemRowNum = lowerRowNum;
      CoinPackedMatrix subproblemMat(rowCoefMatrixByCol);
      subproblemMat.deleteRows(upperRowNum, upperRowInd);
      subproblemMat.deleteCols(upperColNum, upperColInd);

      //Matrix of lower level row coeffs for upper level cols
      CoinPackedMatrix lowerMatOfUpperCols(rowCoefMatrixByCol);
      lowerMatOfUpperCols.deleteRows(upperRowNum, upperRowInd);
      lowerMatOfUpperCols.deleteCols(lowerColNum, lowerColInd);

      double *subproblemColLb = new double[lowerColNum];
      double *subproblemColUb = new double[lowerColNum];
      char *subproblemColType = new char[lowerColNum];
      double *subproblemUpperObjCoef = new double[lowerColNum];
      memcpy(subproblemColLb, &origColLb[upperColNum], sizeof(double)*lowerColNum);
      memcpy(subproblemColUb, &origColUb[upperColNum], sizeof(double)*lowerColNum);
      memcpy(subproblemColType, &colType[upperColNum], sizeof(char)*lowerColNum);
      memcpy(subproblemUpperObjCoef, &upperObjCoef[upperColNum], sizeof(double)*lowerColNum);

      //Note: row LB, UB, & RHS will be modified in each iteration for a given master solution
      double *subproblemRowLb = new double[lowerRowNum];
      double *subproblemRowUb = new double[lowerRowNum];
      double *subproblemRhs = new double[lowerRowNum];
      memset(subproblemRhs, 0, sizeof(double)*lowerRowNum);
      double *subproblemUpperColRowActivity = new double[lowerRowNum];
      memset(subproblemUpperColRowActivity, 0, sizeof(double)*lowerRowNum);

      char *subproblemRowSense = new char[lowerRowNum];
      int *subproblemLowerColInd = new int[lowerColNum];
      int *subproblemLowerRowInd = new int[lowerRowNum];
      CoinIotaN(subproblemLowerColInd, lowerColNum, 0);
      CoinIotaN(subproblemLowerRowInd, lowerRowNum, 0);

      for (i = upperRowNum; i < (upperRowNum + lowerRowNum); i++) {
          if (rowLb[i] > -infinity && rowUb[i] < infinity) {
              std::cout << 
                  "Error: MibS can handle only <= or >= type constraints at present."
                  << std::endl;
              return 0;
          } else if (rowLb[i] > -infinity) {
              subproblemRowSense[i - upperRowNum] = 'G';
          } else if (rowUb[i] < infinity) {
              subproblemRowSense[i - upperRowNum] = 'L';
          }
      }

      int subproblemStructRowNum = 0;
      std::vector<int> subproblemStructRowIndVec;
      for (i = 0; i < structRowNum; i++) {
          if (structRowInd[i] >= upperRowNum) {
              subproblemStructRowIndVec.push_back((structRowInd[i] - upperRowNum));
              subproblemStructRowNum++;
          }
      }
      int *subproblemStructRowInd = &subproblemStructRowIndVec[0];

      double etol = origMibsModel.getTolerance();
      MibSWarmStart *subproblemWS;
      int leafNodeNum = 0;
      int *leafDepth = NULL;
      BlisLpStatus *leafFeasibilityStatus = NULL;
//      MibSBranchObjectInt **leafBranchPath = NULL;
      CoinPackedMatrix *leafDualByRow = NULL;
      CoinPackedMatrix *leafPosDjByRow = NULL;
      CoinPackedMatrix *leafNegDjByRow = NULL;
      //FIXME: remove following two lines after fixing leafBranchPath
      double **leafLb = NULL;
      double **leafUb = NULL;
      bool *leafDualInfoUsageStatus = NULL;


      /** Initial setup for MILP master problem **/
      //Assumption: no 2nd level cols in 1st level rows ==> 1st level rows are here!
      //Various arrays and matrices for the problem setup
      int masterColNum = upperColNum, masterRowNum = upperRowNum;
      int j, branchId, branchDirection;
      bool masterInfeasible = false;
      double branchValue, BIGM = 1e+7, bilevelVFApproxValue;
      //FIXME: set # of threads appropriately later (as a part of argc/argv?)!
      int masterMaxThreads = 1;
      double *masterBestSolution;
      double *masterBestSolutionUpperCols = new double[upperColNum];

      //Copying original matrix at first and then deleting second level rows & cols
      //FIXME: instead of deleting, write a generic adding logic
      CoinPackedMatrix masterMat(rowCoefMatrixByCol);
      masterMat.deleteRows(lowerRowNum, lowerRowInd);
      masterMat.deleteCols(lowerColNum, lowerColInd);

      std::vector<double> masterColLbVec;
      std::vector<double> masterColUbVec;
      std::vector<double> masterObjCoefVec;
      std::vector<double> masterRowLbVec;
      std::vector<double> masterRowUbVec;
      std::vector<char> masterColTypeVec; 

      for (i = 0; i < upperColNum; i++) {
          masterColLbVec.push_back(origColLb[i]);
          masterColUbVec.push_back(origColUb[i]);
          masterObjCoefVec.push_back(upperObjCoef[i]);
          masterColTypeVec.push_back(colType[i]);
      }
      for (i = 0; i < masterRowNum; i++) {
          masterRowLbVec.push_back(rowLb[i]);
          masterRowUbVec.push_back(rowUb[i]);
      }

      double *masterColLb = &masterColLbVec[0];
      double *masterColUb = &masterColUbVec[0];
      double *masterObjCoef = &masterObjCoefVec[0];
      char *masterColType = &masterColTypeVec[0];
      double *masterRowLb = &masterRowLbVec[0];
      double *masterRowUb = &masterRowUbVec[0];

      //Parameter for identifying which MILP solver needs to be used
      //FIXME: Have this as a user input parameter!
      std::string masterProblemSolver = "CPLEX";


      /** Solving initial master problem without any second level information **/
      masterBestSolution = new double[masterColNum];
      masterInfeasible = solveMasterProblem(masterBestSolution, masterProblemSolver, masterMaxThreads,
              masterColNum, masterRowNum, masterObjCoef, upperObjSense,
              masterColLb, masterColUb, masterColType,
              &masterMat, masterRowLb, masterRowUb);
      memcpy(masterBestSolutionUpperCols, masterBestSolution, sizeof(double)*upperColNum);
      delete [] masterBestSolution;

      if (!masterInfeasible) {
          //FIXME: unbounded case? any other invalid case?
          bilevelVFApproxValue = -infinity;
      }

      // Adding a column to master problem representing bilevel VF approx. value
      masterColNum += 1;
      masterColLbVec.resize(masterColNum, -infinity);
      masterColUbVec.resize(masterColNum, +infinity);
      masterObjCoefVec.resize(masterColNum, 1);
      masterColTypeVec.resize(masterColNum, 'C');
      masterMat.appendCol(0, NULL, NULL);

      masterColLb = &masterColLbVec[0];
      masterColUb = &masterColUbVec[0];
      masterObjCoef = &masterObjCoefVec[0];
      masterColType = &masterColTypeVec[0];

      //Temporary LB and UB to be modified in every iteration of decomposition algorithm
      double *lowerColLb = new double[lowerColNum];
      double *lowerColUb = new double[lowerColNum];
      //Matrix coefficient column to iterate through constraint matrix columns 
      //FIXME: Does declaration initialize by default to zeroes?
      //    If yes, delete manual initialization using memset below.
      double *rowCoefMatCol = new double[lowerRowNum];
      memset(rowCoefMatCol, 0, sizeof(double)*(lowerRowNum));

      //Misc declarations
      bool termFlag = false;
      int newArgc = 1, iterCounter = 0;
      double bilevelVFExactValue;
      char** newArgv = new char* [1];
      newArgv[0] = (char *) "mibs";

      while (!termFlag) {
          /** Setting up MibS subproblem using solution to master problem **/
          if (!masterInfeasible) {
              OsiCpxSolverInterface lpSolver;
//              lpSolver.getModelPtr()->setDualBound(1.0e10);
              lpSolver.messageHandler()->setLogLevel(0);

              subproblem = new MibSModel();
              subproblem->setSolver(&lpSolver);

              subproblem->readParameters(argc, argv);
              subproblem->AlpsPar()->setEntry(AlpsParams::instance, "NONE");
              subproblem->MibSPar()->setEntry(MibSParams::auxiliaryInfoFile, "NONE");
              subproblem->AlpsPar()->setEntry(AlpsParams::msgLevel, -1);
              subproblem->MibSPar()->setEntry(MibSParams::bilevelProblemType, 0);

              //Setting subproblem Row LB, UB, RHS for the given masterBestSolution
              memcpy(subproblemRowLb, &rowLb[upperRowNum], sizeof(double)*lowerRowNum);
              memcpy(subproblemRowUb, &rowUb[upperRowNum], sizeof(double)*lowerRowNum);
              lowerMatOfUpperCols.times(masterBestSolutionUpperCols, subproblemUpperColRowActivity);
              for (i = 0; i < lowerRowNum; i++) {
                  if (subproblemRowLb[i] > -infinity) {
                      subproblemRowLb[i] -= subproblemUpperColRowActivity[i];
                      subproblemRhs[i] = subproblemRowLb[i];
                  } else if (subproblemRowUb[i] < infinity) {
                      subproblemRowUb[i] -= subproblemUpperColRowActivity[i];
                      subproblemRhs[i] = subproblemRowUb[i];
                  }
              }

              subproblem->loadAuxiliaryData(lowerColNum,
                      lowerRowNum,
                      subproblemLowerColInd,
                      subproblemLowerRowInd,
                      lowerObjSense, 
                      lowerObjCoef,
                      0, 0, NULL, NULL,
                      subproblemStructRowNum,
                      subproblemStructRowInd,
                      0, NULL);
              subproblem->loadProblemData(subproblemMat,
                      subproblemColLb, subproblemColUb,
                      subproblemUpperObjCoef,
                      subproblemRowLb, subproblemRowUb,
                      subproblemColType, upperObjSense, infinity,
                      subproblemRowSense);


              /** Solving MibS subproblem **/
#ifdef  COIN_HAS_MPI
              broker = new AlpsKnowledgeBrokerMPI(newArgc, newArgv, *subproblem);
#else
              broker = new AlpsKnowledgeBrokerSerial(newArgc, newArgv, *subproblem);
#endif

              if (1) {
                  lpSolver.writeLp("subproblem");
              }
              broker->search(subproblem);
              //FIXME: if subproblem is not solved to optimality, bestQuality is not
              //    the exact value of bilevel VF
              bilevelVFExactValue = broker->getBestQuality();


              /** Getting dual information to bilevel subproblem **/
              ast = broker->getWorkingSubTree();
              subproblem->generateMibsWarmStart(ast);
              subproblemWS = subproblem->getMibsWarmStart();

              delete broker;
              delete subproblem;
          }


          /** Checking termination criterion **/
          if (masterInfeasible || 
                  (fabs(bilevelVFExactValue - bilevelVFApproxValue) <= etol)) {
              termFlag = true;
          }

          /** Generating Benders' cuts and adding them to the master problem **/
          if (!termFlag) {
              /* Finding various dual information products */
              //START finding products
              leafNodeNum = subproblemWS->getLeafNodeNum(); 
              leafDepth = subproblemWS->getLeafDepths();
              leafFeasibilityStatus = subproblemWS->getLeafFeasibilityStati();
              leafDualInfoUsageStatus = subproblemWS->getLeafDualInfoUsageStati();
//              leafBranchPath = subproblemWS->getLeafBranchPaths();
              leafDualByRow = subproblemWS->getLeafDualsByRow();
              leafPosDjByRow = subproblemWS->getLeafPosDjsByRow();
              leafNegDjByRow = subproblemWS->getLeafNegDjsByRow();
              leafLb = subproblemWS->getLeafLBs();
              leafUb = subproblemWS->getLeafUBs();
              int posDjNonzeroNum = leafPosDjByRow->getNumElements();
              int negDjNonzeroNum = leafNegDjByRow->getNumElements();

              //Product of RHS (subproblemRhs) and dual matrix (leafDualByRow)
              double *rhsDualProduct = new double[leafNodeNum];
              memset(rhsDualProduct, 0, sizeof(double)*leafNodeNum);
              leafDualByRow->times(subproblemRhs, rhsDualProduct);

              //Product of (column LBs and positive reduced costs (leafPosDjByRow))
              //    and (column UBs and negative reduced costs (leafNegDjByRow))
              double *lbPosDjProduct = new double[leafNodeNum];
              double *ubNegDjProduct = new double[leafNodeNum];
              //FIXME: Is following initialization necessary?
              memset(lbPosDjProduct, 0, sizeof(double)*leafNodeNum);
              memset(ubNegDjProduct, 0, sizeof(double)*leafNodeNum);


              CoinShallowPackedVector posDj;
              int posDjNum = 0;
              const int *posDjIndices = NULL;
              const double *posDjElements = NULL;
              CoinShallowPackedVector negDj;
              int negDjNum = 0;
              const int *negDjIndices = NULL;
              const double *negDjElements = NULL;

              for (i = 0; i < leafNodeNum; i++) {
                  //i-th leaf node's positive and negative reduced costs
                  if (posDjNonzeroNum) {
                      //Data of i-th leaf node's positive reduced costs
                      posDj = leafPosDjByRow->getVector(i);
                      posDjNum = posDj.getNumElements();
                      posDjIndices = posDj.getIndices();
                      posDjElements = posDj.getElements();
                  }
                  if (negDjNonzeroNum) {
                      //Data of i-th leaf node's negative reduced costs
                      negDj = leafNegDjByRow->getVector(i);
                      negDjNum = negDj.getNumElements();
                      negDjIndices = negDj.getIndices();
                      negDjElements = negDj.getElements();
                  }

                  //Copy of original bounds of lower level variables 
//                  memcpy(lowerColLb, subproblemColLb, sizeof(double)*(lowerColNum));
//                  memcpy(lowerColUb, subproblemColUb, sizeof(double)*(lowerColNum));
                  memcpy(lowerColLb, leafLb[i], sizeof(double)*(lowerColNum));
                  memcpy(lowerColUb, leafUb[i], sizeof(double)*(lowerColNum));

                  //Change original bounds to bounds at the leaf node 'i'
                  /*
                  for (j = 0; j < leafDepth[i]; j++) {
                      //TODO: need to check if variable branching?
                      branchId = lowerIntColInd[leafBranchPath[i][j].getObjectIndex()];
                      branchDirection = leafBranchPath[i][j].getDirection();
                      branchValue = leafBranchPath[i][j].getValue();

                      if (branchDirection == 1) {
                          lowerColLb[branchId] = branchValue;
                      } else {
                          lowerColUb[branchId] = branchValue;
                      }
                  }
                  */

                  //Product of bounds and reduced costs
                  for (j = 0; j < posDjNum; j++) {
                      lbPosDjProduct[i] += posDjElements[j]*lowerColLb[posDjIndices[j]];
                  }
                  for (j = 0; j < negDjNum; j++) {
                      ubNegDjProduct[i] += negDjElements[j]*lowerColUb[negDjIndices[j]];
                  }
              }

              //RHS of dual constraints (Benders' cuts)
              double *rhsDual = new double[leafNodeNum];
              for (i = 0; i < leafNodeNum; i++) {
                  rhsDual[i] = rhsDualProduct[i] + lbPosDjProduct[i] + ubNegDjProduct[i];
              }

              //Product of con matrix (subprobleMat) and dual matrix (leafDualByRow)
              double **lhsDualProduct = new double*[upperColNum];

              for (i = 0; i < upperColNum; i++) {
                  lhsDualProduct[i] = new double[leafNodeNum];
                  //FIXME: Can following memset be avoided?
                  //    Does above declaration initialize by default to zeroes?
                  memset(lhsDualProduct[i], 0, sizeof(double)*leafNodeNum);

                  CoinShallowPackedVector matCol = subproblemMat.getVector(i);

                  int matColNum = matCol.getNumElements();
                  const int *matColIndices = matCol.getIndices();
                  const double *matColElements = matCol.getElements();

                  for (j = 0; j < matColNum; j++) {
                      rowCoefMatCol[matColIndices[j]] = matColElements[j];
                  }

                  leafDualByRow->times(rowCoefMatCol, lhsDualProduct[i]);
              }
              //END finding products

              /* Updating master problem's matrices and vectors/arrays */
              int feasibleLeafNodeNum = 0;
              //Note: actual size of following would be feasibleLeafNodeNum
              int *feasibleLeafNodeInd = new int[leafNodeNum];
              memset(feasibleLeafNodeInd, 0, sizeof(int)*leafNodeNum);
              for (i = 0; i < leafNodeNum; i++) {
//                  if (leafFeasibilityStatus[i] != BlisLpStatusPrimalInfeasible) {
                  if (leafDualInfoUsageStatus[i]) {
                      //TODO: All all other stati good to consider here?
                      feasibleLeafNodeInd[feasibleLeafNodeNum] = i;
                      feasibleLeafNodeNum++;
                  }
              }
              // Note: new # of cols = old # of cols + feasibleLeafNodeNum + feasibleLeafNodeNum * upperColNum
              // Note: new # of rows = old # of rows + 4 * feasibleLeafNodeNum * upperColNum + 2

              //Objective coefficient vector
              masterObjCoefVec.resize((masterColNum + (1 + upperColNum)*feasibleLeafNodeNum), 0);
              masterObjCoef = &masterObjCoefVec[0];

              //Column bounds
              masterColLbVec.resize((masterColNum + feasibleLeafNodeNum), 0);
              masterColLbVec.resize((masterColNum + (1 + upperColNum)*feasibleLeafNodeNum), -infinity);
              masterColUbVec.resize((masterColNum + feasibleLeafNodeNum), 1);
              masterColUbVec.resize((masterColNum + (1 + upperColNum)*feasibleLeafNodeNum), +infinity);
              masterColLb = &masterColLbVec[0];
              masterColUb = &masterColUbVec[0];

              //Column types
              masterColTypeVec.resize((masterColNum + feasibleLeafNodeNum), 'B');
              //Note: following 'C' variables correspond to product to integer and binary variables
              //FIXME: Explicit mention of 'I' for following variables required?
              masterColTypeVec.resize((masterColNum + (1 + upperColNum)*feasibleLeafNodeNum), 'C');
              masterColType = &masterColTypeVec[0];

              //Row bounds
              masterRowLbVec.resize((masterRowNum + feasibleLeafNodeNum*upperColNum), 0);
              masterRowLbVec.resize((masterRowNum + 2*feasibleLeafNodeNum*upperColNum), -infinity);
              for (i = 0; i < upperColNum; i++) {
                  masterRowLbVec.resize((masterRowNum + 2*feasibleLeafNodeNum*upperColNum + (i+1)*feasibleLeafNodeNum),
                          -origColUb[i]);
              }
              masterRowLbVec.resize((masterRowNum + 4*feasibleLeafNodeNum*upperColNum), -infinity);
              masterRowLbVec.resize((masterRowNum + 4*feasibleLeafNodeNum*upperColNum + 1), 1);
              masterRowLbVec.resize((masterRowNum + 4*feasibleLeafNodeNum*upperColNum + 2), 0);

              masterRowUbVec.resize((masterRowNum + feasibleLeafNodeNum*upperColNum), +infinity);
              masterRowUbVec.resize((masterRowNum + 2*feasibleLeafNodeNum*upperColNum),0);
              masterRowUbVec.resize((masterRowNum + 3*feasibleLeafNodeNum*upperColNum), +infinity);
              for (i = 0; i < upperColNum; i++) {
                  masterRowUbVec.resize((masterRowNum + 3*feasibleLeafNodeNum*upperColNum + (i+1)*feasibleLeafNodeNum),
                          -origColLb[i]);
              }
              masterRowUbVec.resize((masterRowNum + 4*feasibleLeafNodeNum*upperColNum + 1), 1);
              masterRowUbVec.resize((masterRowNum + 4*feasibleLeafNodeNum*upperColNum + 2), +infinity);
              masterRowLb = &masterRowLbVec[0];
              masterRowUb = &masterRowUbVec[0];

              //Row coefficient matrix
              //First, appending new columns to existing matrix
              int newColNum = (1 + upperColNum)*feasibleLeafNodeNum;
              CoinBigIndex *newColStarts = new CoinBigIndex[newColNum + 1];
              //FIXME: Is following memset required, or does declaration take care of it?
              memset(newColStarts, 0, sizeof(CoinBigIndex)*(newColNum + 1));
              int errorNum = masterMat.appendCols(newColNum, newColStarts, 
                      NULL, NULL, masterRowNum);
              assert(errorNum == 0);


              //Now, appending new rows to the matrix
              for (i = 0; i < upperColNum; i++) {
                  for (j = 0; j < feasibleLeafNodeNum; j++) {
                      CoinPackedVector row;
                      row.insert(masterColNum + j, -origColLb[i]);
                      row.insert(masterColNum + feasibleLeafNodeNum + i*feasibleLeafNodeNum + j, 1);
                      masterMat.appendRow(row);
                  }
              }
              for (i = 0; i < upperColNum; i++) {
                  for (j = 0; j < feasibleLeafNodeNum; j++) {
                      CoinPackedVector row;
                      row.insert(masterColNum + j, -origColUb[i]);
                      row.insert(masterColNum + feasibleLeafNodeNum + i*feasibleLeafNodeNum + j, 1);
                      masterMat.appendRow(row);
                  }
              }
              for (i = 0; i < upperColNum; i++) {
                  for (j = 0; j < feasibleLeafNodeNum; j++) {
                      CoinPackedVector row;
                      row.insert(i, -1);
                      row.insert(masterColNum + j, -origColUb[i]);
                      row.insert(masterColNum + feasibleLeafNodeNum + i*feasibleLeafNodeNum + j, 1);
                      masterMat.appendRow(row);
                  }
              }
              for (i = 0; i < upperColNum; i++) {
                  for (j = 0; j < feasibleLeafNodeNum; j++) {
                      CoinPackedVector row;
                      row.insert(i, -1);
                      row.insert(masterColNum + j, -origColLb[i]);
                      row.insert(masterColNum + feasibleLeafNodeNum + i*feasibleLeafNodeNum + j, 1);
                      masterMat.appendRow(row);
                  }
              }
              CoinPackedVector oneMoreRow;
              for (i = 0; i < feasibleLeafNodeNum; i++) {
                  oneMoreRow.insert(masterColNum + i, 1);
              }
              masterMat.appendRow(oneMoreRow);
              CoinPackedVector oneLastRow;
              oneLastRow.insert(upperColNum, 1);
              for (i = 0; i < feasibleLeafNodeNum; i++) {
                  oneLastRow.insert(masterColNum + i, -rhsDual[feasibleLeafNodeInd[i]]);
              }
              for (i = 0; i < upperColNum; i++) {
                  for (j = 0; j < feasibleLeafNodeNum; j++) {
                      oneLastRow.insert(masterColNum + feasibleLeafNodeNum + 
                              i*feasibleLeafNodeNum + j, lhsDualProduct[i][feasibleLeafNodeInd[j]]);
                  }
              }
              masterMat.appendRow(oneLastRow);

              //Number of rows and columns
              masterColNum += (1 + upperColNum)*feasibleLeafNodeNum;
              masterRowNum += 4*feasibleLeafNodeNum*upperColNum + 2;


              /** Setting and solving the master problem **/
              masterBestSolution = new double[masterColNum];
              masterInfeasible = solveMasterProblem(masterBestSolution, masterProblemSolver, masterMaxThreads,
                      masterColNum, masterRowNum, masterObjCoef, upperObjSense,
                      masterColLb, masterColUb, masterColType,
                      &masterMat, masterRowLb, masterRowUb);
              memcpy(masterBestSolutionUpperCols, masterBestSolution, sizeof(double)*upperColNum);

              /** Getting solution to master problem **/
              if (!masterInfeasible) {
                  //FIXME: unbounded case? any other invalid case?
                  bilevelVFApproxValue = masterBestSolution[upperColNum];
              }

              delete [] masterBestSolution;
              delete [] newColStarts;
              delete [] lhsDualProduct;
              delete [] rhsDual;
              delete [] ubNegDjProduct;
              delete [] lbPosDjProduct;
              delete [] rhsDualProduct;
              delete [] feasibleLeafNodeInd;
          }
          std::cout << std::endl;
          std::cout << "Iter-" << iterCounter+1 << std::endl;
          std::cout << "masterInfeas = " << masterInfeasible << std::endl;
          std::cout << "VF Exact = " << bilevelVFExactValue << ", VF Approx = " << bilevelVFApproxValue << std::endl;
          std::cout << std::endl;
          iterCounter++;
      }

      delete [] newArgv;
      delete [] rowCoefMatCol;
      delete [] lowerColUb;
      delete [] lowerColLb;
      delete [] masterBestSolutionUpperCols;
      delete [] subproblemLowerRowInd;
      delete [] subproblemLowerColInd;
      delete [] subproblemRowSense;
      delete [] subproblemUpperColRowActivity;
      delete [] subproblemRhs;
      delete [] subproblemRowUb;
      delete [] subproblemRowLb;
      delete [] subproblemUpperObjCoef;
      delete [] subproblemColType;
      delete [] subproblemColUb;
      delete [] subproblemColLb;
//      delete [] lowerObjCoef;
      delete [] lowerIntColInd;
//      delete [] lowerColInd;
    }
    catch(CoinError& er) {
	std::cerr << "ERROR:" << er.message() << std::endl
		  << " from function " << er.methodName() << std::endl
		  << " from class " << er.className() << std::endl;
    }
    catch(...) {
	std::cerr << "Something went wrong!" << std::endl;
    }

    return 0;
}

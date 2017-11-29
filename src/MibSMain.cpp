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
#define COIN_HAS_SYMPHONY 1

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


OsiSolverInterface* getSolver(std::string problemSolver, int maxThreads,
        bool dualInfoRequired) {
    //Check for solver when dual information is required
    if (dualInfoRequired && problemSolver != "SYMPHONY") {
        throw CoinError("Choose SYMPHONY if dual information is required",
                "getSolver",
                "MibSMain");
    }

    //Problem solver
    OsiSolverInterface *solver;
    if (problemSolver == "Cbc") {
        solver = new OsiCbcSolverInterface();

        dynamic_cast<OsiCbcSolverInterface *> 
            (solver)->getModelPtr()->messageHandler()->setLogLevel(0);
    } else if (problemSolver == "SYMPHONY") {
#ifdef COIN_HAS_SYMPHONY
        solver = new OsiSymSolverInterface();

        sym_environment *env = dynamic_cast<OsiSymSolverInterface *> 
            (solver)->getSymphonyEnvironment();

        sym_set_int_param(env, "verbosity", -2);
        sym_set_int_param(env, "max_active_nodes", maxThreads);

        if (dualInfoRequired) {
            sym_set_int_param(env, "keep_warm_start", TRUE);
            sym_set_int_param(env, "warm_start_type", 1);
            sym_set_int_param(env, "sensitivity_analysis", TRUE);
            sym_set_int_param(env, "sensitivity_rhs", TRUE);
            sym_set_int_param(env, "sensitivity_bounds", TRUE);
            sym_set_int_param(env, "do_reduced_cost_fixing", FALSE);
            sym_set_int_param(env, "do_primal_heuristic", FALSE);
            sym_set_int_param(env, "generate_cgl_cuts", FALSE);
            sym_set_int_param(env, "should_use_rel_br", FALSE);
            sym_set_int_param(env, "prep_level", 0);
        }
#else
        throw CoinError("SYMPHONY chosen as solver but it has not been enabled",
                "getSolver",
                "MibSMain");
#endif

    } else if (problemSolver == "CPLEX") {
#ifdef COIN_HAS_CPLEX
        solver = new OsiCpxSolverInterface();

        solver->setHintParam(OsiDoReducePrint);
        solver->messageHandler()->setLogLevel(0);
        CPXENVptr cpxEnv = 
            dynamic_cast<OsiCpxSolverInterface*>(solver)->getEnvironmentPtr();
        assert(cpxEnv);
        CPXsetintparam(cpxEnv, CPX_PARAM_SCRIND, CPX_OFF);
        CPXsetintparam(cpxEnv, CPX_PARAM_THREADS, maxThreads);
#else
        throw CoinError("CPLEX chosen as solver but it has not been enabled",
                "getSolver",
                "MibSMain");
#endif
    } else {
        throw CoinError("Unknown solver chosen",
                "getSolver",
                "MibSMain");
    }
    solver->setHintParam(OsiDoReducePrint, true, OsiHintDo);

    return solver;
}



//#############################################################################
//#############################################################################


bool solve(OsiSolverInterface *solver,
        int colNum, double *objCoef, double objSense,
        double *colLb, double *colUb, char *colType,
        CoinPackedMatrix *mat, double *rowLb, double *rowUb,
        double *objVal, double *bestSolution) {
    int i;
    bool isLp = true;
    // Loading problem data and related information
    solver->loadProblem(*mat, colLb, colUb,
            objCoef, rowLb, rowUb);
    for (i = 0; i < colNum; i++) {
        if (colType[i] == 'C') {
            //do nothing because variables are continuous by default
        } else if (colType[i] == 'I' or colType[i] == 'B') {
            isLp = false;
            solver->setInteger(i);
        } else {
            throw CoinError("Unknown column type in the problem",
                    "solve",
                    "MibSMain");
        }
    }
    solver->setObjSense(objSense);

    if (0) {
        solver->writeLp("problemAtHand");
    }

    // Solving the problem
    if (isLp) {
        solver->initialSolve();
    } else {
        solver->branchAndBound();
    }

    // Getting solution to the problem
    bool isInfeasible = solver->isProvenPrimalInfeasible();
    if (isInfeasible) {
        bestSolution = NULL;
    } else {
        memcpy(bestSolution, solver->getColSolution(), sizeof(double)*colNum);
        *objVal = solver->getObjValue();
    }

    return isInfeasible;
}


//#############################################################################
//#############################################################################


void getDualData(OsiSolverInterface *solver, int *nodeNum, int **feasibilityStatus,
        CoinPackedMatrix *dual, CoinPackedMatrix *posDj, CoinPackedMatrix *negDj,
        int **lbCnt, int ***lbInd, double ***lbVal,
        int **ubCnt, int ***ubInd, double ***ubVal) {
#ifdef COIN_HAS_SYMPHONY
        sym_environment *env = dynamic_cast<OsiSymSolverInterface *> 
            (solver)->getSymphonyEnvironment();

        sym_get_num_leaves(env, nodeNum);

        *feasibilityStatus = new int[*nodeNum];
        sym_get_leaf_feas_stats(env, *feasibilityStatus);

        *lbCnt = new int[*nodeNum];
        *lbInd = new int*[*nodeNum];
        *lbVal = new double*[*nodeNum];
        *ubCnt = new int[*nodeNum];
        *ubInd = new int*[*nodeNum];
        *ubVal = new double*[*nodeNum];
        sym_get_branchdesc_bounds(env, *lbCnt, *lbInd, *lbVal, *ubCnt, *ubInd, *ubVal);

        sym_get_leaf_duals_by_row(env, dual);
        sym_get_leaf_pos_djs_by_row(env, posDj);
        sym_get_leaf_neg_djs_by_row(env, negDj);
#else
        throw CoinError("SYMPHONY chosen as solver but it has not been enabled",
                "getDualData",
                "MibSMain");
#endif
}


//#############################################################################
//#############################################################################

int main(int argc, char* argv[])
{

    //FIXME: use paramters such as subproblemRowNum/ColNum, subproblemLowerRowNum/
    //  LowerColNum, masterRowNum/ColNum, etc to make the implementation generic!

    //FIXME: It is assumed that lower rows follow strictly after upper rows.
    //  Is it OK? Does MibS do this reformulation internally or not?

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

      int *upperRowInd = origMibsModel.getOrigUpperRowInd();
      int *lowerRowInd = origMibsModel.getLowerRowInd();
      double *rowLb = origMibsModel.getOrigRowLb();
      double *rowUb = origMibsModel.getOrigRowUb();

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

      //FIXME: since intColIndices_ is being set when broker is initialized,
      //    following line doesn't work. Hence, manual code below.
      //int *intColInd = origMibsModel.getIntColIndices();
      //Other data is used in continuous restriction of second level MILP.
      int i, lowerContColNum = 0, lowerIntColNum = 0;
      int *lowerIntColInd = new int[lowerColNum];
      int *lowerContColInd = new int[lowerColNum];
      double *lowerContObjCoef = new double[lowerColNum];
      double *lowerContColLb = new double[lowerColNum];
      double *lowerContColUb = new double[lowerColNum];
      memset(lowerIntColInd, 0, sizeof(int)*lowerColNum);
      memset(lowerContColInd, 0, sizeof(int)*lowerColNum);
      for (i = upperColNum; i < (upperColNum + lowerColNum); i++) {
          if (colType[i] == 'I' || colType[i] == 'B') {
              lowerIntColInd[lowerIntColNum++] = i - upperColNum;
          } else {
              lowerContColInd[lowerContColNum] = i - upperColNum;
              lowerContObjCoef[lowerContColNum] = lowerObjCoef[i - upperColNum];
              lowerContColLb[lowerContColNum] = origColLb[i];
              lowerContColUb[lowerContColNum++] = origColUb[i];
          }
      }

      //FIXME: Ideally, infinity should be retrieved from MibSModel itself. 
      CoinMpsIO *mps = new CoinMpsIO;
//      double infinity = origLpSolver.getInfinity();
      double infinity = mps->getInfinity();
      double etol = origMibsModel.getTolerance();


      /** Initial setup for the bilevel subproblem as an MILP **/
      //Various arrays and matrices for the problem setup
      int subproblemColNum = lowerColNum, subproblemRowNum = lowerRowNum + 1;
      CoinPackedMatrix subproblemMat(rowCoefMatrixByCol);
      subproblemMat.deleteRows(upperRowNum, upperRowInd);
      subproblemMat.deleteCols(upperColNum, upperColInd);
      //One more row for objective-bound type constraint
      int *subproblemLowerColInd = new int[subproblemColNum];
      CoinIotaN(subproblemLowerColInd, subproblemColNum, 0);
      subproblemMat.appendRow(subproblemColNum, subproblemLowerColInd, lowerObjCoef);

      //Matrix of lower level row coeffs for upper level cols (A^2 x)
      CoinPackedMatrix lowerMatOfUpperCols(rowCoefMatrixByCol);
      lowerMatOfUpperCols.deleteRows(upperRowNum, upperRowInd);
      lowerMatOfUpperCols.deleteCols(lowerColNum, lowerColInd);

      double *subproblemColLb = new double[subproblemColNum];
      double *subproblemColUb = new double[subproblemColNum];
      char *subproblemColType = new char[subproblemColNum];
      double *subproblemObjCoef = new double[subproblemColNum];
      memcpy(subproblemColLb, &origColLb[upperColNum], sizeof(double)*subproblemColNum);
      memcpy(subproblemColUb, &origColUb[upperColNum], sizeof(double)*subproblemColNum);
      memcpy(subproblemColType, &colType[upperColNum], sizeof(char)*subproblemColNum);
      memcpy(subproblemObjCoef, &upperObjCoef[upperColNum], sizeof(double)*subproblemColNum);

      //Note: row LB, UB, RHS & upperColRowActivity are modified in each 
      //    iteration for a given master solution
      double *subproblemRowLb = new double[subproblemRowNum];
      double *subproblemRowUb = new double[subproblemRowNum];
      double *subproblemRhs = new double[subproblemRowNum];
      memset(subproblemRhs, 0, sizeof(double)*subproblemRowNum);
      double *subproblemUpperColRowActivity = new double[lowerRowNum];
      memset(subproblemUpperColRowActivity, 0, sizeof(double)*lowerRowNum);

      char *subproblemRowSense = new char[subproblemRowNum];
      //lowerIneqRowNum is for continuous restriction later!
      int lowerIneqRowNum = 0;
      //lowerRowRhs will be used while updating master problem in every iteration
      double *lowerRowRhs = new double[lowerRowNum];
      for (i = upperRowNum; i < (upperRowNum + lowerRowNum); i++) {
          if (rowLb[i] > -infinity && rowUb[i] < infinity) {
              //FIXME: following cout should not exist because this algo. does
              //    not require MibS as such..
              std::cout << 
                  "Error: MibS can handle only <= or >= type constraints at present."
                  << std::endl;
              return 0;
          } else if (rowLb[i] > -infinity) {
              subproblemRowSense[i - upperRowNum] = 'G';
              lowerRowRhs[i - upperRowNum] = rowLb[i];
              lowerIneqRowNum++;
          } else if (rowUb[i] < infinity) {
              subproblemRowSense[i - upperRowNum] = 'L';
              lowerRowRhs[i - upperRowNum] = rowUb[i];
              lowerIneqRowNum++;
          }
      }
      //Row sense for objective-bound type constraint
      subproblemRowSense[lowerRowNum] = ((lowerObjSense == 1) ? 'L' : 'G');

      //Temporary LB and UB to be modified in every iteration of decomposition algorithm
      double *tempSubproblemColLb = new double[subproblemColNum];
      double *tempSubproblemColUb = new double[subproblemColNum];

      //Parameter for identifying which MILP solver needs to be used
      //FIXME: Have this as a user input parameter!
      std::string subproblemSolver = "SYMPHONY";
      bool subproblemInfeasible = false;
      double *subproblemBestSolution = new double[subproblemColNum];
      double subproblemObjVal;
      //FIXME: set # of threads appropriately later (as a part of argc/argv?)!
      int subproblemMaxThreads = 1;

      //For gathering dual information
      int leafNodeNum = 0;
      int *leafFeasibilityStatus = NULL;
      CoinPackedMatrix *leafDualByRow = new CoinPackedMatrix();
      CoinPackedMatrix *leafPosDjByRow = new CoinPackedMatrix();
      CoinPackedMatrix *leafNegDjByRow = new CoinPackedMatrix();
      int *leafLbCnt = NULL;
      int *leafUbCnt = NULL;
      int **leafLbInd = NULL;
      int **leafUbInd = NULL;
      double **leafLbVal = NULL;
      double **leafUbVal = NULL;


      /** Initial setup for MILP master problem **/
      //Assumption: no 2nd level cols in 1st level rows ==> 1st level rows are here!
      //FIXME: generalize this code later to remove above assumption
      //Various arrays and matrices for the problem setup
      int masterColNum = upperColNum, masterRowNum = upperRowNum;
      int j;
      bool masterInfeasible = false;
      double bilevelVFApproxValue = -infinity;
      //FIXME: set # of threads appropriately later (as a part of argc/argv?)!
      int masterMaxThreads = 1;
      double *masterBestSolution;
      double masterObjVal;
      double *masterBestSolutionUpperCols = new double[upperColNum];

      //Copying original matrix at first and then deleting second level rows & cols
      //FIXME: instead of deleting, better to write a generic adding logic?
      CoinPackedMatrix masterMat(rowCoefMatrixByCol);
      masterMat.deleteRows(lowerRowNum, lowerRowInd);
      masterMat.deleteCols(lowerColNum, lowerColInd);

      std::vector<double> masterColLbVec;
      std::vector<double> masterColUbVec;
      std::vector<char> masterColTypeVec; 
      std::vector<double> masterObjCoefVec;
      std::vector<double> masterRowLbVec;
      std::vector<double> masterRowUbVec;

      for (i = 0; i < masterColNum; i++) {
          masterColLbVec.push_back(origColLb[i]);
          masterColUbVec.push_back(origColUb[i]);
          masterColTypeVec.push_back(colType[i]);
          masterObjCoefVec.push_back(upperObjCoef[i]);
      }
      for (i = 0; i < masterRowNum; i++) {
          masterRowLbVec.push_back(rowLb[i]);
          masterRowUbVec.push_back(rowUb[i]);
      }

      double *masterColLb = &masterColLbVec[0];
      double *masterColUb = &masterColUbVec[0];
      char *masterColType = &masterColTypeVec[0];
      double *masterObjCoef = &masterObjCoefVec[0];
      double *masterRowLb = &masterRowLbVec[0];
      double *masterRowUb = &masterRowUbVec[0];

      //Parameter for identifying which MILP solver needs to be used
      //FIXME: Have this as a user input parameter!
      std::string masterProblemSolver = "CPLEX";


      /** Solving initial master problem without any second level information **/
      masterBestSolution = new double[masterColNum];
      OsiSolverInterface *solver = getSolver(masterProblemSolver, masterMaxThreads, false);
      masterInfeasible = solve(solver,
              masterColNum, masterObjCoef, upperObjSense,
              masterColLb, masterColUb, masterColType,
              &masterMat, masterRowLb, masterRowUb,
              &masterObjVal, masterBestSolution);
      delete solver;
      memcpy(masterBestSolutionUpperCols, masterBestSolution, sizeof(double)*upperColNum);
      delete [] masterBestSolution;

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


      /** Initial setup for second level problem **/
      bool level2Infeasible = false;
      double *level2BestSolution = new double[lowerColNum];
      double level2ObjVal;
      double *level2IntBestSolution = new double[lowerIntColNum];
      memset(level2IntBestSolution, 0, sizeof(double)*lowerIntColNum);
      CoinPackedMatrix level2Mat(rowCoefMatrixByCol);
      level2Mat.deleteRows(upperRowNum, upperRowInd);
      level2Mat.deleteCols(upperColNum, upperColInd);
      //Parameter for identifying which MILP solver needs to be used
      //FIXME: Have this as a user input parameter!
      std::string level2ProblemSolver = "CPLEX";
      //FIXME: set # of threads appropriately later (as a part of argc/argv?)!
      int level2MaxThreads = 1;


      /** Initial setup for continuous restriction of second level problem **/
      //NOTE: contRestRowNum = lowerRowNum
      int contRestColNum = lowerContColNum + lowerIneqRowNum;
      double *contRestObjCoef = new double[contRestColNum];
      memset(contRestObjCoef, 0, sizeof(double)*contRestColNum);
      memcpy(contRestObjCoef, lowerContObjCoef, sizeof(double)*lowerContColNum);
      double *contRestColLb = new double[contRestColNum];
      memset(contRestColLb, 0, sizeof(double)*contRestColNum);
      memcpy(contRestColLb, lowerContColLb, sizeof(double)*lowerContColNum);
      double *contRestColUb = new double[contRestColNum];
      memcpy(contRestColUb, lowerContColUb, sizeof(double)*lowerContColNum);
      for (i = lowerContColNum; i < contRestColNum; i++) {
          contRestColUb[i] = infinity;
      }
      char *contRestColType = new char[contRestColNum];
      for (i = 0; i < contRestColNum; i++) {
          contRestColType[i] = 'C';
      }
      double *contRestRowLb = new double[lowerRowNum];
      double *contRestRowUb = new double[lowerRowNum];
      CoinPackedMatrix contRestMat(rowCoefMatrixByCol);
      contRestMat.deleteRows(upperRowNum, upperRowInd);
      contRestMat.deleteCols(upperColNum, upperColInd);
      contRestMat.deleteCols(lowerIntColNum, lowerIntColInd);
      int extraColNum = 0;
      for (i = 0; i < lowerRowNum; i++) {
          //FIXME: remove dependency on subproblem data?
          if (subproblemRowSense[i] == 'L') {
              CoinPackedVector col;
              col.insert(extraColNum, 1);
              contRestMat.appendCol(col);
              extraColNum++;
          } else if (subproblemRowSense[i] == 'G') {
              CoinPackedVector col;
              col.insert(extraColNum, -1);
              contRestMat.appendCol(col);
              extraColNum++;
          } else {
              //should not happen due to a check earlier!
          }
      }
      assert(extraColNum == lowerIneqRowNum);
      bool contRestInfeasible = false;
      double *contRestBestSolution = new double[contRestColNum];
      double *contRestDualSolution = new double[lowerRowNum];
      double contRestObjVal = 0;
      double **contRestBasisInverseRow = new double*[lowerRowNum];

      //integer restriction matrix of second level rows for second level cols
      CoinPackedMatrix intRestMat(rowCoefMatrixByCol);
      intRestMat.deleteRows(upperRowNum, upperRowInd);
      intRestMat.deleteCols(upperColNum, upperColInd);
      intRestMat.deleteCols(lowerContColNum, lowerContColInd);
      double *level2IntColRowActivity = new double[lowerRowNum];
      memset(level2IntColRowActivity, 0, sizeof(double)*lowerRowNum);
      //Parameter for identifying which LP solver needs to be used
      //FIXME: Have this as a user input parameter!
      std::string contRestProblemSolver = "CPLEX";
      //FIXME: set # of threads appropriately later (as a part of argc/argv?)!
      int contRestMaxThreads = 1;


      //Misc declarations
      bool termFlag = false;
      int iterCounter = 0;
      double bilevelVFExactValue, bigM = 1e+7;

      /*** while loop for decomposition algorithm ***/
      while (!termFlag) {
          /** Setting up the subproblem in MILP form using master problem's solution **/
          if (!masterInfeasible) {
              //Setting subproblem Row LB, UB, RHS for the given masterBestSolution!
              //NOTE: The data (b^2 - A^2 x) can also be used for solving second level MILP below.
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

              /* Solving second level MILP and gathering required data */
              solver = getSolver(level2ProblemSolver, level2MaxThreads, false);
              level2Infeasible = solve(solver,
                      lowerColNum, lowerObjCoef, lowerObjSense,
                      subproblemColLb, subproblemColUb, subproblemColType,
                      &level2Mat, subproblemRowLb, subproblemRowUb,
                      &level2ObjVal, level2BestSolution);
              delete solver;
              for (i = 0; i < lowerIntColNum; i++) {
                  level2IntBestSolution[i] = level2BestSolution[lowerIntColInd[i]];
              }

              if (!level2Infeasible) {
                  //If second level feasible, solve the continuous restriction
                  //    for the known level2IntBestSolution

                  //Finding product of integer best solution and integer matrix
                  intRestMat.times(level2IntBestSolution, level2IntColRowActivity);

                  //Finding RowLb and RowUb for continuous restriction
                  memcpy(contRestRowLb, subproblemRowLb, sizeof(double)*lowerRowNum);
                  memcpy(contRestRowUb, subproblemRowUb, sizeof(double)*lowerRowNum);
                  for (i = 0; i < lowerRowNum; i++) {
                      if (contRestRowLb[i] > -infinity) {
                          contRestRowLb[i] -= level2IntColRowActivity[i];
                          //Following line because all are rows are '=' type now!
                          contRestRowUb[i] = contRestRowLb[i];
                      } else if (contRestRowUb[i] < infinity) {
                          contRestRowUb[i] -= level2IntColRowActivity[i];
                          //Following line because all are rows are '=' type now!
                          contRestRowLb[i] = contRestRowUb[i];
                      }
                  }

                  //Actual solving
                  solver = getSolver(contRestProblemSolver, contRestMaxThreads, false);
                  contRestInfeasible = solve(solver,
                          contRestColNum, contRestObjCoef, lowerObjSense,
                          contRestColLb, contRestColUb, contRestColType,
                          &contRestMat, contRestRowLb, contRestRowUb,
                          &contRestObjVal, contRestBestSolution);
                  memcpy(contRestDualSolution, solver->getRowPrice(), sizeof(double)*lowerRowNum);
                  for (i = 0; i < lowerRowNum; i++) {
                      contRestBasisInverseRow[i] = new double[lowerRowNum];
                      solver->getBInvRow(i, contRestBasisInverseRow[i]);
                  }
                  delete solver;

                  if (contRestInfeasible) {
                      //Should not happen because level2Infeasible is false!
                      std::cout << 
                          "Error: Restriction is infeasible whereas it should not be!."
                          << std::endl;
                      return 0;
                  }
              } else {
                  //FIXME: What here?
                  //    Option-1: simply say that UBF=infty, solve subproblem and go ahead.
                  //    Option-2: skip solving subproblem and go ahead solving master problem
                  //        by adding a bound cut of type cx > cx^* to remove earlier x solution?
                  contRestObjVal = infinity;
                  throw CoinError("Second level problem is infeasible for the given first level solution",
                          "Main",
                          "MibSMain");
              }


              /* Updating subproblem data with new info. from continuous restriction */
              if (lowerObjSense == 1) {
                  subproblemRowLb[lowerRowNum] = -infinity;
                  subproblemRowUb[lowerRowNum] = contRestObjVal;
              } else {
                  subproblemRowLb[lowerRowNum] = contRestObjVal;
                  subproblemRowUb[lowerRowNum] = infinity;
              }
              subproblemRhs[lowerRowNum] = contRestObjVal;

              /** Solving the subproblem **/
              solver = getSolver(subproblemSolver, subproblemMaxThreads, true);
              subproblemInfeasible = solve(solver,
                      subproblemColNum, subproblemObjCoef, upperObjSense,
                      subproblemColLb, subproblemColUb, subproblemColType,
                      &subproblemMat, subproblemRowLb, subproblemRowUb,
                      &subproblemObjVal, subproblemBestSolution);
              bilevelVFExactValue = subproblemObjVal;


              /** Getting dual information to the subproblem **/
              getDualData(solver, &leafNodeNum, &leafFeasibilityStatus,
                      leafDualByRow, leafPosDjByRow, leafNegDjByRow,
                      &leafLbCnt, &leafLbInd, &leafLbVal,
                      &leafUbCnt, &leafUbInd, &leafUbVal);
              delete solver;
          }


          /** Checking termination criterion **/
          //FIXME: are the criteria correct?
          if (masterInfeasible || 
                  (fabs(bilevelVFExactValue - bilevelVFApproxValue) <= etol)) {
              termFlag = true;
          }

          /** Generating Benders' cuts and adding them to the master problem **/
          if (!termFlag) {
              /* Finding various dual information products */
              //Separating leafDualByRow into two parts
              CoinPackedMatrix leafDualToOrigRows(*leafDualByRow);
              int delColNum = 1, minorDim = leafDualToOrigRows.getMinorDim();
              int *delColInd = new int[delColNum];
              delColInd[0] = subproblemRowNum - 1;
              if ((subproblemRowNum - 1) < minorDim) {
                  leafDualToOrigRows.deleteCols(delColNum, delColInd);
              }
              CoinPackedMatrix copyLeafDualMat(*leafDualByRow);
              copyLeafDualMat.reverseOrdering();
              int majorDim = copyLeafDualMat.getMajorDim();
              double *fullDualOfExtraRow = new double[leafNodeNum];
              memset(fullDualOfExtraRow, 0, sizeof(double)*leafNodeNum);
              if ((subproblemRowNum - 1) < majorDim) {
                  CoinShallowPackedVector leafDualOfExtraRow = copyLeafDualMat.getVector(subproblemRowNum - 1);
                  int leafDualOfExtraRowNum = leafDualOfExtraRow.getNumElements();
                  const int *leafDualOfExtraRowInd = leafDualOfExtraRow.getIndices();
                  const double *leafDualOfExtraRowVal = leafDualOfExtraRow.getElements();
                  for (i = 0; i < leafDualOfExtraRowNum; i++) {
                      fullDualOfExtraRow[leafDualOfExtraRowInd[i]] = leafDualOfExtraRowVal[i];
                  }
              }

              //START finding products
              //Product of constraint matrix (A^2) and dual info.
              double **product1 = new double*[leafNodeNum];
              CoinShallowPackedVector singleDualRow;
              int singleDualNnz, majorDimDual = leafDualToOrigRows.getMajorDim();
              int allDualNnz = leafDualToOrigRows.getNumElements();
              const int *singleDualInd;
              const double *singleDualVal;

              //Products of (column LBs and positive reduced costs (leafPosDjByRow))
              //    and (column UBs and negative reduced costs (leafNegDjByRow))
              double *lbPosDjProduct = new double[leafNodeNum];
              double *ubNegDjProduct = new double[leafNodeNum];
              //FIXME: Is following initialization necessary?
              memset(lbPosDjProduct, 0, sizeof(double)*leafNodeNum);
              memset(ubNegDjProduct, 0, sizeof(double)*leafNodeNum);
              CoinShallowPackedVector posDj;
              int posDjNum, majorDimPosDj = leafPosDjByRow->getMajorDim();
              int allPosDjNnz = leafPosDjByRow->getNumElements();
              const int *posDjIndices;
              const double *posDjElements;
              CoinShallowPackedVector negDj;
              int negDjNum, majorDimNegDj = leafNegDjByRow->getMajorDim();
              int allNegDjNnz = leafNegDjByRow->getNumElements();
              const int *negDjIndices;
              const double *negDjElements;

              for (i = 0; i < leafNodeNum; i++) {
                  //if nonzero dual entries exist, find the product
                  if (allDualNnz && (i < majorDimDual)) {
                      product1[i] = new double[upperColNum];
                      memset(product1[i], 0, sizeof(double)*upperColNum);
                      singleDualRow = leafDualToOrigRows.getVector(i);
                      singleDualNnz = singleDualRow.getNumElements();
                      singleDualInd = singleDualRow.getIndices();
                      singleDualVal = singleDualRow.getElements();
                      //Building a full vector from nonzeroes
                      double *singleDual = new double[subproblemRowNum - 1];
                      memset(singleDual, 0, sizeof(double)*(subproblemRowNum - 1));
                      for (j = 0; j < singleDualNnz; j++) {
                          singleDual[singleDualInd[j]] = singleDualVal[j];
                      }
                      //Product
                      lowerMatOfUpperCols.transposeTimes(singleDual, product1[i]);
                  } else {
                      product1[i] = new double[upperColNum];
                      memset(product1[i], 0, sizeof(double)*upperColNum);
                  }

                  //i-th leaf node's positive and negative reduced costs
                  if (allPosDjNnz && (i < majorDimPosDj)) {
                      //Data of i-th leaf node's positive reduced costs
                      posDj = leafPosDjByRow->getVector(i);
                      posDjNum = posDj.getNumElements();
                      posDjIndices = posDj.getIndices();
                      posDjElements = posDj.getElements();

                      //Copy of original bounds of subproblem 
                      memcpy(tempSubproblemColLb, subproblemColLb, sizeof(double)*(subproblemColNum));

                      //Change original bounds to bounds at the leaf node 'i'
                      for (j = 0; j < leafLbCnt[i]; j++) {
                          if (leafLbVal[i][j] > tempSubproblemColLb[leafLbInd[i][j]]) {
                              tempSubproblemColLb[leafLbInd[i][j]] = leafLbVal[i][j];
                          }
                      }

                      //Product of bounds and reduced costs
                      for (j = 0; j < posDjNum; j++) {
                          lbPosDjProduct[i] += posDjElements[j]*tempSubproblemColLb[posDjIndices[j]];
                      }
                  }
                  if (allNegDjNnz && (i < majorDimNegDj)) {
                      //Data of i-th leaf node's negative reduced costs
                      negDj = leafNegDjByRow->getVector(i);
                      negDjNum = negDj.getNumElements();
                      negDjIndices = negDj.getIndices();
                      negDjElements = negDj.getElements();

                      //Copy of original bounds of subproblem 
                      memcpy(tempSubproblemColUb, subproblemColUb, sizeof(double)*(subproblemColNum));

                      //Change original bounds to bounds at the leaf node 'i'
                      for (j = 0; j < leafUbCnt[i]; j++) {
                          if (leafUbVal[i][j] > tempSubproblemColUb[leafUbInd[i][j]]) {
                              tempSubproblemColUb[leafUbInd[i][j]] = leafUbVal[i][j];
                          }
                      }

                      //Product of bounds and reduced costs
                      for (j = 0; j < negDjNum; j++) {
                          ubNegDjProduct[i] += negDjElements[j]*tempSubproblemColUb[negDjIndices[j]];
                      }
                  }
              }

              //Product of constraint matrix (A^2) and dual of continuous restriction
              double *product2 = new double[upperColNum];
              memset(product2, 0, sizeof(double)*upperColNum);
              lowerMatOfUpperCols.transposeTimes(contRestDualSolution, product2);

              //Product of dual of continuous restriction and lower level's 
              //    row activity of integer restriction
              double product3 = 0;
              for (i = 0; i < lowerRowNum; i++) {
                  product3 += contRestDualSolution[i]*level2IntColRowActivity[i];
              }

              //Product of cont. rest. basis inverse and constraint matrix A^2
              double **product4 = new double*[lowerRowNum];
              for (i = 0; i < lowerRowNum; i++) {
                  product4[i] = new double[upperColNum];
                  memset(product4[i], 0, sizeof(double)*upperColNum);
                  lowerMatOfUpperCols.transposeTimes(contRestBasisInverseRow[i], product4[i]);
              }

              //Product of cont. rest. basis inverse and lower level's row
              //    activity of integer restriction
              double *product5 = new double[lowerRowNum];
              for (i = 0; i < lowerRowNum; i++) {
                  product5[i] = 0;
                  for (j = 0; j < lowerRowNum; j++) {
                      product5[i] += contRestBasisInverseRow[i][j]*
                          (level2IntColRowActivity[j] - lowerRowRhs[j]);
                  }
              }

              //Product of dual of cont. rest. and lowerRowRhs
              double product6 = 0;
              for (i = 0; i < lowerRowNum; i++) {
                  product6 += contRestDualSolution[i]*lowerRowRhs[i];
              }


              //Product of leafDualToOrigRows and lowerRowRhs
              double *product7 = new double[leafNodeNum];
              memset(product7, 0, sizeof(double)*leafNodeNum);
              leafDualToOrigRows.times(lowerRowRhs, product7);
              //END finding products

              /* Updating master problem's matrices and vectors/arrays */
              int feasibleLeafNodeNum = 0;
              //Note: actual size of following would be feasibleLeafNodeNum
              int *feasibleLeafNodeInd = new int[leafNodeNum];
              memset(feasibleLeafNodeInd, 0, sizeof(int)*leafNodeNum);
              for (i = 0; i < leafNodeNum; i++) {
                  //FIXME: this '4' corresponds to INFEASIBLE_PRUNED in SYMPHONY
                  if (leafFeasibilityStatus[i] != 4) {
                      //TODO: All all other stati good to consider here?
                      feasibleLeafNodeInd[feasibleLeafNodeNum] = i;
                      feasibleLeafNodeNum++;
                  }
              }
              // Note: new # of cols = old # of cols + (upperColNum + 2)*feasibleLeafNodeNum + 1
              // Note: new # of rows = old # of rows + (4 * upperColNum + 3)*feasibleLeafNodeNum + 2*lowerRowNum + 1

              //Objective coefficient vector
              masterObjCoefVec.resize((masterColNum + (upperColNum + 2)*feasibleLeafNodeNum + 1), 0);
              masterObjCoef = &masterObjCoefVec[0];

              //Column bounds
              masterColLbVec.resize((masterColNum + upperColNum*feasibleLeafNodeNum), -infinity);
              masterColLbVec.resize((masterColNum + (upperColNum + 2)*feasibleLeafNodeNum + 1), 0);
              masterColUbVec.resize((masterColNum + upperColNum*feasibleLeafNodeNum), infinity);
              masterColUbVec.resize((masterColNum + (upperColNum + 2)*feasibleLeafNodeNum + 1), 1);
              masterColLb = &masterColLbVec[0];
              masterColUb = &masterColUbVec[0];

              //Column types
              masterColTypeVec.resize((masterColNum + upperColNum*feasibleLeafNodeNum), 'C');
              masterColTypeVec.resize((masterColNum + (upperColNum + 2)*feasibleLeafNodeNum + 1), 'B');
              masterColType = &masterColTypeVec[0];

              //Row bounds
              masterRowLbVec.resize((masterRowNum + upperColNum*feasibleLeafNodeNum), -infinity);
              masterRowLbVec.resize((masterRowNum + 2*upperColNum*feasibleLeafNodeNum), 0);
              masterRowLbVec.resize((masterRowNum + 3*upperColNum*feasibleLeafNodeNum), -infinity);
              for (i = 0; i < upperColNum; i++) {
                  masterRowLbVec.resize((masterRowNum + 3*upperColNum*feasibleLeafNodeNum + (i+1)*feasibleLeafNodeNum),
                          -origColUb[i]);
              }
              masterRowLbVec.resize((masterRowNum + (4*upperColNum + 2)*feasibleLeafNodeNum), -infinity);
              masterRowLbVec.resize((masterRowNum + (4*upperColNum + 3)*feasibleLeafNodeNum), -1);
              for (i = 0; i < lowerRowNum; i++) {
                  masterRowLbVec.resize((masterRowNum + (4*upperColNum + 3)*feasibleLeafNodeNum + (i+1)), product5[i]);
              }
              masterRowLbVec.resize((masterRowNum + (4*upperColNum + 3)*feasibleLeafNodeNum + 2*lowerRowNum), -infinity);
              masterRowLbVec.resize((masterRowNum + (4*upperColNum + 3)*feasibleLeafNodeNum + 2*lowerRowNum + 1), 0);

              masterRowUbVec.resize((masterRowNum + upperColNum*feasibleLeafNodeNum), 0);
              masterRowUbVec.resize((masterRowNum + 2*upperColNum*feasibleLeafNodeNum), infinity);
              for (i = 0; i < upperColNum; i++) {
                  masterRowUbVec.resize((masterRowNum + 2*upperColNum*feasibleLeafNodeNum + (i+1)*feasibleLeafNodeNum),
                          -origColLb[i]);
              }
              masterRowUbVec.resize((masterRowNum + 4*upperColNum*feasibleLeafNodeNum), infinity);
              masterRowUbVec.resize((masterRowNum + (4*upperColNum + 2)*feasibleLeafNodeNum), 0);
              masterRowUbVec.resize((masterRowNum + (4*upperColNum + 3)*feasibleLeafNodeNum + lowerRowNum), infinity);
              for (i = 0; i < lowerRowNum; i++) {
                  masterRowUbVec.resize((masterRowNum + (4*upperColNum + 3)*feasibleLeafNodeNum + lowerRowNum + (i+1)), bigM + product5[i]);
              }
              masterRowUbVec.resize((masterRowNum + (4*upperColNum + 3)*feasibleLeafNodeNum + 2*lowerRowNum + 1), infinity);
              masterRowLb = &masterRowLbVec[0];
              masterRowUb = &masterRowUbVec[0];

              //Row coefficient matrix
              //First, appending new columns to existing matrix
              int newColNum = (upperColNum + 2)*feasibleLeafNodeNum + 1;
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
                      row.insert(masterColNum + i*feasibleLeafNodeNum + j, 1);
                      row.insert(masterColNum + (upperColNum + 1)*feasibleLeafNodeNum + j, -origColUb[i]);
                      masterMat.appendRow(row);
                  }
              }
              for (i = 0; i < upperColNum; i++) {
                  for (j = 0; j < feasibleLeafNodeNum; j++) {
                      CoinPackedVector row;
                      row.insert(masterColNum + i*feasibleLeafNodeNum + j, 1);
                      row.insert(masterColNum + (upperColNum + 1)*feasibleLeafNodeNum + j, -origColLb[i]);
                      masterMat.appendRow(row);
                  }
              }
              for (i = 0; i < upperColNum; i++) {
                  for (j = 0; j < feasibleLeafNodeNum; j++) {
                      CoinPackedVector row;
                      row.insert(i, -1);
                      row.insert(masterColNum + i*feasibleLeafNodeNum + j, 1);
                      row.insert(masterColNum + (upperColNum + 1)*feasibleLeafNodeNum + j, -origColLb[i]);
                      masterMat.appendRow(row);
                  }
              }
              for (i = 0; i < upperColNum; i++) {
                  for (j = 0; j < feasibleLeafNodeNum; j++) {
                      CoinPackedVector row;
                      row.insert(i, -1);
                      row.insert(masterColNum + i*feasibleLeafNodeNum + j, 1);
                      row.insert(masterColNum + (upperColNum + 1)*feasibleLeafNodeNum + j, -origColUb[i]);
                      masterMat.appendRow(row);
                  }
              }
              for (i = 0; i < feasibleLeafNodeNum; i++) {
                  CoinPackedVector row;
                  row.insert(masterColNum + upperColNum*feasibleLeafNodeNum + i, 1);
                  row.insert(masterColNum + (upperColNum + 2)*feasibleLeafNodeNum, -1);
                  masterMat.appendRow(row);
              }
              for (i = 0; i < feasibleLeafNodeNum; i++) {
                  CoinPackedVector row;
                  row.insert(masterColNum + upperColNum*feasibleLeafNodeNum + i, 1);
                  row.insert(masterColNum + (upperColNum + 1)*feasibleLeafNodeNum + i, -1);
                  masterMat.appendRow(row);
              }
              for (i = 0; i < feasibleLeafNodeNum; i++) {
                  CoinPackedVector row;
                  row.insert(masterColNum + upperColNum*feasibleLeafNodeNum + i, 1);
                  row.insert(masterColNum + (upperColNum + 1)*feasibleLeafNodeNum + i, -1);
                  row.insert(masterColNum + (upperColNum + 2)*feasibleLeafNodeNum, -1);
                  masterMat.appendRow(row);
              }
              for (i = 0; i < lowerRowNum; i++) {
                  CoinPackedVector row;
                  for (j = 0; j < upperColNum; j++) {
                      row.insert(j, -product4[i][j]);
                  }
                  row.insert(masterColNum + (upperColNum + 2)*feasibleLeafNodeNum, bigM);
                  masterMat.appendRow(row);
              }
              for (i = 0; i < lowerRowNum; i++) {
                  CoinPackedVector row;
                  for (j = 0; j < upperColNum; j++) {
                      row.insert(j, -product4[i][j]);
                  }
                  row.insert(masterColNum + (upperColNum + 2)*feasibleLeafNodeNum, bigM);
                  masterMat.appendRow(row);
              }
              CoinPackedVector oneLastRow;
              oneLastRow.insert(upperColNum, 1);
              for (i = 0; i < upperColNum; i++) {
                  for (j = 0; j < feasibleLeafNodeNum; j++) {
                      oneLastRow.insert(masterColNum + i*feasibleLeafNodeNum + j,
                              (product1[feasibleLeafNodeInd[j]][i] + 
                                  fullDualOfExtraRow[feasibleLeafNodeInd[j]]*product2[i]));
                  }
              }
              for (i = 0; i < feasibleLeafNodeNum; i++) {
                  oneLastRow.insert(masterColNum + upperColNum*feasibleLeafNodeNum + i,
                          -bigM*fullDualOfExtraRow[feasibleLeafNodeInd[i]]);
              }
              for (i = 0; i < feasibleLeafNodeNum; i++) {
                  oneLastRow.insert(masterColNum + (upperColNum + 1)*feasibleLeafNodeNum + i,
                          -(lbPosDjProduct[feasibleLeafNodeInd[i]] + ubNegDjProduct[feasibleLeafNodeInd[i]] - 
                              fullDualOfExtraRow[feasibleLeafNodeInd[i]]*product3 +
                              fullDualOfExtraRow[feasibleLeafNodeInd[i]]*product6 +
                              product7[feasibleLeafNodeInd[i]]));
              }
              masterMat.appendRow(oneLastRow);

              //Number of rows and columns
              masterColNum += (upperColNum + 2)*feasibleLeafNodeNum + 1;
              masterRowNum += (4*upperColNum + 3)*feasibleLeafNodeNum + 2*lowerRowNum + 1;


              /** Setting and solving the master problem **/
              masterBestSolution = new double[masterColNum];
              solver = getSolver(masterProblemSolver, masterMaxThreads, false);
              masterInfeasible = solve(solver,
                      masterColNum, masterObjCoef, upperObjSense,
                      masterColLb, masterColUb, masterColType,
                      &masterMat, masterRowLb, masterRowUb,
                      &masterObjVal, masterBestSolution);
              delete solver;
              memcpy(masterBestSolutionUpperCols, masterBestSolution, sizeof(double)*upperColNum);

              std::cout << std::endl;
              std::cout << "Iter-" << iterCounter << std::endl;
              std::cout << "masterInfeas = " << masterInfeasible << std::endl;
              std::cout << "VF Exact = " << bilevelVFExactValue << ", VF Approx = " << bilevelVFApproxValue << std::endl;
              std::cout << std::endl;

              /** Getting solution to master problem **/
              if (!masterInfeasible) {
                  //FIXME: unbounded case? any other invalid case?
                  bilevelVFApproxValue = masterBestSolution[upperColNum];
              }

              delete [] masterBestSolution;
              delete [] newColStarts;
              delete [] ubNegDjProduct;
              delete [] lbPosDjProduct;
              delete [] feasibleLeafNodeInd;
          } else {
              std::cout << std::endl;
              std::cout << "Iter-" << iterCounter << std::endl;
              std::cout << "masterInfeas = " << masterInfeasible << std::endl;
              std::cout << "VF Exact = " << bilevelVFExactValue << ", VF Approx = " << bilevelVFApproxValue << std::endl;
              std::cout << "Optimal Objective Value = " << masterObjVal << std::endl;
              std::cout << std::endl;
          }
          iterCounter++;
      }

      delete [] tempSubproblemColUb;
      delete [] tempSubproblemColLb;
      delete [] masterBestSolutionUpperCols;
      delete [] subproblemLowerColInd;
      delete [] subproblemRowSense;
      delete [] subproblemUpperColRowActivity;
      delete [] subproblemRhs;
      delete [] subproblemRowUb;
      delete [] subproblemRowLb;
      delete [] subproblemObjCoef;
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

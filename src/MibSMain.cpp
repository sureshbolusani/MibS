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

//FIXME: A bug related to M_P parameter for iter#3 of toy MIBLP in paper.
//FIXED: It is not a bug. It is in fact working better than expected because
//  of better bound on z due to better dualBoundOnLevel2.

#define COIN_HAS_CPLEX 1
#define COIN_HAS_SYMPHONY 1
#define COIN_HAS_SOPLEX 1

#include <iostream>
#include <vector>
#include <ctime>

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

#ifdef COIN_HAS_SOPLEX
#include "soplex.h"
#include "OsiSpxSolverInterface.hpp"
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

        solver->messageHandler()->setLogLevel(0);
        CPXENVptr cpxEnv = 
            dynamic_cast<OsiCpxSolverInterface*>(solver)->getEnvironmentPtr();
        assert(cpxEnv);
        CPXsetintparam(cpxEnv, CPX_PARAM_THREADS, maxThreads);

        double ztol = 1e-12;
        CPXsetdblparam(cpxEnv, CPX_PARAM_EPINT, ztol);
        //std::cout << "\nCPLEX integrality tolerance = " << ztol << std::endl << std::endl;

        solver->setHintParam(OsiDoReducePrint, true, OsiHintDo);
#else
        throw CoinError("CPLEX chosen as solver but it has not been enabled",
                "getSolver",
                "MibSMain");
#endif
    } else if (problemSolver == "SoPlex") {
#ifdef COIN_HAS_SOPLEX
        solver = new OsiSpxSolverInterface();

        solver->messageHandler()->setLogLevel(0);
        soplex::SoPlex *soplex =
            dynamic_cast<OsiSpxSolverInterface*>(solver)->getLpPtr();
        soplex->setIntParam(soplex::SoPlex::VERBOSITY, 0);
        //set various parameter values for exact LP solving
        soplex->setRealParam(soplex::SoPlex::FEASTOL, 0.0);
        soplex->setRealParam(soplex::SoPlex::OPTTOL, 0.0);
        soplex->setIntParam(soplex::SoPlex::SOLVEMODE, 2);
        soplex->setIntParam(soplex::SoPlex::SYNCMODE, 1);
        soplex->setIntParam(soplex::SoPlex::READMODE, 1);
        soplex->setIntParam(soplex::SoPlex::CHECKMODE, 2);
#else
        throw CoinError("SoPlex chosen as solver but it has not been enabled",
                "getSolver",
                "MibSMain");
#endif
    } else {
        throw CoinError("Unknown solver chosen",
                "getSolver",
                "MibSMain");
    }

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
    solver->setObjSense(objSense);
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

    bool indLp = true;
    bool indMps = false;
    if (indLp) {
        solver->writeLp("problemAtHand");
    }
    if (indMps) {
        solver->writeMps("problemAtHand");
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


void getDualData(OsiSolverInterface *solver, double boundOnBigM,
        int *nodeNum, int **feasibilityStatus,
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

    /* Update dual information for infeasible nodes by solving their duals */
    double symInfty = solver->getInfinity();
    int i, j, infeasibleNodeNum = 0;
    for (i = 0; i < *nodeNum; i++) {
        //FIXME: this '4' corresponds to INFEASIBLE_PRUNED in SYMPHONY
        if (feasibilityStatus[0][i] == 4)
            infeasibleNodeNum++;
    }

    if (infeasibleNodeNum) {
        //Obtain problem data from SYMPHONY
        int numCols, numRows;
        double *objCoef, *colLb, *colUb, *rowLb, *rowUb, *rowRange, *rhs;
        double objSense, objVal, bigM = 1e+7, etol = 1e-7;
        CoinPackedMatrix *matByRow;

        numCols = solver->getNumCols();
        numRows = solver->getNumRows();
        objSense = solver->getObjSense();
        //FIXME: Remove this fixme and error output upon making followup code 
        //  of dual problem building robust.
        if (objSense < 0) {
            throw CoinError("Expected a minimization problem! FIXME!",
                    "getDualData",
                    "MibSMain");
        }
        objVal = solver->getObjValue();
        if (solver->isProvenOptimal()) {
            //Note: 2 or 0.5 are random multipliers since we need bigM>objVal
            if (objVal >= etol) {
                bigM = 2*objVal;
            } else {
                bigM = 0.5*objVal;
            }
//            std::cout << bigM << "\t" << boundOnBigM << "\t" << objVal << std::endl;

            //Update bigM further based on max of subproblem's objective function
            if (bigM >= boundOnBigM + etol) {
//                std::cout << "Updating bigM Again!" << std::endl;
                bigM = boundOnBigM;
            }
        }

        objCoef = const_cast <double *> (solver->getObjCoefficients());
        colLb = const_cast <double *> (solver->getColLower());
        colUb = const_cast <double *> (solver->getColUpper());
        rowLb = const_cast <double *> (solver->getRowLower());
        rowUb = const_cast <double *> (solver->getRowUpper());
        rowRange = const_cast <double *> (solver->getRowRange());
        rhs = const_cast <double *> (solver->getRightHandSide());
        matByRow = const_cast <CoinPackedMatrix *> (solver->getMatrixByRow());

        //Data structures for dual problem of current node's LP
        //Assumption: bounds are considered as cons in current node's LP
        int dualNumCols, dualNumRows;
        double *dualObjCoef, dualObjSense, *dualColLb, *dualColUb;
        double *dualRowLb, *dualRowUb, dualObjVal, *dualBestSolution;
        double *unbddRay, *correctFullDual;
        char *dualColType;
        CoinPackedMatrix *dualMatByCol;

        dualNumRows = numCols;
        dualObjSense = -objSense;
        dualRowLb = objCoef;
        dualRowUb = objCoef;
        //Allocating max size to avoid repetitive allocation
        dualNumCols = numRows + 2*numCols; //will be modified in each iteration
        dualObjCoef = new double[dualNumCols];
        CoinZeroN(dualObjCoef, dualNumCols);
        dualColLb = new double[dualNumCols];
        CoinZeroN(dualColLb, dualNumCols);
        dualColUb = new double[dualNumCols];
        CoinZeroN(dualColUb, dualNumCols);
        dualBestSolution = new double[dualNumCols];
        CoinZeroN(dualBestSolution, dualNumCols);
        unbddRay = new double[dualNumCols];
        CoinZeroN(unbddRay, dualNumCols);
        correctFullDual = new double[dualNumCols];
        CoinZeroN(correctFullDual, dualNumCols);
        dualColType = new char[dualNumCols];
        //Setting fixed parts of entries
        memcpy(dualObjCoef, rhs, sizeof(double)*numRows);
        for (i = 0; i < numRows; i++) {
            if ((rowLb[i] > -symInfty) && (rowUb[i] < symInfty)) {
                throw CoinError("Ranged constraint detected! FIXME!",
                        "getDualData",
                        "MibSMain");
            } else if (rowLb[i] > -symInfty) {
                // >= type constraint
                dualColLb[i] = 0;
                dualColUb[i] = symInfty;
            } else if (rowUb[i] < symInfty) {
                // <= type constraint
                dualColLb[i] = -symInfty;
                dualColUb[i] = 0;
            } else {
                throw CoinError("Free constraint detected! FIXME!",
                        "getDualData",
                        "MibSMain");
            }
        }
        for (i = 0; i < dualNumCols; i++) {
            dualColType[i] = 'C';
        }
        matByRow->transpose();

        //Temporary bounds that will be changed in each iteration
        double *tempColLb, *tempColUb;
        tempColLb = new double[numCols];
        tempColUb = new double[numCols];
        //Temporary data structures (finite value indicators, IDs, reverse IDs)
        bool *finiteColLbInd, *finiteColUbInd;
        int *finiteColLbId, *finiteColUbId, *reverseId;
        finiteColLbInd = new bool[numCols];
        CoinZeroN(finiteColLbInd, numCols);
        finiteColUbInd = new bool[numCols];
        CoinZeroN(finiteColUbInd, numCols);
        finiteColLbId = new int[numCols];
        finiteColUbId = new int[numCols];
        reverseId = new int[dualNumCols];
        CoinIotaN(reverseId, numRows, 0);

#ifdef COIN_HAS_CPLEX
        bool isInfeasible = false;
        bool isOptimal = true;
        std::string problemSolver = "CPLEX";
        OsiSolverInterface *dualSolver;
        for (i = 0; i < *nodeNum; i++) {
            //FIXME: this '4' corresponds to INFEASIBLE_PRUNED in SYMPHONY
            if (feasibilityStatus[0][i] == 4) {
                dualSolver = getSolver(problemSolver, 1, false);
                dualSolver->setHintParam(OsiDoPresolveInInitial, false, OsiHintDo);
                dualSolver->setHintParam(OsiDoDualInInitial, false, OsiHintDo);

                dualMatByCol = new CoinPackedMatrix(*matByRow);

                /* Setting variable parts of entries (obj. coeff, matrix, col. bounds) */
                //At first, build current node's bounds into 'temp' structures
                memcpy(tempColLb, colLb, sizeof(double)*numCols);
                memcpy(tempColUb, colUb, sizeof(double)*numCols);
                int cnt = lbCnt[0][i];
                for (j = 0; j < cnt; j++) {
                    if (lbVal[0][i][j] > tempColLb[lbInd[0][i][j]]) {
                        tempColLb[lbInd[0][i][j]] = lbVal[0][i][j];
                    }
                }
                cnt = ubCnt[0][i];
                for (j = 0; j < cnt; j++) {
                    if (ubVal[0][i][j] < tempColUb[ubInd[0][i][j]]) {
                        tempColUb[ubInd[0][i][j]] = ubVal[0][i][j];
                    }
                }
                //Now, set variable parts of entries as required
                dualNumCols = numRows;
                for (j = 0; j < numCols; j++) {
                    if (tempColLb[j] > -symInfty) {
                        dualObjCoef[dualNumCols] = tempColLb[j];
                        dualColLb[dualNumCols] = 0;
                        dualColUb[dualNumCols] = symInfty;

                        CoinPackedVector col;
                        col.insert(j, 1);
                        dualMatByCol->appendCol(col);

                        finiteColLbInd[j] = 1;
                        finiteColLbId[j] = dualNumCols;
                        reverseId[dualNumCols] = j;
                        dualNumCols++;
                    } else {
                        finiteColLbId[j] = -1;
                    }
                }
                for (j = 0; j < numCols; j++) {
                    if (tempColUb[j] < symInfty) {
                        dualObjCoef[dualNumCols] = tempColUb[j];
                        dualColLb[dualNumCols] = -symInfty;
                        dualColUb[dualNumCols] = 0;

                        CoinPackedVector col;
                        col.insert(j, 1);
                        dualMatByCol->appendCol(col);

                        finiteColUbInd[j] = 1;
                        finiteColUbId[j] = dualNumCols;
                        reverseId[dualNumCols] = j;
                        dualNumCols++;
                    } else {
                        finiteColUbId[j] = -1;
                    }
                }

                //Solving the dual of current node's LP
                isInfeasible = solve(dualSolver,
                        dualNumCols, dualObjCoef, dualObjSense,
                        dualColLb, dualColUb, dualColType,
                        dualMatByCol, dualRowLb, dualRowUb,
                        &dualObjVal, dualBestSolution);
                isOptimal = dualSolver->isProvenOptimal();
                if (!isInfeasible && !isOptimal) {
                    memcpy(unbddRay, dualSolver->getPrimalRays(1)[0], sizeof(double)*dualNumCols);
                }
                delete dualSolver;

                if (!isInfeasible && !isOptimal) {
                    //Find appropriate multiplier for unbounded ray
                    //FIXME: Improve following etol later!
                    double lambda = 0, objCoefRayProd = 0;
                    for (j = 0; j < dualNumCols; j++) {
                        objCoefRayProd += dualObjCoef[j]*unbddRay[j];
                    }
                    lambda = (bigM - dualObjVal)/(objCoefRayProd);
                    if (lambda <= etol) {
                        assert(!((objCoefRayProd < -etol) && (bigM - dualObjVal > etol)));
                        //{lambda <= 0} ==> any positive lambda is a satisfactory value
                        lambda = 1;
                    }

                    //Evaluate full correct dual
                    for (j = 0; j < dualNumCols; j++) {
                        correctFullDual[j] = dualBestSolution[j] + lambda * unbddRay[j];
                    }
                } else if (!isInfeasible && isOptimal) {
                    //Evaluate full correct dual which is same as dual solution
                    for (j = 0; j < dualNumCols; j++) {
                        correctFullDual[j] = dualBestSolution[j];
                    }
                } else if (isInfeasible) {
                    throw CoinError("Infeasible dual detected which is supposed to be feasible! Check it!!",
                            "getDualData",
                            "MibSMain");
                }

                //Find appropriate multiplier for unbounded ray
                //FIXME: Improve following etol later!
                int nzDual = 0, nzPosDj = 0, nzNegDj = 0;
                int maxColIndDualNz = -1, maxColIndPosDjNz = -1, maxColIndNegDjNz = -1;

                //Find various paramteres
                for (j = 0; j < dualNumCols; j++) {
                    if (fabs(correctFullDual[j]) > etol) {
                        if (j < numRows) {
                            nzDual++;
                            maxColIndDualNz = j;
                        } else if (correctFullDual[j] > etol) {
                            nzPosDj++;
                            maxColIndPosDjNz = reverseId[j];
                        } else if (correctFullDual[j] < -etol) {
                            nzNegDj++;
                            maxColIndNegDjNz = reverseId[j];
                        }
                    }
                }

                //Update dual matrix
                int dualMajorDim = dual->getMajorDim();
                int dualMinorDim = dual->getMinorDim();
                if ((maxColIndDualNz >= 0) && !dualMajorDim && !dualMinorDim) {
                    //Dual matrix is currently empty
                    int newRowNum = i + 1;
                    CoinBigIndex *newRowStarts = new CoinBigIndex[newRowNum + 1];
                    CoinZeroN(newRowStarts, newRowNum + 1);
                    int errorNum = dual->appendRows(newRowNum, newRowStarts,
                            NULL, NULL, maxColIndDualNz + 1);
                    delete [] newRowStarts;
                    assert(errorNum == 0);
                    dualMajorDim = dual->getMajorDim();
                    dualMinorDim = dual->getMinorDim();
                }
                if ((maxColIndDualNz + 1) > dualMinorDim) {
                    //Minor dimension is not large enough! Append extra columns of all zeroes.
                    int newColNum = maxColIndDualNz + 1 - dualMinorDim;
                    CoinBigIndex *newColStarts = new CoinBigIndex[newColNum + 1];
                    CoinZeroN(newColStarts, newColNum + 1);
                    int errorNum = dual->appendCols(newColNum, newColStarts,
                            NULL, NULL, dualMajorDim);
                    delete [] newColStarts;
                    assert(errorNum == 0);
                    dualMinorDim = maxColIndDualNz + 1;
                }
                if ((maxColIndDualNz >= 0) && (i >= dualMajorDim)) {
                    //Major dimension is not large enough! Append extra rows of all zeroes.
                    int newRowNum = i + 1 - dualMajorDim;
                    CoinBigIndex *newRowStarts = new CoinBigIndex[newRowNum + 1];
                    CoinZeroN(newRowStarts, newRowNum + 1);
                    int errorNum = dual->appendRows(newRowNum, newRowStarts,
                            NULL, NULL, dualMinorDim);
                    delete [] newRowStarts;
                    assert(errorNum == 0);
                }
                //Now, both the minor and major dimensions are large enough
                for (j = 0; j < numRows; j++) {
                    dual->modifyCoefficient(i, j, correctFullDual[j], true);
                }

                //Expand posDj matrix if required
                int posDjMajorDim = posDj->getMajorDim();
                int posDjMinorDim = posDj->getMinorDim();
                if ((maxColIndPosDjNz >= 0) && !posDjMajorDim && !posDjMinorDim) {
                    //posDj matrix is currently empty
                    int newRowNum = i + 1;
                    CoinBigIndex *newRowStarts = new CoinBigIndex[newRowNum + 1];
                    CoinZeroN(newRowStarts, newRowNum + 1);
                    int errorNum = posDj->appendRows(newRowNum, newRowStarts,
                            NULL, NULL, maxColIndPosDjNz + 1);
                    delete [] newRowStarts;
                    assert(errorNum == 0);
                    posDjMajorDim = posDj->getMajorDim();
                    posDjMinorDim = posDj->getMinorDim();
                }
                if ((maxColIndPosDjNz + 1) > posDjMinorDim) {
                    //Minor dimension is not large enough! Append extra columns of all zeroes.
                    int newColNum = maxColIndPosDjNz + 1 - posDjMinorDim;
                    CoinBigIndex *newColStarts = new CoinBigIndex[newColNum + 1];
                    CoinZeroN(newColStarts, newColNum + 1);
                    int errorNum = posDj->appendCols(newColNum, newColStarts,
                            NULL, NULL, posDjMajorDim);
                    delete [] newColStarts;
                    assert(errorNum == 0);
                    posDjMinorDim = maxColIndPosDjNz + 1;
                }
                if ((maxColIndPosDjNz >= 0) && (i >= posDjMajorDim)) {
                    //Major dimension is not large enough! Append extra rows of all zeroes.
                    int newRowNum = i + 1 - posDjMajorDim;
                    CoinBigIndex *newRowStarts = new CoinBigIndex[newRowNum + 1];
                    CoinZeroN(newRowStarts, newRowNum + 1);
                    int errorNum = posDj->appendRows(newRowNum, newRowStarts,
                            NULL, NULL, posDjMinorDim);
                    delete [] newRowStarts;
                    assert(errorNum == 0);
                }

                //Expand negDj matrix if required
                int negDjMajorDim = negDj->getMajorDim();
                int negDjMinorDim = negDj->getMinorDim();
                if ((maxColIndNegDjNz >= 0) && !negDjMajorDim && !negDjMinorDim) {
                    //negDj matrix is currently empty
                    int newRowNum = i + 1;
                    CoinBigIndex *newRowStarts = new CoinBigIndex[newRowNum + 1];
                    CoinZeroN(newRowStarts, newRowNum + 1);
                    int errorNum = negDj->appendRows(newRowNum, newRowStarts,
                            NULL, NULL, maxColIndNegDjNz + 1);
                    delete [] newRowStarts;
                    assert(errorNum == 0);
                    negDjMajorDim = negDj->getMajorDim();
                    negDjMinorDim = negDj->getMinorDim();
                }
                if ((maxColIndNegDjNz + 1) > negDjMinorDim) {
                    //Minor dimension is not large enough! Append extra columns of all zeroes.
                    int newColNum = maxColIndNegDjNz + 1 - negDjMinorDim;
                    CoinBigIndex *newColStarts = new CoinBigIndex[newColNum + 1];
                    CoinZeroN(newColStarts, newColNum + 1);
                    int errorNum = negDj->appendCols(newColNum, newColStarts,
                            NULL, NULL, negDjMajorDim);
                    delete [] newColStarts;
                    assert(errorNum == 0);
                    negDjMinorDim = maxColIndNegDjNz + 1;
                }
                if ((maxColIndNegDjNz >= 0) && (i >= negDjMajorDim)) {
                    //Major dimension is not large enough! Append extra rows of all zeroes.
                    int newRowNum = i + 1 - negDjMajorDim;
                    CoinBigIndex *newRowStarts = new CoinBigIndex[newRowNum + 1];
                    CoinZeroN(newRowStarts, newRowNum + 1);
                    int errorNum = negDj->appendRows(newRowNum, newRowStarts,
                            NULL, NULL, negDjMinorDim);
                    delete [] newRowStarts;
                    assert(errorNum == 0);
                }

                //Now, both posDj and negDj matrices have large enough dimensions
                double djVal = 0;
                for (j = 0; j < numCols; j++) {
                    if (finiteColLbInd[j] && finiteColUbInd[j]) {
                        //Check that only one of two values is nonzero
                        double temp1 = correctFullDual[finiteColLbId[j]];
                        double temp2 = correctFullDual[finiteColUbId[j]];
                        assert(!(temp1 && temp2));
                        if (temp1) {
                            djVal = temp1;
                        } else {
                            djVal = temp2;
                        }
                    } else if (finiteColLbInd[j]) {
                        djVal = correctFullDual[finiteColLbId[j]];
                    } else if (finiteColUbInd[j]) {
                        djVal = correctFullDual[finiteColUbId[j]];
                    }
                    if (djVal > etol) {
                        posDj->modifyCoefficient(i, j, djVal, true);
                        if ((i < negDjMajorDim) && (j < negDjMinorDim)) {
                            negDj->modifyCoefficient(i, j, 0, true);
                        }
                    } else if (djVal < -etol) {
                        negDj->modifyCoefficient(i, j, djVal, true);
                        if ((i < posDjMajorDim) && (j < posDjMinorDim)) {
                            posDj->modifyCoefficient(i, j, 0, true);
                        }
                    } else {
                        if ((i < negDjMajorDim) && (j < negDjMinorDim)) {
                            negDj->modifyCoefficient(i, j, 0, true);
                        }
                        if ((i < posDjMajorDim) && (j < posDjMinorDim)) {
                            posDj->modifyCoefficient(i, j, 0, true);
                        }
                    }
                }

                delete dualMatByCol;

                //Change feasibility status to ON
                feasibilityStatus[0][i] = 1;
            } else {
                //Leaf node is feasible; change the feasibility status to ON
                feasibilityStatus[0][i] = 1;
            }
        }

#else
        throw CoinError("CPLEX is disabled while trying to use it for resolving infeasible nodes",
                "getDualData",
                "MibSMain");
#endif

        delete [] reverseId;
        delete [] finiteColUbId;
        delete [] finiteColLbId;
        delete [] finiteColUbInd;
        delete [] finiteColLbInd;
        delete [] tempColUb;
        delete [] tempColLb;
        delete [] dualColType;
        delete [] correctFullDual;
        delete [] unbddRay;
        delete [] dualBestSolution;
        delete [] dualColUb;
        delete [] dualColLb;
        delete [] dualObjCoef;
    }

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

    //FIXME: There is mix of assumptions. At times, it is assumed that lower
    //  rows follow strictly after upper rows. At times, not!
    //  Make is robust by removing the first assumption.

    //FIXME: objVals are uninitialized before passing to "solve" function.

    clock_t begin = clock();

    try{
       
      /*** Original problem data setup ***/
      // Create a MibS model to read arguments and problem data
      MibSModel origMibsModel;
      origMibsModel.readParameters(argc, argv);
      std::string dataFile = origMibsModel.AlpsPar()->entry(AlpsParams::instance);
      origMibsModel.readInstance(dataFile.c_str());

      //FIXME: Ideally, infinity should be retrieved from MibSModel itself.
      CoinMpsIO *mps = new CoinMpsIO;
      double infinity = mps->getInfinity();
      double etol = origMibsModel.getTolerance();

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
      int i;

#if 1
      //Improving bounds by solving optimization problems
      OsiSolverInterface *boundImprSolver;
      std::string boundImprProbSolver = "CPLEX";
      int boundImprMaxThreads = 1;
      bool boundImprInfeasible = false;
      double *boundImprObjSense = new double[2];
      boundImprObjSense[0] = 1.0;
      boundImprObjSense[1] = -1.0;
      double *boundImprObjCoef = new double[upperColNum + lowerColNum];
      char *boundImprColType = new char[upperColNum + lowerColNum];
      for (i = 0; i < (upperColNum + lowerColNum); i++) {
          boundImprColType[i] = 'C';
      }
      double boundImprObjVal = 0;
      double *boundImprBestSolution = new double[upperColNum + lowerColNum];
      int counter;
      //Performing bound improvement preprocessing two times
      for (counter = 0; counter < 2; counter++) {
          CoinZeroN(boundImprObjCoef, (upperColNum + lowerColNum));
          for (i = 0; i < (upperColNum + lowerColNum); i++) {
              //Setting objective function of only one variable
              boundImprObjCoef[i] = 1.0;
              if (i) {
                  boundImprObjCoef[i-1] = 0.0;
              }

              //Solving minimization problem
              boundImprSolver = getSolver(boundImprProbSolver, boundImprMaxThreads, false);
              boundImprInfeasible = solve(boundImprSolver,
                      (upperColNum + lowerColNum), boundImprObjCoef, boundImprObjSense[0],
                      origColLb, origColUb, boundImprColType,
                      &rowCoefMatrixByCol, rowLb, rowUb,
                      &boundImprObjVal, boundImprBestSolution);
              if (boundImprSolver->isProvenOptimal()) {
                  if (boundImprObjVal >= origColLb[i] + etol) {
                      if ((colType[i] == 'B') || (colType[i] == 'I')) {
                          origColLb[i] = ceil(boundImprObjVal);
                      } else {
                          //NOTE: Technically, following 'floor' is not required.
                          //  Since these bounds are useful in product9 later,
                          //    we are flooring these bounds here.
                          origColLb[i] = floor(boundImprObjVal);
                      }
                  }
              } else {
                  throw CoinError("Bound improvement problem is not solved to optimality.",
                          "main",
                          "MibSMain");
              };
              delete boundImprSolver;

              //Solving maximization problem
              boundImprSolver = getSolver(boundImprProbSolver, boundImprMaxThreads, false);
              boundImprInfeasible = solve(boundImprSolver,
                      (upperColNum + lowerColNum), boundImprObjCoef, boundImprObjSense[1],
                      origColLb, origColUb, boundImprColType,
                      &rowCoefMatrixByCol, rowLb, rowUb,
                      &boundImprObjVal, boundImprBestSolution);
              if (boundImprSolver->isProvenOptimal()) {
                  if (boundImprObjVal <= origColUb[i] - etol) {
                      if ((colType[i] == 'B') || (colType[i] == 'I')) {
                          origColUb[i] = floor(boundImprObjVal);
                      } else {
                          //NOTE: Technically, following 'ceil' is not required.
                          //  Since these bounds are useful in tolProb later for
                          //    finding epsilon values where all integer coefficients
                          //    are assumed, we are ceiling these bounds here.
                          origColUb[i] = ceil(boundImprObjVal);
                      }
                  }
              } else {
                  throw CoinError("Bound improvement problem is not solved to optimality.",
                          "main",
                          "MibSMain");
              };
              delete boundImprSolver;
          }
      }
      delete [] boundImprBestSolution;
      delete [] boundImprColType;
      delete [] boundImprObjCoef;
      delete [] boundImprObjSense;
#endif


      //FIXME: Make "z >= LBF" con. in master prob. robust to account for upperObjSense
      //    i.e., fix the "UBF + infty" or "UBF - infty" variation
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
      int lowerContColNum = 0, lowerIntColNum = 0;
      int *lowerIntColInd = new int[lowerColNum];
      int *lowerContColInd = new int[lowerColNum];
      bool *lowerContColFiniteLbId = new bool[lowerColNum];
      bool *lowerContColFiniteUbId = new bool[lowerColNum];
      double *lowerContObjCoef = new double[lowerColNum];
      double *lowerContColLb = new double[lowerColNum];
      double *lowerContColUb = new double[lowerColNum];
      CoinZeroN(lowerIntColInd, lowerColNum);
      CoinZeroN(lowerContColInd, lowerColNum);
      CoinZeroN(lowerContColFiniteLbId, lowerColNum);
      CoinZeroN(lowerContColFiniteUbId, lowerColNum);
      for (i = upperColNum; i < (upperColNum + lowerColNum); i++) {
          if (colType[i] == 'I' || colType[i] == 'B') {
              lowerIntColInd[lowerIntColNum++] = i - upperColNum;
          } else {
              lowerContColInd[lowerContColNum] = i - upperColNum;
              lowerContObjCoef[lowerContColNum] = lowerObjCoef[i - upperColNum];
              lowerContColLb[lowerContColNum] = origColLb[i];
              lowerContColUb[lowerContColNum] = origColUb[i];
              if (origColLb[i] > -infinity) {
                  lowerContColFiniteLbId[lowerContColNum] = true;
              }
              if (origColUb[i] < infinity) {
                  lowerContColFiniteUbId[lowerContColNum] = true;
              }
              lowerContColNum++;
          }
      }


      /** Determine if upper rows be part of master problem or subproblem **/
      bool upperRowsHaveLowerCols = false;
      CoinPackedMatrix upperMatOfLowerCols(rowCoefMatrixByCol);
      upperMatOfLowerCols.deleteRows(lowerRowNum, lowerRowInd);
      upperMatOfLowerCols.deleteCols(upperColNum, upperColInd);
      if (upperMatOfLowerCols.getNumElements()) {
          //upper rows have lower cols ==> upper rows will be in subproblem
          upperRowsHaveLowerCols = true;
      }


      /** Initial setup for the bilevel subproblem as an MILP **/
      //Various arrays and matrices for the problem setup
      //NOTE: subproblemRowNum may be only lowerRowNum depending on level2Infeasibility
      int subproblemColNum = lowerColNum;
      int subproblemRowNum = lowerRowNum + 1;
      //Subproblem matrix
      CoinPackedMatrix subproblemMat(rowCoefMatrixByCol);
      //Matrix of upper level cols
      CoinPackedMatrix matOfUpperCols(rowCoefMatrixByCol);
      if (!upperRowsHaveLowerCols) {
          subproblemMat.deleteRows(upperRowNum, upperRowInd);
          matOfUpperCols.deleteRows(upperRowNum, upperRowInd);
      } else {
          subproblemRowNum += upperRowNum;
      }
      subproblemMat.deleteCols(upperColNum, upperColInd);
      matOfUpperCols.deleteCols(lowerColNum, lowerColInd);
      //Matrix of upper level cols in lower rows
      CoinPackedMatrix lowerMatOfUpperCols(rowCoefMatrixByCol);
      lowerMatOfUpperCols.deleteRows(upperRowNum, upperRowInd);
      lowerMatOfUpperCols.deleteCols(lowerColNum, lowerColInd);

      int *subproblemLowerColInd = new int[subproblemColNum];
      CoinIotaN(subproblemLowerColInd, subproblemColNum, 0);
      //Indicator for adding one more row for objective-bound type constraint
      bool addRowInd = true;

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
      CoinZeroN(subproblemRhs, subproblemRowNum);
      double *subproblemUpperColRowActivity = new double[subproblemRowNum];
      CoinZeroN(subproblemUpperColRowActivity, subproblemRowNum);

      char *subproblemRowSense = new char[subproblemRowNum];
      //Original RHS (b) to be used while updating master problem in every iteration
      double *subproblemOrigRhs = new double[subproblemRowNum];
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
              lowerRowRhs[i - upperRowNum] = rowLb[i];
              lowerIneqRowNum++;
          } else if (rowUb[i] < infinity) {
              lowerRowRhs[i - upperRowNum] = rowUb[i];
              lowerIneqRowNum++;
          }
      }
      if (!upperRowsHaveLowerCols) {
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
                  subproblemOrigRhs[i - upperRowNum] = rowLb[i];
              } else if (rowUb[i] < infinity) {
                  subproblemRowSense[i - upperRowNum] = 'L';
                  subproblemOrigRhs[i - upperRowNum] = rowUb[i];
              }
          }
      } else {
          for (i = 0; i < (upperRowNum + lowerRowNum); i++) {
              if (rowLb[i] > -infinity && rowUb[i] < infinity) {
                  //FIXME: following cout should not exist because this algo. does
                  //    not require MibS as such..
                  std::cout <<
                      "Error: MibS can handle only <= or >= type constraints at present."
                      << std::endl;
                  return 0;
              } else if (rowLb[i] > -infinity) {
                  subproblemRowSense[i] = 'G';
                  subproblemOrigRhs[i] = rowLb[i];
              } else if (rowUb[i] < infinity) {
                  subproblemRowSense[i] = 'L';
                  subproblemOrigRhs[i] = rowUb[i];
              }
          }
      }
      //Row sense for objective-bound type constraint
      subproblemRowSense[subproblemRowNum - 1] = ((lowerObjSense == 1) ? 'L' : 'G');

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
      int *leafFeasibilityStatusInd = NULL;
      CoinPackedMatrix *leafDualByRow = new CoinPackedMatrix();
      CoinPackedMatrix *leafPosDjByRow = new CoinPackedMatrix();
      CoinPackedMatrix *leafNegDjByRow = new CoinPackedMatrix();
      int *leafLbCnt = NULL;
      int *leafUbCnt = NULL;
      int **leafLbInd = NULL;
      int **leafUbInd = NULL;
      double **leafLbVal = NULL;
      double **leafUbVal = NULL;
      //For detecting infeasible leaf nodes
      int feasibleLeafNodeNum, *feasibleLeafNodeInd;


      /** Initial setup for MILP master problem **/
      //Various arrays and matrices for the problem setup
      int masterColNum = upperColNum;
      int masterRowNum = (upperRowsHaveLowerCols ? 0 : upperRowNum);
      int j;
      bool masterInfeasible = false;
      double rhoApproxValue = -infinity;
      //FIXME: set # of threads appropriately later (as a part of argc/argv?)!
      int masterMaxThreads = 1;
      //Parameter for identifying which MILP solver needs to be used
      //FIXME: Have this as a user input parameter!
      std::string masterProblemSolver = "CPLEX";
      double *masterBestSolution;
      double masterObjVal, masterObjValPrevIter = 0, optObjVal = 0;
      double *masterBestSolutionUpperCols = new double[upperColNum];
      double *masterBestSolutionUpperColsPrevIter = new double[upperColNum];

      //Column and row related data
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

      //Copying original matrix at first and then deleting subproblem's rows & cols
      //FIXME: instead of deleting, better to write a generic adding logic?
      CoinPackedMatrix masterMat(rowCoefMatrixByCol);
      masterMat.deleteRows(lowerRowNum, lowerRowInd);
      masterMat.deleteCols(lowerColNum, lowerColInd);
      if (upperRowsHaveLowerCols || !upperRowNum) {
          masterMat.deleteRows(upperRowNum, upperRowInd);
          //Matrix is empty now! Add a dummy row to avoid error from CPLEX or SYMPHONY
          bool finiteLb = false;
          for (i = 0; i < masterColNum; i++) {
              if (masterColLb[i] > -infinity) {
                  finiteLb = true;
                  break;
              } else if (masterColUb[i] < infinity) {
                  break;
              }
          }
          assert(i < masterColNum);
          assert(masterRowNum == 0);
          masterRowNum = 1;
          CoinPackedVector dummyRow;
          dummyRow.insert(i, 1);
          masterMat.appendRow(dummyRow);
          if (finiteLb) {
              masterRowLbVec.push_back(masterColLb[i]);
              masterRowUbVec.push_back(infinity);
          } else {
              masterRowLbVec.push_back(-infinity);
              masterRowUbVec.push_back(masterColUb[i]);
          }
      }

      double *masterRowLb = &masterRowLbVec[0];
      double *masterRowUb = &masterRowUbVec[0];


      /** Solving initial master problem without any second level information **/
      masterBestSolution = new double[masterColNum];
      OsiSolverInterface *solver = getSolver(masterProblemSolver, masterMaxThreads, false);
      masterInfeasible = solve(solver,
              masterColNum, masterObjCoef, upperObjSense,
              masterColLb, masterColUb, masterColType,
              &masterMat, masterRowLb, masterRowUb,
              &masterObjVal, masterBestSolution);
      delete solver;
      //Rounding integer variables and reevaluating master obj. value
      masterObjVal = 0;
      for (i = 0; i < masterColNum; i++) {
          if (masterColType[i] == 'I' || masterColType[i] == 'B') {
              masterBestSolution[i] = round(masterBestSolution[i]);
          }
          masterObjVal += masterBestSolution[i]*masterObjCoef[i];
      }
      //Save upperColNum part of best solution
      memcpy(masterBestSolutionUpperCols, masterBestSolution, sizeof(double)*upperColNum);
      delete [] masterBestSolution;
      //Evaluating first level part of obj. value
      for (i = 0; i < upperColNum; i++) {
          optObjVal += masterObjCoef[i]*masterBestSolutionUpperCols[i];
      }

      //Adding a column to master problem representing rho approx. value
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
      double level2ObjVal, level2IntObjVal;
      double *level2RowLb = new double[lowerRowNum];
      double *level2RowUb = new double[lowerRowNum];
      CoinZeroN(level2RowLb, lowerRowNum);
      CoinZeroN(level2RowUb, lowerRowNum);
      CoinPackedMatrix level2Mat(rowCoefMatrixByCol);
      level2Mat.deleteRows(upperRowNum, upperRowInd);
      level2Mat.deleteCols(upperColNum, upperColInd);
      //Parameter for identifying which MILP solver needs to be used
      //FIXME: Have this as a user input parameter!
      std::string level2ProblemSolver = "CPLEX";
      //FIXME: set # of threads appropriately later (as a part of argc/argv?)!
      int level2MaxThreads = 1;


      /** Initial setup for continuous restriction of second level problem **/
      //FIXME: A restriction can be used only if first level is minimization.
      //    Make the code robust to account for first level maximization too
      //        by solving a relaxation in place of restriction.
      double *intBestSolution = new double[lowerIntColNum];
      CoinZeroN(intBestSolution, lowerIntColNum);
      //NOTE: contRestRowNum = lowerRowNum
      int contRestColNum = lowerContColNum + lowerIneqRowNum;
      double *contRestObjCoef = new double[contRestColNum];
      CoinZeroN(contRestObjCoef, contRestColNum);
      memcpy(contRestObjCoef, lowerContObjCoef, sizeof(double)*lowerContColNum);
      double *contRestColLb = new double[contRestColNum];
      CoinZeroN(contRestColLb, contRestColNum);
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
      if (!upperRowsHaveLowerCols) {
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
      } else {
          for (i = upperRowNum; i < (upperRowNum + lowerRowNum); i++) {
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
      }
      assert(extraColNum == lowerIneqRowNum);
      bool contRestInfeasible = false;
      double *contRestBestSolution = new double[contRestColNum];
      double *contRestBestSolutionNonBasics = new double[contRestColNum];
      double *contRestDualSolution = new double[lowerRowNum];
      double *contRestDjSolution = new double[contRestColNum];
      double contRestObjVal = 0;
      double **contRestBasisInverseRow = new double*[lowerRowNum];
      int *contRestBasisIndices = new int[lowerRowNum];
      int numBinColsForDomainRest = 0;
      soplex::SSVectorRational contRestBasisInverseRows(lowerRowNum*contRestColNum);
      unsigned long int *contRestBasisInverseRowLcm = new unsigned long int[lowerRowNum];

      //integer restriction matrix of second level rows for second level cols
      CoinPackedMatrix intRestMat(rowCoefMatrixByCol);
      intRestMat.deleteRows(upperRowNum, upperRowInd);
      intRestMat.deleteCols(upperColNum, upperColInd);
      intRestMat.deleteCols(lowerContColNum, lowerContColInd);
      double *level2IntColRowActivity = new double[lowerRowNum];
      CoinZeroN(level2IntColRowActivity, lowerRowNum);
      //Parameter for identifying which LP solver needs to be used
      //FIXME: Have this as a user input parameter!
      std::string contRestProblemSolver = "SoPlex";
      //FIXME: set # of threads appropriately later (as a part of argc/argv?)!
      int contRestMaxThreads = 1;


      /*
      //Data required to identify linking cols
      int *linkingColId = new int[upperColNum + lowerColNum];
      memcpy(linkingColId, origMibsModel.getFixedInd(),
              sizeof(int)*(upperColNum + lowerColNum));
      int linkingColNum = origMibsModel.getSizeFixedInd();
      */

      /** Initial setup for finding tolerance to impose strict inequality while
        restricting domain of 'x' in master problem **/
      //NOTE: 2*lowerRowNum is max size for tol
      double *tol = new double[2*lowerRowNum];
      int tolProbColNum = upperColNum;
      int tolProbRowNum = 1;
      double *tolProbObjCoef = new double[tolProbColNum];
      double tolProbObjSense = 1.0;
      double *tolProbColLb = new double[tolProbColNum];
      double *tolProbColUb = new double[tolProbColNum];
      memcpy(tolProbColLb, origColLb, sizeof(double)*tolProbColNum);
      memcpy(tolProbColUb, origColUb, sizeof(double)*tolProbColNum);
      char *tolProbColType = new char[tolProbColNum];
      memcpy(tolProbColType, colType, sizeof(char)*tolProbColNum);
      double *tolProbRowLb = new double[tolProbRowNum];
      double *tolProbRowUb = new double[tolProbRowNum];
      tolProbRowUb[0] = infinity;
      double tolProbObjVal;
      double *tolProbBestSolution = new double[tolProbColNum];
      //Triplets to initialize tolProbMat
      //Declaring maximum memory
      int maxNumEntries = tolProbRowNum * tolProbColNum;
      int *rowIndex = new int[maxNumEntries];
      CoinFillN(rowIndex, maxNumEntries, 0);
      int *colIndex = new int[maxNumEntries];
      CoinIotaN(colIndex, maxNumEntries, 0);
      double * entryVal = new double[maxNumEntries];
      CoinZeroN(entryVal, maxNumEntries);
      CoinPackedMatrix *tolProbMat = new CoinPackedMatrix(true, rowIndex, colIndex, entryVal, maxNumEntries);
      bool tolProbInfeasible = false;
      //Parameter for identifying which MILP solver needs to be used
      //FIXME: Have this as a user input parameter!
      std::string tolProbSolver = "CPLEX";
      //FIXME: set # of threads appropriately later (as a part of argc/argv?)!
      int tolProbMaxThreads = 1;

      //Misc declarations
      int tempCounter = 0;
      bool termFlag = false, timeUp = false;
      int iterCounter = 0;
      double rho = infinity, boundOnLbf, dualBoundOnLevel2 = 0;
      double bigMForUbf = 1e+7, singleBigMForLbf;
      double dualBoundOnSubproblem = 0, primalBoundOnSubproblem = 0;
      double *maxValForDomainRest = new double[lowerRowNum];
      double *minValForDomainRest = new double[lowerRowNum];
      for (i = 0; i < lowerColNum; i++) {
          double lb = subproblemColLb[i];
          assert(lb > -infinity);
          double ub = subproblemColUb[i];
          assert(ub < infinity);
          double coef1 = subproblemObjCoef[i];
          if (coef1 < -etol) {
              dualBoundOnSubproblem += coef1*lb;
              primalBoundOnSubproblem += coef1*ub;
          } else if (coef1 > etol) {
              dualBoundOnSubproblem += coef1*ub;
              primalBoundOnSubproblem += coef1*lb;
          }
          /*
          //Approach-1: Find dualBoundOnLevel2 w.r.t. simple box constraints on 'y'
          double coef2 = lowerObjCoef[i];
          if (coef2 < -etol) {
              dualBoundOnLevel2 += coef2*lb;
          } else if (coef2 > etol) {
              dualBoundOnLevel2 += coef2*ub;
          }
          */
      }
      //Approach-2: Find dualBoundOnLevel2 by solving a bilevel bounding problem
      {
          //Misc. parameters for bounding problem
          std::string boundProbLpSolverName = "CPLEX";
          int boundProbMaxThreads = 1;
          CoinPackedMatrix boundProbMat(rowCoefMatrixByCol);
          boundProbMat.deleteRows(upperRowNum, upperRowInd);
          double *boundProbUpperObjCoef = new double[upperColNum + lowerColNum];
          CoinZeroN(boundProbUpperObjCoef, upperColNum);
          for (i = upperColNum; i < (upperColNum + lowerColNum); i++) {
              boundProbUpperObjCoef[i] = -1.0 * lowerObjSense * lowerObjCoef[i - upperColNum];
          }
          double *boundProbRowLb = new double[lowerRowNum];
          double *boundProbRowUb = new double[lowerRowNum];
          memcpy(boundProbRowLb, &rowLb[upperRowNum], sizeof(double)*lowerRowNum);
          memcpy(boundProbRowUb, &rowUb[upperRowNum], sizeof(double)*lowerRowNum);
          char *boundProbRowSense = new char[lowerRowNum];
          if (upperRowsHaveLowerCols) {
              memcpy(boundProbRowSense, &subproblemRowSense[upperRowNum], sizeof(char)*lowerRowNum);
          } else {
              memcpy(boundProbRowSense, subproblemRowSense, sizeof(char)*lowerRowNum);
          }
          int *boundProbLowerRowInd = new int[lowerRowNum];
          CoinIotaN(boundProbLowerRowInd, lowerRowNum, 0);
          int boundProbArgc = 1;
          char** boundProbArgv = new char* [1];
          boundProbArgv[0] = (char *) "mibs";

          //Setup an LP solver
          OsiSolverInterface *boundProbLpSolver = getSolver(boundProbLpSolverName, boundProbMaxThreads, false);
//          boundProbLpSolver->getModelPtr()->setDualBound(1.0e10);

          //New MibS model
          MibSModel *boundProbModel = new MibSModel();
          boundProbModel->setSolver(boundProbLpSolver);
          boundProbModel->AlpsPar()->setEntry(AlpsParams::msgLevel, -1);
          boundProbModel->AlpsPar()->setEntry(AlpsParams::timeLimit, 100);
          boundProbModel->MibSPar()->setEntry(MibSParams::cutStrategy, 0);

          boundProbModel->loadAuxiliaryData(lowerColNum,
                  lowerRowNum,
                  lowerColInd,
                  boundProbLowerRowInd,
                  lowerObjSense,
                  lowerObjCoef,
                  upperColNum,
                  0,
                  upperColInd,
                  NULL,
                  0, NULL,
                  0, NULL,
                  NULL, NULL);

          boundProbModel->loadProblemData(boundProbMat,
                  origColLb, origColUb,
                  boundProbUpperObjCoef,
                  boundProbRowLb, boundProbRowUb,
                  colType, 1.0, infinity,
                  boundProbRowSense);

#ifdef  COIN_HAS_MPI
          AlpsKnowledgeBrokerMPI boundProbBroker(boundProbArgc, boundProbArgv, *boundProbModel);
#else
          AlpsKnowledgeBrokerSerial boundProbBroker(boundProbArgc, boundProbArgv, *boundProbModel);
#endif

          boundProbBroker.search(boundProbModel);
          assert(boundProbBroker.getSolStatus() != AlpsExitStatusInfeasible);
          /*
          double *boundProbSolution;
          if (boundProbModel->getNumSolutions() > 0){
              boundProbSolution = boundProbModel->incumbent();
          }
          */

          dualBoundOnLevel2 = -1.0 * boundProbModel->getKnowledgeBroker()->getBestQuality();
          /*
          if (boundProbBroker.getBestNode()) {
              dualBoundOnLevel2 = -1.0 * boundProbBroker.getBestNode()->getQuality();
              assert(dualBoundOnLevel2 >= -etol);
          }
          */
      }
      //Value of boundOnLbf: 2 is a random multiplier
      if (fabs(primalBoundOnSubproblem) >= fabs(dualBoundOnSubproblem)) {
          boundOnLbf = 2*fabs(primalBoundOnSubproblem);
      } else {
          boundOnLbf = 2*fabs(dualBoundOnSubproblem);
      }
      //Also assert that upper level col bounds are finite!
      for (i = 0; i < upperColNum; i++) {
          assert(masterColLb[i] > -infinity);
          assert(masterColUb[i] < infinity);
      }
      //Best solution found so far
      int currentBestIteration = -1;
      double currentBestObjVal = infinity;
      double *currentBestSolution = new double[upperColNum + lowerColNum];
      CoinZeroN(currentBestSolution, upperColNum + lowerColNum);
      /*
      //Known solution
      double *optSol = new double[upperColNum];
      optSol[0] = 1;
      optSol[1] = 0;
      optSol[2] = 0;
      optSol[3] = 0;
      optSol[4] = 1;
      optSol[5] = 1;
      optSol[6] = 0;
      optSol[7] = 0;
      optSol[8] = 0;
      optSol[9] = 1;
//      optSol[0] = masterBestSolutionUpperColsPrevIter[0];
//      optSol[1] = masterBestSolutionUpperColsPrevIter[1];
//      optSol[2] = masterBestSolutionUpperColsPrevIter[2];
//      optSol[3] = masterBestSolutionUpperColsPrevIter[3];
//      optSol[4] = masterBestSolutionUpperColsPrevIter[4];
//      optSol[5] = masterBestSolutionUpperColsPrevIter[5];
//      optSol[6] = masterBestSolutionUpperColsPrevIter[6];
//      optSol[7] = masterBestSolutionUpperColsPrevIter[7];
//      optSol[8] = masterBestSolutionUpperColsPrevIter[8];
//      optSol[9] = masterBestSolutionUpperColsPrevIter[9];
      //      optSol[10] = -5620.75;
      */
      std::cout << std::endl;
      std::cout << "Starting algorithmic iterations..." << std::endl;
      std::cout << std::endl;
      std::cout << "    Iteration";
      std::cout << "  Upper Bound";
      std::cout << "  Lower Bound";
      std::cout << std::endl;



      /*** while loop for decomposition algorithm ***/
      while (!termFlag) {
          /*
          if (iterCounter >= 1) {
              // Checking if the known solution is already cut off from master problem
              // Note: This is for debugging purposes only.
              // Setting and solving the master problem
              masterBestSolution = new double[masterColNum];
              double masterObjValCopy = 0.0;
              double *masterColLbCopy = new double[masterColNum];
              double *masterColUbCopy = new double[masterColNum];
              memcpy(masterColLbCopy, masterColLb, sizeof(double)*masterColNum);
              memcpy(masterColUbCopy, masterColUb, sizeof(double)*masterColNum);
              memcpy(masterColLbCopy, optSol, sizeof(double)*(upperColNum));
              memcpy(masterColUbCopy, optSol, sizeof(double)*(upperColNum));
              solver = getSolver(masterProblemSolver, masterMaxThreads, false);
              masterInfeasible = solve(solver,
                      masterColNum, masterObjCoef, upperObjSense,
                      masterColLbCopy, masterColUbCopy, masterColType,
                      &masterMat, masterRowLb, masterRowUb,
                      &masterObjValCopy, masterBestSolution);
              std::cout << "********" << std::endl;
              std::cout << masterBestSolution[10] << "  " << 2836.33 << std::endl;
              std::cout << "********" << std::endl;
              assert((masterBestSolution[10] - 2836.333333333333) <= etol);
              delete solver;
              delete [] masterColUbCopy;
              delete [] masterColLbCopy;
              delete [] masterBestSolution;
          }
          */


          /** Setting up the subproblem in MILP form using master problem's solution **/
          if (!masterInfeasible) {
              //Setting subproblem Row LB, UB, RHS for the given masterBestSolution!
              matOfUpperCols.times(masterBestSolutionUpperCols, subproblemUpperColRowActivity);
              if (upperRowsHaveLowerCols) {
                  memcpy(subproblemRowLb, rowLb, sizeof(double)*(upperRowNum + lowerRowNum));
                  memcpy(subproblemRowUb, rowUb, sizeof(double)*(upperRowNum + lowerRowNum));
              } else {
                  memcpy(subproblemRowLb, &rowLb[upperRowNum], sizeof(double)*lowerRowNum);
                  memcpy(subproblemRowUb, &rowUb[upperRowNum], sizeof(double)*lowerRowNum);
              }
              for (i = 0; i < subproblemRowNum - 1; i++) {
                  if (subproblemRowLb[i] > -infinity) {
                      subproblemRowLb[i] -= subproblemUpperColRowActivity[i];
                      subproblemRhs[i] = subproblemRowLb[i];
                  } else if (subproblemRowUb[i] < infinity) {
                      subproblemRowUb[i] -= subproblemUpperColRowActivity[i];
                      subproblemRhs[i] = subproblemRowUb[i];
                  }
              }

              /* Solving second level MILP and gathering required data */
              //Building row bounds
              if (upperRowsHaveLowerCols) {
                  memcpy(level2RowLb, &subproblemRowLb[upperRowNum], sizeof(double)*lowerRowNum);
                  memcpy(level2RowUb, &subproblemRowUb[upperRowNum], sizeof(double)*lowerRowNum);
              } else {
                  memcpy(level2RowLb, subproblemRowLb, sizeof(double)*lowerRowNum);
                  memcpy(level2RowUb, subproblemRowUb, sizeof(double)*lowerRowNum);
              }
              //Actual solving
              solver = getSolver(level2ProblemSolver, level2MaxThreads, false);
              level2Infeasible = solve(solver,
                      lowerColNum, lowerObjCoef, lowerObjSense,
                      subproblemColLb, subproblemColUb, subproblemColType,
                      &level2Mat, level2RowLb, level2RowUb,
                      &level2ObjVal, level2BestSolution);
              if (!level2Infeasible) {
                  assert(solver->isProvenOptimal());
                  level2ObjVal = 0;
                  for (i = 0; i < lowerColNum; i++) {
                      if (subproblemColType[i] == 'I' || subproblemColType[i] == 'B') {
                          level2BestSolution[i] = round(level2BestSolution[i]);
                      }
                      level2ObjVal += level2BestSolution[i] * lowerObjCoef[i];
                  }
              }
              delete solver;

              if (!level2Infeasible) {
                  //Updating subproblem data with new info.
                  if (addRowInd) {
                      subproblemMat.appendRow(subproblemColNum, subproblemLowerColInd, lowerObjCoef);
                      addRowInd = false;
                  }
                  if (lowerObjSense == 1) {
                      subproblemRowLb[subproblemRowNum - 1] = -infinity;
                      subproblemRowUb[subproblemRowNum - 1] = level2ObjVal;
                  } else {
                      subproblemRowLb[subproblemRowNum - 1] = level2ObjVal;
                      subproblemRowUb[subproblemRowNum - 1] = infinity;
                  }
                  subproblemRhs[subproblemRowNum - 1] = level2ObjVal;
              } else {
                  if (!addRowInd) {
                      //Delete the objective bound-type row from subproblemMat
                      int rowIndToDelete = subproblemRowNum - 1;
                      subproblemMat.deleteRows(1, &rowIndToDelete);
                      addRowInd = true;
                  }
              }


              /** Solving the subproblem **/
              solver = getSolver(subproblemSolver, subproblemMaxThreads, true);
              subproblemInfeasible = solve(solver,
                      subproblemColNum, subproblemObjCoef, upperObjSense,
                      subproblemColLb, subproblemColUb, subproblemColType,
                      &subproblemMat, subproblemRowLb, subproblemRowUb,
                      &subproblemObjVal, subproblemBestSolution);

              //Updating reaction function value
              if (!subproblemInfeasible) {
                  subproblemObjVal = 0;
                  for (i = 0; i < subproblemColNum; i++) {
                      if (subproblemColType[i] == 'I' || subproblemColType[i] == 'B') {
                          subproblemBestSolution[i] = round(subproblemBestSolution[i]);
                      }
                      subproblemObjVal += subproblemBestSolution[i] * subproblemObjCoef[i];
                  }
                  rho = subproblemObjVal;
              }

              /* Checking termination criterion and gathering dual information if needed */
              //TODO: are the criteria correct?
              clock_t current = clock();
              double timeTillNow = (double) (current - begin) / CLOCKS_PER_SEC;
              timeUp = ((timeTillNow >= 14400) ? true : false);
              if (!subproblemInfeasible &&
                       (fabs(rho - rhoApproxValue) <= etol)
                    || timeUp) {
                 termFlag = true;
              } else {
                 //Getting dual information to the subproblem
                 getDualData(solver, dualBoundOnSubproblem,
                       &leafNodeNum, &leafFeasibilityStatusInd,
                       leafDualByRow, leafPosDjByRow, leafNegDjByRow,
                       &leafLbCnt, &leafLbInd, &leafLbVal,
                       &leafUbCnt, &leafUbInd, &leafUbVal);

                 //Changing sign of leafDualByRow w.r.t. orig. subproblem's row sense
                 //    since SYMPHONY changes all rows internally to 'L' sense
                 //FIXME: Note that there are no 'E' rows in our MIBLPs as of now.
                 for (i = 0; i < subproblemRowNum; i++) {
                    if (subproblemRowSense[i] == 'G') {
                       for (j = 0; j < leafNodeNum; j++) {
                          double element = leafDualByRow->getCoefficient(j, i);
                          if (element) {
                             leafDualByRow->modifyCoefficient(j, i, -element, true);
                          }
                       }
                    }
                 }
              }

              delete solver;

              if (!subproblemInfeasible) {
                 /* Check bilevel feasibility */
                 double level2ObjValForSubproblemSol = 0;
                 for (i = 0; i < lowerColNum; i++) {
                    level2ObjValForSubproblemSol += subproblemBestSolution[i]*lowerObjCoef[i];
                 }
                 assert(fabs(level2ObjValForSubproblemSol - level2ObjVal) <= etol);
                 //Updating original MIBLP objective value with its second level part
                 for (i = 0; i < lowerColNum; i++) {
                    optObjVal += subproblemObjCoef[i]*subproblemBestSolution[i];
                 }

                 /* Updating best solution found so far */
                 if (currentBestObjVal >= optObjVal + etol) {
                    currentBestIteration = iterCounter;
                    currentBestObjVal = optObjVal;
                    memcpy(currentBestSolution, masterBestSolutionUpperCols, sizeof(double)*upperColNum);
                    memcpy(&currentBestSolution[upperColNum], subproblemBestSolution, sizeof(double)*lowerColNum);
                 }
              }
          } else {
             //Master problem infeasible
             termFlag = true;
          }


          /** Generating Benders' cuts and adding them to the master problem **/
          if (!termFlag) {
              /* Finding various dual information products */
              //Partitioning leafDualByRow into two parts (original rows & extra row)
              CoinPackedMatrix leafDualToOrigRows(*leafDualByRow);
              double *fullDualOfExtraRow = new double[leafNodeNum];
              CoinZeroN(fullDualOfExtraRow, leafNodeNum);
              if (!level2Infeasible) {
                  int delColNum = 1, minorDim = leafDualToOrigRows.getMinorDim();
                  int *delColInd = new int[delColNum];
                  delColInd[0] = subproblemRowNum - 1;
                  if ((subproblemRowNum - 1) < minorDim) {
                      leafDualToOrigRows.deleteCols(delColNum, delColInd);
                  }

                  CoinPackedMatrix copyLeafDualMat(*leafDualByRow);
                  copyLeafDualMat.reverseOrdering();
                  int majorDim = copyLeafDualMat.getMajorDim();
                  if ((subproblemRowNum - 1) < majorDim) {
                      CoinShallowPackedVector leafDualOfExtraRow = copyLeafDualMat.getVector(subproblemRowNum - 1);
                      int leafDualOfExtraRowNum = leafDualOfExtraRow.getNumElements();
                      const int *leafDualOfExtraRowInd = leafDualOfExtraRow.getIndices();
                      const double *leafDualOfExtraRowVal = leafDualOfExtraRow.getElements();
                      for (i = 0; i < leafDualOfExtraRowNum; i++) {
                          fullDualOfExtraRow[leafDualOfExtraRowInd[i]] = leafDualOfExtraRowVal[i];
                      }
                  }
                  delete [] delColInd;
              }

              //Setting feasible leaf node number and indices
              feasibleLeafNodeNum = 0;
              //Actual size should be feasibleLeafNodeNum
              feasibleLeafNodeInd = new int[leafNodeNum];
              for (i = 0; i < leafNodeNum; i++) {
                  if (leafFeasibilityStatusInd[i]) {
                      feasibleLeafNodeInd[feasibleLeafNodeNum] = i;
                      feasibleLeafNodeNum++;
                  }
              }

              //START finding products
              //Product of constraint matrix (A^2 or (A^1; A^2)) and dual info.
              double **product1 = new double*[leafNodeNum];
              CoinShallowPackedVector singleDualRow;
              int singleDualNnz, majorDimDual = leafDualToOrigRows.getMajorDim();
              int allDualNnz = leafDualToOrigRows.getNumElements();
              const int *singleDualInd;
              const double *singleDualVal;
              //For building full dual
              //NOTE: Actual size will be 'subproblemRowNum-1' if level2 is infeasible
              double *singleDual = new double[subproblemRowNum];

              //Products of (column LBs and positive reduced costs (leafPosDjByRow))
              //    and (column UBs and negative reduced costs (leafNegDjByRow))
              double *lbPosDjProduct = new double[leafNodeNum];
              double *ubNegDjProduct = new double[leafNodeNum];
              //TODO: Is following initialization necessary?
              CoinZeroN(lbPosDjProduct, leafNodeNum);
              CoinZeroN(ubNegDjProduct, leafNodeNum);
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
                  product1[i] = new double[upperColNum];
                  CoinZeroN(product1[i], upperColNum);
                  //if nonzero dual entries exist, find the product
                  if (allDualNnz && (i < majorDimDual)) {
                      singleDualRow = leafDualToOrigRows.getVector(i);
                      singleDualNnz = singleDualRow.getNumElements();
                      singleDualInd = singleDualRow.getIndices();
                      singleDualVal = singleDualRow.getElements();
                      //Building a full vector from nonzeroes
                      CoinZeroN(singleDual, subproblemRowNum);
                      for (j = 0; j < singleDualNnz; j++) {
                          singleDual[singleDualInd[j]] = singleDualVal[j];
                      }
                      //Product
                      matOfUpperCols.transposeTimes(singleDual, product1[i]);
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
                          if (leafUbVal[i][j] < tempSubproblemColUb[leafUbInd[i][j]]) {
                              tempSubproblemColUb[leafUbInd[i][j]] = leafUbVal[i][j];
                          }
                      }

                      //Product of bounds and reduced costs
                      for (j = 0; j < negDjNum; j++) {
                          ubNegDjProduct[i] += negDjElements[j]*tempSubproblemColUb[negDjIndices[j]];
                      }
                  }
              }

              //Product of leafDualToOrigRows and subproblemOrigRhs
              double *product7 = new double[leafNodeNum];
              CoinZeroN(product7, leafNodeNum);
              leafDualToOrigRows.times(subproblemOrigRhs, product7);

              //Declare remaining products but calculate them only if required
              //Product of constraint matrix (A^2) and dual of continuous restriction
              double *product2 = new double[upperColNum];
              //Product of dual of continuous restriction and lower level's
              //    row activity of integer restriction
              double product3 = 0;
              //Product of cont. rest. basis inverse and constraint matrix A^2
              double **product4 = new double*[lowerRowNum];
              //Product of cont. rest. basis inverse and lower level's row
              //    activity of integer restriction
              double *product5 = new double[lowerRowNum];
              //Product of dual of cont. rest. and lowerRowRhs
              double product6 = 0;
              //Product of cont. rest.'s dj and col. bounds
              double product8 = 0;
              //Product of cont. rest. non-basic part of matrix and non-basic variables
              double *product9 = new double[lowerRowNum];
              //Product of cont. rest. basis inverse and product9
              double *product10 = new double[lowerRowNum];
              if (!level2Infeasible) {
                  /* Continuous restriction building and solving */
                  //If level2 problem feasible, solve the continuous restriction
                  //    for the known integer component of subproblemBestSolution
                  //    (or level2BestSolution if subproblem is infeasible).
                  for (i = 0; i < lowerIntColNum; i++) {
                      intBestSolution[i] = (!subproblemInfeasible ?
                            subproblemBestSolution[lowerIntColInd[i]] :
                            level2BestSolution[lowerIntColInd[i]]);
                  }

                  //Finding product of integer best solution and corresponding obj. vector
                  level2IntObjVal = 0;
                  for (i = 0; i < lowerIntColNum; i++) {
                      level2IntObjVal += intBestSolution[i]*lowerObjCoef[lowerIntColInd[i]];
                  }

                  //Finding product of integer best solution and integer matrix
                  intRestMat.times(intBestSolution, level2IntColRowActivity);

                  //Finding RowLb and RowUb for continuous restriction
                  memcpy(contRestRowLb, level2RowLb, sizeof(double)*lowerRowNum);
                  memcpy(contRestRowUb, level2RowUb, sizeof(double)*lowerRowNum);
                  for (i = 0; i < lowerRowNum; i++) {
                      if (contRestRowLb[i] > -infinity) {
                          contRestRowLb[i] -= level2IntColRowActivity[i];
                          //Following line to avoid soplex related errors
                          contRestRowLb[i] = contRestRowLb[i];
                          //Following line because all are rows are '=' type now!
                          contRestRowUb[i] = contRestRowLb[i];
                      } else if (contRestRowUb[i] < infinity) {
                          contRestRowUb[i] -= level2IntColRowActivity[i];
                          //Following line to avoid soplex related errors
                          contRestRowUb[i] = contRestRowUb[i];
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
                  memcpy(contRestDjSolution, solver->getReducedCost(), sizeof(double)*contRestColNum);
                  soplex::SoPlex *soplex =
                      dynamic_cast<OsiSpxSolverInterface*>(solver)->getLpPtr();
                  for (i = 0; i < lowerRowNum; i++) {
                      soplex::SSVectorRational basisInverseCol(lowerRowNum);
                      soplex->getBasisInverseColRational(i, basisInverseCol);
                      for (j = 0; j < lowerRowNum; j++) {
                          contRestBasisInverseRows.setValue(i + j*lowerRowNum, basisInverseCol[j]);
                      }
                  }
                  for (i = 0; i < lowerRowNum; i++) {
                      soplex::SSVectorRational basisInverseRow(lowerRowNum);
                      for (j = 0; j < lowerRowNum; j++) {
                          basisInverseRow.add(j, contRestBasisInverseRows[i*lowerRowNum + j]);
                      }
                      contRestBasisInverseRowLcm[i] = dlcmRational(basisInverseRow.values(), lowerRowNum);
                  }
                  for (i = 0; i < lowerRowNum; i++) {
                      contRestBasisInverseRow[i] = new double[lowerRowNum];
                      solver->getBInvRow(i, contRestBasisInverseRow[i]);
                  }
                  solver->getBasics(contRestBasisIndices);
                  //Changing contRestBasisIndices to reflect slack/surplus variables based on output from SoPlex
                  //We get the indices as "numcols -1 + rownum" in case of a basic row
                  //    ==> the slack variables corresponding to that row is a basic col
                  //        ==> we can mark our already introduced slack col as basic
                  for (i = 0; i < lowerRowNum; i++) {
                      if (contRestBasisIndices[i] >= contRestColNum) {
                          contRestBasisIndices[i] = contRestBasisIndices[i] - contRestColNum + lowerContColNum;
                      }
                  }
                  //Building non-basic part of opt. sol. by zeroing basic part
                  memcpy(contRestBestSolutionNonBasics, contRestBestSolution, sizeof(double)*contRestColNum);
                  for (i = 0; i < lowerRowNum; i++) {
                      contRestBestSolutionNonBasics[contRestBasisIndices[i]] = 0;
                  }
                  delete solver;

                  if (contRestInfeasible) {
                      //Should not happen because subproblemInfeasible is false!
                      std::cout <<
                          "Error: Restriction is infeasible whereas it should not be!."
                          << std::endl;
                      return 0;
                  }

                  /* Back to product calculations */
                  //Product of constraint matrix (A^2) and dual of continuous restriction
                  CoinZeroN(product2, upperColNum);
                  lowerMatOfUpperCols.transposeTimes(contRestDualSolution, product2);

                  //Product of dual of continuous restriction and lower level's
                  //    row activity of integer restriction
                  for (i = 0; i < lowerRowNum; i++) {
                      product3 += contRestDualSolution[i]*level2IntColRowActivity[i];
                  }

                  //Product of cont. rest. basis inverse and constraint matrix A^2
                  for (i = 0; i < lowerRowNum; i++) {
                      product4[i] = new double[upperColNum];
                      CoinZeroN(product4[i], upperColNum);
                      lowerMatOfUpperCols.transposeTimes(contRestBasisInverseRow[i], product4[i]);

                      //bigM for domain restriction of first-level variable (x)
                      maxValForDomainRest[i] = 0;
                      minValForDomainRest[i] = 0;
                      for (j = 0; j < upperColNum; j++) {
                          //Multiplying product4[i] with i-th row's LCM
                          product4[i][j] *= contRestBasisInverseRowLcm[i];
                          product4[i][j] = round(product4[i][j]);
                          double coef = product4[i][j];
                          if (coef < -etol) {
                              maxValForDomainRest[i] += coef*masterColLb[j];
                              minValForDomainRest[i] += coef*masterColUb[j];
                          } else if (coef > etol) {
                              maxValForDomainRest[i] += coef*masterColUb[j];
                              minValForDomainRest[i] += coef*masterColLb[j];
                          }
                      }
                  }

                  //Product of cont. rest. basis inverse and lower level's row
                  //    activity of integer restriction
                  for (i = 0; i < lowerRowNum; i++) {
                      product5[i] = 0;
                      for (j = 0; j < lowerRowNum; j++) {
                          product5[i] += contRestBasisInverseRowLcm[i]*contRestBasisInverseRow[i][j]*
                              (level2IntColRowActivity[j] - lowerRowRhs[j]);
                      }
                      product5[i] = round(product5[i]);
                  }

                  //Product of dual of cont. rest. and lowerRowRhs
                  for (i = 0; i < lowerRowNum; i++) {
                      product6 += contRestDualSolution[i]*lowerRowRhs[i];
                  }

                  //Product of cont. rest.'s dj and col. bounds
                  for (i = 0; i < contRestColNum; i++) {
                      if (contRestDjSolution[i] > etol) {
                          assert(contRestColLb[i] > -infinity);
                          product8 += contRestDjSolution[i]*contRestColLb[i];
                      } else if (contRestDjSolution[i] < -etol) {
                          assert(contRestColUb[i] < infinity);
                          product8 += contRestDjSolution[i]*contRestColUb[i];
                      }
                  }

                  //Asserting equality of 2nd level objective value and its various components
                  double oneComponent = 0.0;
                  for (i = 0; i < upperColNum; i++) {
                      oneComponent += product2[i]*masterBestSolutionUpperCols[i];
                  }
                  assert(fabs(product6 - oneComponent - product3 + level2IntObjVal + product8 - level2ObjVal) <= etol);

                  //Product of cont. rest. non-basic part of matrix and non-basic variables
                  CoinZeroN(product9, lowerRowNum);
                  contRestMat.times(contRestBestSolutionNonBasics, product9);
                  for (i = 0; i < lowerRowNum; i++) {
                      product9[i] = round(product9[i]);
                  }

                  //Product of cont. rest. basis inverse and product9
                  for (i = 0; i < lowerRowNum; i++) {
                      product10[i] = 0;
                      for (j = 0; j < lowerRowNum; j++) {
                          product10[i] += contRestBasisInverseRowLcm[i]*contRestBasisInverseRow[i][j]*
                              product9[j];
                      }
                      product10[i] = round(product10[i]);
                  }

                  //bigM for UBF of level2 problem
                  double maxVal = 0;
                  double minVal = 0;
                  for (i = 0; i < upperColNum; i++) {
                      double coef = product2[i];
                      if (coef < -etol) {
                          maxVal += coef*masterColLb[i];
                          minVal += coef*masterColUb[i];
                      } else if (coef > etol) {
                          maxVal += coef*masterColUb[i];
                          minVal += coef*masterColLb[i];
                      }
                  }
                  double ubfMin = product6 - product3 + level2IntObjVal + product8 - maxVal;
                  double ubfMax = product6 - product3 + level2IntObjVal + product8 - minVal;
                  /*
                  //Approach-1
                  bigMForUbf = fabs(dualBoundOnLevel2) + fabs(ubfMin);
                  */
                  /*
                  //Approach-2 (incorrect)
                  bigMForUbf = fabs(dualBoundOnLevel2 - ubfMax);
                  */
                  //Approach-3
                  bigMForUbf = fabs(dualBoundOnLevel2 - ubfMin);

                  //bigM for LBF of subproblem
                  singleBigMForLbf = 0;
                  double tempMax = -infinity, tempMin = infinity;
                  for (i = 0; i < feasibleLeafNodeNum; i++) {
                      double bigMForLbf = 0, minVal1 = 0, maxVal1 = 0;
                      double minVal2 = 0, maxVal2 = 0;
                      for (j = 0; j < upperColNum; j++) {
                          double coef = (product1[feasibleLeafNodeInd[i]][j] +
                                      fullDualOfExtraRow[feasibleLeafNodeInd[i]]*product2[j]);
                          if (coef < -etol) {
                              minVal1 += coef*masterColUb[j];
                              maxVal1 += coef*masterColLb[j];
                          } else if (coef > etol) {
                              minVal1 += coef*masterColLb[j];
                              maxVal1 += coef*masterColUb[j];
                          }
                      }
                      if (fullDualOfExtraRow[feasibleLeafNodeInd[i]] < -etol) {
                          minVal2 = bigMForUbf*fullDualOfExtraRow[feasibleLeafNodeInd[i]];
                          maxVal2 = 0;
                      } else if (fullDualOfExtraRow[feasibleLeafNodeInd[i]] > etol) {
                          minVal2 = 0;
                          maxVal2 = bigMForUbf*fullDualOfExtraRow[feasibleLeafNodeInd[i]];
                      }
                      bigMForLbf = product7[feasibleLeafNodeInd[i]] +
                               lbPosDjProduct[feasibleLeafNodeInd[i]] +
                               ubNegDjProduct[feasibleLeafNodeInd[i]] +
                               fullDualOfExtraRow[feasibleLeafNodeInd[i]]*(product6 - product3 + level2IntObjVal + product8);
                      if (tempMax < (bigMForLbf - minVal1 + maxVal2)) {
                          tempMax = (bigMForLbf - minVal1 + maxVal2);
                      }
                      if (tempMin > (bigMForLbf - maxVal1 + minVal2)) {
                          tempMin = (bigMForLbf - maxVal1 + minVal2);
                      }
                      /*
                      //Update Feb. 7, 2019: the following two approaches seem
                      //    to be incorrect.
                      //Approach-1
                      //Note: 4 is a random multiplier
                      bigMForLbf = ((tempMax >= tempMin + etol) ? 4*tempMax : 4*tempMin);
                      //Approach-2: Difference between max and min values
                      bigMForLbf = tempMax - tempMin;
                      assert(bigMForLbf >= -etol);
                      if ((bigMForLbf - singleBigMForLbf) > etol) {
                          singleBigMForLbf = bigMForLbf;
                      }
                      */
                  }
                  if ((tempMax > etol) && (tempMin < -etol)) {
                      //Note: 2 is a random multiplier
                      singleBigMForLbf = 2*(tempMax - tempMin);
                  } else if (tempMax > etol) {
                      singleBigMForLbf = 2*tempMax;
                  } else {
                      singleBigMForLbf = boundOnLbf;
                  }

                  /*
                  //The following case is accounted just above.
                  if (singleBigMForLbf <= etol) {
                      //bigM = 0; make is some nonzero
                      singleBigMForLbf = boundOnLbf;
                  }
                  */
              } else {
                  //bigM for LBF of subproblem
                  singleBigMForLbf = 0;
                  double tempMax = -infinity, tempMin = infinity;
                  for (i = 0; i < feasibleLeafNodeNum; i++) {
                      double bigMForLbf = 0, minVal = 0, maxVal = 0;
                      for (j = 0; j < upperColNum; j++) {
                          double coef = product1[feasibleLeafNodeInd[i]][j];
                          if (coef < -etol) {
                              minVal += coef*masterColUb[j];
                              maxVal += coef*masterColLb[j];
                          } else if (coef > etol) {
                              minVal += coef*masterColLb[j];
                              maxVal += coef*masterColUb[j];
                          }
                      }
                      bigMForLbf = product7[feasibleLeafNodeInd[i]] +
                               lbPosDjProduct[feasibleLeafNodeInd[i]] +
                               ubNegDjProduct[feasibleLeafNodeInd[i]];
                      if (tempMax < (bigMForLbf - minVal)) {
                          tempMax = (bigMForLbf - minVal);
                      }
                      if (tempMin > (bigMForLbf - maxVal)) {
                          tempMin = (bigMForLbf - maxVal);
                      }
                      /*
                      //The following approach seems to be incorrect.
                      //Note: 4 is a random multiplier
                      bigMForLbf = ((tempMax >= tempMin + etol) ? 4*tempMax : 4*tempMin);
                      if ((bigMForLbf - singleBigMForLbf) > etol) {
                          singleBigMForLbf = bigMForLbf;
                      }
                      */
                  }
                  if ((tempMax > etol) && (tempMin < -etol)) {
                      //Note: 2 is a random multiplier
                      singleBigMForLbf = 2*(tempMax - tempMin);
                  } else if (tempMax > etol) {
                      singleBigMForLbf = 2*tempMax;
                  } else {
                      singleBigMForLbf = boundOnLbf;
                  }

                  /*
                  //The following case is accounted just above.
                  if (singleBigMForLbf <= etol) {
                      //bigM = 0; make is some nonzero
                      singleBigMForLbf = boundOnLbf;
                  }
                  */
              }
              //END finding products

              /* Updating master problem's matrices and vectors/arrays */
              //Finding number of extra binary variables for domain restriction on 'x'
              numBinColsForDomainRest = 0;
              //FIXME: combine the code of all for loops of the following kind
              if (!level2Infeasible) {
                  CoinZeroN(tol, 2*lowerRowNum);
                  for (i = 0; i < lowerRowNum; i++) {
                      int ind = contRestBasisIndices[i];
                      if (ind < lowerContColNum) {
                          if (lowerContColFiniteLbId[ind]) {
                              //Finding "tol"
                              //Setting data matrices and vectors
                              memcpy(tolProbObjCoef, product4[i], sizeof(double)*upperColNum);
                              tolProbRowLb[0] = -product5[i] - product10[i] - contRestBasisInverseRowLcm[i]*contRestColLb[ind] + 1;
                              for (j = 0; j < upperColNum; j++) {
                                  tolProbMat->modifyCoefficient(0, j, product4[i][j], true);
                              }
                              tolProbObjVal = 0.0;
                              //Actual solving
                              solver = getSolver(tolProbSolver, tolProbMaxThreads, false);
                              tolProbInfeasible = solve(solver,
                                      tolProbColNum, tolProbObjCoef, tolProbObjSense,
                                      tolProbColLb, tolProbColUb, tolProbColType,
                                      tolProbMat, tolProbRowLb, tolProbRowUb,
                                      &tolProbObjVal, tolProbBestSolution);
                              if (tolProbInfeasible) {
                                  //1 is some random value
                                  tol[numBinColsForDomainRest] = 1.0;
                              } else {
                                  assert(solver->isProvenOptimal());
                                  tolProbObjVal = 0;
                                  for (j = 0; j < tolProbColNum; j++) {
                                      //Following rounding if valid due to the known fact that all cols are binary
                                      tolProbBestSolution[j] = round(tolProbBestSolution[j]);
                                      tolProbObjVal += tolProbBestSolution[j]*tolProbObjCoef[j];
                                  }
                                  double tolTemp = fabs(-product5[i] - product10[i]
                                          - contRestBasisInverseRowLcm[i]*contRestColLb[ind] - tolProbObjVal);
                                  assert(tolTemp > etol); // not equal to zero
                                  //FIXME: 'floor' is used according to 'tol' usage
                                  //  later in the algo. Check again if 'ceil' is reqd.
                                  tol[numBinColsForDomainRest] = floor(tolTemp);
                              }
                              delete solver;

                              //Counting # of binary variables for domain rest.
                              numBinColsForDomainRest++;
                          }
                          if (lowerContColFiniteUbId[ind]) {
                              //Finding "tol"
                              //Setting data matrices and vectors
                              tolProbRowLb[0] = product5[i] + product10[i] + contRestBasisInverseRowLcm[i]*contRestColUb[ind] + 1;
                              for (j = 0; j < upperColNum; j++) {
                                  tolProbObjCoef[j] = -product4[i][j];
                                  tolProbMat->modifyCoefficient(0, j, -product4[i][j], true);
                              }
                              tolProbObjVal = 0.0;
                              //Actual solving
                              solver = getSolver(tolProbSolver, tolProbMaxThreads, false);
                              tolProbInfeasible = solve(solver,
                                      tolProbColNum, tolProbObjCoef, tolProbObjSense,
                                      tolProbColLb, tolProbColUb, tolProbColType,
                                      tolProbMat, tolProbRowLb, tolProbRowUb,
                                      &tolProbObjVal, tolProbBestSolution);
                              if (tolProbInfeasible) {
                                  //1 is some random value
                                  tol[numBinColsForDomainRest] = 1.0;
                              } else {
                                  assert(solver->isProvenOptimal());
                                  tolProbObjVal = 0;
                                  for (j = 0; j < tolProbColNum; j++) {
                                      //Following rounding if valid due to the known fact that all cols are binary
                                      tolProbBestSolution[j] = round(tolProbBestSolution[j]);
                                      tolProbObjVal += tolProbBestSolution[j]*tolProbObjCoef[j];
                                  }
                                  double tolTemp = fabs(product5[i] + product10[i]
                                          + contRestBasisInverseRowLcm[i]*contRestColUb[ind] - tolProbObjVal);
                                  assert(tolTemp > etol); // not equal to zero
                                  //FIXME: 'floor' is used according to 'tol' usage
                                  //  later in the algo. Check again if 'ceil' is reqd.
                                  tol[numBinColsForDomainRest] = floor(tolTemp);
                              }
                              delete solver;

                              //Counting # of binary variables for domain rest.
                              numBinColsForDomainRest++;
                          }
                      } else {
                          //Setting data matrices and vectors
                          memcpy(tolProbObjCoef, product4[i], sizeof(double)*upperColNum);
                          tolProbRowLb[0] = -product5[i] - product10[i] + 1;
                          for (j = 0; j < upperColNum; j++) {
                              tolProbMat->modifyCoefficient(0, j, product4[i][j], true);
                          }
                          tolProbObjVal = 0.0;
                          //Actual solving
                          solver = getSolver(tolProbSolver, tolProbMaxThreads, false);
                          tolProbInfeasible = solve(solver,
                                  tolProbColNum, tolProbObjCoef, tolProbObjSense,
                                  tolProbColLb, tolProbColUb, tolProbColType,
                                  tolProbMat, tolProbRowLb, tolProbRowUb,
                                  &tolProbObjVal, tolProbBestSolution);
                          if (tolProbInfeasible) {
                              //1 is some random value
                              tol[numBinColsForDomainRest] = 1;
                          } else {
                              assert(solver->isProvenOptimal());
                              tolProbObjVal = 0;
                              for (j = 0; j < tolProbColNum; j++) {
                                  //Following rounding if valid due to the known fact that all cols are binary
                                  tolProbBestSolution[j] = round(tolProbBestSolution[j]);
                                  tolProbObjVal += tolProbBestSolution[j]*tolProbObjCoef[j];
                              }
                              double tolTemp = fabs(-product5[i] - product10[i] - tolProbObjVal);
                              assert(tolTemp > etol); // not equal to zero
                              //FIXME: 'floor' is used according to 'tol' usage
                              //  later in the algo. Check again if 'ceil' is reqd.
//                              std::cout << tolTemp << "\t" << floor(tolTemp) << "\t" << ceil(tolTemp) << "\t" << round(tolTemp) << std::endl;
                              tol[numBinColsForDomainRest] = floor(tolTemp);
                          }
                          delete solver;

                          //Counting # of binary variables for domain rest.
                          numBinColsForDomainRest++;
                      }
                  }
              }
              //FIXME: Technically, follwing if-else conditions should depend on presence
              //    of continuous variables in 2nd level prob, i.e., if UBF is generated or not!
              /* NOTE: if level2Infeasible = false:
                         new # of cols = old # of cols + feasibleLeafNodeNum + 1 + numBinColsForDomainRest
                         new # of rows = old # of rows + feasibleLeafNodeNum + 1 + 2*numBinColsForDomainRest + 2
                       else if level2Infeasible = true:
                         new # of cols = old # of cols + feasibleLeafNodeNum
                         new # of rows = old # of rows + feasibleLeafNodeNum + 1 */

              //FIXME: Shall we remove following if-check because else-check is covered WLOG?
              if (!level2Infeasible) {
                  //Objective coefficient vector
                  masterObjCoefVec.resize((masterColNum + feasibleLeafNodeNum + 1 + numBinColsForDomainRest), 0);

                  //Column bounds
                  masterColLbVec.resize((masterColNum + feasibleLeafNodeNum + 1 + numBinColsForDomainRest), 0);
                  masterColUbVec.resize((masterColNum + feasibleLeafNodeNum + 1 + numBinColsForDomainRest), 1);

                  //Column types
                  masterColTypeVec.resize((masterColNum + feasibleLeafNodeNum + 1 + numBinColsForDomainRest), 'B');

                  //Row bounds
                  for (i = 0; i < feasibleLeafNodeNum; i++) {
                      masterRowLbVec.resize((masterRowNum + (i+1)),
                              (product7[feasibleLeafNodeInd[i]] +
                               lbPosDjProduct[feasibleLeafNodeInd[i]] +
                               ubNegDjProduct[feasibleLeafNodeInd[i]] +
                               fullDualOfExtraRow[feasibleLeafNodeInd[i]]*(product6 - product3 + level2IntObjVal + product8) -
                               singleBigMForLbf));
                  }
                  masterRowLbVec.resize((masterRowNum + feasibleLeafNodeNum + 1), 1);
                  int counterTemp = 0;
                  for (i = 0; i < lowerRowNum; i++) {
                      int ind = contRestBasisIndices[i];
                      if (ind < lowerContColNum) {
                          if (lowerContColFiniteLbId[ind]) {
                              masterRowLbVec.resize((masterRowNum + feasibleLeafNodeNum + 1 + (counterTemp+1)),
                                      (product5[i] + product10[i] + contRestBasisInverseRowLcm[i]*contRestColLb[ind]));
                              counterTemp++;
                          }
                          if (lowerContColFiniteUbId[ind]) {
                              masterRowLbVec.resize((masterRowNum + feasibleLeafNodeNum + 1 + (counterTemp+1)),
                                      (-product5[i] - product10[i] - contRestBasisInverseRowLcm[i]*contRestColUb[ind]));
                              counterTemp++;
                          }
                      } else {
                          masterRowLbVec.resize((masterRowNum + feasibleLeafNodeNum + 1 + (counterTemp+1)), product5[i] + product10[i]);
                          counterTemp++;
                      }
                  }
                  assert(counterTemp == numBinColsForDomainRest);
                  masterRowLbVec.resize((masterRowNum + feasibleLeafNodeNum + 1 + 2*numBinColsForDomainRest), -infinity);
                  masterRowLbVec.resize((masterRowNum + feasibleLeafNodeNum + 1 + 2*numBinColsForDomainRest + 1), 0);
                  masterRowLbVec.resize((masterRowNum + feasibleLeafNodeNum + 1 + 2*numBinColsForDomainRest + 2), -infinity);

                  masterRowUbVec.resize((masterRowNum + feasibleLeafNodeNum), infinity);
                  masterRowUbVec.resize((masterRowNum + feasibleLeafNodeNum + 1), 1);
                  masterRowUbVec.resize((masterRowNum + feasibleLeafNodeNum + 1 + numBinColsForDomainRest), infinity);
                  counterTemp = 0;
                  for (i = 0; i < lowerRowNum; i++) {
                      int ind = contRestBasisIndices[i];
                      if (ind < lowerContColNum) {
                          if (lowerContColFiniteLbId[ind]) {
                              double bigMForDomainRest = (-product5[i] - product10[i] - minValForDomainRest[i] -
                                      contRestBasisInverseRowLcm[i]*contRestColLb[ind]);
                              // Offset bigM by epsilon (= tol[counterTemp])
                              bigMForDomainRest += tol[counterTemp];
                              /*
                              if (bigMForDomainRest <= etol) {
                                  // bigM = 0; 10 is a random number!
                                  bigMForDomainRest = 10;
                              }
                              */
                              masterRowUbVec.resize((masterRowNum + feasibleLeafNodeNum + 1 + numBinColsForDomainRest +
                                          (counterTemp+1)), (bigMForDomainRest + product5[i] + product10[i] +
                                              contRestBasisInverseRowLcm[i]*contRestColLb[ind] - tol[counterTemp]));
                              counterTemp++;
                          }
                          if (lowerContColFiniteUbId[ind]) {
                              double bigMForDomainRest = (product5[i] + product10[i] + maxValForDomainRest[i] +
                                      contRestBasisInverseRowLcm[i]*contRestColUb[ind]);
                              // Offset bigM by epsilon (= tol[counterTemp])
                              bigMForDomainRest += tol[counterTemp];
                              /*
                              if (bigMForDomainRest <= etol) {
                                  // bigM = 0; 10 is a random number!
                                  bigMForDomainRest = 10;
                              }
                              */
                              masterRowUbVec.resize((masterRowNum + feasibleLeafNodeNum + 1 + numBinColsForDomainRest +
                                          (counterTemp+1)), (bigMForDomainRest - product5[i] - product10[i] -
                                              contRestBasisInverseRowLcm[i]*contRestColUb[ind] - tol[counterTemp]));
                              counterTemp++;
                          }
                      } else {
                          double bigMForDomainRest = (-product5[i] - product10[i] - minValForDomainRest[i]);
                          // Offset bigM by epsilon (= tol[counterTemp])
                          bigMForDomainRest += tol[counterTemp];
                          /*
                          if (bigMForDomainRest <= etol) {
                              // bigM = 0; 10 is a random number!
                              bigMForDomainRest = 10;
                          }
                          */
                          masterRowUbVec.resize((masterRowNum + feasibleLeafNodeNum + 1 + numBinColsForDomainRest +
                                      (counterTemp+1)), (bigMForDomainRest + product5[i] + product10[i] - tol[counterTemp]));
                          counterTemp++;
                      }
                  }
                  assert(counterTemp == numBinColsForDomainRest);
                  masterRowUbVec.resize((masterRowNum + feasibleLeafNodeNum + 1 + 2*numBinColsForDomainRest + 1), infinity);
                  masterRowUbVec.resize((masterRowNum + feasibleLeafNodeNum + 1 + 2*numBinColsForDomainRest + 2), 0);
              } else {
                  //Objective coefficient vector
                  masterObjCoefVec.resize((masterColNum + feasibleLeafNodeNum), 0);

                  //Column bounds
                  masterColLbVec.resize((masterColNum + feasibleLeafNodeNum), 0);
                  masterColUbVec.resize((masterColNum + feasibleLeafNodeNum), 1);

                  //Column types
                  masterColTypeVec.resize((masterColNum + feasibleLeafNodeNum), 'B');

                  //Row bounds
                  for (i = 0; i < feasibleLeafNodeNum; i++) {
                      masterRowLbVec.resize((masterRowNum + (i+1)),
                              (product7[feasibleLeafNodeInd[i]] +
                               lbPosDjProduct[feasibleLeafNodeInd[i]] +
                               ubNegDjProduct[feasibleLeafNodeInd[i]] -
                               singleBigMForLbf));
                  }
                  masterRowLbVec.resize((masterRowNum + feasibleLeafNodeNum + 1), 1);

                  masterRowUbVec.resize((masterRowNum + feasibleLeafNodeNum), infinity);
                  masterRowUbVec.resize((masterRowNum + feasibleLeafNodeNum + 1), 1);
              }

              masterObjCoef = &masterObjCoefVec[0];
              masterColLb = &masterColLbVec[0];
              masterColUb = &masterColUbVec[0];
              masterColType = &masterColTypeVec[0];
              masterRowLb = &masterRowLbVec[0];
              masterRowUb = &masterRowUbVec[0];

              //Row coefficient matrix
              //First, appending new columns to existing matrix
              int newColNum = (level2Infeasible ? (feasibleLeafNodeNum) :
                      (feasibleLeafNodeNum + 1 + numBinColsForDomainRest));
              CoinBigIndex *newColStarts = new CoinBigIndex[newColNum + 1];
              CoinZeroN(newColStarts, newColNum + 1);
              int errorNum = masterMat.appendCols(newColNum, newColStarts,
                      NULL, NULL, masterRowNum);
              assert(errorNum == 0);

              //Now, appending new rows to the matrix
              //FIXME: Simplify further by removing common rows outside if-else condition.
              if (!level2Infeasible) {
                  for (i = 0; i < feasibleLeafNodeNum; i++) {
                      CoinPackedVector row;
                      for (j = 0; j < upperColNum; j++) {
                          row.insert(j, (product1[feasibleLeafNodeInd[i]][j] +
                                      fullDualOfExtraRow[feasibleLeafNodeInd[i]]*product2[j]));
                      }
                      row.insert(upperColNum, 1);
                      row.insert((masterColNum + i), -singleBigMForLbf);
                      if (fullDualOfExtraRow[feasibleLeafNodeInd[i]] >= -etol) {
                          //dual value = 0 (since we know it is <= 0 to begin with)
                          //   so, following is equivalent to dual = -1 to still have effect to bigMForUbf
                          row.insert((masterColNum + feasibleLeafNodeNum), bigMForUbf);
                      } else {
                          row.insert((masterColNum + feasibleLeafNodeNum), -bigMForUbf*fullDualOfExtraRow[feasibleLeafNodeInd[i]]);
                      }
                      masterMat.appendRow(row);
                  }
                  CoinPackedVector oneMoreRow;
                  for (i = 0; i < feasibleLeafNodeNum; i++) {
                      oneMoreRow.insert((masterColNum + i), 1);
                  }
                  masterMat.appendRow(oneMoreRow);
                  int counterTemp = 0;
                  for (i = 0; i < lowerRowNum; i++) {
                      int ind = contRestBasisIndices[i];
                      if (ind < lowerContColNum) {
                          if (lowerContColFiniteLbId[ind]) {
                              double bigMForDomainRest = -product5[i] - product10[i]- maxValForDomainRest[i] -
                                  contRestBasisInverseRowLcm[i]*contRestColLb[ind];
                              //Note: 10 is a random multiplier
                              if (bigMForDomainRest > etol) {
//                                  bigMForDomainRest = -10*bigMForDomainRest;
                              } else if (bigMForDomainRest < etol) {
//                                  bigMForDomainRest = 10*bigMForDomainRest;
                              } else {
                                  // bigM = 0; 1 is a random positive number to not let v_i = 1 to occur!
                                  bigMForDomainRest = 1;
                              }
                              CoinPackedVector row;
                              for (j = 0; j < upperColNum; j++) {
                                  row.insert(j, -product4[i][j]);
                              }
                              row.insert(masterColNum + feasibleLeafNodeNum + 1 + counterTemp, -bigMForDomainRest);
                              masterMat.appendRow(row);
                              counterTemp++;
                          }
                          if (lowerContColFiniteUbId[ind]) {
                              double bigMForDomainRest = product5[i] + product10[i] + minValForDomainRest[i] +
                                  contRestBasisInverseRowLcm[i]*contRestColUb[ind];
                              if (bigMForDomainRest > etol) {
//                                  bigMForDomainRest = -10*bigMForDomainRest;
                              } else if (bigMForDomainRest < etol) {
//                                  bigMForDomainRest = 10*bigMForDomainRest;
                              } else {
                                  // bigM = 0; 1 is a random positive number to not let v_i = 1 to occur!
                                  bigMForDomainRest = 1;
                              }
                              CoinPackedVector row;
                              for (j = 0; j < upperColNum; j++) {
                                  row.insert(j, product4[i][j]);
                              }
                              row.insert(masterColNum + feasibleLeafNodeNum + 1 + counterTemp, -bigMForDomainRest);
                              masterMat.appendRow(row);
                              counterTemp++;
                          }
                      } else {
                          double bigMForDomainRest = -product5[i] - product10[i] - maxValForDomainRest[i];
                          if (bigMForDomainRest > etol) {
//                              bigMForDomainRest = -10*bigMForDomainRest;
                          } else if (bigMForDomainRest < etol) {
//                              bigMForDomainRest = 10*bigMForDomainRest;
                          } else {
                              // bigM = 0; 1 is a random positive number to not let v_i = 1 to occur!
                              bigMForDomainRest = 1;
                          }
                          CoinPackedVector row;
                          for (j = 0; j < upperColNum; j++) {
                              row.insert(j, -product4[i][j]);
                          }
                          row.insert(masterColNum + feasibleLeafNodeNum + 1 + counterTemp, -bigMForDomainRest);
                          masterMat.appendRow(row);
                          counterTemp++;
                      }
                  }
                  assert(counterTemp == numBinColsForDomainRest);
                  counterTemp = 0;
                  for (i = 0; i < lowerRowNum; i++) {
                      int ind = contRestBasisIndices[i];
                      if (ind < lowerContColNum) {
                          if (lowerContColFiniteLbId[ind]) {
                              double bigMForDomainRest = (-product5[i] - product10[i] - minValForDomainRest[i] -
                                      contRestBasisInverseRowLcm[i]*contRestColLb[ind]);
                              // Offset bigM by epsilon (= tol[counterTemp])
                              bigMForDomainRest += tol[counterTemp];
                              /*
                              if (bigMForDomainRest <= etol) {
                                  // bigM = 0; 10 is a random number!
                                  bigMForDomainRest = 10;
                              }
                              */
                              CoinPackedVector row;
                              for (j = 0; j < upperColNum; j++) {
                                  row.insert(j, -product4[i][j]);
                              }
                              row.insert(masterColNum + feasibleLeafNodeNum + 1 + counterTemp, bigMForDomainRest);
                              masterMat.appendRow(row);
                              counterTemp++;
                          }
                          if (lowerContColFiniteUbId[ind]) {
                              double bigMForDomainRest = (product5[i] + product10[i] + maxValForDomainRest[i] +
                                      contRestBasisInverseRowLcm[i]*contRestColUb[ind]);
                              // Offset bigM by epsilon (= tol[counterTemp])
                              bigMForDomainRest += tol[counterTemp];
                              /*
                              if (bigMForDomainRest <= etol) {
                                  // bigM = 0; 10 is a random number!
                                  bigMForDomainRest = 10;
                              }
                              */
                              CoinPackedVector row;
                              for (j = 0; j < upperColNum; j++) {
                                  row.insert(j, product4[i][j]);
                              }
                              row.insert(masterColNum + feasibleLeafNodeNum + 1 + counterTemp, bigMForDomainRest);
                              masterMat.appendRow(row);
                              counterTemp++;
                          }
                      } else {
                          double bigMForDomainRest = (-product5[i] - product10[i] - minValForDomainRest[i]);
                          // Offset bigM by epsilon (= tol[counterTemp])
                          bigMForDomainRest += tol[counterTemp];
                          /*
                          if (bigMForDomainRest <= etol) {
                              // bigM = 0; 10 is a random number!
                              bigMForDomainRest = 10;
                          }
                          */
                          CoinPackedVector row;
                          for (j = 0; j < upperColNum; j++) {
                              row.insert(j, -product4[i][j]);
                          }
                          row.insert(masterColNum + feasibleLeafNodeNum + 1 + counterTemp, bigMForDomainRest);
                          masterMat.appendRow(row);
                          counterTemp++;
                      }
                  }
                  assert(counterTemp == numBinColsForDomainRest);
                  CoinPackedVector anotherRow;
                  anotherRow.insert((masterColNum + feasibleLeafNodeNum), numBinColsForDomainRest);
                  for (i = 0; i < numBinColsForDomainRest; i++) {
                      anotherRow.insert((masterColNum + feasibleLeafNodeNum + 1 + i), -1);
                  }
                  masterMat.appendRow(anotherRow);
                  CoinPackedVector oneLastRow;
                  oneLastRow.insert((masterColNum + feasibleLeafNodeNum), 1);
                  for (i = 0; i < numBinColsForDomainRest; i++) {
                      oneLastRow.insert((masterColNum + feasibleLeafNodeNum + 1 + i), -1);
                  }
                  masterMat.appendRow(oneLastRow);
              } else {
                  for (i = 0; i < feasibleLeafNodeNum; i++) {
                      CoinPackedVector row;
                      for (j = 0; j < upperColNum; j++) {
                          row.insert(j, product1[feasibleLeafNodeInd[i]][j]);
                      }
                      row.insert(upperColNum, 1);
                      row.insert((masterColNum + i), -singleBigMForLbf);
                      masterMat.appendRow(row);
                  }
                  CoinPackedVector oneMoreRow;
                  for (i = 0; i < feasibleLeafNodeNum; i++) {
                      oneMoreRow.insert((masterColNum + i), 1);
                  }
                  masterMat.appendRow(oneMoreRow);
              }

              //Number of rows and columns
              if (!level2Infeasible) {
                  masterColNum += feasibleLeafNodeNum + 1 + numBinColsForDomainRest;
                  masterRowNum += feasibleLeafNodeNum + 1 + 2*numBinColsForDomainRest + 2;
              } else {
                  masterColNum += feasibleLeafNodeNum;
                  masterRowNum += feasibleLeafNodeNum + 1;
              }

              /*
              std::cout << std::endl;
              std::cout << "Iter-" << iterCounter << std::endl;
              std::cout << "masterInfeas = " << masterInfeasible <<
                  ", level2Infeas = " << level2Infeasible <<
                  ", subproblemInfeas = " << subproblemInfeasible << std::endl;
              std::cout << "RF Exact = " << rho << ", RF Approx = " <<
                  rhoApproxValue  << ", Master ObjVal = " << masterObjVal << std::endl;
              std::cout << std::endl;
              */
              double UB = (iterCounter ? ((masterObjVal - rhoApproxValue) + rho) : infinity);
              double LB = (iterCounter ? masterObjVal : -infinity);
              std::cout << std::right
                        << std::setw(13) << iterCounter
                        << std::setw(13) << UB
                        << std::setw(13) << LB
                        << std::endl;

              //Storing current iteration's data
              masterObjValPrevIter = masterObjVal;
              memcpy(masterBestSolutionUpperColsPrevIter, masterBestSolutionUpperCols, sizeof(double)*upperColNum);


              /** Setting and solving the master problem **/
              masterBestSolution = new double[masterColNum];
              solver = getSolver(masterProblemSolver, masterMaxThreads, false);
              masterInfeasible = solve(solver,
                      masterColNum, masterObjCoef, upperObjSense,
                      masterColLb, masterColUb, masterColType,
                      &masterMat, masterRowLb, masterRowUb,
                      &masterObjVal, masterBestSolution);
              delete solver;

              //Rounding integer variables and reevaluating master obj. value
              masterObjVal = 0;
              for (i = 0; i < masterColNum; i++) {
                  if (masterColType[i] == 'I' || masterColType[i] == 'B') {
                      masterBestSolution[i] = round(masterBestSolution[i]);
                  }
                  masterObjVal += masterBestSolution[i]*masterObjCoef[i];
              }
              if (iterCounter >= 1) {
                  /*
                  std::cout << "********" << std::endl;
                  std::cout << masterObjValPrevIter << "  " << masterObjVal << std::endl;
                  std::cout << "********" << std::endl;
                  */
                  assert(masterObjVal >= masterObjValPrevIter - etol);
              }
              

              //This checks the infinite loop condition which should never happen
              if (fabs(masterBestSolution[upperColNum] - rhoApproxValue) <= etol) {
                  for (i = 0; i < upperColNum; i++) {
                      if (fabs(masterBestSolutionUpperCols[i] - masterBestSolution[i]) > etol) {
                          break;
                      }
                  }
                  if (i == upperColNum) {
                      //NOTE: instead of modifying master problem, we are simply exiting for the time being
                      tempCounter++;
                      if (tempCounter == 1) {
                          termFlag = true;
                      }
//                      termFlag = true;
                      std::cout << "Algoirthm exited because we are about to enter an infinite loop!\n" << std::endl;
                  }
              }

              memcpy(masterBestSolutionUpperCols, masterBestSolution, sizeof(double)*upperColNum);


              /** Getting solution to master problem **/
              if (!masterInfeasible) {
                  //FIXME: unbounded case? any other invalid case?
                  rhoApproxValue = masterBestSolution[upperColNum];
                  //Evaluating first level part of the obj. function
                  optObjVal = 0;
                  for (i = 0; i < upperColNum; i++) {
                      optObjVal += masterObjCoef[i]*masterBestSolutionUpperCols[i];
                  }
              }

              for (i = 0; i < leafNodeNum; i++) {
                  delete [] product1[i];
              }
              delete [] leafUbCnt;
              delete [] leafLbCnt;
              delete [] leafFeasibilityStatusInd;
              delete [] newColStarts;
              delete [] product10;
              delete [] product9;
              delete [] product5;
              delete [] product2;
              delete [] product7;
              delete [] ubNegDjProduct;
              delete [] lbPosDjProduct;
              delete [] singleDual;
              delete [] product1;
              delete [] feasibleLeafNodeInd;
              delete [] fullDualOfExtraRow;
              if (!level2Infeasible) {
                  for (i = 0; i < lowerRowNum; i++) {
                      delete [] contRestBasisInverseRow[i];
                      delete [] product4[i];
                  }
              }
              delete [] product4;
              delete [] masterBestSolution;
          } else {
              /*
              std::cout << std::endl;
              std::cout << "Iter-" << iterCounter << std::endl;
              std::cout << "masterInfeas = " << masterInfeasible <<
                  ", level2Infeas = " << level2Infeasible <<
                  ", subproblemInfeas = " << subproblemInfeasible << std::endl;
              std::cout << "RF Exact = " << rho << ", RF Approx = " <<
                  rhoApproxValue  << ", Master ObjVal = " << masterObjVal << std::endl;
              */
              double UB = (iterCounter ? ((masterObjVal - rhoApproxValue) + rho) : infinity);
              double LB = (iterCounter ? masterObjVal : -infinity);
              std::cout << std::right
                        << std::setw(13) << iterCounter
                        << std::setw(13) << UB
                        << std::setw(13) << LB
                        << std::endl;
              std::cout << std::endl;
              std::cout << "Termination criterion achieved." << std::endl;
              if (!masterInfeasible && !timeUp) {
                  std::cout << "Optimal Objective Value = " << optObjVal << std::endl;
                  std::cout << "Optimal Solution:" << std::endl;
                  for (i = 0; i < upperColNum; i++) {
                      if (fabs(masterBestSolutionUpperCols[i]) > etol) {
                          std::cout << "    x[" << i << "] = " << masterBestSolutionUpperCols[i] << std::endl;
                      }
                  }
                  for (i = 0; i < lowerColNum; i++) {
                      if (fabs(subproblemBestSolution[i]) > etol) {
                          std::cout << "    y[" << i << "] = " << subproblemBestSolution[i] << std::endl;
                      }
                  }
              } else if (masterInfeasible) {
                  std::cout << "Master problem found infeasible!" << std::endl;
              } else if (timeUp) {
                  std::cout << "\nTime limit exceeded!" << std::endl;
                  std::cout << "Current Best Objective Value = " << currentBestObjVal << std::endl;
                  std::cout << " ... found in iteration-" << currentBestIteration << std::endl;
                  std::cout << "Current Best Solution:" << std::endl;
                  for (i = 0; i < upperColNum; i++) {
                      if (fabs(currentBestSolution[i]) > etol) {
                          std::cout << "    x[" << i << "] = " << currentBestSolution[i] << std::endl;
                      }
                  }
                  for (i = 0; i < lowerColNum; i++) {
                      if (fabs(currentBestSolution[upperColNum + i]) > etol) {
                          std::cout << "    y[" << i << "] = " << currentBestSolution[upperColNum + i] << std::endl;
                      }
                  }
              } else {
                  std::cout << "Error: unknown exit criterion triggered!" << std::endl;
              }
              std::cout << std::endl;

              /*
              for (i = 0; i < leafNodeNum; i++) {
                  if (leafUbCnt[i]) {
                      delete [] leafUbVal[i];
                      delete [] leafUbInd[i];
                  }
                  if (leafLbCnt[i]) {
                      delete [] leafLbVal[i];
                      delete [] leafLbInd[i];
                  }
              }
              delete [] leafUbVal;
              delete [] leafLbVal;
              delete [] leafUbInd;
              delete [] leafLbInd;
              delete [] leafUbCnt;
              delete [] leafLbCnt;
              delete [] leafFeasibilityStatusInd;
              if (!level2Infeasible) {
                  for (i = 0; i < lowerRowNum; i++) {
                      delete [] contRestBasisInverseRow[i];
                  }
              }
              */
          }
          iterCounter++;
      }

//      delete [] lowerObjCoef;
//      delete [] lowerColInd;
//      delete [] feasibleLeafNodeInd;
      delete [] currentBestSolution;
      delete [] minValForDomainRest;
      delete [] maxValForDomainRest;
//      delete [] linkingColId;
      delete [] level2IntColRowActivity;
      delete [] contRestBasisInverseRowLcm;
      delete [] contRestBasisIndices;
      delete [] contRestBasisInverseRow;
      delete [] contRestDjSolution;
      delete [] contRestDualSolution;
      delete [] contRestBestSolutionNonBasics;
      delete [] contRestBestSolution;
      delete [] contRestRowUb;
      delete [] contRestRowLb;
      delete [] contRestColType;
      delete [] contRestColUb;
      delete [] contRestColLb;
      delete [] contRestObjCoef;
      delete [] intBestSolution;
      delete [] level2RowUb;
      delete [] level2RowLb;
      delete [] level2BestSolution;
      delete [] masterBestSolutionUpperColsPrevIter;
      delete [] masterBestSolutionUpperCols;
      delete leafNegDjByRow;
      delete leafPosDjByRow;
      delete leafDualByRow;
      delete [] subproblemBestSolution;
      delete [] tempSubproblemColUb;
      delete [] tempSubproblemColLb;
      delete [] lowerRowRhs;
      delete [] subproblemOrigRhs;
      delete [] subproblemRowSense;
      delete [] subproblemUpperColRowActivity;
      delete [] subproblemRhs;
      delete [] subproblemRowUb;
      delete [] subproblemRowLb;
      delete [] subproblemObjCoef;
      delete [] subproblemColType;
      delete [] subproblemColUb;
      delete [] subproblemColLb;
      delete [] subproblemLowerColInd;
      delete [] lowerContColUb;
      delete [] lowerContColLb;
      delete [] lowerContObjCoef;
      delete [] lowerContColFiniteUbId;
      delete [] lowerContColFiniteLbId;
      delete [] lowerContColInd;
      delete [] lowerIntColInd;
    }
    catch(CoinError& er) {
	std::cerr << "ERROR:" << er.message() << std::endl
		  << " from function " << er.methodName() << std::endl
		  << " from class " << er.className() << std::endl;
    }
    catch(...) {
	std::cerr << "Something went wrong!" << std::endl;
    }

    clock_t end = clock();
    double elapsed_time = (double)(end - begin) / CLOCKS_PER_SEC;

    std::cout << "Total elapsed time = " << elapsed_time << std::endl;

    return 0;
}

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
        CPXsetintparam(cpxEnv, CPX_PARAM_SCRIND, CPX_OFF);
        CPXsetintparam(cpxEnv, CPX_PARAM_THREADS, maxThreads);
        //set various tolerances manually to have better control over them
        CPXsetdblparam(cpxEnv, CPX_PARAM_EPINT, 1e-05);
        CPXsetdblparam(cpxEnv, CPX_PARAM_EPOPT, 1e-06);
        CPXsetdblparam(cpxEnv, CPX_PARAM_EPRHS, 1e-06);
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

    bool ind = true;
    if (ind) {
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
        double objSense, objVal, bigM = 1e+8, etol = 1e-7;
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
            bigM = objVal;
        }
        //Update bigM further based on known bound
        /*
        //bound (=subproblemObjVal) can't be used this way (I think)
        if (bigM >= bound + etol) {
            bigM = bound;
        }
        */
        //Update bigM further based on max of subproblem's objective function
        if (bigM >= boundOnBigM + etol) {
            bigM = boundOnBigM;
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
                memcpy(unbddRay, dualSolver->getPrimalRays(1)[0], sizeof(double)*dualNumCols);
                delete dualSolver;

                //Find appropriate multiplier for unbounded ray
                //FIXME: Use isInfeasible as well as improve following etol later!
                double lambda = 0, objCoefRayProd = 0;
                int nzDual = 0, nzPosDj = 0, nzNegDj = 0;
                int maxColIndDualNz = -1, maxColIndPosDjNz = -1, maxColIndNegDjNz = -1;
                for (j = 0; j < dualNumCols; j++) {
                    objCoefRayProd += dualObjCoef[j]*unbddRay[j];
                }
                lambda = (bigM - dualObjVal)/(objCoefRayProd);
                if (lambda < -etol) {
                    //lambda is negative ==> any positive lambda is a satisfactory value
                    lambda = 1;
                }

                //Evaluate full currect dual
                for (j = 0; j < dualNumCols; j++) {
                    correctFullDual[j] = dualBestSolution[j] + lambda * unbddRay[j];
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
                    dual->modifyCoefficient(i, j, correctFullDual[j], false);
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
                        posDj->modifyCoefficient(i, j, djVal, false);
                        if ((i < negDjMajorDim) && (j < negDjMinorDim)) {
                            negDj->modifyCoefficient(i, j, 0, false);
                        }
                    } else if (djVal < -etol) {
                        negDj->modifyCoefficient(i, j, djVal, false);
                        if ((i < posDjMajorDim) && (j < posDjMinorDim)) {
                            posDj->modifyCoefficient(i, j, 0, false);
                        }
                    }
                }

                delete dualMatByCol;
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
      //Making upper bounds finite wherever required
      //FIXME: Find a workaround for following 1500.
      int i;
      for (i = 0; i < (upperColNum + lowerColNum); i++) {
          if (origColUb[i] >= infinity) {
              origColUb[i] = 1500;
          }
      }

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
      //For detecting infeasible leaf nodes
      int feasibleLeafNodeNum, *feasibleLeafNodeInd;


      /** Initial setup for MILP master problem **/
      //Various arrays and matrices for the problem setup
      int masterColNum = upperColNum;
      int masterRowNum = (upperRowsHaveLowerCols ? 0 : upperRowNum);
      int j;
      bool masterInfeasible = false;
      double bilevelVFApproxValue = -infinity;
      //FIXME: set # of threads appropriately later (as a part of argc/argv?)!
      int masterMaxThreads = 1;
      //Parameter for identifying which MILP solver needs to be used
      //FIXME: Have this as a user input parameter!
      std::string masterProblemSolver = "CPLEX";
      double *masterBestSolution;
      double masterObjVal, optObjVal = 0;
      double *masterBestSolutionUpperCols = new double[upperColNum];

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
      if (upperRowsHaveLowerCols) {
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
      //Save upperColNum part of best solution
      memcpy(masterBestSolutionUpperCols, masterBestSolution, sizeof(double)*upperColNum);
      delete [] masterBestSolution;
      //Evaluating first level part of original MIBLP objective value
      for (i = 0; i < upperColNum; i++) {
          optObjVal += masterObjCoef[i]*masterBestSolutionUpperCols[i];
      }

      //Adding a column to master problem representing bilevel VF approx. value
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
      double *level2IntBestSolution = new double[lowerIntColNum];
      CoinZeroN(level2IntBestSolution, lowerIntColNum);
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
      int *contRestBasisIndices = new int[lowerRowNum];
      int numBinColsForDomainRest = 0;

      //integer restriction matrix of second level rows for second level cols
      CoinPackedMatrix intRestMat(rowCoefMatrixByCol);
      intRestMat.deleteRows(upperRowNum, upperRowInd);
      intRestMat.deleteCols(upperColNum, upperColInd);
      intRestMat.deleteCols(lowerContColNum, lowerContColInd);
      double *level2IntColRowActivity = new double[lowerRowNum];
      CoinZeroN(level2IntColRowActivity, lowerRowNum);
      //Parameter for identifying which LP solver needs to be used
      //FIXME: Have this as a user input parameter!
      std::string contRestProblemSolver = "CPLEX";
      //FIXME: set # of threads appropriately later (as a part of argc/argv?)!
      int contRestMaxThreads = 1;


      /** Initial setup for finding tolerance to impose strict inequality while
          restricting domain of 'x' in master problem **/
      //NOTE: 2*lowerRowNum is max size for tol
      double *tol = new double[2*lowerRowNum];
      //Data required to identify linking cols
      int *linkingColId = new int[upperColNum + lowerColNum];
      memcpy(linkingColId, origMibsModel.getFixedInd(),
              sizeof(int)*(upperColNum + lowerColNum));
      int linkingColNum = origMibsModel.getSizeFixedInd();

      //Misc declarations
      bool termFlag = false, timeUp = false;
      int iterCounter = 0;
      double bilevelVFExactValue = infinity, boundOnLbf, dualBoundOnLevel2 = 0;
      double bigMForUbf = 1e+7, *bigMForLbf;
      double dualBoundOnSubproblem = 0, primalBoundOnSubproblem = 0;
      double *maxValForDomainRest = new double[lowerRowNum];
      double *minValForDomainRest = new double[lowerRowNum];
      for (i = 0; i < lowerColNum; i++) {
          double lb = subproblemColLb[i];
          assert(lb > -infinity);
          double ub = subproblemColUb[i];
          assert(ub < infinity);
          double coef1 = subproblemObjCoef[i];
          double coef2 = lowerObjCoef[i];
          if (coef1 < -etol) {
              dualBoundOnSubproblem += coef1*lb;
              primalBoundOnSubproblem += coef1*ub;
          } else if (coef1 > etol) {
              dualBoundOnSubproblem += coef1*ub;
              primalBoundOnSubproblem += coef1*lb;
          }
          if (coef2 < -etol) {
              dualBoundOnLevel2 += coef2*lb;
          } else if (coef2 > etol) {
              dualBoundOnLevel2 += coef2*ub;
          }
      }
      //Value of boundOnLbf: 10 is a random multiplier
      if (fabs(primalBoundOnSubproblem) >= fabs(dualBoundOnSubproblem)) {
          boundOnLbf = 10*fabs(primalBoundOnSubproblem);
      } else {
          boundOnLbf = 10*fabs(dualBoundOnSubproblem);
      }
      //Also assert that upper level col bounds are finite!
      for (i = 0; i < upperColNum; i++) {
          assert(masterColLb[i] > -infinity);
          assert(masterColUb[i] < infinity);
      }



      /*** while loop for decomposition algorithm ***/
      while (!termFlag) {
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
              }
              delete solver;
              for (i = 0; i < lowerIntColNum; i++) {
                  level2IntBestSolution[i] = level2BestSolution[lowerIntColInd[i]];
              }

              if (!level2Infeasible) {
                  //If second level feasible, solve the continuous restriction
                  //    for the known level2IntBestSolution

                  //Finding product of integer best solution and corresponding obj. vector
                  level2IntObjVal = 0;
                  for (i = 0; i < lowerIntColNum; i++) {
                      level2IntObjVal += level2IntBestSolution[i]*lowerObjCoef[lowerIntColInd[i]];
                  }

                  //Finding product of integer best solution and integer matrix
                  intRestMat.times(level2IntBestSolution, level2IntColRowActivity);

                  //Finding RowLb and RowUb for continuous restriction
                  memcpy(contRestRowLb, level2RowLb, sizeof(double)*lowerRowNum);
                  memcpy(contRestRowUb, level2RowUb, sizeof(double)*lowerRowNum);
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
                  solver->getBasics(contRestBasisIndices);
                  delete solver;

                  if (contRestInfeasible) {
                      //Should not happen because level2Infeasible is false!
                      std::cout << 
                          "Error: Restriction is infeasible whereas it should not be!."
                          << std::endl;
                      return 0;
                  }

                  /* Updating subproblem data with new info. from continuous restriction */
                  if (addRowInd) {
                      subproblemMat.appendRow(subproblemColNum, subproblemLowerColInd, lowerObjCoef);
                      addRowInd = false;
                  }
                  if (lowerObjSense == 1) {
                      subproblemRowLb[subproblemRowNum - 1] = -infinity;
                      subproblemRowUb[subproblemRowNum - 1] = contRestObjVal + level2IntObjVal;
                  } else {
                      subproblemRowLb[subproblemRowNum - 1] = contRestObjVal + level2IntObjVal;
                      subproblemRowUb[subproblemRowNum - 1] = infinity;
                  }
                  subproblemRhs[subproblemRowNum - 1] = contRestObjVal + level2IntObjVal;
              } else {
                  contRestObjVal = infinity;
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

              // Updating value function exact value
              if (!subproblemInfeasible) {
                  bilevelVFExactValue = subproblemObjVal;
              }


              /** Getting dual information to the subproblem **/
              getDualData(solver, dualBoundOnSubproblem,
                      &leafNodeNum, &leafFeasibilityStatus,
                      leafDualByRow, leafPosDjByRow, leafNegDjByRow,
                      &leafLbCnt, &leafLbInd, &leafLbVal,
                      &leafUbCnt, &leafUbInd, &leafUbVal);

              delete solver;

              /* Check bilevel feasibility */
              if (!subproblemInfeasible) {
                  double level2ObjValForSubproblemSol = 0;
                  for (i = 0; i < lowerColNum; i++) {
                      level2ObjValForSubproblemSol += subproblemBestSolution[i]*lowerObjCoef[i];
                  }
                  assert(fabs(level2ObjValForSubproblemSol - level2ObjVal) <= etol);
                  //Updating original MIBLP objective value with its second level part
                  for (i = 0; i < lowerColNum; i++) {
                      optObjVal += subproblemObjCoef[i]*subproblemBestSolution[i];
                  }
              }
          }


          /** Checking termination criterion **/
          //FIXME: are the criteria correct?
          clock_t current = clock();
          double timeTillNow = (double) (current - begin) / CLOCKS_PER_SEC;
          timeUp = ((timeTillNow >= 3600) ? true : false);
          if (masterInfeasible || (!subproblemInfeasible &&
                  (fabs(bilevelVFExactValue - bilevelVFApproxValue) <= etol))
                  || timeUp) {
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

              //FIXME: remove feasibleLeafNodeNum and feasibleLeafNodeInd later!
              feasibleLeafNodeNum = leafNodeNum;
              feasibleLeafNodeInd = new int[feasibleLeafNodeNum];
              CoinIotaN(feasibleLeafNodeInd, feasibleLeafNodeNum, 0);

              //bigM for LBF of subproblem
              bigMForLbf = new double[feasibleLeafNodeNum];

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
              //FIXME: Is following initialization necessary?
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
              if (!level2Infeasible) {
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
                          product5[i] += contRestBasisInverseRow[i][j]*
                              (level2IntColRowActivity[j] - lowerRowRhs[j]);
                      }
                  }

                  //Product of dual of cont. rest. and lowerRowRhs
                  for (i = 0; i < lowerRowNum; i++) {
                      product6 += contRestDualSolution[i]*lowerRowRhs[i];
                  }

                  //bigM for UBF of level2 problem
                  double maxVal = 0;
                  for (i = 0; i < upperColNum; i++) {
                      double coef = product2[i];
                      if (coef < -etol) {
                          maxVal += coef*masterColLb[i];
                      } else if (coef > etol) {
                          maxVal += coef*masterColUb[i];
                      }
                  }
                  double ubfMin = product6 - product3 + level2IntObjVal - maxVal;
                  bigMForUbf = fabs(dualBoundOnLevel2) + fabs(ubfMin);

                  //bigM for LBF of subproblem
                  for (i = 0; i < feasibleLeafNodeNum; i++) {
                      bigMForLbf[feasibleLeafNodeInd[i]] = boundOnLbf;
                      double minVal = 0, maxVal = 0;
                      for (j = 0; j < upperColNum; j++) {
                          double coef = (product1[feasibleLeafNodeInd[i]][j] +
                                      fullDualOfExtraRow[feasibleLeafNodeInd[i]]*product2[j]);
                          if (coef < -etol) {
                              minVal += coef*masterColUb[j];
                          } else if (coef > etol) {
                              minVal += coef*masterColLb[j];
                          }
                      }
                      if (fullDualOfExtraRow[feasibleLeafNodeInd[i]] < -etol) {
                          maxVal = 0;
                      } else if (fullDualOfExtraRow[feasibleLeafNodeInd[i]] > etol) {
                          maxVal = bigMForUbf*fullDualOfExtraRow[feasibleLeafNodeInd[i]];
                      }
                      bigMForLbf[feasibleLeafNodeInd[i]] += fabs(product7[feasibleLeafNodeInd[i]] +
                               lbPosDjProduct[feasibleLeafNodeInd[i]] +
                               ubNegDjProduct[feasibleLeafNodeInd[i]] +
                               fullDualOfExtraRow[feasibleLeafNodeInd[i]]*(product6 - product3 + level2IntObjVal)
                               - minVal + maxVal);
                  }
              } else {
                  //bigM for LBF of subproblem
                  for (i = 0; i < feasibleLeafNodeNum; i++) {
                      bigMForLbf[feasibleLeafNodeInd[i]] = boundOnLbf;
                      double minVal = 0;
                      for (j = 0; j < upperColNum; j++) {
                          double coef = product1[feasibleLeafNodeInd[i]][j];
                          if (coef < -etol) {
                              minVal += coef*masterColUb[j];
                          } else if (coef > etol) {
                              minVal += coef*masterColLb[j];
                          }
                      }
                      bigMForLbf[feasibleLeafNodeInd[i]] += fabs(product7[feasibleLeafNodeInd[i]] +
                               lbPosDjProduct[feasibleLeafNodeInd[i]] +
                               ubNegDjProduct[feasibleLeafNodeInd[i]] +
                               - minVal);
                  }
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
                              for (j = 0; j < upperColNum; j++) {
                                  double currCoef = fabs(product4[i][j]);
                                  if (currCoef > etol) {
                                      tol[numBinColsForDomainRest] = currCoef;
                                      break;
                                  }
                              }
                              int k;
                              for (k = j+1; k < upperColNum; k++) {
                                  double currCoef = fabs(product4[i][k]);
                                  if ((currCoef > etol) && (currCoef < tol[numBinColsForDomainRest])) {
                                      tol[numBinColsForDomainRest] = currCoef;
                                  }
                              }

                              //Counting # of binary variables for domain rest.
                              numBinColsForDomainRest++;
                          }
                          if (lowerContColFiniteUbId[ind]) {
                              for (j = 0; j < upperColNum; j++) {
                                  //Finding "tol"
                                  double currCoef = fabs(product4[i][j]);
                                  if (currCoef > etol) {
                                      tol[numBinColsForDomainRest] = currCoef;
                                      break;
                                  }
                              }
                              int k;
                              for (k = j+1; k < upperColNum; k++) {
                                  double currCoef = fabs(product4[i][k]);
                                  if ((currCoef > etol) && (currCoef < tol[numBinColsForDomainRest])) {
                                      tol[numBinColsForDomainRest] = currCoef;
                                  }
                              }

                              //Counting # of binary variables for domain rest.
                              numBinColsForDomainRest++;
                          }
                      } else {
                          for (j = 0; j < upperColNum; j++) {
                              //Finding "tol"
                              double currCoef = fabs(product4[i][j]);
                              if (currCoef > etol) {
                                  tol[numBinColsForDomainRest] = currCoef;
                                  break;
                              }
                          }
                          int k;
                          for (k = j+1; k < upperColNum; k++) {
                              double currCoef = fabs(product4[i][k]);
                              if ((currCoef > etol) && (currCoef < tol[numBinColsForDomainRest])) {
                                  tol[numBinColsForDomainRest] = currCoef;
                              }
                          }

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
                               fullDualOfExtraRow[feasibleLeafNodeInd[i]]*(product6 - product3 + level2IntObjVal) -
                               bigMForLbf[feasibleLeafNodeInd[i]]));
                  }
                  masterRowLbVec.resize((masterRowNum + feasibleLeafNodeNum + 1), 1);
                  int counterTemp = 0;
                  for (i = 0; i < lowerRowNum; i++) {
                      int ind = contRestBasisIndices[i];
                      if (ind < lowerContColNum) {
                          if (lowerContColFiniteLbId[ind]) {
                              masterRowLbVec.resize((masterRowNum + feasibleLeafNodeNum + 1 + (counterTemp+1)),
                                      (product5[i] + contRestColLb[ind]));
                              counterTemp++;
                          }
                          if (lowerContColFiniteUbId[ind]) {
                              masterRowLbVec.resize((masterRowNum + feasibleLeafNodeNum + 1 + (counterTemp+1)),
                                      (-product5[i] - contRestColUb[ind]));
                              counterTemp++;
                          }
                      } else {
                          masterRowLbVec.resize((masterRowNum + feasibleLeafNodeNum + 1 + (counterTemp+1)), product5[i]);
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
                              //Note: 10 is a random multiplier
                              double bigMForDomainRest = 10*fabs(-product5[i] - minValForDomainRest[i] - contRestColLb[ind]);
                              masterRowUbVec.resize((masterRowNum + feasibleLeafNodeNum + 1 + numBinColsForDomainRest +
                                          (counterTemp+1)), (bigMForDomainRest + product5[i] + contRestColLb[ind] - tol[counterTemp]));
                              counterTemp++;
                          }
                          if (lowerContColFiniteUbId[ind]) {
                              double bigMForDomainRest = 10*fabs(product5[i] + maxValForDomainRest[i] + contRestColUb[ind]);
                              masterRowUbVec.resize((masterRowNum + feasibleLeafNodeNum + 1 + numBinColsForDomainRest +
                                          (counterTemp+1)), (bigMForDomainRest - product5[i] - contRestColUb[ind] - tol[counterTemp]));
                              counterTemp++;
                          }
                      } else {
                          double bigMForDomainRest = 10*fabs(-product5[i] - minValForDomainRest[i]);
                          masterRowUbVec.resize((masterRowNum + feasibleLeafNodeNum + 1 + numBinColsForDomainRest +
                                      (counterTemp+1)), (bigMForDomainRest + product5[i] - tol[counterTemp]));
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
                               bigMForLbf[feasibleLeafNodeInd[i]]));
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
                      row.insert((masterColNum + i), -bigMForLbf[feasibleLeafNodeInd[i]]);
                      row.insert((masterColNum + feasibleLeafNodeNum), -bigMForUbf*fullDualOfExtraRow[feasibleLeafNodeInd[i]]);
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
                              double bigMForDomainRest = -product5[i] - maxValForDomainRest[i] - contRestColLb[ind];
                              //Note: 10 is a random multiplier
                              if (bigMForDomainRest > etol) {
                                  bigMForDomainRest = -10*bigMForDomainRest;
                              } else {
                                  bigMForDomainRest = 10*bigMForDomainRest;
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
                              double bigMForDomainRest = product5[i] + minValForDomainRest[i] + contRestColUb[ind];
                              if (bigMForDomainRest > etol) {
                                  bigMForDomainRest = -10*bigMForDomainRest;
                              } else {
                                  bigMForDomainRest = 10*bigMForDomainRest;
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
                          double bigMForDomainRest = -product5[i] - maxValForDomainRest[i];
                          if (bigMForDomainRest > etol) {
                              bigMForDomainRest = -10*bigMForDomainRest;
                          } else {
                              bigMForDomainRest = 10*bigMForDomainRest;
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
                              //Note: 10 is a random multiplier
                              double bigMForDomainRest = 10*fabs(-product5[i] - minValForDomainRest[i] - contRestColLb[ind]);
                              CoinPackedVector row;
                              for (j = 0; j < upperColNum; j++) {
                                  row.insert(j, -product4[i][j]);
                              }
                              row.insert(masterColNum + feasibleLeafNodeNum + 1 + counterTemp, bigMForDomainRest);
                              masterMat.appendRow(row);
                              counterTemp++;
                          }
                          if (lowerContColFiniteUbId[ind]) {
                              double bigMForDomainRest = 10*fabs(product5[i] + maxValForDomainRest[i] + contRestColUb[ind]);
                              CoinPackedVector row;
                              for (j = 0; j < upperColNum; j++) {
                                  row.insert(j, product4[i][j]);
                              }
                              row.insert(masterColNum + feasibleLeafNodeNum + 1 + counterTemp, bigMForDomainRest);
                              masterMat.appendRow(row);
                              counterTemp++;
                          }
                      } else {
                          double bigMForDomainRest = 10*fabs(-product5[i] - minValForDomainRest[i]);
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
                      row.insert((masterColNum + i), -bigMForLbf[feasibleLeafNodeInd[i]]);
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

              std::cout << std::endl;
              std::cout << "Iter-" << iterCounter << std::endl;
              std::cout << "masterInfeas = " << masterInfeasible <<
                  ", level2Infeas = " << level2Infeasible <<
                  ", subproblemInfeas = " << subproblemInfeasible << std::endl;
              std::cout << "VF Exact = " << bilevelVFExactValue << ", VF Approx = " <<
                  bilevelVFApproxValue  << ", Master ObjVal = " << masterObjVal << std::endl;
              std::cout << std::endl;


              /** Setting and solving the master problem **/
              masterBestSolution = new double[masterColNum];
              solver = getSolver(masterProblemSolver, masterMaxThreads, false);
              masterInfeasible = solve(solver,
                      masterColNum, masterObjCoef, upperObjSense,
                      masterColLb, masterColUb, masterColType,
                      &masterMat, masterRowLb, masterRowUb,
                      &masterObjVal, masterBestSolution);
              delete solver;

              //This checks the infinite loop condition which should never happen
              if (masterBestSolution[upperColNum] == bilevelVFApproxValue) {
                  for (i = 0; i < upperColNum; i++) {
                      if (fabs(masterBestSolutionUpperCols[i] - masterBestSolution[i]) >= etol) {
                          break;
                      }
                  }
                  if (i == upperColNum) {
                      //NOTE: instead of modifying master problem, we are simply exiting for the time being
                      termFlag = true;
                      std::cout << "Algoirthm exited because we are about to enter an infinite loop!\n" << std::endl;
                  }
              }

              memcpy(masterBestSolutionUpperCols, masterBestSolution, sizeof(double)*upperColNum);


              /** Getting solution to master problem **/
              if (!masterInfeasible) {
                  //FIXME: unbounded case? any other invalid case?
                  bilevelVFApproxValue = masterBestSolution[upperColNum];
                  //Evaluating first level part of original MIBLP objective function
                  optObjVal = 0;
                  for (i = 0; i < upperColNum; i++) {
                      optObjVal += masterObjCoef[i]*masterBestSolutionUpperCols[i];
                  }
              }

              for (i = 0; i < leafNodeNum; i++) {
                  delete [] product1[i];
//                  delete [] leafUbVal[i];
//                  delete [] leafLbVal[i];
//                  delete [] leafUbInd[i];
//                  delete [] leafLbInd[i];
              }
//              delete [] leafUbVal;
//              delete [] leafLbVal;
//              delete [] leafUbInd;
//              delete [] leafLbInd;
//              delete [] leafUbCnt;
//              delete [] leafLbCnt;
//              delete [] leafFeasibilityStatus;
              delete [] newColStarts;
              delete [] product5;
              delete [] product2;
              delete [] product7;
              delete [] ubNegDjProduct;
              delete [] lbPosDjProduct;
              delete [] singleDual;
              delete [] product1;
              delete [] bigMForLbf;
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
              std::cout << std::endl;
              std::cout << "Iter-" << iterCounter << std::endl;
              std::cout << "masterInfeas = " << masterInfeasible <<
                  ", level2Infeas = " << level2Infeasible <<
                  ", subproblemInfeas = " << subproblemInfeasible << std::endl;
              std::cout << "VF Exact = " << bilevelVFExactValue << ", VF Approx = " <<
                  bilevelVFApproxValue  << ", Master ObjVal = " << masterObjVal << std::endl;
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
              } else {
                  std::cout << "Error: unknown exit criterion triggered!" << std::endl;
              }
              std::cout << std::endl;
          }
          iterCounter++;
      }

//      delete [] lowerObjCoef;
//      delete [] lowerColInd;
//      delete [] bigMForLbf;
//      delete [] feasibleLeafNodeInd;
      delete [] minValForDomainRest;
      delete [] maxValForDomainRest;
      delete [] linkingColId;
      delete [] tol;
      delete [] level2IntColRowActivity;
      delete [] contRestBasisIndices;
      delete [] contRestBasisInverseRow;
      delete [] contRestDualSolution;
      delete [] contRestBestSolution;
      delete [] contRestRowUb;
      delete [] contRestRowLb;
      delete [] contRestColType;
      delete [] contRestColUb;
      delete [] contRestColLb;
      delete [] contRestObjCoef;
      delete [] level2RowUb;
      delete [] level2RowLb;
      delete [] level2IntBestSolution;
      delete [] level2BestSolution;
      delete [] masterBestSolutionUpperCols;
      delete leafNegDjByRow;
      delete leafPosDjByRow;
      delete leafDualByRow;
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

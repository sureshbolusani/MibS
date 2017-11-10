/*===========================================================================*/
/* This file is part of a Mixed Integer Bilevel Solver                       */
/* developed using the BiCePS Linear Integer Solver (BLIS).                  */
/*                                                                           */
/* Authors: Suresh Bolusani, Lehigh University                               */
/*          Scott DeNegre, Lehigh University                                 */
/*          Ted Ralphs, Lehigh University                                    */
/*          Sahar Tahernajad, Lehigh University                              */
/*                                                                           */
/* Copyright (C) 2007-2017 Lehigh University, Scott DeNegre, and Ted Ralphs. */
/* All Rights Reserved.                                                      */
/*                                                                           */
/* This software is licensed under the Eclipse Public License. Please see    */
/* accompanying file for terms.                                              */
/*===========================================================================*/

#ifndef MibSWarmStart_h_
#define MibSWarmStart_h_

#include "MibSTreeNode.hpp"
#include "MibSBranchObjectInt.hpp"

//#############################################################################

class MibSWarmStart {

    friend class MibSModel;

private:

    AlpsTreeNode *root_;

    int leafNodeNum_;
    int *leafDepth_;
    //TODO: Should I use BlisLpStatus or AlpsNodeStatus?
    BlisLpStatus *leafFeasibilityStatus_;
    bool *leafDualInfoUsageStatus_;
    double *leafLowerBound_;
//    MibSBranchObjectInt **leafBranchPath_;
    double **leafLb_;
    double **leafUb_;

    /*
    int leafDualsNonzeroNum_;
    int *leafDualsRowIndex_;
    int *leafDualsColIndex_;
    double *leafDualsVal_;
    */
    CoinPackedMatrix *leafDualsByRow_;

    /*
    int leafDjsNonzeroNum_;
    int *leafDjsRowIndex_;
    int *leafDjsColIndex_;
    double *leafDjsVal_;
    */
    CoinPackedMatrix *leafDjsByRow_;

    /*
    int leafPosDjsNonzeroNum_;
    int *leafPosDjsRowIndex_;
    int *leafPosDjsColIndex_;
    double *leafPosDjsVal_;
    */
    CoinPackedMatrix *leafPosDjsByRow_;

    /*
    int leafNegDjsNonzeroNum_;
    int *leafNegDjsRowIndex_;
    int *leafNegDjsColIndex_;
    double *leafNegDjsVal_;
    */
    CoinPackedMatrix *leafNegDjsByRow_;

public:
   
    MibSWarmStart();
   
    ~MibSWarmStart();

    AlpsTreeNode * getRootNode() {return root_;}
    int getLeafNodeNum() {return leafNodeNum_;}
    int * getLeafDepths() {return leafDepth_;}
    BlisLpStatus * getLeafFeasibilityStati() {return leafFeasibilityStatus_;}
    bool * getLeafDualInfoUsageStati() {return leafDualInfoUsageStatus_;}
    double * getLeafLowerBounds() {return leafLowerBound_;}
//    MibSBranchObjectInt ** getLeafBranchPaths() {return leafBranchPath_;}
    double ** getLeafLBs() {return leafLb_;}
    double ** getLeafUBs() {return leafUb_;}
    CoinPackedMatrix * getLeafDualsByRow() {return leafDualsByRow_;}
    CoinPackedMatrix * getLeafDjsByRow() {return leafDjsByRow_;}
    CoinPackedMatrix * getLeafPosDjsByRow() {return leafPosDjsByRow_;}
    CoinPackedMatrix * getLeafNegDjsByRow() {return leafNegDjsByRow_;}

    void setMibsRootNode(AlpsTreeNode *root) {root_ = root;}
    void setLeafNodeNum(int nodeNum) {leafNodeNum_ = nodeNum;}
    void setLeafDepths(int *depth) {leafDepth_ = depth;}
    void setLeafFeasibilityStati(BlisLpStatus *feasibilityStatus) 
        {leafFeasibilityStatus_ = feasibilityStatus;}
    void setLeafDualInfoUsageStati(bool *val) {leafDualInfoUsageStatus_ = val;};
    void setLeafLowerBounds(double *lowerBound) {leafLowerBound_ = lowerBound;}
//    void setLeafBranchPaths(MibSBranchObjectInt **branchPath) 
//        {leafBranchPath_ = branchPath;}
    void setLeafLBs(double **leafLb) {leafLb_ = leafLb;}
    void setLeafUBs(double **leafUb) {leafUb_ = leafUb;}
    void setLeafDualsByRow(CoinPackedMatrix *duals) {leafDualsByRow_ = duals;}
    void setLeafDjsByRow(CoinPackedMatrix *djs) {leafDjsByRow_ = djs;}
    void setLeafPosDjsByRow(CoinPackedMatrix *posDjs) {leafPosDjsByRow_ = posDjs;}
    void setLeafNegDjsByRow(CoinPackedMatrix *negDjs) {leafNegDjsByRow_ = negDjs;}

   
private:
   
};

#endif

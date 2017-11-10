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

#include "MibSWarmStart.hpp"

//#############################################################################
MibSWarmStart::MibSWarmStart()
{
    root_ = NULL;
    leafNodeNum_ = 0;
    leafDepth_ = NULL;
    leafFeasibilityStatus_ = NULL;
    leafDualInfoUsageStatus_ = NULL;
    leafLowerBound_ = NULL;
//    leafBranchPath_ = NULL;
    leafLb_ = NULL;
    leafUb_ = NULL;

    /*
    leafDualsNonzeroNum_ = 0;
    leafDualsRowIndex_ = NULL;
    leafDualsColIndex_ = NULL;
    leafDualsVal_ = NULL;
    */
    leafDualsByRow_ = NULL;

    /*
    leafDjsNonzeroNum_ = 0;
    leafDjsRowIndex_ = NULL;
    leafDjsColIndex_ = NULL;
    leafDjsVal_ = NULL;
    */
    leafDjsByRow_ = NULL;

    /*
    leafPosDjsNonzeroNum_ = 0;
    leafPosDjsRowIndex_ = NULL;
    leafPosDjsColIndex_ = NULL;
    leafPosDjsVal_ = NULL;
    */
    leafPosDjsByRow_ = NULL;

    /*
    leafNegDjsNonzeroNum_ = 0;
    leafNegDjsRowIndex_ = NULL;
    leafNegDjsColIndex_ = NULL;
    leafNegDjsVal_ = NULL;
    */
    leafNegDjsByRow_ = NULL;
}

//#############################################################################
MibSWarmStart::~MibSWarmStart()
{
    //TODO: Will this destructor prevent from having warm start info at the end
    // of MIBLP solve?
    if (root_) delete root_;
    if (leafDepth_) delete [] leafDepth_;
    if (leafFeasibilityStatus_) delete [] leafFeasibilityStatus_;
    if (leafDualInfoUsageStatus_) delete [] leafDualInfoUsageStatus_;
    if (leafLowerBound_) delete [] leafLowerBound_;
//    if (leafBranchPath_) delete [] leafBranchPath_;
    if (leafLb_) delete [] leafLb_;
    if (leafUb_) delete [] leafUb_;

    /*
    if (leafDualsRowIndex_) delete [] leafDualsRowIndex_;
    if (leafDualsColIndex_) delete [] leafDualsColIndex_;
    if (leafDualsVal_) delete [] leafDualsVal_;
    */
    if (leafDualsByRow_) delete leafDualsByRow_;

    /*
    if (leafDjsRowIndex_) delete [] leafDjsRowIndex_;
    if (leafDjsColIndex_) delete [] leafDjsColIndex_;
    if (leafDjsVal_) delete [] leafDjsVal_;
    */
    if (leafDjsByRow_) delete leafDjsByRow_;

    /*
    if (leafPosDjsRowIndex_) delete [] leafPosDjsRowIndex_;
    if (leafPosDjsColIndex_) delete [] leafPosDjsColIndex_;
    if (leafPosDjsVal_) delete [] leafPosDjsVal_;
    */
    if (leafPosDjsByRow_) delete leafPosDjsByRow_;

    /*
    if (leafNegDjsRowIndex_) delete [] leafNegDjsRowIndex_;
    if (leafNegDjsColIndex_) delete [] leafNegDjsColIndex_;
    if (leafNegDjsVal_) delete [] leafNegDjsVal_;
    */
    if (leafNegDjsByRow_) delete leafNegDjsByRow_;
}

//#############################################################################

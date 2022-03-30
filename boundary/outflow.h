#ifndef _OUTFLOW_H_
#define _OUTFLOW_H_    

#include "../core/cell.h"
#include "../core/localoperator.h"

namespace LBM{

    //文献 "2013 Evaluation of outflow boundary conditions for two-phase lattice Boltzmann equation"

    //边界的方位  X，Y表示边界外法线方向沿着X, Y方向  P表示positive正向  N表示negative负向
    //这个枚举类型并不使用，只是为了说明下面各个类的后缀含义
    //enum class BoundaryOrientation { EdgeXP, EdgeXN, EdgeYP, EdgeYN, CornerPP, CornerPN, CornerNP, CornerNN };
    
    //8种情况
    
    /**************************************** 流出边界OutFlowEdgeXP定义 ****************************************/
    //流出边界用于碰撞中?
    struct OutFlowEdgeXP : public LocalOperator {

        OutFlowEdgeXP(){}

        //后面三个形参的作用是可以使用第1种executeLocalOperator,这是将其用于碰撞中
        void operator()(Cell& cell, const double& a, const double& b, const double& c) {
            for (int k = 0; k < 9; ++k) {
                cell.f[k] = cell.neighbors[3] -> f[k];
            }
        }

    };    
    /**************************************** 流出边界OutFlowEdgeXN定义 ****************************************/
    struct OutFlowEdgeXN : public LocalOperator {

        OutFlowEdgeXN(){}

        //后面三个形参的作用是可以使用第1种executeLocalOperator,这是将其用于碰撞中
        void operator()(Cell& cell, const double& a, const double& b, const double& c) {
            for (int k = 0; k < 9; ++k) {
                cell.f[k] = cell.neighbors[1] -> f[k];
            }
        }

    };    
    /**************************************** 流出边界OutFlowEdgeYP定义 ****************************************/
    struct OutFlowEdgeYP : public LocalOperator {

        OutFlowEdgeYP(){}

        //后面三个形参的作用是可以使用第1种executeLocalOperator,这是将其用于碰撞中
        void operator()(Cell& cell, const double& a, const double& b, const double& c) {
            for (int k = 0; k < 9; ++k) {
                cell.f[k] = cell.neighbors[4] -> f[k];
            }
        }

    };    
    /**************************************** 流出边界OutFlowEdgeYN定义 ****************************************/
    struct OutFlowEdgeYN : public LocalOperator {

        OutFlowEdgeYN(){}

        //后面三个形参的作用是可以使用第1种executeLocalOperator,这是将其用于碰撞中
        void operator()(Cell& cell, const double& a, const double& b, const double& c) {
            for (int k = 0; k < 9; ++k) {
                cell.f[k] = cell.neighbors[2] -> f[k];
            }
        }

    };    
    /**************************************** 流出边界OutFlowCornerPP定义 ****************************************/
    struct OutFlowCornerPP : public LocalOperator {

        OutFlowCornerPP(){}

        //后面三个形参的作用是可以使用第1种executeLocalOperator,这是将其用于碰撞中
        void operator()(Cell& cell, const double& a, const double& b, const double& c) {
            for (int k = 0; k < 9; ++k) {
                cell.f[k] = cell.neighbors[7] -> f[k];
            }
        }

    };    
    /**************************************** 流出边界OutFlowCornerPN定义 ****************************************/
    struct OutFlowCornerPN : public LocalOperator {

        OutFlowCornerPN(){}

        //后面三个形参的作用是可以使用第1种executeLocalOperator,这是将其用于碰撞中
        void operator()(Cell& cell, const double& a, const double& b, const double& c) {
            for (int k = 0; k < 9; ++k) {
                cell.f[k] = cell.neighbors[6] -> f[k];
            }
        }

    };   
    /**************************************** 流出边界OutFlowCornerNP定义 ****************************************/
    struct OutFlowCornerNP : public LocalOperator {

        OutFlowCornerNP(){}

        //后面三个形参的作用是可以使用第1种executeLocalOperator,这是将其用于碰撞中
        void operator()(Cell& cell, const double& a, const double& b, const double& c) {
            for (int k = 0; k < 9; ++k) {
                cell.f[k] = cell.neighbors[8] -> f[k];
            }
        }

    };   
    /**************************************** 流出边界OutFlowCornerNN定义 ****************************************/
    struct OutFlowCornerNN : public LocalOperator {

        OutFlowCornerNN(){}

        //后面三个形参的作用是可以使用第1种executeLocalOperator,这是将其用于碰撞中
        void operator()(Cell& cell, const double& a, const double& b, const double& c) {
            for (int k = 0; k < 9; ++k) {
                cell.f[k] = cell.neighbors[5] -> f[k];
            }
        }

    };

}

#endif
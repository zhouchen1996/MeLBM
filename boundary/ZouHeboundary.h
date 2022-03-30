#ifndef _ZOUHEBOUNDARY_H_
#define _ZOUHEBOUNDARY_H_    

#include "../core/cell.h"
#include "../core/localoperator.h"

namespace LBM{

    //经过检验该方法不适合我的计算模型，因为我是把边界处理放在碰撞这一步的，而一般zouhe边界是在streaming之后的
    //所以只写了一种情况，不建议使用

    //文献 "1997 On pressure and velocity boundary conditions for the lattice Boltzmann BGK model"

    //边界的方位  X，Y表示边界外法线方向沿着X, Y方向  P表示positive正向  N表示negative负向
    //这个枚举类型并不使用，只是为了说明下面各个类的后缀含义
    //enum class BoundaryOrientation { EdgeXP, EdgeXN, EdgeYP, EdgeYN, CornerPP, CornerPN, CornerNP, CornerNN };
    
    //16种情况 速度 压力
    
    /**************************************** 速度边界ZouHeVelocityBoundaryEdgeXN定义 ****************************************/
    struct ZouHeVelocityBoundaryEdgeXN : public LocalOperator {
        
        double velocity[2];

        ZouHeVelocityBoundaryEdgeXN(double ux, double uy){
            velocity[0] = ux;
            velocity[1] = uy;
        }

        //后面三个形参的作用是可以使用第1种executeLocalOperator,这是将其用于碰撞中 
        void operator()(Cell& cell, const double& a, const double& b, const double& c) {
            double u = velocity[0];
            double v = velocity[1];
            //计算密度
            double rho = (cell.f[0] + cell.f[2] + cell.f[4] + 2.0 * (cell.f[3] + cell.f[6] + cell.f[7])) / (1.0 - u);
            cell.f[1] = cell.f[3] + 2.0 * rho * u / 3.0;
			cell.f[5] = cell.f[7] - 0.5 * (cell.f[2] - cell.f[4]) + rho * u / 6.0;
			cell.f[8] = cell.f[6] + 0.5 * (cell.f[2] - cell.f[4]) + rho * u / 6.0;
        }

    };    
    

}

#endif
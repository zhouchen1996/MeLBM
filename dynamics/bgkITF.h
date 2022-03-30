#ifndef _BGKITF_H_
#define _BGKITF_H_

#include "../core/cell.h"
#include "../core/localoperator.h"
#include <cmath>

namespace LBM {

    /**********************************BGKITF*************************************/
    //2006 An adaptive scheme using hierarchical grids for lattice Boltzmann multi-phase flow simulations
    //根据文献中的3DMRT的推导, 推导2DBGK的形式
    //参考 Three-dimensional lattice Boltzmann model for immiscible two-phase flow simulations 式(23 - 24)
    // Dirk Kehrwald 的博士论文 Numerical Analysis of Immiscible Lattice BGK 式 (2.77)
    // 至于表面张力的计算公式可以参考 Liu Haihu 的 2012 Three-dimensional lattice Boltzmann model for immiscible two-phase flow simulations 的式(21-23)
    struct BGKITF : public LocalOperator { // Interfacial Tension Force Jonas Tolke (ITF)

        double invnu_r;
        double invnu_b;
        double rho_0; //constant reference density
        double sigma; //界面张力系数
        double limit; //用于判断相场梯度平方和的上限

        BGKITF(double nu_r, double nu_b, double sigma_, double rho_0_ = 1.0, double limit_ = 1e-8) 
            : invnu_r(1.0 / nu_r), invnu_b(1.0 / nu_b), sigma(sigma_), rho_0(rho_0_), limit(limit_) {}

        void operator()(Cell& cell, const double& rho, const double& u, const double& v, const double& phi, const double& gradx, const double& grady) {
            static double uu;
            static double gg, invggsqrt;
            uu = u * u + v * v;
            gg = gradx * gradx + grady * grady;
            if(gg > limit){
                invggsqrt = 1.0 / sqrt(gg);
            }
            else{
                invggsqrt = 0;
            }
            static double invnu, omega;
            invnu = 0.5 * ((1 + phi) * invnu_r + (1 - phi) * invnu_b);
            omega = 2.0 * invnu / (6.0 + invnu);
            for (int k = 0; k < 9; ++k) {
                static double cu;
                cu = u * Cell::c[k][0] + v * Cell::c[k][1];
                static double cc, cg;
                cc = Cell::c[k][0] * Cell::c[k][0] + Cell::c[k][1] * Cell::c[k][1];
                cg = gradx * Cell::c[k][0] + grady * Cell::c[k][1];
                static double feq;
                //使用相场梯度, 把界面张力直接引入到平衡函数中
                //feq = Cell::w[k] * (rho + rho_0 * (3.0 * cu + 4.5 * cu * cu - 1.5 * uu) - 9.0 / 4.0 * sigma * invggsqrt * (gg * (cc - 1.0 / 3.0) - cg * cg));
                feq = Cell::w[k] * (rho * (1.0 + 3.0 * cu + 4.5 * cu * cu - 1.5 * uu) + 9.0 / 4.0 * sigma * invggsqrt * (cg * cg + gg * (1.0 / 3.0 - cc)));
                cell.f[k] = (1 - omega) * cell.f[k] + omega * feq;
            }
        }

    };

}
#endif
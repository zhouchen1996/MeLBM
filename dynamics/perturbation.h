#ifndef _PERTURBATION_H_
#define _PERTURBATION_H_

#include "../core/cell.h"
#include "../core/localoperator.h"
#include <cmath>

namespace LBM {

    /*********扰动算子Perturbation*********/
    struct Perturbation : public LocalOperator {

        double A;
        double theta; //界面张力系数
        double limit; //用于判断相场梯度平方和的上限

        Perturbation(double _A, double _limit = 1e-8) 
            : A(_A), limit(_limit) {}

        void operator()(Cell& cell, const double& phaseFieldGradient_x, const double& phaseFieldGradient_y) {
            double FF = phaseFieldGradient_x * phaseFieldGradient_x + phaseFieldGradient_y * phaseFieldGradient_y;
            if (FF < limit) {
                return;
            }
            for (int k = 0; k < 9; k++) {
                double Fc = phaseFieldGradient_x * Cell::c[k][0] + phaseFieldGradient_y * Cell::c[k][1];
                cell.f[k] = cell.f[k] + 0.5 * A * sqrt(FF) * (Cell::w[k] * Fc * Fc / FF - Cell::B[k]);
            }
        }

    };




} //namespace

#endif
#ifndef _FORCE_H_
#define _FORCE_H_

#include "../core/cell.h"
#include "../core/localoperator.h"

namespace LBM {

	/**************************************** 局部算子GuoForce定义 ****************************************/

    /*
    使用Guo力项算子，来自论文"Discrete lattice effects on the forcing term in the lattice Boltzmann method"
    与BGKGuoForce不同的是，该算子并不计算bgk碰撞
    作为备用，用于一般或特殊的外力，例如重力，表面张力，拖曳力
    使用Guo力项算子后，注意计算速度时，要进行动量补偿
    */
    struct GuoForce :public LocalOperator {
        
        double omega;
		double coefficientForSourceTerm[9];
		
        GuoForce(double nu) :omega(2.0 / (6.0 * nu + 1.0)) {
			for (int k = 0; k < 9; ++k) {
				coefficientForSourceTerm[k] = (1.0 - 0.5 * omega) * Cell::w[k];
			}
		}

		//对单个晶格的分布函数进行<Guo力项>的碰撞处理
		void operator()(Cell& cell, const double& u, const double& v, const double& Fx, const double& Fy) {
			
            static double Fu;
            Fu = Fx * u + Fy * v;
			
            for (int k = 0; k < 9; ++k) {

                static double cu;
				cu = u * Cell::c[k][0] + v * Cell::c[k][1];

                static double Fc;
				Fc = Fx * Cell::c[k][0] + Fy * Cell::c[k][1];
				
				//源项(力项)
                static double S;
				S = coefficientForSourceTerm[k] * (3.0 * Fc - 3.0 * Fu + 9.0 * Fc * cu);
				
                //碰撞后的分布函数
                cell.f[k] += S;
			}
		}

	};
    
    /**************************************** 局部算子Force定义 ****************************************/

    /*
    该力项算子精度比GuoForce更低，但是更简单，也可以用于力项的描述
    与Guo力项算子不同的是，计算速度时，不要进行动量补偿
    */
    struct Force :public LocalOperator {
        
        double omega;
		double coefficientForSourceTerm[9];
		
        Force(double nu) :omega(2.0 / (6.0 * nu + 1.0)) {
			for (int k = 0; k < 9; ++k) {
				coefficientForSourceTerm[k] = (1.0 - 0.5 * omega) * Cell::w[k];
			}
		}

		//对单个晶格的分布函数进行<力项>的碰撞处理
		void operator()(Cell& cell, const double& Fx, const double& Fy) {
			
            for (int k = 0; k < 9; ++k) {

                static double Fc;
				Fc = Fx * Cell::c[k][0] + Fy * Cell::c[k][1];
				
				//源项(力项)
                static double S;
                S = coefficientForSourceTerm[k] * 3.0 * Fc;
				
                //碰撞后的分布函数
                cell.f[k] += S;
			}
		}

	};


} // namespace LBM

#endif
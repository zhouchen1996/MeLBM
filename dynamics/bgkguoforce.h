#ifndef _BGKGUOFORCE_H_
#define _BGKGUOFORCE_H_

#include "../core/cell.h"
#include "../core/localoperator.h"

namespace LBM {

	/**************************************** 局部算子BGKGuoForce定义 ****************************************/
	
    //使用BGK算子+Guo力项算子，来自论文"Discrete lattice effects on the forcing term in the lattice Boltzmann method"
	//使用Guo力项算子后，注意计算速度时，要进行动量补偿
    struct BGKGuoForce :public LocalOperator {
        
        double omega;
		double coefficientForSourceTerm[9];
		
        BGKGuoForce(double nu) :omega(2.0 / (6.0 * nu + 1.0)) {
			for (int k = 0; k < 9; ++k) {
				coefficientForSourceTerm[k] = (1.0 - 0.5 * omega) * Cell::w[k];
			}
		}

		//对单个晶格的分布函数进行<BGK+Guo力项>的碰撞处理
		void operator()(Cell& cell, const double& rho, const double& u, const double& v, const double& Fx, const double& Fy) {
			
            static double uu;
            uu = u * u + v * v;
			
            static double Fu;
            Fu = Fx * u + Fy * v;
			
            for (int k = 0; k < 9; ++k) {

				static double cu;
				cu = u * Cell::c[k][0] + v * Cell::c[k][1];

                static double Fc;
				Fc = Fx * Cell::c[k][0] + Fy * Cell::c[k][1];

                //平衡分布函数
				static double feq;
				feq = rho * Cell::w[k] * (1 + 3.0 * cu + 4.5 * cu * cu - 1.5 * uu);
				
				//源项(力项)
                static double S;
				S = coefficientForSourceTerm[k] * (3.0 * Fc - 3.0 * Fu + 9.0 * Fc * cu);
				
                //碰撞后的分布函数
                cell.f[k] = (1 - omega) * cell.f[k] + omega * feq + S;
			}
		}


	};

} // namespace LBM

#endif
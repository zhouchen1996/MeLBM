#ifndef _BGK_H_
#define _BGK_H_

#include "../core/cell.h"
#include "../core/localoperator.h"

namespace LBM {

	/**************************************** 局部算子BGK定义 ****************************************/

	//BGK算子，最简单的单相碰撞算子，使用粘滞系数nu构造
	struct BGK :public LocalOperator {
		
		double omega;

		BGK(double nu) :omega(2.0 / (6.0 * nu + 1.0)) {}

		//对单个晶格的分布函数进行<BGK>碰撞处理
		void operator()(Cell& cell, const double& rho, const double& u, const double& v) {
			
			static double uu;
			uu = u * u + v * v;
			
			for (int k = 0; k < 9; ++k) {

				static double cu;
				cu = u * Cell::c[k][0] + v * Cell::c[k][1];

				//平衡分布函数
				static double feq;
				feq = rho * Cell::w[k] * (1 + 3.0 * cu + 4.5 * cu * cu - 1.5 * uu);

				//碰撞后的分布函数
				cell.f[k] = (1 - omega) * cell.f[k] + omega * feq;
			}
		}
		
	};

} // namespace LBM

#endif // _HALFWAYBOUNCEBACK_H_
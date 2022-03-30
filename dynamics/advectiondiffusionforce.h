#ifndef _ADVECTIONDIFFUSIONFORCE_H_
#define _ADVECTIONDIFFUSIONFORCE_H_

#include "../core/cell.h"
#include "../core/localoperator.h"

/*
一般不用
*/

namespace LBM {
	
	/**************************************** 局部算子AdvectionDiffusionFirstOrderForce定义 ****************************************/
	
	struct AdvectionDiffusionFirstOrderForce :public LocalOperator {
		
		double omega;
		double coefficientForSourceTerm[9];

		AdvectionDiffusionFirstOrderForce(double nu) :omega(2.0 / (6.0 * nu + 1.0)) {
			for (int k = 0; k < 9; ++k) {
				coefficientForSourceTerm[k] = (1.0 - 0.5 * omega) * Cell::w[k];
			}
		}

		void operator()(Cell& cell, const double& rho, const double& u, const double& v, const double& Fx, const double& Fy) {
			
			static double uu;
			uu = u * u + v * v;
			
			for (int k = 0; k < 9; ++k) {
				
				static double cu;
				cu = u * Cell::c[k][0] + v * Cell::c[k][1];
				
				static double feq;
				feq = rho * Cell::w[k] * (1 + 3.0 * cu);
				
				static double Fc;
				Fc = Fx * Cell::c[k][0] + Fy * Cell::c[k][1];
				
				static double S;
				S = coefficientForSourceTerm[k] * 3.0 * Fc;
				
				cell.f[k] = (1 - omega) * cell.f[k] + omega * feq + S;
			}
		}

	};


	/**************************************** 局部算子AdvectionDiffusionFirstOrderGuoForce定义 ****************************************/
	
	struct AdvectionDiffusionFirstOrderGuoForce :public LocalOperator {
		
		double omega;
		double coefficientForSourceTerm[9];

		AdvectionDiffusionFirstOrderGuoForce(double nu) :omega(2.0 / (6.0 * nu + 1.0)) {
			for (int k = 0; k < 9; ++k) {
				coefficientForSourceTerm[k] = (1.0 - 0.5 * omega) * Cell::w[k];
			}
		}

		void operator()(Cell& cell, const double& rho, const double& u, const double& v, const double& Fx, const double& Fy) {
			
			static double uu;
			uu = u * u + v * v;
			
			static double Fu;
			Fu = Fx * u + Fy * v;
			
			for (int k = 0; k < 9; ++k) {
				
				static double cu;
				cu = u * Cell::c[k][0] + v * Cell::c[k][1];
				
				static double feq;
				feq = rho * Cell::w[k] * (1 + 3.0 * cu);
				
				static double Fc;
				Fc = Fx * Cell::c[k][0] + Fy * Cell::c[k][1];
				
				static double S;
				S = coefficientForSourceTerm[k] * (3.0 * Fc - 3.0 * Fu + 9.0 * Fc * cu);
				
				cell.f[k] = (1 - omega) * cell.f[k] + omega * feq + S;
			}
		}

	};
}//namespace

#endif

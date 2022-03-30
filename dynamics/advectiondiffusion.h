#ifndef _ADVECTIONDIFFUSION_H_
#define _ADVECTIONDIFFUSION_H_

#include "../core/cell.h"
#include "../core/localoperator.h"

namespace LBM {

	/**************************************** 局部算子AdvectionDiffusionFirstOrder定义 ****************************************/

    /*
    一阶对流扩散，虽然一阶，但是可以适用对流扩散的所有情况
    */
    struct AdvectionDiffusionFirstOrder :public LocalOperator {
        
        double omega;

        /*
        这里的nu与其对应omega与流场的不一样
		如果nu=1.0/6.0，则omega=1.0，这时就与以下文献中特殊的对流扩散情况一致：
		"Elucidating the Role of Interfacial Tension for Hydrological Properties of Two-Phase Flow in Natural Sandstone by an Improved Lattice Boltzmann Method"
        "Pore–scale analysis of supercritical CO2–brine immiscible displacement under fractional–wettability conditions"
        "Scaling of Imbibition Front Dynamics in Heterogeneous Porous Media"
        */
        AdvectionDiffusionFirstOrder(double nu) :omega(2.0 / (6.0 * nu + 1.0)) {}

        //计算在流场中对流扩散的标量场，注意：下面的rho不仅仅可以表示密度，也可以泛指一切标量
		void operator()(Cell& cell, const double& rho, const double& u, const double& v) {
			
            static double uu;
			uu = u * u + v * v;
			
            for (int k = 0; k < 9; ++k) {
				
                static double cu;
				cu = u * Cell::c[k][0] + v * Cell::c[k][1];
				
                static double feq;
				feq = rho * Cell::w[k] * (1 + 3.0 * cu);
				
                //对流扩散后的标量场分布
                cell.f[k] = (1 - omega) * cell.f[k] + omega * feq;
			}
		}

	};

    /**************************************** 局部算子AdvectionDiffusionSecondOrder定义 ****************************************/
	/*
    二阶对流扩散，虽然二阶，但是没必要使用
    参数意义与一阶的相同
    */ 
    struct AdvectionDiffusionSecondOrder :public LocalOperator {

        double omega;

		AdvectionDiffusionSecondOrder(double nu) :omega(2.0 / (6.0 * nu + 1.0)) {}
		
        //计算在流场中对流扩散的标量场，注意：下面的rho不仅仅可以表示密度，也可以泛指一切标量
        void operator()(Cell& cell, const double& rho, const double& u, const double& v) {
			
            static double uu;
			uu = u * u + v * v;
			
            for (int k = 0; k < 9; ++k) {
				
                static double cu;
				cu = u * Cell::c[k][0] + v * Cell::c[k][1];
				
                static double feq;
                //与一阶对流扩散唯一不同的地方
				feq = rho * Cell::w[k] * (1 + 3.0 * cu + 4.5 * cu * cu - 1.5 * uu);
				
                //对流扩散后的标量场分布
                cell.f[k] = (1 - omega) * cell.f[k] + omega * feq;
			}
		}

	};

} // namespace LBM

#endif
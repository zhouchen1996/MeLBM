#ifndef _BGKREGULARIZED_H_
#define _BGKREGULARIZED_H_

#include "../core/cell.h"
#include "../core/localoperator.h"

namespace LBM {

    /**************************************** 局部算子BGKRegularized定义 ****************************************/
	
	/* 参考论文
	1 "Regularized lattice Bhatnagar-Gross-Krook model for two- and three-dimensional cavity flow simulations"
	2 "General regularized boundary condition for multi-speed lattice Boltzmann models"
	3 "Generalized three-dimensional lattice Boltzmann color-gradient method for immiscible two-phase pore-scale imbibition and drainage in porous media"
	*/

	/*
	正则化用于处理边界，但是BGKRegularized有没有处理边界的能力待考察
	*/

    //正则化BGK(不包括力项)，该算子比BGK更加稳定
    struct BGKRegularized :public LocalOperator {

        double omega;
		double** cc;	//用于正则化非平衡分布函数，对于D2Q9，是一个9*9的常量矩阵

		BGKRegularized(double nu) :omega(2.0 / (6.0 * nu + 1.0)) {
			
            cc = new double* [9];
			for (int i = 0; i < 9; ++i) {
				cc[i] = new double[9];
			}
			
            for (int i = 0; i < 9; ++i) {
				for (int j = 0; j < 9; ++j) {
					double ccij = Cell::c[i][0] * Cell::c[j][0] + Cell::c[i][1] * Cell::c[j][1];
					double ccii= Cell::c[i][0] * Cell::c[i][0] + Cell::c[i][1] * Cell::c[i][1];
					double ccjj = Cell::c[j][0] * Cell::c[j][0] + Cell::c[j][1] * Cell::c[j][1];

					
					//这里的cc二维数组是数学推导来的，为了提高计算速度，避免重复计算
					//具体的方程可以参考正则化BGK的相关文献

					//cc[i][j] = ccij * ccij - Cell::Cs2 * (ccii + ccjj);
					cc[i][j] = ccij * ccij - Cell::Cs2 * (ccii + ccjj) + 2 * Cell::Cs4;
				}
			}

			/*
			
			计算 cc[i][j] = ccij * ccij - Cell::Cs2 * (ccii + ccjj) + 2 * Cell::Cs4;

			static const double cc[9][9]{
             {2. / 9, -1. / 9, -1. / 9, -1. / 9, -1. / 9, -4. / 9, -4. / 9, -4. / 9, -4. / 9},
             {-1. / 9, 5. / 9, -4. / 9, 5. / 9, -4. / 9, 2. / 9, 2. / 9, 2. / 9, 2. / 9},
             { -1. / 9, -4. / 9, 5. / 9, -4. / 9, 5. / 9, 2. / 9, 2. / 9, 2. / 9, 2. / 9},
             {-1. / 9, 5. / 9, -4. / 9, 5. / 9, -4. / 9, 2. / 9, 2. / 9, 2. / 9, 2. / 9},
             {-1. / 9, -4. / 9, 5. / 9, -4. / 9, 5. / 9, 2. / 9, 2. / 9, 2. / 9, 2. / 9},
             {-4. / 9, 2. / 9, 2. / 9, 2. / 9, 2. / 9, 26. / 9, -10. / 9, 26. / 9, -10. / 9},
             {-4. / 9, 2. / 9, 2. / 9, 2. / 9, 2. / 9, -10. / 9, 26. / 9, -10. / 9, 26. / 9},
             {-4. / 9, 2. / 9, 2. / 9, 2. / 9, 2. / 9, 26. / 9, -10. / 9, 26. / 9, -10. / 9},
             { -4. / 9, 2. / 9, 2. / 9, 2. / 9, 2. / 9, -10. / 9, 26. / 9, -10. / 9, 26. / 9} };
			
			*/

		}

		~BGKRegularized() {
			for (int i = 0; i < 9; ++i) {
				delete[] cc[i];
			}
			delete[] cc;			
		}

		//对单个晶格的分布函数进行<正则化BGK>碰撞处理
		void operator()(Cell& cell, const double& rho, const double& u, const double& v) {
			
            static double uu;
			uu = u * u + v * v;
			
            //平衡分布函数
            static double feq[9];

            //非平衡分布函数
			static double fneq[9];
			
            for (int k = 0; k < 9; ++k) {
				
                static double cu;
				cu = u * Cell::c[k][0] + v * Cell::c[k][1];
				
                feq[k] = rho * Cell::w[k] * (1 + 3.0 * cu + 4.5 * cu * cu - 1.5 * uu);

                fneq[k] = cell.f[k] - feq[k];

            }
			
            for (int k = 0; k < 9; ++k) {
				
                static double temp;
				temp = 0;
				
                for (int j = 0; j < 9; ++j) {
					temp += fneq[j] * cc[k][j];
				}
                
                //碰撞后的分布函数
				cell.f[k] = feq[k] + 4.5 * (1.0 - omega) * Cell::w[k] * temp;

			}
		}

	};

} // namespace LBM

#endif // _HALFWAYBOUNCEBACK_H_
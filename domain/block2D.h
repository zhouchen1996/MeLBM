#ifndef _BLOCK2D_H_
#define _BLOCK2D_H_

#include "../core/cell.h"
#include "../core/localoperator.h"
#include "../core/utilities.h"
#include "scalarfield2D.h"
#include "vectorfield2D.h"

namespace LBM {

    /**************************************** Block2D ****************************************/
    struct Block2D {
        
        //定义一个二级指针指向每个晶格
        Cell** pc;

        //计算域尺寸
        int nx; 
		int ny;
        //除外轮廓边界外，内部晶格区域尺寸
        int X;
        int Y;

        //不可以定义默认构造函数
        Block2D() = delete;

        //使用X,Y是最后显示的区域，不是实际的计算区域，这是为了处理边界时能有操作空间
        //计算域为0~X+2,0~Y+2,但是一般操作的内部流体晶格为1~X,1~Y,除非需要设置边界
        Block2D(int _X, int _Y) :nx(_X + 2), ny(_Y + 2), X(_X), Y(_Y) { //注意：在边界留有宽度为1的边，这个边用于边界，包括出口入口无滑移周期
			
			//先分配内存

			pc = new Cell * [nx];
			for (int i = 0; i < nx; ++i) {
				pc[i] = new Cell[ny];
			}
			
			//分配相邻节点位置，并将所有晶格内所有分布函数初始化为Cell::w

			for(int i = 0; i < nx; ++i) {
				for (int j = 0; j < ny; ++j) {
					if (i + 1 < nx) pc[i][j].neighbors[1] = &pc[i + 1][j]; 
					if (j + 1 < ny) pc[i][j].neighbors[2] = &pc[i][j + 1]; 
					if (i - 1 >= 0) pc[i][j].neighbors[3] = &pc[i - 1][j]; 
					if (j - 1 >= 0) pc[i][j].neighbors[4] = &pc[i][j - 1]; 
					if (i + 1 < nx && j + 1 < ny) pc[i][j].neighbors[5] = &pc[i + 1][j + 1]; 
					if (i - 1 >= 0 && j + 1 < ny) pc[i][j].neighbors[6] = &pc[i - 1][j + 1];
					if (i - 1 >= 0 && j - 1 >= 0) pc[i][j].neighbors[7] = &pc[i - 1][j - 1];
					if (i + 1 < nx && j - 1 >= 0) pc[i][j].neighbors[8] = &pc[i + 1][j - 1];   
				}
			}

		}

		~Block2D() {
			for (int i = 0; i < nx; ++i) {
				delete[] pc[i];
			}
			delete pc;
		}

		/********** 流 **********/
		Block2D& streaming() {
			for (int j = 0; j < Y + 2; ++j) {
				for (int i = X + 1; i > 0; --i) {
					pc[i][j].f[1] = pc[i - 1][j].f[1];
				}
			}
			for (int i = 0; i < X + 2; ++i) {
				for (int j = Y + 1; j > 0; --j) {
					pc[i][j].f[2] = pc[i][j - 1].f[2];
				}
			}
			for (int j = 0; j < Y + 2; ++j) {
				for (int i = 0; i < X + 1; ++i) {
					pc[i][j].f[3] = pc[i + 1][j].f[3];
				}
			}
			for (int i = 0; i < X + 2; ++i) {
				for (int j = 0; j < Y + 1; ++j) {
					pc[i][j].f[4] = pc[i][j + 1].f[4];
				}
			}
			for (int j = Y + 1; j > 0; --j) {
				for (int i = X + 1; i > 0; --i) {
					pc[i][j].f[5] = pc[i - 1][j - 1].f[5];
				}
			}
			for (int j = Y + 1; j > 0; --j) {
				for (int i = 0; i < X + 1; ++i) {
					pc[i][j].f[6] = pc[i + 1][j - 1].f[6];
				}
			}
			for (int j = 0; j < Y + 1; ++j) {
				for (int i = 0; i < X + 1; ++i) {
					pc[i][j].f[7] = pc[i + 1][j + 1].f[7];
				}
			}
			for (int j = 0; j < Y + 1; ++j) {
				for (int i = X + 1; i > 0; --i) {
					pc[i][j].f[8] = pc[i - 1][j + 1].f[8];
				}
			}
			return *this;
		}

		//重着色
		Block2D& recolor(Block2D& redblock, Block2D& blueblock, 
			ScalarField2D<>& reddensityfield, ScalarField2D<>& bluedensityfield, 
			ScalarField2D<>& blenddensityfield, VectorField2D<>& phaseFieldGradient, 
			int x0, int x1, int y0, int y1, double _beta = 0.99) {
			
			static double beta = _beta;
			for (int i = x0; i <= x1; ++i) {
				for (int j = y0; j <= y1; ++j) {
					static double FF;
					FF = phaseFieldGradient.pv[i][j][0] * phaseFieldGradient.pv[i][j][0] + phaseFieldGradient.pv[i][j][1] * phaseFieldGradient.pv[i][j][1];
					for (int k = 0; k < 9; k++) {
						static double Fc;
						Fc = phaseFieldGradient.pv[i][j][0] * Cell::c[k][0] + phaseFieldGradient.pv[i][j][1] * Cell::c[k][1];
						static double redDivideByblend;
						static double blueDivideByblend;
						redDivideByblend = reddensityfield.ps[i][j] / blenddensityfield.ps[i][j];
						blueDivideByblend = bluedensityfield.ps[i][j] / blenddensityfield.ps[i][j];
						redblock.pc[i][j].f[k] = redDivideByblend * pc[i][j].f[k] + beta * reddensityfield.ps[i][j] * blueDivideByblend * Cell::w[k] * cosphi(k, Fc, FF);
						blueblock.pc[i][j].f[k] = blueDivideByblend * pc[i][j].f[k] - beta *  bluedensityfield.ps[i][j] * redDivideByblend * Cell::w[k] * cosphi(k, Fc, FF);
					}
				}
			}
			return *this;
		}

        /********** 遍历全部晶格 **********/
		template <typename Functor>
		bool traversal(Functor& functor) {
			//遍历一次所有的Cell对象
			//形参建议是函数对象，用于处理Cell对象
			for (int i = 0; i < nx; ++i) {
				for (int j = 0; j < ny; ++j) {
					functor(pc[i][j]);
				}
			}
            return true;
		}

        /********** 遍历一个矩形区域 **********/
		template <typename Functor>
		bool traversal(Functor& functor, const int& x0, const int& x1, const int& y0, const int& y1) {
			//遍历一次指定矩形区域内的Cell对象
			//注意：x0, x1, y0, y1 均会遍历到，即为闭区间
            if (x0 < 0 || x0 >= nx || x1 < 0 || x1 >= nx || y0 < 0 || y0 >= ny || y1 < 0 || y1 >= ny) {
                return false;
            }
            if (x0 > x1 || y0 > y1) {
                return false;
            } 
			for (int i = x0; i <= x1; ++i) {
				for (int j = y0; j <= y1; ++j) {
					functor(pc[i][j]);
				}
			}
            return true;
		}

        /********** 在整个计算域 混合(合并) red 与 blue 的分布函数 **********/
		Block2D& combine(Block2D& redblock, Block2D& blueblock) {
            //注意 必须满足 redblock.nx == blueblock.nx && redblock.ny == blueblock.ny
			for (int i = 0; i < nx; ++i) {
				for (int j = 0; j < ny; ++j) {
					for (int k = 0; k < 9; ++k) {
						pc[i][j].f[k] = redblock.pc[i][j].f[k] + blueblock.pc[i][j].f[k];
					}
				}
			}
			return *this;
		}

        /********** 计算密度 **********/
		Block2D& computeDensityField(ScalarField2D<>& densityField) {
			for (int i = 0; i < X + 2; ++i) {
				for (int j = 0; j < Y + 2; ++j) {
					densityField.ps[i][j] = pc[i][j].calculateDensityOneCell();
				}
			}
			return *this;
		}

		/********** 计算速度 **********/
		//不建议使用
		Block2D& computeVelocityField(VectorField2D<>& velocityField) {
			for (int i = 0; i < nx; ++i) {
				for (int j = 0; j < ny; ++j) {
					velocityField.pv[i][j][0] = pc[i][j].calculateVelocityXOneCell();//会再次计算密度
					velocityField.pv[i][j][1] = pc[i][j].calculateVelocityYOneCell();//会再次计算密度
				}
			}
			return *this;
		}
		//不采用Guo力项时，建议使用
		Block2D& computeVelocityField(VectorField2D<>& velocityField, ScalarField2D<>& densityField) {
			for (int i = 0; i < nx; ++i) {
				for (int j = 0; j < ny; ++j) {
					static double invRho;
					invRho = 1.0 / densityField.ps[i][j];
					velocityField.pv[i][j][0] = invRho * pc[i][j].calculateMomentumXOneCell();
					velocityField.pv[i][j][1] = invRho * pc[i][j].calculateMomentumYOneCell();
				}
			}
			return *this;
		}
		//采用Guo力项时，必须使用
		Block2D& computeVelocityFieldGuoForce(VectorField2D<>& velocityField, ScalarField2D<>& densityField, VectorField2D<>& forceField) {
			for (int i = 0; i < nx; ++i) {
				for (int j = 0; j < ny; ++j) {
					static double invRho;
					invRho = 1.0 / densityField.ps[i][j];
					velocityField.pv[i][j][0] = invRho * (pc[i][j].calculateMomentumXOneCell() + 0.5 * forceField.pv[i][j][0]);
					velocityField.pv[i][j][1] = invRho * (pc[i][j].calculateMomentumYOneCell() + 0.5 * forceField.pv[i][j][1]);
				}
			}
			return *this;
		}

		/****************** 执行每个晶格内的算子localOperators[0] *********************/
		//第1种 SV 执行每个晶格内的算子localOperators[0] 专门用于处理简单局部的碰撞，包括bgk,固体晶格的反弹
		Block2D& executeLocalOperator(ScalarField2D<double>& densityField, VectorField2D<double>& velocityField, int level = 0) {
			for (int i = 0; i < nx; ++i) {
				for (int j = 0; j < ny; ++j) {
					(*pc[i][j].localOperators[level])(pc[i][j], densityField.ps[i][j], velocityField.pv[i][j][0], velocityField.pv[i][j][1]);
				}
			}
			return *this;
		}
		//第2种 SVV 执行每个晶格内的算子localOperators[0] 用于执行BGKGuoForce局部操作
		Block2D& executeLocalOperator(ScalarField2D<>& densityField, VectorField2D<>& velocityField, VectorField2D<>& forceField, int level = 0) {
			for (int i = 0; i < X + 2; ++i) {
				for (int j = 0; j < Y + 2; ++j) {
					(*pc[i][j].localOperators[level])(pc[i][j], densityField.ps[i][j], velocityField.pv[i][j][0], velocityField.pv[i][j][1], forceField.pv[i][j][0], forceField.pv[i][j][1]);
				}
			}
			return *this;
		}
		//第3种 SVSV 执行每个晶格内的算子localOperators[0] 用于执行BGKITF局部操作
		Block2D& executeLocalOperator(ScalarField2D<>& densityField, VectorField2D<>& velocityField, ScalarField2D<>& phaseField, VectorField2D<>& phaseFieldGradient, int level = 0) {
			for (int i = 0; i < X + 2; ++i) {
				for (int j = 0; j < Y + 2; ++j) {
					(*pc[i][j].localOperators[level])(pc[i][j], densityField.ps[i][j], velocityField.pv[i][j][0], velocityField.pv[i][j][1], phaseField.ps[i][j], phaseFieldGradient.pv[i][j][0], phaseFieldGradient.pv[i][j][1]);
				}
			}
			return *this;
		}

		/********** 执行每个晶格内的算子localOperators[1] **********/
		//第4种 V 执行每个晶格内的算子localOperators[1] 用于执行Perturbation局部操作
		Block2D& executeLocalOperator(VectorField2D<>& phaseFieldGradient, int level = 1) {
			for (int i = 1; i <= X; ++i) {
				for (int j = 1; j <= Y; ++j) {
					(*pc[i][j].localOperators[level])(pc[i][j], phaseFieldGradient.pv[i][j][0], phaseFieldGradient.pv[i][j][1]);
				}
			}
			return *this;
		}

    };

}//namespace

#endif

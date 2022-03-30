#ifndef _VECTORFIELD2D_H_
#define _VECTORFIELD2D_H_

#include "../core/cell.h"
#include "scalarfield2D.h"
#include <cmath>

namespace LBM{

	/**************************************** 向量场 ****************************************/
	template <typename T = double>
	struct VectorField2D {

        //使用三级指针指向向量场各个向量的分量
		T*** pv;

		//nx,ny含义与Block2D中的相同
		int nx;
		int ny;

		//X,Y含义与Block2D中的相同
		int X;
		int Y;

        //不可以定义默认构造函数
		VectorField2D() = delete;

        //含义 参考Block2D的构造函数
		VectorField2D(int _X, int _Y) :nx(_X + 2), ny(_Y + 2), X(_X), Y(_Y) {
			pv = new T** [nx];
			for (int i = 0; i < nx; ++i) {
				pv[i] = new T* [ny];
			}
            //每个向量只有两个分量
			for (int i = 0; i < nx; ++i) {
				for (int j = 0; j < ny; ++j) {
					pv[i][j] = new T[2];
				}
			}
		}

		~VectorField2D() {
			for (int i = 0; i < nx; ++i) {
				for (int j = 0; j < ny; ++j) {
					delete pv[i][j];
				}
			}
			for (int i = 0; i < nx; ++i) {
				delete pv[i];
			}
			delete pv;
		}

        /********** 遍历全部向量 **********/
		template <typename Functor>
		bool traversal(Functor& functor) {
			//遍历一次
			//形参为可以是函数，也可以是函数对象
			for (int i = 0; i < X + 2; ++i) {
				for (int j = 0; j < Y + 2; ++j) {
					functor(pv[i][j]);
				}
			}
            return true;
		}

        /********** 遍历一个矩形区域 **********/
		template <typename Functor>
		bool traversal(Functor& functor, int x0, int x1, int y0, int y1) {
			if (x0 < 0 || x0 >= nx || x1 < 0 || x1 >= nx || y0 < 0 || y0 >= ny || y1 < 0 || y1 >= ny) {
                return false;
            }
            if (x0 > x1 || y0 > y1) {
                return false;
            } 
            //遍历一次指定矩形区域内的数据
			//注意：x0, x1, y0, y1 均会遍历到，即为闭区间
			for (int i = x0; i <= x1; ++i) {
				for (int j = y0; j <= y1; ++j) {
					functor(pv[i][j]);
				}
			}
            return true;
		}

		/*********** 各向同性差分计算标量场的梯度 **************/
        //计算1~X 1~Y区域范围内的标量场梯度
		VectorField2D& gradient(ScalarField2D<>& scalarField, const int& x0, const int& x1, const int& y0, const int& y1) {
            //注意 尺寸必须相同
            for (int i = x0; i <= x1; ++i) {
				for (int j = y0; j <= y1; ++j) {
					//X方向
					pv[i][j][0] =
						Cell::derivativeCoefficient[0][0] * scalarField.ps[i][j] +
						Cell::derivativeCoefficient[1][0] * scalarField.ps[i + 1][j] +
						Cell::derivativeCoefficient[2][0] * scalarField.ps[i][j + 1] +
						Cell::derivativeCoefficient[3][0] * scalarField.ps[i - 1][j] +
						Cell::derivativeCoefficient[4][0] * scalarField.ps[i][j - 1] +
						Cell::derivativeCoefficient[5][0] * scalarField.ps[i + 1][j + 1] +
						Cell::derivativeCoefficient[6][0] * scalarField.ps[i - 1][j + 1] +
						Cell::derivativeCoefficient[7][0] * scalarField.ps[i - 1][j - 1] +
						Cell::derivativeCoefficient[8][0] * scalarField.ps[i + 1][j - 1];
					//Y方向
					pv[i][j][1] =
						Cell::derivativeCoefficient[0][1] * scalarField.ps[i][j] +
						Cell::derivativeCoefficient[1][1] * scalarField.ps[i + 1][j] +
						Cell::derivativeCoefficient[2][1] * scalarField.ps[i][j + 1] +
						Cell::derivativeCoefficient[3][1] * scalarField.ps[i - 1][j] +
						Cell::derivativeCoefficient[4][1] * scalarField.ps[i][j - 1] +
						Cell::derivativeCoefficient[5][1] * scalarField.ps[i + 1][j + 1] +
						Cell::derivativeCoefficient[6][1] * scalarField.ps[i - 1][j + 1] +
						Cell::derivativeCoefficient[7][1] * scalarField.ps[i - 1][j - 1] +
						Cell::derivativeCoefficient[8][1] * scalarField.ps[i + 1][j - 1];
				}
			}
			return *this;
		}

		/*********** 各向同性差分计算标量场在FB与FL的梯度并且估算SB的梯度 **************/
		//注意：这与文献中的预估SB的相场不一样
		VectorField2D& gradient(ScalarField2D<>& scalarField, ScalarField2D<AreaType>& area) {
			//注意 尺寸必须相同
            for (int i = 1; i <= X; ++i) {
				for (int j = 1; j <= Y; ++j) {
					if (area.ps[i][j] == AreaType::FB || area.ps[i][j] == AreaType::FL) {
						//X方向
						pv[i][j][0] =
						Cell::derivativeCoefficient[0][0] * scalarField.ps[i][j] +
						Cell::derivativeCoefficient[1][0] * scalarField.ps[i + 1][j] +
						Cell::derivativeCoefficient[2][0] * scalarField.ps[i][j + 1] +
						Cell::derivativeCoefficient[3][0] * scalarField.ps[i - 1][j] +
						Cell::derivativeCoefficient[4][0] * scalarField.ps[i][j - 1] +
						Cell::derivativeCoefficient[5][0] * scalarField.ps[i + 1][j + 1] +
						Cell::derivativeCoefficient[6][0] * scalarField.ps[i - 1][j + 1] +
						Cell::derivativeCoefficient[7][0] * scalarField.ps[i - 1][j - 1] +
						Cell::derivativeCoefficient[8][0] * scalarField.ps[i + 1][j - 1];
					    //Y方向
                        pv[i][j][1] =
                            Cell::derivativeCoefficient[0][1] * scalarField.ps[i][j] +
                            Cell::derivativeCoefficient[1][1] * scalarField.ps[i + 1][j] +
                            Cell::derivativeCoefficient[2][1] * scalarField.ps[i][j + 1] +
                            Cell::derivativeCoefficient[3][1] * scalarField.ps[i - 1][j] +
                            Cell::derivativeCoefficient[4][1] * scalarField.ps[i][j - 1] +
                            Cell::derivativeCoefficient[5][1] * scalarField.ps[i + 1][j + 1] +
                            Cell::derivativeCoefficient[6][1] * scalarField.ps[i - 1][j + 1] +
                            Cell::derivativeCoefficient[7][1] * scalarField.ps[i - 1][j - 1] +
                            Cell::derivativeCoefficient[8][1] * scalarField.ps[i + 1][j - 1];
					}
				}
			}
			//估算SB处的相场梯度
			for (int i = 0; i < nx; ++i) {
				for (int j = 0; j < ny; ++j) {
					if (area.ps[i][j] == AreaType::SB) {
						pv[i][j][0] = 0;
						pv[i][j][1] = 0;
						static double temp0;
						static double temp1;
						static double tempw;
						temp0 = 0;
						temp1 = 0;
						tempw = 0;
						if (i + 1 < nx && area.ps[i + 1][j] == AreaType::FB) {
							temp0 += Cell::w[1] * pv[i + 1][j][0];
							temp1 += Cell::w[1] * pv[i + 1][j][1];
							tempw += Cell::w[1];
						}						
						if (j + 1 < ny && area.ps[i][j + 1] == AreaType::FB) {
							temp0 += Cell::w[2] * pv[i][j + 1][0]; 
							temp1 += Cell::w[2] * pv[i][j + 1][1]; 
							tempw += Cell::w[2];
						}						
						if (i - 1 >= 0 && area.ps[i - 1][j] == AreaType::FB) {
							temp0 += Cell::w[3] * pv[i - 1][j][0]; 
							temp1 += Cell::w[3] * pv[i - 1][j][1]; 
							tempw += Cell::w[3];
						}
						if (j - 1 >= 0 && area.ps[i][j - 1] == AreaType::FB) {
							temp0 += Cell::w[4] * pv[i][j - 1][0]; 
							temp1 += Cell::w[4] * pv[i][j - 1][1]; 
							tempw += Cell::w[4];
						}
						if (i + 1 < nx && j + 1 < ny && area.ps[i + 1][j + 1] == AreaType::FB) {
							temp0 += Cell::w[5] * pv[i + 1][j + 1][0];
							temp1 += Cell::w[5] * pv[i + 1][j + 1][1];
							tempw += Cell::w[5];
						}
						if (i - 1 >= 0 && j + 1 < ny && area.ps[i - 1][j + 1] == AreaType::FB) {
							temp0 += Cell::w[6] * pv[i - 1][j + 1][0];
							temp1 += Cell::w[6] * pv[i - 1][j + 1][1];
							tempw += Cell::w[6];
						}
						if (i - 1 >= 0 && j - 1 >=0 && area.ps[i - 1][j - 1] == AreaType::FB) {
							temp0 += Cell::w[7] * pv[i - 1][j - 1][0];
							temp1 += Cell::w[7] * pv[i - 1][j - 1][1];
							tempw += Cell::w[7];
						}
						if (i + 1 < nx && j - 1 >= 0 && area.ps[i + 1][j - 1] == AreaType::FB) {
							temp0 += Cell::w[8] * pv[i + 1][j - 1][0];
							temp1 += Cell::w[8] * pv[i + 1][j - 1][1];
							tempw += Cell::w[8];
						}
						pv[i][j][0] = temp0 / tempw;
						pv[i][j][1] = temp1 / tempw;
					}
				}
			}
			return *this;
		}
        
		/***************** 归一化 *****************/
		VectorField2D& normalize(VectorField2D& phaseFieldGradient, const double& limit = 1e-8) {
			for (int i = 1; i <= X; ++i) {
				for (int j = 1; j <= Y; ++j) {
					double norm = sqrt(phaseFieldGradient.pv[i][j][0] * phaseFieldGradient.pv[i][j][0] + phaseFieldGradient.pv[i][j][1] * phaseFieldGradient.pv[i][j][1]);
					if (norm < limit) {
						pv[i][j][0] = pv[i][j][1] = 0;
						continue;
					}
					pv[i][j][0] = phaseFieldGradient.pv[i][j][0] / norm;
					pv[i][j][1] = phaseFieldGradient.pv[i][j][1] / norm;
				}
			}
			return *this;
		}

		/*************** 归一化在FB与FL的梯度 *****************/
		VectorField2D& normalize(VectorField2D& phaseFieldGradient, ScalarField2D<AreaType>& area, const double& limit = 1e-8) {
			for (int i = 1; i <= X; ++i) {
				for (int j = 1; j <= Y; ++j) {
					if (area.ps[i][j] == AreaType::FB || area.ps[i][j] == AreaType::FL) {
						double norm = sqrt(phaseFieldGradient.pv[i][j][0] * phaseFieldGradient.pv[i][j][0] + phaseFieldGradient.pv[i][j][1] * phaseFieldGradient.pv[i][j][1]);
						if (norm < limit) {
							pv[i][j][0] = pv[i][j][1] = 0;
							continue;
						}
						pv[i][j][0] = phaseFieldGradient.pv[i][j][0] / norm;
						pv[i][j][1] = phaseFieldGradient.pv[i][j][1] / norm;
					}
				}
			}
			return *this;
		}

		/************* 表面张力力 ****************/
		VectorField2D& SurfaceTensionForce(VectorField2D& phaseFieldGradient, const double& theta) {
			
			return *this;
		}
	};

}//namespace

#endif
#ifndef _SCALARFIELD2D_H_
#define _SCALARFIELD2D_H_

#include "../core/cell.h"

namespace LBM{

	/**************************************** 标量场 ****************************************/
	template <typename T = double>
	struct ScalarField2D{

		//定义一个二级指针指向每个标量
		T** ps;

		//nx,ny含义与Block2D中的相同
		int nx;
		int ny;

		//X,Y含义与Block2D中的相同
		int X;
		int Y;

		//不可以定义默认构造函数
		ScalarField2D() = delete;

		//含义 参考Block2D的构造函数
		ScalarField2D(int _X, int _Y) :nx(_X + 2), ny(_Y + 2), X(_X), Y(_Y) {
			ps = new T* [nx];
			for (int i = 0; i < nx; ++i) {
				ps[i] = new T[ny];
			}
		}

		~ScalarField2D(){
			for (int i = 0; i < nx; ++i) {
				delete ps[i];
			}
			delete ps;
		}

		/********** 遍历全部标量 **********/
		template <typename Functor>
		bool traversal(Functor& functor) {
			//遍历一次
			//形参为可以是函数，也可以是函数对象
			for (int i = 0; i < nx; ++i) {
				for (int j = 0; j < ny; ++j) {
					functor(ps[i][j]);
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
					functor(ps[i][j]);
				}
			}
			return true;
		}

		/********** 标量场相减 **********/
		// this -> ps[i][j] = scalarField1.ps[i][j] - scalarField2.ps[i][j]
		ScalarField2D& minus(ScalarField2D& scalarField1, ScalarField2D& scalarField2) {
			//注意 必须满足 scalarField1.nx == scalarField2.nx && scalarField1.ny == scalarField2.ny
			for (int i = 0; i < nx; ++i) {
				for (int j = 0; j < ny; ++j) {
					ps[i][j] = scalarField1.ps[i][j] - scalarField2.ps[i][j];
				}
			}
			return *this;
		}
		
		/********** 矩形区域内标量场相减 **********/
		// this -> ps[i][j] = scalarField1.ps[i][j] - scalarField2.ps[i][j]
		ScalarField2D& minus(ScalarField2D& scalarField1, ScalarField2D& scalarField2, const int& x0, const int& x1, const int& y0, const int& y1) {
			//注意 必须满足 scalarField1.nx == scalarField2.nx && scalarField1.ny == scalarField2.ny
			for (int i = x0; i <= x1; ++i) {
				for (int j = y0; j <= y1; ++j) {
					ps[i][j] = scalarField1.ps[i][j] - scalarField2.ps[i][j];
				}
			}
			return *this;
		}
		
		/********** 标量场相加 **********/
		// this -> ps[i][j] = scalarField1.ps[i][j] + scalarField2.ps[i][j]
		ScalarField2D& plus(ScalarField2D& scalarField1, ScalarField2D& scalarField2) {
			//注意 必须满足 scalarField1.nx == scalarField2.nx && scalarField1.ny == scalarField2.ny
			for (int i = 0; i < nx; ++i) {
				for (int j = 0; j < ny; ++j) {
					ps[i][j] = scalarField1.ps[i][j] + scalarField2.ps[i][j];
				}
			}
			return *this;
		}

		/********** 矩形区域内标量场相加 **********/
		// this -> ps[i][j] = scalarField1.ps[i][j] + scalarField2.ps[i][j]
		ScalarField2D& plus(ScalarField2D& scalarField1, ScalarField2D& scalarField2, const int& x0, const int& x1, const int& y0, const int& y1) {
			//注意 必须满足 scalarField1.nx == scalarField2.nx && scalarField1.ny == scalarField2.ny
			for (int i = x0; i <= x1; ++i) {
				for (int j = y0; j <= y1; ++j) {
					ps[i][j] = scalarField1.ps[i][j] + scalarField2.ps[i][j];
				}
			}
			return *this;
		}

		/************* 求相场 *************/
		ScalarField2D& computePhaseField(ScalarField2D& reddensityfield, ScalarField2D& bluedensityfield, 
			ScalarField2D& blenddensityfield, int x0, int x1, int y0, int y1) {
			//注意 尺寸必须相同
			for (int i = x0; i <= x1; ++i) {
				for (int j = y0; j <= y1; ++j) {
					ps[i][j] = (reddensityfield.ps[i][j] - bluedensityfield.ps[i][j]) / blenddensityfield.ps[i][j];
				}
			}
			return *this;
		}

		//只在流体晶格内求相场
		ScalarField2D& computePhaseField(ScalarField2D& reddensityfield, ScalarField2D& bluedensityfield, 
			ScalarField2D& blenddensityfield, ScalarField2D<AreaType>& area) {
			//注意 尺寸必须相同
			for (int i = 0; i < nx; ++i) {
				for (int j = 0; j < ny; ++j) {
					if (area.ps[i][j] == AreaType::FL || area.ps[i][j] == AreaType::FB) {
						ps[i][j] = (reddensityfield.ps[i][j] - bluedensityfield.ps[i][j]) / blenddensityfield.ps[i][j];
					}
				}
			}
			return *this;
		}

	};

}//namespace

#endif
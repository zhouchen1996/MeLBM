#ifndef _CELL_H_
#define _CELL_H_

namespace LBM {

    struct LocalOperator;

	enum class AreaType { F, FB, FL, S, SB, SL};

	//D2Q9晶格,包含分布函数与局部算子
	struct Cell {

		//分布函数
		double f[9]{ 0 };

		//晶格所在的区域,默认为F,即流体晶格
		AreaType areatype = AreaType::F;

		//访问相邻的晶格,按下标访问的晶格与f[1...8]方向对应
		Cell* neighbors[8]{nullptr}; 

        /*
        //局部算子
		//预设2个，一般 localOperator[0] 为单相碰撞算子，localOperator[1] 为扰动算子
        */
		LocalOperator* localOperators[2]{nullptr};

        //分布函数的各个分量的权重
		static const double w[9];

        //D2Q9的速度集
		static const double c[9][2];

        //D2Q9颜色梯度模型的扰动算子中的固定参数
		static const double B[9];
		
        //晶格声速的平方
        static const double Cs2;

        //晶格声速的四次方
		static const double Cs4;

        //晶格声速的平方的倒数
		static const double InvCs2;

		//使用各项同性差分求解一个标量的导数或者梯度的系数
		static const double derivativeCoefficient[9][2];

		/*
        //以下成员函数均是用于局部操作，即，只操作单个晶格的分布函数，不涉及邻近的晶格
        */

        //计算密度 rho
		double calculateDensityOneCell() {
			return f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7] + f[8];
		}

        //计算x方向上的动量 rho * u_x
		double calculateMomentumXOneCell() {
			return f[1] + f[5] + f[8] - f[3] - f[6] - f[7];
		}

        //计算y方向上的动量 rho * u_y
		double calculateMomentumYOneCell() {
			return f[2] + f[5] + f[6] - f[4] - f[7] - f[8];
		}

        //计算x方向上的速度 u_x
		double calculateVelocityXOneCell() {
			return (f[1] + f[5] + f[8] - f[3] - f[6] - f[7]) / (f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7] + f[8]);
		}

        //计算y方向上的速度 u_y
		double calculateVelocityYOneCell() {
			return (f[2] + f[5] + f[6] - f[4] - f[7] - f[8]) / (f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7] + f[8]);
		}

	};

	//分布函数的各个分量的权重
	const double Cell::w[9]{ 4. / 9,1. / 9,1. / 9 ,1. / 9 ,1. / 9 ,1. / 36,1. / 36 ,1. / 36 ,1. / 36 };

	//D2Q9的速度集
	const double Cell::c[9][2]{ { 0, 0 }, { 1, 0 }, { 0, 1 },{ -1, 0 }, { 0, -1 }, { 1, 1 },{ -1, 1 } ,{ -1, -1 }, { 1, -1 } };

	//D2Q9颜色梯度模型的扰动算子中的固定参数
	const double Cell::B[9]{ -4.0 / 27.0, 2.0 / 27.0, 2.0 / 27.0, 2.0 / 27.0, 2.0 / 27.0, 5.0 / 108.0, 5.0 / 108.0, 5.0 / 108.0, 5.0 / 108.0 };
	
	//晶格声速的平方
	const double Cell::Cs2 = 1.0 / 3.0;

	//晶格声速的四次方
	const double Cell::Cs4 = 1.0 / 9.0;

	//晶格声速的平方的倒数
	const double Cell::InvCs2 = 3.0;

	//使用各项同性差分求解一个标量的导数或者梯度的系数
	const double Cell::derivativeCoefficient[9][2] = {	
			{Cell::InvCs2 * Cell::w[0] * Cell::c[0][0], Cell::InvCs2 * Cell::w[0] * Cell::c[0][1]},
			{Cell::InvCs2 * Cell::w[1] * Cell::c[1][0], Cell::InvCs2 * Cell::w[1] * Cell::c[1][1]},
			{Cell::InvCs2 * Cell::w[2] * Cell::c[2][0], Cell::InvCs2 * Cell::w[2] * Cell::c[2][1]},
			{Cell::InvCs2 * Cell::w[3] * Cell::c[3][0], Cell::InvCs2 * Cell::w[3] * Cell::c[3][1]},
			{Cell::InvCs2 * Cell::w[4] * Cell::c[4][0], Cell::InvCs2 * Cell::w[4] * Cell::c[4][1]},
			{Cell::InvCs2 * Cell::w[5] * Cell::c[5][0], Cell::InvCs2 * Cell::w[5] * Cell::c[5][1]},
			{Cell::InvCs2 * Cell::w[6] * Cell::c[6][0], Cell::InvCs2 * Cell::w[6] * Cell::c[6][1]},
			{Cell::InvCs2 * Cell::w[7] * Cell::c[7][0], Cell::InvCs2 * Cell::w[7] * Cell::c[7][1]},
			{Cell::InvCs2 * Cell::w[8] * Cell::c[8][0], Cell::InvCs2 * Cell::w[8] * Cell::c[8][1]} };

}//namespace

#endif
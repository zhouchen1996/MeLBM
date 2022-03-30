#ifndef _REGULARIZEBOUNDARY_H_
#define _REGULARIZEBOUNDARY_H_

#include "../core/cell.h"
#include "../core/localoperator.h"

namespace LBM {

    //边界的方位  X，Y表示边界外法线方向沿着X, Y方向  P表示positive正向  N表示negative负向
    //这个枚举类型并不使用，只是为了说明下面各个类的后缀含义
    enum class BoundaryOrientation { EdgeXP, EdgeXN, EdgeYP, EdgeYN, CornerPP, CornerPN, CornerNP, CornerNN };
    
    //列举了16种情况，包含了所有压力速度边界情况

    /**************************************** 派生类RegularizedBoundary定义 ****************************************/
    class RegularizedBoundary :public LocalOperator {
    public:
        static const double cc[9][9];   //用于正则化非平衡分布函数
    };

    const double RegularizedBoundary::cc[9][9] = {   //用于正则化非平衡分布函数
             {2. / 9, -1. / 9, -1. / 9, -1. / 9, -1. / 9, -4. / 9, -4. / 9, -4. / 9, -4. / 9},
             {-1. / 9, 5. / 9, -4. / 9, 5. / 9, -4. / 9, 2. / 9, 2. / 9, 2. / 9, 2. / 9},
             { -1. / 9, -4. / 9, 5. / 9, -4. / 9, 5. / 9, 2. / 9, 2. / 9, 2. / 9, 2. / 9},
             {-1. / 9, 5. / 9, -4. / 9, 5. / 9, -4. / 9, 2. / 9, 2. / 9, 2. / 9, 2. / 9},
             {-1. / 9, -4. / 9, 5. / 9, -4. / 9, 5. / 9, 2. / 9, 2. / 9, 2. / 9, 2. / 9},
             {-4. / 9, 2. / 9, 2. / 9, 2. / 9, 2. / 9, 26. / 9, -10. / 9, 26. / 9, -10. / 9},
             {-4. / 9, 2. / 9, 2. / 9, 2. / 9, 2. / 9, -10. / 9, 26. / 9, -10. / 9, 26. / 9},
             {-4. / 9, 2. / 9, 2. / 9, 2. / 9, 2. / 9, 26. / 9, -10. / 9, 26. / 9, -10. / 9},
             { -4. / 9, 2. / 9, 2. / 9, 2. / 9, 2. / 9, -10. / 9, 26. / 9, -10. / 9, 26. / 9} };

    /**************************************** 派生类RegularizedVelocityBoundaryEdgeXN定义 ****************************************/
    class RegularizedVelocityBoundaryEdgeXN :public RegularizedBoundary {
    public:
        RegularizedVelocityBoundaryEdgeXN(double nu, double ux, double uy)
            :omega(2.0 / (6.0 * nu + 1.0)) {
            velocity[0] = ux;
            velocity[1] = uy;
        }
        //后面三个形参的作用是可以使用第1种executeLocalOperator 
        void operator()(Cell& cell, const double& a, const double& b, const double& c) {
            double u = velocity[0];
            double v = velocity[1];
            //根据方向修改
            double rho = (cell.f[0] + cell.f[2] + cell.f[4] + 2.0 * (cell.f[3] + cell.f[6] + cell.f[7])) / (1.0 - u);
            static double uu;
            uu = u * u + v * v;
            static double feq[9];
            for (int k = 0; k < 9; ++k) {
                static double cu;
                cu = u * Cell::c[k][0] + v * Cell::c[k][1];
                feq[k] = rho * Cell::w[k] * (1 + 3.0 * cu + 4.5 * cu * cu - 1.5 * uu);
            }
            //根据方向修改
            static double fneq[9];
            fneq[0] = cell.f[0] - feq[0];
            fneq[2] = cell.f[2] - feq[2];
            fneq[3] = cell.f[3] - feq[3];
            fneq[4] = cell.f[4] - feq[4];
            fneq[6] = cell.f[6] - feq[6];
            fneq[7] = cell.f[7] - feq[7];
            //未知的分布函数均采用非平衡反弹得到；这是文献中采取的方法，
            //但是个人认为，没有考虑rho的守恒性，因为非平衡函数的总和应该为0
            //所以，这里计算未知分布函数的非平衡部分属于预处理，后面的正则化是否确保了rho守恒还有待验证
            //不过，从模拟的结果来看，这么做反而确保了相当程度的稳定
            fneq[1] = fneq[3];
            fneq[5] = fneq[7];
            fneq[8] = fneq[6];
            for (int k = 0; k < 9; ++k) {
                static double temp;
                temp = 0;
                for (int j = 0; j < 9; ++j) {
                    temp += fneq[j] * cc[k][j];
                }
                cell.f[k] = feq[k] + 4.5 * (1.0 - omega) * Cell::w[k] * temp;
            }
        }
    private:
        double omega;
        double velocity[2];
    };
    /**************************************** 派生类RegularizedVelocityBoundaryEdgeXP定义 ****************************************/
    class RegularizedVelocityBoundaryEdgeXP :public RegularizedBoundary {
    public:

        RegularizedVelocityBoundaryEdgeXP(double nu, double ux, double uy)
            :omega(2.0 / (6.0 * nu + 1.0)) {
            velocity[0] = ux;
            velocity[1] = uy;
        }
        //后面三个形参的作用是可以使用 第1种 executeLocalOperator 
        void operator()(Cell& cell, const double& a, const double& b, const double& c) {
            double u = velocity[0];
            double v = velocity[1];
            //根据方向修改
            double rho = (cell.f[0] + cell.f[2] + cell.f[4] + 2.0 * (cell.f[1] + cell.f[5] + cell.f[8])) / (1.0 + u);
            static double uu;
            uu = u * u + v * v;
            static double feq[9];
            for (int k = 0; k < 9; ++k) {
                static double cu;
                cu = u * Cell::c[k][0] + v * Cell::c[k][1];
                feq[k] = rho * Cell::w[k] * (1 + 3.0 * cu + 4.5 * cu * cu - 1.5 * uu);
            }
            //根据方向修改
            static double fneq[9];
            fneq[0] = cell.f[0] - feq[0];
            fneq[2] = cell.f[2] - feq[2];
            fneq[4] = cell.f[4] - feq[4];
            fneq[1] = cell.f[1] - feq[1];
            fneq[5] = cell.f[5] - feq[5];
            fneq[8] = cell.f[8] - feq[8];
            //未知的分布函数均采用非平衡反弹得到
            fneq[3] = fneq[1];
            fneq[6] = fneq[8];
            fneq[7] = fneq[5];
            for (int k = 0; k < 9; ++k) {
                static double temp;
                temp = 0;
                for (int j = 0; j < 9; ++j) {
                    temp += fneq[j] * cc[k][j];
                }
                cell.f[k] = feq[k] + 4.5 * (1.0 - omega) * Cell::w[k] * temp;
            }
        }
    private:
        double omega;
        double velocity[2];
    };
    /**************************************** 派生类RegularizedVelocityBoundaryEdgeYP定义 ****************************************/
    class RegularizedVelocityBoundaryEdgeYP :public RegularizedBoundary {
    public:
        RegularizedVelocityBoundaryEdgeYP(double nu, double ux, double uy)
            :omega(2.0 / (6.0 * nu + 1.0)) {
            velocity[0] = ux;
            velocity[1] = uy;
        }
        //后面三个形参的作用是可以使用 第1种 executeLocalOperator 
        void operator()(Cell& cell, const double& a, const double& b, const double& c) {
            double u = velocity[0];
            double v = velocity[1];
            //根据方向修改
            double rho = (cell.f[0] + cell.f[1] + cell.f[3] + 2.0 * (cell.f[2] + cell.f[5] + cell.f[6])) / (1.0 + v);
            static double uu;
            uu = u * u + v * v;
            static double feq[9];
            for (int k = 0; k < 9; ++k) {
                static double cu;
                cu = u * Cell::c[k][0] + v * Cell::c[k][1];
                feq[k] = rho * Cell::w[k] * (1 + 3.0 * cu + 4.5 * cu * cu - 1.5 * uu);
            }
            //根据方向修改
            static double fneq[9];
            fneq[0] = cell.f[0] - feq[0];
            fneq[1] = cell.f[1] - feq[1];
            fneq[2] = cell.f[2] - feq[2];
            fneq[3] = cell.f[3] - feq[3];
            fneq[5] = cell.f[5] - feq[5];
            fneq[6] = cell.f[6] - feq[6];
            //未知的分布函数均采用非平衡反弹得到
            fneq[4] = fneq[4];
            fneq[7] = fneq[7];
            fneq[8] = fneq[8];
            for (int k = 0; k < 9; ++k) {
                static double temp;
                temp = 0;
                for (int j = 0; j < 9; ++j) {
                    temp += fneq[j] * cc[k][j];
                }
                cell.f[k] = feq[k] + 4.5 * (1.0 - omega) * Cell::w[k] * temp;
            }
        }
    private:
        double omega;
        double velocity[2];
    };
    /**************************************** 派生类RegularizedVelocityBoundaryEdgeYN定义 ****************************************/
    class RegularizedVelocityBoundaryEdgeYN :public RegularizedBoundary {
    public:
        RegularizedVelocityBoundaryEdgeYN(double nu, double ux, double uy)
            :omega(2.0 / (6.0 * nu + 1.0)) {
            velocity[0] = ux;
            velocity[1] = uy;
        }
        //后面三个形参的作用是可以使用 第1种 executeLocalOperator 
        void operator()(Cell& cell, const double& a, const double& b, const double& c) {
            double u = velocity[0];
            double v = velocity[1];
            //根据方向修改
            double rho = (cell.f[0] + cell.f[1] + cell.f[3] + 2.0 * (cell.f[4] + cell.f[7] + cell.f[8])) / (1.0 - v);
            static double uu;
            uu = u * u + v * v;
            static double feq[9];
            for (int k = 0; k < 9; ++k) {
                static double cu;
                cu = u * Cell::c[k][0] + v * Cell::c[k][1];
                feq[k] = rho * Cell::w[k] * (1 + 3.0 * cu + 4.5 * cu * cu - 1.5 * uu);
            }
            //根据方向修改
            static double fneq[9];
            fneq[0] = cell.f[0] - feq[0];
            fneq[1] = cell.f[1] - feq[1];
            fneq[2] = cell.f[2] - feq[2];
            fneq[4] = cell.f[4] - feq[4];
            fneq[7] = cell.f[7] - feq[7];
            fneq[8] = cell.f[8] - feq[8];
            //未知的分布函数均采用非平衡反弹得到
            fneq[2] = fneq[2];
            fneq[5] = fneq[5];
            fneq[6] = fneq[6];
            for (int k = 0; k < 9; ++k) {
                static double temp;
                temp = 0;
                for (int j = 0; j < 9; ++j) {
                    temp += fneq[j] * cc[k][j];
                }
                cell.f[k] = feq[k] + 4.5 * (1.0 - omega) * Cell::w[k] * temp;
            }
        }
    private:
        double omega;
        double velocity[2];
    };
    /**************************************** 派生类RegularizedVelocityBoundaryCornerPP定义 ****************************************/
    class RegularizedVelocityBoundaryCornerPP :public RegularizedBoundary {
    public:
        RegularizedVelocityBoundaryCornerPP(double nu, double ux, double uy)
            :omega(2.0 / (6.0 * nu + 1.0)) {
            velocity[0] = ux;
            velocity[1] = uy;
        }
        //后面三个形参的作用是可以使用 第1种 executeLocalOperator 
        void operator()(Cell& cell, const double& a, const double& b, const double& c) {
            double u = velocity[0];
            double v = velocity[1];
            //根据方向修改
            double rho = cell.f[0] + cell.f[6] + cell.f[8] + 2.0 * (cell.f[1] + cell.f[2] + cell.f[5]);
            static double uu;
            uu = u * u + v * v;
            static double feq[9];
            for (int k = 0; k < 9; ++k) {
                static double cu;
                cu = u * Cell::c[k][0] + v * Cell::c[k][1];
                feq[k] = rho * Cell::w[k] * (1 + 3.0 * cu + 4.5 * cu * cu - 1.5 * uu);
            }
            //根据方向修改
            static double fneq[9];
            fneq[0] = cell.f[0] - feq[0];
            fneq[1] = cell.f[1] - feq[1];
            fneq[2] = cell.f[2] - feq[2];
            fneq[5] = cell.f[5] - feq[5];
            fneq[6] = cell.f[6] - feq[6];
            fneq[8] = cell.f[8] - feq[8];
            //未知的分布函数均采用非平衡反弹得到
            fneq[3] = fneq[3];
            fneq[4] = fneq[4];
            fneq[7] = fneq[7];
            for (int k = 0; k < 9; ++k) {
                static double temp;
                temp = 0;
                for (int j = 0; j < 9; ++j) {
                    temp += fneq[j] * cc[k][j];
                }
                cell.f[k] = feq[k] + 4.5 * (1.0 - omega) * Cell::w[k] * temp;
            }
        }
    private:
        double omega;
        double velocity[2];
    };
    /**************************************** 派生类RegularizedVelocityBoundaryCornerPN定义 ****************************************/
    class RegularizedVelocityBoundaryCornerPN :public RegularizedBoundary {
    public:
        RegularizedVelocityBoundaryCornerPN(double nu, double ux, double uy)
            :omega(2.0 / (6.0 * nu + 1.0)) {
            velocity[0] = ux;
            velocity[1] = uy;
        }
        //后面三个形参的作用是可以使用 第1种 executeLocalOperator 
        void operator()(Cell& cell, const double& a, const double& b, const double& c) {
            double u = velocity[0];
            double v = velocity[1];
            //根据方向修改
            double rho = cell.f[0] + cell.f[5] + cell.f[7] + 2.0 * (cell.f[1] + cell.f[4] + cell.f[8]);
            static double uu;
            uu = u * u + v * v;
            static double feq[9];
            for (int k = 0; k < 9; ++k) {
                static double cu;
                cu = u * Cell::c[k][0] + v * Cell::c[k][1];
                feq[k] = rho * Cell::w[k] * (1 + 3.0 * cu + 4.5 * cu * cu - 1.5 * uu);
            }
            //根据方向修改
            static double fneq[9];
            fneq[0] = cell.f[0] - feq[0];
            fneq[1] = cell.f[1] - feq[1];
            fneq[4] = cell.f[4] - feq[4];
            fneq[5] = cell.f[5] - feq[5];
            fneq[7] = cell.f[7] - feq[7];
            fneq[8] = cell.f[8] - feq[8];
            //未知的分布函数均采用非平衡反弹得到
            fneq[2] = fneq[2];
            fneq[3] = fneq[3];
            fneq[6] = fneq[6];
            for (int k = 0; k < 9; ++k) {
                static double temp;
                temp = 0;
                for (int j = 0; j < 9; ++j) {
                    temp += fneq[j] * cc[k][j];
                }
                cell.f[k] = feq[k] + 4.5 * (1.0 - omega) * Cell::w[k] * temp;
            }
        }
    private:
        double omega;
        double velocity[2];
    };
    /**************************************** 派生类RegularizedVelocityBoundaryCornerNP定义 ****************************************/
    class RegularizedVelocityBoundaryCornerNP :public RegularizedBoundary {
    public:
        RegularizedVelocityBoundaryCornerNP(double nu, double ux, double uy)
            :omega(2.0 / (6.0 * nu + 1.0)) {
            velocity[0] = ux;
            velocity[1] = uy;
        }
        //后面三个形参的作用是可以使用 第1种 executeLocalOperator 
        void operator()(Cell& cell, const double& a, const double& b, const double& c) {
            double u = velocity[0];
            double v = velocity[1];
            //根据方向修改
            double rho = cell.f[0] + cell.f[5] + cell.f[7] + 2.0 * (cell.f[2] + cell.f[3] + cell.f[6]);
            static double uu;
            uu = u * u + v * v;
            static double feq[9];
            for (int k = 0; k < 9; ++k) {
                static double cu;
                cu = u * Cell::c[k][0] + v * Cell::c[k][1];
                feq[k] = rho * Cell::w[k] * (1 + 3.0 * cu + 4.5 * cu * cu - 1.5 * uu);
            }
            //根据方向修改
            static double fneq[9];
            fneq[0] = cell.f[0] - feq[0];
            fneq[2] = cell.f[2] - feq[2];
            fneq[3] = cell.f[3] - feq[3];
            fneq[5] = cell.f[5] - feq[5];
            fneq[6] = cell.f[6] - feq[6];
            fneq[7] = cell.f[7] - feq[7];
            //未知的分布函数均采用非平衡反弹得到
            fneq[1] = fneq[1];
            fneq[4] = fneq[4];
            fneq[8] = fneq[8];
            for (int k = 0; k < 9; ++k) {
                static double temp;
                temp = 0;
                for (int j = 0; j < 9; ++j) {
                    temp += fneq[j] * cc[k][j];
                }
                cell.f[k] = feq[k] + 4.5 * (1.0 - omega) * Cell::w[k] * temp;
            }
        }
    private:
        double omega;
        double velocity[2];
    };
    /**************************************** 派生类RegularizedVelocityBoundaryCornerNN定义 ****************************************/
    class RegularizedVelocityBoundaryCornerNN :public RegularizedBoundary {
    public:
        RegularizedVelocityBoundaryCornerNN(double nu, double ux, double uy)
            :omega(2.0 / (6.0 * nu + 1.0)) {
            velocity[0] = ux;
            velocity[1] = uy;
        }
        //后面三个形参的作用是可以使用 第1种 executeLocalOperator 
        void operator()(Cell& cell, const double& a, const double& b, const double& c) {
            double u = velocity[0];
            double v = velocity[1];
            //根据方向修改
            double rho = cell.f[0] + cell.f[6] + cell.f[8] + 2.0 * (cell.f[3] + cell.f[4] + cell.f[7]);
            static double uu;
            uu = u * u + v * v;
            static double feq[9];
            for (int k = 0; k < 9; ++k) {
                static double cu;
                cu = u * Cell::c[k][0] + v * Cell::c[k][1];
                feq[k] = rho * Cell::w[k] * (1 + 3.0 * cu + 4.5 * cu * cu - 1.5 * uu);
            }
            //根据方向修改
            static double fneq[9];
            fneq[0] = cell.f[0] - feq[0];
            fneq[3] = cell.f[3] - feq[3];
            fneq[4] = cell.f[4] - feq[4];
            fneq[6] = cell.f[6] - feq[6];
            fneq[7] = cell.f[7] - feq[7];
            fneq[8] = cell.f[8] - feq[8];
            //未知的分布函数均采用非平衡反弹得到
            fneq[1] = fneq[1];
            fneq[2] = fneq[2];
            fneq[5] = fneq[5];
            for (int k = 0; k < 9; ++k) {
                static double temp;
                temp = 0;
                for (int j = 0; j < 9; ++j) {
                    temp += fneq[j] * cc[k][j];
                }
                cell.f[k] = feq[k] + 4.5 * (1.0 - omega) * Cell::w[k] * temp;
            }
        }
    private:
        double omega;
        double velocity[2];
    };
    /**************************************** 派生类RegularizedPressureBoundaryEdgeXN定义 ****************************************/
    class RegularizedPressureBoundaryEdgeXN :public RegularizedBoundary {
    public:
        RegularizedPressureBoundaryEdgeXN(double nu, double pressure)
            :omega(2.0 / (6.0 * nu + 1.0)) {
            rho = 3.0 * pressure;
        }
        //后面三个形参的作用是可以使用 第1种 executeLocalOperator 
        void operator()(Cell& cell, const double& a, const double& b, const double& c) {
            //根据方向修改
            double u = 1.0 - (cell.f[0] + cell.f[2] + cell.f[4] + 2.0 * (cell.f[3] + cell.f[6] + cell.f[7])) / rho;
            double v = cell.f[2] + cell.f[6] - cell.f[4] - cell.f[7];   //压力边界下，平行于边界的流速可以是0，但是这里我自己并不把它直接设为0
            static double uu;
            uu = u * u + v * v;
            static double feq[9];
            for (int k = 0; k < 9; ++k) {
                static double cu;
                cu = u * Cell::c[k][0] + v * Cell::c[k][1];
                feq[k] = rho * Cell::w[k] * (1 + 3.0 * cu + 4.5 * cu * cu - 1.5 * uu);
            }
            //根据方向修改
            static double fneq[9];
            fneq[0] = cell.f[0] - feq[0];
            fneq[2] = cell.f[2] - feq[2];
            fneq[3] = cell.f[3] - feq[3];
            fneq[4] = cell.f[4] - feq[4];
            fneq[6] = cell.f[6] - feq[6];
            fneq[7] = cell.f[7] - feq[7];
            fneq[1] = fneq[3];
            fneq[5] = fneq[7];
            fneq[8] = fneq[6];
            for (int k = 0; k < 9; ++k) {
                static double temp;
                temp = 0;
                for (int j = 0; j < 9; ++j) {
                    temp += fneq[j] * cc[k][j];
                }
                cell.f[k] = feq[k] + 4.5 * (1.0 - omega) * Cell::w[k] * temp;
            }
        }
    private:
        double omega;
        double rho;
    };
    /**************************************** 派生类RegularizedPressureBoundaryEdgeXP定义 ****************************************/
    class RegularizedPressureBoundaryEdgeXP :public RegularizedBoundary {
    public:
        RegularizedPressureBoundaryEdgeXP(double nu, double pressure)
            :omega(2.0 / (6.0 * nu + 1.0)) {
            rho = 3.0 * pressure;
        }
        //后面三个形参的作用是可以使用 第1种 executeLocalOperator 
        void operator()(Cell& cell, const double& a, const double& b, const double& c) {
            //根据方向修改
            double u = (cell.f[0] + cell.f[2] + cell.f[4] + 2.0 * (cell.f[1] + cell.f[5] + cell.f[8])) / rho - 1.0;
            double v = cell.f[2] + cell.f[5] - cell.f[4] - cell.f[8];
            static double uu;
            uu = u * u + v * v;
            static double feq[9];
            for (int k = 0; k < 9; ++k) {
                static double cu;
                cu = u * Cell::c[k][0] + v * Cell::c[k][1];
                feq[k] = rho * Cell::w[k] * (1 + 3.0 * cu + 4.5 * cu * cu - 1.5 * uu);
            }
            //根据方向修改
            static double fneq[9];
            fneq[0] = cell.f[0] - feq[0];
            fneq[2] = cell.f[2] - feq[2];
            fneq[4] = cell.f[4] - feq[4];
            fneq[1] = cell.f[1] - feq[1];
            fneq[5] = cell.f[5] - feq[5];
            fneq[8] = cell.f[8] - feq[8];
            //未知的分布函数均采用非平衡反弹得到
            fneq[3] = fneq[1];
            fneq[6] = fneq[8];
            fneq[7] = fneq[5];
            for (int k = 0; k < 9; ++k) {
                static double temp;
                temp = 0;
                for (int j = 0; j < 9; ++j) {
                    temp += fneq[j] * cc[k][j];
                }
                cell.f[k] = feq[k] + 4.5 * (1.0 - omega) * Cell::w[k] * temp;
            }
        }
    private:
        double omega;
        double rho;
    };
    /**************************************** 派生类RegularizedPressureBoundaryEdgeYP定义 ****************************************/
    class RegularizedPressureBoundaryEdgeYP :public RegularizedBoundary {
    public:
        RegularizedPressureBoundaryEdgeYP(double nu, double pressure)
            :omega(2.0 / (6.0 * nu + 1.0)) {
            rho = 3.0 * pressure;
        }
        //后面三个形参的作用是可以使用 第1种 executeLocalOperator 
        void operator()(Cell& cell, const double& a, const double& b, const double& c) {
            //根据方向修改
            double u = cell.f[1] + cell.f[5] - cell.f[3] - cell.f[6];
            double v = (cell.f[0] + cell.f[1] + cell.f[3] + 2.0 * (cell.f[2] + cell.f[5] + cell.f[6])) / rho - 1.0;
            static double uu;
            uu = u * u + v * v;
            static double feq[9];
            for (int k = 0; k < 9; ++k) {
                static double cu;
                cu = u * Cell::c[k][0] + v * Cell::c[k][1];
                feq[k] = rho * Cell::w[k] * (1 + 3.0 * cu + 4.5 * cu * cu - 1.5 * uu);
            }
            //根据方向修改
            static double fneq[9];
            fneq[0] = cell.f[0] - feq[0];
            fneq[1] = cell.f[1] - feq[1];
            fneq[2] = cell.f[2] - feq[2];
            fneq[3] = cell.f[3] - feq[3];
            fneq[5] = cell.f[5] - feq[5];
            fneq[6] = cell.f[6] - feq[6];
            //未知的分布函数均采用非平衡反弹得到
            fneq[4] = fneq[4];
            fneq[7] = fneq[7];
            fneq[8] = fneq[8];
            for (int k = 0; k < 9; ++k) {
                static double temp;
                temp = 0;
                for (int j = 0; j < 9; ++j) {
                    temp += fneq[j] * cc[k][j];
                }
                cell.f[k] = feq[k] + 4.5 * (1.0 - omega) * Cell::w[k] * temp;
            }
        }
    private:
        double omega;
        double rho;
    };
    /**************************************** 派生类RegularizedPressureBoundaryEdgeYN定义 ****************************************/
    class RegularizedPressureBoundaryEdgeYN :public RegularizedBoundary {
    public:
        RegularizedPressureBoundaryEdgeYN(double nu, double pressure)
            :omega(2.0 / (6.0 * nu + 1.0)) {
            rho = 3.0 * pressure;
        }
        //后面三个形参的作用是可以使用 第1种 executeLocalOperator 
        void operator()(Cell& cell, const double& a, const double& b, const double& c) {
            //根据方向修改
            double u = cell.f[1] + cell.f[8] - cell.f[3] - cell.f[7];
            double v = 1.0 - (cell.f[0] + cell.f[1] + cell.f[3] + 2.0 * (cell.f[4] + cell.f[7] + cell.f[8])) / rho;
            static double uu;
            uu = u * u + v * v;
            static double feq[9];
            for (int k = 0; k < 9; ++k) {
                static double cu;
                cu = u * Cell::c[k][0] + v * Cell::c[k][1];
                feq[k] = rho * Cell::w[k] * (1 + 3.0 * cu + 4.5 * cu * cu - 1.5 * uu);
            }
            //根据方向修改
            static double fneq[9];
            fneq[0] = cell.f[0] - feq[0];
            fneq[1] = cell.f[1] - feq[1];
            fneq[2] = cell.f[2] - feq[2];
            fneq[4] = cell.f[4] - feq[4];
            fneq[7] = cell.f[7] - feq[7];
            fneq[8] = cell.f[8] - feq[8];
            //未知的分布函数均采用非平衡反弹得到
            fneq[2] = fneq[2];
            fneq[5] = fneq[5];
            fneq[6] = fneq[6];
            for (int k = 0; k < 9; ++k) {
                static double temp;
                temp = 0;
                for (int j = 0; j < 9; ++j) {
                    temp += fneq[j] * cc[k][j];
                }
                cell.f[k] = feq[k] + 4.5 * (1.0 - omega) * Cell::w[k] * temp;
            }
        }
    private:
        double omega;
        double rho;
    };
    /**************************************** 派生类RegularizedPressureBoundaryCornerPP定义 ****************************************/
    class RegularizedPressureBoundaryCornerPP :public RegularizedBoundary {
    public:
        RegularizedPressureBoundaryCornerPP(double nu, double pressure)
            :omega(2.0 / (6.0 * nu + 1.0)) {
            rho = 3.0 * pressure;
        }
        //后面三个形参的作用是可以使用 第1种 executeLocalOperator 
        void operator()(Cell& cell, const double& a, const double& b, const double& c) {
            double u = 0;   //在边界角处的速度为0，这算是正常的设置方案
            double v = 0;
            static double uu;
            uu = u * u + v * v;
            static double feq[9];
            for (int k = 0; k < 9; ++k) {
                static double cu;
                cu = u * Cell::c[k][0] + v * Cell::c[k][1];
                feq[k] = rho * Cell::w[k] * (1 + 3.0 * cu + 4.5 * cu * cu - 1.5 * uu);
            }
            //根据方向修改
            static double fneq[9];
            fneq[0] = cell.f[0] - feq[0];
            fneq[1] = cell.f[1] - feq[1];
            fneq[2] = cell.f[2] - feq[2];
            fneq[5] = cell.f[5] - feq[5];
            fneq[6] = cell.f[6] - feq[6];
            fneq[8] = cell.f[8] - feq[8];
            //未知的分布函数均采用非平衡反弹得到
            fneq[3] = fneq[3];
            fneq[4] = fneq[4];
            fneq[7] = fneq[7];
            for (int k = 0; k < 9; ++k) {
                static double temp;
                temp = 0;
                for (int j = 0; j < 9; ++j) {
                    temp += fneq[j] * cc[k][j];
                }
                cell.f[k] = feq[k] + 4.5 * (1.0 - omega) * Cell::w[k] * temp;
            }
        }
    private:
        double omega;
        double rho;
    };
    /**************************************** 派生类RegularizedPressureBoundaryCornerPN定义 ****************************************/
    class RegularizedPressureBoundaryCornerPN :public RegularizedBoundary {
    public:
        RegularizedPressureBoundaryCornerPN(double nu, double pressure)
            :omega(2.0 / (6.0 * nu + 1.0)) {
            rho = 3.0 * pressure;
        }
        //后面三个形参的作用是可以使用 第1种 executeLocalOperator 
        void operator()(Cell& cell, const double& a, const double& b, const double& c) {
            double u = 0;   //在边界角处的速度为0，这算是正常的设置方案
            double v = 0;
            static double uu;
            uu = u * u + v * v;
            static double feq[9];
            for (int k = 0; k < 9; ++k) {
                static double cu;
                cu = u * Cell::c[k][0] + v * Cell::c[k][1];
                feq[k] = rho * Cell::w[k] * (1 + 3.0 * cu + 4.5 * cu * cu - 1.5 * uu);
            }
            //根据方向修改
            static double fneq[9];
            fneq[0] = cell.f[0] - feq[0];
            fneq[1] = cell.f[1] - feq[1];
            fneq[4] = cell.f[4] - feq[4];
            fneq[5] = cell.f[5] - feq[5];
            fneq[7] = cell.f[7] - feq[7];
            fneq[8] = cell.f[8] - feq[8];
            //未知的分布函数均采用非平衡反弹得到
            fneq[2] = fneq[2];
            fneq[3] = fneq[3];
            fneq[6] = fneq[6];
            for (int k = 0; k < 9; ++k) {
                static double temp;
                temp = 0;
                for (int j = 0; j < 9; ++j) {
                    temp += fneq[j] * cc[k][j];
                }
                cell.f[k] = feq[k] + 4.5 * (1.0 - omega) * Cell::w[k] * temp;
            }
        }
        void operator()(Cell& cell, const double& a, const double& b, const double& c, const double& d, const double& e) {
            //用于ForceBGK时Block的成员函数executeLocalOperator,实际与上面的一模一样
            this->operator()(cell, a, b, c);
        }
    private:
        double omega;
        double rho;
    };
    /**************************************** 派生类RegularizedPressureBoundaryCornerNP定义 ****************************************/
    class RegularizedPressureBoundaryCornerNP :public RegularizedBoundary {
    public:
        RegularizedPressureBoundaryCornerNP(double nu, double pressure)
            :omega(2.0 / (6.0 * nu + 1.0)) {
            rho = 3.0 * pressure;
        }
        //后面三个形参的作用是可以使用 第1种 executeLocalOperator 
        void operator()(Cell& cell, const double& a, const double& b, const double& c) {
            double u = 0;   //在边界角处的速度为0，这算是正常的设置方案
            double v = 0;
            static double uu;
            uu = u * u + v * v;
            static double feq[9];
            for (int k = 0; k < 9; ++k) {
                static double cu;
                cu = u * Cell::c[k][0] + v * Cell::c[k][1];
                feq[k] = rho * Cell::w[k] * (1 + 3.0 * cu + 4.5 * cu * cu - 1.5 * uu);
            }
            //根据方向修改
            static double fneq[9];
            fneq[0] = cell.f[0] - feq[0];
            fneq[2] = cell.f[2] - feq[2];
            fneq[3] = cell.f[3] - feq[3];
            fneq[5] = cell.f[5] - feq[5];
            fneq[6] = cell.f[6] - feq[6];
            fneq[7] = cell.f[7] - feq[7];
            //未知的分布函数均采用非平衡反弹得到
            fneq[1] = fneq[1];
            fneq[4] = fneq[4];
            fneq[8] = fneq[8];
            for (int k = 0; k < 9; ++k) {
                static double temp;
                temp = 0;
                for (int j = 0; j < 9; ++j) {
                    temp += fneq[j] * cc[k][j];
                }
                cell.f[k] = feq[k] + 4.5 * (1.0 - omega) * Cell::w[k] * temp;
            }
        }
    private:
        double omega;
        double rho;
    };
    /**************************************** 派生类RegularizedPressureBoundaryCornerNN定义 ****************************************/
    class RegularizedPressureBoundaryCornerNN :public RegularizedBoundary {
    public:
        RegularizedPressureBoundaryCornerNN(double nu, double pressure)
            :omega(2.0 / (6.0 * nu + 1.0)) {
            rho = 3.0 * pressure;
        }
        //后面三个形参的作用是可以使用 第1种 executeLocalOperator 
        void operator()(Cell& cell, const double& a, const double& b, const double& c) {
            double u = 0;
            double v = 0;
            static double uu;
            uu = u * u + v * v;
            static double feq[9];
            for (int k = 0; k < 9; ++k) {
                static double cu;
                cu = u * Cell::c[k][0] + v * Cell::c[k][1];
                feq[k] = rho * Cell::w[k] * (1 + 3.0 * cu + 4.5 * cu * cu - 1.5 * uu);
            }
            //根据方向修改
            static double fneq[9];
            fneq[0] = cell.f[0] - feq[0];
            fneq[3] = cell.f[3] - feq[3];
            fneq[4] = cell.f[4] - feq[4];
            fneq[6] = cell.f[6] - feq[6];
            fneq[7] = cell.f[7] - feq[7];
            fneq[8] = cell.f[8] - feq[8];
            //未知的分布函数均采用非平衡反弹得到
            fneq[1] = fneq[1];
            fneq[2] = fneq[2];
            fneq[5] = fneq[5];
            for (int k = 0; k < 9; ++k) {
                static double temp;
                temp = 0;
                for (int j = 0; j < 9; ++j) {
                    temp += fneq[j] * cc[k][j];
                }
                cell.f[k] = feq[k] + 4.5 * (1.0 - omega) * Cell::w[k] * temp;
            }
        }
    private:
        double omega;
        double rho;
    };

}//namespace

#endif
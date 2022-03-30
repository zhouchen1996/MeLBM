#ifndef _UTILITIES_H_
#define _UTILITIES_H_

#include "cell.h"
#include <cmath>

namespace LBM {

    /**************************************** 定义局部操作 ****************************************/
	struct SetLocalOperator {

        LocalOperator* localOperator;	//数据成员为一个指向LocalOperator的指针,用于多态处理不同局部操作

        int level; // 0 表示 localOperators[0]; 1 表示 localOperators[1]

		SetLocalOperator(LocalOperator* localOperator_, int _level) 
            :localOperator(localOperator_), level(_level) {}
	    
        SetLocalOperator& changeTo(LocalOperator* localOperatorNew, int levelNew) {
			localOperator = localOperatorNew;
            level = levelNew;
			return *this;
		}

		void operator()(Cell& cell) {
			cell.localOperators[level] = localOperator;
		}

	};

	/**************************************** cosphi ****************************************/
	double cosphi(int& k, double& Fc, double& FF) {
		
        //这一步对于质量的守恒至关重要,一定要是整个流体区域内做重新着色
		
        if(k > 0){
			
            double temp = Fc / sqrt(FF);
			
            return ((std::isnormal(temp) == 1) ? temp : 0);
		}
		
        return 0;
	}
	/**************************************** Omega ****************************************/
	struct Omega {	//用于构造多相流中会随相分布动态变化的弛豫率
		
        double invnu_r, invnu_b;	//使用粘滞系数的倒数，减少除法的使用
        
        Omega(double nu_r_, double  nu_b_) : invnu_r(1.0 / nu_r_), invnu_b(1.0 / nu_b_) {}
		
        const double& operator()(const double& Phi) {
			
            //使用的是刘海湖的方法，我通过使用倒数，将除法减少到只使用一次
			
            static double invnu_;
			invnu_ = 0.5 * ((1.0 + Phi) * invnu_r + (1.0 - Phi) * invnu_b);
			
            static double omegaEff;
			omegaEff = (2.0 * invnu_) / (6.0 + invnu_);
			
            return omegaEff;
		}
		
	};

	/**************************************** 初始化器 ****************************************/
	template <typename T = double>
	struct ScalarInitializer {	//用于初始化标量场

        T rho;

		ScalarInitializer(T rho_ = 0) :rho(rho_) {}
	    
        void operator()(T& a) {	//形参为引用类型
			a = rho;
		}
		
        ScalarInitializer<T>& changeTo(T rhoNew) {
			
            rho = rhoNew;
			
            return *this;
		}

	};

	template <typename T = double>
	struct VectorInitializer {	//用于初始化向量场

		T ux;
		T uy;

		VectorInitializer(T ux_ = 0, T uy_ = 0) :ux(ux_), uy(uy_) {}
	    
        void operator()(T* a) {	//形参为指针
			
            a[0] = ux;
			
            a[1] = uy;

		}

		VectorInitializer<T>& changeTo(T uxNew, T uyNew) {
			
            ux = uxNew;
			
            uy = uyNew;
			
            return *this;
		}

	};

	struct DistributionFunctionInitializer {	//用于初始化分布函数

        double f[9];

		DistributionFunctionInitializer(const double(&f_)[9]) {
		    
            for (int i = 0; i < 9; ++i)
				f[i] = f_[i];
		
        }
		
        void operator()(LBM::Cell& cell) {	//形参为Cell的引用类型
			
            for (int i = 0; i < 9; ++i)
				cell.f[i] = f[i];
		
        }

		
	};
	
}//namespace

#endif // !_UTILITIES_H_

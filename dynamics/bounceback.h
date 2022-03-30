#ifndef _BOUNCEBACK_H_
#define _BOUNCEBACK_H_

#include "../core/cell.h"
#include "../core/localoperator.h"

namespace LBM {

	/**************************************** 局部算子BounceBack定义 ****************************************/
	
	//固体晶格在碰撞步骤执行的是反弹
	struct BounceBack :public LocalOperator {

		//full way bounce back

		//交换a与b的值
		void swp(double& a, double& b) {
			double tmp = a;
			a = b;
			b = tmp;
		}

		void BB(Cell& cell) {
			swp(cell.f[1], cell.f[3]);
			swp(cell.f[2], cell.f[4]);
			swp(cell.f[5], cell.f[7]);
			swp(cell.f[6], cell.f[8]);
		}

        virtual void operator()(Cell& cell, const double& a1) {
			BB(cell);
		}

		virtual void operator()(Cell& cell, const double& a1, const double& a2) {
			BB(cell);
		}

		virtual void operator()(Cell& cell, const double& a1, const double& a2, const double& a3) {
			BB(cell);
		}

		virtual void operator()(Cell& cell, const double& a1, const double& a2, const double& a3, const double& a4) {
			BB(cell);
		}

		virtual void operator()(Cell& cell, const double& a1, const double& a2, const double& a3, const double& a4, 
        	const double& a5) {
			BB(cell);
		}

		virtual void operator()(Cell& cell, const double& a1, const double& a2, const double& a3, const double& a4, 
        	const double& a5, const double& a6) {
			BB(cell);
		}

		virtual void operator()(Cell& cell, const double& a1, const double& a2, const double& a3, const double& a4, 
        	const double& a5, const double& a6, const double &a7) {
			BB(cell);
		}

        virtual void operator()(Cell& cell, const double& a1, const double& a2, const double& a3, const double& a4, 
			const double& a5, const double& a6, const double &a7, const double& a8) {
			BB(cell);
		}

	};

} // namespace LBM

#endif
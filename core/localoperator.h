#ifndef _LOCALOPERATOR_H_
#define _LOCALOPERATOR_H_

namespace LBM {

    struct Cell;

    //局部算子，用于局部地处理分布函数，例如: 碰撞 或 固体晶格的反弹
	struct LocalOperator {	
		
        virtual void operator()(Cell& cell, const double& a1) {}

		virtual void operator()(Cell& cell, const double& a1, const double& a2) {}

		virtual void operator()(Cell& cell, const double& a1, const double& a2, const double& a3) {}

		virtual void operator()(Cell& cell, const double& a1, const double& a2, const double& a3, const double& a4) {}

		virtual void operator()(Cell& cell, const double& a1, const double& a2, const double& a3, const double& a4, 
            const double& a5) {}

		virtual void operator()(Cell& cell, const double& a1, const double& a2, const double& a3, const double& a4, 
            const double& a5, const double& a6) {}

		virtual void operator()(Cell& cell, const double& a1, const double& a2, const double& a3, const double& a4, 
            const double& a5, const double& a6, const double &a7) {}

        virtual void operator()(Cell& cell, const double& a1, const double& a2, const double& a3, const double& a4, 
            const double& a5, const double& a6, const double &a7, const double& a8) {}

    };

}//namespace

#endif
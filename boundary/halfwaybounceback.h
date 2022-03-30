#ifndef _HALFWAYBOUNCEBACK_H_
#define _HALFWAYBOUNCEBACK_H_

#include "../core/cell.h"
#include "../core/localoperator.h"

namespace LBM {

	//不建议使用

	/**************************************** 半步回弹边界HalfWayBounceBack定义 ****************************************/
	
	//用于流之后,半步回弹边界需要相邻晶格的分布函数,即将流出去的分布函数补给与之方向相反的自己的未知分布函数
	struct HalfWayBounceBack :public LocalOperator {
        
		//half way bounce back

	};

} // namespace LBM

#endif
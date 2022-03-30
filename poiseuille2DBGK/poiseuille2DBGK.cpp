#include "../domain/block2D.h"
#include "../domain/scalarfield2D.h"
#include "../domain/vectorfield2D.h"
#include "../dynamics/bgk.h"
#include "../dynamics/bounceback.h"
#include "../boundary/regularizedboundary.h"
#include "../boundary/outflow.h"
#include "../io/vtkWriter.h"

#include <ctime>

int main() {
    using namespace LBM;

    double nu = 0.1; //粘滞系数
    double ux = 0.001; //流速
    double uy = 0;
    double rho0 = 1; //初始密度
    double velx0 = 0; //初始速度
    double vely0 = 0;
    int X = 200;
    int Y = 50;
    int nstep = 1000; //运行时间
    int interval = 10; //文件输出间隔

    std::cout << "X Y =";
    std::cin >> X >> Y;

    std::cout <<"Viscosity =";
    std::cin >> nu;

    std::cout << "Velocity =";
    std::cin >> ux;

    std::cout << "All timestep =";
    std::cin >> nstep;

    std::cout << "Interval of outfilevtk =";
    std::cin >> interval;

    //定义 X*Y的计算域
    Block2D block(X, Y);
    ScalarField2D<double> densityfield(X, Y);
    VectorField2D<double> velocityfield(X, Y);

    //定义输出文件
    VtkWriter outvtk;
    outvtk.add(densityfield, "density");
    outvtk.add(velocityfield, "velocity");

    //定义碰撞算子
    BGK* bgk = new BGK(nu);
    BounceBack* bb = new BounceBack();
    RegularizedVelocityBoundaryEdgeXN* velboundaryXN = new RegularizedVelocityBoundaryEdgeXN(nu, ux, uy);
    OutFlowEdgeXP* outflow = new OutFlowEdgeXP();
    
    //定义一个算子设置器
    SetLocalOperator setOp(bb, 0);

    //定义算子
    block.traversal(setOp); //先全部定义为bb
    block.traversal(setOp.changeTo(bgk, 0), 1, X, 1, Y); //内部流体晶格定义为bgk
    block.traversal(setOp.changeTo(velboundaryXN, 0), 0, 0, 1, Y); //左边界晶格定义一个速度
    block.traversal(setOp.changeTo(outflow, 0), X + 1, X + 1, 1, Y);

    //初始化密度\速度\分布函数
    ScalarInitializer<double> inidensity(rho0);
    VectorInitializer<double> inivelocity(velx0, vely0);
    DistributionFunctionInitializer inidisfunction(Cell::w);
    densityfield.traversal(inidensity);
    velocityfield.traversal(inivelocity);
    block.traversal(inidisfunction);

    //开始计算
    clock_t start, stop;
	start = clock();

	/**************************************** 开始运算 ****************************************/
	for (int t = 1; t <= nstep; ++t) {

		block.executeLocalOperator(densityfield, velocityfield);
		
        block.streaming();

		block.computeDensityField(densityfield); //先计算密度再计算速度
		block.computeVelocityField(velocityfield, densityfield);

		if (t % interval == 0) {
            std::cout << "step = " << t << std::endl;
			//std::cout << densityfield.ps[X / 2][Y / 2] << std::endl;
			outvtk.ascii("./outfile/poiseuille2dbgk" + std::to_string(t) + std::string(".vti"), 1, X, 1, Y);
		}
	}
	stop = clock();
	std::cout << "Time: " << ((double)stop - (double)start) / CLK_TCK << std::endl;

    system("pause");

    return 0;
}
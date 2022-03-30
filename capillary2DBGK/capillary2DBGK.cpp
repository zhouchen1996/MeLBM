#include "../domain/block2D.h"
#include "../domain/scalarfield2D.h"
#include "../domain/vectorfield2D.h"
#include "../dynamics/bgk.h"
#include "../dynamics/bgkITF.h"
#include "../dynamics/bounceback.h"
#include "../boundary/regularizedboundary.h"
#include "../boundary/outflow.h"
#include "../io/vtkWriter.h"

#include <ctime>

int main() {

    using namespace LBM;

    double nu_r = 0.1; //粘滞系数
    double nu_b = 0.1;

    double rho0 = 1;  //初始密度
    double rho_r0 = 1; 
    double rho_b0 = 1;

    double theta = 45.0 / 180.0 * 3.1415926; //接触角
    double sigma = 0.1; //表面张力系数

    int X = 100;
    int Y = 30;
    int nstep = 200000; //运行时间
    int interval = 500; //文件输出间隔

    //定义 X*Y的计算域
    Block2D block(X, Y);
    ScalarField2D<double> densityfield(X, Y);
    VectorField2D<double> velocityfield(X, Y);

    Block2D redblock(X, Y);
    ScalarField2D<double> reddensityfield(X, Y);
    
    Block2D blueblock(X, Y);
    ScalarField2D<double> bluedensityfield(X, Y);

    ScalarField2D<double> phasefield(X, Y);
    VectorField2D<double> phasefieldgradient(X, Y);

    //定义输出文件
    VtkWriter outvtk;
    outvtk.add(densityfield, "density");
    outvtk.add(phasefield, "phasefield");
    outvtk.add(velocityfield, "velocity");

    //定义碰撞算子
    BGKITF* bgkITF = new BGKITF(nu_r, nu_b, sigma);
    BounceBack* bb = new BounceBack();

    OutFlowEdgeXP* openboundaryleft = new OutFlowEdgeXP();
    OutFlowEdgeXN* openboundaryright = new OutFlowEdgeXN();
    
    //定义一个算子设置器
    SetLocalOperator setOp(bb, 0);

    //定义算子
    block.traversal(setOp); //先全部定义为bb
    block.traversal(setOp.changeTo(bgkITF, 0), 1, X, 1, Y); //内部流体晶格定义为bgkITF
    block.traversal(setOp.changeTo(openboundaryleft, 0), 0, 0, 1, Y); //左边界为开口
    block.traversal(setOp.changeTo(openboundaryright, 0), X + 1, X + 1, 1, Y);//右边界为开口

    redblock.traversal(setOp.changeTo(bb, 0));
    blueblock.traversal(setOp.changeTo(bb, 0));

    //初始化密度\速度\分布函数
    ScalarInitializer<double> inidensity(rho0);
    ScalarInitializer<double> iniphasefield(cos(theta));
    VectorInitializer<double> inivelocity(0, 0);
    VectorInitializer<double> iniphasefieldgradient(0, 0);
    DistributionFunctionInitializer inidisfunction(Cell::w);

    block.traversal(inidisfunction);
    densityfield.traversal(inidensity);
    velocityfield.traversal(inivelocity);
    
    redblock.traversal(inidisfunction);
    reddensityfield.traversal(inidensity.changeTo(0));
    reddensityfield.traversal(inidensity.changeTo(rho_r0), 0, X / 2, 1, Y);

    blueblock.traversal(inidisfunction);
    bluedensityfield.traversal(inidensity.changeTo(rho_b0));
    bluedensityfield.minus(bluedensityfield, reddensityfield);

    phasefield.traversal(iniphasefield);
    phasefield.traversal(iniphasefield.changeTo(0), 0, X + 1, 1, Y);

    phasefieldgradient.traversal(iniphasefieldgradient);

    //开始计算
    clock_t start, stop;
	start = clock();

	/**************************************** 开始运算 ****************************************/
	for (int t = 0; t <= nstep; ++t) {

        phasefield.computePhaseField(reddensityfield, bluedensityfield, densityfield, 1, X, 1, Y);
        
        phasefieldgradient.gradient(phasefield, 1, X, 1, Y);
		
        redblock.executeLocalOperator(reddensityfield, velocityfield);
        blueblock.executeLocalOperator(bluedensityfield, velocityfield);
        block.executeLocalOperator(densityfield, velocityfield, phasefield, phasefieldgradient);

        block.recolor(redblock, blueblock, reddensityfield, bluedensityfield, 
            densityfield, phasefieldgradient, 1, X, 1, Y);

        redblock.streaming();
        blueblock.streaming();

		redblock.computeDensityField(reddensityfield);
        blueblock.computeDensityField(bluedensityfield);
        densityfield.plus(reddensityfield, bluedensityfield);
		block.computeVelocityField(velocityfield, densityfield);

		if (t % interval == 0) {
            std::cout << "step = " << t << std::endl;
			outvtk.ascii("./outfile/capillary2dbgk" + std::to_string(t) + std::string(".vti"), 1, X, 1, Y);
		}
	}
	stop = clock();
	std::cout << "Time: " << ((double)stop - (double)start) / CLK_TCK << std::endl;

    system("pause");

    return 0;
}
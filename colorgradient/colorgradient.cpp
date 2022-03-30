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

//经检验满足拉普拉斯定律

int main() {

    using namespace LBM;

    double nu_r = 0.1; //粘滞系数
    double nu_b = 0.1;

    double rho0 = 1;  //初始密度
    double rho_r0 = 1; 
    double rho_b0 = 1;

    double theta = 0.001; //表面张力系数

    int X = 50;
    int Y = 50;
    int nstep = 50000; //运行时间
    int interval = 1000; //文件输出间隔

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
    BGKITF* bgkITF = new BGKITF(nu_r, nu_b, theta);
    BounceBack* bb = new BounceBack();
    
    //定义一个算子设置器
    SetLocalOperator setOp(bb, 0);

    //定义算子
    block.traversal(setOp); //先全部定义为bb
    block.traversal(setOp.changeTo(bgkITF, 0), 1, X, 1, Y); //内部流体晶格定义为bgkITF

    redblock.traversal(setOp.changeTo(bb, 0));
    blueblock.traversal(setOp.changeTo(bb, 0));

    //初始化密度\速度\分布函数
    ScalarInitializer<double> inidensity(rho0);
    ScalarInitializer<double> iniphasefield(0);
    VectorInitializer<double> inivelocity(0, 0);
    VectorInitializer<double> iniphasefieldgradient(0, 0);
    DistributionFunctionInitializer inidisfunction(Cell::w);

    block.traversal(inidisfunction);
    densityfield.traversal(inidensity);
    velocityfield.traversal(inivelocity);
    
    redblock.traversal(inidisfunction);
    reddensityfield.traversal(inidensity.changeTo(0));
    reddensityfield.traversal(inidensity.changeTo(rho_r0), 1 * X / 4, 3 * X / 4, 1 * Y / 4, 3 * Y / 4);

    blueblock.traversal(inidisfunction);
    bluedensityfield.traversal(inidensity.changeTo(rho_b0));
    bluedensityfield.minus(bluedensityfield, reddensityfield);

    phasefield.traversal(iniphasefield);
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
			outvtk.ascii("./outfile/cgm2dbgk" + std::to_string(t) + std::string(".vti"), 1, X, 1, Y);
		}
	}
	stop = clock();
	std::cout << "Time: " << ((double)stop - (double)start) / CLK_TCK << std::endl;

    system("pause");

    return 0;
}
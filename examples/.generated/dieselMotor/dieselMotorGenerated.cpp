
// CODEGEN FOR MODEL "dieselMotor"

// includes
#define _USE_MATH_DEFINES
#include "dieselMotorGeneratedParams.h"
#include <cmath>
#include <string>
#include "constants.h"
#include <problem.h>
#include "integrator.h"
#include "mesh.h"
#include "gdop.h"
#include "solver.h"


// runtime parameters and global constants


// mayer term
class MayerdieselMotor : public Expression {
public:
	static std::unique_ptr<MayerdieselMotor> create() {
		Adjacency adj{{0, 1, 2, 3}, {}, {}};
		AdjacencyDiff adjDiff{{{0, 0}, {1, 1}, {2, 2}, {3, 3}}, {}, {}, {}, {}, {}};
		return std::unique_ptr<MayerdieselMotor>(new MayerdieselMotor(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return pow(x[0] - 0.51530917068559601, 2) + pow(x[1] - 0.54705585422599101, 2) + pow(x[2] - 0.381048005791294, 2) + pow(x[3] - 0.27144300053768, 2);
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{2*x[0] - 1.030618341371192, 2*x[1] - 1.094111708451982, 2*x[2] - 0.762096011582588, 2*x[3] - 0.54288600107536}, {}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{2, 2, 2, 2}, {}, {}, {}, {}, {}};
	}
private:
	MayerdieselMotor(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


// lagrange term
class LagrangedieselMotor : public Expression {
public:
	static std::unique_ptr<LagrangedieselMotor> create() {
		Adjacency adj{{0}, {0}, {}};
		AdjacencyDiff adjDiff{{}, {{0, 0}}, {}, {}, {}, {}};
		return std::unique_ptr<LagrangedieselMotor>(new LagrangedieselMotor(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return 0.08249999999999999*u[0]*x[0]/M_PI;
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        const double x0 = 0.08249999999999999/M_PI;
		return {std::vector<double>{u[0]*x0}, {x0*x[0]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        const double x0 = 0.08249999999999999/M_PI;
		return {std::vector<double>{}, {x0}, {}, {}, {}, {}};
	}
private:
	LagrangedieselMotor(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


// dynamic constraints
class F0dieselMotor : public Expression {
public:
	static std::unique_ptr<F0dieselMotor> create() {
		Adjacency adj{{0, 1, 2}, {0}, {}};
		AdjacencyDiff adjDiff{{{0, 0}}, {}, {}, {}, {}, {}};
		return std::unique_ptr<F0dieselMotor>(new F0dieselMotor(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
        const double x0 = M_1_PI;
		return 0.000324675324675325*x0*(27936.417547812125*u[0] + 1185.5800082209614*x0*x[0] - 39807.940708436552*pow(x[0], 2)/pow(M_PI, 2) + 2540.0*x[1] - 3810.0*x[2] - 455.98330550589981);
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        const double x0 = M_1_PI;
		return {std::vector<double>{0.000324675324675325*x0*(1185.5800082209614*x0 - 79615.881416873104*x[0]/pow(M_PI, 2)), 0.82467532467532545*x0, -1.2370129870129882*x0}, {9.0702654376013481*x0}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-25.84931214833545/pow(M_PI, 3)}, {}, {}, {}, {}, {}};
	}
private:
	F0dieselMotor(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F1dieselMotor : public Expression {
public:
	static std::unique_ptr<F1dieselMotor> create() {
		Adjacency adj{{0, 1, 3}, {}, {}};
		AdjacencyDiff adjDiff{{{1, 0}, {1, 1}, {3, 1}, {3, 3}}, {}, {}, {}, {}, {}};
		return std::unique_ptr<F1dieselMotor>(new F1dieselMotor(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return -29.202218227585707*x[0]*x[1]/M_PI + 21.105690996021139*sqrt(-pow(x[1], 2)*pow(0.38110788300674586*pow(x[3], 2) + 1, -7.0452961672473862) + 0.25558758431328077);
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        const double x0 = M_1_PI;
        const double x1 = pow(x[1], 2);
        const double x2 = 0.38110788300674586*pow(x[3], 2) + 1;
        const double x3 = pow(x2, -7.0452961672473862);
        const double x4 = 0.50555670731707314/sqrt(-x1*x3 + 0.25558758431328077);
		return {std::vector<double>{-29.202218227585707*x0*x[1], -29.202218227585707*x0*x[0] - 41.747425542084933*x3*x4*x[1], 112.09258517065034*x1*pow(x2, -8.0452961672473862)*x4*x[3]}, {}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        const double x0 = pow(x[1], 2);
        const double x1 = pow(x[3], 2);
        const double x2 = 0.38110788300674586*x1 + 1;
        const double x3 = pow(x2, -7.0452961672473862);
        const double x4 = -3.9125531182855591*x0*x3 + 1;
        const double x5 = pow(x4, -1.0/2.0);
        const double x6 = pow(x4, -3.0/2.0);
        const double x7 = pow(x2, -8.0452961672473862)*x5;
		return {std::vector<double>{-29.202218227585707/M_PI, -163.33901998507861*x0*pow(x2, -14.090592334494772)*x6 - 41.747425542084933*x3*x5, 438.56819364611761*pow(x2, -15.090592334494772)*x6*pow(x[1], 3)*x[3] + 224.18517034130068*x7*x[1]*x[3], -687.37993262256464*x0*x1*pow(x2, -9.0452961672473862)*x5 + 112.09258517065034*x0*x7 - 1177.5634535801021*x1*pow(x2, -16.090592334494772)*x6*pow(x[1], 4)}, {}, {}, {}, {}, {}};
	}
private:
	F1dieselMotor(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F2dieselMotor : public Expression {
public:
	static std::unique_ptr<F2dieselMotor> create() {
		Adjacency adj{{0, 1, 2}, {0, 1}, {}};
		AdjacencyDiff adjDiff{{{0, 0}, {1, 0}, {1, 1}, {2, 0}, {2, 1}, {2, 2}}, {{0, 0}, {0, 1}, {0, 2}, {1, 0}, {1, 1}, {1, 2}}, {{0, 0}, {1, 0}}, {}, {}, {}};
		return std::unique_ptr<F2dieselMotor>(new F2dieselMotor(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
        const double x0 = 1.0/x[2];
        const double x1 = x[0]/M_PI;
        const double x2 = u[0]*x1;
        const double x3 = 1.4420483946140072*x1*x[1] + 0.08249999999999999*x2;
        const double x4 = x2/x3;
        const double x5 = 0.72900098370469424*pow(x[2]/x[1], 0.28387734915924823)*pow(x4 + 0.32842745509493809, -0.28387734915924823)*(4305.1556436685933*x4 + 930.65395757192096);
        const double x6 = x[2]/sqrt(x5);
        const double x7 = sqrt(x0);
		return 0.018144065864949988*x5*(-32.986003833080048*u[1]*x6*sqrt(pow(x0, 1.5705705705705706) - 0.79174558284615604*pow(x0, 1.7852852852852854)) + x3 - 57.149103182867073*x6*sqrt(pow(x7, 1.5705705705705706) - 0.88980086696190397*pow(x7, 1.7852852852852854)));
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        const double x0 = M_1_PI;
        const double x1 = 1.4420483946140072*x0;
        const double x2 = x1*x[0];
        const double x3 = 0.08249999999999999*x0;
        const double x4 = x3*x[0];
        const double x5 = u[0]*x4 + x2*x[1];
        const double x6 = x0/x5;
        const double x7 = 4305.1556436685933*x6;
        const double x8 = u[0]*x3 + x1*x[1];
        const double x9 = u[0]/pow(x5, 2);
        const double x10 = x0*x8*x9*x[0];
        const double x11 = 3.0448124372273666*x6;
        const double x12 = x11*x[0];
        const double x13 = u[0]*x12 + 1;
        const double x14 = pow(x13, -0.28387734915924823);
        const double x15 = 1.0/x[1];
        const double x16 = x15*x[2];
        const double x17 = pow(x16, 0.28387734915924823);
        const double x18 = x14*x17;
        const double x19 = x18*(u[0]*x7 - 4305.1556436685933*x10);
        const double x20 = u[0]*x11 - 3.0448124372273666*x10;
        const double x21 = x7*x[0];
        const double x22 = u[0]*x21 + 930.65395757192096;
        const double x23 = 0.28387734915924823*x22;
        const double x24 = pow(x13, -1.2838773491592481)*x17;
        const double x25 = x23*x24;
        const double x26 = x19 - x20*x25;
        const double x27 = 43.795847987340309*x[2];
        const double x28 = 1.0/x[2];
        const double x29 = sqrt(x28);
        const double x30 = 0.65244887140199537*sqrt(pow(x29, 1.5705705705705706) - 0.88980086696190397*pow(x29, 1.7852852852852854));
        const double x31 = x18*x22;
        const double x32 = pow(x31, -3.0/2.0);
        const double x33 = x30*x32;
        const double x34 = x27*x33;
        const double x35 = 38.744203843894169*x[2];
        const double x36 = 0.42568952979373742*sqrt(pow(x28, 1.5705705705705706) - 0.79174558284615604*pow(x28, 1.7852852852852854));
        const double x37 = u[1]*x36;
        const double x38 = x32*x37;
        const double x39 = x35*x38;
        const double x40 = 0.018144065864949988*x31;
        const double x41 = pow(x31, -1.0/2.0);
        const double x42 = 77.488407687788339*x37*x41;
        const double x43 = 87.591695974680619*x30*x41;
        const double x44 = -x42*x[2] - x43*x[2] + x5;
        const double x45 = 0.018144065864949988*x44;
        const double x46 = 0.0051506893207128044*x44;
        const double x47 = x22*x46;
        const double x48 = x24*x47;
        const double x49 = x14*pow(x16, -0.71612265084075177);
        const double x50 = x49*x[2]/pow(x[1], 2);
        const double x51 = x9*pow(x[0], 2)/pow(M_PI, 2);
        const double x52 = x18*x51;
        const double x53 = x22*x24*x51;
        const double x54 = -x23*x50 - 6208.2427845157281*x52 + 1.2464392646590381*x53;
        const double x55 = pow(x[2], -3.0/2.0);
        const double x56 = x22*x49;
        const double x57 = x16*x56;
        const double x58 = x18*(x21 - 355.17534060265888*x51);
        const double x59 = x12 - 0.25119702607125771*x51;
        const double x60 = -x25*x59 + x58;
		return {std::vector<double>{x19*x45 - x20*x48 + x40*(x26*x34 + x26*x39 + x8), x40*(x2 + x34*x54 + x39*x54) - 112.64276598785388*x44*x52 + 0.022615476114633416*x44*x53 - x47*x50, x15*x46*x56 + x40*(-u[1]*x35*x41*(0.25614106532433373*pow(x[2], -2.7852852852852852) - 0.28460556796052533*pow(x[2], -2.5705705705705704))/x36 - x27*x41*(-0.33428772384703404*pow(x29, 0.57057057057057059)*x55 + 0.33811420958044264*pow(x29, 0.78528528528528541)*x55)/x30 + 12.432649230827563*x33*x57 + 10.998601882490233*x38*x57 - x42 - x43)}, {x40*(x34*x60 + x39*x60 + x4) + x45*x58 - x48*x59, -1.4059547728573285*x31*x36*x41*x[2]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        const double x0 = M_1_PI;
        const double x1 = 0.08249999999999999*x0;
        const double x2 = u[0]*x1;
        const double x3 = 1.4420483946140072*x0;
        const double x4 = x3*x[1];
        const double x5 = x2 + x4;
        const double x6 = x2*x[0] + x4*x[0];
        const double x7 = pow(x6, -3);
        const double x8 = u[0]*x7;
        const double x9 = x0*pow(x5, 2)*x8*x[0];
        const double x10 = u[0]*x5;
        const double x11 = pow(x6, -2);
        const double x12 = x0*x11;
        const double x13 = x10*x12;
        const double x14 = -6.0896248744547332*x13 + 6.0896248744547332*x9;
        const double x15 = 1.0/x[2];
        const double x16 = 0.18121157577601327*pow(x15, 1.5705705705705706) - 0.14347346468125*pow(x15, 1.7852852852852854);
        const double x17 = sqrt(x16);
        const double x18 = x0/x6;
        const double x19 = 4305.1556436685933*x18;
        const double x20 = x19*x[0];
        const double x21 = u[0]*x20 + 930.65395757192096;
        const double x22 = 3.0448124372273666*x18;
        const double x23 = u[0]*x22;
        const double x24 = x23*x[0] + 1;
        const double x25 = pow(x24, -0.28387734915924823);
        const double x26 = 1.0/x[1];
        const double x27 = x26*x[2];
        const double x28 = pow(x27, 0.28387734915924823);
        const double x29 = x25*x28;
        const double x30 = x21*x29;
        const double x31 = pow(x30, -1.0/2.0);
        const double x32 = 77.488407687788339*x31;
        const double x33 = x17*x32;
        const double x34 = u[1]*x33;
        const double x35 = sqrt(x15);
        const double x36 = 0.42568952979373742*pow(x35, 1.5705705705705706) - 0.37877891266707281*pow(x35, 1.7852852852852854);
        const double x37 = sqrt(x36);
        const double x38 = 87.591695974680619*x31;
        const double x39 = x37*x38;
        const double x40 = -x34*x[2] - x39*x[2] + x6;
        const double x41 = pow(x24, -1.2838773491592481);
        const double x42 = x21*x28;
        const double x43 = 0.0051506893207128044*x42;
        const double x44 = x41*x43;
        const double x45 = x40*x44;
        const double x46 = 4305.1556436685933*x[0];
        const double x47 = u[0]*x19 - x13*x46;
        const double x48 = 3.0448124372273666*x[0];
        const double x49 = -x13*x48 + x23;
        const double x50 = x28*x41;
        const double x51 = x49*x50;
        const double x52 = x47*x51;
        const double x53 = 0.010301378641425609*x40;
        const double x54 = x29*(-8610.3112873371865*x13 + 8610.3112873371865*x9);
        const double x55 = 0.018144065864949988*x40;
        const double x56 = pow(x49, 2);
        const double x57 = pow(x24, -2.2838773491592481)*x42;
        const double x58 = x40*x57;
        const double x59 = 0.0066128533514196036*x58;
        const double x60 = x29*x47;
        const double x61 = 0.28387734915924823*x41;
        const double x62 = x42*x61;
        const double x63 = -x49*x62 + x60;
        const double x64 = pow(x30, -3.0/2.0);
        const double x65 = 43.795847987340309*x37;
        const double x66 = x64*x65;
        const double x67 = x63*x66;
        const double x68 = x63*x64;
        const double x69 = 38.744203843894169*u[1];
        const double x70 = x17*x69;
        const double x71 = x68*x70;
        const double x72 = x5 + x67*x[2] + x71*x[2];
        const double x73 = x41*x49;
        const double x74 = 0.010301378641425609*x42;
        const double x75 = 0.36446369852492994*x57;
        const double x76 = -x14*x62 - 0.56775469831849645*x52 + x54 + x56*x75;
        const double x77 = x64*x70;
        const double x78 = x77*x[2];
        const double x79 = pow(x63, 2);
        const double x80 = pow(x30, -5.0/2.0);
        const double x81 = x80*x[2];
        const double x82 = u[1]*x17;
        const double x83 = 58.116305765841254*x82;
        const double x84 = x81*x83;
        const double x85 = x66*x[2];
        const double x86 = 65.693771981010457*x37;
        const double x87 = x81*x86;
        const double x88 = 0.018144065864949988*x30;
        const double x89 = pow(x[0], 2);
        const double x90 = pow(M_PI, -2);
        const double x91 = x11*x90;
        const double x92 = x89*x91;
        const double x93 = u[0]*x92;
        const double x94 = x20 - 355.17534060265888*x93;
        const double x95 = 0.0051506893207128044*x40;
        const double x96 = x94*x95;
        const double x97 = x12*x5;
        const double x98 = x10*x7*x89*x90;
        const double x99 = u[0]*x91*x[0];
        const double x100 = x22 - x48*x97 + 0.50239405214251542*x98 - 0.50239405214251542*x99;
        const double x101 = x22*x[0] - 0.25119702607125771*x93;
        const double x102 = x101*x49;
        const double x103 = x101*x47;
        const double x104 = x29*x94;
        const double x105 = -x101*x62 + x104;
        const double x106 = x105*x66;
        const double x107 = x105*x77;
        const double x108 = x1*x[0] + x106*x[2] + x107*x[2];
        const double x109 = 0.018144065864949988*x60;
        const double x110 = x101*x44;
        const double x111 = x43*x73;
        const double x112 = x29*(x19 - x46*x97 + 710.35068120531776*x98 - 710.35068120531776*x99);
        const double x113 = 0.018144065864949988*x104;
        const double x114 = 0.28387734915924823*x94;
        const double x115 = -x100*x62 + x102*x75 - x103*x28*x61 + x112 - x114*x51;
        const double x116 = x105*x84;
        const double x117 = x105*x87;
        const double x118 = -x100*x45 + x102*x59 - x103*x50*x95 + x108*x109 - x108*x111 - x110*x72 + x112*x55 + x113*x72 - x51*x96 + x88*(x1 + x115*x78 + x115*x85 - x116*x63 - x117*x63);
        const double x119 = x31*x[2];
        const double x120 = 0.39911871395653137*x17;
        const double x121 = x119*x120*x42;
        const double x122 = 1.4059547728573285*x17;
        const double x123 = x119*x122;
        const double x124 = x68*x[2];
        const double x125 = 0.70297738642866425*x30;
        const double x126 = x125*x17;
        const double x127 = x121*x73 - x123*x60 + x124*x126;
        const double x128 = x41*x42;
        const double x129 = x128*x99;
        const double x130 = 0.045230952229266833*x40;
        const double x131 = pow(x[1], -2);
        const double x132 = pow(x27, -0.71612265084075177);
        const double x133 = x132*x25;
        const double x134 = x131*x133;
        const double x135 = x134*x[2];
        const double x136 = 0.28387734915924823*x135;
        const double x137 = 6208.2427845157281*x92;
        const double x138 = x137*x29;
        const double x139 = 1.2464392646590381*x128;
        const double x140 = -u[0]*x138 - x136*x21 + x139*x93;
        const double x141 = x140*x64;
        const double x142 = x141*x65;
        const double x143 = x141*x70;
        const double x144 = x142*x[2] + x143*x[2] + x3*x[0];
        const double x145 = 0.018144065864949988*x144;
        const double x146 = x134*x40;
        const double x147 = x146*x[2];
        const double x148 = 0.0051506893207128044*x147;
        const double x149 = x40*x93;
        const double x150 = 31.976729810597487*x149;
        const double x151 = x128*x98;
        const double x152 = 225.28553197570776*x29;
        const double x153 = x152*x40;
        const double x154 = x140*x63;
        const double x155 = 1.2464392646590381*x93;
        const double x156 = x47*x50;
        const double x157 = 12416.485569031456*x29;
        const double x158 = 1.6002751389984484*x57;
        const double x159 = x49*x93;
        const double x160 = 1762.3795046053549*x93;
        const double x161 = x132*x21;
        const double x162 = x161*x73;
        const double x163 = x131*x[2];
        const double x164 = 0.080586349365681736*x163;
        const double x165 = 2.4928785293180762*x129 - x136*x47 - 2.4928785293180762*x151 + x155*x156 + x157*x98 - x157*x99 - x158*x159 + x160*x51 + x162*x164;
        const double x166 = 0.0051506893207128044*x21;
        const double x167 = x135*x166;
        const double x168 = 0.0014621640307067999*x40;
        const double x169 = x168*x73;
        const double x170 = x161*x163;
        const double x171 = 112.64276598785388*x29;
        const double x172 = x171*x93;
        const double x173 = 0.022615476114633416*x128;
        const double x174 = x173*x93;
        const double x175 = 0.022615476114633416*x149;
        const double x176 = x8*pow(x[0], 3)/pow(M_PI, 3);
        const double x177 = x128*x176;
        const double x178 = x177*x40;
        const double x179 = pow(u[0], 2)*pow(x[0], 4)/(pow(M_PI, 4)*pow(x6, 4));
        const double x180 = x179*x50;
        const double x181 = x149*x41;
        const double x182 = x21*pow(x[2], 2);
        const double x183 = x25*pow(x27, -1.7161226508407519);
        const double x184 = x182*x183;
        const double x185 = x184/pow(x[1], 4);
        const double x186 = 0.0036885252900060047*x40;
        const double x187 = x144*x93;
        const double x188 = x133*x21;
        const double x189 = pow(x[1], -3);
        const double x190 = x189*x[2];
        const double x191 = x188*x190;
        const double x192 = x176*x29;
        const double x193 = x41*x93;
        const double x194 = 3524.7590092107098*x135*x93 - 0.70767174867882077*x170*x193 - 3.594851481170859*x177 + 7.0264350904108897*x179*x57 - 15476.395142313126*x180 - 0.20329099979356649*x185 + 0.56775469831849645*x191 + 17905.173081569799*x192;
        const double x195 = pow(x140, 2)*x81;
        const double x196 = x192*x40;
        const double x197 = x135*x21;
        const double x198 = x40*x92;
        const double x199 = 1024.360059445095*x176;
        const double x200 = x101*x41;
        const double x201 = x161*x200;
        const double x202 = x101*x50;
        const double x203 = x50*x94;
        const double x204 = x101*x93;
        const double x205 = -x114*x135 + x155*x203 - x158*x204 + x160*x202 + x164*x201;
        const double x206 = -x138 + x139*x92 - 0.20566247866874127*x177 + x199*x29 + x205;
        const double x207 = -x116*x140 - x117*x140;
        const double x208 = x204*x58;
        const double x209 = x104*x145 - x108*x167 - x108*x172 + x108*x174 - x110*x144 - x148*x94 + x163*x168*x201 + x175*x203;
        const double x210 = x119*x17*x93;
        const double x211 = x120*x31;
        const double x212 = x141*x[2];
        const double x213 = x126*x212 - 1.7524372332241536*x128*x210 + x134*x182*x211 + 8728.5085739469596*x210*x29;
        const double x214 = x161*x26;
        const double x215 = x133*x26;
        const double x216 = x166*x215;
        const double x217 = 1.0/x37;
        const double x218 = pow(x[2], -3.0/2.0);
        const double x219 = pow(x35, 0.57057057057057059);
        const double x220 = pow(x35, 0.78528528528528541);
        const double x221 = -0.33428772384703404*x218*x219 + 0.33811420958044264*x218*x220;
        const double x222 = x217*x221;
        const double x223 = 43.795847987340309*x119;
        const double x224 = x37*x64;
        const double x225 = 12.432649230827563*x224;
        const double x226 = x188*x27;
        const double x227 = 0.25614106532433373*pow(x[2], -2.7852852852852852) - 0.28460556796052533*pow(x[2], -2.5705705705705704);
        const double x228 = 1.0/x17;
        const double x229 = x227*x228;
        const double x230 = u[1]*x229;
        const double x231 = 38.744203843894169*x119;
        const double x232 = 10.998601882490233*x226;
        const double x233 = x64*x82;
        const double x234 = -x222*x223 + x225*x226 - x230*x231 + x232*x233 - x34 - x39;
        const double x235 = x162*x27;
        const double x236 = 3.1222539468592445*x233;
        const double x237 = x226*x80;
        const double x238 = 16.497902823735348*x82;
        const double x239 = x237*x238;
        const double x240 = 19.372101921947085*x230;
        const double x241 = x133*x27;
        const double x242 = 10.998601882490233*x233;
        const double x243 = x241*x242;
        const double x244 = 3.5293475066740951*x224;
        const double x245 = x225*x241;
        const double x246 = 18.648973846241343*x37;
        const double x247 = x237*x246;
        const double x248 = 21.897923993670155*x222;
        const double x249 = x248*x[2];
        const double x250 = x241*x93;
        const double x251 = x161*x193*x27;
        const double x252 = x184*x189;
        const double x253 = 31.976729810597483*x149;
        const double x254 = x131*x183*x21;
        const double x255 = x21*x215;
        const double x256 = pow(x21, 2)*pow(x24, -0.56775469831849645);
        const double x257 = x131*x256*pow(x27, -1.4322453016815035)*x81;
        const double x258 = x64*x[2];
        const double x259 = x254*x258;
        const double x260 = pow(x[2], -5.0/2.0);
        const double x261 = pow(x[2], -3);
        const double x262 = x201*x27;
        const double x263 = x105*x237;
        const double x264 = x105*x258;
        const double x265 = x106 + x107 - x238*x263 + x240*x264 + x248*x264;
        const double x266 = x200*x214;
        const double x267 = x108*x216 - x110*x234 + x113*x234 - x168*x266 + x215*x96;
        const double x268 = -x211*x226;
        const double x269 = x17*x64;
        const double x270 = 0.72447653635574016*x176 - 4.3907668870044869*x92;
        const double x271 = x29*(-x137 + x199);
        const double x272 = x205 - x270*x62 + x271;
        const double x273 = x114*x215 - 0.080586349365681736*x266;
        const double x274 = x101*x203;
        const double x275 = x29*(58.603931199438705*x176 - 710.35068120531776*x92);
        const double x276 = 0.041447509301757515*x176 - 0.50239405214251542*x92;
        const double x277 = pow(x101, 2);
        const double x278 = -0.56775469831849645*x274 + x275 - x276*x62 + x277*x75;
        const double x279 = pow(x105, 2);
		return {std::vector<double>{-x14*x45 - x52*x53 + x54*x55 + x56*x59 + 0.036288131729899975*x60*x72 - x72*x73*x74 + x88*(x76*x78 + x76*x85 - x79*x84 - x79*x87), -x111*x144 + x129*x130 - x130*x151 + x145*x60 - x148*x47 + x150*x51 + x153*x98 - x153*x99 + x156*x175 - 0.029035497524029844*x159*x58 - x167*x72 + x169*x170 - x172*x72 + x174*x72 + x88*(-x154*x84 - x154*x87 + x165*x78 + x165*x85 + x3), 0.045230952229266833*x128*x187 - 0.010301378641425609*x144*x197 + 63.953459621194966*x147*x93 - x152*x187 - 0.012840042818792856*x170*x181 - 0.065225222049077089*x178 + 0.127488101076211*x179*x58 - 280.80473281412139*x180*x40 - x185*x186 + x191*x53 + 324.87263971533196*x196 + x88*(x194*x78 + x194*x85 - x195*x83 - x195*x86), x109*x234 - x111*x234 - x169*x214 + x215*x47*x95 + x216*x72 + x88*(x124*x240 - x235*x236 - x235*x244 - x239*x63 + x243*x47 + x245*x47 - x247*x63 + x249*x68 + x67 + x71), x144*x216 - x146*x166 - x167*x234 - x172*x234 + x174*x234 + 0.0064200214093964279*x181*x214 + x183*x186*x190*x21 - x215*x253 + x88*(-x140*x239 - x140*x247 + x141*x249 + x142 + x143 - x197*x225 - x197*x242 + x212*x240 - 77184.904879700232*x224*x250 + 15.496542165036463*x224*x251 + 8.9033017241534687*x224*x252 - 68281.990776731094*x233*x250 + 13.709089242688638*x233*x251 + 7.8763479356309878*x233*x252), -x186*x254 + 0.010301378641425609*x234*x255 + x88*(19.372101921947085*u[1]*x119*pow(x227, 2)/pow(x16, 3.0/2.0) + 21.897923993670155*x119*pow(x221, 2)/pow(x36, 3.0/2.0) - x119*x228*x69*(-0.7134259402051637*pow(x[2], -3.7852852852852852) + 0.73159869721984883*pow(x[2], -3.5705705705705704)) - x217*x223*(0.50143158577055102*x219*x260 - 0.50717131437066398*x220*x260 + 0.095367368665069771*x261*pow(x35, -0.42942942942942941) - 0.13275805676469335*x261*pow(x35, -0.21471471471471459)) + 12.432649230827563*x222*x226*x64 - x222*x38 + 24.865298461655126*x224*x255 + x230*x232*x64 - x230*x32 + 21.997203764980465*x233*x255 - 5.2940212600111423*x257*x37 - 4.683380920288867*x257*x82 - 8.9033017241534687*x259*x37 - 7.8763479356309878*x259*x82)}, {x118, x202*x253 - 0.029035497524029841*x208 + x209 - x270*x45 + x271*x55 + x88*(x207 + x272*x78 + x272*x85), x267 + x88*(-18.648973846241347*x263*x37 + x265 + x273*x78 + x273*x85), x127, x213, -x119*x125*x229 - x122*x30*x31 + 0.19955935697826568*x256*x269*pow(x27, 0.56775469831849645) + x268}, {0.036288131729899975*x104*x108 - x108*x200*x74 - x274*x53 + x275*x55 - x276*x45 + x277*x59 + x88*(x278*x78 + x278*x85 - x279*x84 - x279*x87), -x104*x123 + x121*x200 + x126*x264}, {}, {}, {}};
	}
private:
	F2dieselMotor(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F3dieselMotor : public Expression {
public:
	static std::unique_ptr<F3dieselMotor> create() {
		Adjacency adj{{0, 1, 2, 3}, {0}, {}};
		AdjacencyDiff adjDiff{{{0, 0}, {1, 0}, {1, 1}, {2, 0}, {2, 1}, {2, 2}, {3, 0}, {3, 1}, {3, 2}, {3, 3}}, {{0, 0}, {0, 1}, {0, 2}, {0, 3}}, {{0, 0}}, {}, {}, {}};
		return std::unique_ptr<F3dieselMotor>(new F3dieselMotor(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
        const double x0 = pow(x[3], 2);
        const double x1 = 1.0/x[2];
        const double x2 = sqrt(x1);
        const double x3 = x[0]/M_PI;
        const double x4 = u[0]*x3;
        const double x5 = x4/(1.4420483946140072*x3*x[1] + 0.08249999999999999*x4);
        const double x6 = 0.72900098370469424*pow(x[2]/x[1], 0.28387734915924823)*pow(x5 + 0.32842745509493809, -0.28387734915924823)*(4305.1556436685933*x5 + 930.65395757192096);
		return -0.24723010996875103*x0 + 5.0561342312332735e-5*(19879.505327143732*sqrt(x6)*x[2]*(1 - 0.79174558284615615*pow(x1, 0.21471471471471473))*sqrt(pow(x2, 1.5705705705705706) - 0.88980086696190397*pow(x2, 1.7852852852852854)) - 602154.98178743932*(1.2136487758138137*pow(x[1], 0.28387734915924828) - 1)*sqrt(-pow(x[1], 2)*pow(0.38110788300674586*x0 + 1, -7.0452961672473862) + 0.25558758431328077))/x[3];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        const double x0 = M_1_PI;
        const double x1 = u[0]*x0;
        const double x2 = 1.4420483946140072*x0*x[1];
        const double x3 = 0.08249999999999999*x1;
        const double x4 = x2*x[0] + x3*x[0];
        const double x5 = 1.0/x4;
        const double x6 = 3.0448124372273666*x5;
        const double x7 = x1*x6;
        const double x8 = pow(x4, -2);
        const double x9 = x1*x8*x[0]*(x2 + x3);
        const double x10 = x7 - 3.0448124372273666*x9;
        const double x11 = 1.0/x[2];
        const double x12 = sqrt(x11);
        const double x13 = 0.65244887140199537*sqrt(pow(x12, 1.5705705705705706) - 0.88980086696190397*pow(x12, 1.7852852852852854));
        const double x14 = 4305.1556436685933*x5;
        const double x15 = x1*x14;
        const double x16 = x15*x[0] + 930.65395757192096;
        const double x17 = x7*x[0] + 1;
        const double x18 = pow(x17, -0.28387734915924823);
        const double x19 = x[2]/x[1];
        const double x20 = pow(x19, 0.28387734915924823);
        const double x21 = x18*x20;
        const double x22 = x16*x21;
        const double x23 = pow(x22, -1.0/2.0);
        const double x24 = x13*x23;
        const double x25 = x24*x[2];
        const double x26 = pow(x17, -1.2838773491592481)*x20;
        const double x27 = x16*x26;
        const double x28 = 1 - 0.79174558284615615*pow(x11, 0.21471471471471473);
        const double x29 = 8649.4766444153556*x28;
        const double x30 = x25*x27*x29;
        const double x31 = x21*(x15 - 4305.1556436685933*x9);
        const double x32 = 0.28387734915924823*x16;
        const double x33 = x26*x32;
        const double x34 = 15234.531162898826*x22*x28*x[2];
        const double x35 = x13/pow(x22, 3.0/2.0);
        const double x36 = x34*x35;
        const double x37 = 30469.062325797651*x28;
        const double x38 = x25*x37;
        const double x39 = 5.0561342312332735e-5/x[3];
        const double x40 = u[0]*x8*pow(x[0], 2)/pow(M_PI, 2);
        const double x41 = x21*x40;
        const double x42 = x25*x28;
        const double x43 = pow(x[1], 2);
        const double x44 = x18*pow(x19, -0.71612265084075177);
        const double x45 = x44/x43;
        const double x46 = x16*x24*x29;
        const double x47 = pow(x[3], 2);
        const double x48 = 0.38110788300674586*x47 + 1;
        const double x49 = pow(x48, -7.0452961672473862);
        const double x50 = 1.2136487758138137*pow(x[1], 0.28387734915924828) - 1;
        const double x51 = 1.9780174716835943*sqrt(-x43*x49 + 0.25558758431328077);
        const double x52 = x50/x51;
        const double x53 = x27*x40;
        const double x54 = x22*x24;
        const double x55 = x37*x54;
        const double x56 = pow(x[2], -3.0/2.0);
        const double x57 = x0*x[0];
        const double x58 = -0.25119702607125771*x40 + x57*x6;
        const double x59 = x21*(x14*x57 - 355.17534060265888*x40);
		return {std::vector<double>{x39*(-x10*x30 + x31*x38 - x36*(-x10*x33 + x31)), x39*(-x36*(-x32*x45*x[2] - 6208.2427845157281*x41 + 1.2464392646590381*x53) - 189159336.33509329*x41*x42 + 37977.835640217629*x42*x53 - x45*x46*pow(x[2], 2) + 1191073.0746368715*x49*x52*x[1] - 104882.23264422762*x51*pow(x[1], -0.71612265084075166)), x39*(-4324.7383222076778*pow(x16, 2)*pow(x17, -0.56775469831849645)*pow(x19, 0.56775469831849645)*x28*x35 + x19*x44*x46 + 5179.7231350116763*x54*pow(x[2], -0.21471471471471482) + x55 + x23*x34*(-0.33428772384703404*pow(x12, 0.57057057057057059)*x56 + 0.33811420958044264*pow(x12, 0.78528528528528541)*x56)/x13), -161.69782892905661*x43*pow(x48, -8.0452961672473862)*x52 - 0.49446021993750205*x[3] - 5.0561342312332735e-5*(-304423.48988702998*x50*x51 + x55*x[2])/x47}, {x39*(-x30*x58 - x36*(-x33*x58 + x59) + x38*x59)}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        const double x0 = M_1_PI;
        const double x1 = u[0]*x0;
        const double x2 = 8610.3112873371865*x1;
        const double x3 = 1.4420483946140072*x0*x[1];
        const double x4 = 0.08249999999999999*x1;
        const double x5 = x3*x[0] + x4*x[0];
        const double x6 = pow(x5, -3);
        const double x7 = x3 + x4;
        const double x8 = x6*pow(x7, 2)*x[0];
        const double x9 = pow(x5, -2);
        const double x10 = x7*x9;
        const double x11 = x0/x5;
        const double x12 = 3.0448124372273666*x11;
        const double x13 = x12*x[0];
        const double x14 = u[0]*x13 + 1;
        const double x15 = pow(x14, -0.28387734915924823);
        const double x16 = 1.0/x[1];
        const double x17 = x16*x[2];
        const double x18 = pow(x17, 0.28387734915924823);
        const double x19 = x15*x18;
        const double x20 = x19*(-x10*x2 + x2*x8);
        const double x21 = 1.0/x[2];
        const double x22 = 1 - 0.79174558284615615*pow(x21, 0.21471471471471473);
        const double x23 = 30469.062325797651*x22;
        const double x24 = sqrt(x21);
        const double x25 = 0.42568952979373742*pow(x24, 1.5705705705705706) - 0.37877891266707281*pow(x24, 1.7852852852852854);
        const double x26 = sqrt(x25);
        const double x27 = 4305.1556436685933*x11;
        const double x28 = u[0]*x27;
        const double x29 = x28*x[0] + 930.65395757192096;
        const double x30 = x19*x29;
        const double x31 = pow(x30, -1.0/2.0);
        const double x32 = x26*x31;
        const double x33 = x32*x[2];
        const double x34 = x23*x33;
        const double x35 = x10*x[0];
        const double x36 = x1*x35;
        const double x37 = x28 - 4305.1556436685933*x36;
        const double x38 = x19*x37;
        const double x39 = u[0]*x12 - 3.0448124372273666*x36;
        const double x40 = x18*x29;
        const double x41 = pow(x14, -1.2838773491592481);
        const double x42 = 0.28387734915924823*x41;
        const double x43 = x40*x42;
        const double x44 = x38 - x39*x43;
        const double x45 = pow(x30, -3.0/2.0);
        const double x46 = x26*x45;
        const double x47 = x38*x46;
        const double x48 = x23*x[2];
        const double x49 = pow(x39, 2);
        const double x50 = pow(x14, -2.2838773491592481)*x40;
        const double x51 = 0.36446369852492994*x50;
        const double x52 = 6.0896248744547332*x1;
        const double x53 = -x10*x52 + x52*x8;
        const double x54 = 0.56775469831849645*x41;
        const double x55 = x18*x37;
        const double x56 = x39*x55;
        const double x57 = x30*x46;
        const double x58 = 15234.531162898826*x22;
        const double x59 = x57*x58;
        const double x60 = x59*x[2];
        const double x61 = x26/pow(x30, 5.0/2.0);
        const double x62 = x22*x[2];
        const double x63 = x30*x62;
        const double x64 = 22851.796744348238*x61*x63;
        const double x65 = x22*x33;
        const double x66 = 11104.867145846816*x50*x65;
        const double x67 = 17298.953288830711*x22;
        const double x68 = x33*x41;
        const double x69 = x67*x68;
        const double x70 = x40*x41;
        const double x71 = 8649.4766444153556*x22;
        const double x72 = x33*x71;
        const double x73 = x70*x72;
        const double x74 = x39*x41;
        const double x75 = x40*x74;
        const double x76 = x71*x75;
        const double x77 = x44*x46;
        const double x78 = 1.0/x[3];
        const double x79 = 5.0561342312332735e-5*x78;
        const double x80 = x58*x[2];
        const double x81 = pow(x[0], 2);
        const double x82 = pow(M_PI, -2);
        const double x83 = x82*x9;
        const double x84 = x81*x83;
        const double x85 = u[0]*x84;
        const double x86 = x27*x[0] - 355.17534060265888*x85;
        const double x87 = x19*x86;
        const double x88 = x46*x87;
        const double x89 = x13 - 0.25119702607125771*x85;
        const double x90 = x39*x89;
        const double x91 = x18*x86;
        const double x92 = x39*x91;
        const double x93 = x0*x35;
        const double x94 = 0.50239405214251542*u[0];
        const double x95 = x7*x81*x82;
        const double x96 = x6*x95;
        const double x97 = x83*x[0];
        const double x98 = x12 - 3.0448124372273666*x93 + x94*x96 - x94*x97;
        const double x99 = x55*x89;
        const double x100 = 710.35068120531776*u[0];
        const double x101 = x19*(x100*x96 - x100*x97 + x27 - 4305.1556436685933*x93);
        const double x102 = x38*x80;
        const double x103 = -x43*x89 + x87;
        const double x104 = x103*x46;
        const double x105 = x44*x64;
        const double x106 = x41*x89;
        const double x107 = x106*x40;
        const double x108 = 4324.7383222076778*x62;
        const double x109 = x107*x108;
        const double x110 = x108*x75;
        const double x111 = x68*x71;
        const double x112 = x79*(x101*x34 - x102*x104 + x103*x105 + x104*x110 + x109*x77 - x111*x92 - x111*x99 - x44*x80*x88 - x60*(x101 - x42*x92 - x42*x99 - x43*x98 + x51*x90) + x66*x90 - x73*x98);
        const double x113 = 0.28387734915924823*x15;
        const double x114 = pow(x[1], 2);
        const double x115 = 1.0/x114;
        const double x116 = pow(x17, -0.71612265084075177);
        const double x117 = x116*x[2];
        const double x118 = x115*x117;
        const double x119 = x113*x118;
        const double x120 = 6208.2427845157281*x84;
        const double x121 = x120*x19;
        const double x122 = 1.2464392646590381*x41;
        const double x123 = x122*x40;
        const double x124 = -u[0]*x121 - x119*x29 + x123*x85;
        const double x125 = x124*x46;
        const double x126 = x124*x80;
        const double x127 = x19*x85;
        const double x128 = 94579668.167546645*x127;
        const double x129 = x128*x62;
        const double x130 = x18*x85;
        const double x131 = 53698050.967528947*x130*x65;
        const double x132 = x15*x29;
        const double x133 = x116*x132;
        const double x134 = x115*x133;
        const double x135 = pow(x[2], 2);
        const double x136 = x135*x22;
        const double x137 = 4324.7383222076778*x134*x136;
        const double x138 = u[0]*x97;
        const double x139 = x138*x65;
        const double x140 = 378318672.67018658*x19;
        const double x141 = x15*x32;
        const double x142 = x116*x141;
        const double x143 = x135*x71;
        const double x144 = x115*x143;
        const double x145 = 75955.671280435257*x70;
        const double x146 = u[0]*x6;
        const double x147 = x146*x95;
        const double x148 = x147*x65;
        const double x149 = x116*x135;
        const double x150 = x22*x32;
        const double x151 = x115*x149*x150;
        const double x152 = 2455.3905014314605*x151;
        const double x153 = x29*x74;
        const double x154 = 37977.835640217629*x22*x68*x85;
        const double x155 = x70*x85;
        const double x156 = 18988.917820108814*x155;
        const double x157 = x156*x62;
        const double x158 = x50*x85;
        const double x159 = x158*x65;
        const double x160 = 48758.882948568229*x159;
        const double x161 = 2.4928785293180762*x70;
        const double x162 = x122*x85;
        const double x163 = 12416.485569031456*x19;
        const double x164 = 1.6002751389984484*x158;
        const double x165 = 1762.3795046053549*x130;
        const double x166 = x118*x29;
        const double x167 = pow(u[0], 2)*pow(x[0], 4)/(pow(M_PI, 4)*pow(x5, 4));
        const double x168 = x167*x50;
        const double x169 = x167*x18*x41;
        const double x170 = pow(x[1], 3);
        const double x171 = 1.0/x170;
        const double x172 = x132*x171;
        const double x173 = x146*pow(x[0], 3)/pow(M_PI, 3);
        const double x174 = x173*x19;
        const double x175 = x41*x85;
        const double x176 = x173*x70;
        const double x177 = pow(x[1], 4);
        const double x178 = 1.0/x177;
        const double x179 = x132*pow(x17, -1.7161226508407519);
        const double x180 = x135*x179;
        const double x181 = 189159336.33509329*x19;
        const double x182 = x181*x85;
        const double x183 = x125*x62;
        const double x184 = pow(x[3], 2);
        const double x185 = 0.38110788300674586*x184 + 1;
        const double x186 = pow(x185, -7.0452961672473862);
        const double x187 = -3.9125531182855591*x114*x186 + 1;
        const double x188 = sqrt(x187);
        const double x189 = x176*x65;
        const double x190 = x142*x85;
        const double x191 = x32*x67;
        const double x192 = 1.0/x188;
        const double x193 = x186*x192;
        const double x194 = x174*x65;
        const double x195 = x134*x143;
        const double x196 = pow(x[2], 3);
        const double x197 = 6194.0861429838951*x22;
        const double x198 = x179*x197;
        const double x199 = 1.2136487758138137*pow(x[1], 0.28387734915924828) - 1;
        const double x200 = x199/pow(x187, 3.0/2.0);
        const double x201 = 1191073.0746368715*x193*x199;
        const double x202 = x175*x29;
        const double x203 = 37977.835640217629*x155;
        const double x204 = 1024.360059445095*x173;
        const double x205 = 0.080586349365681736*x106;
        const double x206 = x113*x86;
        const double x207 = x106*x165 - x118*x206 + x162*x91 - x164*x89 + x166*x205;
        const double x208 = x65*x84;
        const double x209 = x106*x29;
        const double x210 = x142*x86;
        const double x211 = x103*x124*x64 + x104*x129 + x104*x137 - x104*x157 + x106*x131 + x109*x125 - x126*x88 - x144*x210 + x152*x209 + x154*x91;
        const double x212 = pow(x14, -1.5677546983184965);
        const double x213 = x17*x22;
        const double x214 = 2455.3905014314605*x213;
        const double x215 = x212*x214;
        const double x216 = pow(x29, 2);
        const double x217 = pow(x17, -0.43224530168150355);
        const double x218 = x216*x217*x46;
        const double x219 = x17*x71;
        const double x220 = x219*x37;
        const double x221 = pow(x14, -0.56775469831849645);
        const double x222 = x217*x221;
        const double x223 = x222*x29*x46;
        const double x224 = x44*x59;
        const double x225 = pow(x[2], -0.21471471471471482);
        const double x226 = 1470.4060729459454*x32;
        const double x227 = x225*x226;
        const double x228 = 4324.7383222076778*x213;
        const double x229 = x133*x228;
        const double x230 = x32*x38;
        const double x231 = x23*x230;
        const double x232 = 1.0/x26;
        const double x233 = pow(x[2], -3.0/2.0);
        const double x234 = pow(x24, 0.57057057057057059);
        const double x235 = pow(x24, 0.78528528528528541);
        const double x236 = -0.33428772384703404*x233*x234 + 0.33811420958044264*x233*x235;
        const double x237 = x232*x236;
        const double x238 = x237*x31;
        const double x239 = x116*x32;
        const double x240 = x214*x239;
        const double x241 = x237*x45;
        const double x242 = 7617.2655814494128*x241*x63;
        const double x243 = 2589.8615675058381*x225;
        const double x244 = x243*x57;
        const double x245 = 5179.7231350116763*x225;
        const double x246 = x216*x222;
        const double x247 = x213*x246*x61;
        const double x248 = 6487.1074833115163*x247;
        const double x249 = x216*x221*x46;
        const double x250 = x217*x249;
        const double x251 = 4324.7383222076778*x250;
        const double x252 = x115*x62;
        const double x253 = x225*x32;
        const double x254 = 10781.0473083506*x213;
        const double x255 = pow(x17, -1.4322453016815035)*x249;
        const double x256 = x124*x57;
        const double x257 = x256*x58;
        const double x258 = x238*x62;
        const double x259 = 53698050.967528947*x213;
        const double x260 = x16*x250;
        const double x261 = x30*x32;
        const double x262 = x30*x31;
        const double x263 = x237*x262;
        const double x264 = x116*x16;
        const double x265 = x264*x29;
        const double x266 = x141*x265;
        const double x267 = x133*x219;
        const double x268 = pow(x[2], -5.0/2.0);
        const double x269 = 1.0/x196;
        const double x270 = 1470.4060729459457*x225;
        const double x271 = x216*x217*x46*x89;
        const double x272 = x223*x86;
        const double x273 = x103*x59;
        const double x274 = x107*x32;
        const double x275 = x32*x87;
        const double x276 = x23*x275;
        const double x277 = -x103*x242 - x103*x244 - x104*x229 - x109*x238 - x209*x240 + x210*x219 + x238*x80*x87 + x245*x275 - x273 - x274*x71 + x276;
        const double x278 = 1.0/x184;
        const double x279 = 5.0561342312332735e-5*x278;
        const double x280 = pow(x185, -8.0452961672473862)*x192;
        const double x281 = x199*x280;
        const double x282 = x23*x261;
        const double x283 = -x279*(-x107*x72 - x273*x[2] + x276*x[2]);
        const double x284 = 0.72447653635574016*x173 - 4.3907668870044869*x84;
        const double x285 = x19*(-x120 + x204);
        const double x286 = pow(x89, 2);
        const double x287 = 0.041447509301757515*x173 - 0.50239405214251542*x84;
        const double x288 = x89*x91;
        const double x289 = x19*(58.603931199438705*x173 - 710.35068120531776*x84);
		return {std::vector<double>{x79*(x20*x34 + pow(x44, 2)*x64 - x44*x47*x48 + x49*x66 - x53*x73 - x56*x69 - x60*(x20 - x43*x53 + x49*x51 - x54*x56) + x76*x77*x[2]), x79*(x105*x124 + x110*x125 - x126*x47 + x129*x77 + x131*x74 + x137*x77 - x139*x140 + x139*x145 + x140*x148 - x142*x144*x37 - x145*x148 + x152*x153 + x154*x55 - x157*x77 - x160*x39 - x60*(-x119*x37 + x138*x161 - x138*x163 - x147*x161 + x147*x163 + x162*x55 - x164*x39 + x165*x74 + 0.080586349365681736*x166*x74)), x79*(4660136.6722764596*x114*pow(x185, -14.090592334494772)*x200 + 107396101.93505789*x115*x136*x190 + pow(x124, 2)*x64 + x125*x195 + x149*x172*x191 - 21562.094616701201*x151*x202 + 214088.88869790107*x168*x65 - 471551248.16981071*x169*x65 - x178*x196*x198*x32 + x182*x183 - x183*x203 + 75108.542467280698*x188*pow(x[1], -1.7161226508407517) - 109531.75383178092*x189 + 820714.61276984843*x193*pow(x[1], 0.28387734915924834) + 545553834.57654464*x194 + x201 - x60*(0.56775469831849645*x117*x172 + 3524.7590092107098*x118*x15*x85 - 0.70767174867882077*x166*x175 + 7.0264350904108897*x168 - 15476.395142313126*x169 + 17905.173081569799*x174 - 3.594851481170859*x176 - 0.20329099979356649*x178*x180)), x79*(x102*x238 - x110*x238 + x142*x220 - x153*x240 + x215*x218*x39 - x220*x223 - x224 - x227*x75 - x229*x77 + x230*x245 + x231 - x242*x44 - x244*x44 + x248*x44 - x32*x76), x79*(-x118*x132*x191 - x124*x242 + x124*x248 - x125*x229 - 32156978.778725427*x127*x253 - x128*x258 - x134*x226*pow(x[2], 0.78528528528528518) - 1869.3478207762171*x136*x171*x255 - x137*x238 - x150*x182 + x150*x203 + 6456.2102955413611*x155*x253 + x156*x258 + x171*x180*x197*x32 - x190*x259 + x202*x239*x254 - x212*x218*x254*x85 + x223*x259*x85 - x243*x256 + x251*x252 - x257), x79*(-x115*x198*x33 + 1841.5428760735954*pow(x14, -0.85163204747774468)*pow(x17, -1.1483679525222552)*x252*pow(x29, 3)*x61 + 2940.8121458918913*x225*x266 - x228*x241*x246 + x23*x263 + x232*x262*x80*(0.50143158577055102*x234*x268 - 0.50717131437066398*x235*x268 + 0.095367368665069771*pow(x24, -0.42942942942942941)*x269 - 0.13275805676469335*pow(x24, -0.21471471471471459)*x269) - 7617.2655814494128*pow(x236, 2)*x262*x62/pow(x25, 3.0/2.0) + x238*x267 + x245*x263 + 641.65257006048682*x252*x255 - x260*x270 - x260*x71 + 4067.5603597764357*x261*pow(x[2], -1.2147147147147148) + x266*x67), -x279*(-x224*x[2] + x231*x[2] - x33*x76), -632.65134479638527*x170*pow(x185, -15.090592334494772)*x200 - 5.0561342312332735e-5*x278*(-x182*x65 - 104882.23264422762*x188*pow(x[1], -0.71612265084075166) - x195*x32 + x201*x[1] + x203*x65 - x257*x[2]) - 55.709332148111571*x280*pow(x[1], 1.2838773491592483) - 323.39565785811322*x281*x[1], -x279*(-x213*x251 + x245*x261 + x263*x80 + x267*x32 + x282), 991.57176708216582*x114*pow(x185, -9.0452961672473862)*x192*x199*x[3] + 161.69782892905661*x114*x281*x78 + 1698.6801899539037*x177*pow(x185, -16.090592334494772)*x200*x[3] - 0.49446021993750205 + 0.00010112268462466547*(-304423.48988702998*x188*x199 + x282*x[2])/pow(x[3], 3)}, {x112, x79*(-48758.882948568222*x159*x89 + x211 - x284*x73 + x285*x34 - x60*(x207 - x284*x43 + x285)), x79*(6487.1074833115172*x103*x247 + 1227.6952507157303*pow(x14, -1.5677546983184962)*x213*x271 - x228*x272 - x270*x274 + x277 - x60*(-x205*x265 + x206*x264)), x283}, {x79*(pow(x103, 2)*x64 - x103*x48*x88 + x104*x107*x71*x[2] + x286*x66 - x287*x73 - x288*x69 + x289*x34 - x60*(x286*x51 - x287*x43 - x288*x54 + x289))}, {}, {}, {}};
	}
private:
	F3dieselMotor(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


Problem createProblem_dieselMotor() {

    std::vector<std::unique_ptr<Expression>> F;
    F.push_back(F0dieselMotor::create());
    F.push_back(F1dieselMotor::create());
    F.push_back(F2dieselMotor::create());
    F.push_back(F3dieselMotor::create());
    
    std::vector<std::unique_ptr<Constraint>> G;
    
    
    std::vector<std::unique_ptr<Constraint>> R;
    
    
    std::vector<std::unique_ptr<ParamConstraint>> A;
    

    Problem problem(
            4, 2, 0,  // #vars
            {0.24989941562646081, 0.50614999999999999, 0.33926666666666666, 0.068099999999999994},  // x0
            {0.018181818181818181, 0.40444536585365859, 0.3370378048780488, 0.029999999999999999},  // lb x
            {1.0, 1.0111134146341463, 1.0111134146341465, 1.0},  // ub x
            {0, 0},  // lb u
            {1.0, 1},  // ub u
            {},  // lb p
            {},  // ub p
            MayerdieselMotor::create(),
            LagrangedieselMotor::create(),
            std::move(F),
            std::move(G),
            std::move(R),
            std::move(A),
            "dieselMotor");
    return problem;
};

int main() {
    auto problem = std::make_shared<const Problem>(createProblem_dieselMotor());
    InitVars initVars = INIT_VARS;
    Integrator rk = Integrator::radauIIA(RADAU_INTEGRATOR);
    Mesh mesh = Mesh::createEquidistantMesh(INTERVALS, FINAL_TIME);
    LinearSolver linearSolver = LINEAR_SOLVER;
    MeshAlgorithm meshAlgorithm = MESH_ALGORITHM;
    int meshIterations = MESH_ITERATIONS;

    Solver solver = Solver(create_gdop(problem, mesh, rk, initVars), meshIterations, linearSolver, meshAlgorithm);

    // set solver flags
    #ifdef EXPORT_OPTIMUM_PATH
    solver.setExportOptimumPath(EXPORT_OPTIMUM_PATH);
    #endif
    
    #ifdef EXPORT_HESSIAN_PATH
    solver.setExportHessianPath(EXPORT_HESSIAN_PATH);
    #endif
    
    #ifdef EXPORT_JACOBIAN_PATH
    solver.setExportJacobianPath(EXPORT_JACOBIAN_PATH);
    #endif
    
    #ifdef TOLERANCE
    solver.setTolerance(TOLERANCE);
    #endif
    
    // set solver mesh parameters
    #ifdef LEVEL
    solver.setMeshParameter("level", LEVEL);
    #endif
    
    #ifdef C_TOL
    solver.setMeshParameter("ctol", C_TOL);
    #endif
    
    #ifdef SIGMA
    solver.setMeshParameter("sigma", SIGMA);
    #endif
    
    // optimize
    int status = solver.solve();
    return status;
}        
        
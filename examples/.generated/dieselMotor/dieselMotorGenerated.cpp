
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
		return pow(-0.547055854225991 + x[1], 2) + pow(-0.515309170685596 + x[0], 2) + pow(-0.381048005791294 + x[2], 2) + pow(-0.27144300053768 + x[3], 2);
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{2*(-0.515309170685596 + x[0]), 2*(-0.547055854225991 + x[1]), 2*(-0.381048005791294 + x[2]), 2*(-0.27144300053768 + x[3])}, {}, {}};
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
		return 0.0825*u[0]*x[0]/acos(-1);
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        const double x0 = 0.0825/acos(-1);
		return {std::vector<double>{x0*u[0]}, {x0*x[0]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        const double x0 = 0.0825/acos(-1);
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
        const double x0 = pow(acos(-1), -1);
		return 0.000324675324675325*x0*(27936.4175478121*u[0] - 0.0127*(-200000.0*x[1] + 300000.0*x[2]) - 1270.0*(0.359041972839291 - 0.933527565528316*x0*x[0] + 31.3448352034933*pow(x[0], 2)/pow(acos(-1), 2)));
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        const double x0 = pow(acos(-1), -1);
		return {std::vector<double>{0.000324675324675325*x0*(1185.58000822096*x0 - 79615.8814168731*x[0]/pow(acos(-1), 2)), 0.824675324675325*x0, -1.23701298701299*x0}, {9.07026543760135*x0}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-25.8493121483354/pow(acos(-1), 3)}, {}, {}, {}, {}, {}};
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
		return 20.2505119361145*(-1.44204839461401*x[0]*x[1]/acos(-1) + 0.526906365590249*sqrt(1 - 3.91255311828556*pow(x[1], 2)*pow(1 + 0.381107883006746*pow(x[3], 2), -7.04529616724739)));
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        const double x0 = pow(acos(-1), -1);
        const double x1 = pow(x[1], 2);
        const double x2 = 1 + 0.381107883006746*pow(x[3], 2);
        const double x3 = pow(x2, -7.04529616724739);
        const double x4 = pow(1 - 3.91255311828556*x1*x3, -1.0/2.0);
		return {std::vector<double>{-29.2022182275857*x0*x[1], 20.2505119361145*(-1.44204839461401*x0*x[0] - 2.06154914373464*x4*x3*x[1]), 112.09258517065*pow(x2, -8.04529616724739)*x1*x4*x[3]}, {}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        const double x0 = pow(x[1], 2);
        const double x1 = pow(x[3], 2);
        const double x2 = 1 + 0.381107883006746*x1;
        const double x3 = pow(x2, -7.04529616724739);
        const double x4 = 1 - 3.91255311828556*x0*x3;
        const double x5 = pow(x4, -1.0/2.0);
        const double x6 = pow(x4, -3.0/2.0);
        const double x7 = pow(x2, -8.04529616724739)*x5;
		return {std::vector<double>{-29.2022182275857/acos(-1), 20.2505119361145*(-2.06154914373464*x3*x5 - 8.06592053081789*x0*pow(x2, -14.0905923344948)*x6), 224.185170341301*x7*x[3]*x[1] + 438.568193646118*pow(x2, -15.0905923344948)*x6*x[3]*pow(x[1], 3), 112.09258517065*x0*x7 - 687.379932622565*x0*pow(x2, -9.04529616724739)*x1*x5 - 1177.5634535801*pow(x2, -16.0905923344948)*x1*x6*pow(x[1], 4)}, {}, {}, {}, {}, {}};
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
        const double x0 = pow(x[2], -1);
        const double x1 = x[0]/acos(-1);
        const double x2 = x1*u[0];
        const double x3 = 0.0825*x2 + 1.44204839461401*x1*x[1];
        const double x4 = x2/x3;
        const double x5 = (930.653957571921 + 4305.15564366859*x4)*pow(1 + 3.04481243722737*x4, -0.283877349159248)*pow(x[2]/x[1], 0.283877349159248);
        const double x6 = x[2]/sqrt(x5);
        const double x7 = sqrt(x0);
		return 0.01814406586495*x5*(x3 - 87.5916959746806*x6*sqrt(0.425689529793737*pow(x7, 1.57057057057057) - 0.378778912667073*pow(x7, 1.78528528528529)) - 77.4884076877883*x6*u[1]*sqrt(0.181211575776013*pow(x0, 1.57057057057057) - 0.14347346468125*pow(x0, 1.78528528528529)));
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        const double x0 = pow(acos(-1), -1);
        const double x1 = 1.44204839461401*x0;
        const double x2 = x1*x[0];
        const double x3 = 0.0825*x0;
        const double x4 = x3*x[0];
        const double x5 = x2*x[1] + x4*u[0];
        const double x6 = x0/x5;
        const double x7 = 4305.15564366859*x6;
        const double x8 = x1*x[1] + x3*u[0];
        const double x9 = u[0]/pow(x5, 2);
        const double x10 = x0*x8*x9*x[0];
        const double x11 = 3.04481243722737*x6;
        const double x12 = x11*x[0];
        const double x13 = 1 + x12*u[0];
        const double x14 = pow(x13, -0.283877349159248);
        const double x15 = pow(x[1], -1);
        const double x16 = x15*x[2];
        const double x17 = pow(x16, 0.283877349159248);
        const double x18 = x14*x17;
        const double x19 = x18*(-4305.15564366859*x10 + x7*u[0]);
        const double x20 = -3.04481243722737*x10 + x11*u[0];
        const double x21 = x7*x[0];
        const double x22 = 930.653957571921 + x21*u[0];
        const double x23 = 0.283877349159248*x22;
        const double x24 = pow(x13, -1.28387734915925)*x17;
        const double x25 = x24*x23;
        const double x26 = x19 - x25*x20;
        const double x27 = 43.7958479873403*x[2];
        const double x28 = pow(x[2], -1);
        const double x29 = sqrt(x28);
        const double x30 = sqrt(0.425689529793737*pow(x29, 1.57057057057057) - 0.378778912667073*pow(x29, 1.78528528528529));
        const double x31 = x22*x18;
        const double x32 = pow(x31, -3.0/2.0);
        const double x33 = x30*x32;
        const double x34 = x33*x27;
        const double x35 = 38.7442038438942*x[2];
        const double x36 = sqrt(0.181211575776013*pow(x28, 1.57057057057057) - 0.14347346468125*pow(x28, 1.78528528528529));
        const double x37 = x36*u[1];
        const double x38 = x32*x37;
        const double x39 = x35*x38;
        const double x40 = 0.01814406586495*x31;
        const double x41 = pow(x31, -1.0/2.0);
        const double x42 = 77.4884076877883*x41*x37;
        const double x43 = 87.5916959746806*x41*x30;
        const double x44 = x5 - x42*x[2] - x43*x[2];
        const double x45 = 0.01814406586495*x44;
        const double x46 = 0.0051506893207128*x44;
        const double x47 = x46*x22;
        const double x48 = x47*x24;
        const double x49 = x14*pow(x16, -0.716122650840752);
        const double x50 = x49*x[2]/pow(x[1], 2);
        const double x51 = x9*pow(x[0], 2)/pow(acos(-1), 2);
        const double x52 = x51*x18;
        const double x53 = x51*x24*x22;
        const double x54 = -6208.24278451573*x52 + 1.24643926465904*x53 - x50*x23;
        const double x55 = pow(x[2], -3.0/2.0);
        const double x56 = x49*x22;
        const double x57 = x56*x16;
        const double x58 = x18*(x21 - 355.175340602659*x51);
        const double x59 = x12 - 0.251197026071258*x51;
        const double x60 = x58 - x59*x25;
		return {std::vector<double>{x45*x19 - x48*x20 + (x8 + x34*x26 + x39*x26)*x40, -x50*x47 - 112.642765987854*x52*x44 + 0.0226154761146334*x53*x44 + (x2 + x54*x34 + x54*x39)*x40, x40*(-x42 - x43 + 12.4326492308276*x57*x33 + 10.9986018824902*x57*x38 - x41*x27*(-0.334287723847034*x55*pow(x29, 0.570570570570571) + 0.338114209580443*x55*pow(x29, 0.785285285285285))/x30 - x41*x35*u[1]*(0.256141065324334*pow(x[2], -2.78528528528529) - 0.284605567960525*pow(x[2], -2.57057057057057))/x36) + x56*x46*x15}, {x58*x45 - x59*x48 + (x4 + x60*x34 + x60*x39)*x40, -1.40595477285733*x41*x31*x36*x[2]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        const double x0 = pow(acos(-1), -1);
        const double x1 = 0.0825*x0;
        const double x2 = x1*u[0];
        const double x3 = 1.44204839461401*x0;
        const double x4 = x3*x[1];
        const double x5 = x2 + x4;
        const double x6 = x2*x[0] + x4*x[0];
        const double x7 = pow(x6, -3);
        const double x8 = x7*u[0];
        const double x9 = x0*pow(x5, 2)*x8*x[0];
        const double x10 = x5*u[0];
        const double x11 = pow(x6, -2);
        const double x12 = x0*x11;
        const double x13 = x12*x10;
        const double x14 = -6.08962487445473*x13 + 6.08962487445473*x9;
        const double x15 = pow(x[2], -1);
        const double x16 = 0.181211575776013*pow(x15, 1.57057057057057) - 0.14347346468125*pow(x15, 1.78528528528529);
        const double x17 = sqrt(x16);
        const double x18 = x0/x6;
        const double x19 = 4305.15564366859*x18;
        const double x20 = x19*x[0];
        const double x21 = 930.653957571921 + x20*u[0];
        const double x22 = 3.04481243722737*x18;
        const double x23 = x22*u[0];
        const double x24 = 1 + x23*x[0];
        const double x25 = pow(x24, -0.283877349159248);
        const double x26 = pow(x[1], -1);
        const double x27 = x26*x[2];
        const double x28 = pow(x27, 0.283877349159248);
        const double x29 = x25*x28;
        const double x30 = x21*x29;
        const double x31 = pow(x30, -1.0/2.0);
        const double x32 = 77.4884076877883*x31;
        const double x33 = x32*x17;
        const double x34 = x33*u[1];
        const double x35 = sqrt(x15);
        const double x36 = 0.425689529793737*pow(x35, 1.57057057057057) - 0.378778912667073*pow(x35, 1.78528528528529);
        const double x37 = sqrt(x36);
        const double x38 = 87.5916959746806*x31;
        const double x39 = x38*x37;
        const double x40 = x6 - x34*x[2] - x39*x[2];
        const double x41 = pow(x24, -1.28387734915925);
        const double x42 = x21*x28;
        const double x43 = 0.0051506893207128*x42;
        const double x44 = x41*x43;
        const double x45 = x40*x44;
        const double x46 = 4305.15564366859*x[0];
        const double x47 = x19*u[0] - x46*x13;
        const double x48 = 3.04481243722737*x[0];
        const double x49 = x23 - x48*x13;
        const double x50 = x41*x28;
        const double x51 = x50*x49;
        const double x52 = x51*x47;
        const double x53 = 0.0103013786414256*x40;
        const double x54 = x29*(-8610.31128733719*x13 + 8610.31128733719*x9);
        const double x55 = 0.01814406586495*x40;
        const double x56 = pow(x49, 2);
        const double x57 = x42*pow(x24, -2.28387734915925);
        const double x58 = x57*x40;
        const double x59 = 0.0066128533514196*x58;
        const double x60 = x47*x29;
        const double x61 = 0.283877349159248*x41;
        const double x62 = x61*x42;
        const double x63 = x60 - x62*x49;
        const double x64 = pow(x30, -3.0/2.0);
        const double x65 = 43.7958479873403*x37;
        const double x66 = x64*x65;
        const double x67 = x63*x66;
        const double x68 = x63*x64;
        const double x69 = 38.7442038438942*u[1];
        const double x70 = x69*x17;
        const double x71 = x70*x68;
        const double x72 = x5 + x67*x[2] + x71*x[2];
        const double x73 = x41*x49;
        const double x74 = 0.0103013786414256*x42;
        const double x75 = 0.36446369852493*x57;
        const double x76 = -0.567754698318496*x52 + x54 - x62*x14 + x75*x56;
        const double x77 = x70*x64;
        const double x78 = x77*x[2];
        const double x79 = pow(x63, 2);
        const double x80 = pow(x30, -5.0/2.0);
        const double x81 = x80*x[2];
        const double x82 = x17*u[1];
        const double x83 = 58.1163057658413*x82;
        const double x84 = x81*x83;
        const double x85 = x66*x[2];
        const double x86 = 65.6937719810105*x37;
        const double x87 = x81*x86;
        const double x88 = 0.01814406586495*x30;
        const double x89 = pow(x[0], 2);
        const double x90 = pow(acos(-1), -2);
        const double x91 = x90*x11;
        const double x92 = x89*x91;
        const double x93 = x92*u[0];
        const double x94 = x20 - 355.175340602659*x93;
        const double x95 = 0.0051506893207128*x40;
        const double x96 = x95*x94;
        const double x97 = x5*x12;
        const double x98 = x7*x89*x90*x10;
        const double x99 = x91*u[0]*x[0];
        const double x100 = x22 + 0.502394052142515*x98 - 0.502394052142515*x99 - x97*x48;
        const double x101 = -0.251197026071258*x93 + x22*x[0];
        const double x102 = x49*x101;
        const double x103 = x47*x101;
        const double x104 = x94*x29;
        const double x105 = x104 - x62*x101;
        const double x106 = x66*x105;
        const double x107 = x77*x105;
        const double x108 = x1*x[0] + x106*x[2] + x107*x[2];
        const double x109 = 0.01814406586495*x60;
        const double x110 = x44*x101;
        const double x111 = x73*x43;
        const double x112 = x29*(x19 + 710.350681205318*x98 - 710.350681205318*x99 - x97*x46);
        const double x113 = 0.01814406586495*x104;
        const double x114 = 0.283877349159248*x94;
        const double x115 = x112 - x51*x114 - x62*x100 + x75*x102 - x61*x28*x103;
        const double x116 = x84*x105;
        const double x117 = x87*x105;
        const double x118 = -x108*x111 + x109*x108 - x45*x100 - x51*x96 + x55*x112 + x59*x102 - x72*x110 + x72*x113 + (x1 - x63*x116 - x63*x117 + x78*x115 + x85*x115)*x88 - x50*x95*x103;
        const double x119 = x31*x[2];
        const double x120 = 0.399118713956531*x17;
        const double x121 = x42*x119*x120;
        const double x122 = 1.40595477285733*x17;
        const double x123 = x119*x122;
        const double x124 = x68*x[2];
        const double x125 = 0.702977386428664*x30;
        const double x126 = x17*x125;
        const double x127 = x124*x126 - x60*x123 + x73*x121;
        const double x128 = x41*x42;
        const double x129 = x99*x128;
        const double x130 = 0.0452309522292668*x40;
        const double x131 = pow(x[1], -2);
        const double x132 = pow(x27, -0.716122650840752);
        const double x133 = x25*x132;
        const double x134 = x133*x131;
        const double x135 = x134*x[2];
        const double x136 = 0.283877349159248*x135;
        const double x137 = 6208.24278451573*x92;
        const double x138 = x29*x137;
        const double x139 = 1.24643926465904*x128;
        const double x140 = -x138*u[0] - x21*x136 + x93*x139;
        const double x141 = x64*x140;
        const double x142 = x65*x141;
        const double x143 = x70*x141;
        const double x144 = x142*x[2] + x143*x[2] + x3*x[0];
        const double x145 = 0.01814406586495*x144;
        const double x146 = x40*x134;
        const double x147 = x146*x[2];
        const double x148 = 0.0051506893207128*x147;
        const double x149 = x93*x40;
        const double x150 = 31.9767298105975*x149;
        const double x151 = x98*x128;
        const double x152 = 225.285531975708*x29;
        const double x153 = x40*x152;
        const double x154 = x63*x140;
        const double x155 = 1.24643926465904*x93;
        const double x156 = x50*x47;
        const double x157 = 12416.4855690315*x29;
        const double x158 = 1.60027513899845*x57;
        const double x159 = x93*x49;
        const double x160 = 1762.37950460535*x93;
        const double x161 = x21*x132;
        const double x162 = x73*x161;
        const double x163 = x131*x[2];
        const double x164 = 0.0805863493656817*x163;
        const double x165 = 2.49287852931808*x129 - 2.49287852931808*x151 + x155*x156 - x158*x159 + x164*x162 - x47*x136 + x51*x160 + x98*x157 - x99*x157;
        const double x166 = 0.0051506893207128*x21;
        const double x167 = x166*x135;
        const double x168 = 0.0014621640307068*x40;
        const double x169 = x73*x168;
        const double x170 = x161*x163;
        const double x171 = 112.642765987854*x29;
        const double x172 = x93*x171;
        const double x173 = 0.0226154761146334*x128;
        const double x174 = x93*x173;
        const double x175 = 0.0226154761146334*x149;
        const double x176 = x8*pow(x[0], 3)/pow(acos(-1), 3);
        const double x177 = x128*x176;
        const double x178 = x40*x177;
        const double x179 = pow(u[0], 2)*pow(x[0], 4)/(pow(acos(-1), 4)*pow(x6, 4));
        const double x180 = x50*x179;
        const double x181 = x41*x149;
        const double x182 = x21*pow(x[2], 2);
        const double x183 = x25*pow(x27, -1.71612265084075);
        const double x184 = x183*x182;
        const double x185 = x184/pow(x[1], 4);
        const double x186 = 0.003688525290006*x40;
        const double x187 = x93*x144;
        const double x188 = x21*x133;
        const double x189 = pow(x[1], -3);
        const double x190 = x189*x[2];
        const double x191 = x188*x190;
        const double x192 = x29*x176;
        const double x193 = x93*x41;
        const double x194 = -3.59485148117086*x177 - 15476.3951423131*x180 - 0.203290999793566*x185 + 0.567754698318496*x191 + 17905.1730815698*x192 - 0.707671748678821*x170*x193 + 7.02643509041089*x57*x179 + 3524.75900921071*x93*x135;
        const double x195 = x81*pow(x140, 2);
        const double x196 = x40*x192;
        const double x197 = x21*x135;
        const double x198 = x92*x40;
        const double x199 = 1024.36005944509*x176;
        const double x200 = x41*x101;
        const double x201 = x200*x161;
        const double x202 = x50*x101;
        const double x203 = x50*x94;
        const double x204 = x93*x101;
        const double x205 = -x114*x135 + x201*x164 + x202*x160 + x203*x155 - x204*x158;
        const double x206 = -x138 - 0.205662478668741*x177 + x205 + x29*x199 + x92*x139;
        const double x207 = -x116*x140 - x117*x140;
        const double x208 = x58*x204;
        const double x209 = x104*x145 - x108*x167 - x108*x172 + x108*x174 - x110*x144 + x203*x175 - x94*x148 + x201*x168*x163;
        const double x210 = x93*x17*x119;
        const double x211 = x31*x120;
        const double x212 = x141*x[2];
        const double x213 = -1.75243723322415*x210*x128 + x212*x126 + 8728.50857394696*x29*x210 + x211*x182*x134;
        const double x214 = x26*x161;
        const double x215 = x26*x133;
        const double x216 = x215*x166;
        const double x217 = pow(x37, -1);
        const double x218 = pow(x[2], -3.0/2.0);
        const double x219 = pow(x35, 0.570570570570571);
        const double x220 = pow(x35, 0.785285285285285);
        const double x221 = 0.338114209580443*x218*x220 - 0.334287723847034*x219*x218;
        const double x222 = x217*x221;
        const double x223 = 43.7958479873403*x119;
        const double x224 = x64*x37;
        const double x225 = 12.4326492308276*x224;
        const double x226 = x27*x188;
        const double x227 = 0.256141065324334*pow(x[2], -2.78528528528529) - 0.284605567960525*pow(x[2], -2.57057057057057);
        const double x228 = pow(x17, -1);
        const double x229 = x227*x228;
        const double x230 = x229*u[1];
        const double x231 = 38.7442038438942*x119;
        const double x232 = 10.9986018824902*x226;
        const double x233 = x82*x64;
        const double x234 = -x34 - x39 - x223*x222 + x225*x226 - x230*x231 + x233*x232;
        const double x235 = x27*x162;
        const double x236 = 3.12225394685924*x233;
        const double x237 = x80*x226;
        const double x238 = 16.4979028237353*x82;
        const double x239 = x238*x237;
        const double x240 = 19.3721019219471*x230;
        const double x241 = x27*x133;
        const double x242 = 10.9986018824902*x233;
        const double x243 = x241*x242;
        const double x244 = 3.5293475066741*x224;
        const double x245 = x225*x241;
        const double x246 = 18.6489738462413*x37;
        const double x247 = x237*x246;
        const double x248 = 21.8979239936702*x222;
        const double x249 = x248*x[2];
        const double x250 = x93*x241;
        const double x251 = x27*x161*x193;
        const double x252 = x189*x184;
        const double x253 = 31.9767298105975*x149;
        const double x254 = x21*x183*x131;
        const double x255 = x21*x215;
        const double x256 = pow(x24, -0.567754698318496)*pow(x21, 2);
        const double x257 = x81*pow(x27, -1.4322453016815)*x256*x131;
        const double x258 = x64*x[2];
        const double x259 = x254*x258;
        const double x260 = pow(x[2], -5.0/2.0);
        const double x261 = pow(x[2], -3);
        const double x262 = x27*x201;
        const double x263 = x237*x105;
        const double x264 = x258*x105;
        const double x265 = x106 + x107 - x238*x263 + x264*x240 + x264*x248;
        const double x266 = x214*x200;
        const double x267 = x216*x108 - x234*x110 + x234*x113 - x266*x168 + x96*x215;
        const double x268 = -x211*x226;
        const double x269 = x64*x17;
        const double x270 = 0.72447653635574*x176 - 4.39076688700449*x92;
        const double x271 = (-x137 + x199)*x29;
        const double x272 = x205 + x271 - x62*x270;
        const double x273 = -0.0805863493656817*x266 + x215*x114;
        const double x274 = x203*x101;
        const double x275 = x29*(58.6039311994387*x176 - 710.350681205318*x92);
        const double x276 = 0.0414475093017575*x176 - 0.502394052142515*x92;
        const double x277 = pow(x101, 2);
        const double x278 = -0.567754698318496*x274 + x275 - x62*x276 + x75*x277;
        const double x279 = pow(x105, 2);
		return {std::vector<double>{-x45*x14 - x53*x52 + x54*x55 + x56*x59 + 0.0362881317299*x72*x60 + (x78*x76 - x84*x79 + x85*x76 - x87*x79)*x88 - x73*x72*x74, -x111*x144 + x129*x130 - x130*x151 + x169*x170 + x175*x156 - x47*x148 + x51*x150 - 0.0290354975240298*x58*x159 + x60*x145 - x72*x167 - x72*x172 + x72*x174 + x98*x153 - x99*x153 + (x3 + x78*x165 - x84*x154 + x85*x165 - x87*x154)*x88, -0.0652252220490771*x178 + 324.872639715332*x196 + 0.0452309522292668*x128*x187 - 0.0128400428187929*x170*x181 - x186*x185 - x187*x152 - 0.0103013786414256*x197*x144 - 280.804732814121*x40*x180 + x53*x191 + 0.127488101076211*x58*x179 + 63.953459621195*x93*x147 + (x78*x194 - x83*x195 + x85*x194 - x86*x195)*x88, -x214*x169 + x234*x109 - x234*x111 + x72*x216 + x88*(x67 + x71 - x235*x244 - x236*x235 + x240*x124 + x47*x243 + x47*x245 - x63*x239 - x63*x247 + x68*x249) + x95*x47*x215, -x166*x146 + 0.00642002140939643*x214*x181 - x215*x253 + x216*x144 - x234*x167 - x234*x172 + x234*x174 + x88*(x142 + x143 + x212*x240 - x225*x197 - 68281.9907767311*x233*x250 + 13.7090892426886*x233*x251 + 7.87634793563099*x233*x252 - x239*x140 - x242*x197 - x247*x140 + x249*x141 - 77184.9048797002*x250*x224 + 15.4965421650365*x251*x224 + 8.90330172415347*x252*x224) + x21*x186*x183*x190, 0.0103013786414256*x234*x255 - x254*x186 + x88*(21.9972037649805*x233*x255 + 24.8652984616551*x255*x224 - x32*x230 - 5.29402126001114*x37*x257 - 8.90330172415347*x37*x259 - x38*x222 - 4.68338092028887*x82*x257 - 7.87634793563099*x82*x259 - x217*x223*(0.501431585770551*x219*x260 - 0.507171314370664*x260*x220 + 0.0953673686650698*pow(x35, -0.429429429429429)*x261 - 0.132758056764693*pow(x35, -0.214714714714715)*x261) + 21.8979239936702*pow(x221, 2)*x119/pow(x36, 3.0/2.0) + 12.4326492308276*x64*x222*x226 + x64*x230*x232 + 19.3721019219471*pow(x227, 2)*x119*u[1]/pow(x16, 3.0/2.0) - x69*x228*x119*(-0.713425940205164*pow(x[2], -3.78528528528529) + 0.731598697219849*pow(x[2], -3.57057057057057)))}, {x118, -0.0290354975240298*x208 + x209 + x202*x253 - x45*x270 + x55*x271 + x88*(x207 + x78*x272 + x85*x272), x267 + x88*(x265 - 18.6489738462413*x37*x263 + x78*x273 + x85*x273), x127, x213, x268 - x229*x119*x125 + 0.199559356978266*pow(x27, 0.567754698318496)*x269*x256 - x30*x31*x122}, {0.0362881317299*x108*x104 - x45*x276 - x53*x274 + x55*x275 + x59*x277 + (x78*x278 - x84*x279 + x85*x278 - x87*x279)*x88 - x74*x200*x108, -x104*x123 + x200*x121 + x264*x126}, {}, {}, {}};
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
        const double x1 = pow(x[2], -1);
        const double x2 = sqrt(x1);
        const double x3 = x[0]/acos(-1);
        const double x4 = x3*u[0];
        const double x5 = x4/(0.0825*x4 + 1.44204839461401*x3*x[1]);
        const double x6 = (930.653957571921 + 4305.15564366859*x5)*pow(1 + 3.04481243722737*x5, -0.283877349159248)*pow(x[2]/x[1], 0.283877349159248);
		return 0.0001*(-2472.30109968751*x0 + 0.505613423123327*(-304423.48988703*(-1 + 1.21364877581381*pow(x[1], 0.283877349159248))*sqrt(1 - 3.91255311828556*pow(x[1], 2)*pow(1 + 0.381107883006746*x0, -7.04529616724739)) + 30469.0623257977*sqrt(x6)*x[2]*sqrt(0.425689529793737*pow(x2, 1.57057057057057) - 0.378778912667073*pow(x2, 1.78528528528529))*(1 - 0.791745582846156*pow(x1, 0.214714714714715)))/x[3]);
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        const double x0 = pow(acos(-1), -1);
        const double x1 = x0*u[0];
        const double x2 = 1.44204839461401*x0*x[1];
        const double x3 = 0.0825*x1;
        const double x4 = x2*x[0] + x3*x[0];
        const double x5 = pow(x4, -1);
        const double x6 = 3.04481243722737*x5;
        const double x7 = x1*x6;
        const double x8 = pow(x4, -2);
        const double x9 = (x2 + x3)*x1*x8*x[0];
        const double x10 = x7 - 3.04481243722737*x9;
        const double x11 = pow(x[2], -1);
        const double x12 = sqrt(x11);
        const double x13 = sqrt(0.425689529793737*pow(x12, 1.57057057057057) - 0.378778912667073*pow(x12, 1.78528528528529));
        const double x14 = 4305.15564366859*x5;
        const double x15 = x1*x14;
        const double x16 = 930.653957571921 + x15*x[0];
        const double x17 = 1 + x7*x[0];
        const double x18 = pow(x17, -0.283877349159248);
        const double x19 = x[2]/x[1];
        const double x20 = pow(x19, 0.283877349159248);
        const double x21 = x20*x18;
        const double x22 = x21*x16;
        const double x23 = pow(x22, -1.0/2.0);
        const double x24 = x23*x13;
        const double x25 = x24*x[2];
        const double x26 = x20*pow(x17, -1.28387734915925);
        const double x27 = x26*x16;
        const double x28 = 1 - 0.791745582846156*pow(x11, 0.214714714714715);
        const double x29 = 8649.47664441536*x28;
        const double x30 = x25*x29*x27;
        const double x31 = x21*(x15 - 4305.15564366859*x9);
        const double x32 = 0.283877349159248*x16;
        const double x33 = x32*x26;
        const double x34 = 15234.5311628988*x22*x28*x[2];
        const double x35 = x13/pow(x22, 3.0/2.0);
        const double x36 = x34*x35;
        const double x37 = 30469.0623257977*x28;
        const double x38 = x37*x25;
        const double x39 = 5.05613423123327e-05/x[3];
        const double x40 = x8*u[0]*pow(x[0], 2)/pow(acos(-1), 2);
        const double x41 = x40*x21;
        const double x42 = x25*x28;
        const double x43 = pow(x[1], 2);
        const double x44 = pow(x19, -0.716122650840752)*x18;
        const double x45 = x44/x43;
        const double x46 = x24*x29*x16;
        const double x47 = pow(x[3], 2);
        const double x48 = 1 + 0.381107883006746*x47;
        const double x49 = pow(x48, -7.04529616724739);
        const double x50 = -1 + 1.21364877581381*pow(x[1], 0.283877349159248);
        const double x51 = sqrt(1 - 3.91255311828556*x43*x49);
        const double x52 = x50/x51;
        const double x53 = x40*x27;
        const double x54 = x24*x22;
        const double x55 = x54*x37;
        const double x56 = pow(x[2], -3.0/2.0);
        const double x57 = x0*x[0];
        const double x58 = -0.251197026071258*x40 + x6*x57;
        const double x59 = x21*(-355.175340602659*x40 + x57*x14);
		return {std::vector<double>{x39*(-x30*x10 + x31*x38 - x36*(x31 - x33*x10)), x39*(-x36*(-6208.24278451573*x41 + 1.24643926465904*x53 - x45*x32*x[2]) - 189159336.335093*x41*x42 - 104882.232644228*x51*pow(x[1], -0.716122650840752) + 37977.8356402176*x53*x42 - x45*x46*pow(x[2], 2) + 1191073.07463687*x52*x49*x[1]), x39*(x55 + 5179.72313501168*x54*pow(x[2], -0.214714714714715) + x44*x46*x19 + x34*x23*(-0.334287723847034*x56*pow(x12, 0.570570570570571) + 0.338114209580443*x56*pow(x12, 0.785285285285285))/x13 - 4324.73832220768*x35*x28*pow(x19, 0.567754698318496)*pow(x17, -0.567754698318496)*pow(x16, 2)), 0.0001*(-4944.60219937502*x[3] - 0.505613423123327*(-304423.48988703*x50*x51 + x55*x[2])/x47 - 1616978.28929057*x52*x43*pow(x48, -8.04529616724739))}, {x39*(-x36*(x59 - x58*x33) - x58*x30 + x59*x38)}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        const double x0 = pow(acos(-1), -1);
        const double x1 = x0*u[0];
        const double x2 = 8610.31128733719*x1;
        const double x3 = 1.44204839461401*x0*x[1];
        const double x4 = 0.0825*x1;
        const double x5 = x3*x[0] + x4*x[0];
        const double x6 = pow(x5, -3);
        const double x7 = x3 + x4;
        const double x8 = x6*pow(x7, 2)*x[0];
        const double x9 = pow(x5, -2);
        const double x10 = x7*x9;
        const double x11 = x0/x5;
        const double x12 = 3.04481243722737*x11;
        const double x13 = x12*x[0];
        const double x14 = 1 + x13*u[0];
        const double x15 = pow(x14, -0.283877349159248);
        const double x16 = pow(x[1], -1);
        const double x17 = x16*x[2];
        const double x18 = pow(x17, 0.283877349159248);
        const double x19 = x15*x18;
        const double x20 = x19*(-x2*x10 + x2*x8);
        const double x21 = pow(x[2], -1);
        const double x22 = 1 - 0.791745582846156*pow(x21, 0.214714714714715);
        const double x23 = 30469.0623257977*x22;
        const double x24 = sqrt(x21);
        const double x25 = 0.425689529793737*pow(x24, 1.57057057057057) - 0.378778912667073*pow(x24, 1.78528528528529);
        const double x26 = sqrt(x25);
        const double x27 = 4305.15564366859*x11;
        const double x28 = x27*u[0];
        const double x29 = 930.653957571921 + x28*x[0];
        const double x30 = x29*x19;
        const double x31 = pow(x30, -1.0/2.0);
        const double x32 = x31*x26;
        const double x33 = x32*x[2];
        const double x34 = x33*x23;
        const double x35 = x10*x[0];
        const double x36 = x1*x35;
        const double x37 = x28 - 4305.15564366859*x36;
        const double x38 = x37*x19;
        const double x39 = -3.04481243722737*x36 + x12*u[0];
        const double x40 = x29*x18;
        const double x41 = pow(x14, -1.28387734915925);
        const double x42 = 0.283877349159248*x41;
        const double x43 = x40*x42;
        const double x44 = x38 - x43*x39;
        const double x45 = pow(x30, -3.0/2.0);
        const double x46 = x45*x26;
        const double x47 = x46*x38;
        const double x48 = x23*x[2];
        const double x49 = pow(x39, 2);
        const double x50 = x40*pow(x14, -2.28387734915925);
        const double x51 = 0.36446369852493*x50;
        const double x52 = 6.08962487445473*x1;
        const double x53 = -x52*x10 + x8*x52;
        const double x54 = 0.567754698318496*x41;
        const double x55 = x37*x18;
        const double x56 = x55*x39;
        const double x57 = x46*x30;
        const double x58 = 15234.5311628988*x22;
        const double x59 = x58*x57;
        const double x60 = x59*x[2];
        const double x61 = x26/pow(x30, 5.0/2.0);
        const double x62 = x22*x[2];
        const double x63 = x62*x30;
        const double x64 = 22851.7967443482*x63*x61;
        const double x65 = x33*x22;
        const double x66 = 11104.8671458468*x65*x50;
        const double x67 = 17298.9532888307*x22;
        const double x68 = x41*x33;
        const double x69 = x67*x68;
        const double x70 = x40*x41;
        const double x71 = 8649.47664441536*x22;
        const double x72 = x71*x33;
        const double x73 = x70*x72;
        const double x74 = x41*x39;
        const double x75 = x74*x40;
        const double x76 = x71*x75;
        const double x77 = x44*x46;
        const double x78 = pow(x[3], -1);
        const double x79 = 5.05613423123327e-05*x78;
        const double x80 = x58*x[2];
        const double x81 = pow(x[0], 2);
        const double x82 = pow(acos(-1), -2);
        const double x83 = x9*x82;
        const double x84 = x81*x83;
        const double x85 = x84*u[0];
        const double x86 = -355.175340602659*x85 + x27*x[0];
        const double x87 = x86*x19;
        const double x88 = x87*x46;
        const double x89 = x13 - 0.251197026071258*x85;
        const double x90 = x89*x39;
        const double x91 = x86*x18;
        const double x92 = x91*x39;
        const double x93 = x0*x35;
        const double x94 = 0.502394052142515*u[0];
        const double x95 = x7*x81*x82;
        const double x96 = x6*x95;
        const double x97 = x83*x[0];
        const double x98 = x12 - 3.04481243722737*x93 + x96*x94 - x97*x94;
        const double x99 = x89*x55;
        const double x100 = 710.350681205318*u[0];
        const double x101 = x19*(x27 - 4305.15564366859*x93 + x96*x100 - x97*x100);
        const double x102 = x80*x38;
        const double x103 = x87 - x89*x43;
        const double x104 = x46*x103;
        const double x105 = x64*x44;
        const double x106 = x89*x41;
        const double x107 = x40*x106;
        const double x108 = 4324.73832220768*x62;
        const double x109 = x108*x107;
        const double x110 = x75*x108;
        const double x111 = x71*x68;
        const double x112 = x79*(-x102*x104 + x103*x105 + x104*x110 + x34*x101 - x60*(x101 + x51*x90 - x92*x42 - x98*x43 - x99*x42) + x66*x90 - x73*x98 + x77*x109 - x92*x111 - x99*x111 - x80*x88*x44);
        const double x113 = 0.283877349159248*x15;
        const double x114 = pow(x[1], 2);
        const double x115 = pow(x114, -1);
        const double x116 = pow(x17, -0.716122650840752);
        const double x117 = x116*x[2];
        const double x118 = x115*x117;
        const double x119 = x113*x118;
        const double x120 = 6208.24278451573*x84;
        const double x121 = x19*x120;
        const double x122 = 1.24643926465904*x41;
        const double x123 = x40*x122;
        const double x124 = -x121*u[0] - x29*x119 + x85*x123;
        const double x125 = x46*x124;
        const double x126 = x80*x124;
        const double x127 = x85*x19;
        const double x128 = 94579668.1675466*x127;
        const double x129 = x62*x128;
        const double x130 = x85*x18;
        const double x131 = 53698050.9675289*x65*x130;
        const double x132 = x29*x15;
        const double x133 = x116*x132;
        const double x134 = x115*x133;
        const double x135 = pow(x[2], 2);
        const double x136 = x22*x135;
        const double x137 = 4324.73832220768*x134*x136;
        const double x138 = x97*u[0];
        const double x139 = x65*x138;
        const double x140 = 378318672.670187*x19;
        const double x141 = x32*x15;
        const double x142 = x116*x141;
        const double x143 = x71*x135;
        const double x144 = x115*x143;
        const double x145 = 75955.6712804353*x70;
        const double x146 = x6*u[0];
        const double x147 = x95*x146;
        const double x148 = x65*x147;
        const double x149 = x116*x135;
        const double x150 = x32*x22;
        const double x151 = x115*x149*x150;
        const double x152 = 2455.39050143146*x151;
        const double x153 = x74*x29;
        const double x154 = 37977.8356402176*x85*x68*x22;
        const double x155 = x85*x70;
        const double x156 = 18988.9178201088*x155;
        const double x157 = x62*x156;
        const double x158 = x85*x50;
        const double x159 = x65*x158;
        const double x160 = 48758.8829485682*x159;
        const double x161 = 2.49287852931808*x70;
        const double x162 = x85*x122;
        const double x163 = 12416.4855690315*x19;
        const double x164 = 1.60027513899845*x158;
        const double x165 = 1762.37950460535*x130;
        const double x166 = x29*x118;
        const double x167 = pow(u[0], 2)*pow(x[0], 4)/(pow(acos(-1), 4)*pow(x5, 4));
        const double x168 = x50*x167;
        const double x169 = x41*x18*x167;
        const double x170 = pow(x[1], 3);
        const double x171 = pow(x170, -1);
        const double x172 = x171*x132;
        const double x173 = x146*pow(x[0], 3)/pow(acos(-1), 3);
        const double x174 = x19*x173;
        const double x175 = x85*x41;
        const double x176 = x70*x173;
        const double x177 = pow(x[1], 4);
        const double x178 = pow(x177, -1);
        const double x179 = pow(x17, -1.71612265084075)*x132;
        const double x180 = x179*x135;
        const double x181 = 189159336.335093*x19;
        const double x182 = x85*x181;
        const double x183 = x62*x125;
        const double x184 = pow(x[3], 2);
        const double x185 = 1 + 0.381107883006746*x184;
        const double x186 = pow(x185, -7.04529616724739);
        const double x187 = 1 - 3.91255311828556*x114*x186;
        const double x188 = sqrt(x187);
        const double x189 = x65*x176;
        const double x190 = x85*x142;
        const double x191 = x67*x32;
        const double x192 = pow(x188, -1);
        const double x193 = x186*x192;
        const double x194 = x65*x174;
        const double x195 = x134*x143;
        const double x196 = pow(x[2], 3);
        const double x197 = 6194.0861429839*x22;
        const double x198 = x179*x197;
        const double x199 = -1 + 1.21364877581381*pow(x[1], 0.283877349159248);
        const double x200 = x199/pow(x187, 3.0/2.0);
        const double x201 = 1191073.07463687*x199*x193;
        const double x202 = x29*x175;
        const double x203 = 37977.8356402176*x155;
        const double x204 = 1024.36005944509*x173;
        const double x205 = 0.0805863493656817*x106;
        const double x206 = x86*x113;
        const double x207 = x106*x165 + x205*x166 - x206*x118 - x89*x164 + x91*x162;
        const double x208 = x84*x65;
        const double x209 = x29*x106;
        const double x210 = x86*x142;
        const double x211 = x104*x129 + x104*x137 - x104*x157 + x106*x131 + x109*x125 + x209*x152 - x210*x144 - x88*x126 + x91*x154 + x64*x103*x124;
        const double x212 = pow(x14, -1.5677546983185);
        const double x213 = x22*x17;
        const double x214 = 2455.39050143146*x213;
        const double x215 = x212*x214;
        const double x216 = pow(x29, 2);
        const double x217 = pow(x17, -0.432245301681504);
        const double x218 = x46*x217*x216;
        const double x219 = x71*x17;
        const double x220 = x37*x219;
        const double x221 = pow(x14, -0.567754698318496);
        const double x222 = x217*x221;
        const double x223 = x46*x29*x222;
        const double x224 = x59*x44;
        const double x225 = pow(x[2], -0.214714714714715);
        const double x226 = 1470.40607294595*x32;
        const double x227 = x225*x226;
        const double x228 = 4324.73832220768*x213;
        const double x229 = x228*x133;
        const double x230 = x32*x38;
        const double x231 = x23*x230;
        const double x232 = pow(x26, -1);
        const double x233 = pow(x[2], -3.0/2.0);
        const double x234 = pow(x24, 0.570570570570571);
        const double x235 = pow(x24, 0.785285285285285);
        const double x236 = 0.338114209580443*x233*x235 - 0.334287723847034*x234*x233;
        const double x237 = x232*x236;
        const double x238 = x31*x237;
        const double x239 = x32*x116;
        const double x240 = x214*x239;
        const double x241 = x45*x237;
        const double x242 = 7617.26558144941*x63*x241;
        const double x243 = 2589.86156750584*x225;
        const double x244 = x57*x243;
        const double x245 = 5179.72313501168*x225;
        const double x246 = x216*x222;
        const double x247 = x61*x213*x246;
        const double x248 = 6487.10748331152*x247;
        const double x249 = x46*x216*x221;
        const double x250 = x217*x249;
        const double x251 = 4324.73832220768*x250;
        const double x252 = x62*x115;
        const double x253 = x32*x225;
        const double x254 = 10781.0473083506*x213;
        const double x255 = pow(x17, -1.4322453016815)*x249;
        const double x256 = x57*x124;
        const double x257 = x58*x256;
        const double x258 = x62*x238;
        const double x259 = 53698050.9675289*x213;
        const double x260 = x16*x250;
        const double x261 = x30*x32;
        const double x262 = x30*x31;
        const double x263 = x237*x262;
        const double x264 = x16*x116;
        const double x265 = x29*x264;
        const double x266 = x265*x141;
        const double x267 = x219*x133;
        const double x268 = pow(x[2], -5.0/2.0);
        const double x269 = pow(x196, -1);
        const double x270 = 1470.40607294595*x225;
        const double x271 = x89*x46*x217*x216;
        const double x272 = x86*x223;
        const double x273 = x59*x103;
        const double x274 = x32*x107;
        const double x275 = x87*x32;
        const double x276 = x23*x275;
        const double x277 = -x273 + x276 - x209*x240 + x210*x219 - x229*x104 - x238*x109 - x242*x103 - x244*x103 + x275*x245 - x71*x274 + x80*x87*x238;
        const double x278 = pow(x184, -1);
        const double x279 = 5.05613423123327e-05*x278;
        const double x280 = pow(x185, -8.04529616724739)*x192;
        const double x281 = x280*x199;
        const double x282 = x23*x261;
        const double x283 = -x279*(-x273*x[2] + x276*x[2] - x72*x107);
        const double x284 = 0.72447653635574*x173 - 4.39076688700449*x84;
        const double x285 = (-x120 + x204)*x19;
        const double x286 = pow(x89, 2);
        const double x287 = 0.0414475093017575*x173 - 0.502394052142515*x84;
        const double x288 = x89*x91;
        const double x289 = x19*(58.6039311994387*x173 - 710.350681205318*x84);
		return {std::vector<double>{x79*(x34*x20 - x60*(x20 + x51*x49 - x53*x43 - x54*x56) + x64*pow(x44, 2) + x66*x49 - x69*x56 - x73*x53 - x44*x47*x48 + x77*x76*x[2]), x79*(x105*x124 + x110*x125 - x139*x140 + x139*x145 + x140*x148 - x145*x148 + x152*x153 - x39*x160 - x47*x126 + x55*x154 - x60*(x161*x138 - x161*x147 - x163*x138 + x163*x147 - x37*x119 - x39*x164 + x55*x162 + x74*x165 + 0.0805863493656817*x74*x166) + x74*x131 + x77*x129 + x77*x137 - x77*x157 - x37*x142*x144), x79*(-109531.753831781*x189 + 545553834.576545*x194 + x201 + x125*x195 + x183*x182 + 75108.5424672807*x188*pow(x[1], -1.71612265084075) + 820714.612769848*x193*pow(x[1], 0.283877349159248) - 21562.0946167012*x202*x151 - x203*x183 - x60*(7.02643509041089*x168 - 15476.3951423131*x169 + 17905.1730815698*x174 - 3.59485148117086*x176 + 0.567754698318496*x117*x172 - 0.707671748678821*x166*x175 - 0.203290999793566*x178*x180 + 3524.75900921071*x85*x15*x118) + x64*pow(x124, 2) + 214088.888697901*x65*x168 - 471551248.169811*x65*x169 + 107396101.935058*x115*x190*x136 + x172*x191*x149 + 4660136.67227646*x200*x114*pow(x185, -14.0905923344948) - x32*x178*x196*x198), x79*(-x224 + x231 + x220*x142 - x220*x223 + x230*x245 + x238*x102 - x238*x110 - x240*x153 - x44*x242 - x44*x244 + x44*x248 - x75*x227 - x76*x32 - x77*x229 + x39*x218*x215), x79*(-x257 - x182*x150 + x203*x150 - x229*x125 - x238*x137 - x242*x124 + x248*x124 + x252*x251 - 32156978.7787254*x253*x127 + 6456.21029554136*x253*x155 - x256*x243 - x258*x128 + x258*x156 - x259*x190 - x118*x191*x132 - x226*x134*pow(x[2], 0.785285285285285) + x239*x202*x254 - 1869.34782077622*x255*x171*x136 + x85*x259*x223 + x32*x171*x180*x197 - x85*x212*x218*x254), x79*(x23*x263 + x238*x267 + 641.652570060487*x252*x255 + 4067.56035977644*x261*pow(x[2], -1.21471471471471) + x263*x245 + 2940.81214589189*x266*x225 - x270*x260 + x67*x266 - x71*x260 - x228*x241*x246 - x33*x115*x198 - 7617.26558144941*x62*pow(x236, 2)*x262/pow(x25, 3.0/2.0) + x80*x232*x262*(0.501431585770551*x234*x268 - 0.507171314370664*x235*x268 + 0.0953673686650698*pow(x24, -0.429429429429429)*x269 - 0.132758056764693*pow(x24, -0.214714714714715)*x269) + 1841.5428760736*x61*pow(x29, 3)*pow(x14, -0.851632047477745)*pow(x17, -1.14836795252226)*x252), -x279*(-x224*x[2] + x231*x[2] - x76*x33), 0.0001*(-0.505613423123327*x278*(-104882.232644228*x188*pow(x[1], -0.716122650840752) + x201*x[1] - x257*x[2] - x32*x195 - x65*x182 + x65*x203) - 557093.321481116*x280*pow(x[1], 1.28387734915925) - 3233956.57858113*x281*x[1] - 6326513.44796385*x200*x170*pow(x185, -15.0905923344948)), -x279*(x282 - x213*x251 + x261*x245 + x32*x267 + x80*x263), 0.0001*(-4944.60219937502 + 1.01122684624665*(-304423.48988703*x188*x199 + x282*x[2])/pow(x[3], 3) + 1616978.28929057*x78*x281*x114 + 16986801.899539*x200*x177*pow(x185, -16.0905923344948)*x[3] + 9915717.67082166*x114*pow(x185, -9.04529616724739)*x199*x192*x[3])}, {x112, x79*(x211 + x34*x285 - x60*(x207 + x285 - x43*x284) - x73*x284 - 48758.8829485682*x89*x159), x79*(x277 + 6487.10748331152*x247*x103 - x272*x228 - x274*x270 - (-x205*x265 + x206*x264)*x60 + 1227.69525071573*pow(x14, -1.5677546983185)*x213*x271), x283}, {x79*(x34*x289 - x60*(x289 - x43*x287 + x51*x286 - x54*x288) + x64*pow(x103, 2) + x66*x286 - x69*x288 - x73*x287 - x88*x48*x103 + x71*x104*x107*x[2])}, {}, {}, {}};
	}
private:
	F3dieselMotor(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


std::vector<double> uInitialGuess(double t) {
	 return {0.1, 0.1};
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
            {0.249899415626461, 0.50615, 0.339266666666667, 0.0681},  // x0
            {0.0181818181818182, 0.404445365853659, 0.337037804878049, 0.03},  // lb x
            {1.0, 1.01111341463415, 1.01111341463415, 1.0},  // ub x
            &uInitialGuess,  // u0 initial guesses for optimization
            {0, 0},  // lb u
            {1.0, 1},  // ub u
            {},  // p0 initial guesses for optimization
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
        
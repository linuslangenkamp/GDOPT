
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
		return (0.384928574097715*acos(-1)*x[0] + pow(acos(-1), 2)*(-0.148046527761656 + 9.07026543760135*u[0] + 0.824675324675325*x[1] - 1.23701298701299*x[2]) - 12.9246560741677*pow(x[0], 2))/pow(acos(-1), 3);
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        const double x0 = pow(acos(-1), -1);
		return {std::vector<double>{(-25.8493121483354*x[0] + 0.384928574097715*acos(-1))/pow(acos(-1), 3), 0.824675324675325*x0, -1.23701298701299*x0}, {9.07026543760135*x0}, {}};
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
		return -29.2022182275857*x[0]*x[1]/acos(-1) + 21.1056909960211*sqrt(0.255587584313281 - pow(x[1], 2)*pow(1 + 0.381107883006746*pow(x[3], 2), -7.04529616724739));
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        const double x0 = 29.2022182275857/acos(-1);
        const double x1 = 1 + 0.381107883006746*pow(x[3], 2);
        const double x2 = pow(x1, -7.04529616724739);
        const double x3 = pow(x[1], 2);
        const double x4 = pow(0.255587584313281 - x2*x3, -1.0/2.0);
		return {std::vector<double>{-x0*x[1], -x0*x[0] - 21.1056909960211*x2*x4*x[1], 56.6691582735326*pow(x1, -8.04529616724739)*x4*x3*x[3]}, {}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        const double x0 = pow(x[3], 2);
        const double x1 = 1 + 0.381107883006746*x0;
        const double x2 = pow(x1, -7.04529616724739);
        const double x3 = pow(x[1], 2);
        const double x4 = 0.255587584313281 - x2*x3;
        const double x5 = pow(x4, -1.0/2.0);
        const double x6 = pow(x4, -3.0/2.0);
        const double x7 = pow(x1, -8.04529616724739);
        const double x8 = x3*x5;
		return {std::vector<double>{-29.2022182275857/acos(-1), -21.1056909960211*x2*x5 - 21.1056909960211*pow(x1, -14.0905923344948)*x3*x6, 56.6691582735326*pow(x1, -15.0905923344948)*x6*x[3]*pow(x[1], 3) + 113.338316547065*x5*x7*x[3]*x[1], 56.6691582735326*x8*x7 - 347.509535412495*x0*pow(x1, -9.04529616724739)*x8 - 152.157704764847*x0*pow(x1, -16.0905923344948)*x6*pow(x[1], 4)}, {}, {}, {}, {}, {}};
	}
private:
	F1dieselMotor(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F2dieselMotor : public Expression {
public:
	static std::unique_ptr<F2dieselMotor> create() {
		Adjacency adj{{0, 1, 2}, {0, 1}, {}};
		AdjacencyDiff adjDiff{{{1, 0}, {1, 1}, {2, 0}, {2, 1}, {2, 2}}, {{0, 0}, {0, 1}, {0, 2}, {1, 1}, {1, 2}}, {{0, 0}, {1, 0}}, {}, {}, {}};
		return std::unique_ptr<F2dieselMotor>(new F2dieselMotor(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
        const double x0 = 0.0825*u[0] + 1.44204839461401*x[1];
        const double x1 = pow(x0, -1);
        const double x2 = sqrt(x1*(4381.93459516828*u[0] + 1342.04804545776*x[1])*pow(x1*(1.02709526504533*u[0] + 0.473608284366819*x[1]), -0.283877349159248)*pow(x[2]/x[1], 0.283877349159248));
		return 0.0132270418639513*x2*(x0*x2*x[0] - acos(-1)*x[2]*(38.6336417191289*u[1]*sqrt(-0.791745582846156*pow(x[2], -1.78528528528529) + pow(x[2], -1.57057057057057)) + 66.9337816156513*sqrt(-0.889800866961904*pow(x[2], -0.892642642642643) + pow(x[2], -0.785285285285285))))/acos(-1);
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        const double x0 = pow(acos(-1), -1);
        const double x1 = 0.0132270418639513*x0;
        const double x2 = 4381.93459516828*u[0] + 1342.04804545776*x[1];
        const double x3 = 0.0825*u[0] + 1.44204839461401*x[1];
        const double x4 = pow(x3, -1);
        const double x5 = 1.02709526504533*u[0] + 0.473608284366819*x[1];
        const double x6 = x4*x5;
        const double x7 = pow(x6, -0.283877349159248);
        const double x8 = pow(x[1], -1);
        const double x9 = x8*x[2];
        const double x10 = pow(x9, 0.283877349159248);
        const double x11 = x7*x10;
        const double x12 = x2*x11;
        const double x13 = sqrt(-0.791745582846156*pow(x[2], -1.78528528528529) + pow(x[2], -1.57057057057057));
        const double x14 = sqrt(-0.889800866961904*pow(x[2], -0.892642642642643) + pow(x[2], -0.785285285285285));
        const double x15 = acos(-1)*(66.9337816156513*x14 + 38.6336417191289*x13*u[1]);
        const double x16 = x4*x11;
        const double x17 = sqrt(x2*x16);
        const double x18 = x17*x[0];
        const double x19 = x0*(-x15*x[2] + x3*x18);
        const double x20 = 0.00661352093197565*x19;
        const double x21 = pow(x17, -1);
        const double x22 = x2*x4;
        const double x23 = 0.283877349159248*x22;
        const double x24 = x7*pow(x9, -0.716122650840752);
        const double x25 = pow(x3, -2);
        const double x26 = x5*x25;
        const double x27 = pow(x6, -1.28387734915925)*x23*x10;
        const double x28 = x25*x12;
        const double x29 = x21*(1342.04804545776*x16 - 1.44204839461401*x28 - x27*(-1.44204839461401*x26 + 0.473608284366819*x4) - x24*x23*x[2]/pow(x[1], 2));
        const double x30 = (1.0/2.0)*x3*x[0];
        const double x31 = x1*x17;
        const double x32 = x8*x24*x21;
        const double x33 = x21*(4381.93459516828*x16 - 0.0825*x28 - x27*(-0.0825*x26 + 1.02709526504533*x4));
		return {std::vector<double>{x1*x12, x20*x29 + x31*(1.44204839461401*x18 + x30*x29), x31*(-x15 + 0.141938674579624*x2*x32*x[0] - acos(-1)*x[2]*(33.4668908078257*(0.794274197310588*pow(x[2], -1.89264264264264) - 0.785285285285285*pow(x[2], -1.78528528528529))/x14 + 19.3168208595645*u[1]*(1.41349173874486*pow(x[2], -2.78528528528529) - 1.57057057057057*pow(x[2], -2.57057057057057))/x13)) + 0.00187742879077845*x32*x22*x19}, {x31*(0.0825*x18 + x30*x33) + x33*x20, -0.511008796375814*x13*x17*x[2]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        const double x0 = pow(x[1], -1);
        const double x1 = x0*x[2];
        const double x2 = pow(x1, 0.283877349159248);
        const double x3 = 0.0825*u[0] + 1.44204839461401*x[1];
        const double x4 = pow(x3, -1);
        const double x5 = 1.02709526504533*u[0] + 0.473608284366819*x[1];
        const double x6 = x4*x5;
        const double x7 = pow(x6, -1.28387734915925);
        const double x8 = x2*x7;
        const double x9 = 4381.93459516828*u[0] + 1342.04804545776*x[1];
        const double x10 = pow(x3, -2);
        const double x11 = x5*x10;
        const double x12 = -0.0825*x11 + 1.02709526504533*x4;
        const double x13 = x9*x12;
        const double x14 = pow(acos(-1), -1);
        const double x15 = 0.0037548575815569*x14;
        const double x16 = pow(x6, -0.283877349159248);
        const double x17 = x2*x16;
        const double x18 = (1.0/2.0)*x3;
        const double x19 = x4*x17;
        const double x20 = x9*x19;
        const double x21 = sqrt(x20);
        const double x22 = pow(x21, -1);
        const double x23 = x4*x9;
        const double x24 = pow(x[1], -2);
        const double x25 = pow(x1, -0.716122650840752);
        const double x26 = x25*x16;
        const double x27 = x24*x26;
        const double x28 = x27*x[2];
        const double x29 = -1.44204839461401*x11 + 0.473608284366819*x4;
        const double x30 = x4*x8;
        const double x31 = x30*x29;
        const double x32 = 0.283877349159248*x9;
        const double x33 = x10*x17;
        const double x34 = x9*x33;
        const double x35 = 1342.04804545776*x19 - 1.44204839461401*x34 - 0.283877349159248*x23*x28 - x32*x31;
        const double x36 = x35*x22;
        const double x37 = 1.44204839461401*x21;
        const double x38 = 0.0132270418639513*x14;
        const double x39 = x38*x21;
        const double x40 = 0.00661352093197565*x14;
        const double x41 = x3*x40;
        const double x42 = pow(x3, -3);
        const double x43 = x9*x42*x17;
        const double x44 = x5*x42;
        const double x45 = x30*x32;
        const double x46 = x9*x10;
        const double x47 = x8*x46;
        const double x48 = pow(x[1], -3);
        const double x49 = x23*x[2];
        const double x50 = 0.36446369852493*x2*pow(x6, -2.28387734915925)*x23;
        const double x51 = x7*x25;
        const double x52 = x51*x49*x24;
        const double x53 = pow(x1, -1.71612265084075)*x16;
        const double x54 = x53*x23;
        const double x55 = -761.954083177799*x31 - 3870.59645889446*x33 + 4.15900714481767*x43 - 761.954083177799*x4*x28 - x45*(-1.36593213229413*x10 + 4.15900714481767*x44) + 0.818729751244748*x46*x28 + 0.818729751244748*x47*x29 + x50*pow(x29, 2) + 0.161172698731363*x52*x29 + 0.567754698318496*x48*x49*x26 - 0.203290999793566*x54*pow(x[2], 2)/pow(x[1], 4);
        const double x56 = x3*x[0];
        const double x57 = (1.0/2.0)*x56;
        const double x58 = x57*x22;
        const double x59 = pow(x20, -3.0/2.0);
        const double x60 = x59*pow(x35, 2);
        const double x61 = (1.0/4.0)*x56;
        const double x62 = x36*x[0];
        const double x63 = x37*x[0] + x57*x36;
        const double x64 = -0.791745582846156*pow(x[2], -1.78528528528529) + pow(x[2], -1.57057057057057);
        const double x65 = sqrt(x64);
        const double x66 = 38.6336417191289*x65;
        const double x67 = -0.889800866961904*pow(x[2], -0.892642642642643) + pow(x[2], -0.785285285285285);
        const double x68 = sqrt(x67);
        const double x69 = acos(-1)*(66.9337816156513*x68 + x66*u[1]);
        const double x70 = x56*x21 - x69*x[2];
        const double x71 = x70*x40;
        const double x72 = x71*x22;
        const double x73 = x70*x14;
        const double x74 = 0.00330676046598783*x73;
        const double x75 = x30*x12;
        const double x76 = 4381.93459516828*x19 - 0.0825*x34 - x75*x32;
        const double x77 = x76*x22;
        const double x78 = x77*x40;
        const double x79 = 0.0825*x21;
        const double x80 = x77*x57 + x79*x[0];
        const double x81 = x40*x36;
        const double x82 = x47*x12;
        const double x83 = 0.023419881305638*x46;
        const double x84 = 1243.93197706557*x4;
        const double x85 = 0.0805863493656817*x12;
        const double x86 = -1243.93197706557*x31 - 6429.68071201626*x33 + 0.237937985111311*x43 - 380.977041588899*x75 + 0.409364875622374*x82 - x45*(-1.52019376153453*x10 + 0.237937985111311*x44) + x83*x28 - x84*x28 + x85*x52 + x50*x29*x12 + x8*x83*x29;
        const double x87 = x76*x59;
        const double x88 = x87*x35;
        const double x89 = x77*x[0];
        const double x90 = x39*(0.04125*x62 + 0.721024197307004*x89 + x86*x58 - x88*x61) + x78*x63 + x80*x81 + x86*x72 - x88*x74;
        const double x91 = 0.255504398187907*x[2];
        const double x92 = x65*x91;
        const double x93 = -x92*x36;
        const double x94 = x0*x26;
        const double x95 = x9*x94;
        const double x96 = x22*x[0];
        const double x97 = 0.794274197310588*pow(x[2], -1.89264264264264) - 0.785285285285285*pow(x[2], -1.78528528528529);
        const double x98 = 33.4668908078257/x68;
        const double x99 = pow(x[2], -2.78528528528529);
        const double x100 = 1.41349173874486*x99 - 1.57057057057057*pow(x[2], -2.57057057057057);
        const double x101 = pow(x65, -1);
        const double x102 = 19.3168208595645*x101;
        const double x103 = x102*u[1];
        const double x104 = acos(-1)*(x100*x103 + x98*x97);
        const double x105 = -x69 - x104*x[2] + 0.141938674579624*x96*x95;
        const double x106 = x22*x23;
        const double x107 = x73*x106;
        const double x108 = x94*x22;
        const double x109 = x73*x108;
        const double x110 = x4*x109;
        const double x111 = x9*x96;
        const double x112 = x96*x94;
        const double x113 = x0*x51;
        const double x114 = 0.0402931746828409*x113;
        const double x115 = 0.101645499896783*x111;
        const double x116 = x48*x[2];
        const double x117 = x59*x[0];
        const double x118 = 0.0709693372898121*x95;
        const double x119 = x46*x109;
        const double x120 = x73*x59;
        const double x121 = 0.000938714395389225*x94*x23;
        const double x122 = x23*x108;
        const double x123 = 0.00187742879077845*x14*x122;
        const double x124 = x23*x113;
        const double x125 = x73*x22;
        const double x126 = 0.000532959508361439*x124*x125;
        const double x127 = acos(-1)*x[2];
        const double x128 = x53*x24;
        const double x129 = pow(x1, -1.4322453016815)*pow(x6, -0.567754698318496)*pow(x9, 2)*x24;
        const double x130 = -x87*x118*x[0];
        const double x131 = x78*x105 + x80*x123 - x87*x73*x121;
        const double x132 = -0.0725319112561121*x1*x65*x26*x106;
        const double x133 = x22*(-x83*x94 + x84*x94 - x85*x124);
        const double x134 = -723.019208202766*x33 + 0.0136125*x43 - 2487.86395413115*x75 + 0.046839762611276*x82 - x45*(-0.16947071873248*x10 + 0.0136125*x44) + x50*pow(x12, 2);
        const double x135 = pow(x76, 2)*x59;
		return {std::vector<double>{x39*(x37 + x36*x18) + x41*x35, x39*(1.44204839461401*x62 + x58*x55 - x60*x61) + x72*x55 - x74*x60 + x63*x36*x38, x95*x15, 2.51959963915035*x110 - 0.00270734317374418*x119 - 0.00187742879077845*x27*x107 - x29*x126 + x39*(190.48852079445*x112 - 0.141938674579624*x27*x111 - x29*x111*x114 - x35*x118*x117 + x53*x115*x116) + x63*x123 + x81*x105 - x35*x120*x121 + 0.00134446928241701*x54*x116*x125, -0.00134446928241701*x107*x128 + x39*(-2*x104 - x115*x128 - x127*(x103*(-3.93697774079838*pow(x[2], -3.78528528528529) + 4.03726248771294*pow(x[2], -3.57057057057057)) - 16.7334454039128*pow(x97, 2)/pow(x67, 3.0/2.0) + x98*(1.40195826457088*x99 - 1.50327721578078*pow(x[2], -2.89264264264264)) - 9.65841042978223*pow(x100, 2)*u[1]/pow(x64, 3.0/2.0)) - 0.0201465873414204*x4*x117*x129) - 0.000266479754180719*x10*x120*x129 + x15*x105*x122}, {x39*(x79 + x77*x18) + x76*x41, x90, x131 + x39*(x130 + 0.011709940652819*x23*x112 + x57*x133) + x71*x133, x93, x132 - 0.511008796375814*x65*x21 - x91*x21*x101*x100}, {x39*(0.0825*x89 + x58*x134 - x61*x135) + x72*x134 - x74*x135 + x80*x77*x38, -x77*x92}, {}, {}, {}};
	}
private:
	F2dieselMotor(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F3dieselMotor : public Expression {
public:
	static std::unique_ptr<F3dieselMotor> create() {
		Adjacency adj{{1, 2, 3}, {0}, {}};
		AdjacencyDiff adjDiff{{{1, 1}, {2, 1}, {2, 2}, {3, 1}, {3, 2}, {3, 3}}, {{0, 1}, {0, 2}, {0, 3}}, {{0, 0}}, {}, {}, {}};
		return std::unique_ptr<F3dieselMotor>(new F3dieselMotor(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
        const double x0 = (4381.93459516828*u[0] + 1342.04804545776*x[1])*pow(x[2]/x[1], 0.283877349159248);
        const double x1 = 0.0825*u[0] + 1.44204839461401*x[1];
        const double x2 = pow(1 + 0.381107883006746*pow(x[3], 2), 7.04529616724739);
        const double x3 = pow(x1, -1);
        const double x4 = pow(x3*(1.02709526504533*u[0] + 0.473608284366819*x[1]), 0.283877349159248);
        const double x5 = x3/x4;
        const double x6 = sqrt(x0*x5);
        const double x7 = x4*x6;
		return x5*(x7*pow(x[3], 3)*(-0.020396484072422*u[0] - 0.356517783180682*x[1]) - 0.858199437890784*x0*x[2]*sqrt(-0.889800866961904*pow(x[2], -0.892642642642643) + pow(x[2], -0.785285285285285))*(-1 + 0.791745582846156*pow(x[2], -0.214714714714715)) - 30.4457641592312*x1*x7*(-1 + 1.21364877581381*pow(x[1], 0.283877349159248))*sqrt((0.255587584313281*x2 - pow(x[1], 2))/x2))/(x6*x[3]);
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        const double x0 = 0.0825*u[0] + 1.44204839461401*x[1];
        const double x1 = pow(x0, -1);
        const double x2 = 1.02709526504533*u[0] + 0.473608284366819*x[1];
        const double x3 = x2*x1;
        const double x4 = pow(x3, 0.283877349159248);
        const double x5 = pow(x4, -1);
        const double x6 = x1*x5;
        const double x7 = pow(x[1], 2);
        const double x8 = 4381.93459516828*u[0] + 1342.04804545776*x[1];
        const double x9 = pow(x[1], -1);
        const double x10 = x9*x[2];
        const double x11 = x8*pow(x10, -0.716122650840752);
        const double x12 = x11/x7;
        const double x13 = pow(x0, -2);
        const double x14 = x2*x13;
        const double x15 = 0.473608284366819*x1 - 1.44204839461401*x14;
        const double x16 = 0.283877349159248*x15;
        const double x17 = pow(x10, 0.283877349159248);
        const double x18 = x8*x17;
        const double x19 = x1*pow(x3, -1.28387734915925);
        const double x20 = x19*x18;
        const double x21 = x5*x13;
        const double x22 = x21*x18;
        const double x23 = x6*x17;
        const double x24 = -1.44204839461401*x22 + 1342.04804545776*x23 - x20*x16 - 0.283877349159248*x6*x12*x[2];
        const double x25 = (1.0/2.0)*x24;
        const double x26 = x6*x18;
        const double x27 = pow(x[3], -1);
        const double x28 = sqrt(-0.889800866961904*pow(x[2], -0.892642642642643) + pow(x[2], -0.785285285285285));
        const double x29 = -1 + 0.791745582846156*pow(x[2], -0.214714714714715);
        const double x30 = x28*x29;
        const double x31 = x30*x[2];
        const double x32 = 0.858199437890784*x18;
        const double x33 = -1 + 1.21364877581381*pow(x[1], 0.283877349159248);
        const double x34 = pow(x[3], 2);
        const double x35 = 1 + 0.381107883006746*x34;
        const double x36 = pow(x35, 7.04529616724739);
        const double x37 = pow(x36, -1);
        const double x38 = 0.255587584313281*x36 - x7;
        const double x39 = sqrt(x38*x37);
        const double x40 = sqrt(x26);
        const double x41 = x4*x40;
        const double x42 = x41*x39;
        const double x43 = x42*x33;
        const double x44 = 30.4457641592312*x0;
        const double x45 = -0.020396484072422*u[0] - 0.356517783180682*x[1];
        const double x46 = pow(x[3], 3);
        const double x47 = x41*x46;
        const double x48 = -x32*x31 - x43*x44 + x45*x47;
        const double x49 = x48*x27;
        const double x50 = x49/pow(x26, 3.0/2.0);
        const double x51 = x6*x50;
        const double x52 = x45*x46;
        const double x53 = pow(x3, -0.716122650840752)*x40;
        const double x54 = x53*x52;
        const double x55 = pow(x40, -1);
        const double x56 = x4*x55;
        const double x57 = x52*x56;
        const double x58 = x41*x33/x39;
        const double x59 = x33*x39;
        const double x60 = 15.2228820796156*x0;
        const double x61 = x60*x56*x59;
        const double x62 = 8.6428628226502*x0*x53*x59;
        const double x63 = x31*x17;
        const double x64 = 0.243623381478393*x30;
        const double x65 = x6*x55;
        const double x66 = x65*x27;
        const double x67 = x55*x49;
        const double x68 = x67*x19;
        const double x69 = x67*x21;
        const double x70 = x9*x11;
        const double x71 = x70*x55;
        const double x72 = 1.02709526504533*x1 - 0.0825*x14;
        const double x73 = 0.283877349159248*x72;
        const double x74 = -0.0825*x22 + 4381.93459516828*x23 - x73*x20;
        const double x75 = (1.0/2.0)*x74;
		return {std::vector<double>{-1.44204839461401*x69 - x51*x25 + x66*(-43.904265328616*x43 - 0.356517783180682*x47 - 1151.74487823428*x63 + x54*x16 + x57*x25 - x61*x24 - x62*x15 - 10.4893998842361*x0*x42*pow(x[1], -0.716122650840752) + x64*x12*pow(x[2], 2) + x58*x44*x37*x[1]) - x68*x16, x66*(-x30*x32 - 4.3214314113251*x71*x59 + 0.141938674579624*x1*x71*x52 + 0.145893412648055*x28*x18*pow(x[2], -0.214714714714715) - x64*x11*x10 - 0.429099718945392*x29*x18*x[2]*(0.794274197310588*pow(x[2], -1.89264264264264) - 0.785285285285285*pow(x[2], -1.78528528528529))/x28) - 0.141938674579624*pow(x3, -0.567754698318496)*x70*x50*x13, x66*(3*x41*x45*x34 - x60*x58*(1.37251448160875*pow(x35, -1.0)*x[3] - 5.37003581491038*pow(x35, -8.04529616724739)*x38*x[3])) - x65*x48/x34}, {-0.0825*x69 + x66*(-2.51177554313657*x43 - 0.020396484072422*x47 - 3760.5738064476*x63 - x72*x62 + x73*x54 - x74*x61 + x75*x57) - x73*x68 - x75*x51}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        const double x0 = 0.0825*u[0] + 1.44204839461401*x[1];
        const double x1 = pow(x0, -3);
        const double x2 = pow(x0, -1);
        const double x3 = 1.02709526504533*u[0] + 0.473608284366819*x[1];
        const double x4 = x2*x3;
        const double x5 = pow(x4, 0.283877349159248);
        const double x6 = pow(x5, -1);
        const double x7 = 4381.93459516828*u[0] + 1342.04804545776*x[1];
        const double x8 = pow(x[1], -1);
        const double x9 = x8*x[2];
        const double x10 = pow(x9, 0.283877349159248);
        const double x11 = x6*x10;
        const double x12 = x2*x11;
        const double x13 = x7*x12;
        const double x14 = sqrt(x13);
        const double x15 = pow(x14, -1);
        const double x16 = pow(x[3], -1);
        const double x17 = -1 + 0.791745582846156*pow(x[2], -0.214714714714715);
        const double x18 = -0.889800866961904*pow(x[2], -0.892642642642643) + pow(x[2], -0.785285285285285);
        const double x19 = sqrt(x18);
        const double x20 = x10*x19;
        const double x21 = x20*x17;
        const double x22 = 0.858199437890784*x7;
        const double x23 = x22*x21;
        const double x24 = -1 + 1.21364877581381*pow(x[1], 0.283877349159248);
        const double x25 = pow(x[3], 2);
        const double x26 = 1 + 0.381107883006746*x25;
        const double x27 = pow(x26, 7.04529616724739);
        const double x28 = pow(x27, -1);
        const double x29 = pow(x[1], 2);
        const double x30 = 0.255587584313281*x27 - x29;
        const double x31 = x30*x28;
        const double x32 = sqrt(x31);
        const double x33 = x5*x14;
        const double x34 = x32*x33;
        const double x35 = x34*x24;
        const double x36 = 30.4457641592312*x0;
        const double x37 = -0.020396484072422*u[0] - 0.356517783180682*x[1];
        const double x38 = pow(x[3], 3);
        const double x39 = x33*x38;
        const double x40 = -x23*x[2] - x36*x35 + x37*x39;
        const double x41 = x40*x16;
        const double x42 = x41*x15;
        const double x43 = x6*x42;
        const double x44 = x1*x43;
        const double x45 = x1*x3;
        const double x46 = pow(x0, -2);
        const double x47 = 4.15900714481767*x45 - 1.36593213229413*x46;
        const double x48 = pow(x4, -1.28387734915925);
        const double x49 = x2*x42;
        const double x50 = 0.283877349159248*x49;
        const double x51 = x50*x48;
        const double x52 = x3*x46;
        const double x53 = 0.473608284366819*x2 - 1.44204839461401*x52;
        const double x54 = pow(x53, 2);
        const double x55 = 0.36446369852493*pow(x4, -2.28387734915925);
        const double x56 = x54*x55;
        const double x57 = x6*x46;
        const double x58 = x57*x15;
        const double x59 = pow(x4, -0.716122650840752);
        const double x60 = x38*x14;
        const double x61 = x60*x59;
        const double x62 = x61*x53;
        const double x63 = 0.283877349159248*x37;
        const double x64 = x5*x38;
        const double x65 = 0.283877349159248*x7;
        const double x66 = x2*x6;
        const double x67 = pow(x9, -0.716122650840752);
        const double x68 = pow(x29, -1);
        const double x69 = x67*x68;
        const double x70 = x69*x[2];
        const double x71 = x70*x66;
        const double x72 = x2*x10;
        const double x73 = x53*x48;
        const double x74 = x73*x72;
        const double x75 = x46*x11;
        const double x76 = x7*x75;
        const double x77 = 1342.04804545776*x12 - 1.44204839461401*x76 - x71*x65 - x74*x65;
        const double x78 = x77*x15;
        const double x79 = x78*x64;
        const double x80 = (1.0/2.0)*x37;
        const double x81 = pow(x[1], -0.716122650840752);
        const double x82 = x81*x34;
        const double x83 = 10.4893998842361*x0;
        const double x84 = pow(x32, -1);
        const double x85 = x84*x33;
        const double x86 = x85*x24;
        const double x87 = x86*x28;
        const double x88 = x87*x[1];
        const double x89 = x32*x24;
        const double x90 = x5*x78;
        const double x91 = x89*x90;
        const double x92 = 15.2228820796156*x0;
        const double x93 = x89*x14;
        const double x94 = x53*x59;
        const double x95 = x93*x94;
        const double x96 = 8.6428628226502*x0;
        const double x97 = 1151.74487823428*x21;
        const double x98 = x19*x17;
        const double x99 = pow(x[2], 2);
        const double x100 = x69*x99;
        const double x101 = x98*x100;
        const double x102 = 0.243623381478393*x7;
        const double x103 = -43.904265328616*x35 - 0.356517783180682*x39 + x101*x102 + x63*x62 + x80*x79 - x82*x83 + x88*x36 - x92*x91 - x96*x95 - x97*x[2];
        const double x104 = x16*x103;
        const double x105 = x58*x104;
        const double x106 = pow(x[1], -3);
        const double x107 = x7*x67;
        const double x108 = x107*x106;
        const double x109 = 0.487246762956786*x98;
        const double x110 = pow(x4, -1.71612265084075);
        const double x111 = x54*x110;
        const double x112 = x0*x93;
        const double x113 = 6.18934983540924*x112;
        const double x114 = 8.6428628226502*x59*x112;
        const double x115 = pow(x13, -3.0/2.0);
        const double x116 = pow(x77, 2);
        const double x117 = 7.6114410398078*x0;
        const double x118 = x89*x117;
        const double x119 = x5*x118;
        const double x120 = x78*x94;
        const double x121 = x33*x24/pow(x31, 3.0/2.0);
        const double x122 = x84*x24;
        const double x123 = x28*x[1];
        const double x124 = x123*x122;
        const double x125 = pow(x[1], -1.71612265084075);
        const double x126 = x0*x85;
        const double x127 = x7*x11;
        const double x128 = x1*x127;
        const double x129 = x72*x65;
        const double x130 = x48*x129;
        const double x131 = x7*x10;
        const double x132 = x73*x46;
        const double x133 = 0.818729751244748*x132;
        const double x134 = x7*x72;
        const double x135 = x7*x70;
        const double x136 = x2*x135;
        const double x137 = x7*pow(x9, -1.71612265084075);
        const double x138 = x137/pow(x[1], 4);
        const double x139 = 4.15900714481767*x128 - 761.954083177799*x71 - 761.954083177799*x74 - 3870.59645889446*x75 + x133*x131 - x47*x130 + x56*x134 + 0.818729751244748*x57*x135 + 0.161172698731363*x73*x136 + 0.567754698318496*x66*x108*x[2] - 0.203290999793566*x66*x99*x138;
        const double x140 = x89*x15;
        const double x141 = x5*x92*x140;
        const double x142 = x94*x14;
        const double x143 = x0*x142;
        const double x144 = x63*x61;
        const double x145 = 0.203290999793566*x60*x37;
        const double x146 = (1.0/4.0)*x64*x37;
        const double x147 = x115*x146;
        const double x148 = x63*x38;
        const double x149 = x64*x15;
        const double x150 = x80*x149;
        const double x151 = x0*x32;
        const double x152 = x81*x151;
        const double x153 = 0.174464221751094*x98;
        const double x154 = x66*x15;
        const double x155 = x16*x154;
        const double x156 = x41*x115;
        const double x157 = x57*x156;
        const double x158 = x77*x157;
        const double x159 = x2*x73*x156;
        const double x160 = x66*x115;
        const double x161 = (1.0/2.0)*x41*x160;
        const double x162 = x77*x160;
        const double x163 = x2*x15;
        const double x164 = x73*x163;
        const double x165 = x41/pow(x13, 5.0/2.0);
        const double x166 = (3.0/4.0)*x66*x165;
        const double x167 = 1.02709526504533*x2 - 0.0825*x52;
        const double x168 = x59*x167;
        const double x169 = x60*x168;
        const double x170 = x78*x168;
        const double x171 = x38*x37;
        const double x172 = 0.141938674579624*x171;
        const double x173 = x48*x167;
        const double x174 = x72*x173;
        const double x175 = 4381.93459516828*x12 - 0.0825*x76 - x65*x174;
        const double x176 = x5*x175;
        const double x177 = x15*x176;
        const double x178 = 5.24469994211807*x81;
        const double x179 = 4.3214314113251*x140;
        const double x180 = x110*x145;
        const double x181 = x53*x167;
        const double x182 = x14*x168;
        const double x183 = x89*x182;
        const double x184 = 0.237937985111311*x45 - 1.52019376153453*x46;
        const double x185 = 1067.54172349205*x98;
        const double x186 = x46*x173;
        const double x187 = x186*x131;
        const double x188 = 0.023419881305638*x132;
        const double x189 = x48*x184;
        const double x190 = 1243.93197706557*x66;
        const double x191 = 0.023419881305638*x57;
        const double x192 = x55*x181;
        const double x193 = 0.0805863493656817*x173;
        const double x194 = 0.237937985111311*x128 - 380.977041588899*x174 + 0.409364875622374*x187 - 1243.93197706557*x74 - 6429.68071201626*x75 - x129*x189 + x188*x131 + x191*x135 + x192*x134 + x193*x136 - x70*x190;
        const double x195 = x77*x175;
        const double x196 = x15*x175;
        const double x197 = x175*x149;
        const double x198 = x176*x140;
        const double x199 = x110*x113;
        const double x200 = x77*x115;
        const double x201 = -0.101207323217464*x169 - 12.4634264582718*x183 - 0.178258891590341*x197 - 21.952132664308*x198 - 0.00579009983064797*x62 - 0.010198242036211*x79 + 2.51177554313657*x88 - 1.25588777156829*x91 - 0.713036182868642*x95 + x100*x185 - x114*x184 + x170*x172 - x181*x180 + x181*x199 - 2.97770303340828*x182*x152 + x184*x144 - x194*x141 + x194*x150 - x195*x147 - 4.3214314113251*x0*x89*x170 - x178*x177*x151 + x200*x118*x176 + x92*x124*x177 + x94*x172*x196 + x96*x124*x182 - x0*x94*x179*x175;
        const double x202 = x42*x186;
        const double x203 = x2*x173*x156;
        const double x204 = x160*x175;
        const double x205 = (1.0/2.0)*x204;
        const double x206 = 3760.5738064476*x21;
        const double x207 = -2.51177554313657*x35 - 0.020396484072422*x39 - x206*x[2] + x63*x169 + x80*x197 - x92*x198 - x96*x183;
        const double x208 = x16*x207;
        const double x209 = x58*x208;
        const double x210 = x175*x157;
        const double x211 = (1.0/2.0)*x162;
        const double x212 = 0.283877349159248*x164;
        const double x213 = x163*x173;
        const double x214 = 0.283877349159248*x213;
        const double x215 = -0.0825*x105 + 0.04125*x158 + 0.409364875622374*x202 - 1.44204839461401*x209 + 0.721024197307004*x210 + 0.237937985111311*x44 - x161*x194 + x166*x195 + 0.141938674579624*x175*x159 - x205*x104 - x211*x208 - x212*x208 - x214*x104 + x42*x188 + x49*x192 - x50*x189 + 0.141938674579624*x77*x203;
        const double x216 = x8*x67;
        const double x217 = x7*x216;
        const double x218 = x217*x163;
        const double x219 = pow(x19, -1);
        const double x220 = 0.794274197310588*pow(x[2], -1.89264264264264) - 0.785285285285285*pow(x[2], -1.78528528528529);
        const double x221 = x219*x220;
        const double x222 = x17*x[2];
        const double x223 = x222*x131;
        const double x224 = 0.429099718945392*x223;
        const double x225 = pow(x[2], -0.214714714714715);
        const double x226 = x20*x225;
        const double x227 = 0.145893412648055*x7;
        const double x228 = x9*x67;
        const double x229 = x98*x228;
        const double x230 = -x23 - x217*x179 + x218*x172 - x224*x221 + x227*x226 - x229*x102;
        const double x231 = x16*x230;
        const double x232 = x7*x69;
        const double x233 = x106*x137;
        const double x234 = x233*x[2];
        const double x235 = x163*x171;
        const double x236 = 0.101645499896783*x235;
        const double x237 = x2*x217;
        const double x238 = 0.0709693372898121*x237*x171;
        const double x239 = x216*x140;
        const double x240 = 2.16071570566255*x89*x217;
        const double x241 = x15*x217;
        const double x242 = x46*x241*x171;
        const double x243 = x10*x219*x220;
        const double x244 = x222*x243;
        const double x245 = 190.48852079445*x216;
        const double x246 = x17*x221;
        const double x247 = x38*x218;
        const double x248 = x15*x107;
        const double x249 = 4.3214314113251*x122;
        const double x250 = 3.09467491770462*x140;
        const double x251 = pow(x4, -0.567754698318496);
        const double x252 = x46*x251;
        const double x253 = x252*x156;
        const double x254 = x217*x156;
        const double x255 = x1*x251*x254;
        const double x256 = x217*x252;
        const double x257 = x256*x115;
        const double x258 = 0.141938674579624*x257;
        const double x259 = 0.0805863493656817*pow(x4, -1.5677546983185);
        const double x260 = 0.101645499896783*x253;
        const double x261 = 1.44204839461401*x58;
        const double x262 = 0.212908011869436*x256*x165;
        const double x263 = x68*x137;
        const double x264 = x19*x217;
        const double x265 = pow(x7, 2)*pow(x9, -1.4322453016815)*x68;
        const double x266 = 621.965988532787*x216;
        const double x267 = x115*x175;
        const double x268 = -x206 - 1880.2869032238*x244 - 0.00289504991532399*x247 - x228*x185 - x238*x267 + x267*x240;
        const double x269 = x46*x254*x167;
        const double x270 = 0.0825*x58;
        const double x271 = -x208*x258 - x214*x231 - x231*x205 - x231*x270 + x262*x175;
        const double x272 = pow(x25, -1);
        const double x273 = x40*x272;
        const double x274 = x33*x25;
        const double x275 = 1.37251448160875*pow(x26, -1.0);
        const double x276 = pow(x26, -8.04529616724739);
        const double x277 = 5.37003581491038*x30*x276;
        const double x278 = x275*x[3] - x277*x[3];
        const double x279 = x86*x278;
        const double x280 = 3*x37*x274 - x92*x279;
        const double x281 = x16*x280;
        const double x282 = x33*x[3];
        const double x283 = x278*x249;
        const double x284 = x278*x122;
        const double x285 = x284*x117;
        const double x286 = x37*x25;
        const double x287 = 0.851632047477745*x286;
        const double x288 = (3.0/2.0)*x286;
        const double x289 = x272*x154;
        const double x290 = x155*(-0.0611894522172659*x274 - 1.25588777156829*x279 - x285*x177 + x287*x182 + x288*x177 - x0*x283*x182) + x205*x273 - x205*x281 - x207*x289 + x214*x273 - x214*x281 - x270*x281 + x273*x270;
        const double x291 = x216*x190 - x217*x191 - x237*x193;
        const double x292 = pow(x4, -1.0)*x167;
        const double x293 = 0.0136125*x45 - 0.16947071873248*x46;
        const double x294 = 0.0136125*x1;
        const double x295 = pow(x167, 2);
        const double x296 = x55*x295;
        const double x297 = -2487.86395413115*x174 + 0.046839762611276*x187 - 723.019208202766*x75 - x293*x130 + x294*x127 + x296*x134;
        const double x298 = pow(x175, 2);
        const double x299 = x298*x115;
		return {std::vector<double>{-2.88409678922801*x105 + 1.44204839461401*x158 + 4.15900714481767*x44 - x104*x162 - 0.567754698318496*x104*x164 + x116*x166 + x155*(653.908565881775*x101 - 0.202414646434927*x62 - 0.356517783180682*x79 - 30.2524445270542*x82 + 87.8085306572321*x88 - 43.904265328616*x91 - 24.9268529165436*x95 + x111*x113 - x111*x145 - x116*x147 + x120*x148 + 17.2857256453004*x124*x143 - x139*x141 + x139*x150 - 5.95540606681656*x142*x152 - x47*x114 + x47*x144 + x87*x36 + 7.51169685082786*x0*x34*x125 + x119*x115*x116 + x138*x153*pow(x[2], 3) + 20.9787997684723*x28*x126*pow(x[1], 0.283877349159248) - x89*x96*x120 + x90*x36*x124 - x99*x109*x108 + x36*x29*pow(x26, -14.0905923344948)*x121 - x81*x83*x90*x32) - x161*x139 + x42*x133 - x51*x47 + x56*x49 + 0.283877349159248*x77*x159, 0.409364875622374*x255 + x155*(195.795969289485*x226 - 326.954282940887*x229 - 5799.56857914862*x239 - 0.204682437811187*x242 - 575.872439117138*x244 - 0.0506036616087318*x247 - x97 + x109*x135 + x200*x240 + x232*x179 + x234*x236 - x234*x250 + x235*x245 - x238*x200 - 0.0414158352423263*x19*x232*pow(x[2], 0.785285285285285) - x232*x163*x172 + x28*x248*x249 - 1.48885151670414*x32*x248*x125 + 0.121811690739196*x7*x246*x100 - x99*x233*x153) - x211*x231 - x212*x231 - x231*x261 + 0.141938674579624*x232*x253 - x234*x260 - x253*x245 - x258*x104 + x77*x262 + x53*x46*x254*x259, x155*(-0.487246762956786*x17*x264 - x236*x263 + x263*x250 + 0.0828316704846526*x264*x225 + 0.214549859472696*pow(x220, 2)*x223/pow(x18, 3.0/2.0) - x219*x224*(-1.50327721578078*pow(x[2], -2.89264264264264) + 1.40195826457088*pow(x[2], -2.78528528528529)) - x22*x17*x243 + x227*x225*x243 - x228*x246*x102 + x263*x153*x[2] + 0.114567950172572*x7*x20*pow(x[2], -1.21471471471471) + 0.613378246810239*x89*x265*x160 - 0.0201465873414204*x57*x265*x115*x171) - 0.283877349159248*x231*x257 + x263*x260 + 0.0604397620242613*x1*pow(x4, -0.851632047477745)*x265*x165, x155*(-1.06955334954205*x274 - 21.952132664308*x279 - x283*x143 + x287*x142 - x90*x285 + x90*x288 - x278*x126*x178 - x92*x278*x123*x121 - 163.494843947386*x0*x276*x282*x122*x[1]) + x211*x273 - x211*x281 + x212*x273 - x212*x281 - x261*x281 + x273*x261 - x289*x103, x155*(0.425816023738872*x218*x286 - 2.16071570566255*x284*x241) - x230*x289 + x273*x258 - x281*x258, x155*(6*x37*x282 + pow(x278, 2)*x117*x121 - x86*x92*(x275 - x277 - 8.41660409968614*x25*pow(x26, -2.0) + 32.930410615602*x30*x25*pow(x26, -9.04529616724739))) - 2*x289*x280 + 2*x40*x154/x38}, {x215 + x155*(x201 - 0.865375490449482*x82), 0.011709940652819*x255 + x271 + x155*(639.295392089675*x226 + x268 - 0.356518091434321*x237*x140 - x291*x141 + x291*x150 - 1.22675649362048*x217*x292*x140 + 0.0402931746828409*x218*x292*x171) - x291*x161 + 0.0402931746828409*pow(x4, -1.5677546983185)*x269, x290}, {0.046839762611276*x202 - 0.165*x209 + 0.0825*x210 + x155*(-0.0115801996612959*x169 - 1.42607236573728*x183 - 0.020396484072422*x197 - 2.51177554313657*x198 - x293*x114 + x293*x144 - x295*x180 + x295*x199 - x297*x141 + x297*x150 + x299*x119 - x299*x146 + x168*x196*x148 - x96*x168*x175*x140) + 0.283877349159248*x203*x175 - x208*x204 - 0.567754698318496*x213*x208 - x297*x161 + x298*x166 + x43*x294 + x49*x296 - x51*x293}, {}, {}, {}};
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
            {0.249899415626461, 0.50615, 0.339266666666667, 0.0681},  // x0
            {0.0181818181818182, 0.404445365853659, 0.337037804878049, 0.03},  // lb x
            {1.0, 1.01111341463415, 1.01111341463415, 1.0},  // ub x
            {0, 0},  // u0 initial guesses for optimization
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
        

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
        const double s1 = 0.08249999999999999/M_PI;
		return {std::vector<double>{u[0]*s1}, {x[0]*s1}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        const double s0 = 0.08249999999999999/M_PI;
		return {std::vector<double>{}, {s0}, {}, {}, {}, {}};
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
		return (-12.924656074167725*pow(x[0], 2) + 0.38492857409771514*M_PI*x[0] + pow(M_PI, 2)*(9.0702654376013481*u[0] + 0.82467532467532545*x[1] - 1.2370129870129882*x[2] - 0.14804652776165592))/pow(M_PI, 3);
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        const double s0 = M_1_PI;
		return {std::vector<double>{(-25.84931214833545*x[0] + 0.38492857409771514*M_PI)/pow(M_PI, 3), 0.82467532467532545*s0, -1.2370129870129882*s0}, {9.0702654376013481*s0}, {}};
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
        const double s8 = 29.202218227585707/M_PI;
        const double s9 = 0.38110788300674586*pow(x[3], 2) + 1;
        const double s10 = pow(s9, -7.0452961672473862);
        const double s11 = pow(x[1], 2);
        const double s12 = pow(-s10*s11 + 0.25558758431328077, -1.0/2.0);
		return {std::vector<double>{-x[1]*s8, -x[0]*s8 - 21.105690996021139*x[1]*s10*s12, 56.669158273532567*x[3]*s11*s12*pow(s9, -8.0452961672473862)}, {}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        const double s0 = pow(x[1], 2);
        const double s1 = pow(x[3], 2);
        const double s2 = 0.38110788300674586*s1 + 1;
        const double s3 = pow(s2, -7.0452961672473862);
        const double s4 = -s0*s3 + 0.25558758431328077;
        const double s5 = pow(s4, -1.0/2.0);
        const double s6 = s0/s4;
        const double s7 = pow(s2, -8.0452961672473862);
		return {std::vector<double>{-29.202218227585707/M_PI, -21.105690996021139*s5*(pow(s2, -14.090592334494772)*s6 + s3), x[1]*x[3]*s5*(56.669158273532567*pow(s2, -15.090592334494772)*s6 + 113.33831654706513*s7), s0*s5*(-152.15770476484749*s1*pow(s2, -16.090592334494772)*s6 - 347.50953541249538*s1*pow(s2, -9.0452961672473862) + 56.669158273532567*s7)}, {}, {}, {}, {}, {}};
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
        const double s0 = 0.08249999999999999*u[0] + 1.4420483946140072*x[1];
        const double s1 = 1.0/s0;
        const double s2 = sqrt(s1*pow(x[2]/x[1], 0.28387734915924823)*pow(s1*(1.0270952650453324*u[0] + 0.47360828436681945*x[1]), -0.28387734915924823)*(4381.9345951682772*u[0] + 1342.048045457761*x[1]));
        const double s3 = 1.0/x[2];
		return 0.013227041863951305*s2*(x[0]*s0*s2 - M_PI*x[2]*(38.633641719128939*u[1]*sqrt(pow(s3, 1.5705705705705706) - 0.79174558284615604*pow(s3, 1.7852852852852854)) + 66.933781615651327*sqrt(pow(s3, 0.78528528528528529) - 0.88980086696190397*pow(s3, 0.8926426426426427))))/M_PI;
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        const double s119 = M_1_PI;
        const double s120 = 0.08249999999999999*u[0] + 1.4420483946140072*x[1];
        const double s121 = 4381.9345951682772*u[0] + 1342.048045457761*x[1];
        const double s122 = 1.0/s120;
        const double s123 = 1.0270952650453324*u[0] + 0.47360828436681945*x[1];
        const double s124 = pow(s122*s123, 0.28387734915924823);
        const double s125 = 1.0/s124;
        const double s126 = 1.0/x[1];
        const double s127 = pow(x[2]*s126, 0.28387734915924823);
        const double s128 = s122*s125*s127;
        const double s129 = s121*s128;
        const double s130 = pow(0.057210285249880775*u[0] + x[1], 2);
        const double s131 = 1.0/s130;
        const double s132 = s121*s125*s127;
        const double s133 = s131*s132;
        const double s134 = (1.0/2.0)*s132/s123;
        const double s135 = -0.14193867457962411*s126*s129 + 671.02402272888048*s128 - 0.34672900151442898*s133 + s134*(-0.13444666430591212*s122 + 0.19685701965309813*s123*s131);
        const double s136 = 1.0/s127;
        const double s137 = 1.0/s121;
        const double s138 = s124*s136*s137;
        const double s139 = 1.0/x[2];
        const double s140 = pow(s139, 0.78528528528528529);
        const double s141 = pow(s139, 0.8926426426426427);
        const double s142 = sqrt(s140 - 0.88980086696190397*s141);
        const double s143 = pow(s139, 1.5705705705705706);
        const double s144 = pow(s139, 1.7852852852852854);
        const double s145 = sqrt(s143 - 0.79174558284615604*s144);
        const double s146 = 38.633641719128939*u[1];
        const double s147 = M_PI*(66.933781615651327*s142 + s145*s146);
        const double s148 = sqrt(s129);
        const double s149 = x[0]*s148;
        const double s150 = -x[2]*s147 + s120*s149;
        const double s151 = 0.013227041863951305*s119*s148;
        const double s152 = s120*s150*s151;
        const double s153 = 2.0795035724088358*s130;
        const double s154 = 2190.9672975841386*s128 - 0.019836465081046827*s133 + s134*(-0.29156908117508445*s122 + 0.01126224624779513*s123*s131);
		return {std::vector<double>{0.013227041863951305*s119*s120*s129, s135*s138*s152 + s151*(x[0]*s124*s135*s136*s137*s148*s153 + 1.4420483946140072*s149), 0.001877428790778449*s119*s139*s148*s150 + s151*(0.14193867457962411*x[0]*s120*s139*s148 - M_PI*x[2]*(s146*(-0.78528528528528529*s139*s143 + 0.70674586937243211*s139*s144)/s145 + 66.933781615651327*(-0.39264264264264265*s139*s140 + 0.39713709865529423*s139*s141)/s142) - s147)}, {s138*s152*s154 + s151*(x[0]*s124*s136*s137*s148*s153*s154 + 0.08249999999999999*s149), -0.51100879637581409*x[2]*s145*s148}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        const double s4 = M_1_PI;
        const double s5 = 0.08249999999999999*u[0] + 1.4420483946140072*x[1];
        const double s6 = 1.0/s5;
        const double s7 = 1.0270952650453324*u[0] + 0.47360828436681945*x[1];
        const double s8 = pow(s6*s7, -0.28387734915924823);
        const double s9 = 1.0/x[1];
        const double s10 = pow(x[2]*s9, 0.28387734915924823);
        const double s11 = s10*s8;
        const double s12 = s11*s6;
        const double s13 = 4381.9345951682772*s6;
        const double s14 = 0.057210285249880775*u[0] + x[1];
        const double s15 = pow(s14, 2);
        const double s16 = 1.0/s15;
        const double s17 = s16*(173.84418516897318*u[0] + 53.242978381620041*x[1]);
        const double s18 = s16*(0.01156739979488494*u[0] + 0.0053338931235349012*x[1]);
        const double s19 = s18 - 0.29156908117508445*s6;
        const double s20 = -s19;
        const double s21 = 4381.9345951682772*u[0] + 1342.048045457761*x[1];
        const double s22 = 1.0/s7;
        const double s23 = s21*s22;
        const double s24 = s20*s23;
        const double s25 = -s13 + s17 + s24;
        const double s26 = s15*s25;
        const double s27 = s12*s21;
        const double s28 = 1.0/s21;
        const double s29 = s15*s28;
        const double s30 = 1.0397517862044179*s29;
        const double s31 = s25*s30 - 0.08249999999999999;
        const double s32 = 0.013227041863951305*s31;
        const double s33 = -s4*(0.013752840404243983*s12*s26 + s27*s32);
        const double s34 = s6*s9;
        const double s35 = s16*(0.20219091277663304*u[0] + 0.093233115343469064*x[1]);
        const double s36 = -s35 + 0.13444666430591212*s6;
        const double s37 = s23*s36;
        const double s38 = s16*(3038.6876137684608*u[0] + 930.65395757192096*x[1]) + s34*(1243.9319770655741*u[0] + 380.97704158889945*x[1]) + s37 - 1342.048045457761*s6;
        const double s39 = s38*s5;
        const double s40 = s11*s4;
        const double s41 = s30*s38 - 1.4420483946140072;
        const double s42 = s16*s21;
        const double s43 = pow(u[0] + 0.30626838815384533*x[1], -2);
        const double s44 = s39*s43;
        const double s45 = -s6;
        const double s46 = s15*s43;
        const double s47 = 0.097529056070809*s46*pow(s16*(2.2642167126974879*u[0] + 0.69345800302885796*x[1]) + s34*(0.92689079297550758*u[0] + 0.28387734915924823*x[1]) + 0.00074512980618283953*s37 + s45, 2);
        const double s48 = s28*s39;
        const double s49 = 0.28387734915924823*s9;
        const double s50 = 2.0795035724088358*s29;
        const double s51 = s22*s36;
        const double s52 = s50*s51;
        const double s53 = pow(s14, -3);
        const double s54 = s16*s9;
        const double s55 = 1.0/s14;
        const double s56 = s22*s42;
        const double s57 = s21/pow(u[0] + 0.46111427097847252*x[1], 2);
        const double s58 = 0.44894985564763606*s57;
        const double s59 = s5*s57;
        const double s60 = s5*s56;
        const double s61 = s36*s60;
        const double s62 = -1861.3079151438419*s16 - s22*(-s16*(542.699838602402*u[0] + 250.24664043728129*x[1]) + 360.86776610013015*s6) - 761.95408317779891*s34 + s36*s58 + 0.56775469831849645*s37*s9 + s53*(6077.3752275369216*u[0] + 1861.3079151438419*x[1]) + s54*(1725.2291694392643*u[0] + 528.38315692016079*x[1]) + s56*(-s55*(0.40438182555326607*u[0] + 0.18646623068693813*x[1]) + 0.18646623068693813) + 0.036735026484465166*s59*pow(-s16*s7 + 0.68296606614706623*s6, 2) + 0.69345800302885796*s61 + s6*(1597.056089249372*u[0] + 489.1277942456889*x[1])/pow(x[1], 2);
        const double s63 = s28*s5;
        const double s64 = sqrt(s27);
        const double s65 = 1.0/x[2];
        const double s66 = pow(s65, 0.78528528528528529);
        const double s67 = pow(s65, 0.8926426426426427);
        const double s68 = s66 - 0.88980086696190397*s67;
        const double s69 = sqrt(s68);
        const double s70 = pow(s65, 1.5705705705705706);
        const double s71 = pow(s65, 1.7852852852852854);
        const double s72 = s70 - 0.79174558284615604*s71;
        const double s73 = sqrt(s72);
        const double s74 = 38.633641719128939*u[1];
        const double s75 = M_PI*(66.933781615651327*s69 + s73*s74);
        const double s76 = x[2]*s75;
        const double s77 = s5*s64;
        const double s78 = x[0]*s77;
        const double s79 = -s76 + s78;
        const double s80 = s64*s79;
        const double s81 = 0.0066135209319756524*s80;
        const double s82 = 0.0001453436132944744*s46;
        const double s83 = s28*(0.2379379851113112*u[0] + 4.1590071448176715*x[1]);
        const double s84 = 0.59032396170260715*s29*s9;
        const double s85 = x[0]*s27;
        const double s86 = 0.0066135209319756524*s85;
        const double s87 = x[0]*s11;
        const double s88 = s25*s87;
        const double s89 = s25*s28;
        const double s90 = s77*s79;
        const double s91 = s25*s80;
        const double s92 = -3091.9305921500809*s16 + 0.9479349717200366*s20*s36*s59 + s20*s58 - s22*(-1342.048045457761*s18 + 391.29971550693733*s6) - s22*(-4381.9345951682772*s35 + 589.13648952705228*s6) + s24*s49 - 1243.9319770655741*s34 + s53*(347.68837033794637*u[0] + 106.48595676324008*x[1]) + s54*(49.350426452517603*u[0] + 15.114475564317456*x[1]) + s56*(-s55*(0.02313479958976988*u[0] + 0.010667786247069802*x[1]) + 0.20752480590016797) + 0.039672930162093653*s61;
        const double s93 = s28*s90;
        const double s94 = s4*(4.6224158573505945e-7*s25*s43*s90 - 0.013752840404243983*s29*s51*s91 + 0.0066135209319756524*s31*s38*s87 + 3.5812169024251026e-10*s38*s46*s91 + 0.0066135209319756524*s41*s88 - 0.0095370172427016228*s80*s89 - s86*(-s25*s82 + s25*s83 + s25*s84 - 5.4149929202012645e-8*s26*s44 + 2.0795035724088358*s26*s51*s63 + 0.08249999999999999*s48 - s50*s92) - 0.001877428790778449*s89*s9*s90 + 0.0066135209319756524*s92*s93);
        const double s95 = s48*s64;
        const double s96 = 0.25550439818790704*x[2]*s73;
        const double s97 = s95*s96;
        const double s98 = 0.001877428790778449*s65*s85;
        const double s99 = s16*(431.30729235981607*u[0] + 132.0957892300402*x[1]) + s34*(176.56205609189897*u[0] + 54.075376328394725*x[1]) + 0.14193867457962411*s37 - 190.48852079444973*s6;
        const double s100 = 0.013227041863951305*s65*s93;
        const double s101 = 0.14758099042565179*s29;
        const double s102 = s65*s78;
        const double s103 = 1.0/s69;
        const double s104 = 1.0/s73;
        const double s105 = s104*s74;
        const double s106 = s103*(26.281056895634869*s66 - 26.58188783286684*s67) + s105*(0.78528528528528529*s70 - 0.70674586937243211*s71);
        const double s107 = 0.00093871439538922451*s102 + 0.0066135209319756524*M_PI*s106 - 0.0066135209319756524*s75;
        const double s108 = s16*(24.675213226258801*u[0] + 7.5572377821587278*x[1]) + 0.14193867457962411*s24 - 621.96598853278704*s6;
        const double s109 = s77*s89;
        const double s110 = s4*(0.013227041863951305*x[0]*s10*s21*s6*s65*s8*(s101*s25 - s108*s50 + 0.011709940652818988) - s100*s108 - s107*s109 + 0.00093871439538922451*s25*s28*s5*s64*s65*s79 - s31*s98);
        const double s111 = s64*(s104*(0.40128768844527141*s70 - 0.36115335605158488*s71) - 0.58354070763192611*s73);
        const double s112 = s19*s23;
        const double s113 = s112 + s13 - s17;
        const double s114 = 1.0397517862044179*s46*pow(0.00022820970470500542*s112 - s16*(0.039672930162093653*u[0] + 0.012150564374084497*x[1]) + s6, 2);
        const double s115 = s113*s28;
        const double s116 = 2.0795035724088358*s15*s19*s22;
        const double s117 = 347.68837033794637*s16 + 0.97361952102453075*s19*s57 + 0.039672930162093653*s19*s60 - s22*(s16*(101.37517867469751*u[0] + 46.745541609895533*x[1]) - 2555.2732873650607*s6) - s53*(19.891350845100096*u[0] + 6.0920919615314357*x[1]) + s56*(s55*(0.0013235484837295596*u[0] + 0.00061030709417961857*x[1]) - 0.023134799589769883) - 0.080586349365681736*s59*pow(s16*(0.039672930162093646*u[0] + 0.018293754269273665*x[1]) + s45, 2);
        const double s118 = s113*s63;
		return {std::vector<double>{s40*(-0.013227041863951305*s39 + s6*(83.581171581021223*u[0] + 25.598270700129355*x[1])), s4*(0.013227041863951305*x[0]*s10*s38*s41*s8 - s81*(s28*(s34*(1793.8101105364392*u[0] + 549.38733120806637*x[1]) + 1.4420483946140072*s37 + s42 - 1935.2982294472304*s6) + s38*s52 - 6.9893418421067017e-5*s44 - s47 + s48*s49 - s62*s63) - s86*(-s38*s82 + s38*s83 + s38*s84 + s39*s52 - s47*s5 + 1.4420483946140072*s48 - s50*s62)), 0.003754857581556898*s21*s40*s65, s4*(0.013227041863951305*x[0]*s10*s21*s6*s65*s8*(s101*s38 - s50*s99 + 0.20468243781118695) - s100*s99 - s107*s95 + 0.00093871439538922451*s28*s38*s5*s64*s65*s79 - s41*s98), s4*s64*s65*(-0.0010779895282362903*s102 + 0.030208941309459506*M_PI*s106 - s65*(-0.0016109490365977294*s76 + 0.0016109490365977294*s78) - 0.003754857581556898*s75 + 0.013227041863951305*M_PI*(23.824322939587383*u[1]*pow(s70 - 0.89998613575915831*s71, 2)/pow(s72, 3.0/2.0) - s103*(46.919184157522309*s66 - 50.310014434427416*s67) - s105*(2.0186312438564693*s70 - 1.9684888703991916*s71) + 10.556653810725203*pow(0.98868286033243979*s66 - s67, 2)/pow(s68, 3.0/2.0)))}, {s33, s94, s110, s97, s111}, {s4*(s32*s88 - s81*(0.00022820970470500542*s113*s43*s5 - s114 + s115*s116 + s117*s63 - s28*(0.08249999999999999*s112 - s16*(14.342145276440286*u[0] + 4.392545716483653*x[1]) + 361.50960410138282*s6)) + s86*(-0.00047456289619242429*s113*s46 + s114*s5 + s115*(0.013612499999999998*u[0] + 0.2379379851113112*x[1]) - s116*s118 - s117*s50 + 0.08249999999999999*s118)), s109*s96}, {}, {}, {}};
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
        const double s0 = 0.08249999999999999*u[0] + 1.4420483946140072*x[1];
        const double s1 = 1.0/s0;
        const double s2 = pow(s1*(1.0270952650453324*u[0] + 0.47360828436681945*x[1]), 0.28387734915924823);
        const double s3 = s1/s2;
        const double s4 = pow(x[2]/x[1], 0.28387734915924823)*(4381.9345951682772*u[0] + 1342.048045457761*x[1]);
        const double s5 = sqrt(s3*s4);
        const double s6 = 1.0/x[2];
        const double s7 = pow(0.38110788300674586*pow(x[3], 2) + 1, 7.0452961672473862);
		return s3*(-0.85819943789078446*x[2]*s4*(0.79174558284615615*pow(s6, 0.21471471471471473) - 1)*sqrt(pow(s6, 0.78528528528528529) - 0.88980086696190397*pow(s6, 0.8926426426426427)) + pow(x[3], 3)*s2*s5*(-0.020396484072421957*u[0] - 0.35651778318068189*x[1]) - 30.445764159231203*s0*s2*s5*sqrt((-pow(x[1], 2) + 0.25558758431328077*s7)/s7)*(1.2136487758138137*pow(x[1], 0.28387734915924828) - 1))/(x[3]*s5);
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        const double s257 = pow(0.057210285249880775*u[0] + x[1], 2);
        const double s258 = 1.0/s257;
        const double s259 = 1.0270952650453324*u[0] + 0.47360828436681945*x[1];
        const double s260 = 0.08249999999999999*u[0] + 1.4420483946140072*x[1];
        const double s261 = 1.0/s260;
        const double s262 = s259*s261;
        const double s263 = pow(s262, 0.28387734915924823);
        const double s264 = 1.0/s263;
        const double s265 = s261*s264;
        const double s266 = 1.0/x[1];
        const double s267 = pow(x[2]*s266, 0.28387734915924823);
        const double s268 = 4381.9345951682772*u[0] + 1342.048045457761*x[1];
        const double s269 = s267*s268;
        const double s270 = s265*s269;
        const double s271 = sqrt(s270);
        const double s272 = 1.0/s271;
        const double s273 = 1.0/x[3];
        const double s274 = 1.0/x[2];
        const double s275 = pow(s274, 0.78528528528528529);
        const double s276 = pow(s274, 0.8926426426426427);
        const double s277 = sqrt(s275 - 0.88980086696190397*s276);
        const double s278 = pow(s274, 0.21471471471471473);
        const double s279 = 0.79174558284615615*s278 - 1;
        const double s280 = s277*s279;
        const double s281 = x[2]*s280;
        const double s282 = 0.85819943789078446*s269;
        const double s283 = pow(x[3], 3);
        const double s284 = -0.020396484072421957*u[0] - 0.35651778318068189*x[1];
        const double s285 = 1.2136487758138137*pow(x[1], 0.28387734915924828) - 1;
        const double s286 = pow(x[3], 2);
        const double s287 = 0.38110788300674586*s286 + 1;
        const double s288 = pow(s287, 7.0452961672473862);
        const double s289 = -pow(x[1], 2) + 0.25558758431328077*s288;
        const double s290 = sqrt(s289/s288);
        const double s291 = s263*s271;
        const double s292 = s290*s291;
        const double s293 = s285*s292;
        const double s294 = s260*s293;
        const double s295 = 30.445764159231203*s294;
        const double s296 = s263*s271*s283*s284 - s281*s282 - s295;
        const double s297 = s272*s273*s296;
        const double s298 = s264*s297;
        const double s299 = s258*s298;
        const double s300 = 1.0/s259;
        const double s301 = s258*s259;
        const double s302 = 0.13444666430591212*s261 - 0.19685701965309813*s301;
        const double s303 = -s300*s302;
        const double s304 = s265*s267;
        const double s305 = s264*s269;
        const double s306 = s258*s305;
        const double s307 = (1.0/2.0)*s305;
        const double s308 = -0.14193867457962411*s266*s270 + s303*s307 + 671.02402272888048*s304 - 0.34672900151442898*s306;
        const double s309 = 1.0/s267;
        const double s310 = 1.0/s268;
        const double s311 = s297*s309*s310;
        const double s312 = s267*s281;
        const double s313 = s283*s291;
        const double s314 = 1.0/s289;
        const double s315 = 2.0795035724088358*s257;
        const double s316 = pow(s262, 0.56775469831849645);
        const double s317 = s265*s272;
        const double s318 = s273*s317;
        const double s319 = 0.29156908117508445*s261 - 0.01126224624779513*s301;
        const double s320 = -s300*s319;
        const double s321 = 2190.9672975841386*s304 - 0.019836465081046827*s306 + s307*s320;
		return {std::vector<double>{s298*s303 - 0.69345800302885796*s299 - s308*s311 + s318*(-10.489399884236141*pow(x[1], -0.71612265084075166)*s260*s292 + 30.445764159231203*x[1]*s260*s263*s271*s285*s290*s314 + 0.24362338147839277*x[2]*s266*s267*s268*s277*s279 + s260*s263*s271*s283*s284*s300*s302 + s260*s271*s283*s284*s308*s309*s310*s316 - 30.445764159231203*s263*s271*s285*s290*s300*s302*s315 - 30.445764159231203*s271*s285*s290*s308*s309*s310*s315*s316 - 43.904265328616034*s293 - 1151.7448782342765*s312 - 0.35651778318068189*s313), -0.14193867457962411*s274*s296*s318 + s318*(-x[2]*s279*s282*(-0.39264264264264265*s274*s275 + 0.39713709865529423*s274*s276)/s277 + 0.14193867457962411*s263*s271*s274*s283*s284 + 0.1458934126480553*s267*s268*s277*s278 - 1.1018228193691773*s269*s280 - 4.3214314113251007*s274*s294), s318*(3*s284*s286*s291 - s288*s295*s314*(-2.6850179074551921*x[3]*pow(s287, -8.0452961672473862)*s289 + 0.68625724080437256*x[3]*1.0/s287)) - s296*s317/s286}, {s298*s320 - 0.039672930162093653*s299 - s311*s321 + s318*(s260*s263*s271*s283*s284*s300*s319 + s260*s271*s283*s284*s309*s310*s316*s321 - 30.445764159231203*s263*s271*s285*s290*s300*s315*s319 - 30.445764159231203*s271*s285*s290*s309*s310*s315*s316*s321 - 2.511775543136574*s293 - 3760.5738064475977*s312 - 0.020396484072421957*s313)}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        const double s8 = 0.057210285249880775*u[0] + x[1];
        const double s9 = pow(s8, -3);
        const double s10 = 1.0/x[2];
        const double s11 = pow(s10, 0.21471471471471473);
        const double s12 = 0.79174558284615615*s11 - 1;
        const double s13 = pow(s10, 0.78528528528528529);
        const double s14 = pow(s10, 0.8926426426426427);
        const double s15 = s13 - 0.88980086696190397*s14;
        const double s16 = sqrt(s15);
        const double s17 = 4381.9345951682772*u[0] + 1342.048045457761*x[1];
        const double s18 = 1.0/x[1];
        const double s19 = x[2]*s18;
        const double s20 = pow(s19, 0.28387734915924823);
        const double s21 = s17*s20;
        const double s22 = s16*s21;
        const double s23 = s12*s22;
        const double s24 = x[2]*s23;
        const double s25 = 0.020396484072421957*u[0] + 0.35651778318068189*x[1];
        const double s26 = pow(x[3], 3);
        const double s27 = 1.0270952650453324*u[0] + 0.47360828436681945*x[1];
        const double s28 = 0.08249999999999999*u[0] + 1.4420483946140072*x[1];
        const double s29 = 1.0/s28;
        const double s30 = pow(s27*s29, 0.28387734915924823);
        const double s31 = 1.0/s30;
        const double s32 = sqrt(s21*s29*s31);
        const double s33 = s30*s32;
        const double s34 = s26*s33;
        const double s35 = s25*s34;
        const double s36 = 1.3869160060577159*s35;
        const double s37 = 1.2136487758138137*pow(x[1], 0.28387734915924828) - 1;
        const double s38 = pow(x[3], 2);
        const double s39 = 0.38110788300674586*s38 + 1;
        const double s40 = pow(s39, 7.0452961672473862);
        const double s41 = pow(x[1], 2);
        const double s42 = -0.25558758431328077*s40 + s41;
        const double s43 = sqrt(-s42/s40);
        const double s44 = s33*s43;
        const double s45 = s37*s44;
        const double s46 = s28*s45;
        const double s47 = pow(s8, 2);
        const double s48 = 1.0/s47;
        const double s49 = s28*s48;
        const double s50 = 0.85819943789078446*s24 + s35 + 30.445764159231203*s46;
        const double s51 = 1.0/s27;
        const double s52 = s48*(0.20219091277663304*u[0] + 0.093233115343469064*x[1]);
        const double s53 = 0.13444666430591212*s29 - s52;
        const double s54 = s51*s53;
        const double s55 = s50*s54;
        const double s56 = s49*s55;
        const double s57 = s48*s50;
        const double s58 = s18*s29;
        const double s59 = s17*s54;
        const double s60 = -1342.048045457761*s29 + s48*(3038.6876137684608*u[0] + 930.65395757192096*x[1]) + s58*(1243.9319770655741*u[0] + 380.97704158889945*x[1]) + s59;
        const double s61 = 1.0/s17;
        const double s62 = s28*s61;
        const double s63 = s60*s62;
        const double s64 = 0.40438182555326607*u[0] + 0.18646623068693813*x[1];
        const double s65 = 1.0/s8;
        const double s66 = s51*(-s64*s65 + 0.18646623068693813);
        const double s67 = s49*s66;
        const double s68 = pow(u[0] + 0.46111427097847252*x[1], -2);
        const double s69 = 0.44894985564763606*s68;
        const double s70 = s53*s69;
        const double s71 = pow(-s27*s48 + 0.68296606614706623*s29, 2);
        const double s72 = s47*s68;
        const double s73 = s71*s72;
        const double s74 = s29*s50;
        const double s75 = s16*s20;
        const double s76 = s12*s75;
        const double s77 = x[2]*s76;
        const double s78 = s19*s23;
        const double s79 = pow(x[1], -0.71612265084075166);
        const double s80 = s44*s79;
        const double s81 = s28*s80;
        const double s82 = s28*s54;
        const double s83 = s35*s82;
        const double s84 = 1.0/s42;
        const double s85 = x[1]*s84;
        const double s86 = s46*s85;
        const double s87 = s45*s47;
        const double s88 = s54*s87;
        const double s89 = s60*s61;
        const double s90 = 87.808530657232083*s87;
        const double s91 = 0.71303556636136378*s34;
        const double s92 = 87.808530657232069*s45;
        const double s93 = 60.891528318462406*s46;
        const double s94 = s84*s93;
        const double s95 = 126.62415066767636*s45;
        const double s96 = s47*s95;
        const double s97 = s35*s63;
        const double s98 = 63.31207533383818*s47;
        const double s99 = s45*s98;
        const double s100 = x[1]*s94 + s54*s96 + 2303.489756468553*s77 - 0.48724676295678554*s78 + 20.978799768472282*s81 + 2*s83 - s89*s99 + s91 + s92 - s97;
        const double s101 = s17*s48;
        const double s102 = pow(u[0] + 0.30626838815384533*x[1], -2);
        const double s103 = s102*s28;
        const double s104 = s103*s60;
        const double s105 = 6.9893418421067017e-5*s104;
        const double s106 = -s29;
        const double s107 = s102*s47;
        const double s108 = s107*pow(s106 + s48*(2.2642167126974879*u[0] + 0.69345800302885796*x[1]) + s58*(0.92689079297550758*u[0] + 0.28387734915924823*x[1]) + 0.00074512980618283953*s59, 2);
        const double s109 = 0.097529056070809*s108;
        const double s110 = 0.28387734915924823*s18;
        const double s111 = s110*s63;
        const double s112 = 2.0795035724088358*s47;
        const double s113 = s112*s89;
        const double s114 = s113*s54;
        const double s115 = s18*s48;
        const double s116 = 1.0/s41;
        const double s117 = s28*s68;
        const double s118 = s117*s17;
        const double s119 = s101*s82;
        const double s120 = s101*s66 + s115*(1725.2291694392643*u[0] + 528.38315692016079*x[1]) + s116*s29*(1597.056089249372*u[0] + 489.1277942456889*x[1]) + 0.036735026484465166*s118*s71 + 0.69345800302885796*s119 + s17*s70 + 0.56775469831849645*s18*s59 - 1861.3079151438419*s48 - s51*(360.86776610013015*s29 - s48*(542.699838602402*u[0] + 250.24664043728129*x[1])) - 761.95408317779891*s58 + s9*(6077.3752275369216*u[0] + 1861.3079151438419*x[1]);
        const double s121 = s120*s62;
        const double s122 = (1.0/2.0)*s29;
        const double s123 = s122*s50;
        const double s124 = s19*s76;
        const double s125 = s28*s44;
        const double s126 = 2*s35;
        const double s127 = 0.89789971129527213*s117*s35;
        const double s128 = s45*s85;
        const double s129 = pow(s42, -2);
        const double s130 = s47*s80;
        const double s131 = 56.847894163757815*s45*s72;
        const double s132 = 0.2379379851113112*u[0] + 4.1590071448176715*x[1];
        const double s133 = s132*s45;
        const double s134 = s47*s54;
        const double s135 = 1.4420483946140072*s35;
        const double s136 = s45*s63;
        const double s137 = 30.445764159231203*s133;
        const double s138 = s87*s89;
        const double s139 = 17.972864115540606*s18;
        const double s140 = 0.26889332861182424*s29 - s48*s64;
        const double s141 = s140*s35;
        const double s142 = s51*s98;
        const double s143 = s61*s98;
        const double s144 = 1.0/x[3];
        const double s145 = 1.0/s32;
        const double s146 = s145*s31;
        const double s147 = s144*s146;
        const double s148 = 0.079345860324187306*s35;
        const double s149 = 2.4157453504415072*s46;
        const double s150 = 0.02313479958976988*u[0] + 0.010667786247069802*x[1];
        const double s151 = s51*(-s150*s65 + 0.20752480590016797);
        const double s152 = s48*(0.01156739979488494*u[0] + 0.0053338931235349012*x[1]);
        const double s153 = s152 - 0.29156908117508445*s29;
        const double s154 = -s153;
        const double s155 = s154*s69;
        const double s156 = 0.9479349717200366*s154*s53;
        const double s157 = 4381.9345951682772*s29;
        const double s158 = s48*(173.84418516897318*u[0] + 53.242978381620041*x[1]);
        const double s159 = s154*s51;
        const double s160 = s159*s17;
        const double s161 = -s157 + s158 + s160;
        const double s162 = s50*s61;
        const double s163 = 0.14193867457962411*s161;
        const double s164 = s161*s50;
        const double s165 = s101*s151 + s110*s160 + s115*(49.350426452517603*u[0] + 15.114475564317456*x[1]) + s118*s156 + 0.039672930162093653*s119 + s155*s17 - 3091.9305921500809*s48 - s51*(-1342.048045457761*s152 + 391.29971550693733*s29) - s51*(589.13648952705228*s29 - 4381.9345951682772*s52) - 1243.9319770655741*s58 + s9*(347.68837033794637*u[0] + 106.48595676324008*x[1]);
        const double s166 = s159*s28;
        const double s167 = s161*s62;
        const double s168 = s167*s35;
        const double s169 = s161*s61;
        const double s170 = s169*s87;
        const double s171 = s166*s35;
        const double s172 = s159*s96 - s168 - s169*s99 + 2*s171 + 0.040792968144843914*s34 + 5.0235510862731481*s45 + 7521.1476128951954*s77;
        const double s173 = s112*s51;
        const double s174 = s147*((1.0/2.0)*s100*s154*s51 - 1.0/4.0*s100*s169 + 3.4946709210533508e-5*s102*s161*s50 - 1.3019917330386444e-8*s104*s164 - s117*s156*s50 - s122*(126.62415066767636*x[1]*s154*s30*s32*s37*s43*s47*s51*s84 + 5.0235510862731481*x[1]*s30*s32*s37*s43*s84 + 6.9893418421067017e-5*s102*s161*s25*s26*s28*s30*s32 + 5.4149929202012645e-8*s102*s161*s25*s26*s30*s32*s47*s60 + 1.6486359737235435e-6*s102*s161*s28*s30*s32*s37*s43*s47*s60 + 0.0044250973724140686*s102*s161*s30*s32*s37*s43*s47 - s110*s168 - s113*s159*s35 - 2135.0834469840956*s124 - s126*s151*s49 - s127*s154 - s128*s169*s98 - 21.812744531693884*s130*s169 - s131*s154 + 60.891528318462406*s132*s154*s30*s32*s37*s43*s51 - s135*s169 - 2.511775543136574*s136 - s137*s169 - s139*s170 - s140*s167*s51*s99 - s141*s169*s173 - s151*s95 + 3.9424683202061694*s154*s25*s26*s30*s32*s47*s53*s68 + 2.8840967892280145*s154*s25*s26*s30*s32*s51 + 0.71303556636136378*s154*s26*s28*s30*s32*s51 + 120.03146068223745*s154*s28*s30*s32*s37*s43*s47*s53*s68 + 43.625489063387768*s154*s30*s32*s43*s47*s51*s79 - s159*s63*s99 + s165*s25*s26*s28*s30*s32*s61 + 63.31207533383818*s165*s30*s32*s37*s43*s47*s61 - 0.35651778318068189*s167*s34 + 0.040792968144843914*s26*s28*s30*s32*s51*s53 + 5.0235510862731481*s28*s30*s32*s37*s43*s51*s53 + 1.7307509808989632*s30*s32*s43*s79 - 0.020396484072421957*s34*s63) - s151*s57 + (1.0/2.0)*s154*s28*s50*s51*s60*s61 - s155*s50 - s162*s163*s18 + (1.0/2.0)*s165*s50*s61 + (1.0/2.0)*s172*s51*s53 - 1.0/4.0*s172*s89 + 0.019836465081046827*s28*s48*s50*s60*s61 + (1.0/2.0)*s48*(s159*s90 + s166*s36 - 0.69345800302885796*s168 - 43.904265328616042*s170 + 0.028288210227343277*s34 + 3.4836217044004276*s45 + 5215.6000041235638*s77) + (1.0/2.0)*s48*(-2.5117755431365745*s138 + s148*s82 + s149*s85 + 0.028288210227343281*s34 + 3.4836217044004272*s45 + 91.386188237475025*s77 - 0.019330506798490753*s78 + 0.83229045809914737*s81 + 5.0235510862731489*s88 - 0.039672930162093653*s97) - 0.039672930162093653*s56 - s9*(s148 + s149 + 0.068094572729178238*s24));
        const double s175 = 0.098428509826549063*s35;
        const double s176 = 2.9967311967236836*s46;
        const double s177 = s10*s48;
        const double s178 = 0.14193867457962411*s10*s50;
        const double s179 = s11*s22;
        const double s180 = 0.39264264264264265*s13 - 0.39713709865529423*s14;
        const double s181 = 1.0/s16;
        const double s182 = s12*s181;
        const double s183 = s180*s182;
        const double s184 = s183*s21;
        const double s185 = -190.48852079444973*s29 + s48*(431.30729235981607*u[0] + 132.0957892300402*x[1]) + s58*(176.56205609189897*u[0] + 54.075376328394725*x[1]) + 0.14193867457962411*s59;
        const double s186 = s10*s162;
        const double s187 = s10*s35;
        const double s188 = 0.14193867457962411*s187;
        const double s189 = s10*s46;
        const double s190 = 4.3214314113251007*s189;
        const double s191 = -0.1458934126480553*s179 - 0.85819943789078446*s184 + s188 + s190 + 1.1018228193691773*s23;
        const double s192 = (1.0/2.0)*s191;
        const double s193 = s11*s75;
        const double s194 = s183*s20;
        const double s195 = s10*s34;
        const double s196 = s10*s45;
        const double s197 = s196*s47;
        const double s198 = 8.986432057770303*s197;
        const double s199 = 0.070969337289812057*s187;
        const double s200 = s187*s62;
        const double s201 = 4.4932160288851515*s197;
        const double s202 = s143*s196;
        const double s203 = 0.16208526192104455*s35;
        const double s204 = s144*s29;
        const double s205 = 0.005631123123897565*s35;
        const double s206 = 0.17144384658177853*s46;
        const double s207 = 0.14193867457962411*s160 - 621.96598853278704*s29 + s48*(24.675213226258801*u[0] + 7.5572377821587278*x[1]);
        const double s208 = s147*(0.21290801186943617*s10*s161*s50*s61 + (1.0/2.0)*s10*s29*(17.972864115540609*s159*s87 - s163*s35*s62 - 8.9864320577703047*s170 + 0.28387734915924823*s171 + 0.0057900998306479731*s34 + 0.71303618286864157*s45 + 1067.5417234920478*s77) + s154*s191*s51 - s159*s178 - s169*s192 - s177*(s205 + s206 + 0.0048326266996226884*s24) - s186*s207 - s29*(s159*s198 + s166*s188 + s167*s199 + s169*s201 - 639.29539208967469*s193 - 3760.5738064475977*s194 + 0.0028950499153239865*s195 + 0.35651809143432078*s196 - s200*s207 - s202*s207 + 4828.1155299396451*s76) + s48*(s10*s205 + s10*s206 - 0.0057880191710958087*s179 - 0.034047286364589119*s184 + 0.043712539763834497*s23));
        const double s209 = x[3]*s25;
        const double s210 = 1.0/s39;
        const double s211 = 0.68625724080437256*s210;
        const double s212 = pow(s39, -8.0452961672473862)*s42;
        const double s213 = 2.6850179074551921*s212;
        const double s214 = s211 + s213;
        const double s215 = s37*s43;
        const double s216 = s214*s215;
        const double s217 = s40*s84;
        const double s218 = s216*s217;
        const double s219 = s218*s28;
        const double s220 = 3*s209;
        const double s221 = -30.445764159231203*s219 + s220;
        const double s222 = (1.0/2.0)*s221;
        const double s223 = 1.0/s38;
        const double s224 = s146*s223;
        const double s225 = s224*s57;
        const double s226 = s224*s50;
        const double s227 = pow(s39, 6.0452961672473862);
        const double s228 = s214*s217;
        const double s229 = 126.62415066767636*s215*s228;
        const double s230 = s215*s28;
        const double s231 = -s122*(6*x[3]*s154*s25*s28*s51 + 0.12237890443453174*x[3] - s159*s229*s47 + 63.31207533383818*s161*s214*s37*s40*s43*s47*s61*s84 - s167*s220 - 5.0235510862731481*s218) + (1.0/2.0)*s145*s161*s223*s31*s50*s61 + (1.0/2.0)*s145*s172*s223*s29*s31 + s154*s221*s51 - s159*s226 - s169*s222 - 0.039672930162093653*s225 + s48*(0.11901879048628096*s209 - 1.2078726752207536*s219);
        const double s232 = s146*s50;
        const double s233 = 0.039672930162093653*s48;
        const double s234 = s146*s164*s62;
        const double s235 = s65*(0.0013235484837295596*u[0] + 0.00061030709417961857*x[1]) - 0.023134799589769883;
        const double s236 = s235*s51;
        const double s237 = 0.97361952102453075*s68;
        const double s238 = s153*s28;
        const double s239 = pow(s106 + s48*(0.039672930162093646*u[0] + 0.018293754269273665*x[1]), 2);
        const double s240 = s146*s172;
        const double s241 = s153*s17;
        const double s242 = s241*s51;
        const double s243 = s157 - s158 + s242;
        const double s244 = 0.00022820970470500542*s103*s243;
        const double s245 = pow(0.00022820970470500542*s242 + s29 - s48*(0.039672930162093653*u[0] + 0.012150564374084497*x[1]), 2);
        const double s246 = s153*s173;
        const double s247 = s238*s51;
        const double s248 = s101*s236 + 0.039672930162093653*s101*s247 - 0.080586349365681736*s118*s239 + s237*s241 + 347.68837033794637*s48 - s51*(-2555.2732873650607*s29 + s48*(101.37517867469751*u[0] + 46.745541609895533*x[1])) - s9*(19.891350845100096*u[0] + 6.0920919615314357*x[1]);
        const double s249 = s248*s62;
        const double s250 = s26*s51;
        const double s251 = 0.013612499999999998*u[0] + 0.2379379851113112*x[1];
        const double s252 = s25*s26;
        const double s253 = s215*s243;
        const double s254 = s243*s25*s26*s61;
        const double s255 = s150*s48 - 0.58313816235016891*s29;
        const double s256 = s142*s253*s62;
		return {std::vector<double>{s147*(s100*s54 - 1.0/2.0*s100*s89 + s122*(15.023393701655717*pow(x[1], -1.7161226508407517)*s125 - 41.957599536944564*pow(x[1], 0.28387734915924834)*s125*s84 - s105*s35 - 0.0044250973724140686*s107*s45*s60 - 2.9693466398042871*s108*s46 - s109*s35 + s111*s35 + s113*s141*s51 + s114*s35 - 0.62556508241138242*s116*s24 - s120*s143*s45 - s121*s35 + 1307.8171317635502*s124 + s126*s67 + s127*s53 - 253.24830133535272*s128*s134 - 175.61706131446414*s128 + s129*s41*s93 - 87.250978126775536*s130*s54 + 43.625489063387768*s130*s89 + s131*s53 - 60.891528318462406*s133*s54 + s135*s89 + s136*s140*s142 + s136*s54*s98 + 43.904265328616034*s136 + s137*s89 + s138*s139 - 1.4260711327227276*s34*s82 - 2.8840967892280145*s35*s54 - 0.15278123761395701*s35*s73 - 4.6515415283499992*s46*s73 + s63*s91 + s66*s95 - 60.504889054108318*s80 - s82*s92 + s85*s89*s96 - s94) - s123*(-s105 + s109 + s111 + s114 - s121 + s61*(s101 - 1935.2982294472304*s29 + s58*(1793.8101105364392*u[0] + 549.38733120806637*x[1]) + 1.4420483946140072*s59)) + (1.0/2.0)*s48*(0.98892043987500411*s34 - s36*s63 + 121.78305663692481*s45 + 3194.7468130362263*s77 - 0.67577033444457568*s78 + 29.095833186774112*s81 + 2.7738320121154318*s83 + 84.451435258192191*s86 + 175.61706131446417*s88 - s89*s90) + s55*s63 - 1.3869160060577159*s56 + 0.69345800302885796*s57*s63 - s74*(s28*s70 - s51*(0.19387859642354893*s29 - s48*(0.29156908117508445*u[0] + 0.13444666430591212*x[1])) + s67 + 0.076390618806978505*s73) - s9*(1.1902505368004637*s24 + s36 + 42.225717629096096*s46)), s147*((1.0/2.0)*s10*s29*(-8.9864320577703047*s138 + 0.10120732321746359*s34 + 12.46342645827181*s45 + 326.95428294088754*s77 - 0.069159159727298428*s78 + 2.977703033408281*s81 + 0.28387734915924823*s83 + 8.6428628226502013*s86 + 17.972864115540609*s88 - 0.14193867457962411*s97) + 0.21290801186943617*s10*s50*s60*s61 - s177*(s175 + s176 + 0.08447129180557196*s24) - s178*s54 - s185*s186 + s191*s51*s53 - s192*s89 - s29*(1.4888515167041405*s10*s81 + 0.041415835242326282*s179*s18 + 0.24362338147839277*s18*s184 - 0.31278254120569121*s18*s23 - s185*s200 - s185*s202 + s188*s82 + s190*s85 - 195.79596928948524*s193 - 1151.7448782342765*s194 + 0.050603661608731795*s195 + 6.231713229135905*s196 + s198*s54 + s199*s63 + s201*s89 + 1478.6991611751641*s76) + s48*(s10*s175 + s10*s176 - 0.10117095458998555*s179 - 0.59512526840023183*s184 + 0.76406785201137573*s23)), s10*s146*s204*(s10*s203 - s10*(s203 + 0.13910148067102102*s24 + 4.9348096581353404*s46) - 0.29178682529611061*s11*s180*s181*s21 + 0.13535339178047454*s12*s21*pow(0.98868286033243979*s13 - s14, 2)/pow(s15, 3.0/2.0) + 0.15598378541489827*s179 - 0.85819943789078446*s182*s21*(0.70097913228543862*s13 - 0.75163860789038806*s14) + 1.9600222572599617*s184 + 4.9348096581353396*s189), (1.0/2.0)*s100*s145*s223*s29*s31 - s122*(-60.891528318462406*x[1]*s129*s216*s28*s40 + 83.574504424378645*x[1]*s129*s227*s28*s37*s43 + 6*x[3]*s25*s28*s51*s53 + 2.1391066990840915*x[3] - s134*s229 + 63.31207533383818*s214*s37*s40*s43*s47*s60*s61*s84 - 87.808530657232069*s218 - s220*s63 - 20.978799768472282*s228*s28*s43*s79) + (1.0/2.0)*s145*s223*s31*s50*s60*s61 + s221*s51*s53 - s222*s89 - 0.69345800302885796*s225 - s226*s54 + s48*(2.0803740090865741*s209 - 21.112858814548048*s219), s29*(-s178*s224 + s191*s224), s204*(41.787252212189323*s129*s214*s28*s37*s38*pow(s39, 13.090592334494772)*s43 - 219.49329188766245*s129*s230*s38*pow(s39, 14.090592334494772)*pow(0.25558758431328077*s210 + s212, 2) + 163.4948439473865*s214*s227*s28*s37*s38*s43*s84 - 30.445764159231203*s217*s230*(-s211 - s213 + 16.465205307801018*s38*pow(s39, -9.0452961672473862)*s42 + 4.2083020498430708*s38*pow(s39, -2.0)) - 60.891528318462406*s219 - 2*s226)}, {s174, s208, s231}, {s144*(-s122*(1.0397517862044179*s102*s245*s25*s26*s47 + 31.656037666919087*s102*s245*s28*s37*s43*s47 - 0.01444843001619627*s107*s253 - s143*s215*s248 - 60.891528318462406*s153*s215*s251*s51 - 0.16499999999999998*s153*s25*s250 + 1.9472390420490615*s153*s25*s26*s28*s68 - s153*s256 + 123.28374492320107*s153*s37*s43*s47*s68 - s173*s254*s255 - 5.0235510862731481*s215*s247 + 2*s235*s25*s26*s28*s48*s51 + 126.62415066767636*s235*s37*s43*s51 - 0.081585936289687827*s238*s250 + 0.33515920278664341*s239*s25*s26*s47*s68 + 10.204178043838089*s239*s28*s37*s43*s47*s68 + 0.08249999999999999*s243*s25*s26*s61 + 30.445764159231203*s243*s251*s37*s43*s61 + 0.040792968144843914*s243*s26*s28*s61 + 2.511775543136574*s243*s28*s37*s43*s61 - s244*s252 - s246*s254 - s249*s252 - s255*s256) + s123*s146*(-1.0397517862044179*s107*s245 - s243*s246*s61 - s244 - s249 + s61*(0.08249999999999999*s242 + 361.50960410138282*s29 - s48*(14.342145276440286*u[0] + 4.392545716483653*x[1]))) - 1.0/2.0*s146*s169*s172 - s146*s74*(-s236*s49 - s237*s238 + 0.16757960139332168*s239*s72 + s51*(-0.024054449196944466*s29 + s48*(0.00095431048307800738*u[0] + 0.00044004618269162928*x[1]))) - 0.079345860324187306*s159*s232*s49 + s159*s234 + s159*s240 - 0.0045393993025439529*s232*s9 + s233*s234 + s233*s240)}, {}, {}, {}};
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
        

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


// runtime parameters


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
		return 0.026260565610162729*u[0]*x[0];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{0.026260565610162729*u[0]}, {0.026260565610162729*x[0]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {0.026260565610162729}, {}, {}, {}, {}};
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
		return 2.887155159099656*u[0] - 0.41683999041328129*pow(x[0], 2) + 0.039001418745337815*x[0] + 0.26250230873598346*x[1] - 0.39375346310397524*x[2] - 0.047124673401718109;
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{0.039001418745337815 - 0.83367998082656258*x[0], 0.26250230873598346, -0.39375346310397524}, {2.887155159099656}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-0.83367998082656258}, {}, {}, {}, {}, {}};
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
		return -9.295354760337025*x[0]*x[1] + 21.105690996021139*sqrt(-pow(x[1], 2)*pow(0.38110788300674586*pow(x[3], 2) + 1, -7.0452961672473862) + 0.25558758431328077);
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        const double s8 = 0.38110788300674586*pow(x[3], 2) + 1;
        const double s9 = pow(s8, -7.0452961672473862);
        const double s10 = pow(x[1], 2);
        const double s11 = pow(-s10*s9 + 0.25558758431328077, -1.0/2.0);
		return {std::vector<double>{-9.295354760337025*x[1], -9.295354760337025*x[0] - 21.105690996021139*x[1]*s11*s9, 56.669158273532567*x[3]*s10*s11*pow(s8, -8.0452961672473862)}, {}, {}};
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
		return {std::vector<double>{-9.295354760337025, -21.105690996021139*s5*(pow(s2, -14.090592334494772)*s6 + s3), x[1]*x[3]*s5*(56.669158273532567*pow(s2, -15.090592334494772)*s6 + 113.33831654706513*s7), s0*s5*(-152.15770476484749*s1*pow(s2, -16.090592334494772)*s6 - 347.50953541249538*s1*pow(s2, -9.0452961672473862) + 56.669158273532567*s7)}, {}, {}, {}, {}, {}};
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
        const double s0 = 0.45901826036110277*x[1];
        const double s1 = 0.026260565610162729*u[0] + s0;
        const double s2 = 1.0/s1;
        const double s3 = sqrt(s2*pow(x[2]/x[1], 0.28387734915924823)*pow(s2*(0.99545446595499631*u[0] + s0), -0.28387734915924823)*(1394.8131022528294*u[0] + 427.18716060283873*x[1]));
        const double s4 = 1.0/x[2];
		return s3*(-0.59850022616889709*u[1]*x[2]*sqrt(pow(s4, 1.5705705705705706) - 0.79174558284615604*pow(s4, 1.7852852852852854)) + 0.018144065864949988*x[0]*s1*s3 - 1.0369170922727629*x[2]*sqrt(pow(s4, 0.78528528528528529) - 0.88980086696190397*pow(s4, 0.8926426426426427)));
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        const double s117 = 0.45901826036110277*x[1];
        const double s118 = 0.026260565610162729*u[0] + s117;
        const double s119 = 1.0/s118;
        const double s120 = 0.99545446595499631*u[0] + s117;
        const double s121 = pow(s119*s120, 0.28387734915924823);
        const double s122 = 1.0/s121;
        const double s123 = 1.0/x[1];
        const double s124 = pow(x[2]*s123, 0.28387734915924823);
        const double s125 = 1394.8131022528294*u[0] + 427.18716060283873*x[1];
        const double s126 = s122*s124*s125;
        const double s127 = s119*s126;
        const double s128 = s119*s122*s124;
        const double s129 = pow(0.057210285249880768*u[0] + x[1], 2);
        const double s130 = 1.0/s129;
        const double s131 = s126*s130;
        const double s132 = s120*s130;
        const double s133 = (1.0/2.0)*s126/s120;
        const double s134 = -0.14193867457962411*s123*s127 + 213.59358030141937*s128 - 1.0892812839442543*s131 + s133*(-0.13030488696699949*s119 + 0.61844456674975457*s132);
        const double s135 = 1.0/s124;
        const double s136 = 1.0/s125;
        const double s137 = s121*s135*s136;
        const double s138 = sqrt(s127);
        const double s139 = 1.0/x[2];
        const double s140 = pow(s139, 0.78528528528528529);
        const double s141 = pow(s139, 0.8926426426426427);
        const double s142 = sqrt(s140 - 0.88980086696190397*s141);
        const double s143 = 1.0369170922727629*s142;
        const double s144 = pow(s139, 1.5705705705705706);
        const double s145 = pow(s139, 1.7852852852852854);
        const double s146 = sqrt(s144 - 0.79174558284615604*s145);
        const double s147 = 0.59850022616889709*s146;
        const double s148 = u[1]*s147;
        const double s149 = s138*(0.018144065864949988*x[0]*s118*s138 - x[2]*s143 - x[2]*s148);
        const double s150 = s118*s149;
        const double s151 = x[0]*s138;
        const double s152 = 0.21069776334493312*s129;
        const double s153 = 697.40655112641468*s128 - 0.062318092971807154*s131 + s133*(-0.28258697500403945*s119 + 0.035381390074992386*s132);
		return {std::vector<double>{0.018144065864949988*s118*s127, s134*s137*s150 + s138*(0.018144065864949988*x[0]*s121*s134*s135*s136*s138*s152 + 0.0083284575492066099*s151), s138*(-0.59850022616889709*u[1]*x[2]*(-0.78528528528528529*s139*s144 + 0.70674586937243211*s139*s145)/s146 + 0.0025753446603564022*x[0]*s118*s138*s139 - 1.0369170922727629*x[2]*(-0.39264264264264265*s139*s140 + 0.39713709865529423*s139*s141)/s142 - s143 - s148) + 0.14193867457962411*s139*s149}, {s137*s150*s153 + s138*(0.018144065864949988*x[0]*s121*s135*s136*s138*s152*s153 + 0.0004764734320816331*s151), -x[2]*s138*s147}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        const double s5 = 0.45901826036110277*x[1];
        const double s6 = 0.026260565610162729*u[0] + s5;
        const double s7 = 1.0/s6;
        const double s8 = 0.057210285249880768*u[0] + x[1];
        const double s9 = pow(s8, 2);
        const double s10 = 1.0/s9;
        const double s11 = s10*(0.035220562761846951*u[0] + 0.016240704121380591*x[1]);
        const double s12 = s11 - 0.28258697500403945*s7;
        const double s13 = 1394.8131022528294*u[0] + 427.18716060283873*x[1];
        const double s14 = 0.99545446595499631*u[0] + s5;
        const double s15 = 1.0/s14;
        const double s16 = s13*s15;
        const double s17 = s12*s16;
        const double s18 = -s10*(173.84418516897315*u[0] + 53.242978381620034*x[1]) + s17 + 1394.8131022528294*s7;
        const double s19 = 1.0/x[1];
        const double s20 = pow(x[2]*s19, 0.28387734915924823)*pow(s14*s7, -0.28387734915924823);
        const double s21 = 0.0019114570478640559*s9;
        const double s22 = 1.0/s13;
        const double s23 = s18*s22;
        const double s24 = s21*s23 + 0.0004764734320816331;
        const double s25 = s13*s20;
        const double s26 = s25*s7;
        const double s27 = s24*s26;
        const double s28 = s18*s20*s21*s7 + s27;
        const double s29 = s19*s7;
        const double s30 = 0.28387734915924823*x[1];
        const double s31 = s10*(0.61563340591664595*u[0] + s30);
        const double s32 = s31 - 0.13030488696699949*s7;
        const double s33 = s16*s32;
        const double s34 = -s10*(3038.6876137684608*u[0] + 930.65395757192107*x[1]) - s29*(395.95584604012066*u[0] + 121.2687587467999*x[1]) + s33 + 427.18716060283873*s7;
        const double s35 = 0.018144065864949988*s6;
        const double s36 = s22*s34;
        const double s37 = s21*s36 + 0.0083284575492066099;
        const double s38 = x[0]*s20;
        const double s39 = s34*s38;
        const double s40 = s10*s13;
        const double s41 = pow(u[0] + 0.30626838815384533*x[1], -2);
        const double s42 = s41*s6;
        const double s43 = s41*s9;
        const double s44 = s43*pow(-s10*(7.1132465907456588*u[0] + 2.1785625678885085*x[1]) - s29*(0.92689079297550736*u[0] + s30) + 0.0023408943250747944*s33 + s7, 2);
        const double s45 = s22*s6;
        const double s46 = s34*s45;
        const double s47 = 0.28387734915924823*s19;
        const double s48 = s36*s9;
        const double s49 = 0.21069776334493312*s15;
        const double s50 = pow(s8, -3);
        const double s51 = s10*s19;
        const double s52 = 1.0/s8;
        const double s53 = s15*s40;
        const double s54 = s13/pow(u[0] + 0.46111427097847246*x[1], 2);
        const double s55 = 0.46321985258873616*s54;
        const double s56 = s54*s6;
        const double s57 = s32*s6;
        const double s58 = s53*s57;
        const double s59 = -1861.3079151438421*s10 + s15*(s10*(525.98137329157373*u[0] + 242.53751749359981*x[1]) - 111.32914935221271*s7) - 0.56775469831849645*s19*s33 - 242.53751749359981*s29 - s32*s55 + s50*(6077.3752275369216*u[0] + 1861.3079151438421*x[1]) + s51*(1725.2291694392643*u[0] + 528.38315692016079*x[1]) + s53*(-s52*(1.2312668118332919*u[0] + 0.56775469831849645*x[1]) + 0.56775469831849645) + 0.38597462881247235*s56*pow(s10*s14 - 0.21069776334493312*s7, 2) - 2.1785625678885085*s58 + s7*(508.35874199809757*u[0] + 155.69421251567385*x[1])/pow(x[1], 2);
        const double s60 = sqrt(s26);
        const double s61 = 1.0/x[2];
        const double s62 = pow(s61, 0.78528528528528529);
        const double s63 = pow(s61, 0.8926426426426427);
        const double s64 = s62 - 0.88980086696190397*s63;
        const double s65 = sqrt(s64);
        const double s66 = 1.0369170922727629*s65;
        const double s67 = pow(s61, 1.5705705705705706);
        const double s68 = pow(s61, 1.7852852852852854);
        const double s69 = s67 - 0.79174558284615604*s68;
        const double s70 = sqrt(s69);
        const double s71 = u[1]*s70;
        const double s72 = 0.59850022616889709*s71;
        const double s73 = x[0]*s60;
        const double s74 = x[2]*s66 + x[2]*s72 - s35*s73;
        const double s75 = s60*s74;
        const double s76 = (1.0/2.0)*s75;
        const double s77 = 4.1971133489432667e-7*s43;
        const double s78 = 0.00021871000590239524*u[0] + 0.0038229140957281118*x[1];
        const double s79 = 0.00054261935977941041*s19;
        const double s80 = s15*s21;
        const double s81 = s32*s80;
        const double s82 = s21*s22;
        const double s83 = x[0]*s26;
        const double s84 = s18*s38;
        const double s85 = s23*s75;
        const double s86 = s18*s42;
        const double s87 = s18*s45;
        const double s88 = s18*s43;
        const double s89 = s34*s88;
        const double s90 = s12*s54;
        const double s91 = -3091.9305921500804*s10 - s12*s55 + s15*(427.18716060283873*s11 - 120.71752747532098*s7) + s15*(1394.8131022528294*s31 - 181.75096362914482*s7) - s17*s47 - 395.95584604012066*s29 + s50*(347.68837033794631*u[0] + 106.48595676324007*x[1]) + s51*(49.350426452517596*u[0] + 15.114475564317457*x[1]) + s53*(-s52*(0.070441125523693901*u[0] + 0.032481408242761182*x[1]) + 0.63187411003802652) + 1.0091534315526534*s57*s90 - 0.12463618594361431*s58;
        const double s92 = s23*s9;
        const double s93 = 0.10534888167246656*s15*s32*s85*s9 - 0.14193867457962411*s19*s75*s87 + (1.0/2.0)*s24*s39 + (1.0/2.0)*s37*s84 - s45*s76*s91 + 0.00010978832492295084*s75*s86 - 2.7074964601006316e-8*s75*s89 + s83*(-s18*s77 + s23*s78 + 0.00023823671604081655*s46 + 4.9124994101184804e-10*s6*s89 + s79*s92 - s81*s87 + s82*s91) - 0.22950913018055139*s85;
        const double s94 = s46*s60;
        const double s95 = 0.29925011308444854*x[2]*s70;
        const double s96 = -s94*s95;
        const double s97 = s61*s83;
        const double s98 = s60*s61;
        const double s99 = s74*s98;
        const double s100 = 0.070969337289812057*s99;
        const double s101 = -s10*(431.30729235981607*u[0] + 132.0957892300402*x[1]) - s29*(56.201447978988433*u[0] + 17.212726884436979*x[1]) + 0.14193867457962411*s33 + 60.634379373399952*s7;
        const double s102 = s45*s99;
        const double s103 = 0.0002713096798897052*s9;
        const double s104 = 0.0038229140957281118*s22*s9;
        const double s105 = 1.0/s65;
        const double s106 = s105*(0.40713786731130253*s62 - 0.41179824557128902*s63);
        const double s107 = 1.0/s70;
        const double s108 = u[1]*s107*(0.78528528528528529*s67 - 0.70674586937243211*s68);
        const double s109 = x[0]*s6*s98;
        const double s110 = (1.0/2.0)*s106 + 0.29925011308444854*s108 + 0.0012876723301782011*s109 - 1.0/2.0*s66 - 1.0/2.0*s72;
        const double s111 = -s10*(24.675213226258798*u[0] + 7.5572377821587287*x[1]) + 0.14193867457962411*s17 + 197.97792302006033*s7;
        const double s112 = s60*s87;
        const double s113 = 0.14193867457962411*x[0]*s27*s61 + s100*s87 - s102*s111 + s110*s112 + s97*(-s103*s23 + s104*s111 + 6.7630007422071558e-5);
        const double s114 = s60*(s107*(0.46999342085035012*s67 - 0.42298756266333443*s68) - 0.68345055500691565*s70);
        const double s115 = s43*pow(-s10*(0.12463618594361431*u[0] + 0.038172123774593711*x[1]) + 0.00071694193177914095*s17 + s7, 2);
        const double s116 = 347.68837033794631*s10 + 0.12463618594361431*s12*s53*s6 - s15*(s10*(98.252204817884447*u[0] + 45.305493796626351*x[1]) - 788.31203052325407*s7) - s50*(19.891350845100092*u[0] + 6.0920919615314348*x[1]) + s53*(s52*(0.0040299568845331849*u[0] + 0.0018582706308861959*x[1]) - 0.070441125523693901) - 0.080586349365681723*s56*pow(s10*(0.12463618594361432*u[0] + 0.057471524018927057*x[1]) - s7, 2) + 1.0045662902728985*s90;
		return {std::vector<double>{s20*(s34*s35 + s7*(11.616641711189867*u[0] + 3.5578101326468485*x[1])), s37*s39 - s76*(s22*(-s29*(181.75096362914482*u[0] + 55.664574676106355*x[1]) + 0.45901826036110277*s33 - s40 + 196.08670730851406*s7) - s32*s48*s49 - 0.00021957664984590167*s34*s42 + 0.0098817594006143008*s44 + s45*s59 + s46*s47) + s83*(-s34*s77 + s36*s78 + 8.9647646713167302e-5*s44*s6 - s46*s81 + 0.004164228774603305*s46 + s48*s79 + s59*s82), 0.0051506893207128044*s25*s61, s100*s46 - s101*s102 + s110*s94 + 0.14193867457962411*s37*s97 + s97*(s101*s104 - s103*s36 + 0.0011821302258270508), s98*(-0.59850022616889709*u[1]*s107*(2.0186312438564693*s67 - 1.9684888703991916*s68) + 0.36907891757467437*u[1]*pow(s67 - 0.89998613575915831*s68, 2)/pow(s69, 3.0/2.0) - s105*(0.72685724359330139*s62 - 0.77938691973364849*s63) + s105*(0.81427573462260505*s62 - 0.82359649114257805*s63) + 0.28387734915924823*s106 + 1.3669011100138313*s108 - 0.0014787216373263021*s109 + s61*(0.12628829696086882*x[2]*s65 + 0.072892591757646938*x[2]*s71 - 0.0022098036526797022*s6*s73) - 0.2943572754523075*s65 - 0.16990065767603704*s71 + 0.1635403604775221*pow(0.98868286033243979*s62 - s63, 2)/pow(s64, 3.0/2.0))}, {s28, s93, s113, s96, s114}, {s24*s84 + s76*(-0.10534888167246656*s115 + s116*s45 + s12*s49*s92 - s22*(-s10*(4.5652466305750981*u[0] + 1.3981907270710088*x[1]) + 0.026260565610162729*s17 + 36.628580985625042*s7) + 0.00071694193177914106*s86) + s83*(0.00095572852393202817*s115*s6 - s116*s82 - s12*s80*s87 + 0.0090720329324749938*s23*(0.0013792346123313226*u[0] + 0.024108158284950991*x[1]) + 0.00023823671604081655*s87 - 1.3704037084085102e-6*s88), -s112*s95}, {}, {}, {}};
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
        const double s0 = pow(x[3], 2);
        const double s1 = u[0]*x[0];
        const double s2 = s1/(0.45901826036110277*x[0]*x[1] + 0.026260565610162729*s1);
        const double s3 = pow(x[2]/x[1], 0.28387734915924823)*pow(0.96919390034483355*s2 + 1, -0.28387734915924823)*(1370.3736029396543*s2 + 930.65395757192096);
        const double s4 = 1.0/x[2];
		return -0.24723010996875103*s0 - 5.0561342312332735e-5*(19879.505327143728*x[2]*sqrt(s3)*(0.79174558284615615*pow(s4, 0.21471471471471473) - 1)*sqrt(pow(s4, 0.78528528528528529) - 0.88980086696190397*pow(s4, 0.8926426426426427)) + 602154.98178743932*(1.2136487758138137*pow(x[1], 0.28387734915924828) - 1)*sqrt(-pow(x[1], 2)*pow(0.38110788300674586*s0 + 1, -7.0452961672473862) + 0.25558758431328077))/x[3];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        const double s163 = 0.026260565610162729*u[0];
        const double s164 = x[0]*x[1];
        const double s165 = 1.0/(x[0]*s163 + 0.45901826036110277*s164);
        const double s166 = u[0]*s165;
        const double s167 = 1370.3736029396543*s166;
        const double s168 = u[0]*x[0];
        const double s169 = pow(s164 + 0.057210285249880768*s168, -2);
        const double s170 = s168*s169*(-0.45901826036110277*x[1] - s163);
        const double s171 = s167 + 6503.9779311573275*s170;
        const double s172 = 1.0/x[1];
        const double s173 = x[2]*s172;
        const double s174 = pow(s173, 0.28387734915924823);
        const double s175 = 0.96919390034483355*x[0]*s166 + 1;
        const double s176 = s174*pow(s175, -0.28387734915924823);
        const double s177 = 1.0/x[2];
        const double s178 = pow(s177, 0.21471471471471473);
        const double s179 = 0.79174558284615615*s178 - 1;
        const double s180 = x[0]*s167 + 930.65395757192096;
        const double s181 = s176*s180;
        const double s182 = pow(s181, -1.0/2.0);
        const double s183 = pow(s177, 0.78528528528528529);
        const double s184 = pow(s177, 0.8926426426426427);
        const double s185 = sqrt(s183 - 0.88980086696190397*s184);
        const double s186 = s182*s185;
        const double s187 = s179*s186;
        const double s188 = 19879.505327143728*x[2];
        const double s189 = s187*s188;
        const double s190 = s176*s189;
        const double s191 = -0.27513219525120397*s166 - 1.3058145035967244*s170;
        const double s192 = s174*pow(s175, -1.2838773491592481)*s180;
        const double s193 = s189*s192;
        const double s194 = (1.0/2.0)*s176;
        const double s195 = (1.0/2.0)*s192;
        const double s196 = 5.0561342312332735e-5/x[3];
        const double s197 = pow(x[1], 2);
        const double s198 = pow(x[3], 2);
        const double s199 = 0.38110788300674586*s198 + 1;
        const double s200 = pow(s199, -7.0452961672473862);
        const double s201 = sqrt(-s197*s200 + 0.25558758431328077);
        const double s202 = 1.0/s201;
        const double s203 = 1.2136487758138137*pow(x[1], 0.28387734915924828) - 1;
        const double s204 = 602154.98178743932*s203;
        const double s205 = u[0]*pow(x[0], 2)*s169;
        const double s206 = s176*s205;
        const double s207 = x[2]*s187;
        const double s208 = s181*s187;
        const double s209 = s192*s205;
        const double s210 = 1370.3736029396543*x[0]*s165 - 170.79813918820744*s205;
        const double s211 = -0.27513219525120397*x[0]*s165 + 0.034291427446403858*s205;
		return {std::vector<double>{-s196*(s171*s190 + s189*(-s171*s194 - s191*s195) + s191*s193), -s196*(207458.88863946567*pow(x[1], -0.71612265084075166)*s201 - x[1]*s200*s202*s204 - 5643.3412748667151*s173*s208 + s189*(0.14193867457962411*s172*s181 + 1492.7223176934203*s206 - 0.29969635089763269*s209) - 59349162.533065364*s206*s207 + 11915.63040839005*s207*s209), -s196*(-3379.504513613173*s178*s181*s186 + s179*s181*s182*s188*(-0.39264264264264265*s177*s183 + 0.39713709865529423*s177*s184)/s185 + 22701.175964577087*s208), -0.49446021993750205*x[3] - 81.74742197369325*s197*pow(s199, -8.0452961672473862)*s202*s203 + 5.0561342312332735e-5*(s181*s189 + s201*s204)/s198}, {-s196*(s189*(-s194*s210 - s195*s211) + s190*s210 + s193*s211)}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        const double s5 = 0.057210285249880768*u[0] + x[1];
        const double s6 = 1.0/s5;
        const double s7 = 0.026260565610162729*u[0] + 0.45901826036110277*x[1];
        const double s8 = pow(s5, -2);
        const double s9 = 1.0/s7;
        const double s10 = 0.96919390034483355*u[0]*s9 + 1;
        const double s11 = pow(s10, 0.28387734915924823);
        const double s12 = 1.0/s11;
        const double s13 = s12*s8;
        const double s14 = 19879.505327143728*s13;
        const double s15 = 1370.3736029396543*s9;
        const double s16 = u[0]*s15 + 930.65395757192096;
        const double s17 = pow(s10, -1.2838773491592481);
        const double s18 = s16*s17;
        const double s19 = 19879.505327143728*s18;
        const double s20 = s19*s8;
        const double s21 = -s15;
        const double s22 = 170.79813918820744*u[0];
        const double s23 = s21 + s8*(2985.4446353868407*x[1] + s22);
        const double s24 = -0.27513219525120397*s9;
        const double s25 = 0.034291427446403858*u[0];
        const double s26 = s24 + s8*(0.59939270179526538*x[1] + s25);
        const double s27 = -s26;
        const double s28 = u[0]*s27;
        const double s29 = s17*s28;
        const double s30 = -s23;
        const double s31 = 39759.010654287456*s17;
        const double s32 = -1.2443260955960374*s9;
        const double s33 = 0.15508805862519931*u[0];
        const double s34 = -s32 - s8*(2.7108422541123849*x[1] + s33);
        const double s35 = pow(s10, -2.2838773491592481)*s16;
        const double s36 = 19879.505327143728*s35;
        const double s37 = s18*s27;
        const double s38 = s12*s30 - s37;
        const double s39 = 1.0/s16;
        const double s40 = 9939.7526635718641*s11*s39;
        const double s41 = 19879.505327143728*s38;
        const double s42 = u[0]*s23;
        const double s43 = pow(s10, -0.99999999999999989);
        const double s44 = 1.0/x[3];
        const double s45 = 1.0/x[2];
        const double s46 = pow(s45, 0.78528528528528529);
        const double s47 = pow(s45, 0.8926426426426427);
        const double s48 = s46 - 0.88980086696190397*s47;
        const double s49 = sqrt(s48);
        const double s50 = pow(s45, 0.21471471471471473);
        const double s51 = 0.79174558284615615*s50 - 1;
        const double s52 = s49*s51;
        const double s53 = 1.0/x[1];
        const double s54 = x[2]*s53;
        const double s55 = pow(s54, 0.28387734915924823);
        const double s56 = s12*s16;
        const double s57 = s55*s56;
        const double s58 = pow(s57, -1.0/2.0);
        const double s59 = s55*s58;
        const double s60 = x[2]*s59;
        const double s61 = s52*s60;
        const double s62 = 2.5280671156166367e-5*s44*s61;
        const double s63 = s24 + s25*s8;
        const double s64 = s17*s63;
        const double s65 = u[0]*s30;
        const double s66 = s21 + s22*s8;
        const double s67 = 341.59627837641489*u[0];
        const double s68 = pow(s5, -3);
        const double s69 = s68*(0.11442057049976154*u[0] + 2*x[1]);
        const double s70 = 9939.7526635718641*s12;
        const double s71 = s35*s63;
        const double s72 = 0.068582854892807715*u[0];
        const double s73 = -s25*s69 + s26 + s72*s8;
        const double s74 = 9939.7526635718641*s18;
        const double s75 = s43*s63;
        const double s76 = s12*s66;
        const double s77 = s18*s63;
        const double s78 = -s76 + s77;
        const double s79 = 1.0/x[0];
        const double s80 = 5.0561342312332735e-5*s44;
        const double s81 = s59*s80;
        const double s82 = s79*s81;
        const double s83 = x[2]*s52;
        const double s84 = -s82*s83*(4969.8763317859321*u[0]*s11*s38*s39*s78 + 9939.7526635718641*u[0]*s17*s27*s66 - 9939.7526635718641*u[0]*s34*s71 + 9939.7526635718641*u[0]*s38*s39*s66 - 9939.7526635718641*u[0]*s38*s75 + s19*s73 - 19879.505327143728*s42*s64 - 9939.7526635718641*s64*s65 - s70*(-s22*s69 + s23 + s67*s8) - s73*s74);
        const double s85 = -s6*(5970.8892707736813*x[1] + s67) + 5970.8892707736813;
        const double s86 = s70*s8;
        const double s87 = s12*s30;
        const double s88 = 2821.6706374333576*s53;
        const double s89 = u[0]*s8;
        const double s90 = s17*s89;
        const double s91 = 5957.8152041950252*s17;
        const double s92 = s74*s8;
        const double s93 = s18*s26;
        const double s94 = 5643.3412748667151*s53;
        const double s95 = 26945.101515836737*s89;
        const double s96 = u[0]*s13;
        const double s97 = s53*s56;
        const double s98 = s18*s89;
        const double s99 = 1492.7223176934203*s96 + 0.14193867457962411*s97 - 0.29969635089763269*s98;
        const double s100 = 19879.505327143728*s99;
        const double s101 = s100*s43;
        const double s102 = s100*s39;
        const double s103 = s40*s99;
        const double s104 = u[0]*s82;
        const double s105 = pow(x[1], 2);
        const double s106 = pow(x[3], 2);
        const double s107 = 0.38110788300674586*s106 + 1;
        const double s108 = pow(s107, -7.0452961672473862);
        const double s109 = -s105*s108 + 0.25558758431328077;
        const double s110 = sqrt(s109);
        const double s111 = 1.0/s110;
        const double s112 = s108*s111;
        const double s113 = 1.2136487758138137*pow(x[1], 0.28387734915924828) - 1;
        const double s114 = 602154.98178743932*s113;
        const double s115 = s112*s114;
        const double s116 = pow(s109, -3.0/2.0);
        const double s117 = u[0]*s68;
        const double s118 = s117*s12;
        const double s119 = pow(u[0], 2)/pow(s5, 4);
        const double s120 = s119*s17;
        const double s121 = s52*s96;
        const double s122 = s54*s59;
        const double s123 = 1.0/s105;
        const double s124 = s57*s58;
        const double s125 = s124*s83;
        const double s126 = s117*s18;
        const double s127 = s119*s35;
        const double s128 = s52*s98;
        const double s129 = s52*s99;
        const double s130 = 5643.3412748667151*s54;
        const double s131 = s129*s89;
        const double s132 = 11915.63040839005*s60;
        const double s133 = 59349162.533065364*s60;
        const double s134 = 19879.505327143728*s52;
        const double s135 = -s66;
        const double s136 = s12*s135;
        const double s137 = s6*s72 - 0.59939270179526549;
        const double s138 = s136 + s77;
        const double s139 = s81*s83*(-s101*s63 + s102*s66 + s103*s138 + s135*s89*s91 - s136*s88 + s137*s20 - s137*s92 + 29674581.266532682*s64*s89 + 11915.63040839005*s66*s90 - s71*s95 - s76*s94 + s77*s88 - s86*(s6*s67 - 2985.4446353868407));
        const double s140 = s49*s50;
        const double s141 = 3379.504513613173*s140;
        const double s142 = s12*s23;
        const double s143 = 22701.175964577087*s52;
        const double s144 = 0.39264264264264265*s46 - 0.39713709865529423*s47;
        const double s145 = 1.0/s49;
        const double s146 = s145*s51;
        const double s147 = 19879.505327143728*s144*s146;
        const double s148 = 1689.7522568065865*s140;
        const double s149 = s144*s146;
        const double s150 = 9939.7526635718641*s149;
        const double s151 = 19879.505327143728*s76;
        const double s152 = s81*(22701.175964577087*s12*s49*s51*s66 - s138*s148 - s138*s150 + 11350.587982288542*s138*s49*s51 - s141*s76 - s143*s77 + 19879.505327143728*s144*s145*s16*s17*s51*s63 - s149*s151 + 3379.504513613173*s16*s17*s49*s50*s63);
        const double s153 = 1.0/s106;
        const double s154 = 5.0561342312332735e-5*s153;
        const double s155 = s154*s61;
        const double s156 = pow(s107, -8.0452961672473862);
        const double s157 = s113*s156;
        const double s158 = x[3]*s113;
        const double s159 = s105*s111;
        const double s160 = -s155*(9939.7526635718641*s136 + s151 - 9939.7526635718641*s77);
        const double s161 = u[0]*s6;
        const double s162 = 19879.505327143728*s78;
		return {std::vector<double>{-u[0]*s62*(u[0]*pow(s38, 2)*s40 + s14*s7*(s6*(744.18886540088022*u[0] + 13007.955862314655*x[1]) - 13007.955862314655) - s20*s7*(s6*(0.14941204046840012*u[0] + 2.6116290071934487*x[1]) - 2.6116290071934487) + 79518.021308574913*s23*s29 + s28*s30*s31 + s28*s34*s36 + s28*s41*s43 + s39*s41*s42)/pow(x[0], 2), s104*s83*(s101*s27 - s102*s30 + s103*(s87 + s93) + s14*s85 + s26*s35*s95 - 29674581.266532682*s26*s90 + 53890.203031673474*s28*s35*s8 - 59349162.533065364*s29*s8 - s37*s94 - s65*s8*s91 - s85*s86 + s87*s88 - s88*s93 - s92*(-s6*(1.1987854035905308*x[1] + s72) + 1.1987854035905308)), s80*(148566.00927297046*pow(x[1], -1.7161226508407517)*s110 + 414917.77727893135*pow(x[1], 0.28387734915924834)*s112 + s105*pow(s107, -14.090592334494772)*s114*s116 - 44295909724.760422*s11*s39*s61*pow(s96 + 9.5087125647689213e-5*s97 - 0.00020077166887993512*s98, 2) + s115 - 118698325.06613073*s118*s61 + 71146909.759960771*s120*s61 - 33695765.869415939*s121*s122 + 6765.1551477901958*s122*s128 - 7245.3580363768506*s123*s125 + 23831.260816780101*s126*s61 - 32301.394395450166*s127*s61 + s129*s130*s59 - s131*s132*s43 + s131*s133*s39 + s134*s60*(2985.4446353868407*s118 - 1789.4537260646994*s120 + 0.182231849262465*s123*s56 - 0.59939270179526538*s126 + 0.81242953141659502*s127 + 847.5001091553147*s53*s96 - 0.17015401129103971*s53*s98)), -s104*(s141*s142 + s141*s37 - s142*s143 + s142*s147 - s143*s37 + s147*s37 + s148*s38 + s150*s38 - 11350.587982288542*s38*s52), -s81*(-s100*s149 - 67773104.000419348*s121 + 13606.919195337599*s128 + 17057.834689710369*s129 + s134*(423.75005457765735*s96 + 0.040293174682840868*s97 - 0.085077005645519857*s98) + 10089323.62043206*s140*s96 + 959.36478279622202*s140*s97 - 2025.6503411438939*s140*s98 - s141*s99 + 59349162.533065364*s149*s96 + 5643.3412748667151*s149*s97 - 11915.63040839005*s149*s98 - 6444.3496556217833*s52*s97), -s124*s45*s80*(-3613.2399488918518*s140 + 6759.009027226346*s144*s145*s50 + 19879.505327143728*s145*s51*(0.70097913228543862*s46 - 0.75163860789038806*s47) - 45402.351929154174*s149 + 3222.1748278108917*s49*s51 - 3135.3533387997218*s51*pow(0.98868286033243979*s46 - s47, 2)/pow(s48, 3.0/2.0)), -u[0]*s155*s79*(19879.505327143728*s37 - 9939.7526635718641*s87 + 9939.7526635718641*s93), -5.0561342312332735e-5*s111*(1616796.9091626296*pow(x[1], 3)*pow(s107, -15.090592334494772)*s113/s109 + 3233593.8183252593*x[1]*s157 + 557030.83105771779*pow(x[1], 1.2838773491592483)*s156) + 5.0561342312332735e-5*s153*(207458.88863946567*pow(x[1], -0.71612265084075166)*s110 - x[1]*s115 + s100*s61 - s121*s133 - s124*s130*s52 + s128*s132), -s124*s154*(s141 - s143 + s147), 219.49329188766245*pow(x[1], 4)*pow(s107, -16.090592334494772)*s116*s158 + 501.29575763463157*pow(s107, -9.0452961672473862)*s158*s159 + 81.74742197369325*s157*s159*s44 - 0.49446021993750205 - (60.891528318462406*s110*s113 + 2.0102689476911122*s125)/pow(x[3], 3)}, {s84, s139, s152, s160}, {-s62*(s14*(19.542820526212374*s161 - 341.59627837641489) + s162*s39*s66 - s162*s75 - s20*(0.0039236446916687106*s161 - 0.068582854892807715) - s31*s63*s66 + s36*s63*(s32 + s33*s8) + s40*pow(s78, 2))}, {}, {}, {}};
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
            {0.2498994156264608, 0.50615, 0.33926666666666666, 0.0681},  // x0
            {0.01818181818181818, 0.4044453658536586, 0.3370378048780488, 0.03},  // lb x
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
    
    // optimize
    int status = solver.solve();
    return status;
}        
        
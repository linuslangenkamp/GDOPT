
// CODEGEN FOR MODEL "dieselMotor"

#define _USE_MATH_DEFINES
#include <cmath>
#include <string>
#include "dieselMotor.h"
#include "constants.h"




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
        return sqrt(pow(x[2]/x[1], 0.28387734915924823)*pow((0.99545446595499631*u[0] + 0.45901826036110277*x[1])/(0.026260565610162729*u[0] + 0.45901826036110277*x[1]), -0.28387734915924823)*(1394.8131022528294*u[0] + 427.18716060283873*x[1])/(0.026260565610162729*u[0] + 0.45901826036110277*x[1]))*(-0.59850022616889709*u[1]*x[2]*sqrt(pow(1.0/x[2], 1.5705705705705706) - 0.79174558284615604*pow(1.0/x[2], 1.7852852852852854)) + 0.018144065864949988*x[0]*sqrt(pow(x[2]/x[1], 0.28387734915924823)*pow((0.99545446595499631*u[0] + 0.45901826036110277*x[1])/(0.026260565610162729*u[0] + 0.45901826036110277*x[1]), -0.28387734915924823)*(1394.8131022528294*u[0] + 427.18716060283873*x[1])/(0.026260565610162729*u[0] + 0.45901826036110277*x[1]))*(0.026260565610162729*u[0] + 0.45901826036110277*x[1]) - 1.0369170922727629*x[2]*sqrt(pow(1.0/x[2], 0.78528528528528529) - 0.88980086696190397*pow(1.0/x[2], 0.8926426426426427)));
    }

    std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        const double s112 = 0.45901826036110277*x[1];
        const double s113 = 0.026260565610162729*u[0] + s112;
        const double s114 = 1.0/s113;
        const double s115 = 0.99545446595499631*u[0] + s112;
        const double s116 = pow(s114*s115, 0.28387734915924823);
        const double s117 = 1.0/s116;
        const double s118 = 1.0/x[1];
        const double s119 = pow(x[2]*s118, 0.28387734915924823);
        const double s120 = 1394.8131022528294*u[0] + 427.18716060283873*x[1];
        const double s121 = s117*s119*s120;
        const double s122 = s114*s121;
        const double s123 = s114*s117*s119;
        const double s124 = pow(0.057210285249880768*u[0] + x[1], 2);
        const double s125 = 1.0/s124;
        const double s126 = s121*s125;
        const double s127 = s115*s125;
        const double s128 = (1.0/2.0)*s121/s115;
        const double s129 = -0.14193867457962411*s118*s122 + 213.59358030141937*s123 - 1.0892812839442543*s126 + s128*(-0.13030488696699949*s114 + 0.61844456674975457*s127);
        const double s130 = 1.0/s119;
        const double s131 = 1.0/s120;
        const double s132 = s116*s130*s131;
        const double s133 = sqrt(s122);
        const double s134 = 1.0/x[2];
        const double s135 = pow(s134, 0.78528528528528529);
        const double s136 = pow(s134, 0.8926426426426427);
        const double s137 = sqrt(s135 - 0.88980086696190397*s136);
        const double s138 = 1.0369170922727629*s137;
        const double s139 = pow(s134, 1.5705705705705706);
        const double s140 = pow(s134, 1.7852852852852854);
        const double s141 = sqrt(s139 - 0.79174558284615604*s140);
        const double s142 = 0.59850022616889709*s141;
        const double s143 = u[1]*s142;
        const double s144 = s133*(0.018144065864949988*x[0]*s113*s133 - x[2]*s138 - x[2]*s143);
        const double s145 = s113*s144;
        const double s146 = x[0]*s133;
        const double s147 = 0.21069776334493312*s124;
        const double s148 = 697.40655112641468*s123 - 0.062318092971807154*s126 + s128*(-0.28258697500403945*s114 + 0.035381390074992386*s127);
        return {std::vector<double>{0.018144065864949988*s113*s122, s129*s132*s145 + s133*(0.018144065864949988*x[0]*s116*s129*s130*s131*s133*s147 + 0.0083284575492066099*s146), s133*(-0.59850022616889709*u[1]*x[2]*(-0.78528528528528529*s134*s139 + 0.70674586937243211*s134*s140)/s141 + 0.0025753446603564022*x[0]*s113*s133*s134 - 1.0369170922727629*x[2]*(-0.39264264264264265*s134*s135 + 0.39713709865529423*s134*s136)/s137 - s138 - s143) + 0.14193867457962411*s134*s144}, {s132*s145*s148 + s133*(0.018144065864949988*x[0]*s116*s130*s131*s133*s147*s148 + 0.0004764734320816331*s146), -x[2]*s133*s142}, {}};
    }

    std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        const double s0 = 0.45901826036110277*x[1];
        const double s1 = 0.026260565610162729*u[0] + s0;
        const double s2 = 1.0/s1;
        const double s3 = 0.057210285249880768*u[0] + x[1];
        const double s4 = pow(s3, 2);
        const double s5 = 1.0/s4;
        const double s6 = s5*(0.035220562761846951*u[0] + 0.016240704121380591*x[1]);
        const double s7 = -0.28258697500403945*s2 + s6;
        const double s8 = 1394.8131022528294*u[0] + 427.18716060283873*x[1];
        const double s9 = 0.99545446595499631*u[0] + s0;
        const double s10 = 1.0/s9;
        const double s11 = s10*s8;
        const double s12 = s11*s7;
        const double s13 = s12 + 1394.8131022528294*s2 - s5*(173.84418516897315*u[0] + 53.242978381620034*x[1]);
        const double s14 = 1.0/x[1];
        const double s15 = pow(x[2]*s14, 0.28387734915924823)*pow(s2*s9, -0.28387734915924823);
        const double s16 = 0.0019114570478640559*s4;
        const double s17 = 1.0/s8;
        const double s18 = s13*s17;
        const double s19 = s16*s18 + 0.0004764734320816331;
        const double s20 = s15*s8;
        const double s21 = s2*s20;
        const double s22 = s19*s21;
        const double s23 = s13*s15*s16*s2 + s22;
        const double s24 = s14*s2;
        const double s25 = 0.28387734915924823*x[1];
        const double s26 = s5*(0.61563340591664595*u[0] + s25);
        const double s27 = -0.13030488696699949*s2 + s26;
        const double s28 = s11*s27;
        const double s29 = 427.18716060283873*s2 - s24*(395.95584604012066*u[0] + 121.2687587467999*x[1]) + s28 - s5*(3038.6876137684608*u[0] + 930.65395757192107*x[1]);
        const double s30 = 0.018144065864949988*s1;
        const double s31 = s17*s29;
        const double s32 = s16*s31 + 0.0083284575492066099;
        const double s33 = x[0]*s15;
        const double s34 = s29*s33;
        const double s35 = s5*s8;
        const double s36 = pow(u[0] + 0.30626838815384533*x[1], -2);
        const double s37 = s1*s36;
        const double s38 = s36*s4;
        const double s39 = s38*pow(s2 - s24*(0.92689079297550736*u[0] + s25) + 0.0023408943250747944*s28 - s5*(7.1132465907456588*u[0] + 2.1785625678885085*x[1]), 2);
        const double s40 = s1*s17;
        const double s41 = s29*s40;
        const double s42 = 0.28387734915924823*s14;
        const double s43 = s31*s4;
        const double s44 = 0.21069776334493312*s10;
        const double s45 = pow(s3, -3);
        const double s46 = s14*s5;
        const double s47 = 1.0/s3;
        const double s48 = s10*s35;
        const double s49 = s8/pow(u[0] + 0.46111427097847246*x[1], 2);
        const double s50 = 0.46321985258873616*s49;
        const double s51 = s1*s49;
        const double s52 = s1*s27;
        const double s53 = s48*s52;
        const double s54 = s10*(-111.32914935221271*s2 + s5*(525.98137329157373*u[0] + 242.53751749359981*x[1])) - 0.56775469831849645*s14*s28 - 242.53751749359981*s24 - s27*s50 + s45*(6077.3752275369216*u[0] + 1861.3079151438421*x[1]) + s46*(1725.2291694392643*u[0] + 528.38315692016079*x[1]) + s48*(-s47*(1.2312668118332919*u[0] + 0.56775469831849645*x[1]) + 0.56775469831849645) - 1861.3079151438421*s5 + 0.38597462881247235*s51*pow(-0.21069776334493312*s2 + s5*s9, 2) - 2.1785625678885085*s53 + s2*(508.35874199809757*u[0] + 155.69421251567385*x[1])/pow(x[1], 2);
        const double s55 = sqrt(s21);
        const double s56 = 1.0/x[2];
        const double s57 = pow(s56, 0.78528528528528529);
        const double s58 = pow(s56, 0.8926426426426427);
        const double s59 = s57 - 0.88980086696190397*s58;
        const double s60 = sqrt(s59);
        const double s61 = 1.0369170922727629*s60;
        const double s62 = pow(s56, 1.5705705705705706);
        const double s63 = pow(s56, 1.7852852852852854);
        const double s64 = s62 - 0.79174558284615604*s63;
        const double s65 = sqrt(s64);
        const double s66 = u[1]*s65;
        const double s67 = 0.59850022616889709*s66;
        const double s68 = x[0]*s55;
        const double s69 = x[2]*s61 + x[2]*s67 - s30*s68;
        const double s70 = s55*s69;
        const double s71 = (1.0/2.0)*s70;
        const double s72 = 4.1971133489432667e-7*s38;
        const double s73 = 0.00021871000590239524*u[0] + 0.0038229140957281118*x[1];
        const double s74 = 0.00054261935977941041*s14;
        const double s75 = s10*s16;
        const double s76 = s27*s75;
        const double s77 = s16*s17;
        const double s78 = x[0]*s21;
        const double s79 = s13*s33;
        const double s80 = s18*s70;
        const double s81 = s13*s37;
        const double s82 = s13*s40;
        const double s83 = s13*s38;
        const double s84 = s29*s83;
        const double s85 = s49*s7;
        const double s86 = s10*(-181.75096362914482*s2 + 1394.8131022528294*s26) + s10*(-120.71752747532098*s2 + 427.18716060283873*s6) - s12*s42 - 395.95584604012066*s24 + s45*(347.68837033794631*u[0] + 106.48595676324007*x[1]) + s46*(49.350426452517596*u[0] + 15.114475564317457*x[1]) + s48*(-s47*(0.070441125523693901*u[0] + 0.032481408242761182*x[1]) + 0.63187411003802652) - 3091.9305921500804*s5 - s50*s7 + 1.0091534315526534*s52*s85 - 0.12463618594361431*s53;
        const double s87 = s18*s4;
        const double s88 = 0.10534888167246656*s10*s27*s4*s80 - 0.14193867457962411*s14*s70*s82 + (1.0/2.0)*s19*s34 + (1.0/2.0)*s32*s79 - s40*s71*s86 + 0.00010978832492295084*s70*s81 - 2.7074964601006316e-8*s70*s84 + s78*(4.9124994101184804e-10*s1*s84 - s13*s72 + s18*s73 + 0.00023823671604081655*s41 + s74*s87 - s76*s82 + s77*s86) - 0.22950913018055139*s80;
        const double s89 = s41*s55;
        const double s90 = 0.29925011308444854*x[2]*s65;
        const double s91 = -s89*s90;
        const double s92 = s56*s78;
        const double s93 = s55*s56;
        const double s94 = s69*s93;
        const double s95 = 0.070969337289812057*s94;
        const double s96 = 60.634379373399952*s2 - s24*(56.201447978988433*u[0] + 17.212726884436979*x[1]) + 0.14193867457962411*s28 - s5*(431.30729235981607*u[0] + 132.0957892300402*x[1]);
        const double s97 = s40*s94;
        const double s98 = 0.0002713096798897052*s4;
        const double s99 = 0.0038229140957281118*s17*s4;
        const double s100 = 1.0/s60;
        const double s101 = s100*(0.40713786731130253*s57 - 0.41179824557128902*s58);
        const double s102 = 1.0/s65;
        const double s103 = u[1]*s102*(0.78528528528528529*s62 - 0.70674586937243211*s63);
        const double s104 = x[0]*s1*s93;
        const double s105 = (1.0/2.0)*s101 + 0.29925011308444854*s103 + 0.0012876723301782011*s104 - 1.0/2.0*s61 - 1.0/2.0*s67;
        const double s106 = 0.14193867457962411*s12 + 197.97792302006033*s2 - s5*(24.675213226258798*u[0] + 7.5572377821587287*x[1]);
        const double s107 = s55*s82;
        const double s108 = 0.14193867457962411*x[0]*s22*s56 + s105*s107 - s106*s97 + s82*s95 + s92*(s106*s99 - s18*s98 + 6.7630007422071558e-5);
        const double s109 = s55*(s102*(0.46999342085035012*s62 - 0.42298756266333443*s63) - 0.68345055500691565*s65);
        const double s110 = s38*pow(0.00071694193177914095*s12 + s2 - s5*(0.12463618594361431*u[0] + 0.038172123774593711*x[1]), 2);
        const double s111 = 0.12463618594361431*s1*s48*s7 - s10*(-788.31203052325407*s2 + s5*(98.252204817884447*u[0] + 45.305493796626351*x[1])) - s45*(19.891350845100092*u[0] + 6.0920919615314348*x[1]) + s48*(s47*(0.0040299568845331849*u[0] + 0.0018582706308861959*x[1]) - 0.070441125523693901) + 347.68837033794631*s5 - 0.080586349365681723*s51*pow(-s2 + s5*(0.12463618594361432*u[0] + 0.057471524018927057*x[1]), 2) + 1.0045662902728985*s85;
        return {std::vector<double>{s15*(s2*(11.616641711189867*u[0] + 3.5578101326468485*x[1]) + s29*s30), s32*s34 - s71*(s17*(196.08670730851406*s2 - s24*(181.75096362914482*u[0] + 55.664574676106355*x[1]) + 0.45901826036110277*s28 - s35) - s27*s43*s44 - 0.00021957664984590167*s29*s37 + 0.0098817594006143008*s39 + s40*s54 + s41*s42) + s78*(8.9647646713167302e-5*s1*s39 - s29*s72 + s31*s73 - s41*s76 + 0.004164228774603305*s41 + s43*s74 + s54*s77), 0.0051506893207128044*s20*s56, s105*s89 + 0.14193867457962411*s32*s92 + s41*s95 + s92*(-s31*s98 + s96*s99 + 0.0011821302258270508) - s96*s97, s93*(-0.59850022616889709*u[1]*s102*(2.0186312438564693*s62 - 1.9684888703991916*s63) + 0.36907891757467437*u[1]*pow(s62 - 0.89998613575915831*s63, 2)/pow(s64, 3.0/2.0) - s100*(0.72685724359330139*s57 - 0.77938691973364849*s58) + s100*(0.81427573462260505*s57 - 0.82359649114257805*s58) + 0.28387734915924823*s101 + 1.3669011100138313*s103 - 0.0014787216373263021*s104 + s56*(0.12628829696086882*x[2]*s60 + 0.072892591757646938*x[2]*s66 - 0.0022098036526797022*s1*s68) - 0.2943572754523075*s60 - 0.16990065767603704*s66 + 0.1635403604775221*pow(0.98868286033243979*s57 - s58, 2)/pow(s59, 3.0/2.0))}, {s23, s88, s108, s91, s109}, {s19*s79 + s71*(-0.10534888167246656*s110 + s111*s40 - s17*(0.026260565610162729*s12 + 36.628580985625042*s2 - s5*(4.5652466305750981*u[0] + 1.3981907270710088*x[1])) + s44*s7*s87 + 0.00071694193177914106*s81) + s78*(0.00095572852393202817*s1*s110 - s111*s77 + 0.0090720329324749938*s18*(0.0013792346123313226*u[0] + 0.024108158284950991*x[1]) - s7*s75*s82 + 0.00023823671604081655*s82 - 1.3704037084085102e-6*s83), -s107*s90}, {}, {}, {}};
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
        return -0.24723010996875103*pow(x[3], 2) - 5.0561342312332735e-5*(19879.505327143728*x[2]*pow(x[2]/x[1], 0.28387734915924823)*pow(0.96919390034483355*u[0]*x[0]/(0.026260565610162729*u[0]*x[0] + 0.45901826036110277*x[0]*x[1]) + 1, -0.28387734915924823)*(1370.3736029396543*u[0]*x[0]/(0.026260565610162729*u[0]*x[0] + 0.45901826036110277*x[0]*x[1]) + 930.65395757192096)*(0.79174558284615615*pow(1.0/x[2], 0.21471471471471473) - 1)*sqrt(pow(1.0/x[2], 0.78528528528528529) - 0.88980086696190397*pow(1.0/x[2], 0.8926426426426427))/sqrt(pow(x[2]/x[1], 0.28387734915924823)*pow(0.96919390034483355*u[0]*x[0]/(0.026260565610162729*u[0]*x[0] + 0.45901826036110277*x[0]*x[1]) + 1, -0.28387734915924823)*(1370.3736029396543*u[0]*x[0]/(0.026260565610162729*u[0]*x[0] + 0.45901826036110277*x[0]*x[1]) + 930.65395757192096)) + 602154.98178743932*(1.2136487758138137*pow(x[1], 0.28387734915924828) - 1)*sqrt(-pow(x[1], 2)*pow(0.38110788300674586*pow(x[3], 2) + 1, -7.0452961672473862) + 0.25558758431328077))/x[3];
    }

    std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        const double s158 = 0.026260565610162729*u[0];
        const double s159 = x[0]*x[1];
        const double s160 = 1.0/(x[0]*s158 + 0.45901826036110277*s159);
        const double s161 = u[0]*s160;
        const double s162 = 1370.3736029396543*s161;
        const double s163 = u[0]*x[0];
        const double s164 = pow(s159 + 0.057210285249880768*s163, -2);
        const double s165 = s163*s164*(-0.45901826036110277*x[1] - s158);
        const double s166 = s162 + 6503.9779311573275*s165;
        const double s167 = 1.0/x[1];
        const double s168 = x[2]*s167;
        const double s169 = pow(s168, 0.28387734915924823);
        const double s170 = 0.96919390034483355*x[0]*s161 + 1;
        const double s171 = s169*pow(s170, -0.28387734915924823);
        const double s172 = 1.0/x[2];
        const double s173 = pow(s172, 0.21471471471471473);
        const double s174 = 0.79174558284615615*s173 - 1;
        const double s175 = x[0]*s162 + 930.65395757192096;
        const double s176 = s171*s175;
        const double s177 = pow(s176, -1.0/2.0);
        const double s178 = pow(s172, 0.78528528528528529);
        const double s179 = pow(s172, 0.8926426426426427);
        const double s180 = sqrt(s178 - 0.88980086696190397*s179);
        const double s181 = s177*s180;
        const double s182 = s174*s181;
        const double s183 = 19879.505327143728*x[2];
        const double s184 = s182*s183;
        const double s185 = s171*s184;
        const double s186 = -0.27513219525120397*s161 - 1.3058145035967244*s165;
        const double s187 = s169*pow(s170, -1.2838773491592481)*s175;
        const double s188 = s184*s187;
        const double s189 = (1.0/2.0)*s171;
        const double s190 = (1.0/2.0)*s187;
        const double s191 = 5.0561342312332735e-5/x[3];
        const double s192 = pow(x[1], 2);
        const double s193 = pow(x[3], 2);
        const double s194 = 0.38110788300674586*s193 + 1;
        const double s195 = pow(s194, -7.0452961672473862);
        const double s196 = sqrt(-s192*s195 + 0.25558758431328077);
        const double s197 = 1.0/s196;
        const double s198 = 1.2136487758138137*pow(x[1], 0.28387734915924828) - 1;
        const double s199 = 602154.98178743932*s198;
        const double s200 = u[0]*pow(x[0], 2)*s164;
        const double s201 = s171*s200;
        const double s202 = x[2]*s182;
        const double s203 = s176*s182;
        const double s204 = s187*s200;
        const double s205 = 1370.3736029396543*x[0]*s160 - 170.79813918820744*s200;
        const double s206 = -0.27513219525120397*x[0]*s160 + 0.034291427446403858*s200;
        return {std::vector<double>{-s191*(s166*s185 + s184*(-s166*s189 - s186*s190) + s186*s188), -s191*(207458.88863946567*pow(x[1], -0.71612265084075166)*s196 - x[1]*s195*s197*s199 - 5643.3412748667151*s168*s203 + s184*(0.14193867457962411*s167*s176 + 1492.7223176934203*s201 - 0.29969635089763269*s204) - 59349162.533065364*s201*s202 + 11915.63040839005*s202*s204), -s191*(-3379.504513613173*s173*s176*s181 + s174*s176*s177*s183*(-0.39264264264264265*s172*s178 + 0.39713709865529423*s172*s179)/s180 + 22701.175964577087*s203), -0.49446021993750205*x[3] - 81.74742197369325*s192*pow(s194, -8.0452961672473862)*s197*s198 + 5.0561342312332735e-5*(s176*s184 + s196*s199)/s193}, {-s191*(s184*(-s189*s205 - s190*s206) + s185*s205 + s188*s206)}, {}};
    }

    std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        const double s0 = 0.057210285249880768*u[0] + x[1];
        const double s1 = 1.0/s0;
        const double s2 = 0.026260565610162729*u[0] + 0.45901826036110277*x[1];
        const double s3 = pow(s0, -2);
        const double s4 = 1.0/s2;
        const double s5 = 0.96919390034483355*u[0]*s4 + 1;
        const double s6 = pow(s5, 0.28387734915924823);
        const double s7 = 1.0/s6;
        const double s8 = s3*s7;
        const double s9 = 19879.505327143728*s8;
        const double s10 = 1370.3736029396543*s4;
        const double s11 = u[0]*s10 + 930.65395757192096;
        const double s12 = pow(s5, -1.2838773491592481);
        const double s13 = s11*s12;
        const double s14 = 19879.505327143728*s13;
        const double s15 = s14*s3;
        const double s16 = -s10;
        const double s17 = 170.79813918820744*u[0];
        const double s18 = s16 + s3*(2985.4446353868407*x[1] + s17);
        const double s19 = -0.27513219525120397*s4;
        const double s20 = 0.034291427446403858*u[0];
        const double s21 = s19 + s3*(0.59939270179526538*x[1] + s20);
        const double s22 = -s21;
        const double s23 = u[0]*s22;
        const double s24 = s12*s23;
        const double s25 = -s18;
        const double s26 = 39759.010654287456*s12;
        const double s27 = -1.2443260955960374*s4;
        const double s28 = 0.15508805862519931*u[0];
        const double s29 = -s27 - s3*(2.7108422541123849*x[1] + s28);
        const double s30 = s11*pow(s5, -2.2838773491592481);
        const double s31 = 19879.505327143728*s30;
        const double s32 = s13*s22;
        const double s33 = s25*s7 - s32;
        const double s34 = 1.0/s11;
        const double s35 = 9939.7526635718641*s34*s6;
        const double s36 = 19879.505327143728*s33;
        const double s37 = u[0]*s18;
        const double s38 = pow(s5, -0.99999999999999989);
        const double s39 = 1.0/x[3];
        const double s40 = 1.0/x[2];
        const double s41 = pow(s40, 0.78528528528528529);
        const double s42 = pow(s40, 0.8926426426426427);
        const double s43 = s41 - 0.88980086696190397*s42;
        const double s44 = sqrt(s43);
        const double s45 = pow(s40, 0.21471471471471473);
        const double s46 = 0.79174558284615615*s45 - 1;
        const double s47 = s44*s46;
        const double s48 = 1.0/x[1];
        const double s49 = x[2]*s48;
        const double s50 = pow(s49, 0.28387734915924823);
        const double s51 = s11*s7;
        const double s52 = s50*s51;
        const double s53 = pow(s52, -1.0/2.0);
        const double s54 = s50*s53;
        const double s55 = x[2]*s54;
        const double s56 = s47*s55;
        const double s57 = 2.5280671156166367e-5*s39*s56;
        const double s58 = s19 + s20*s3;
        const double s59 = s12*s58;
        const double s60 = u[0]*s25;
        const double s61 = s16 + s17*s3;
        const double s62 = 341.59627837641489*u[0];
        const double s63 = pow(s0, -3);
        const double s64 = s63*(0.11442057049976154*u[0] + 2*x[1]);
        const double s65 = 9939.7526635718641*s7;
        const double s66 = s30*s58;
        const double s67 = 0.068582854892807715*u[0];
        const double s68 = -s20*s64 + s21 + s3*s67;
        const double s69 = 9939.7526635718641*s13;
        const double s70 = s38*s58;
        const double s71 = s61*s7;
        const double s72 = s13*s58;
        const double s73 = -s71 + s72;
        const double s74 = 1.0/x[0];
        const double s75 = 5.0561342312332735e-5*s39;
        const double s76 = s54*s75;
        const double s77 = s74*s76;
        const double s78 = x[2]*s47;
        const double s79 = -s77*s78*(9939.7526635718641*u[0]*s12*s22*s61 - 9939.7526635718641*u[0]*s29*s66 + 4969.8763317859321*u[0]*s33*s34*s6*s73 + 9939.7526635718641*u[0]*s33*s34*s61 - 9939.7526635718641*u[0]*s33*s70 + s14*s68 - 19879.505327143728*s37*s59 - 9939.7526635718641*s59*s60 - s65*(-s17*s64 + s18 + s3*s62) - s68*s69);
        const double s80 = -s1*(5970.8892707736813*x[1] + s62) + 5970.8892707736813;
        const double s81 = s3*s65;
        const double s82 = s25*s7;
        const double s83 = 2821.6706374333576*s48;
        const double s84 = u[0]*s3;
        const double s85 = s12*s84;
        const double s86 = 5957.8152041950252*s12;
        const double s87 = s3*s69;
        const double s88 = s13*s21;
        const double s89 = 5643.3412748667151*s48;
        const double s90 = 26945.101515836737*s84;
        const double s91 = u[0]*s8;
        const double s92 = s48*s51;
        const double s93 = s13*s84;
        const double s94 = 1492.7223176934203*s91 + 0.14193867457962411*s92 - 0.29969635089763269*s93;
        const double s95 = 19879.505327143728*s94;
        const double s96 = s38*s95;
        const double s97 = s34*s95;
        const double s98 = s35*s94;
        const double s99 = u[0]*s77;
        const double s100 = pow(x[1], 2);
        const double s101 = pow(x[3], 2);
        const double s102 = 0.38110788300674586*s101 + 1;
        const double s103 = pow(s102, -7.0452961672473862);
        const double s104 = -s100*s103 + 0.25558758431328077;
        const double s105 = sqrt(s104);
        const double s106 = 1.0/s105;
        const double s107 = s103*s106;
        const double s108 = 1.2136487758138137*pow(x[1], 0.28387734915924828) - 1;
        const double s109 = 602154.98178743932*s108;
        const double s110 = s107*s109;
        const double s111 = pow(s104, -3.0/2.0);
        const double s112 = u[0]*s63;
        const double s113 = s112*s7;
        const double s114 = pow(u[0], 2)/pow(s0, 4);
        const double s115 = s114*s12;
        const double s116 = s47*s91;
        const double s117 = s49*s54;
        const double s118 = 1.0/s100;
        const double s119 = s52*s53;
        const double s120 = s119*s78;
        const double s121 = s112*s13;
        const double s122 = s114*s30;
        const double s123 = s47*s93;
        const double s124 = s47*s94;
        const double s125 = 5643.3412748667151*s49;
        const double s126 = s124*s84;
        const double s127 = 11915.63040839005*s55;
        const double s128 = 59349162.533065364*s55;
        const double s129 = 19879.505327143728*s47;
        const double s130 = -s61;
        const double s131 = s130*s7;
        const double s132 = s1*s67 - 0.59939270179526549;
        const double s133 = s131 + s72;
        const double s134 = s76*s78*(s130*s84*s86 - s131*s83 + s132*s15 - s132*s87 + s133*s98 - s58*s96 + 29674581.266532682*s59*s84 + 11915.63040839005*s61*s85 + s61*s97 - s66*s90 - s71*s89 + s72*s83 - s81*(s1*s62 - 2985.4446353868407));
        const double s135 = s44*s45;
        const double s136 = 3379.504513613173*s135;
        const double s137 = s18*s7;
        const double s138 = 22701.175964577087*s47;
        const double s139 = 0.39264264264264265*s41 - 0.39713709865529423*s42;
        const double s140 = 1.0/s44;
        const double s141 = s140*s46;
        const double s142 = 19879.505327143728*s139*s141;
        const double s143 = 1689.7522568065865*s135;
        const double s144 = s139*s141;
        const double s145 = 9939.7526635718641*s144;
        const double s146 = 19879.505327143728*s71;
        const double s147 = s76*(19879.505327143728*s11*s12*s139*s140*s46*s58 + 3379.504513613173*s11*s12*s44*s45*s58 - s133*s143 - s133*s145 + 11350.587982288542*s133*s44*s46 - s136*s71 - s138*s72 - s144*s146 + 22701.175964577087*s44*s46*s61*s7);
        const double s148 = 1.0/s101;
        const double s149 = 5.0561342312332735e-5*s148;
        const double s150 = s149*s56;
        const double s151 = pow(s102, -8.0452961672473862);
        const double s152 = s108*s151;
        const double s153 = x[3]*s108;
        const double s154 = s100*s106;
        const double s155 = -s150*(9939.7526635718641*s131 + s146 - 9939.7526635718641*s72);
        const double s156 = u[0]*s1;
        const double s157 = 19879.505327143728*s73;
        return {std::vector<double>{-u[0]*s57*(u[0]*pow(s33, 2)*s35 - s15*s2*(s1*(0.14941204046840012*u[0] + 2.6116290071934487*x[1]) - 2.6116290071934487) + 79518.021308574913*s18*s24 + s2*s9*(s1*(744.18886540088022*u[0] + 13007.955862314655*x[1]) - 13007.955862314655) + s23*s25*s26 + s23*s29*s31 + s23*s36*s38 + s34*s36*s37)/pow(x[0], 2), s78*s99*(s21*s30*s90 - 29674581.266532682*s21*s85 + s22*s96 + 53890.203031673474*s23*s3*s30 - 59349162.533065364*s24*s3 - s25*s97 - s3*s60*s86 - s32*s89 - s80*s81 + s80*s9 + s82*s83 - s83*s88 - s87*(-s1*(1.1987854035905308*x[1] + s67) + 1.1987854035905308) + s98*(s82 + s88)), s75*(148566.00927297046*pow(x[1], -1.7161226508407517)*s105 + 414917.77727893135*pow(x[1], 0.28387734915924834)*s107 + s100*pow(s102, -14.090592334494772)*s109*s111 + s110 - 118698325.06613073*s113*s56 + 71146909.759960771*s115*s56 - 33695765.869415939*s116*s117 + 6765.1551477901958*s117*s123 - 7245.3580363768506*s118*s120 + 23831.260816780101*s121*s56 - 32301.394395450166*s122*s56 + s124*s125*s54 - s126*s127*s38 + s126*s128*s34 + s129*s55*(2985.4446353868407*s113 - 1789.4537260646994*s115 + 0.182231849262465*s118*s51 - 0.59939270179526538*s121 + 0.81242953141659502*s122 + 847.5001091553147*s48*s91 - 0.17015401129103971*s48*s93) - 44295909724.760422*s34*s56*s6*pow(s91 + 9.5087125647689213e-5*s92 - 0.00020077166887993512*s93, 2)), -s99*(s136*s137 + s136*s32 - s137*s138 + s137*s142 - s138*s32 + s142*s32 + s143*s33 + s145*s33 - 11350.587982288542*s33*s47), -s76*(-67773104.000419348*s116 + 13606.919195337599*s123 + 17057.834689710369*s124 + s129*(423.75005457765735*s91 + 0.040293174682840868*s92 - 0.085077005645519857*s93) + 10089323.62043206*s135*s91 + 959.36478279622202*s135*s92 - 2025.6503411438939*s135*s93 - s136*s94 + 59349162.533065364*s144*s91 + 5643.3412748667151*s144*s92 - 11915.63040839005*s144*s93 - s144*s95 - 6444.3496556217833*s47*s92), -s119*s40*s75*(-3613.2399488918518*s135 + 6759.009027226346*s139*s140*s45 + 19879.505327143728*s140*s46*(0.70097913228543862*s41 - 0.75163860789038806*s42) - 45402.351929154174*s144 + 3222.1748278108917*s44*s46 - 3135.3533387997218*s46*pow(0.98868286033243979*s41 - s42, 2)/pow(s43, 3.0/2.0)), -u[0]*s150*s74*(19879.505327143728*s32 - 9939.7526635718641*s82 + 9939.7526635718641*s88), -5.0561342312332735e-5*s106*(1616796.9091626296*pow(x[1], 3)*pow(s102, -15.090592334494772)*s108/s104 + 3233593.8183252593*x[1]*s152 + 557030.83105771779*pow(x[1], 1.2838773491592483)*s151) + 5.0561342312332735e-5*s148*(207458.88863946567*pow(x[1], -0.71612265084075166)*s105 - x[1]*s110 - s116*s128 - s119*s125*s47 + s123*s127 + s56*s95), -s119*s149*(s136 - s138 + s142), 219.49329188766245*pow(x[1], 4)*pow(s102, -16.090592334494772)*s111*s153 + 501.29575763463157*pow(s102, -9.0452961672473862)*s153*s154 + 81.74742197369325*s152*s154*s39 - 0.49446021993750205 - (60.891528318462406*s105*s108 + 2.0102689476911122*s120)/pow(x[3], 3)}, {s79, s134, s147, s155}, {-s57*(-s15*(0.0039236446916687106*s156 - 0.068582854892807715) + s157*s34*s61 - s157*s70 - s26*s58*s61 + s31*s58*(s27 + s28*s3) + s35*pow(s73, 2) + s9*(19.542820526212374*s156 - 341.59627837641489))}, {}, {}, {}};
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

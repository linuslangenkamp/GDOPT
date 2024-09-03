
// CODEGEN FOR MODEL "generalizedBatchReactor"

// includes
#define _USE_MATH_DEFINES
#include "generalizedBatchReactorGeneratedParams.h"
#include <cmath>
#include <string>
#include "constants.h"
#include <problem.h>
#include "integrator.h"
#include "mesh.h"
#include "gdop.h"
#include "solver.h"


// runtime parameters and global constants
const double EXPONENT_ENERGY = EXPONENT_ENERGY_VALUE;
const double DEPLETION_COEFF = DEPLETION_COEFF_VALUE;


// mayer term
class MayergeneralizedBatchReactor : public Expression {
public:
	static std::unique_ptr<MayergeneralizedBatchReactor> create() {
		Adjacency adj{{0}, {}, {}};
		AdjacencyDiff adjDiff{{}, {}, {}, {}, {}, {}};
		return std::unique_ptr<MayergeneralizedBatchReactor>(new MayergeneralizedBatchReactor(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return -x[0];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-1}, {}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {}, {}, {}, {}, {}};
	}
private:
	MayergeneralizedBatchReactor(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


// dynamic constraints
class F0generalizedBatchReactor : public Expression {
public:
	static std::unique_ptr<F0generalizedBatchReactor> create() {
		Adjacency adj{{0, 1, 2, 3, 4, 5}, {0}, {}};
		AdjacencyDiff adjDiff{{}, {{0, 0}, {0, 1}, {0, 2}, {0, 3}, {0, 4}, {0, 5}}, {{0, 0}}, {}, {}, {}};
		return std::unique_ptr<F0generalizedBatchReactor>(new F0generalizedBatchReactor(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return 0.8552241295636*u[0]*x[1] + 0.292864672862223*u[0]*x[2] + 0.342951321558795*u[0]*x[3] + 0.11005832421052*u[0]*x[4] + 0.311601109024375*u[0]*x[5] - pow(u[0], 2)*x[0]*DEPLETION_COEFF_VALUE;
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-pow(u[0], 2)*DEPLETION_COEFF_VALUE, 0.8552241295636*u[0], 0.292864672862223*u[0], 0.342951321558795*u[0], 0.11005832421052*u[0], 0.311601109024375*u[0]}, {0.8552241295636*x[1] + 0.292864672862223*x[2] + 0.342951321558795*x[3] + 0.11005832421052*x[4] + 0.311601109024375*x[5] - 2*u[0]*x[0]*DEPLETION_COEFF_VALUE}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        const double x0 = 2*DEPLETION_COEFF_VALUE;
        const double x1 = -x0*u[0];
		return {std::vector<double>{}, {x1, 0.8552241295636, 0.292864672862223, 0.342951321558795, 0.11005832421052, 0.311601109024375}, {-x0*x[0]}, {}, {}, {}};
	}
private:
	F0generalizedBatchReactor(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F1generalizedBatchReactor : public Expression {
public:
	static std::unique_ptr<F1generalizedBatchReactor> create() {
		Adjacency adj{{1, 2, 3, 4, 5, 6}, {0}, {}};
		AdjacencyDiff adjDiff{{}, {{0, 1}, {0, 2}, {0, 3}, {0, 4}, {0, 5}, {0, 6}}, {}, {}, {}, {}};
		return std::unique_ptr<F1generalizedBatchReactor>(new F1generalizedBatchReactor(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return -0.8552241295636*u[0]*x[1] + 0.453291327685434*u[0]*x[2] + 0.962012695622935*u[0]*x[3] + 0.659136923379733*u[0]*x[4] + 0.616813662433078*u[0]*x[5] + 0.490811454960079*u[0]*x[6];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-0.8552241295636*u[0], 0.453291327685434*u[0], 0.962012695622935*u[0], 0.659136923379733*u[0], 0.616813662433078*u[0], 0.490811454960079*u[0]}, {-0.8552241295636*x[1] + 0.453291327685434*x[2] + 0.962012695622935*x[3] + 0.659136923379733*x[4] + 0.616813662433078*x[5] + 0.490811454960079*x[6]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-0.8552241295636, 0.453291327685434, 0.962012695622935, 0.659136923379733, 0.616813662433078, 0.490811454960079}, {}, {}, {}, {}};
	}
private:
	F1generalizedBatchReactor(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F2generalizedBatchReactor : public Expression {
public:
	static std::unique_ptr<F2generalizedBatchReactor> create() {
		Adjacency adj{{2, 3, 4, 5, 6, 7}, {0}, {}};
		AdjacencyDiff adjDiff{{}, {{0, 2}, {0, 3}, {0, 4}, {0, 5}, {0, 6}, {0, 7}}, {}, {}, {}, {}};
		return std::unique_ptr<F2generalizedBatchReactor>(new F2generalizedBatchReactor(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return -0.746156000547657*u[0]*x[2] + 0.456953645030287*u[0]*x[3] + 0.705827793744213*u[0]*x[4] + 0.372838217499846*u[0]*x[5] + 0.960078238377733*u[0]*x[6] + 0.938710138954149*u[0]*x[7];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-0.746156000547657*u[0], 0.456953645030287*u[0], 0.705827793744213*u[0], 0.372838217499846*u[0], 0.960078238377733*u[0], 0.938710138954149*u[0]}, {-0.746156000547657*x[2] + 0.456953645030287*x[3] + 0.705827793744213*x[4] + 0.372838217499846*x[5] + 0.960078238377733*x[6] + 0.938710138954149*x[7]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-0.746156000547657, 0.456953645030287, 0.705827793744213, 0.372838217499846, 0.960078238377733, 0.938710138954149}, {}, {}, {}, {}};
	}
private:
	F2generalizedBatchReactor(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F3generalizedBatchReactor : public Expression {
public:
	static std::unique_ptr<F3generalizedBatchReactor> create() {
		Adjacency adj{{3, 4, 5, 6, 7, 8}, {0}, {}};
		AdjacencyDiff adjDiff{{}, {{0, 3}, {0, 4}, {0, 5}, {0, 6}, {0, 7}, {0, 8}}, {}, {}, {}, {}};
		return std::unique_ptr<F3generalizedBatchReactor>(new F3generalizedBatchReactor(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return -1.76191766221202*u[0]*x[3] + 0.203505927974201*u[0]*x[4] + 0.671139779294799*u[0]*x[5] + 0.166092254776028*u[0]*x[6] + 0.653581553338871*u[0]*x[7] + 0.657998985282944*u[0]*x[8];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-1.76191766221202*u[0], 0.203505927974201*u[0], 0.671139779294799*u[0], 0.166092254776028*u[0], 0.653581553338871*u[0], 0.657998985282944*u[0]}, {-1.76191766221202*x[3] + 0.203505927974201*x[4] + 0.671139779294799*x[5] + 0.166092254776028*x[6] + 0.653581553338871*x[7] + 0.657998985282944*x[8]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-1.76191766221202, 0.203505927974201, 0.671139779294799, 0.166092254776028, 0.653581553338871, 0.657998985282944}, {}, {}, {}, {}};
	}
private:
	F3generalizedBatchReactor(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F4generalizedBatchReactor : public Expression {
public:
	static std::unique_ptr<F4generalizedBatchReactor> create() {
		Adjacency adj{{4, 5, 6, 7, 8, 9}, {0}, {}};
		AdjacencyDiff adjDiff{{}, {{0, 4}, {0, 5}, {0, 6}, {0, 7}, {0, 8}, {0, 9}}, {}, {}, {}, {}};
		return std::unique_ptr<F4generalizedBatchReactor>(new F4generalizedBatchReactor(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return -1.67852896930867*u[0]*x[4] + 0.926087218355965*u[0]*x[5] + 0.525113317807733*u[0]*x[6] + 0.0270696841378247*u[0]*x[7] + 0.693649834784741*u[0]*x[8] + 0.476861308321223*u[0]*x[9];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-1.67852896930867*u[0], 0.926087218355965*u[0], 0.525113317807733*u[0], 0.0270696841378247*u[0], 0.693649834784741*u[0], 0.476861308321223*u[0]}, {-1.67852896930867*x[4] + 0.926087218355965*x[5] + 0.525113317807733*x[6] + 0.0270696841378247*x[7] + 0.693649834784741*x[8] + 0.476861308321223*x[9]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-1.67852896930867, 0.926087218355965, 0.525113317807733, 0.0270696841378247, 0.693649834784741, 0.476861308321223}, {}, {}, {}, {}};
	}
private:
	F4generalizedBatchReactor(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F5generalizedBatchReactor : public Expression {
public:
	static std::unique_ptr<F5generalizedBatchReactor> create() {
		Adjacency adj{{5, 6, 7, 8, 9, 10}, {0}, {}};
		AdjacencyDiff adjDiff{{}, {{0, 5}, {0, 6}, {0, 7}, {0, 8}, {0, 9}, {0, 10}}, {}, {}, {}, {}};
		return std::unique_ptr<F5generalizedBatchReactor>(new F5generalizedBatchReactor(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return 0.886782915330311*u[0]*x[10] - 2.89847998660806*u[0]*x[5] + 0.662749072714428*u[0]*x[6] + 0.828640761528211*u[0]*x[7] + 0.758165838937764*u[0]*x[8] + 0.732629293026911*u[0]*x[9];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-2.89847998660806*u[0], 0.662749072714428*u[0], 0.828640761528211*u[0], 0.758165838937764*u[0], 0.732629293026911*u[0], 0.886782915330311*u[0]}, {0.886782915330311*x[10] - 2.89847998660806*x[5] + 0.662749072714428*x[6] + 0.828640761528211*x[7] + 0.758165838937764*x[8] + 0.732629293026911*x[9]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-2.89847998660806, 0.662749072714428, 0.828640761528211, 0.758165838937764, 0.732629293026911, 0.886782915330311}, {}, {}, {}, {}};
	}
private:
	F5generalizedBatchReactor(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F6generalizedBatchReactor : public Expression {
public:
	static std::unique_ptr<F6generalizedBatchReactor> create() {
		Adjacency adj{{6, 7, 8, 9, 10, 11}, {0}, {}};
		AdjacencyDiff adjDiff{{}, {{0, 6}, {0, 7}, {0, 8}, {0, 9}, {0, 10}, {0, 11}}, {}, {}, {}, {}};
		return std::unique_ptr<F6generalizedBatchReactor>(new F6generalizedBatchReactor(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return 0.595095939995323*u[0]*x[10] + 0.00452751637940219*u[0]*x[11] - 2.804844338636*u[0]*x[6] + 0.837096950885916*u[0]*x[7] + 0.658096821640587*u[0]*x[8] + 0.538460202958348*u[0]*x[9];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-2.804844338636*u[0], 0.837096950885916*u[0], 0.658096821640587*u[0], 0.538460202958348*u[0], 0.595095939995323*u[0], 0.00452751637940219*u[0]}, {0.595095939995323*x[10] + 0.00452751637940219*x[11] - 2.804844338636*x[6] + 0.837096950885916*x[7] + 0.658096821640587*x[8] + 0.538460202958348*x[9]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-2.804844338636, 0.837096950885916, 0.658096821640587, 0.538460202958348, 0.595095939995323, 0.00452751637940219}, {}, {}, {}, {}};
	}
private:
	F6generalizedBatchReactor(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F7generalizedBatchReactor : public Expression {
public:
	static std::unique_ptr<F7generalizedBatchReactor> create() {
		Adjacency adj{{7, 8, 9, 10, 11, 12}, {0}, {}};
		AdjacencyDiff adjDiff{{}, {{0, 7}, {0, 8}, {0, 9}, {0, 10}, {0, 11}, {0, 12}}, {}, {}, {}, {}};
		return std::unique_ptr<F7generalizedBatchReactor>(new F7generalizedBatchReactor(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return 0.892504982402579*u[0]*x[10] + 0.0898400313811519*u[0]*x[11] + 0.456327961909505*u[0]*x[12] - 3.28509908884497*u[0]*x[7] + 0.0323502921590433*u[0]*x[8] + 0.499867503130388*u[0]*x[9];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-3.28509908884497*u[0], 0.0323502921590433*u[0], 0.499867503130388*u[0], 0.892504982402579*u[0], 0.0898400313811519*u[0], 0.456327961909505*u[0]}, {0.892504982402579*x[10] + 0.0898400313811519*x[11] + 0.456327961909505*x[12] - 3.28509908884497*x[7] + 0.0323502921590433*x[8] + 0.499867503130388*x[9]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-3.28509908884497, 0.0323502921590433, 0.499867503130388, 0.892504982402579, 0.0898400313811519, 0.456327961909505}, {}, {}, {}, {}};
	}
private:
	F7generalizedBatchReactor(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F8generalizedBatchReactor : public Expression {
public:
	static std::unique_ptr<F8generalizedBatchReactor> create() {
		Adjacency adj{{8, 9, 10, 11, 12, 13}, {0}, {}};
		AdjacencyDiff adjDiff{{}, {{0, 8}, {0, 9}, {0, 10}, {0, 11}, {0, 12}, {0, 13}}, {}, {}, {}, {}};
		return std::unique_ptr<F8generalizedBatchReactor>(new F8generalizedBatchReactor(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return 0.0166139035604809*u[0]*x[10] + 0.352274002993723*u[0]*x[11] + 0.960845220231029*u[0]*x[12] + 0.742077528972963*u[0]*x[13] - 2.80026177280508*u[0]*x[8] + 0.889290646797178*u[0]*x[9];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-2.80026177280508*u[0], 0.889290646797178*u[0], 0.0166139035604809*u[0], 0.352274002993723*u[0], 0.960845220231029*u[0], 0.742077528972963*u[0]}, {0.0166139035604809*x[10] + 0.352274002993723*x[11] + 0.960845220231029*x[12] + 0.742077528972963*x[13] - 2.80026177280508*x[8] + 0.889290646797178*x[9]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-2.80026177280508, 0.889290646797178, 0.0166139035604809, 0.352274002993723, 0.960845220231029, 0.742077528972963}, {}, {}, {}, {}};
	}
private:
	F8generalizedBatchReactor(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F9generalizedBatchReactor : public Expression {
public:
	static std::unique_ptr<F9generalizedBatchReactor> create() {
		Adjacency adj{{9, 10, 11, 12, 13, 14}, {0}, {}};
		AdjacencyDiff adjDiff{{}, {{0, 9}, {0, 10}, {0, 11}, {0, 12}, {0, 13}, {0, 14}}, {}, {}, {}, {}};
		return std::unique_ptr<F9generalizedBatchReactor>(new F9generalizedBatchReactor(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return 0.302107579323058*u[0]*x[10] + 0.387416630575024*u[0]*x[11] + 0.655872767844423*u[0]*x[12] + 0.164857802428697*u[0]*x[13] + 0.82847780190645*u[0]*x[14] - 3.13710895423405*u[0]*x[9];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-3.13710895423405*u[0], 0.302107579323058*u[0], 0.387416630575024*u[0], 0.655872767844423*u[0], 0.164857802428697*u[0], 0.82847780190645*u[0]}, {0.302107579323058*x[10] + 0.387416630575024*x[11] + 0.655872767844423*x[12] + 0.164857802428697*x[13] + 0.82847780190645*x[14] - 3.13710895423405*x[9]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-3.13710895423405, 0.302107579323058, 0.387416630575024, 0.655872767844423, 0.164857802428697, 0.82847780190645}, {}, {}, {}, {}};
	}
private:
	F9generalizedBatchReactor(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F10generalizedBatchReactor : public Expression {
public:
	static std::unique_ptr<F10generalizedBatchReactor> create() {
		Adjacency adj{{10, 11, 12, 13, 14}, {0}, {}};
		AdjacencyDiff adjDiff{{}, {{0, 10}, {0, 11}, {0, 12}, {0, 13}, {0, 14}}, {}, {}, {}, {}};
		return std::unique_ptr<F10generalizedBatchReactor>(new F10generalizedBatchReactor(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return -2.69310532061175*u[0]*x[10] + 0.985372000161767*u[0]*x[11] + 0.766013677675046*u[0]*x[12] + 0.482682241153235*u[0]*x[13] + 0.643629829888353*u[0]*x[14];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-2.69310532061175*u[0], 0.985372000161767*u[0], 0.766013677675046*u[0], 0.482682241153235*u[0], 0.643629829888353*u[0]}, {-2.69310532061175*x[10] + 0.985372000161767*x[11] + 0.766013677675046*x[12] + 0.482682241153235*x[13] + 0.643629829888353*x[14]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-2.69310532061175, 0.985372000161767, 0.766013677675046, 0.482682241153235, 0.643629829888353}, {}, {}, {}, {}};
	}
private:
	F10generalizedBatchReactor(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F11generalizedBatchReactor : public Expression {
public:
	static std::unique_ptr<F11generalizedBatchReactor> create() {
		Adjacency adj{{11, 12, 13, 14}, {0}, {}};
		AdjacencyDiff adjDiff{{}, {{0, 11}, {0, 12}, {0, 13}, {0, 14}}, {}, {}, {}, {}};
		return std::unique_ptr<F11generalizedBatchReactor>(new F11generalizedBatchReactor(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return -1.81943018149107*u[0]*x[11] + 0.31177572374159*u[0]*x[12] + 0.917145607057463*u[0]*x[13] + 0.380363765970804*u[0]*x[14];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-1.81943018149107*u[0], 0.31177572374159*u[0], 0.917145607057463*u[0], 0.380363765970804*u[0]}, {-1.81943018149107*x[11] + 0.31177572374159*x[12] + 0.917145607057463*x[13] + 0.380363765970804*x[14]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-1.81943018149107, 0.31177572374159, 0.917145607057463, 0.380363765970804}, {}, {}, {}, {}};
	}
private:
	F11generalizedBatchReactor(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F12generalizedBatchReactor : public Expression {
public:
	static std::unique_ptr<F12generalizedBatchReactor> create() {
		Adjacency adj{{12, 13, 14}, {0}, {}};
		AdjacencyDiff adjDiff{{}, {{0, 12}, {0, 13}, {0, 14}}, {}, {}, {}, {}};
		return std::unique_ptr<F12generalizedBatchReactor>(new F12generalizedBatchReactor(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return -3.15083535140159*u[0]*x[12] + 0.607927627506543*u[0]*x[13] + 0.238500319786097*u[0]*x[14];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-3.15083535140159*u[0], 0.607927627506543*u[0], 0.238500319786097*u[0]}, {-3.15083535140159*x[12] + 0.607927627506543*x[13] + 0.238500319786097*x[14]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-3.15083535140159, 0.607927627506543, 0.238500319786097}, {}, {}, {}, {}};
	}
private:
	F12generalizedBatchReactor(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F13generalizedBatchReactor : public Expression {
public:
	static std::unique_ptr<F13generalizedBatchReactor> create() {
		Adjacency adj{{13, 14}, {0}, {}};
		AdjacencyDiff adjDiff{{}, {{0, 13}, {0, 14}}, {}, {}, {}, {}};
		return std::unique_ptr<F13generalizedBatchReactor>(new F13generalizedBatchReactor(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return -2.9146908071189*u[0]*x[13] + 0.85772698137391*u[0]*x[14];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-2.9146908071189*u[0], 0.85772698137391*u[0]}, {-2.9146908071189*x[13] + 0.85772698137391*x[14]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-2.9146908071189, 0.85772698137391}, {}, {}, {}, {}};
	}
private:
	F13generalizedBatchReactor(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F14generalizedBatchReactor : public Expression {
public:
	static std::unique_ptr<F14generalizedBatchReactor> create() {
		Adjacency adj{{0, 14}, {0}, {}};
		AdjacencyDiff adjDiff{{}, {{0, 0}, {0, 14}}, {{0, 0}}, {}, {}, {}};
		return std::unique_ptr<F14generalizedBatchReactor>(new F14generalizedBatchReactor(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return -2.94869869892561*u[0]*x[14] + pow(u[0], 2)*x[0]*DEPLETION_COEFF_VALUE;
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{pow(u[0], 2)*DEPLETION_COEFF_VALUE, -2.94869869892561*u[0]}, {-2.94869869892561*x[14] + 2*u[0]*x[0]*DEPLETION_COEFF_VALUE}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        const double x0 = 2*DEPLETION_COEFF_VALUE;
        const double x1 = x0*u[0];
		return {std::vector<double>{}, {x1, -2.94869869892561}, {x0*x[0]}, {}, {}, {}};
	}
private:
	F14generalizedBatchReactor(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F15generalizedBatchReactor : public Expression {
public:
	static std::unique_ptr<F15generalizedBatchReactor> create() {
		Adjacency adj{{}, {0}, {}};
		AdjacencyDiff adjDiff{{}, {}, {{0, 0}}, {}, {}, {}};
		return std::unique_ptr<F15generalizedBatchReactor>(new F15generalizedBatchReactor(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return pow(u[0], EXPONENT_ENERGY_VALUE);
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {pow(u[0], -1 + EXPONENT_ENERGY_VALUE)*EXPONENT_ENERGY_VALUE}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {}, {pow(u[0], -2 + EXPONENT_ENERGY_VALUE)*(-1 + EXPONENT_ENERGY_VALUE)*EXPONENT_ENERGY_VALUE}, {}, {}, {}};
	}
private:
	F15generalizedBatchReactor(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


// final constraints
class R0generalizedBatchReactor : public Constraint {
public:
	static std::unique_ptr<R0generalizedBatchReactor> create() {
		Adjacency adj{{15}, {}, {}};
		AdjacencyDiff adjDiff{{}, {}, {}, {}, {}, {}};
		return std::unique_ptr<R0generalizedBatchReactor>(new R0generalizedBatchReactor(std::move(adj), std::move(adjDiff), MINUS_INFINITY, 0.5));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return x[15];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{1}, {}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {}, {}, {}, {}, {}};
	}
private:
	R0generalizedBatchReactor(Adjacency adj, AdjacencyDiff adjDiff, double lb, double ub) : Constraint(std::move(adj), std::move(adjDiff), lb, ub) {}
};


std::vector<double> uInitialGuess(double t) {
	 return {1.5 - t};
};

Problem createProblem_generalizedBatchReactor() {

    std::vector<std::unique_ptr<Expression>> F;
    F.push_back(F0generalizedBatchReactor::create());
    F.push_back(F1generalizedBatchReactor::create());
    F.push_back(F2generalizedBatchReactor::create());
    F.push_back(F3generalizedBatchReactor::create());
    F.push_back(F4generalizedBatchReactor::create());
    F.push_back(F5generalizedBatchReactor::create());
    F.push_back(F6generalizedBatchReactor::create());
    F.push_back(F7generalizedBatchReactor::create());
    F.push_back(F8generalizedBatchReactor::create());
    F.push_back(F9generalizedBatchReactor::create());
    F.push_back(F10generalizedBatchReactor::create());
    F.push_back(F11generalizedBatchReactor::create());
    F.push_back(F12generalizedBatchReactor::create());
    F.push_back(F13generalizedBatchReactor::create());
    F.push_back(F14generalizedBatchReactor::create());
    F.push_back(F15generalizedBatchReactor::create());
    
    std::vector<std::unique_ptr<Constraint>> G;
    
    
    std::vector<std::unique_ptr<Constraint>> R;
    R.push_back(R0generalizedBatchReactor::create());
    
    std::vector<std::unique_ptr<ParamConstraint>> A;
    

    Problem problem(
            16, 1, 0,  // #vars
            {0, 0.0714285714285714, 0.0714285714285714, 0.0714285714285714, 0.0714285714285714, 0.0714285714285714, 0.0714285714285714, 0.0714285714285714, 0.0714285714285714, 0.0714285714285714, 0.0714285714285714, 0.0714285714285714, 0.0714285714285714, 0.0714285714285714, 0.0714285714285714, 0},  // x0
            {MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY},  // lb x
            {PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY},  // ub x
            &uInitialGuess,  // u0 initial guesses for optimization
            {0},  // lb u
            {5},  // ub u
            {},  // p0 initial guesses for optimization
            {},  // lb p
            {},  // ub p
            MayergeneralizedBatchReactor::create(),
            {},
            std::move(F),
            std::move(G),
            std::move(R),
            std::move(A),
            "generalizedBatchReactor");
    return problem;
};

int main() {
    auto problem = std::make_shared<const Problem>(createProblem_generalizedBatchReactor());
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
        
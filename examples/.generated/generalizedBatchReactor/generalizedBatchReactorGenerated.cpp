
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
		return 0.597544669755805*u[0]*x[1] + 0.0539728740676567*u[0]*x[2] + 0.363897340580728*u[0]*x[3] + 0.979879288221738*u[0]*x[4] + 0.442227952315806*u[0]*x[5] - pow(u[0], 2)*x[0]*DEPLETION_COEFF_VALUE;
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-pow(u[0], 2)*DEPLETION_COEFF_VALUE, 0.597544669755805*u[0], 0.0539728740676567*u[0], 0.363897340580728*u[0], 0.979879288221738*u[0], 0.442227952315806*u[0]}, {0.597544669755805*x[1] + 0.0539728740676567*x[2] + 0.363897340580728*x[3] + 0.979879288221738*x[4] + 0.442227952315806*x[5] - 2*u[0]*x[0]*DEPLETION_COEFF_VALUE}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        const double x0 = 2*DEPLETION_COEFF_VALUE;
        const double x1 = -x0*u[0];
		return {std::vector<double>{}, {x1, 0.597544669755805, 0.0539728740676567, 0.363897340580728, 0.979879288221738, 0.442227952315806}, {-x0*x[0]}, {}, {}, {}};
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
		return -0.597544669755805*u[0]*x[1] + 0.13828443233923*u[0]*x[2] + 0.945391323401778*u[0]*x[3] + 0.683980514484315*u[0]*x[4] + 0.316182218750967*u[0]*x[5] + 0.412292622086736*u[0]*x[6];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-0.597544669755805*u[0], 0.13828443233923*u[0], 0.945391323401778*u[0], 0.683980514484315*u[0], 0.316182218750967*u[0], 0.412292622086736*u[0]}, {-0.597544669755805*x[1] + 0.13828443233923*x[2] + 0.945391323401778*x[3] + 0.683980514484315*x[4] + 0.316182218750967*x[5] + 0.412292622086736*x[6]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-0.597544669755805, 0.13828443233923, 0.945391323401778, 0.683980514484315, 0.316182218750967, 0.412292622086736}, {}, {}, {}, {}};
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
		return -0.192257306406886*u[0]*x[2] + 0.159558059433433*u[0]*x[3] + 0.749357552207486*u[0]*x[4] + 0.830644009629674*u[0]*x[5] + 0.104642955039831*u[0]*x[6] + 0.238479139550985*u[0]*x[7];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-0.192257306406886*u[0], 0.159558059433433*u[0], 0.749357552207486*u[0], 0.830644009629674*u[0], 0.104642955039831*u[0], 0.238479139550985*u[0]}, {-0.192257306406886*x[2] + 0.159558059433433*x[3] + 0.749357552207486*x[4] + 0.830644009629674*x[5] + 0.104642955039831*x[6] + 0.238479139550985*x[7]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-0.192257306406886, 0.159558059433433, 0.749357552207486, 0.830644009629674, 0.104642955039831, 0.238479139550985}, {}, {}, {}, {}};
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
		return -1.46884672341594*u[0]*x[3] + 0.291499791657283*u[0]*x[4] + 0.666889673662558*u[0]*x[5] + 0.627081489598424*u[0]*x[6] + 0.663246254989886*u[0]*x[7] + 0.53042021884027*u[0]*x[8];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-1.46884672341594*u[0], 0.291499791657283*u[0], 0.666889673662558*u[0], 0.627081489598424*u[0], 0.663246254989886*u[0], 0.53042021884027*u[0]}, {-1.46884672341594*x[3] + 0.291499791657283*x[4] + 0.666889673662558*x[5] + 0.627081489598424*x[6] + 0.663246254989886*x[7] + 0.53042021884027*x[8]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-1.46884672341594, 0.291499791657283, 0.666889673662558, 0.627081489598424, 0.663246254989886, 0.53042021884027}, {}, {}, {}, {}};
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
		return -2.70471714657082*u[0]*x[4] + 0.955816677980712*u[0]*x[5] + 0.397781229192407*u[0]*x[6] + 0.0914480644083759*u[0]*x[7] + 0.797846461495176*u[0]*x[8] + 0.988387833198566*u[0]*x[9];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-2.70471714657082*u[0], 0.955816677980712*u[0], 0.397781229192407*u[0], 0.0914480644083759*u[0], 0.797846461495176*u[0], 0.988387833198566*u[0]}, {-2.70471714657082*x[4] + 0.955816677980712*x[5] + 0.397781229192407*x[6] + 0.0914480644083759*x[7] + 0.797846461495176*x[8] + 0.988387833198566*x[9]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-2.70471714657082, 0.955816677980712, 0.397781229192407, 0.0914480644083759, 0.797846461495176, 0.988387833198566}, {}, {}, {}, {}};
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
		return 0.383018941154278*u[0]*x[10] - 3.21176053233972*u[0]*x[5] + 0.118977221524535*u[0]*x[6] + 0.976804710229121*u[0]*x[7] + 0.378305193315893*u[0]*x[8] + 0.962345716497645*u[0]*x[9];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-3.21176053233972*u[0], 0.118977221524535*u[0], 0.976804710229121*u[0], 0.378305193315893*u[0], 0.962345716497645*u[0], 0.383018941154278*u[0]}, {0.383018941154278*x[10] - 3.21176053233972*x[5] + 0.118977221524535*x[6] + 0.976804710229121*x[7] + 0.378305193315893*x[8] + 0.962345716497645*x[9]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-3.21176053233972, 0.118977221524535, 0.976804710229121, 0.378305193315893, 0.962345716497645, 0.383018941154278}, {}, {}, {}, {}};
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
		return 0.206764362053459*u[0]*x[10] + 0.678749653703692*u[0]*x[11] - 1.66077551744193*u[0]*x[6] + 0.723137583375984*u[0]*x[7] + 0.428296556292115*u[0]*x[8] + 0.767954713905873*u[0]*x[9];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-1.66077551744193*u[0], 0.723137583375984*u[0], 0.428296556292115*u[0], 0.767954713905873*u[0], 0.206764362053459*u[0], 0.678749653703692*u[0]}, {0.206764362053459*x[10] + 0.678749653703692*x[11] - 1.66077551744193*x[6] + 0.723137583375984*x[7] + 0.428296556292115*x[8] + 0.767954713905873*x[9]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-1.66077551744193, 0.723137583375984, 0.428296556292115, 0.767954713905873, 0.206764362053459, 0.678749653703692}, {}, {}, {}, {}};
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
		return 0.658281413349725*u[0]*x[10] + 0.858938910650542*u[0]*x[11] + 0.376833295326455*u[0]*x[12] - 2.69311575255435*u[0]*x[7] + 0.397399352762824*u[0]*x[8] + 0.440688857475468*u[0]*x[9];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-2.69311575255435*u[0], 0.397399352762824*u[0], 0.440688857475468*u[0], 0.658281413349725*u[0], 0.858938910650542*u[0], 0.376833295326455*u[0]}, {0.658281413349725*x[10] + 0.858938910650542*x[11] + 0.376833295326455*x[12] - 2.69311575255435*x[7] + 0.397399352762824*x[8] + 0.440688857475468*x[9]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-2.69311575255435, 0.397399352762824, 0.440688857475468, 0.658281413349725, 0.858938910650542, 0.376833295326455}, {}, {}, {}, {}};
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
		return 0.399294370890498*u[0]*x[10] + 0.870353049376835*u[0]*x[11] + 0.216112863243079*u[0]*x[12] + 0.695174275914929*u[0]*x[13] - 2.53226778270628*u[0]*x[8] + 0.922595191897896*u[0]*x[9];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-2.53226778270628*u[0], 0.922595191897896*u[0], 0.399294370890498*u[0], 0.870353049376835*u[0], 0.216112863243079*u[0], 0.695174275914929*u[0]}, {0.399294370890498*x[10] + 0.870353049376835*x[11] + 0.216112863243079*x[12] + 0.695174275914929*x[13] - 2.53226778270628*x[8] + 0.922595191897896*x[9]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-2.53226778270628, 0.922595191897896, 0.399294370890498, 0.870353049376835, 0.216112863243079, 0.695174275914929}, {}, {}, {}, {}};
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
		return 0.180925108233109*u[0]*x[10] + 0.0306010916517576*u[0]*x[11] + 0.384187340891354*u[0]*x[12] + 0.647864059948507*u[0]*x[13] + 0.891478392601203*u[0]*x[14] - 4.08197231297545*u[0]*x[9];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-4.08197231297545*u[0], 0.180925108233109*u[0], 0.0306010916517576*u[0], 0.384187340891354*u[0], 0.647864059948507*u[0], 0.891478392601203*u[0]}, {0.180925108233109*x[10] + 0.0306010916517576*x[11] + 0.384187340891354*x[12] + 0.647864059948507*x[13] + 0.891478392601203*x[14] - 4.08197231297545*x[9]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-4.08197231297545, 0.180925108233109, 0.0306010916517576, 0.384187340891354, 0.647864059948507, 0.891478392601203}, {}, {}, {}, {}};
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
		return -1.82828419568107*u[0]*x[10] + 0.476186888977727*u[0]*x[11] + 0.853574240067793*u[0]*x[12] + 0.769646747361565*u[0]*x[13] + 0.00434551243564307*u[0]*x[14];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-1.82828419568107*u[0], 0.476186888977727*u[0], 0.853574240067793*u[0], 0.769646747361565*u[0], 0.00434551243564307*u[0]}, {-1.82828419568107*x[10] + 0.476186888977727*x[11] + 0.853574240067793*x[12] + 0.769646747361565*x[13] + 0.00434551243564307*x[14]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-1.82828419568107, 0.476186888977727, 0.853574240067793, 0.769646747361565, 0.00434551243564307}, {}, {}, {}, {}};
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
		return -2.91482959436055*u[0]*x[11] + 0.0773441802267341*u[0]*x[12] + 0.720638601429823*u[0]*x[13] + 0.893965265926955*u[0]*x[14];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-2.91482959436055*u[0], 0.0773441802267341*u[0], 0.720638601429823*u[0], 0.893965265926955*u[0]}, {-2.91482959436055*x[11] + 0.0773441802267341*x[12] + 0.720638601429823*x[13] + 0.893965265926955*x[14]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-2.91482959436055, 0.0773441802267341, 0.720638601429823, 0.893965265926955}, {}, {}, {}, {}};
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
		return -1.90805191975541*u[0]*x[12] + 0.180309799136384*u[0]*x[13] + 0.521795675868917*u[0]*x[14];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-1.90805191975541*u[0], 0.180309799136384*u[0], 0.521795675868917*u[0]}, {-1.90805191975541*x[12] + 0.180309799136384*x[13] + 0.521795675868917*x[14]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-1.90805191975541, 0.180309799136384, 0.521795675868917}, {}, {}, {}, {}};
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
		return -3.01363348379121*u[0]*x[13] + 0.111921111187218*u[0]*x[14];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-3.01363348379121*u[0], 0.111921111187218*u[0]}, {-3.01363348379121*x[13] + 0.111921111187218*x[14]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-3.01363348379121, 0.111921111187218}, {}, {}, {}, {}};
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
		return -2.42350595801994*u[0]*x[14] + pow(u[0], 2)*x[0]*DEPLETION_COEFF_VALUE;
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{pow(u[0], 2)*DEPLETION_COEFF_VALUE, -2.42350595801994*u[0]}, {-2.42350595801994*x[14] + 2*u[0]*x[0]*DEPLETION_COEFF_VALUE}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        const double x0 = 2*DEPLETION_COEFF_VALUE;
        const double x1 = x0*u[0];
		return {std::vector<double>{}, {x1, -2.42350595801994}, {x0*x[0]}, {}, {}, {}};
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
            {1},  // u0 initial guesses for optimization
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
        
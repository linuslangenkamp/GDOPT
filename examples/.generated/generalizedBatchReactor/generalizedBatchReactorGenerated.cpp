
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
		return 0.0555486799671117*u[0]*x[1] + 0.82077766243196*u[0]*x[2] + 0.544332040896934*u[0]*x[3] + 0.364633248507526*u[0]*x[4] + 0.501779570968709*u[0]*x[5] - pow(u[0], 2)*x[0]*DEPLETION_COEFF_VALUE;
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-pow(u[0], 2)*DEPLETION_COEFF_VALUE, 0.0555486799671117*u[0], 0.82077766243196*u[0], 0.544332040896934*u[0], 0.364633248507526*u[0], 0.501779570968709*u[0]}, {0.0555486799671117*x[1] + 0.82077766243196*x[2] + 0.544332040896934*x[3] + 0.364633248507526*x[4] + 0.501779570968709*x[5] - 2*u[0]*x[0]*DEPLETION_COEFF_VALUE}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        const double x0 = 2*DEPLETION_COEFF_VALUE;
        const double x1 = -x0*u[0];
		return {std::vector<double>{}, {x1, 0.0555486799671117, 0.82077766243196, 0.544332040896934, 0.364633248507526, 0.501779570968709}, {-x0*x[0]}, {}, {}, {}};
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
		return -0.0555486799671117*u[0]*x[1] + 0.523930148921627*u[0]*x[2] + 0.627396493238211*u[0]*x[3] + 0.592629537669276*u[0]*x[4] + 0.260345009513549*u[0]*x[5] + 0.417762866524527*u[0]*x[6];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-0.0555486799671117*u[0], 0.523930148921627*u[0], 0.627396493238211*u[0], 0.592629537669276*u[0], 0.260345009513549*u[0], 0.417762866524527*u[0]}, {-0.0555486799671117*x[1] + 0.523930148921627*x[2] + 0.627396493238211*x[3] + 0.592629537669276*x[4] + 0.260345009513549*x[5] + 0.417762866524527*x[6]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-0.0555486799671117, 0.523930148921627, 0.627396493238211, 0.592629537669276, 0.260345009513549, 0.417762866524527}, {}, {}, {}, {}};
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
		return -1.34470781135359*u[0]*x[2] + 0.35521352345907*u[0]*x[3] + 0.551515730487981*u[0]*x[4] + 0.327289601374793*u[0]*x[5] + 0.721705963650284*u[0]*x[6] + 0.544712200765281*u[0]*x[7];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-1.34470781135359*u[0], 0.35521352345907*u[0], 0.551515730487981*u[0], 0.327289601374793*u[0], 0.721705963650284*u[0], 0.544712200765281*u[0]}, {-1.34470781135359*x[2] + 0.35521352345907*x[3] + 0.551515730487981*x[4] + 0.327289601374793*x[5] + 0.721705963650284*x[6] + 0.544712200765281*x[7]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-1.34470781135359, 0.35521352345907, 0.551515730487981, 0.327289601374793, 0.721705963650284, 0.544712200765281}, {}, {}, {}, {}};
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
		return -1.52694205759422*u[0]*x[3] + 0.784272010532514*u[0]*x[4] + 0.965493934452077*u[0]*x[5] + 0.75736235430309*u[0]*x[6] + 0.23079217859386*u[0]*x[7] + 0.316949318593201*u[0]*x[8];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-1.52694205759422*u[0], 0.784272010532514*u[0], 0.965493934452077*u[0], 0.75736235430309*u[0], 0.23079217859386*u[0], 0.316949318593201*u[0]}, {-1.52694205759422*x[3] + 0.784272010532514*x[4] + 0.965493934452077*x[5] + 0.75736235430309*x[6] + 0.23079217859386*x[7] + 0.316949318593201*x[8]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-1.52694205759422, 0.784272010532514, 0.965493934452077, 0.75736235430309, 0.23079217859386, 0.316949318593201}, {}, {}, {}, {}};
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
		return -2.2930505271973*u[0]*x[4] + 0.203522627779163*u[0]*x[5] + 0.761206059747744*u[0]*x[6] + 0.00493147540516248*u[0]*x[7] + 0.141516532940813*u[0]*x[8] + 0.808315375885871*u[0]*x[9];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-2.2930505271973*u[0], 0.203522627779163*u[0], 0.761206059747744*u[0], 0.00493147540516248*u[0], 0.141516532940813*u[0], 0.808315375885871*u[0]}, {-2.2930505271973*x[4] + 0.203522627779163*x[5] + 0.761206059747744*x[6] + 0.00493147540516248*x[7] + 0.141516532940813*x[8] + 0.808315375885871*x[9]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-2.2930505271973, 0.203522627779163, 0.761206059747744, 0.00493147540516248, 0.141516532940813, 0.808315375885871}, {}, {}, {}, {}};
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
		return 0.749380594554308*u[0]*x[10] - 2.25843074408829*u[0]*x[5] + 0.55766819497798*u[0]*x[6] + 0.142334867787304*u[0]*x[7] + 0.338685301366266*u[0]*x[8] + 0.694640752786613*u[0]*x[9];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-2.25843074408829*u[0], 0.55766819497798*u[0], 0.142334867787304*u[0], 0.338685301366266*u[0], 0.694640752786613*u[0], 0.749380594554308*u[0]}, {0.749380594554308*x[10] - 2.25843074408829*x[5] + 0.55766819497798*x[6] + 0.142334867787304*x[7] + 0.338685301366266*x[8] + 0.694640752786613*x[9]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-2.25843074408829, 0.55766819497798, 0.142334867787304, 0.338685301366266, 0.694640752786613, 0.749380594554308}, {}, {}, {}, {}};
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
		return 0.0557107760764743*u[0]*x[10] + 0.85721725724932*u[0]*x[11] - 3.21570543920362*u[0]*x[6] + 0.592377366390014*u[0]*x[7] + 0.762824427499181*u[0]*x[8] + 0.89796445945797*u[0]*x[9];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-3.21570543920362*u[0], 0.592377366390014*u[0], 0.762824427499181*u[0], 0.89796445945797*u[0], 0.0557107760764743*u[0], 0.85721725724932*u[0]}, {0.0557107760764743*x[10] + 0.85721725724932*x[11] - 3.21570543920362*x[6] + 0.592377366390014*x[7] + 0.762824427499181*x[8] + 0.89796445945797*x[9]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-3.21570543920362, 0.592377366390014, 0.762824427499181, 0.89796445945797, 0.0557107760764743, 0.85721725724932}, {}, {}, {}, {}};
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
		return 0.837712348374487*u[0]*x[10] + 0.883131715103596*u[0]*x[11] + 0.90397672001838*u[0]*x[12] - 1.51514808894162*u[0]*x[7] + 0.70826613492235*u[0]*x[8] + 0.770357053802267*u[0]*x[9];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-1.51514808894162*u[0], 0.70826613492235*u[0], 0.770357053802267*u[0], 0.837712348374487*u[0], 0.883131715103596*u[0], 0.90397672001838*u[0]}, {0.837712348374487*x[10] + 0.883131715103596*x[11] + 0.90397672001838*x[12] - 1.51514808894162*x[7] + 0.70826613492235*x[8] + 0.770357053802267*x[9]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-1.51514808894162, 0.70826613492235, 0.770357053802267, 0.837712348374487, 0.883131715103596, 0.90397672001838}, {}, {}, {}, {}};
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
		return 0.558517897457815*u[0]*x[10] + 0.700947463516176*u[0]*x[11] + 0.567441873181581*u[0]*x[12] + 0.102845726663574*u[0]*x[13] - 2.26824171532181*u[0]*x[8] + 0.968001262486662*u[0]*x[9];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-2.26824171532181*u[0], 0.968001262486662*u[0], 0.558517897457815*u[0], 0.700947463516176*u[0], 0.567441873181581*u[0], 0.102845726663574*u[0]}, {0.558517897457815*x[10] + 0.700947463516176*x[11] + 0.567441873181581*x[12] + 0.102845726663574*x[13] - 2.26824171532181*x[8] + 0.968001262486662*x[9]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-2.26824171532181, 0.968001262486662, 0.558517897457815, 0.700947463516176, 0.567441873181581, 0.102845726663574}, {}, {}, {}, {}};
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
		return 0.972993214046657*u[0]*x[10] + 0.0215362970675818*u[0]*x[11] + 0.586334772501055*u[0]*x[12] + 0.34336583980797*u[0]*x[13] + 0.920360371807376*u[0]*x[14] - 4.13927890441938*u[0]*x[9];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-4.13927890441938*u[0], 0.972993214046657*u[0], 0.0215362970675818*u[0], 0.586334772501055*u[0], 0.34336583980797*u[0], 0.920360371807376*u[0]}, {0.972993214046657*x[10] + 0.0215362970675818*x[11] + 0.586334772501055*x[12] + 0.34336583980797*x[13] + 0.920360371807376*x[14] - 4.13927890441938*x[9]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-4.13927890441938, 0.972993214046657, 0.0215362970675818, 0.586334772501055, 0.34336583980797, 0.920360371807376}, {}, {}, {}, {}};
	}
private:
	F9generalizedBatchReactor(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F10generalizedBatchReactor : public Expression {
public:
	static std::unique_ptr<F10generalizedBatchReactor> create() {
		Adjacency adj{{10, 11, 12, 13, 14, 15}, {0}, {}};
		AdjacencyDiff adjDiff{{}, {{0, 10}, {0, 11}, {0, 12}, {0, 13}, {0, 14}, {0, 15}}, {}, {}, {}, {}};
		return std::unique_ptr<F10generalizedBatchReactor>(new F10generalizedBatchReactor(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return -3.17431483050974*u[0]*x[10] + 0.364111487826355*u[0]*x[11] + 0.425427003889252*u[0]*x[12] + 0.0682597591389495*u[0]*x[13] + 0.445411765224119*u[0]*x[14] + 0.23495404593911*u[0]*x[15];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-3.17431483050974*u[0], 0.364111487826355*u[0], 0.425427003889252*u[0], 0.0682597591389495*u[0], 0.445411765224119*u[0], 0.23495404593911*u[0]}, {-3.17431483050974*x[10] + 0.364111487826355*x[11] + 0.425427003889252*x[12] + 0.0682597591389495*x[13] + 0.445411765224119*x[14] + 0.23495404593911*x[15]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-3.17431483050974, 0.364111487826355, 0.425427003889252, 0.0682597591389495, 0.445411765224119, 0.23495404593911}, {}, {}, {}, {}};
	}
private:
	F10generalizedBatchReactor(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F11generalizedBatchReactor : public Expression {
public:
	static std::unique_ptr<F11generalizedBatchReactor> create() {
		Adjacency adj{{11, 12, 13, 14, 15, 16}, {0}, {}};
		AdjacencyDiff adjDiff{{}, {{0, 11}, {0, 12}, {0, 13}, {0, 14}, {0, 15}, {0, 16}}, {}, {}, {}, {}};
		return std::unique_ptr<F11generalizedBatchReactor>(new F11generalizedBatchReactor(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return -2.82694422076303*u[0]*x[11] + 0.465995384617403*u[0]*x[12] + 0.160635419265447*u[0]*x[13] + 0.382627314486551*u[0]*x[14] + 0.086335577081366*u[0]*x[15] + 0.834784630014535*u[0]*x[16];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-2.82694422076303*u[0], 0.465995384617403*u[0], 0.160635419265447*u[0], 0.382627314486551*u[0], 0.086335577081366*u[0], 0.834784630014535*u[0]}, {-2.82694422076303*x[11] + 0.465995384617403*x[12] + 0.160635419265447*x[13] + 0.382627314486551*x[14] + 0.086335577081366*x[15] + 0.834784630014535*x[16]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-2.82694422076303, 0.465995384617403, 0.160635419265447, 0.382627314486551, 0.086335577081366, 0.834784630014535}, {}, {}, {}, {}};
	}
private:
	F11generalizedBatchReactor(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F12generalizedBatchReactor : public Expression {
public:
	static std::unique_ptr<F12generalizedBatchReactor> create() {
		Adjacency adj{{12, 13, 14, 15, 16, 17}, {0}, {}};
		AdjacencyDiff adjDiff{{}, {{0, 12}, {0, 13}, {0, 14}, {0, 15}, {0, 16}, {0, 17}}, {}, {}, {}, {}};
		return std::unique_ptr<F12generalizedBatchReactor>(new F12generalizedBatchReactor(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return -2.94917575420767*u[0]*x[12] + 0.0557490381233168*u[0]*x[13] + 0.470027445959506*u[0]*x[14] + 0.325371636792972*u[0]*x[15] + 0.981938832299811*u[0]*x[16] + 0.205097837509714*u[0]*x[17];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-2.94917575420767*u[0], 0.0557490381233168*u[0], 0.470027445959506*u[0], 0.325371636792972*u[0], 0.981938832299811*u[0], 0.205097837509714*u[0]}, {-2.94917575420767*x[12] + 0.0557490381233168*x[13] + 0.470027445959506*x[14] + 0.325371636792972*x[15] + 0.981938832299811*x[16] + 0.205097837509714*x[17]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-2.94917575420767, 0.0557490381233168, 0.470027445959506, 0.325371636792972, 0.981938832299811, 0.205097837509714}, {}, {}, {}, {}};
	}
private:
	F12generalizedBatchReactor(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F13generalizedBatchReactor : public Expression {
public:
	static std::unique_ptr<F13generalizedBatchReactor> create() {
		Adjacency adj{{13, 14, 15, 16, 17, 18}, {0}, {}};
		AdjacencyDiff adjDiff{{}, {{0, 13}, {0, 14}, {0, 15}, {0, 16}, {0, 17}, {0, 18}}, {}, {}, {}, {}};
		return std::unique_ptr<F13generalizedBatchReactor>(new F13generalizedBatchReactor(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return -0.730855782999258*u[0]*x[13] + 0.116559511682509*u[0]*x[14] + 0.348861322371123*u[0]*x[15] + 0.800824992841427*u[0]*x[16] + 0.0863137579021616*u[0]*x[17] + 0.0446177497615747*u[0]*x[18];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-0.730855782999258*u[0], 0.116559511682509*u[0], 0.348861322371123*u[0], 0.800824992841427*u[0], 0.0863137579021616*u[0], 0.0446177497615747*u[0]}, {-0.730855782999258*x[13] + 0.116559511682509*x[14] + 0.348861322371123*x[15] + 0.800824992841427*x[16] + 0.0863137579021616*x[17] + 0.0446177497615747*x[18]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-0.730855782999258, 0.116559511682509, 0.348861322371123, 0.800824992841427, 0.0863137579021616, 0.0446177497615747}, {}, {}, {}, {}};
	}
private:
	F13generalizedBatchReactor(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F14generalizedBatchReactor : public Expression {
public:
	static std::unique_ptr<F14generalizedBatchReactor> create() {
		Adjacency adj{{14, 15, 16, 17, 18, 19}, {0}, {}};
		AdjacencyDiff adjDiff{{}, {{0, 14}, {0, 15}, {0, 16}, {0, 17}, {0, 18}, {0, 19}}, {}, {}, {}, {}};
		return std::unique_ptr<F14generalizedBatchReactor>(new F14generalizedBatchReactor(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return -2.33498640916006*u[0]*x[14] + 0.554588839951002*u[0]*x[15] + 0.166567442239406*u[0]*x[16] + 0.591351054795758*u[0]*x[17] + 0.00796201210813197*u[0]*x[18] + 0.449894181692955*u[0]*x[19];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-2.33498640916006*u[0], 0.554588839951002*u[0], 0.166567442239406*u[0], 0.591351054795758*u[0], 0.00796201210813197*u[0], 0.449894181692955*u[0]}, {-2.33498640916006*x[14] + 0.554588839951002*x[15] + 0.166567442239406*x[16] + 0.591351054795758*x[17] + 0.00796201210813197*x[18] + 0.449894181692955*x[19]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-2.33498640916006, 0.554588839951002, 0.166567442239406, 0.591351054795758, 0.00796201210813197, 0.449894181692955}, {}, {}, {}, {}};
	}
private:
	F14generalizedBatchReactor(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F15generalizedBatchReactor : public Expression {
public:
	static std::unique_ptr<F15generalizedBatchReactor> create() {
		Adjacency adj{{15, 16, 17, 18, 19, 20}, {0}, {}};
		AdjacencyDiff adjDiff{{}, {{0, 15}, {0, 16}, {0, 17}, {0, 18}, {0, 19}, {0, 20}}, {}, {}, {}, {}};
		return std::unique_ptr<F15generalizedBatchReactor>(new F15generalizedBatchReactor(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return -1.55011142213557*u[0]*x[15] + 0.381720065686155*u[0]*x[16] + 0.801812996241705*u[0]*x[17] + 0.455956437816518*u[0]*x[18] + 0.560562494515818*u[0]*x[19] + 0.698790547794545*u[0]*x[20];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-1.55011142213557*u[0], 0.381720065686155*u[0], 0.801812996241705*u[0], 0.455956437816518*u[0], 0.560562494515818*u[0], 0.698790547794545*u[0]}, {-1.55011142213557*x[15] + 0.381720065686155*x[16] + 0.801812996241705*x[17] + 0.455956437816518*x[18] + 0.560562494515818*x[19] + 0.698790547794545*x[20]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-1.55011142213557, 0.381720065686155, 0.801812996241705, 0.455956437816518, 0.560562494515818, 0.698790547794545}, {}, {}, {}, {}};
	}
private:
	F15generalizedBatchReactor(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F16generalizedBatchReactor : public Expression {
public:
	static std::unique_ptr<F16generalizedBatchReactor> create() {
		Adjacency adj{{16, 17, 18, 19, 20, 21}, {0}, {}};
		AdjacencyDiff adjDiff{{}, {{0, 16}, {0, 17}, {0, 18}, {0, 19}, {0, 20}, {0, 21}}, {}, {}, {}, {}};
		return std::unique_ptr<F16generalizedBatchReactor>(new F16generalizedBatchReactor(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return -3.16583596308133*u[0]*x[16] + 0.5669782879309*u[0]*x[17] + 0.101710949998167*u[0]*x[18] + 0.996569128665827*u[0]*x[19] + 0.348159743141782*u[0]*x[20] + 0.589604060335158*u[0]*x[21];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-3.16583596308133*u[0], 0.5669782879309*u[0], 0.101710949998167*u[0], 0.996569128665827*u[0], 0.348159743141782*u[0], 0.589604060335158*u[0]}, {-3.16583596308133*x[16] + 0.5669782879309*x[17] + 0.101710949998167*x[18] + 0.996569128665827*x[19] + 0.348159743141782*x[20] + 0.589604060335158*x[21]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-3.16583596308133, 0.5669782879309, 0.101710949998167, 0.996569128665827, 0.348159743141782, 0.589604060335158}, {}, {}, {}, {}};
	}
private:
	F16generalizedBatchReactor(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F17generalizedBatchReactor : public Expression {
public:
	static std::unique_ptr<F17generalizedBatchReactor> create() {
		Adjacency adj{{17, 18, 19, 20, 21, 22}, {0}, {}};
		AdjacencyDiff adjDiff{{}, {{0, 17}, {0, 18}, {0, 19}, {0, 20}, {0, 21}, {0, 22}}, {}, {}, {}, {}};
		return std::unique_ptr<F17generalizedBatchReactor>(new F17generalizedBatchReactor(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return -2.25155393438024*u[0]*x[17] + 0.935609604750318*u[0]*x[18] + 0.770663908354191*u[0]*x[19] + 0.534093066599977*u[0]*x[20] + 0.964998250647389*u[0]*x[21] + 0.13896039650688*u[0]*x[22];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-2.25155393438024*u[0], 0.935609604750318*u[0], 0.770663908354191*u[0], 0.534093066599977*u[0], 0.964998250647389*u[0], 0.13896039650688*u[0]}, {-2.25155393438024*x[17] + 0.935609604750318*x[18] + 0.770663908354191*x[19] + 0.534093066599977*x[20] + 0.964998250647389*x[21] + 0.13896039650688*x[22]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-2.25155393438024, 0.935609604750318, 0.770663908354191, 0.534093066599977, 0.964998250647389, 0.13896039650688}, {}, {}, {}, {}};
	}
private:
	F17generalizedBatchReactor(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F18generalizedBatchReactor : public Expression {
public:
	static std::unique_ptr<F18generalizedBatchReactor> create() {
		Adjacency adj{{18, 19, 20, 21, 22, 23}, {0}, {}};
		AdjacencyDiff adjDiff{{}, {{0, 18}, {0, 19}, {0, 20}, {0, 21}, {0, 22}, {0, 23}}, {}, {}, {}, {}};
		return std::unique_ptr<F18generalizedBatchReactor>(new F18generalizedBatchReactor(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return -1.54585675443471*u[0]*x[18] + 0.610860779442975*u[0]*x[19] + 0.839738174418102*u[0]*x[20] + 0.0809196524338009*u[0]*x[21] + 0.922820647555851*u[0]*x[22] + 0.482053542960284*u[0]*x[23];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-1.54585675443471*u[0], 0.610860779442975*u[0], 0.839738174418102*u[0], 0.0809196524338009*u[0], 0.922820647555851*u[0], 0.482053542960284*u[0]}, {-1.54585675443471*x[18] + 0.610860779442975*x[19] + 0.839738174418102*x[20] + 0.0809196524338009*x[21] + 0.922820647555851*x[22] + 0.482053542960284*x[23]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-1.54585675443471, 0.610860779442975, 0.839738174418102, 0.0809196524338009, 0.922820647555851, 0.482053542960284}, {}, {}, {}, {}};
	}
private:
	F18generalizedBatchReactor(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F19generalizedBatchReactor : public Expression {
public:
	static std::unique_ptr<F19generalizedBatchReactor> create() {
		Adjacency adj{{19, 20, 21, 22, 23, 24}, {0}, {}};
		AdjacencyDiff adjDiff{{}, {{0, 19}, {0, 20}, {0, 21}, {0, 22}, {0, 23}, {0, 24}}, {}, {}, {}, {}};
		return std::unique_ptr<F19generalizedBatchReactor>(new F19generalizedBatchReactor(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return -3.38855049267177*u[0]*x[19] + 0.220673632292271*u[0]*x[20] + 0.956185736385657*u[0]*x[21] + 0.580728020120306*u[0]*x[22] + 0.627467587379673*u[0]*x[23] + 0.623358747752024*u[0]*x[24];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-3.38855049267177*u[0], 0.220673632292271*u[0], 0.956185736385657*u[0], 0.580728020120306*u[0], 0.627467587379673*u[0], 0.623358747752024*u[0]}, {-3.38855049267177*x[19] + 0.220673632292271*x[20] + 0.956185736385657*x[21] + 0.580728020120306*x[22] + 0.627467587379673*x[23] + 0.623358747752024*x[24]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-3.38855049267177, 0.220673632292271, 0.956185736385657, 0.580728020120306, 0.627467587379673, 0.623358747752024}, {}, {}, {}, {}};
	}
private:
	F19generalizedBatchReactor(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F20generalizedBatchReactor : public Expression {
public:
	static std::unique_ptr<F20generalizedBatchReactor> create() {
		Adjacency adj{{20, 21, 22, 23, 24, 25}, {0}, {}};
		AdjacencyDiff adjDiff{{}, {{0, 20}, {0, 21}, {0, 22}, {0, 23}, {0, 24}, {0, 25}}, {}, {}, {}, {}};
		return std::unique_ptr<F20generalizedBatchReactor>(new F20generalizedBatchReactor(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return -2.64145516424668*u[0]*x[20] + 0.189782161967765*u[0]*x[21] + 0.669650858320891*u[0]*x[22] + 0.660697482989425*u[0]*x[23] + 0.939023479531182*u[0]*x[24] + 0.179197694611771*u[0]*x[25];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-2.64145516424668*u[0], 0.189782161967765*u[0], 0.669650858320891*u[0], 0.660697482989425*u[0], 0.939023479531182*u[0], 0.179197694611771*u[0]}, {-2.64145516424668*x[20] + 0.189782161967765*x[21] + 0.669650858320891*x[22] + 0.660697482989425*x[23] + 0.939023479531182*x[24] + 0.179197694611771*x[25]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-2.64145516424668, 0.189782161967765, 0.669650858320891, 0.660697482989425, 0.939023479531182, 0.179197694611771}, {}, {}, {}, {}};
	}
private:
	F20generalizedBatchReactor(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F21generalizedBatchReactor : public Expression {
public:
	static std::unique_ptr<F21generalizedBatchReactor> create() {
		Adjacency adj{{21, 22, 23, 24, 25, 26}, {0}, {}};
		AdjacencyDiff adjDiff{{}, {{0, 21}, {0, 22}, {0, 23}, {0, 24}, {0, 25}, {0, 26}}, {}, {}, {}, {}};
		return std::unique_ptr<F21generalizedBatchReactor>(new F21generalizedBatchReactor(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return -2.78148986176977*u[0]*x[21] + 0.517255959376217*u[0]*x[22] + 0.71641748645842*u[0]*x[23] + 0.343315804490311*u[0]*x[24] + 0.869794720038914*u[0]*x[25] + 0.976011515338833*u[0]*x[26];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-2.78148986176977*u[0], 0.517255959376217*u[0], 0.71641748645842*u[0], 0.343315804490311*u[0], 0.869794720038914*u[0], 0.976011515338833*u[0]}, {-2.78148986176977*x[21] + 0.517255959376217*x[22] + 0.71641748645842*x[23] + 0.343315804490311*x[24] + 0.869794720038914*x[25] + 0.976011515338833*x[26]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-2.78148986176977, 0.517255959376217, 0.71641748645842, 0.343315804490311, 0.869794720038914, 0.976011515338833}, {}, {}, {}, {}};
	}
private:
	F21generalizedBatchReactor(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F22generalizedBatchReactor : public Expression {
public:
	static std::unique_ptr<F22generalizedBatchReactor> create() {
		Adjacency adj{{22, 23, 24, 25, 26, 27}, {0}, {}};
		AdjacencyDiff adjDiff{{}, {{0, 22}, {0, 23}, {0, 24}, {0, 25}, {0, 26}, {0, 27}}, {}, {}, {}, {}};
		return std::unique_ptr<F22generalizedBatchReactor>(new F22generalizedBatchReactor(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return -2.82941588188015*u[0]*x[22] + 0.201495132985613*u[0]*x[23] + 0.298017793415622*u[0]*x[24] + 0.0934947735910401*u[0]*x[25] + 0.548588115971628*u[0]*x[26] + 0.431368751868528*u[0]*x[27];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-2.82941588188015*u[0], 0.201495132985613*u[0], 0.298017793415622*u[0], 0.0934947735910401*u[0], 0.548588115971628*u[0], 0.431368751868528*u[0]}, {-2.82941588188015*x[22] + 0.201495132985613*x[23] + 0.298017793415622*x[24] + 0.0934947735910401*x[25] + 0.548588115971628*x[26] + 0.431368751868528*x[27]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-2.82941588188015, 0.201495132985613, 0.298017793415622, 0.0934947735910401, 0.548588115971628, 0.431368751868528}, {}, {}, {}, {}};
	}
private:
	F22generalizedBatchReactor(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F23generalizedBatchReactor : public Expression {
public:
	static std::unique_ptr<F23generalizedBatchReactor> create() {
		Adjacency adj{{23, 24, 25, 26, 27, 28}, {0}, {}};
		AdjacencyDiff adjDiff{{}, {{0, 23}, {0, 24}, {0, 25}, {0, 26}, {0, 27}, {0, 28}}, {}, {}, {}, {}};
		return std::unique_ptr<F23generalizedBatchReactor>(new F23generalizedBatchReactor(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return -2.68813123277342*u[0]*x[23] + 0.188531246624573*u[0]*x[24] + 0.0807810262701528*u[0]*x[25] + 0.00849333263382646*u[0]*x[26] + 0.100111425109776*u[0]*x[27] + 0.661813341920466*u[0]*x[28];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-2.68813123277342*u[0], 0.188531246624573*u[0], 0.0807810262701528*u[0], 0.00849333263382646*u[0], 0.100111425109776*u[0], 0.661813341920466*u[0]}, {-2.68813123277342*x[23] + 0.188531246624573*x[24] + 0.0807810262701528*x[25] + 0.00849333263382646*x[26] + 0.100111425109776*x[27] + 0.661813341920466*x[28]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-2.68813123277342, 0.188531246624573, 0.0807810262701528, 0.00849333263382646, 0.100111425109776, 0.661813341920466}, {}, {}, {}, {}};
	}
private:
	F23generalizedBatchReactor(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F24generalizedBatchReactor : public Expression {
public:
	static std::unique_ptr<F24generalizedBatchReactor> create() {
		Adjacency adj{{24, 25, 26, 27, 28, 29}, {0}, {}};
		AdjacencyDiff adjDiff{{}, {{0, 24}, {0, 25}, {0, 26}, {0, 27}, {0, 28}, {0, 29}}, {}, {}, {}, {}};
		return std::unique_ptr<F24generalizedBatchReactor>(new F24generalizedBatchReactor(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return -2.39224707181371*u[0]*x[24] + 0.892449797147899*u[0]*x[25] + 0.38434299521615*u[0]*x[26] + 0.0892601855464117*u[0]*x[27] + 0.233459605301198*u[0]*x[28] + 0.420124597627528*u[0]*x[29];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-2.39224707181371*u[0], 0.892449797147899*u[0], 0.38434299521615*u[0], 0.0892601855464117*u[0], 0.233459605301198*u[0], 0.420124597627528*u[0]}, {-2.39224707181371*x[24] + 0.892449797147899*x[25] + 0.38434299521615*x[26] + 0.0892601855464117*x[27] + 0.233459605301198*x[28] + 0.420124597627528*x[29]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-2.39224707181371, 0.892449797147899, 0.38434299521615, 0.0892601855464117, 0.233459605301198, 0.420124597627528}, {}, {}, {}, {}};
	}
private:
	F24generalizedBatchReactor(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F25generalizedBatchReactor : public Expression {
public:
	static std::unique_ptr<F25generalizedBatchReactor> create() {
		Adjacency adj{{25, 26, 27, 28, 29, 30}, {0}, {}};
		AdjacencyDiff adjDiff{{}, {{0, 25}, {0, 26}, {0, 27}, {0, 28}, {0, 29}, {0, 30}}, {}, {}, {}, {}};
		return std::unique_ptr<F25generalizedBatchReactor>(new F25generalizedBatchReactor(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return -2.11571801165978*u[0]*x[25] + 0.103729866387657*u[0]*x[26] + 0.547804551921661*u[0]*x[27] + 0.566944967754185*u[0]*x[28] + 0.622266327274425*u[0]*x[29] + 0.651701061053077*u[0]*x[30];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-2.11571801165978*u[0], 0.103729866387657*u[0], 0.547804551921661*u[0], 0.566944967754185*u[0], 0.622266327274425*u[0], 0.651701061053077*u[0]}, {-2.11571801165978*x[25] + 0.103729866387657*x[26] + 0.547804551921661*x[27] + 0.566944967754185*x[28] + 0.622266327274425*x[29] + 0.651701061053077*x[30]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-2.11571801165978, 0.103729866387657, 0.547804551921661, 0.566944967754185, 0.622266327274425, 0.651701061053077}, {}, {}, {}, {}};
	}
private:
	F25generalizedBatchReactor(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F26generalizedBatchReactor : public Expression {
public:
	static std::unique_ptr<F26generalizedBatchReactor> create() {
		Adjacency adj{{26, 27, 28, 29, 30, 31}, {0}, {}};
		AdjacencyDiff adjDiff{{}, {{0, 26}, {0, 27}, {0, 28}, {0, 29}, {0, 30}, {0, 31}}, {}, {}, {}, {}};
		return std::unique_ptr<F26generalizedBatchReactor>(new F26generalizedBatchReactor(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return -2.02116582554809*u[0]*x[26] + 0.0764794130707378*u[0]*x[27] + 0.222963701581816*u[0]*x[28] + 0.786046429891347*u[0]*x[29] + 0.239360900524536*u[0]*x[30] + 0.981064186437453*u[0]*x[31];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-2.02116582554809*u[0], 0.0764794130707378*u[0], 0.222963701581816*u[0], 0.786046429891347*u[0], 0.239360900524536*u[0], 0.981064186437453*u[0]}, {-2.02116582554809*x[26] + 0.0764794130707378*x[27] + 0.222963701581816*x[28] + 0.786046429891347*x[29] + 0.239360900524536*x[30] + 0.981064186437453*x[31]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-2.02116582554809, 0.0764794130707378, 0.222963701581816, 0.786046429891347, 0.239360900524536, 0.981064186437453}, {}, {}, {}, {}};
	}
private:
	F26generalizedBatchReactor(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F27generalizedBatchReactor : public Expression {
public:
	static std::unique_ptr<F27generalizedBatchReactor> create() {
		Adjacency adj{{27, 28, 29, 30, 31, 32}, {0}, {}};
		AdjacencyDiff adjDiff{{}, {{0, 27}, {0, 28}, {0, 29}, {0, 30}, {0, 31}, {0, 32}}, {}, {}, {}, {}};
		return std::unique_ptr<F27generalizedBatchReactor>(new F27generalizedBatchReactor(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return -1.24502432751711*u[0]*x[27] + 0.892431913256511*u[0]*x[28] + 0.873868769578652*u[0]*x[29] + 0.629332181014187*u[0]*x[30] + 0.0486899741579899*u[0]*x[31] + 0.424002644552821*u[0]*x[32];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-1.24502432751711*u[0], 0.892431913256511*u[0], 0.873868769578652*u[0], 0.629332181014187*u[0], 0.0486899741579899*u[0], 0.424002644552821*u[0]}, {-1.24502432751711*x[27] + 0.892431913256511*x[28] + 0.873868769578652*x[29] + 0.629332181014187*x[30] + 0.0486899741579899*x[31] + 0.424002644552821*x[32]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-1.24502432751711, 0.892431913256511, 0.873868769578652, 0.629332181014187, 0.0486899741579899, 0.424002644552821}, {}, {}, {}, {}};
	}
private:
	F27generalizedBatchReactor(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F28generalizedBatchReactor : public Expression {
public:
	static std::unique_ptr<F28generalizedBatchReactor> create() {
		Adjacency adj{{28, 29, 30, 31, 32, 33}, {0}, {}};
		AdjacencyDiff adjDiff{{}, {{0, 28}, {0, 29}, {0, 30}, {0, 31}, {0, 32}, {0, 33}}, {}, {}, {}, {}};
		return std::unique_ptr<F28generalizedBatchReactor>(new F28generalizedBatchReactor(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return -2.57761352981418*u[0]*x[28] + 0.940093955086718*u[0]*x[29] + 0.364250226048989*u[0]*x[30] + 0.749875013324617*u[0]*x[31] + 0.469993686035266*u[0]*x[32] + 0.90822890346605*u[0]*x[33];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-2.57761352981418*u[0], 0.940093955086718*u[0], 0.364250226048989*u[0], 0.749875013324617*u[0], 0.469993686035266*u[0], 0.90822890346605*u[0]}, {-2.57761352981418*x[28] + 0.940093955086718*x[29] + 0.364250226048989*x[30] + 0.749875013324617*x[31] + 0.469993686035266*x[32] + 0.90822890346605*x[33]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-2.57761352981418, 0.940093955086718, 0.364250226048989, 0.749875013324617, 0.469993686035266, 0.90822890346605}, {}, {}, {}, {}};
	}
private:
	F28generalizedBatchReactor(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F29generalizedBatchReactor : public Expression {
public:
	static std::unique_ptr<F29generalizedBatchReactor> create() {
		Adjacency adj{{29, 30, 31, 32, 33, 34}, {0}, {}};
		AdjacencyDiff adjDiff{{}, {{0, 29}, {0, 30}, {0, 31}, {0, 32}, {0, 33}, {0, 34}}, {}, {}, {}, {}};
		return std::unique_ptr<F29generalizedBatchReactor>(new F29generalizedBatchReactor(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return -3.64240007945867*u[0]*x[29] + 0.920734417341879*u[0]*x[30] + 0.10275559707294*u[0]*x[31] + 0.0699162531307997*u[0]*x[32] + 0.804224799514677*u[0]*x[33] + 0.114219317514066*u[0]*x[34];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-3.64240007945867*u[0], 0.920734417341879*u[0], 0.10275559707294*u[0], 0.0699162531307997*u[0], 0.804224799514677*u[0], 0.114219317514066*u[0]}, {-3.64240007945867*x[29] + 0.920734417341879*x[30] + 0.10275559707294*x[31] + 0.0699162531307997*x[32] + 0.804224799514677*x[33] + 0.114219317514066*x[34]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-3.64240007945867, 0.920734417341879, 0.10275559707294, 0.0699162531307997, 0.804224799514677, 0.114219317514066}, {}, {}, {}, {}};
	}
private:
	F29generalizedBatchReactor(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F30generalizedBatchReactor : public Expression {
public:
	static std::unique_ptr<F30generalizedBatchReactor> create() {
		Adjacency adj{{30, 31, 32, 33, 34, 35}, {0}, {}};
		AdjacencyDiff adjDiff{{}, {{0, 30}, {0, 31}, {0, 32}, {0, 33}, {0, 34}, {0, 35}}, {}, {}, {}, {}};
		return std::unique_ptr<F30generalizedBatchReactor>(new F30generalizedBatchReactor(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return -2.80537878598267*u[0]*x[30] + 0.890753072043354*u[0]*x[31] + 0.235086273077543*u[0]*x[32] + 0.482894017486791*u[0]*x[33] + 0.29319610335881*u[0]*x[34] + 0.148162798608384*u[0]*x[35];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-2.80537878598267*u[0], 0.890753072043354*u[0], 0.235086273077543*u[0], 0.482894017486791*u[0], 0.29319610335881*u[0], 0.148162798608384*u[0]}, {-2.80537878598267*x[30] + 0.890753072043354*x[31] + 0.235086273077543*x[32] + 0.482894017486791*x[33] + 0.29319610335881*x[34] + 0.148162798608384*x[35]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-2.80537878598267, 0.890753072043354, 0.235086273077543, 0.482894017486791, 0.29319610335881, 0.148162798608384}, {}, {}, {}, {}};
	}
private:
	F30generalizedBatchReactor(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F31generalizedBatchReactor : public Expression {
public:
	static std::unique_ptr<F31generalizedBatchReactor> create() {
		Adjacency adj{{31, 32, 33, 34, 35, 36}, {0}, {}};
		AdjacencyDiff adjDiff{{}, {{0, 31}, {0, 32}, {0, 33}, {0, 34}, {0, 35}, {0, 36}}, {}, {}, {}, {}};
		return std::unique_ptr<F31generalizedBatchReactor>(new F31generalizedBatchReactor(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return -2.77313784303635*u[0]*x[31] + 0.504674801013974*u[0]*x[32] + 0.742144936686235*u[0]*x[33] + 0.298329045335977*u[0]*x[34] + 0.301839638913647*u[0]*x[35] + 0.115608418140323*u[0]*x[36];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-2.77313784303635*u[0], 0.504674801013974*u[0], 0.742144936686235*u[0], 0.298329045335977*u[0], 0.301839638913647*u[0], 0.115608418140323*u[0]}, {-2.77313784303635*x[31] + 0.504674801013974*x[32] + 0.742144936686235*x[33] + 0.298329045335977*x[34] + 0.301839638913647*x[35] + 0.115608418140323*x[36]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-2.77313784303635, 0.504674801013974, 0.742144936686235, 0.298329045335977, 0.301839638913647, 0.115608418140323}, {}, {}, {}, {}};
	}
private:
	F31generalizedBatchReactor(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F32generalizedBatchReactor : public Expression {
public:
	static std::unique_ptr<F32generalizedBatchReactor> create() {
		Adjacency adj{{32, 33, 34, 35, 36, 37}, {0}, {}};
		AdjacencyDiff adjDiff{{}, {{0, 32}, {0, 33}, {0, 34}, {0, 35}, {0, 36}, {0, 37}}, {}, {}, {}, {}};
		return std::unique_ptr<F32generalizedBatchReactor>(new F32generalizedBatchReactor(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return -1.7036736578104*u[0]*x[32] + 0.537686749957713*u[0]*x[33] + 0.447158297883282*u[0]*x[34] + 0.865091678747471*u[0]*x[35] + 0.275336609405286*u[0]*x[36] + 0.374285053608546*u[0]*x[37];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-1.7036736578104*u[0], 0.537686749957713*u[0], 0.447158297883282*u[0], 0.865091678747471*u[0], 0.275336609405286*u[0], 0.374285053608546*u[0]}, {-1.7036736578104*x[32] + 0.537686749957713*x[33] + 0.447158297883282*x[34] + 0.865091678747471*x[35] + 0.275336609405286*x[36] + 0.374285053608546*x[37]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-1.7036736578104, 0.537686749957713, 0.447158297883282, 0.865091678747471, 0.275336609405286, 0.374285053608546}, {}, {}, {}, {}};
	}
private:
	F32generalizedBatchReactor(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F33generalizedBatchReactor : public Expression {
public:
	static std::unique_ptr<F33generalizedBatchReactor> create() {
		Adjacency adj{{33, 34, 35, 36, 37, 38}, {0}, {}};
		AdjacencyDiff adjDiff{{}, {{0, 33}, {0, 34}, {0, 35}, {0, 36}, {0, 37}, {0, 38}}, {}, {}, {}, {}};
		return std::unique_ptr<F33generalizedBatchReactor>(new F33generalizedBatchReactor(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return -3.47517940711147*u[0]*x[33] + 0.432510890088799*u[0]*x[34] + 0.206869416752813*u[0]*x[35] + 0.405461858431165*u[0]*x[36] + 0.536970907675484*u[0]*x[37] + 0.942826248013742*u[0]*x[38];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-3.47517940711147*u[0], 0.432510890088799*u[0], 0.206869416752813*u[0], 0.405461858431165*u[0], 0.536970907675484*u[0], 0.942826248013742*u[0]}, {-3.47517940711147*x[33] + 0.432510890088799*x[34] + 0.206869416752813*x[35] + 0.405461858431165*x[36] + 0.536970907675484*x[37] + 0.942826248013742*x[38]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-3.47517940711147, 0.432510890088799, 0.206869416752813, 0.405461858431165, 0.536970907675484, 0.942826248013742}, {}, {}, {}, {}};
	}
private:
	F33generalizedBatchReactor(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F34generalizedBatchReactor : public Expression {
public:
	static std::unique_ptr<F34generalizedBatchReactor> create() {
		Adjacency adj{{34, 35, 36, 37, 38, 39}, {0}, {}};
		AdjacencyDiff adjDiff{{}, {{0, 34}, {0, 35}, {0, 36}, {0, 37}, {0, 38}, {0, 39}}, {}, {}, {}, {}};
		return std::unique_ptr<F34generalizedBatchReactor>(new F34generalizedBatchReactor(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return -1.58541365418093*u[0]*x[34] + 0.49427601054591*u[0]*x[35] + 0.148661140204073*u[0]*x[36] + 0.971711827161295*u[0]*x[37] + 0.93660533290218*u[0]*x[38] + 0.258743100147279*u[0]*x[39];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-1.58541365418093*u[0], 0.49427601054591*u[0], 0.148661140204073*u[0], 0.971711827161295*u[0], 0.93660533290218*u[0], 0.258743100147279*u[0]}, {-1.58541365418093*x[34] + 0.49427601054591*x[35] + 0.148661140204073*x[36] + 0.971711827161295*x[37] + 0.93660533290218*x[38] + 0.258743100147279*x[39]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-1.58541365418093, 0.49427601054591, 0.148661140204073, 0.971711827161295, 0.93660533290218, 0.258743100147279}, {}, {}, {}, {}};
	}
private:
	F34generalizedBatchReactor(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F35generalizedBatchReactor : public Expression {
public:
	static std::unique_ptr<F35generalizedBatchReactor> create() {
		Adjacency adj{{35, 36, 37, 38, 39, 40}, {0}, {}};
		AdjacencyDiff adjDiff{{}, {{0, 35}, {0, 36}, {0, 37}, {0, 38}, {0, 39}, {0, 40}}, {}, {}, {}, {}};
		return std::unique_ptr<F35generalizedBatchReactor>(new F35generalizedBatchReactor(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return -2.01623954356822*u[0]*x[35] + 0.92136604189966*u[0]*x[36] + 0.638444356678476*u[0]*x[37] + 0.918712462954774*u[0]*x[38] + 0.606502101183564*u[0]*x[39] + 0.90310424866317*u[0]*x[40];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-2.01623954356822*u[0], 0.92136604189966*u[0], 0.638444356678476*u[0], 0.918712462954774*u[0], 0.606502101183564*u[0], 0.90310424866317*u[0]}, {-2.01623954356822*x[35] + 0.92136604189966*x[36] + 0.638444356678476*x[37] + 0.918712462954774*x[38] + 0.606502101183564*x[39] + 0.90310424866317*x[40]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-2.01623954356822, 0.92136604189966, 0.638444356678476, 0.918712462954774, 0.606502101183564, 0.90310424866317}, {}, {}, {}, {}};
	}
private:
	F35generalizedBatchReactor(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F36generalizedBatchReactor : public Expression {
public:
	static std::unique_ptr<F36generalizedBatchReactor> create() {
		Adjacency adj{{36, 37, 38, 39, 40, 41}, {0}, {}};
		AdjacencyDiff adjDiff{{}, {{0, 36}, {0, 37}, {0, 38}, {0, 39}, {0, 40}, {0, 41}}, {}, {}, {}, {}};
		return std::unique_ptr<F36generalizedBatchReactor>(new F36generalizedBatchReactor(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return -1.86643406808051*u[0]*x[36] + 0.904612204169701*u[0]*x[37] + 0.431568773980985*u[0]*x[38] + 0.149503560060005*u[0]*x[39] + 0.8882603465276*u[0]*x[40] + 0.577967350985635*u[0]*x[41];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-1.86643406808051*u[0], 0.904612204169701*u[0], 0.431568773980985*u[0], 0.149503560060005*u[0], 0.8882603465276*u[0], 0.577967350985635*u[0]}, {-1.86643406808051*x[36] + 0.904612204169701*x[37] + 0.431568773980985*x[38] + 0.149503560060005*x[39] + 0.8882603465276*x[40] + 0.577967350985635*x[41]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-1.86643406808051, 0.904612204169701, 0.431568773980985, 0.149503560060005, 0.8882603465276, 0.577967350985635}, {}, {}, {}, {}};
	}
private:
	F36generalizedBatchReactor(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F37generalizedBatchReactor : public Expression {
public:
	static std::unique_ptr<F37generalizedBatchReactor> create() {
		Adjacency adj{{37, 38, 39, 40, 41, 42}, {0}, {}};
		AdjacencyDiff adjDiff{{}, {{0, 37}, {0, 38}, {0, 39}, {0, 40}, {0, 41}, {0, 42}}, {}, {}, {}, {}};
		return std::unique_ptr<F37generalizedBatchReactor>(new F37generalizedBatchReactor(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return -3.4260243492935*u[0]*x[37] + 0.196398871711529*u[0]*x[38] + 0.551083739146974*u[0]*x[39] + 0.389321156033436*u[0]*x[40] + 0.32536646427409*u[0]*x[41] + 0.998998402425469*u[0]*x[42];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-3.4260243492935*u[0], 0.196398871711529*u[0], 0.551083739146974*u[0], 0.389321156033436*u[0], 0.32536646427409*u[0], 0.998998402425469*u[0]}, {-3.4260243492935*x[37] + 0.196398871711529*x[38] + 0.551083739146974*x[39] + 0.389321156033436*x[40] + 0.32536646427409*x[41] + 0.998998402425469*x[42]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-3.4260243492935, 0.196398871711529, 0.551083739146974, 0.389321156033436, 0.32536646427409, 0.998998402425469}, {}, {}, {}, {}};
	}
private:
	F37generalizedBatchReactor(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F38generalizedBatchReactor : public Expression {
public:
	static std::unique_ptr<F38generalizedBatchReactor> create() {
		Adjacency adj{{38, 39, 40, 41, 42, 43}, {0}, {}};
		AdjacencyDiff adjDiff{{}, {{0, 38}, {0, 39}, {0, 40}, {0, 41}, {0, 42}, {0, 43}}, {}, {}, {}, {}};
		return std::unique_ptr<F38generalizedBatchReactor>(new F38generalizedBatchReactor(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return -3.42611168956321*u[0]*x[38] + 0.13657956997807*u[0]*x[39] + 0.882726827639407*u[0]*x[40] + 0.763882813740335*u[0]*x[41] + 0.04157112688015*u[0]*x[42] + 0.442937951144497*u[0]*x[43];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-3.42611168956321*u[0], 0.13657956997807*u[0], 0.882726827639407*u[0], 0.763882813740335*u[0], 0.04157112688015*u[0], 0.442937951144497*u[0]}, {-3.42611168956321*x[38] + 0.13657956997807*x[39] + 0.882726827639407*x[40] + 0.763882813740335*x[41] + 0.04157112688015*x[42] + 0.442937951144497*x[43]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-3.42611168956321, 0.13657956997807, 0.882726827639407, 0.763882813740335, 0.04157112688015, 0.442937951144497}, {}, {}, {}, {}};
	}
private:
	F38generalizedBatchReactor(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F39generalizedBatchReactor : public Expression {
public:
	static std::unique_ptr<F39generalizedBatchReactor> create() {
		Adjacency adj{{39, 40, 41, 42, 43, 44}, {0}, {}};
		AdjacencyDiff adjDiff{{}, {{0, 39}, {0, 40}, {0, 41}, {0, 42}, {0, 43}, {0, 44}}, {}, {}, {}, {}};
		return std::unique_ptr<F39generalizedBatchReactor>(new F39generalizedBatchReactor(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return -1.70241207051589*u[0]*x[39] + 0.859340822608575*u[0]*x[40] + 0.988494243628834*u[0]*x[41] + 0.926947266663991*u[0]*x[42] + 0.966506951028157*u[0]*x[43] + 0.175177508291523*u[0]*x[44];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-1.70241207051589*u[0], 0.859340822608575*u[0], 0.988494243628834*u[0], 0.926947266663991*u[0], 0.966506951028157*u[0], 0.175177508291523*u[0]}, {-1.70241207051589*x[39] + 0.859340822608575*x[40] + 0.988494243628834*x[41] + 0.926947266663991*x[42] + 0.966506951028157*x[43] + 0.175177508291523*x[44]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-1.70241207051589, 0.859340822608575, 0.988494243628834, 0.926947266663991, 0.966506951028157, 0.175177508291523}, {}, {}, {}, {}};
	}
private:
	F39generalizedBatchReactor(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F40generalizedBatchReactor : public Expression {
public:
	static std::unique_ptr<F40generalizedBatchReactor> create() {
		Adjacency adj{{40, 41, 42, 43, 44, 45}, {0}, {}};
		AdjacencyDiff adjDiff{{}, {{0, 40}, {0, 41}, {0, 42}, {0, 43}, {0, 44}, {0, 45}}, {}, {}, {}, {}};
		return std::unique_ptr<F40generalizedBatchReactor>(new F40generalizedBatchReactor(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return -3.92275340147219*u[0]*x[40] + 0.213841649340535*u[0]*x[41] + 0.149526378225927*u[0]*x[42] + 0.16797682459807*u[0]*x[43] + 0.895334603930632*u[0]*x[44] + 0.0890122517205226*u[0]*x[45];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-3.92275340147219*u[0], 0.213841649340535*u[0], 0.149526378225927*u[0], 0.16797682459807*u[0], 0.895334603930632*u[0], 0.0890122517205226*u[0]}, {-3.92275340147219*x[40] + 0.213841649340535*x[41] + 0.149526378225927*x[42] + 0.16797682459807*x[43] + 0.895334603930632*x[44] + 0.0890122517205226*x[45]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-3.92275340147219, 0.213841649340535, 0.149526378225927, 0.16797682459807, 0.895334603930632, 0.0890122517205226}, {}, {}, {}, {}};
	}
private:
	F40generalizedBatchReactor(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F41generalizedBatchReactor : public Expression {
public:
	static std::unique_ptr<F41generalizedBatchReactor> create() {
		Adjacency adj{{41, 42, 43, 44, 45, 46}, {0}, {}};
		AdjacencyDiff adjDiff{{}, {{0, 41}, {0, 42}, {0, 43}, {0, 44}, {0, 45}, {0, 46}}, {}, {}, {}, {}};
		return std::unique_ptr<F41generalizedBatchReactor>(new F41generalizedBatchReactor(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return -2.86955252196943*u[0]*x[41] + 0.450482622507049*u[0]*x[42] + 0.207591873499632*u[0]*x[43] + 0.0726866017564705*u[0]*x[44] + 0.831571264172694*u[0]*x[45] + 0.0515347793257441*u[0]*x[46];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-2.86955252196943*u[0], 0.450482622507049*u[0], 0.207591873499632*u[0], 0.0726866017564705*u[0], 0.831571264172694*u[0], 0.0515347793257441*u[0]}, {-2.86955252196943*x[41] + 0.450482622507049*x[42] + 0.207591873499632*x[43] + 0.0726866017564705*x[44] + 0.831571264172694*x[45] + 0.0515347793257441*x[46]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-2.86955252196943, 0.450482622507049, 0.207591873499632, 0.0726866017564705, 0.831571264172694, 0.0515347793257441}, {}, {}, {}, {}};
	}
private:
	F41generalizedBatchReactor(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F42generalizedBatchReactor : public Expression {
public:
	static std::unique_ptr<F42generalizedBatchReactor> create() {
		Adjacency adj{{42, 43, 44, 45, 46, 47}, {0}, {}};
		AdjacencyDiff adjDiff{{}, {{0, 42}, {0, 43}, {0, 44}, {0, 45}, {0, 46}, {0, 47}}, {}, {}, {}, {}};
		return std::unique_ptr<F42generalizedBatchReactor>(new F42generalizedBatchReactor(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return -2.56752579670259*u[0]*x[42] + 0.339647362068514*u[0]*x[43] + 0.459866300843307*u[0]*x[44] + 0.204932215198737*u[0]*x[45] + 0.950915589449735*u[0]*x[46] + 0.420572744283715*u[0]*x[47];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-2.56752579670259*u[0], 0.339647362068514*u[0], 0.459866300843307*u[0], 0.204932215198737*u[0], 0.950915589449735*u[0], 0.420572744283715*u[0]}, {-2.56752579670259*x[42] + 0.339647362068514*x[43] + 0.459866300843307*x[44] + 0.204932215198737*x[45] + 0.950915589449735*x[46] + 0.420572744283715*x[47]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-2.56752579670259, 0.339647362068514, 0.459866300843307, 0.204932215198737, 0.950915589449735, 0.420572744283715}, {}, {}, {}, {}};
	}
private:
	F42generalizedBatchReactor(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F43generalizedBatchReactor : public Expression {
public:
	static std::unique_ptr<F43generalizedBatchReactor> create() {
		Adjacency adj{{43, 44, 45, 46, 47, 48}, {0}, {}};
		AdjacencyDiff adjDiff{{}, {{0, 43}, {0, 44}, {0, 45}, {0, 46}, {0, 47}, {0, 48}}, {}, {}, {}, {}};
		return std::unique_ptr<F43generalizedBatchReactor>(new F43generalizedBatchReactor(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return -2.12466096233887*u[0]*x[43] + 0.111359011152279*u[0]*x[44] + 0.0559432230104395*u[0]*x[45] + 0.0671790552687289*u[0]*x[46] + 0.0542056192017173*u[0]*x[47] + 0.991543739634836*u[0]*x[48];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-2.12466096233887*u[0], 0.111359011152279*u[0], 0.0559432230104395*u[0], 0.0671790552687289*u[0], 0.0542056192017173*u[0], 0.991543739634836*u[0]}, {-2.12466096233887*x[43] + 0.111359011152279*x[44] + 0.0559432230104395*x[45] + 0.0671790552687289*x[46] + 0.0542056192017173*x[47] + 0.991543739634836*x[48]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-2.12466096233887, 0.111359011152279, 0.0559432230104395, 0.0671790552687289, 0.0542056192017173, 0.991543739634836}, {}, {}, {}, {}};
	}
private:
	F43generalizedBatchReactor(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F44generalizedBatchReactor : public Expression {
public:
	static std::unique_ptr<F44generalizedBatchReactor> create() {
		Adjacency adj{{44, 45, 46, 47, 48, 49}, {0}, {}};
		AdjacencyDiff adjDiff{{}, {{0, 44}, {0, 45}, {0, 46}, {0, 47}, {0, 48}, {0, 49}}, {}, {}, {}, {}};
		return std::unique_ptr<F44generalizedBatchReactor>(new F44generalizedBatchReactor(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return -1.71442402597421*u[0]*x[44] + 0.842502159577236*u[0]*x[45] + 0.45520023747691*u[0]*x[46] + 0.370036645505251*u[0]*x[47] + 0.789057795242135*u[0]*x[48] + 0.16149963701895*u[0]*x[49];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-1.71442402597421*u[0], 0.842502159577236*u[0], 0.45520023747691*u[0], 0.370036645505251*u[0], 0.789057795242135*u[0], 0.16149963701895*u[0]}, {-1.71442402597421*x[44] + 0.842502159577236*x[45] + 0.45520023747691*x[46] + 0.370036645505251*x[47] + 0.789057795242135*x[48] + 0.16149963701895*x[49]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-1.71442402597421, 0.842502159577236, 0.45520023747691, 0.370036645505251, 0.789057795242135, 0.16149963701895}, {}, {}, {}, {}};
	}
private:
	F44generalizedBatchReactor(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F45generalizedBatchReactor : public Expression {
public:
	static std::unique_ptr<F45generalizedBatchReactor> create() {
		Adjacency adj{{45, 46, 47, 48, 49}, {0}, {}};
		AdjacencyDiff adjDiff{{}, {{0, 45}, {0, 46}, {0, 47}, {0, 48}, {0, 49}}, {}, {}, {}, {}};
		return std::unique_ptr<F45generalizedBatchReactor>(new F45generalizedBatchReactor(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return -2.02396111367963*u[0]*x[45] + 0.862393814179248*u[0]*x[46] + 0.813714602929293*u[0]*x[47] + 0.504026444563673*u[0]*x[48] + 0.49150149918123*u[0]*x[49];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-2.02396111367963*u[0], 0.862393814179248*u[0], 0.813714602929293*u[0], 0.504026444563673*u[0], 0.49150149918123*u[0]}, {-2.02396111367963*x[45] + 0.862393814179248*x[46] + 0.813714602929293*x[47] + 0.504026444563673*x[48] + 0.49150149918123*x[49]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-2.02396111367963, 0.862393814179248, 0.813714602929293, 0.504026444563673, 0.49150149918123}, {}, {}, {}, {}};
	}
private:
	F45generalizedBatchReactor(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F46generalizedBatchReactor : public Expression {
public:
	static std::unique_ptr<F46generalizedBatchReactor> create() {
		Adjacency adj{{46, 47, 48, 49}, {0}, {}};
		AdjacencyDiff adjDiff{{}, {{0, 46}, {0, 47}, {0, 48}, {0, 49}}, {}, {}, {}, {}};
		return std::unique_ptr<F46generalizedBatchReactor>(new F46generalizedBatchReactor(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return -2.38722347570037*u[0]*x[46] + 0.775818307290483*u[0]*x[47] + 0.0563080592826334*u[0]*x[48] + 0.298844388313121*u[0]*x[49];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-2.38722347570037*u[0], 0.775818307290483*u[0], 0.0563080592826334*u[0], 0.298844388313121*u[0]}, {-2.38722347570037*x[46] + 0.775818307290483*x[47] + 0.0563080592826334*x[48] + 0.298844388313121*x[49]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-2.38722347570037, 0.775818307290483, 0.0563080592826334, 0.298844388313121}, {}, {}, {}, {}};
	}
private:
	F46generalizedBatchReactor(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F47generalizedBatchReactor : public Expression {
public:
	static std::unique_ptr<F47generalizedBatchReactor> create() {
		Adjacency adj{{47, 48, 49}, {0}, {}};
		AdjacencyDiff adjDiff{{}, {{0, 47}, {0, 48}, {0, 49}}, {}, {}, {}, {}};
		return std::unique_ptr<F47generalizedBatchReactor>(new F47generalizedBatchReactor(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return -2.43434791921046*u[0]*x[47] + 0.475398609629483*u[0]*x[48] + 0.21850345064235*u[0]*x[49];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-2.43434791921046*u[0], 0.475398609629483*u[0], 0.21850345064235*u[0]}, {-2.43434791921046*x[47] + 0.475398609629483*x[48] + 0.21850345064235*x[49]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-2.43434791921046, 0.475398609629483, 0.21850345064235}, {}, {}, {}, {}};
	}
private:
	F47generalizedBatchReactor(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F48generalizedBatchReactor : public Expression {
public:
	static std::unique_ptr<F48generalizedBatchReactor> create() {
		Adjacency adj{{48, 49}, {0}, {}};
		AdjacencyDiff adjDiff{{}, {{0, 48}, {0, 49}}, {}, {}, {}, {}};
		return std::unique_ptr<F48generalizedBatchReactor>(new F48generalizedBatchReactor(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return -2.81633464835276*u[0]*x[48] + 0.482818063200095*u[0]*x[49];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-2.81633464835276*u[0], 0.482818063200095*u[0]}, {-2.81633464835276*x[48] + 0.482818063200095*x[49]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-2.81633464835276, 0.482818063200095}, {}, {}, {}, {}};
	}
private:
	F48generalizedBatchReactor(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F49generalizedBatchReactor : public Expression {
public:
	static std::unique_ptr<F49generalizedBatchReactor> create() {
		Adjacency adj{{0, 49}, {0}, {}};
		AdjacencyDiff adjDiff{{}, {{0, 0}, {0, 49}}, {{0, 0}}, {}, {}, {}};
		return std::unique_ptr<F49generalizedBatchReactor>(new F49generalizedBatchReactor(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return -1.65316703835575*u[0]*x[49] + pow(u[0], 2)*x[0]*DEPLETION_COEFF_VALUE;
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{pow(u[0], 2)*DEPLETION_COEFF_VALUE, -1.65316703835575*u[0]}, {-1.65316703835575*x[49] + 2*u[0]*x[0]*DEPLETION_COEFF_VALUE}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        const double x0 = 2*DEPLETION_COEFF_VALUE;
        const double x1 = x0*u[0];
		return {std::vector<double>{}, {x1, -1.65316703835575}, {x0*x[0]}, {}, {}, {}};
	}
private:
	F49generalizedBatchReactor(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F50generalizedBatchReactor : public Expression {
public:
	static std::unique_ptr<F50generalizedBatchReactor> create() {
		Adjacency adj{{}, {0}, {}};
		AdjacencyDiff adjDiff{{}, {}, {{0, 0}}, {}, {}, {}};
		return std::unique_ptr<F50generalizedBatchReactor>(new F50generalizedBatchReactor(std::move(adj), std::move(adjDiff)));
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
	F50generalizedBatchReactor(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


// final constraints
class R0generalizedBatchReactor : public Constraint {
public:
	static std::unique_ptr<R0generalizedBatchReactor> create() {
		Adjacency adj{{50}, {}, {}};
		AdjacencyDiff adjDiff{{}, {}, {}, {}, {}, {}};
		return std::unique_ptr<R0generalizedBatchReactor>(new R0generalizedBatchReactor(std::move(adj), std::move(adjDiff), MINUS_INFINITY, 0.5));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return x[50];
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
    F.push_back(F16generalizedBatchReactor::create());
    F.push_back(F17generalizedBatchReactor::create());
    F.push_back(F18generalizedBatchReactor::create());
    F.push_back(F19generalizedBatchReactor::create());
    F.push_back(F20generalizedBatchReactor::create());
    F.push_back(F21generalizedBatchReactor::create());
    F.push_back(F22generalizedBatchReactor::create());
    F.push_back(F23generalizedBatchReactor::create());
    F.push_back(F24generalizedBatchReactor::create());
    F.push_back(F25generalizedBatchReactor::create());
    F.push_back(F26generalizedBatchReactor::create());
    F.push_back(F27generalizedBatchReactor::create());
    F.push_back(F28generalizedBatchReactor::create());
    F.push_back(F29generalizedBatchReactor::create());
    F.push_back(F30generalizedBatchReactor::create());
    F.push_back(F31generalizedBatchReactor::create());
    F.push_back(F32generalizedBatchReactor::create());
    F.push_back(F33generalizedBatchReactor::create());
    F.push_back(F34generalizedBatchReactor::create());
    F.push_back(F35generalizedBatchReactor::create());
    F.push_back(F36generalizedBatchReactor::create());
    F.push_back(F37generalizedBatchReactor::create());
    F.push_back(F38generalizedBatchReactor::create());
    F.push_back(F39generalizedBatchReactor::create());
    F.push_back(F40generalizedBatchReactor::create());
    F.push_back(F41generalizedBatchReactor::create());
    F.push_back(F42generalizedBatchReactor::create());
    F.push_back(F43generalizedBatchReactor::create());
    F.push_back(F44generalizedBatchReactor::create());
    F.push_back(F45generalizedBatchReactor::create());
    F.push_back(F46generalizedBatchReactor::create());
    F.push_back(F47generalizedBatchReactor::create());
    F.push_back(F48generalizedBatchReactor::create());
    F.push_back(F49generalizedBatchReactor::create());
    F.push_back(F50generalizedBatchReactor::create());
    
    std::vector<std::unique_ptr<Constraint>> G;
    
    
    std::vector<std::unique_ptr<Constraint>> R;
    R.push_back(R0generalizedBatchReactor::create());
    
    std::vector<std::unique_ptr<ParamConstraint>> A;
    

    Problem problem(
            51, 1, 0,  // #vars
            {0, 0.0204081632653061, 0.0204081632653061, 0.0204081632653061, 0.0204081632653061, 0.0204081632653061, 0.0204081632653061, 0.0204081632653061, 0.0204081632653061, 0.0204081632653061, 0.0204081632653061, 0.0204081632653061, 0.0204081632653061, 0.0204081632653061, 0.0204081632653061, 0.0204081632653061, 0.0204081632653061, 0.0204081632653061, 0.0204081632653061, 0.0204081632653061, 0.0204081632653061, 0.0204081632653061, 0.0204081632653061, 0.0204081632653061, 0.0204081632653061, 0.0204081632653061, 0.0204081632653061, 0.0204081632653061, 0.0204081632653061, 0.0204081632653061, 0.0204081632653061, 0.0204081632653061, 0.0204081632653061, 0.0204081632653061, 0.0204081632653061, 0.0204081632653061, 0.0204081632653061, 0.0204081632653061, 0.0204081632653061, 0.0204081632653061, 0.0204081632653061, 0.0204081632653061, 0.0204081632653061, 0.0204081632653061, 0.0204081632653061, 0.0204081632653061, 0.0204081632653061, 0.0204081632653061, 0.0204081632653061, 0.0204081632653061, 0},  // x0
            {MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY},  // lb x
            {PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY},  // ub x
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
            
    #ifdef INITIAL_STATES_PATH
    problem.initialStatesPath = INITIAL_STATES_PATH "/initialValues.csv";
    #endif
    
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
        

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
		return 0.262837954846131*u[0]*x[1] + 0.112364109510204*u[0]*x[2] + 0.242955829895105*u[0]*x[3] + 0.622942338739852*u[0]*x[4] + 0.716087186158134*u[0]*x[5] - pow(u[0], 2)*x[0]*DEPLETION_COEFF_VALUE;
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-pow(u[0], 2)*DEPLETION_COEFF_VALUE, 0.262837954846131*u[0], 0.112364109510204*u[0], 0.242955829895105*u[0], 0.622942338739852*u[0], 0.716087186158134*u[0]}, {0.262837954846131*x[1] + 0.112364109510204*x[2] + 0.242955829895105*x[3] + 0.622942338739852*x[4] + 0.716087186158134*x[5] - 2*u[0]*x[0]*DEPLETION_COEFF_VALUE}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        const double x0 = 2*DEPLETION_COEFF_VALUE;
        const double x1 = -x0*u[0];
		return {std::vector<double>{}, {x1, 0.262837954846131, 0.112364109510204, 0.242955829895105, 0.622942338739852, 0.716087186158134}, {-x0*x[0]}, {}, {}, {}};
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
		return -0.262837954846131*u[0]*x[1] + 0.868971411358313*u[0]*x[2] + 0.661014601059129*u[0]*x[3] + 0.187232302195506*u[0]*x[4] + 0.6150528991228*u[0]*x[5] + 0.265313650155053*u[0]*x[6];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-0.262837954846131*u[0], 0.868971411358313*u[0], 0.661014601059129*u[0], 0.187232302195506*u[0], 0.6150528991228*u[0], 0.265313650155053*u[0]}, {-0.262837954846131*x[1] + 0.868971411358313*x[2] + 0.661014601059129*x[3] + 0.187232302195506*x[4] + 0.6150528991228*x[5] + 0.265313650155053*x[6]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-0.262837954846131, 0.868971411358313, 0.661014601059129, 0.187232302195506, 0.6150528991228, 0.265313650155053}, {}, {}, {}, {}};
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
		return -0.981335520868517*u[0]*x[2] + 0.469443088818075*u[0]*x[3] + 0.137381947162097*u[0]*x[4] + 0.998902564119753*u[0]*x[5] + 0.549146732319529*u[0]*x[6] + 0.58921650148467*u[0]*x[7];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-0.981335520868517*u[0], 0.469443088818075*u[0], 0.137381947162097*u[0], 0.998902564119753*u[0], 0.549146732319529*u[0], 0.58921650148467*u[0]}, {-0.981335520868517*x[2] + 0.469443088818075*x[3] + 0.137381947162097*x[4] + 0.998902564119753*x[5] + 0.549146732319529*x[6] + 0.58921650148467*x[7]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-0.981335520868517, 0.469443088818075, 0.137381947162097, 0.998902564119753, 0.549146732319529, 0.58921650148467}, {}, {}, {}, {}};
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
		return -1.37341351977231*u[0]*x[3] + 0.808722970205673*u[0]*x[4] + 0.505711638831961*u[0]*x[5] + 0.599024378269749*u[0]*x[6] + 0.45089609562194*u[0]*x[7] + 0.390304717922356*u[0]*x[8];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-1.37341351977231*u[0], 0.808722970205673*u[0], 0.505711638831961*u[0], 0.599024378269749*u[0], 0.45089609562194*u[0], 0.390304717922356*u[0]}, {-1.37341351977231*x[3] + 0.808722970205673*x[4] + 0.505711638831961*x[5] + 0.599024378269749*x[6] + 0.45089609562194*x[7] + 0.390304717922356*x[8]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-1.37341351977231, 0.808722970205673, 0.505711638831961, 0.599024378269749, 0.45089609562194, 0.390304717922356}, {}, {}, {}, {}};
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
		return -1.75627955830313*u[0]*x[4] + 0.209442433556091*u[0]*x[5] + 0.124000153334855*u[0]*x[6] + 0.213579687145035*u[0]*x[7] + 0.220736955986666*u[0]*x[8] + 0.426421214426915*u[0]*x[9];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-1.75627955830313*u[0], 0.209442433556091*u[0], 0.124000153334855*u[0], 0.213579687145035*u[0], 0.220736955986666*u[0], 0.426421214426915*u[0]}, {-1.75627955830313*x[4] + 0.209442433556091*x[5] + 0.124000153334855*x[6] + 0.213579687145035*x[7] + 0.220736955986666*x[8] + 0.426421214426915*x[9]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-1.75627955830313, 0.209442433556091, 0.124000153334855, 0.213579687145035, 0.220736955986666, 0.426421214426915}, {}, {}, {}, {}};
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
		return 0.572015914528953*u[0]*x[10] - 3.04519672178874*u[0]*x[5] + 0.411434459459512*u[0]*x[6] + 0.996756294627156*u[0]*x[7] + 0.279875253855369*u[0]*x[8] + 0.450847649395036*u[0]*x[9];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-3.04519672178874*u[0], 0.411434459459512*u[0], 0.996756294627156*u[0], 0.279875253855369*u[0], 0.450847649395036*u[0], 0.572015914528953*u[0]}, {0.572015914528953*x[10] - 3.04519672178874*x[5] + 0.411434459459512*x[6] + 0.996756294627156*x[7] + 0.279875253855369*x[8] + 0.450847649395036*x[9]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-3.04519672178874, 0.411434459459512, 0.996756294627156, 0.279875253855369, 0.450847649395036, 0.572015914528953}, {}, {}, {}, {}};
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
		return 0.522911070516841*u[0]*x[10] + 0.797942071820214*u[0]*x[11] - 1.9489193735387*u[0]*x[6] + 0.877452948308734*u[0]*x[7] + 0.0178836775722142*u[0]*x[8] + 0.439856216898264*u[0]*x[9];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-1.9489193735387*u[0], 0.877452948308734*u[0], 0.0178836775722142*u[0], 0.439856216898264*u[0], 0.522911070516841*u[0], 0.797942071820214*u[0]}, {0.522911070516841*x[10] + 0.797942071820214*x[11] - 1.9489193735387*x[6] + 0.877452948308734*x[7] + 0.0178836775722142*x[8] + 0.439856216898264*x[9]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-1.9489193735387, 0.877452948308734, 0.0178836775722142, 0.439856216898264, 0.522911070516841, 0.797942071820214}, {}, {}, {}, {}};
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
		return 0.896006224442746*u[0]*x[10] + 0.469635148025699*u[0]*x[11] + 0.0298159579527932*u[0]*x[12] - 3.12790152718754*u[0]*x[7] + 0.467425323365412*u[0]*x[8] + 0.783396135185322*u[0]*x[9];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-3.12790152718754*u[0], 0.467425323365412*u[0], 0.783396135185322*u[0], 0.896006224442746*u[0], 0.469635148025699*u[0], 0.0298159579527932*u[0]}, {0.896006224442746*x[10] + 0.469635148025699*x[11] + 0.0298159579527932*x[12] - 3.12790152718754*x[7] + 0.467425323365412*x[8] + 0.783396135185322*x[9]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-3.12790152718754, 0.467425323365412, 0.783396135185322, 0.896006224442746, 0.469635148025699, 0.0298159579527932}, {}, {}, {}, {}};
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
		return 0.974956747013097*u[0]*x[10] + 0.989015045004831*u[0]*x[11] + 0.878427784085761*u[0]*x[12] + 0.968292871725957*u[0]*x[13] - 1.37622592870202*u[0]*x[8] + 0.977992879098449*u[0]*x[9];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-1.37622592870202*u[0], 0.977992879098449*u[0], 0.974956747013097*u[0], 0.989015045004831*u[0], 0.878427784085761*u[0], 0.968292871725957*u[0]}, {0.974956747013097*x[10] + 0.989015045004831*x[11] + 0.878427784085761*x[12] + 0.968292871725957*x[13] - 1.37622592870202*x[8] + 0.977992879098449*x[9]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-1.37622592870202, 0.977992879098449, 0.974956747013097, 0.989015045004831, 0.878427784085761, 0.968292871725957}, {}, {}, {}, {}};
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
		return 0.920590812717961*u[0]*x[10] + 0.00269678297809628*u[0]*x[11] + 0.0386749962186403*u[0]*x[12] + 0.703439767346247*u[0]*x[13] + 0.416017559854925*u[0]*x[14] - 3.07851409500399*u[0]*x[9];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-3.07851409500399*u[0], 0.920590812717961*u[0], 0.00269678297809628*u[0], 0.0386749962186403*u[0], 0.703439767346247*u[0], 0.416017559854925*u[0]}, {0.920590812717961*x[10] + 0.00269678297809628*x[11] + 0.0386749962186403*x[12] + 0.703439767346247*x[13] + 0.416017559854925*x[14] - 3.07851409500399*x[9]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-3.07851409500399, 0.920590812717961, 0.00269678297809628, 0.0386749962186403, 0.703439767346247, 0.416017559854925}, {}, {}, {}, {}};
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
		return -3.8864807692196*u[0]*x[10] + 0.113984715960886*u[0]*x[11] + 0.90484184612424*u[0]*x[12] + 0.978420947419502*u[0]*x[13] + 0.526719910416832*u[0]*x[14] + 0.867803591739988*u[0]*x[15];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-3.8864807692196*u[0], 0.113984715960886*u[0], 0.90484184612424*u[0], 0.978420947419502*u[0], 0.526719910416832*u[0], 0.867803591739988*u[0]}, {-3.8864807692196*x[10] + 0.113984715960886*x[11] + 0.90484184612424*x[12] + 0.978420947419502*x[13] + 0.526719910416832*x[14] + 0.867803591739988*x[15]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-3.8864807692196, 0.113984715960886, 0.90484184612424, 0.978420947419502, 0.526719910416832, 0.867803591739988}, {}, {}, {}, {}};
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
		return -2.37327376378973*u[0]*x[11] + 0.216009898248265*u[0]*x[12] + 0.390474902221211*u[0]*x[13] + 0.466103699608691*u[0]*x[14] + 0.526363414491039*u[0]*x[15] + 0.684715008234527*u[0]*x[16];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-2.37327376378973*u[0], 0.216009898248265*u[0], 0.390474902221211*u[0], 0.466103699608691*u[0], 0.526363414491039*u[0], 0.684715008234527*u[0]}, {-2.37327376378973*x[11] + 0.216009898248265*x[12] + 0.390474902221211*x[13] + 0.466103699608691*x[14] + 0.526363414491039*x[15] + 0.684715008234527*x[16]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-2.37327376378973, 0.216009898248265, 0.390474902221211, 0.466103699608691, 0.526363414491039, 0.684715008234527}, {}, {}, {}, {}};
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
		return -2.0677704826297*u[0]*x[12] + 0.0479075005400251*u[0]*x[13] + 0.71833420240236*u[0]*x[14] + 0.462263802114936*u[0]*x[15] + 0.021586430504078*u[0]*x[16] + 0.604692485855283*u[0]*x[17];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-2.0677704826297*u[0], 0.0479075005400251*u[0], 0.71833420240236*u[0], 0.462263802114936*u[0], 0.021586430504078*u[0], 0.604692485855283*u[0]}, {-2.0677704826297*x[12] + 0.0479075005400251*x[13] + 0.71833420240236*x[14] + 0.462263802114936*x[15] + 0.021586430504078*x[16] + 0.604692485855283*x[17]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-2.0677704826297, 0.0479075005400251, 0.71833420240236, 0.462263802114936, 0.021586430504078, 0.604692485855283}, {}, {}, {}, {}};
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
		return -3.08853598925294*u[0]*x[13] + 0.716237655883514*u[0]*x[14] + 0.948998541772577*u[0]*x[15] + 0.18957632638383*u[0]*x[16] + 0.585928584363033*u[0]*x[17] + 0.287158146740105*u[0]*x[18];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-3.08853598925294*u[0], 0.716237655883514*u[0], 0.948998541772577*u[0], 0.18957632638383*u[0], 0.585928584363033*u[0], 0.287158146740105*u[0]}, {-3.08853598925294*x[13] + 0.716237655883514*x[14] + 0.948998541772577*x[15] + 0.18957632638383*x[16] + 0.585928584363033*x[17] + 0.287158146740105*x[18]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-3.08853598925294, 0.716237655883514, 0.948998541772577, 0.18957632638383, 0.585928584363033, 0.287158146740105}, {}, {}, {}, {}};
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
		return -2.84341302816632*u[0]*x[14] + 0.999171823100972*u[0]*x[15] + 0.138640185494095*u[0]*x[16] + 0.466988639517415*u[0]*x[17] + 0.684158804360738*u[0]*x[18] + 0.258252207064527*u[0]*x[19];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-2.84341302816632*u[0], 0.999171823100972*u[0], 0.138640185494095*u[0], 0.466988639517415*u[0], 0.684158804360738*u[0], 0.258252207064527*u[0]}, {-2.84341302816632*x[14] + 0.999171823100972*x[15] + 0.138640185494095*x[16] + 0.466988639517415*x[17] + 0.684158804360738*x[18] + 0.258252207064527*x[19]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-2.84341302816632, 0.999171823100972, 0.138640185494095, 0.466988639517415, 0.684158804360738, 0.258252207064527}, {}, {}, {}, {}};
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
		return -3.80460117321951*u[0]*x[15] + 0.82776996713313*u[0]*x[16] + 0.984783520247337*u[0]*x[17] + 0.0329191349684085*u[0]*x[18] + 0.524081937600162*u[0]*x[19] + 0.0207814534358743*u[0]*x[20];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-3.80460117321951*u[0], 0.82776996713313*u[0], 0.984783520247337*u[0], 0.0329191349684085*u[0], 0.524081937600162*u[0], 0.0207814534358743*u[0]}, {-3.80460117321951*x[15] + 0.82776996713313*x[16] + 0.984783520247337*x[17] + 0.0329191349684085*x[18] + 0.524081937600162*x[19] + 0.0207814534358743*x[20]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-3.80460117321951, 0.82776996713313, 0.984783520247337, 0.0329191349684085, 0.524081937600162, 0.0207814534358743}, {}, {}, {}, {}};
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
		return -1.86228791774966*u[0]*x[16] + 0.750291068731939*u[0]*x[17] + 0.00551507419663422*u[0]*x[18] + 0.499748935921602*u[0]*x[19] + 0.52589292274709*u[0]*x[20] + 0.338904991052394*u[0]*x[21];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-1.86228791774966*u[0], 0.750291068731939*u[0], 0.00551507419663422*u[0], 0.499748935921602*u[0], 0.52589292274709*u[0], 0.338904991052394*u[0]}, {-1.86228791774966*x[16] + 0.750291068731939*x[17] + 0.00551507419663422*x[18] + 0.499748935921602*x[19] + 0.52589292274709*x[20] + 0.338904991052394*x[21]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-1.86228791774966, 0.750291068731939, 0.00551507419663422, 0.499748935921602, 0.52589292274709, 0.338904991052394}, {}, {}, {}, {}};
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
		return -3.39268429871501*u[0]*x[17] + 0.842622169389521*u[0]*x[18] + 0.921508410314994*u[0]*x[19] + 0.550057218863126*u[0]*x[20] + 0.387427471760572*u[0]*x[21] + 0.58677830963635*u[0]*x[22];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-3.39268429871501*u[0], 0.842622169389521*u[0], 0.921508410314994*u[0], 0.550057218863126*u[0], 0.387427471760572*u[0], 0.58677830963635*u[0]}, {-3.39268429871501*x[17] + 0.842622169389521*x[18] + 0.921508410314994*x[19] + 0.550057218863126*x[20] + 0.387427471760572*x[21] + 0.58677830963635*x[22]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-3.39268429871501, 0.842622169389521, 0.921508410314994, 0.550057218863126, 0.387427471760572, 0.58677830963635}, {}, {}, {}, {}};
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
		return -1.85237332965541*u[0]*x[18] + 0.568103511165832*u[0]*x[19] + 0.953121626986024*u[0]*x[20] + 0.862990089285466*u[0]*x[21] + 0.177076753474592*u[0]*x[22] + 0.585347383583152*u[0]*x[23];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-1.85237332965541*u[0], 0.568103511165832*u[0], 0.953121626986024*u[0], 0.862990089285466*u[0], 0.177076753474592*u[0], 0.585347383583152*u[0]}, {-1.85237332965541*x[18] + 0.568103511165832*x[19] + 0.953121626986024*x[20] + 0.862990089285466*x[21] + 0.177076753474592*x[22] + 0.585347383583152*x[23]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-1.85237332965541, 0.568103511165832, 0.953121626986024, 0.862990089285466, 0.177076753474592, 0.585347383583152}, {}, {}, {}, {}};
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
		return -2.77169500206712*u[0]*x[19] + 0.962618292454326*u[0]*x[20] + 0.54255068456849*u[0]*x[21] + 0.64260760535423*u[0]*x[22] + 0.426599945198362*u[0]*x[23] + 0.507288262614702*u[0]*x[24];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-2.77169500206712*u[0], 0.962618292454326*u[0], 0.54255068456849*u[0], 0.64260760535423*u[0], 0.426599945198362*u[0], 0.507288262614702*u[0]}, {-2.77169500206712*x[19] + 0.962618292454326*x[20] + 0.54255068456849*x[21] + 0.64260760535423*x[22] + 0.426599945198362*x[23] + 0.507288262614702*x[24]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-2.77169500206712, 0.962618292454326, 0.54255068456849, 0.64260760535423, 0.426599945198362, 0.507288262614702}, {}, {}, {}, {}};
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
		return -3.01247151448644*u[0]*x[20] + 0.425014421029746*u[0]*x[21] + 0.633372611686713*u[0]*x[22] + 0.785404547293533*u[0]*x[23] + 0.99378340303531*u[0]*x[24] + 0.0671539683108534*u[0]*x[25];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-3.01247151448644*u[0], 0.425014421029746*u[0], 0.633372611686713*u[0], 0.785404547293533*u[0], 0.99378340303531*u[0], 0.0671539683108534*u[0]}, {-3.01247151448644*x[20] + 0.425014421029746*x[21] + 0.633372611686713*x[22] + 0.785404547293533*x[23] + 0.99378340303531*x[24] + 0.0671539683108534*x[25]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-3.01247151448644, 0.425014421029746, 0.633372611686713, 0.785404547293533, 0.99378340303531, 0.0671539683108534}, {}, {}, {}, {}};
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
		return -2.55688765769667*u[0]*x[21] + 0.929003405600682*u[0]*x[22] + 0.860824560131906*u[0]*x[23] + 0.324707124572482*u[0]*x[24] + 0.669877259556439*u[0]*x[25] + 0.226653040608637*u[0]*x[26];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-2.55688765769667*u[0], 0.929003405600682*u[0], 0.860824560131906*u[0], 0.324707124572482*u[0], 0.669877259556439*u[0], 0.226653040608637*u[0]}, {-2.55688765769667*x[21] + 0.929003405600682*x[22] + 0.860824560131906*x[23] + 0.324707124572482*x[24] + 0.669877259556439*x[25] + 0.226653040608637*x[26]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-2.55688765769667, 0.929003405600682, 0.860824560131906, 0.324707124572482, 0.669877259556439, 0.226653040608637}, {}, {}, {}, {}};
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
		return -2.96883868575257*u[0]*x[22] + 0.253842542804843*u[0]*x[23] + 0.987351853537923*u[0]*x[24] + 0.23219132632934*u[0]*x[25] + 0.886656756198478*u[0]*x[26] + 0.301798077361254*u[0]*x[27];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-2.96883868575257*u[0], 0.253842542804843*u[0], 0.987351853537923*u[0], 0.23219132632934*u[0], 0.886656756198478*u[0], 0.301798077361254*u[0]}, {-2.96883868575257*x[22] + 0.253842542804843*x[23] + 0.987351853537923*x[24] + 0.23219132632934*x[25] + 0.886656756198478*x[26] + 0.301798077361254*x[27]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-2.96883868575257, 0.253842542804843, 0.987351853537923, 0.23219132632934, 0.886656756198478, 0.301798077361254}, {}, {}, {}, {}};
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
		return -2.9120189790118*u[0]*x[23] + 0.452580447772077*u[0]*x[24] + 0.940195177543766*u[0]*x[25] + 0.0290119434783812*u[0]*x[26] + 0.0220053229696739*u[0]*x[27] + 0.40496103785304*u[0]*x[28];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-2.9120189790118*u[0], 0.452580447772077*u[0], 0.940195177543766*u[0], 0.0290119434783812*u[0], 0.0220053229696739*u[0], 0.40496103785304*u[0]}, {-2.9120189790118*x[23] + 0.452580447772077*x[24] + 0.940195177543766*x[25] + 0.0290119434783812*x[26] + 0.0220053229696739*x[27] + 0.40496103785304*x[28]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-2.9120189790118, 0.452580447772077, 0.940195177543766, 0.0290119434783812, 0.0220053229696739, 0.40496103785304}, {}, {}, {}, {}};
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
		return -3.26571109153249*u[0]*x[24] + 0.535508667376141*u[0]*x[25] + 0.979920338559838*u[0]*x[26] + 0.601910783419637*u[0]*x[27] + 0.377755140860824*u[0]*x[28] + 0.29240810451974*u[0]*x[29];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-3.26571109153249*u[0], 0.535508667376141*u[0], 0.979920338559838*u[0], 0.601910783419637*u[0], 0.377755140860824*u[0], 0.29240810451974*u[0]}, {-3.26571109153249*x[24] + 0.535508667376141*x[25] + 0.979920338559838*x[26] + 0.601910783419637*x[27] + 0.377755140860824*x[28] + 0.29240810451974*x[29]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-3.26571109153249, 0.535508667376141, 0.979920338559838, 0.601910783419637, 0.377755140860824, 0.29240810451974}, {}, {}, {}, {}};
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
		return -2.44492639911654*u[0]*x[25] + 0.502878622104055*u[0]*x[26] + 0.0123321407571459*u[0]*x[27] + 0.543368298942586*u[0]*x[28] + 0.873861850040888*u[0]*x[29] + 0.223390770776782*u[0]*x[30];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-2.44492639911654*u[0], 0.502878622104055*u[0], 0.0123321407571459*u[0], 0.543368298942586*u[0], 0.873861850040888*u[0], 0.223390770776782*u[0]}, {-2.44492639911654*x[25] + 0.502878622104055*x[26] + 0.0123321407571459*x[27] + 0.543368298942586*x[28] + 0.873861850040888*x[29] + 0.223390770776782*x[30]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-2.44492639911654, 0.502878622104055, 0.0123321407571459, 0.543368298942586, 0.873861850040888, 0.223390770776782}, {}, {}, {}, {}};
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
		return -2.62512070094939*u[0]*x[26] + 0.623236602647147*u[0]*x[27] + 0.711088734434572*u[0]*x[28] + 0.00541756504573609*u[0]*x[29] + 0.641638748486164*u[0]*x[30] + 0.994941064030392*u[0]*x[31];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-2.62512070094939*u[0], 0.623236602647147*u[0], 0.711088734434572*u[0], 0.00541756504573609*u[0], 0.641638748486164*u[0], 0.994941064030392*u[0]}, {-2.62512070094939*x[26] + 0.623236602647147*x[27] + 0.711088734434572*x[28] + 0.00541756504573609*x[29] + 0.641638748486164*x[30] + 0.994941064030392*x[31]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-2.62512070094939, 0.623236602647147, 0.711088734434572, 0.00541756504573609, 0.641638748486164, 0.994941064030392}, {}, {}, {}, {}};
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
		return -1.56128292715486*u[0]*x[27] + 0.467327187822675*u[0]*x[28] + 0.0473254600265938*u[0]*x[29] + 0.598153192449048*u[0]*x[30] + 0.778612018302338*u[0]*x[31] + 0.737222343048654*u[0]*x[32];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-1.56128292715486*u[0], 0.467327187822675*u[0], 0.0473254600265938*u[0], 0.598153192449048*u[0], 0.778612018302338*u[0], 0.737222343048654*u[0]}, {-1.56128292715486*x[27] + 0.467327187822675*x[28] + 0.0473254600265938*x[29] + 0.598153192449048*x[30] + 0.778612018302338*x[31] + 0.737222343048654*x[32]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-1.56128292715486, 0.467327187822675, 0.0473254600265938, 0.598153192449048, 0.778612018302338, 0.737222343048654}, {}, {}, {}, {}};
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
		return -2.5045003999137*u[0]*x[28] + 0.336385380909935*u[0]*x[29] + 0.223723489811174*u[0]*x[30] + 0.566207863242435*u[0]*x[31] + 0.0281711342982055*u[0]*x[32] + 0.909853802830485*u[0]*x[33];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-2.5045003999137*u[0], 0.336385380909935*u[0], 0.223723489811174*u[0], 0.566207863242435*u[0], 0.0281711342982055*u[0], 0.909853802830485*u[0]}, {-2.5045003999137*x[28] + 0.336385380909935*x[29] + 0.223723489811174*x[30] + 0.566207863242435*x[31] + 0.0281711342982055*x[32] + 0.909853802830485*x[33]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-2.5045003999137, 0.336385380909935, 0.223723489811174, 0.566207863242435, 0.0281711342982055, 0.909853802830485}, {}, {}, {}, {}};
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
		return -1.55539836054289*u[0]*x[29] + 0.643951602420491*u[0]*x[30] + 0.946434807589009*u[0]*x[31] + 0.675713859647608*u[0]*x[32] + 0.100788670306937*u[0]*x[33] + 0.778970392565905*u[0]*x[34];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-1.55539836054289*u[0], 0.643951602420491*u[0], 0.946434807589009*u[0], 0.675713859647608*u[0], 0.100788670306937*u[0], 0.778970392565905*u[0]}, {-1.55539836054289*x[29] + 0.643951602420491*x[30] + 0.946434807589009*x[31] + 0.675713859647608*x[32] + 0.100788670306937*x[33] + 0.778970392565905*x[34]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-1.55539836054289, 0.643951602420491, 0.946434807589009, 0.675713859647608, 0.100788670306937, 0.778970392565905}, {}, {}, {}, {}};
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
		return -2.33085780394366*u[0]*x[30] + 0.219794325177768*u[0]*x[31] + 0.379470624360871*u[0]*x[32] + 0.359089532853284*u[0]*x[33] + 0.24318136925114*u[0]*x[34] + 0.84555852935173*u[0]*x[35];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-2.33085780394366*u[0], 0.219794325177768*u[0], 0.379470624360871*u[0], 0.359089532853284*u[0], 0.24318136925114*u[0], 0.84555852935173*u[0]}, {-2.33085780394366*x[30] + 0.219794325177768*x[31] + 0.379470624360871*x[32] + 0.359089532853284*x[33] + 0.24318136925114*x[34] + 0.84555852935173*x[35]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-2.33085780394366, 0.219794325177768, 0.379470624360871, 0.359089532853284, 0.24318136925114, 0.84555852935173}, {}, {}, {}, {}};
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
		return -3.50599007834194*u[0]*x[31] + 0.689283294952872*u[0]*x[32] + 0.531690873090384*u[0]*x[33] + 0.395667855647068*u[0]*x[34] + 0.593310542253176*u[0]*x[35] + 0.389330886806684*u[0]*x[36];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-3.50599007834194*u[0], 0.689283294952872*u[0], 0.531690873090384*u[0], 0.395667855647068*u[0], 0.593310542253176*u[0], 0.389330886806684*u[0]}, {-3.50599007834194*x[31] + 0.689283294952872*x[32] + 0.531690873090384*x[33] + 0.395667855647068*x[34] + 0.593310542253176*x[35] + 0.389330886806684*x[36]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-3.50599007834194, 0.689283294952872, 0.531690873090384, 0.395667855647068, 0.593310542253176, 0.389330886806684}, {}, {}, {}, {}};
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
		return -2.50986125630821*u[0]*x[32] + 0.754638610961782*u[0]*x[33] + 0.18592173947098*u[0]*x[34] + 0.511976331263828*u[0]*x[35] + 0.267008360951234*u[0]*x[36] + 0.846280395699658*u[0]*x[37];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-2.50986125630821*u[0], 0.754638610961782*u[0], 0.18592173947098*u[0], 0.511976331263828*u[0], 0.267008360951234*u[0], 0.846280395699658*u[0]}, {-2.50986125630821*x[32] + 0.754638610961782*x[33] + 0.18592173947098*x[34] + 0.511976331263828*x[35] + 0.267008360951234*x[36] + 0.846280395699658*x[37]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-2.50986125630821, 0.754638610961782, 0.18592173947098, 0.511976331263828, 0.267008360951234, 0.846280395699658}, {}, {}, {}, {}};
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
		return -2.65606149004287*u[0]*x[33] + 0.487997303932508*u[0]*x[34] + 0.554227584986253*u[0]*x[35] + 0.604336575180743*u[0]*x[36] + 0.522932435754969*u[0]*x[37] + 0.934890315732395*u[0]*x[38];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-2.65606149004287*u[0], 0.487997303932508*u[0], 0.554227584986253*u[0], 0.604336575180743*u[0], 0.522932435754969*u[0], 0.934890315732395*u[0]}, {-2.65606149004287*x[33] + 0.487997303932508*x[34] + 0.554227584986253*x[35] + 0.604336575180743*x[36] + 0.522932435754969*x[37] + 0.934890315732395*x[38]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-2.65606149004287, 0.487997303932508, 0.554227584986253, 0.604336575180743, 0.522932435754969, 0.934890315732395}, {}, {}, {}, {}};
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
		return -2.0917386608676*u[0]*x[34] + 0.0624376197492708*u[0]*x[35] + 0.752028209954402*u[0]*x[36] + 0.843896934152684*u[0]*x[37] + 0.927291669498414*u[0]*x[38] + 0.802268789440801*u[0]*x[39];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-2.0917386608676*u[0], 0.0624376197492708*u[0], 0.752028209954402*u[0], 0.843896934152684*u[0], 0.927291669498414*u[0], 0.802268789440801*u[0]}, {-2.0917386608676*x[34] + 0.0624376197492708*x[35] + 0.752028209954402*x[36] + 0.843896934152684*x[37] + 0.927291669498414*x[38] + 0.802268789440801*x[39]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-2.0917386608676, 0.0624376197492708, 0.752028209954402, 0.843896934152684, 0.927291669498414, 0.802268789440801}, {}, {}, {}, {}};
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
		return -2.56751060760426*u[0]*x[35] + 0.346848229688394*u[0]*x[36] + 0.958475139959143*u[0]*x[37] + 0.391598980834086*u[0]*x[38] + 0.917987901068765*u[0]*x[39] + 0.582120280527124*u[0]*x[40];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-2.56751060760426*u[0], 0.346848229688394*u[0], 0.958475139959143*u[0], 0.391598980834086*u[0], 0.917987901068765*u[0], 0.582120280527124*u[0]}, {-2.56751060760426*x[35] + 0.346848229688394*x[36] + 0.958475139959143*x[37] + 0.391598980834086*x[38] + 0.917987901068765*x[39] + 0.582120280527124*x[40]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-2.56751060760426, 0.346848229688394, 0.958475139959143, 0.391598980834086, 0.917987901068765, 0.582120280527124}, {}, {}, {}, {}};
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
		return -2.35955226258146*u[0]*x[36] + 0.786280369824922*u[0]*x[37] + 0.425582125725915*u[0]*x[38] + 0.558937026389117*u[0]*x[39] + 0.670252662129815*u[0]*x[40] + 0.639345708288811*u[0]*x[41];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-2.35955226258146*u[0], 0.786280369824922*u[0], 0.425582125725915*u[0], 0.558937026389117*u[0], 0.670252662129815*u[0], 0.639345708288811*u[0]}, {-2.35955226258146*x[36] + 0.786280369824922*x[37] + 0.425582125725915*x[38] + 0.558937026389117*x[39] + 0.670252662129815*x[40] + 0.639345708288811*x[41]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-2.35955226258146, 0.786280369824922, 0.425582125725915, 0.558937026389117, 0.670252662129815, 0.639345708288811}, {}, {}, {}, {}};
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
		return -3.95786527539138*u[0]*x[37] + 0.875024885349743*u[0]*x[38] + 0.491662216522932*u[0]*x[39] + 0.234043253861775*u[0]*x[40] + 0.831805167011595*u[0]*x[41] + 0.102893058073924*u[0]*x[42];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-3.95786527539138*u[0], 0.875024885349743*u[0], 0.491662216522932*u[0], 0.234043253861775*u[0], 0.831805167011595*u[0], 0.102893058073924*u[0]}, {-3.95786527539138*x[37] + 0.875024885349743*x[38] + 0.491662216522932*x[39] + 0.234043253861775*x[40] + 0.831805167011595*x[41] + 0.102893058073924*x[42]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-3.95786527539138, 0.875024885349743, 0.491662216522932, 0.234043253861775, 0.831805167011595, 0.102893058073924}, {}, {}, {}, {}};
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
		return -3.55438797714055*u[0]*x[38] + 0.481802880546061*u[0]*x[39] + 0.618571583043053*u[0]*x[40] + 0.895473492220389*u[0]*x[41] + 0.839840496507866*u[0]*x[42] + 0.752745394759449*u[0]*x[43];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-3.55438797714055*u[0], 0.481802880546061*u[0], 0.618571583043053*u[0], 0.895473492220389*u[0], 0.839840496507866*u[0], 0.752745394759449*u[0]}, {-3.55438797714055*x[38] + 0.481802880546061*x[39] + 0.618571583043053*x[40] + 0.895473492220389*x[41] + 0.839840496507866*x[42] + 0.752745394759449*x[43]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-3.55438797714055, 0.481802880546061, 0.618571583043053, 0.895473492220389, 0.839840496507866, 0.752745394759449}, {}, {}, {}, {}};
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
		return -3.25265881396768*u[0]*x[39] + 0.94984910924729*u[0]*x[40] + 0.267650029345338*u[0]*x[41] + 0.954007908729337*u[0]*x[42] + 0.542885229063261*u[0]*x[43] + 0.135670150072439*u[0]*x[44];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-3.25265881396768*u[0], 0.94984910924729*u[0], 0.267650029345338*u[0], 0.954007908729337*u[0], 0.542885229063261*u[0], 0.135670150072439*u[0]}, {-3.25265881396768*x[39] + 0.94984910924729*x[40] + 0.267650029345338*x[41] + 0.954007908729337*x[42] + 0.542885229063261*x[43] + 0.135670150072439*x[44]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-3.25265881396768, 0.94984910924729, 0.267650029345338, 0.954007908729337, 0.542885229063261, 0.135670150072439}, {}, {}, {}, {}};
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
		return -3.05483688880906*u[0]*x[40] + 0.892431859205601*u[0]*x[41] + 0.0528374640872615*u[0]*x[42] + 0.632775085209332*u[0]*x[43] + 0.187879215268798*u[0]*x[44] + 0.504748824792223*u[0]*x[45];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-3.05483688880906*u[0], 0.892431859205601*u[0], 0.0528374640872615*u[0], 0.632775085209332*u[0], 0.187879215268798*u[0], 0.504748824792223*u[0]}, {-3.05483688880906*x[40] + 0.892431859205601*x[41] + 0.0528374640872615*x[42] + 0.632775085209332*x[43] + 0.187879215268798*x[44] + 0.504748824792223*x[45]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-3.05483688880906, 0.892431859205601, 0.0528374640872615, 0.632775085209332, 0.187879215268798, 0.504748824792223}, {}, {}, {}, {}};
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
		return -3.52670625607173*u[0]*x[41] + 0.167993096003758*u[0]*x[42] + 0.713247801151574*u[0]*x[43] + 0.523850846290651*u[0]*x[44] + 0.1776396288159*u[0]*x[45] + 0.701074960352789*u[0]*x[46];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-3.52670625607173*u[0], 0.167993096003758*u[0], 0.713247801151574*u[0], 0.523850846290651*u[0], 0.1776396288159*u[0], 0.701074960352789*u[0]}, {-3.52670625607173*x[41] + 0.167993096003758*x[42] + 0.713247801151574*x[43] + 0.523850846290651*x[44] + 0.1776396288159*x[45] + 0.701074960352789*x[46]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-3.52670625607173, 0.167993096003758, 0.713247801151574, 0.523850846290651, 0.1776396288159, 0.701074960352789}, {}, {}, {}, {}};
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
		return -2.11757202340215*u[0]*x[42] + 0.977285576116882*u[0]*x[43] + 0.515441472328528*u[0]*x[44] + 0.954886730382985*u[0]*x[45] + 0.386967980766608*u[0]*x[46] + 0.499534112104472*u[0]*x[47];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-2.11757202340215*u[0], 0.977285576116882*u[0], 0.515441472328528*u[0], 0.954886730382985*u[0], 0.386967980766608*u[0], 0.499534112104472*u[0]}, {-2.11757202340215*x[42] + 0.977285576116882*x[43] + 0.515441472328528*x[44] + 0.954886730382985*x[45] + 0.386967980766608*x[46] + 0.499534112104472*x[47]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-2.11757202340215, 0.977285576116882, 0.515441472328528, 0.954886730382985, 0.386967980766608, 0.499534112104472}, {}, {}, {}, {}};
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
		return -3.6189390863005*u[0]*x[43] + 0.265601479823131*u[0]*x[44] + 0.906602643729759*u[0]*x[45] + 0.957178595368757*u[0]*x[46] + 0.193376030748537*u[0]*x[47] + 0.936527838667765*u[0]*x[48];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-3.6189390863005*u[0], 0.265601479823131*u[0], 0.906602643729759*u[0], 0.957178595368757*u[0], 0.193376030748537*u[0], 0.936527838667765*u[0]}, {-3.6189390863005*x[43] + 0.265601479823131*x[44] + 0.906602643729759*x[45] + 0.957178595368757*x[46] + 0.193376030748537*x[47] + 0.936527838667765*x[48]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-3.6189390863005, 0.265601479823131, 0.906602643729759, 0.957178595368757, 0.193376030748537, 0.936527838667765}, {}, {}, {}, {}};
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
		return -1.62844316378355*u[0]*x[44] + 0.519858769674969*u[0]*x[45] + 0.636526043363135*u[0]*x[46] + 0.63956768992169*u[0]*x[47] + 0.833888855090984*u[0]*x[48] + 0.272850130040555*u[0]*x[49];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-1.62844316378355*u[0], 0.519858769674969*u[0], 0.636526043363135*u[0], 0.63956768992169*u[0], 0.833888855090984*u[0], 0.272850130040555*u[0]}, {-1.62844316378355*x[44] + 0.519858769674969*x[45] + 0.636526043363135*x[46] + 0.63956768992169*x[47] + 0.833888855090984*x[48] + 0.272850130040555*x[49]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-1.62844316378355, 0.519858769674969, 0.636526043363135, 0.63956768992169, 0.833888855090984, 0.272850130040555}, {}, {}, {}, {}};
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
		return -3.06373659739584*u[0]*x[45] + 0.956185707987931*u[0]*x[46] + 0.467861544853982*u[0]*x[47] + 0.694426499970916*u[0]*x[48] + 0.156454543033952*u[0]*x[49];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-3.06373659739584*u[0], 0.956185707987931*u[0], 0.467861544853982*u[0], 0.694426499970916*u[0], 0.156454543033952*u[0]}, {-3.06373659739584*x[45] + 0.956185707987931*x[46] + 0.467861544853982*x[47] + 0.694426499970916*x[48] + 0.156454543033952*x[49]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-3.06373659739584, 0.956185707987931, 0.467861544853982, 0.694426499970916, 0.156454543033952}, {}, {}, {}, {}};
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
		return -3.63793328783922*u[0]*x[46] + 0.103367585181001*u[0]*x[47] + 0.862729113134381*u[0]*x[48] + 0.579848462800823*u[0]*x[49];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-3.63793328783922*u[0], 0.103367585181001*u[0], 0.862729113134381*u[0], 0.579848462800823*u[0]}, {-3.63793328783922*x[46] + 0.103367585181001*x[47] + 0.862729113134381*x[48] + 0.579848462800823*x[49]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-3.63793328783922, 0.103367585181001, 0.862729113134381, 0.579848462800823}, {}, {}, {}, {}};
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
		return -1.90370696280968*u[0]*x[47] + 0.912740194270146*u[0]*x[48] + 0.734600643736961*u[0]*x[49];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-1.90370696280968*u[0], 0.912740194270146*u[0], 0.734600643736961*u[0]}, {-1.90370696280968*x[47] + 0.912740194270146*x[48] + 0.734600643736961*x[49]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-1.90370696280968, 0.912740194270146, 0.734600643736961}, {}, {}, {}, {}};
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
		return -4.24031250113419*u[0]*x[48] + 0.556455839315072*u[0]*x[49];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-4.24031250113419*u[0], 0.556455839315072*u[0]}, {-4.24031250113419*x[48] + 0.556455839315072*x[49]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-4.24031250113419, 0.556455839315072}, {}, {}, {}, {}};
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
		return -2.30020961892736*u[0]*x[49] + pow(u[0], 2)*x[0]*DEPLETION_COEFF_VALUE;
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{pow(u[0], 2)*DEPLETION_COEFF_VALUE, -2.30020961892736*u[0]}, {-2.30020961892736*x[49] + 2*u[0]*x[0]*DEPLETION_COEFF_VALUE}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        const double x0 = 2*DEPLETION_COEFF_VALUE;
        const double x1 = x0*u[0];
		return {std::vector<double>{}, {x1, -2.30020961892736}, {x0*x[0]}, {}, {}, {}};
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
        

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
		return 0.468626387712207*u[0]*x[1] + 0.136548656087503*u[0]*x[2] + 0.0139353590234034*u[0]*x[3] + 0.356477747573193*u[0]*x[4] + 0.515097894100223*u[0]*x[5] - pow(u[0], 2)*x[0]*DEPLETION_COEFF_VALUE;
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-pow(u[0], 2)*DEPLETION_COEFF_VALUE, 0.468626387712207*u[0], 0.136548656087503*u[0], 0.0139353590234034*u[0], 0.356477747573193*u[0], 0.515097894100223*u[0]}, {0.468626387712207*x[1] + 0.136548656087503*x[2] + 0.0139353590234034*x[3] + 0.356477747573193*x[4] + 0.515097894100223*x[5] - 2*u[0]*x[0]*DEPLETION_COEFF_VALUE}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        const double x0 = 2*DEPLETION_COEFF_VALUE;
        const double x1 = -x0*u[0];
		return {std::vector<double>{}, {x1, 0.468626387712207, 0.136548656087503, 0.0139353590234034, 0.356477747573193, 0.515097894100223}, {-x0*x[0]}, {}, {}, {}};
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
		return -0.468626387712207*u[0]*x[1] + 0.00701438597494297*u[0]*x[2] + 0.253673105694702*u[0]*x[3] + 0.339182778727429*u[0]*x[4] + 0.901804331590947*u[0]*x[5] + 0.634022260001245*u[0]*x[6];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-0.468626387712207*u[0], 0.00701438597494297*u[0], 0.253673105694702*u[0], 0.339182778727429*u[0], 0.901804331590947*u[0], 0.634022260001245*u[0]}, {-0.468626387712207*x[1] + 0.00701438597494297*x[2] + 0.253673105694702*x[3] + 0.339182778727429*x[4] + 0.901804331590947*x[5] + 0.634022260001245*x[6]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-0.468626387712207, 0.00701438597494297, 0.253673105694702, 0.339182778727429, 0.901804331590947, 0.634022260001245}, {}, {}, {}, {}};
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
		return -0.143563042062446*u[0]*x[2] + 0.902993845468682*u[0]*x[3] + 0.630259341968257*u[0]*x[4] + 0.339487012156792*u[0]*x[5] + 0.505216709540402*u[0]*x[6] + 0.72263736822388*u[0]*x[7];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-0.143563042062446*u[0], 0.902993845468682*u[0], 0.630259341968257*u[0], 0.339487012156792*u[0], 0.505216709540402*u[0], 0.72263736822388*u[0]}, {-0.143563042062446*x[2] + 0.902993845468682*x[3] + 0.630259341968257*x[4] + 0.339487012156792*x[5] + 0.505216709540402*x[6] + 0.72263736822388*x[7]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-0.143563042062446, 0.902993845468682, 0.630259341968257, 0.339487012156792, 0.505216709540402, 0.72263736822388}, {}, {}, {}, {}};
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
		return -1.17060231018679*u[0]*x[3] + 0.650678950681845*u[0]*x[4] + 0.138731058137815*u[0]*x[5] + 0.214337498464053*u[0]*x[6] + 0.310717650919698*u[0]*x[7] + 0.0604732888684647*u[0]*x[8];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-1.17060231018679*u[0], 0.650678950681845*u[0], 0.138731058137815*u[0], 0.214337498464053*u[0], 0.310717650919698*u[0], 0.0604732888684647*u[0]}, {-1.17060231018679*x[3] + 0.650678950681845*x[4] + 0.138731058137815*x[5] + 0.214337498464053*x[6] + 0.310717650919698*x[7] + 0.0604732888684647*x[8]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-1.17060231018679, 0.650678950681845, 0.138731058137815, 0.214337498464053, 0.310717650919698, 0.0604732888684647}, {}, {}, {}, {}};
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
		return -1.97659881895072*u[0]*x[4] + 0.569712563480458*u[0]*x[5] + 0.610981910215908*u[0]*x[6] + 0.230121239588827*u[0]*x[7] + 0.688226293032254*u[0]*x[8] + 0.0778287945787205*u[0]*x[9];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-1.97659881895072*u[0], 0.569712563480458*u[0], 0.610981910215908*u[0], 0.230121239588827*u[0], 0.688226293032254*u[0], 0.0778287945787205*u[0]}, {-1.97659881895072*x[4] + 0.569712563480458*x[5] + 0.610981910215908*x[6] + 0.230121239588827*x[7] + 0.688226293032254*x[8] + 0.0778287945787205*x[9]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-1.97659881895072, 0.569712563480458, 0.610981910215908, 0.230121239588827, 0.688226293032254, 0.0778287945787205}, {}, {}, {}, {}};
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
		return 0.606201074129241*u[0]*x[10] - 2.46483285946623*u[0]*x[5] + 0.0581594677766095*u[0]*x[6] + 0.848721089452619*u[0]*x[7] + 0.3014397181145*u[0]*x[8] + 0.875160059335475*u[0]*x[9];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-2.46483285946623*u[0], 0.0581594677766095*u[0], 0.848721089452619*u[0], 0.3014397181145*u[0], 0.875160059335475*u[0], 0.606201074129241*u[0]}, {0.606201074129241*x[10] - 2.46483285946623*x[5] + 0.0581594677766095*x[6] + 0.848721089452619*x[7] + 0.3014397181145*x[8] + 0.875160059335475*x[9]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-2.46483285946623, 0.0581594677766095, 0.848721089452619, 0.3014397181145, 0.875160059335475, 0.606201074129241}, {}, {}, {}, {}};
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
		return 0.835561269311991*u[0]*x[10] + 0.690469791423505*u[0]*x[11] - 2.02271784599822*u[0]*x[6] + 0.977156835213313*u[0]*x[7] + 0.461396259354151*u[0]*x[8] + 0.177943652377825*u[0]*x[9];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-2.02271784599822*u[0], 0.977156835213313*u[0], 0.461396259354151*u[0], 0.177943652377825*u[0], 0.835561269311991*u[0], 0.690469791423505*u[0]}, {0.835561269311991*x[10] + 0.690469791423505*x[11] - 2.02271784599822*x[6] + 0.977156835213313*x[7] + 0.461396259354151*x[8] + 0.177943652377825*x[9]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-2.02271784599822, 0.977156835213313, 0.461396259354151, 0.177943652377825, 0.835561269311991, 0.690469791423505}, {}, {}, {}, {}};
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
		return 0.148841260522142*u[0]*x[10] + 0.486515266758492*u[0]*x[11] + 0.907201671316369*u[0]*x[12] - 3.08935418339834*u[0]*x[7] + 0.999308576469669*u[0]*x[8] + 0.552293164508093*u[0]*x[9];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-3.08935418339834*u[0], 0.999308576469669*u[0], 0.552293164508093*u[0], 0.148841260522142*u[0], 0.486515266758492*u[0], 0.907201671316369*u[0]}, {0.148841260522142*x[10] + 0.486515266758492*x[11] + 0.907201671316369*x[12] - 3.08935418339834*x[7] + 0.999308576469669*x[8] + 0.552293164508093*x[9]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-3.08935418339834, 0.999308576469669, 0.552293164508093, 0.148841260522142, 0.486515266758492, 0.907201671316369}, {}, {}, {}, {}};
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
		return 0.352068089300232*u[0]*x[10] + 0.709214561130748*u[0]*x[11] + 0.72286994391235*u[0]*x[12] + 0.755925612287257*u[0]*x[13] - 2.51084413583904*u[0]*x[8] + 0.638905629452865*u[0]*x[9];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-2.51084413583904*u[0], 0.638905629452865*u[0], 0.352068089300232*u[0], 0.709214561130748*u[0], 0.72286994391235*u[0], 0.755925612287257*u[0]}, {0.352068089300232*x[10] + 0.709214561130748*x[11] + 0.72286994391235*x[12] + 0.755925612287257*x[13] - 2.51084413583904*x[8] + 0.638905629452865*x[9]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-2.51084413583904, 0.638905629452865, 0.352068089300232, 0.709214561130748, 0.72286994391235, 0.755925612287257}, {}, {}, {}, {}};
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
		return 0.926368764953382*u[0]*x[10] + 0.87103924839005*u[0]*x[11] + 0.25537115819334*u[0]*x[12] + 0.691679981048153*u[0]*x[13] + 0.383857512926589*u[0]*x[14] - 2.32213130025298*u[0]*x[9];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-2.32213130025298*u[0], 0.926368764953382*u[0], 0.87103924839005*u[0], 0.25537115819334*u[0], 0.691679981048153*u[0], 0.383857512926589*u[0]}, {0.926368764953382*x[10] + 0.87103924839005*x[11] + 0.25537115819334*x[12] + 0.691679981048153*x[13] + 0.383857512926589*x[14] - 2.32213130025298*x[9]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-2.32213130025298, 0.926368764953382, 0.87103924839005, 0.25537115819334, 0.691679981048153, 0.383857512926589}, {}, {}, {}, {}};
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
		return -2.86904045821699*u[0]*x[10] + 0.149120491032765*u[0]*x[11] + 0.484456198031775*u[0]*x[12] + 0.955945156972338*u[0]*x[13] + 0.559974704400233*u[0]*x[14];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-2.86904045821699*u[0], 0.149120491032765*u[0], 0.484456198031775*u[0], 0.955945156972338*u[0], 0.559974704400233*u[0]}, {-2.86904045821699*x[10] + 0.149120491032765*x[11] + 0.484456198031775*x[12] + 0.955945156972338*x[13] + 0.559974704400233*x[14]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-2.86904045821699, 0.149120491032765, 0.484456198031775, 0.955945156972338, 0.559974704400233}, {}, {}, {}, {}};
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
		return -2.90635935873556*u[0]*x[11] + 0.518678830913885*u[0]*x[12] + 0.89471883491316*u[0]*x[13] + 0.214765545524678*u[0]*x[14];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-2.90635935873556*u[0], 0.518678830913885*u[0], 0.89471883491316*u[0], 0.214765545524678*u[0]}, {-2.90635935873556*x[11] + 0.518678830913885*x[12] + 0.89471883491316*x[13] + 0.214765545524678*x[14]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-2.90635935873556, 0.518678830913885, 0.89471883491316, 0.214765545524678}, {}, {}, {}, {}};
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
		return -2.88857780236772*u[0]*x[12] + 0.896095933755998*u[0]*x[13] + 0.28741826294945*u[0]*x[14];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-2.88857780236772*u[0], 0.896095933755998*u[0], 0.28741826294945*u[0]}, {-2.88857780236772*x[12] + 0.896095933755998*x[13] + 0.28741826294945*x[14]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-2.88857780236772, 0.896095933755998, 0.28741826294945}, {}, {}, {}, {}};
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
		return -4.19436551897691*u[0]*x[13] + 0.72385370267454*u[0]*x[14];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-4.19436551897691*u[0], 0.72385370267454*u[0]}, {-4.19436551897691*x[13] + 0.72385370267454*x[14]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-4.19436551897691, 0.72385370267454}, {}, {}, {}, {}};
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
		return -2.16986972847549*u[0]*x[14] + pow(u[0], 2)*x[0]*DEPLETION_COEFF_VALUE;
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{pow(u[0], 2)*DEPLETION_COEFF_VALUE, -2.16986972847549*u[0]}, {-2.16986972847549*x[14] + 2*u[0]*x[0]*DEPLETION_COEFF_VALUE}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        const double x0 = 2*DEPLETION_COEFF_VALUE;
        const double x1 = x0*u[0];
		return {std::vector<double>{}, {x1, -2.16986972847549}, {x0*x[0]}, {}, {}, {}};
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
        
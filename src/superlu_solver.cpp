#include "quake/superlu_solver.hpp"

#include <slu_ddefs.h>
#include <stdexcept>
#include <utility>

namespace quake {

struct SuperLUFactor::Impl {
    int n{};
    SuperMatrix L{};
    SuperMatrix U{};
    std::vector<int> perm_r;
    std::vector<int> perm_c;
    bool factored{false};

    ~Impl() {
        if (L.Store) Destroy_SuperNode_Matrix(&L);
        if (U.Store) Destroy_CompCol_Matrix(&U);
    }
};

static SuperMatrix make_csc_view(const SparseMatrixCSC& matrix) {
    SuperMatrix A{};
    // SuperLU's creation API is not const-correct; dgssv does not require us to
    // mutate the caller's numerical values for this use.
    dCreate_CompCol_Matrix(&A, matrix.rows(), matrix.cols(), matrix.nnz(),
                           const_cast<double*>(matrix.values().data()),
                           const_cast<int*>(matrix.row_ind().data()),
                           const_cast<int*>(matrix.col_ptr().data()),
                           SLU_NC, SLU_D, SLU_GE);
    return A;
}

SuperLUFactor::SuperLUFactor(const SparseMatrixCSC& matrix)
    : impl_(std::make_unique<Impl>()) {
    if (matrix.rows() != matrix.cols()) throw std::invalid_argument("Factor requires square matrix");
    impl_->n = matrix.rows();
    impl_->perm_r.resize(static_cast<std::size_t>(impl_->n));
    impl_->perm_c.resize(static_cast<std::size_t>(impl_->n));

    SuperMatrix A = make_csc_view(matrix);
    std::vector<double> dummy(static_cast<std::size_t>(impl_->n), 0.0);
    SuperMatrix B{};
    dCreate_Dense_Matrix(&B, impl_->n, 1, dummy.data(), impl_->n,
                         SLU_DN, SLU_D, SLU_GE);

    superlu_options_t options;
    set_default_options(&options);
    options.ColPerm = COLAMD;
    SuperLUStat_t stat;
    StatInit(&stat);
    int info = 0;
    dgssv(&options, &A, impl_->perm_c.data(), impl_->perm_r.data(),
          &impl_->L, &impl_->U, &B, &stat, &info);
    StatFree(&stat);
    Destroy_SuperMatrix_Store(&A);
    Destroy_SuperMatrix_Store(&B);
    if (info != 0) throw std::runtime_error("SuperLU factorization failed, info=" + std::to_string(info));
    impl_->factored = true;
}

SuperLUFactor::~SuperLUFactor() = default;
SuperLUFactor::SuperLUFactor(SuperLUFactor&&) noexcept = default;
SuperLUFactor& SuperLUFactor::operator=(SuperLUFactor&&) noexcept = default;

int SuperLUFactor::size() const { return impl_->n; }

std::vector<double> SuperLUFactor::solve(const std::vector<double>& rhs) const {
    return solve_multiple(rhs, 1);
}

std::vector<double> SuperLUFactor::solve_multiple(const std::vector<double>& rhs, int nrhs) const {
    if (nrhs <= 0 || static_cast<int>(rhs.size()) != impl_->n * nrhs) {
        throw std::invalid_argument("SuperLU RHS dimension mismatch");
    }
    std::vector<double> x = rhs;
    SuperMatrix B{};
    dCreate_Dense_Matrix(&B, impl_->n, nrhs, x.data(), impl_->n,
                         SLU_DN, SLU_D, SLU_GE);
    SuperLUStat_t stat;
    StatInit(&stat);
    int info = 0;
    dgstrs(NOTRANS, &impl_->L, &impl_->U, impl_->perm_c.data(), impl_->perm_r.data(),
           &B, &stat, &info);
    StatFree(&stat);
    Destroy_SuperMatrix_Store(&B);
    if (info != 0) throw std::runtime_error("SuperLU triangular solve failed, info=" + std::to_string(info));
    return x;
}


struct SuperLUSamePatternSolver::Impl {
    int n{};
    std::vector<int> pattern_col_ptr;
    std::vector<int> pattern_row_ind;
    std::vector<int> perm_r;
    std::vector<int> perm_c;
    std::vector<int> etree;
    std::vector<double> R;
    std::vector<double> C;
    char equed[1]{'N'};
    SuperMatrix L{};
    SuperMatrix U{};
    GlobalLU_t Glu{};
    bool factored{false};
    std::size_t factorization_count{0};

    void destroy_factors() noexcept {
        if (L.Store) Destroy_SuperNode_Matrix(&L);
        if (U.Store) Destroy_CompCol_Matrix(&U);
        L = {};
        U = {};
        Glu = {};
        factored = false;
    }

    ~Impl() { destroy_factors(); }
};

static void validate_same_pattern(const SparseMatrixCSC& matrix,
                                  const SuperLUSamePatternSolver::Impl& impl) {
    if (matrix.rows()!=impl.n || matrix.cols()!=impl.n ||
        matrix.col_ptr()!=impl.pattern_col_ptr || matrix.row_ind()!=impl.pattern_row_ind) {
        throw std::invalid_argument("SuperLU same-pattern solver received different sparsity pattern");
    }
}

static std::vector<double> dgssvx_factor_solve(SuperLUSamePatternSolver::Impl& impl,
                                               const SparseMatrixCSC& matrix,
                                               const std::vector<double>& rhs,
                                               fact_t fact_mode) {
    if (static_cast<int>(rhs.size()) != impl.n) throw std::invalid_argument("SuperLU RHS dimension mismatch");

    // SuperLU's SamePattern mode reuses only the column permutation and
    // elimination tree. It produces a new L/U factorization. Destroy the
    // previous factors before the call so their stores are not overwritten
    // and leaked. SamePattern_SameRowPerm is the mode that reuses L/U storage.
    if (fact_mode == SamePattern && impl.factored) impl.destroy_factors();

    std::vector<double> avals=matrix.values();
    std::vector<int> rows=matrix.row_ind();
    std::vector<int> cols=matrix.col_ptr();
    SuperMatrix A{};
    dCreate_CompCol_Matrix(&A,impl.n,impl.n,matrix.nnz(),avals.data(),rows.data(),cols.data(),SLU_NC,SLU_D,SLU_GE);
    std::vector<double> b=rhs, x(static_cast<std::size_t>(impl.n),0.0);
    SuperMatrix B{}, X{};
    dCreate_Dense_Matrix(&B,impl.n,1,b.data(),impl.n,SLU_DN,SLU_D,SLU_GE);
    dCreate_Dense_Matrix(&X,impl.n,1,x.data(),impl.n,SLU_DN,SLU_D,SLU_GE);

    superlu_options_t options; set_default_options(&options);
    options.Fact=fact_mode;
    options.PrintStat=NO;
    // Disable equilibration so solve_current() can use dgstrs directly on the
    // retained factors without an additional scaling/unscaling stage.
    options.Equil=NO;
    SuperLUStat_t stat; StatInit(&stat);
    int_t info=0;
    double rpg=0.0,rcond=0.0,ferr=0.0,berr=0.0;
    mem_usage_t mem_usage{};
    dgssvx(&options,&A,impl.perm_c.data(),impl.perm_r.data(),impl.etree.data(),
           impl.equed,impl.R.data(),impl.C.data(),&impl.L,&impl.U,
           nullptr,0,&B,&X,&rpg,&rcond,&ferr,&berr,&impl.Glu,&mem_usage,&stat,&info);
    StatFree(&stat);
    Destroy_SuperMatrix_Store(&A);
    Destroy_SuperMatrix_Store(&B);
    Destroy_SuperMatrix_Store(&X);
    if(info!=0 && info!=impl.n+1){
        impl.destroy_factors();
        throw std::runtime_error("SuperLU same-pattern factor/solve failed, info="+std::to_string(info));
    }
    impl.factored=true;
    ++impl.factorization_count;
    return x;
}

SuperLUSamePatternSolver::SuperLUSamePatternSolver(const SparseMatrixCSC& initial_matrix)
    : impl_(std::make_unique<Impl>()) {
    if(initial_matrix.rows()!=initial_matrix.cols()) throw std::invalid_argument("Same-pattern factor requires square matrix");
    impl_->n=initial_matrix.rows();
    impl_->pattern_col_ptr=initial_matrix.col_ptr();
    impl_->pattern_row_ind=initial_matrix.row_ind();
    impl_->perm_r.resize(static_cast<std::size_t>(impl_->n));
    impl_->perm_c.resize(static_cast<std::size_t>(impl_->n));
    impl_->etree.resize(static_cast<std::size_t>(impl_->n));
    impl_->R.resize(static_cast<std::size_t>(impl_->n));
    impl_->C.resize(static_cast<std::size_t>(impl_->n));
    std::vector<double> zero(static_cast<std::size_t>(impl_->n),0.0);
    (void)dgssvx_factor_solve(*impl_,initial_matrix,zero,DOFACT);
}
SuperLUSamePatternSolver::~SuperLUSamePatternSolver()=default;
SuperLUSamePatternSolver::SuperLUSamePatternSolver(SuperLUSamePatternSolver&&) noexcept=default;
SuperLUSamePatternSolver& SuperLUSamePatternSolver::operator=(SuperLUSamePatternSolver&&) noexcept=default;
int SuperLUSamePatternSolver::size() const{return impl_->n;}
std::size_t SuperLUSamePatternSolver::factorizations() const{return impl_->factorization_count;}

std::vector<double> SuperLUSamePatternSolver::solve_current(const std::vector<double>& rhs) const {
    if(!impl_->factored)throw std::runtime_error("SuperLU factors invalid after failed refactorization");
    if(static_cast<int>(rhs.size())!=impl_->n) throw std::invalid_argument("SuperLU RHS dimension mismatch");
    std::vector<double> x=rhs;
    SuperMatrix B{};
    dCreate_Dense_Matrix(&B,impl_->n,1,x.data(),impl_->n,SLU_DN,SLU_D,SLU_GE);
    SuperLUStat_t stat; StatInit(&stat); int info=0;
    dgstrs(NOTRANS,&impl_->L,&impl_->U,impl_->perm_c.data(),impl_->perm_r.data(),&B,&stat,&info);
    StatFree(&stat); Destroy_SuperMatrix_Store(&B);
    if(info!=0) throw std::runtime_error("SuperLU same-pattern triangular solve failed, info="+std::to_string(info));
    return x;
}

std::vector<double> SuperLUSamePatternSolver::refactor_and_solve(const SparseMatrixCSC& matrix,
                                                                  const std::vector<double>& rhs) {
    validate_same_pattern(matrix,*impl_);
    return dgssvx_factor_solve(*impl_,matrix,rhs,impl_->factored?SamePattern:DOFACT);
}

std::vector<double> superlu_solve_once(const SparseMatrixCSC& matrix,
                                       const std::vector<double>& rhs) {
    if (matrix.rows() != matrix.cols() || static_cast<int>(rhs.size()) != matrix.rows())
        throw std::invalid_argument("SuperLU one-shot solve dimension mismatch");
    SuperMatrix A = make_csc_view(matrix);
    std::vector<double> x = rhs;
    SuperMatrix B{};
    dCreate_Dense_Matrix(&B, matrix.rows(), 1, x.data(), matrix.rows(),
                         SLU_DN, SLU_D, SLU_GE);
    std::vector<int> perm_r(static_cast<std::size_t>(matrix.rows()));
    std::vector<int> perm_c(static_cast<std::size_t>(matrix.rows()));
    SuperMatrix L{}, U{};
    superlu_options_t options;
    set_default_options(&options);
    options.ColPerm = COLAMD;
    SuperLUStat_t stat;
    StatInit(&stat);
    int info = 0;
    dgssv(&options, &A, perm_c.data(), perm_r.data(), &L, &U, &B, &stat, &info);
    StatFree(&stat);
    Destroy_SuperMatrix_Store(&A);
    Destroy_SuperMatrix_Store(&B);
    if (L.Store) Destroy_SuperNode_Matrix(&L);
    if (U.Store) Destroy_CompCol_Matrix(&U);
    if (info != 0) throw std::runtime_error("SuperLU one-shot solve failed, info=" + std::to_string(info));
    return x;
}

} // namespace quake

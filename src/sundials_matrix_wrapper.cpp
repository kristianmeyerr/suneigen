#include "sundials_matrix_wrapper.h"
#include <sundials/sundials_matrix.h> // return codes

#include <new> // bad_alloc
#include <utility>
#include <stdexcept> // invalid_argument and domain_error

namespace suneigen {

    SUNMatrixWrapper::SUNMatrixWrapper(sunindextype M, sunindextype N,
                                       sunindextype NNZ, int sparsetype)
            : matrix_(SUNSparseMatrix(M, N, NNZ, sparsetype)), id_(SUNMATRIX_SPARSE),
              sparsetype_(sparsetype) {

        if (sparsetype != CSC_MAT && sparsetype != CSR_MAT)
            throw std::invalid_argument("Invalid sparsetype. Must be CSC_MAT or "
                                        "CSR_MAT");

        if (NNZ && M && N && !matrix_)
            throw std::bad_alloc();

        finish_init();
        assert(num_nonzeros() == 0);
        assert(NNZ == capacity() || !matrix_);
        assert(M == rows() || !matrix_);
        assert(N == columns() || !matrix_);
    }

    SUNMatrixWrapper::SUNMatrixWrapper(sunindextype M, sunindextype N)
            : matrix_(SUNDenseMatrix(M, N)), id_(SUNMATRIX_DENSE) {
        if (M && N && !matrix_)
            throw std::bad_alloc();

        finish_init();
        assert(M == rows());
        assert(N == columns());
    }

    SUNMatrixWrapper::SUNMatrixWrapper(sunindextype M, sunindextype ubw,
                                       sunindextype lbw)
            : matrix_(SUNBandMatrix(M, ubw, lbw)), id_(SUNMATRIX_BAND) {
        if (M && !matrix_)
            throw std::bad_alloc();
        finish_init();
    }

    SUNMatrixWrapper::SUNMatrixWrapper(const SUNMatrixWrapper &A, realtype droptol,
                                       int sparsetype)
            : id_(SUNMATRIX_SPARSE), sparsetype_(sparsetype) {
        if (sparsetype != CSC_MAT && sparsetype != CSR_MAT)
            throw std::invalid_argument("Invalid sparsetype. Must be CSC_MAT or "
                                        "CSR_MAT");

        switch (A.matrix_id()) {
            case SUNMATRIX_DENSE:
                matrix_ = SUNSparseFromDenseMatrix(A.get(), droptol, sparsetype);
                break;
            case SUNMATRIX_BAND:
                matrix_ = SUNSparseFromBandMatrix(A.get(), droptol, sparsetype);
                break;
            default:
                throw std::invalid_argument("Invalid Matrix. Must be SUNMATRIX_DENSE or"
                                            " SUNMATRIX_BAND");
        }
        if (!matrix_)
            throw std::bad_alloc();
        finish_init();
        num_nonzeros_ = indexptrs_[num_indexptrs()];
    }

    static inline SUNMatrix_ID get_sparse_id_w_default(SUNMatrix mat) {
        if (mat)
            return SUNMatGetID(mat);
        return SUNMATRIX_CUSTOM;
    }

    static inline int get_sparse_type_w_default(SUNMatrix mat) {
        if (mat && SUNMatGetID(mat) == SUNMATRIX_SPARSE)
            return SM_SPARSETYPE_S(mat);
        return CSC_MAT;
    }

    SUNMatrixWrapper::SUNMatrixWrapper(SUNMatrix mat)
            : matrix_(mat),  id_(get_sparse_id_w_default(mat)),
              sparsetype_(get_sparse_type_w_default(mat)), ownmat(false) {
        finish_init();
    }

    SUNMatrixWrapper::~SUNMatrixWrapper() {
        if (matrix_ && ownmat)
            SUNMatDestroy(matrix_);
    }

    SUNMatrixWrapper::SUNMatrixWrapper(const SUNMatrixWrapper& other)
            : id_(get_sparse_id_w_default(other.matrix_)),
              sparsetype_(get_sparse_type_w_default(other.matrix_))  {
        if (!other.matrix_)
            return;

        matrix_ = SUNMatClone(other.matrix_);
        if (!matrix_)
            throw std::bad_alloc();

        SUNMatCopy(other.matrix_, matrix_);
        finish_init();
    }

    SUNMatrixWrapper::SUNMatrixWrapper(SUNMatrixWrapper &&other) noexcept
            : id_(get_sparse_id_w_default(other.matrix_)),
              sparsetype_(get_sparse_type_w_default(other.matrix_)) {
        std::swap(matrix_, other.matrix_);
        finish_init();
    }

    SUNMatrixWrapper &SUNMatrixWrapper::operator=(const SUNMatrixWrapper &other) {
        if(&other == this)
            return *this;
        return *this = SUNMatrixWrapper(other);
    }

    SUNMatrixWrapper &SUNMatrixWrapper::operator=(SUNMatrixWrapper &&other) noexcept {
        std::swap(matrix_, other.matrix_);
        id_ = other.id_;
        sparsetype_ = other.sparsetype_;
        finish_init();
        return *this;
    }

    void SUNMatrixWrapper::reallocate(sunindextype NNZ) {
        if (sparsetype() != CSC_MAT && sparsetype() != CSR_MAT)
            throw std::invalid_argument("Invalid sparsetype. Must be CSC_MAT or "
                                        "CSR_MAT.");

        if (int ret = SUNSparseMatrix_Reallocate(matrix_, NNZ) != SUNMAT_SUCCESS)
            throw std::runtime_error("SUNSparseMatrix_Reallocate failed with "
                                     "error code " + std::to_string(ret) + ".");

        update_ptrs();
        capacity_ = NNZ;
        assert((NNZ && columns() && rows()) || !matrix_);
    }

    void SUNMatrixWrapper::realloc() {
        if (sparsetype() != CSC_MAT && sparsetype() != CSR_MAT)
            throw std::invalid_argument("Invalid sparsetype. Must be CSC_MAT or "
                                        "CSR_MAT.");
        if (int ret = SUNSparseMatrix_Realloc(matrix_) != SUNMAT_SUCCESS)
            throw std::runtime_error("SUNSparseMatrix_Realloc failed with "
                                     "error code " + std::to_string(ret) + ".");

        update_ptrs();
        capacity_ = num_nonzeros_;
        assert(capacity() || !matrix_);
    }

    sunindextype SUNMatrixWrapper::num_indexptrs() const {
        assert(matrix_id() == SUNMATRIX_SPARSE);
        assert(!matrix_ ||
               (sparsetype() == CSC_MAT ?
                num_indexptrs_ == num_columns_ :
                num_indexptrs_ == num_rows_));
        assert(!matrix_ || num_indexptrs_ == SM_NP_S(matrix_));
        return num_indexptrs_;
    }

    sunindextype SUNMatrixWrapper::capacity() const {
        assert(matrix_id() == SUNMATRIX_SPARSE);
        assert(!matrix_ || capacity_ == SM_NNZ_S(matrix_));
        return capacity_;
    }

    sunindextype SUNMatrixWrapper::num_nonzeros() const {
        assert(matrix_id() == SUNMATRIX_SPARSE);
        assert(!matrix_ ||
               num_nonzeros_ == SM_INDEXPTRS_S(matrix_)[SM_NP_S(matrix_)]);
        return num_nonzeros_;
    }

    const realtype *SUNMatrixWrapper::data() const {
        return data_;
    }

    realtype *SUNMatrixWrapper::data() {
        return data_;
    }

    int SUNMatrixWrapper::sparsetype() const {
        assert(matrix_);
        assert(matrix_id() == SUNMATRIX_SPARSE);
        return sparsetype_;
    }


#ifndef NDEBUG
    static inline void check_csc(const SUNMatrixWrapper *mat) {
        assert(mat->matrix_id() == SUNMATRIX_SPARSE);
        assert(mat->sparsetype() == CSC_MAT);
    }
#else
    // avoid "unused parameter" warning
static inline void check_csc(const SUNMatrixWrapper */*mat*/) {}
#endif

    void SUNMatrixWrapper::to_dense(SUNMatrixWrapper &D) const {
        if (!matrix_ || !D.matrix_)
            return;
        check_csc(this);
        assert(rows() == D.rows());
        assert(columns() == D.columns());

        D.zero();
        if (!num_nonzeros())
            return;

        sunindextype icol;
        sunindextype idx;
        for (icol = 0; icol < columns(); ++icol)
            for (idx = get_indexptr(icol); idx < get_indexptr(icol+1); ++idx) {
                D.set_data(get_indexval(idx), icol, get_data(idx));
            }
    }

    void SUNMatrixWrapper::to_diag(N_Vector v) const {
        if (!matrix_ || !v)
            return;
        check_csc(this);
        assert(rows() == columns());
        assert(rows() == NV_LENGTH_S(v));

        N_VConst(0.0, v);
        if (!num_nonzeros())
            return;

        sunindextype icol;
        sunindextype idx;
        for (icol = 0; icol < columns(); ++icol)
            for (idx = get_indexptr(icol); idx < get_indexptr(icol+1); ++idx)
                if (get_indexval(idx) == icol)
                    NV_Ith_S(v, icol) = get_data(idx);
    }


    void SUNMatrixWrapper::zero()
    {
        if (!matrix_)
            return;
        if(int res = SUNMatZero(matrix_))
            throw std::runtime_error("SUNMatrixWrapper::zero() failed with "
                                     + std::to_string(res) + ".");
    }

    void SUNMatrixWrapper::finish_init() {
        update_ptrs();
        update_size();
    }

    void SUNMatrixWrapper::update_ptrs() {
        if (!matrix_) {
            data_ = nullptr;
            indexptrs_ = nullptr;
            indexvals_ = nullptr;
            return;
        }

        switch (matrix_id()) {
            case SUNMATRIX_DENSE:
                data_ = SM_DATA_D(matrix_);
                break;
            case SUNMATRIX_SPARSE:
                data_ = SM_DATA_S(matrix_);
                indexptrs_ = SM_INDEXPTRS_S(matrix_);
                indexvals_ = SM_INDEXVALS_S(matrix_);
                break;
            default:
                throw std::domain_error("Not Implemented.");
        }
    }

    void SUNMatrixWrapper::update_size() {
        num_indexptrs_ = 0;
        if (!matrix_) {
            num_rows_ = 0;
            num_columns_ = 0;
            return;
        }

        switch (matrix_id()) {
            case SUNMATRIX_DENSE:
                num_rows_ = SM_ROWS_D(matrix_);
                num_columns_ = SM_COLUMNS_D(matrix_);
                capacity_ = num_rows_ * num_columns_;
                break;
            case SUNMATRIX_SPARSE:
                num_rows_ = SM_ROWS_S(matrix_);
                num_columns_ = SM_COLUMNS_S(matrix_);
                capacity_ = SM_NNZ_S(matrix_);
                num_indexptrs_ = SM_NP_S(matrix_);
                break;
            default:
                throw std::domain_error("Not Implemented.");
        }
    }

    void SUNMatrixWrapper::refresh() {
        update_ptrs();
        update_size();
        if (matrix_id() == SUNMATRIX_SPARSE)
            num_nonzeros_ = SM_INDEXPTRS_S(matrix_)[SM_NP_S(matrix_)];
    }

    SUNMatrix SUNMatrixWrapper::get() const { return matrix_; }

} // namespace suneigen

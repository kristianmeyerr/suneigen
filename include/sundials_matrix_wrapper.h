#ifndef SUNEIGEN_SUNDIALS_MATRIX_WRAPPER_H
#define SUNEIGEN_SUNDIALS_MATRIX_WRAPPER_H

#include <gsl/gsl-lite.hpp>

#include <sunmatrix/sunmatrix_band.h>   // SUNMatrix_Band
#include <sunmatrix/sunmatrix_dense.h>  // SUNMatrix_Dense
#include <sunmatrix/sunmatrix_sparse.h>  // SUNMatrix_Dense

#include <gsl/gsl-lite.hpp>

#include <vector>

#include <assert.h>

#include "vector.h"

namespace suneigen {

    /**
     * @brief A RAII wrapper for SUNMatrix structs.
     *
     * This can create dense, sparse, or banded matrices using the respective
     * constructor.
     */
    class SUNMatrixWrapper {
    public:
        SUNMatrixWrapper() = default;

        /**
         * @brief Create sparse matrix. See SUNSparseMatrix in sunmatrix_sparse.h
         * @param M Number of rows
         * @param N Number of columns
         * @param NNZ Number of nonzeros
         * @param sparsetype Sparse type
         */
        SUNMatrixWrapper(sunindextype M, sunindextype N, sunindextype NNZ,
                         int sparsetype);

        /**
         * @brief Create dense matrix. See SUNDenseMatrix in sunmatrix_dense.h
         * @param M Number of rows
         * @param N Number of columns
         */
        SUNMatrixWrapper(sunindextype M, sunindextype N);

        /**
         * @brief Create banded matrix. See SUNBandMatrix in sunmatrix_band.h
         * @param M Number of rows and columns
         * @param ubw Upper bandwidth
         * @param lbw Lower bandwidth
         */
        SUNMatrixWrapper(sunindextype M, sunindextype ubw, sunindextype lbw);

        /**
         * @brief Create sparse matrix from dense or banded matrix. See
         * SUNSparseFromDenseMatrix and SUNSparseFromBandMatrix in
         * sunmatrix_sparse.h
         * @param A Wrapper for dense matrix
         * @param droptol tolerance for dropping entries
         * @param sparsetype Sparse type
         */
        SUNMatrixWrapper(const SUNMatrixWrapper &A, realtype droptol,
                         int sparsetype);

        /**
         * @brief Wrap existing SUNMatrix
         * @param mat input matrix
         */
        explicit SUNMatrixWrapper(SUNMatrix mat);

        ~SUNMatrixWrapper();

        /**
         * @brief Copy constructor
         * @param other Input matrix
         */
        SUNMatrixWrapper(const SUNMatrixWrapper &other);

        /**
         * @brief Move constructor
         * @param other Input matrix
         */
        SUNMatrixWrapper(SUNMatrixWrapper &&other) noexcept;

        /**
         * @brief Copy assignment
         * @param other Input matrix
         * @return matrix
         */
        SUNMatrixWrapper &operator=(const SUNMatrixWrapper &other);

        /**
         * @brief Move assignment
         * @param other Input matrix
         * @return matrix
         */
        SUNMatrixWrapper &operator=(SUNMatrixWrapper &&other) noexcept;

        /**
         * @brief Reallocate space for sparse matrix according to specified nnz
         * @param nnz new number of nonzero entries
         */
        void reallocate(sunindextype nnz);

        /**
         * @brief Reallocate space for sparse matrix to used space according to last entry in indexptrs
         */
        void realloc();

        /**
         * @brief Get the wrapped SUNMatrix
         * @return raw SunMatrix object
         * @note Even though the returned matrix_ pointer is const qualified, matrix_->content will not be const.
         * This is a shortcoming in the underlying C library, which we cannot address and it is not intended that
         * any of those values are modified externally. If matrix_->content is manipulated,
         * cpp:meth:SUNMatrixWrapper:`refresh` needs to be called.
         */
        [[nodiscard]] SUNMatrix get() const;

        /**
         * @brief Get the number of rows
         * @return number of rows
         */
        [[nodiscard]] sunindextype rows() const {
            assert(!matrix_ ||
                   (matrix_id() == SUNMATRIX_SPARSE ?
                    num_rows_ == SM_ROWS_S(matrix_) :
                    num_rows_ == SM_ROWS_D(matrix_)));
            return num_rows_;
        }

        /**
         * @brief Get the number of columns
         * @return number of columns
         */
        [[nodiscard]] sunindextype columns() const {
            assert(!matrix_ ||
                   (matrix_id() == SUNMATRIX_SPARSE ?
                    num_columns_ == SM_COLUMNS_S(matrix_) :
                    num_columns_ == SM_COLUMNS_D(matrix_)));
            return num_columns_;
        }

        /**
         * @brief Get the number of specified non-zero elements (sparse matrices only)
         * @note value will be 0 before indexptrs are set.
         * @return number of nonzero entries
         */
        [[nodiscard]] sunindextype num_nonzeros() const;

        /**
         * @brief Get the number of indexptrs that can be specified (sparse matrices only)
         * @return number of indexptrs
         */
        [[nodiscard]] sunindextype num_indexptrs() const;

        /**
         * @brief Get the number of allocated data elements
         * @return number of allocated entries
         */
        [[nodiscard]] sunindextype capacity() const;

        /**
         * @brief Get  raw data of a sparse matrix
         * @return pointer to first data entry
         */
        realtype *data();

        /**
         * @brief Get const raw data of a sparse matrix
         * @return pointer to first data entry
         */
        [[nodiscard]] const realtype *data() const;

        /**
         * @brief Get data of a sparse matrix
         * @param idx data index
         * @return idx-th data entry
         */
        [[nodiscard]] realtype get_data(sunindextype idx) const{
            assert(matrix_);
            assert(matrix_id() == SUNMATRIX_SPARSE);
            assert(idx < capacity());
            assert(SM_DATA_S(matrix_) == data_);
            return data_[idx];
        }

        /**
         * @brief Get data entry for a dense matrix
         * @param irow row
         * @param icol col
         * @return A(irow,icol)
         */
        [[nodiscard]] realtype get_data(sunindextype irow, sunindextype icol) const{
            assert(matrix_);
            assert(matrix_id() == SUNMATRIX_DENSE);
            assert(irow < rows());
            assert(icol < columns());
            return SM_ELEMENT_D(matrix_, irow, icol);
        }

        /**
         * @brief Set data entry for a sparse matrix
         * @param idx data index
         * @param data data for idx-th entry
         */
        void set_data(sunindextype idx, realtype data) {
            assert(matrix_);
            assert(matrix_id() == SUNMATRIX_SPARSE);
            assert(idx < capacity());
            assert(SM_DATA_S(matrix_) == data_);
            data_[idx] = data;
        }

        /**
         * @brief Set data entry for a dense matrix
         * @param irow row
         * @param icol col
         * @param data data for idx-th entry
         */
        void set_data(sunindextype irow, sunindextype icol, realtype data) {
            assert(matrix_);
            assert(matrix_id() == SUNMATRIX_DENSE);
            assert(irow < rows());
            assert(icol < columns());
            SM_ELEMENT_D(matrix_, irow, icol) = data;
        }

        /**
         * @brief Get the index value of a sparse matrix
         * @param idx data index
         * @return row (CSC) or column (CSR) for idx-th data entry
         */
        sunindextype get_indexval(sunindextype idx) const {
            assert(matrix_);
            assert(matrix_id() == SUNMATRIX_SPARSE);
            assert(idx < capacity());
            assert(indexvals_ == SM_INDEXVALS_S(matrix_));
            return indexvals_[idx];
        }

        /**
         * @brief Set the index value of a sparse matrix
         * @param idx data index
         * @param val row (CSC) or column (CSR) for idx-th data entry
         */
        void set_indexval(sunindextype idx, sunindextype val) {
            assert(matrix_);
            assert(matrix_id() == SUNMATRIX_SPARSE);
            assert(idx < capacity());
            assert(indexvals_ == SM_INDEXVALS_S(matrix_));
            indexvals_[idx] = val;
        }

        /**
         * @brief Set the index values of a sparse matrix
         * @param vals rows (CSC) or columns (CSR) for data entries
         */
        void set_indexvals(const gsl::span<const sunindextype> vals) {
            assert(matrix_);
            assert(matrix_id() == SUNMATRIX_SPARSE);
            assert(static_cast<sunindextype>(vals.size()) == capacity());
            assert(indexvals_ == SM_INDEXVALS_S(matrix_));
            std::copy_n(vals.begin(), capacity(), indexvals_);
        }

        /**
         * @brief Get the index pointer of a sparse matrix
         * @param ptr_idx pointer index
         * @return index where the ptr_idx-th column (CSC) or row (CSR) starts
         */
        sunindextype get_indexptr(sunindextype ptr_idx) const {
            assert(matrix_);
            assert(matrix_id() == SUNMATRIX_SPARSE);
            assert(ptr_idx <= num_indexptrs());
            assert(indexptrs_ == SM_INDEXPTRS_S(matrix_));
            return indexptrs_[ptr_idx];
        }

        /**
         * @brief Set the index pointer of a sparse matrix
         * @param ptr_idx pointer index
         * @param ptr data-index where the ptr_idx-th column (CSC) or row (CSR) starts
         */
        void set_indexptr(sunindextype ptr_idx, sunindextype ptr) {
            assert(matrix_);
            assert(matrix_id() == SUNMATRIX_SPARSE);
            assert(ptr_idx <= num_indexptrs());
            assert(ptr <= capacity());
            assert(indexptrs_ == SM_INDEXPTRS_S(matrix_));
            indexptrs_[ptr_idx] = ptr;
            if (ptr_idx == num_indexptrs())
                num_nonzeros_ = ptr;
        }

        /**
         * @brief Set the index pointers of a sparse matrix
         * @param ptrs starting data-indices where the columns (CSC) or rows (CSR) start
         */
        void set_indexptrs(const gsl::span<const sunindextype> ptrs) {
            assert(matrix_);
            assert(matrix_id() == SUNMATRIX_SPARSE);
            assert(static_cast<sunindextype>(ptrs.size()) == num_indexptrs() + 1);
            assert(indexptrs_ == SM_INDEXPTRS_S(matrix_));
            std::copy_n(ptrs.begin(), num_indexptrs() + 1, indexptrs_);
            num_nonzeros_ = indexptrs_[num_indexptrs()];
        }

        /**
         * @brief Get the type of sparse matrix
         * @return matrix type
         */
        [[nodiscard]] int sparsetype() const;

        /**
         * @brief Writes a sparse matrix A to a dense matrix D.
         *
         * @param D dense output matrix
         */
        void to_dense(SUNMatrixWrapper &D) const;

        /**
         * @brief Writes the diagonal of sparse matrix A to a dense vector v.
         *
         * @param v dense outut vector
         */
        void to_diag(N_Vector v) const;

        /**
         * @brief Set to 0.0, for sparse matrices also resets indexptr/indexvals
         */
        void zero();

        /**
         * @brief Get matrix id
         * @return SUNMatrix_ID
         */
        [[nodiscard]] SUNMatrix_ID matrix_id() const {return id_;}

        /**
         * @brief Update internal cache, needs to be called after external manipulation of matrix_->content
         */
        void refresh();

    private:

        /**
         * @brief SUNMatrix to which all methods are applied
         */
        SUNMatrix matrix_ {nullptr};

        /**
         * @brief cache for SUNMatrixGetId(matrix_)
         */
        SUNMatrix_ID id_ {SUNMATRIX_CUSTOM};

        /**
         * @brief cache for SUNMatrixGetId(matrix_)
         */
        int sparsetype_ {CSC_MAT};

        /**
         * @brief cache for SM_INDEXPTRS_S(matrix_)[SM_NP_S(matrix_)]
         */
        sunindextype num_nonzeros_ {0};
        /**
         * @brief cache for SM_NNZ_S(matrix_)
         */
        sunindextype capacity_ {0};

        /**
         * @brief cache for SM_DATA_S(matrix_)
         */
        realtype *data_ {nullptr};
        /**
         * @brief cache for SM_INDEXPTRS_S(matrix_)
         */
        sunindextype *indexptrs_ {nullptr};
        /**
         * @brief cache for SM_INDEXVALS_S(matrix_)
         */
        sunindextype *indexvals_ {nullptr};

        /**
         * @brief cache for SM_ROWS_X(matrix_)
         */
        sunindextype num_rows_ {0};
        /**
         * @brief cache for SM_COLUMS_X(matrix_)
         */
        sunindextype num_columns_ {0};

        /**
         * @brief cache for SM_NP_S(matrix_)
         */
        sunindextype num_indexptrs_ {0};

        /**
         * @brief call update_ptrs & update_size
         */
        void finish_init();
        /**
         * @brief update data_, indexptrs_, indexvals_ if applicable
         */
        void update_ptrs();
        /**
         * @brief update num_rows_, num_columns_, num_indexptrs if applicable
         */
        void update_size();
        /**
         * @brief indicator whether this wrapper allocated matrix_ and is responsible for deallocation
         */
        bool ownmat = true;


    };

}

namespace gsl {
    /**
     * @brief Create span from SUNMatrix
     * @param m SUNMatrix
     * @return Created span
     */
    inline span<realtype> make_span(SUNMatrix m)
    {
        switch (SUNMatGetID(m)) {
            case SUNMATRIX_DENSE:
                return span<realtype>(SM_DATA_D(m), static_cast<std::size_t>(SM_LDATA_D(m)));
            case SUNMATRIX_SPARSE:
                return span<realtype>(SM_DATA_S(m), static_cast<std::size_t>(SM_NNZ_S(m)));
            default:
                throw suneigen::SunException("Unimplemented SUNMatrix type for make_span");
        }
    }
} // namespace gsl

#endif //SUNEIGEN_SUNDIALS_MATRIX_WRAPPER_H

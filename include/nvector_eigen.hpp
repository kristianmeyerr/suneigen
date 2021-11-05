//
// Created by Kristian Meyer on 2019-09-01.
//

#ifndef SUNEIGENMATRIX_NVECTOREIGEN_HPP
#define SUNEIGENMATRIX_NVECTOREIGEN_HPP

#include <sundials/sundials_nvector.h>
#include "Eigen/Dense"


struct _N_VectorContent_Eigen {
    sunindextype length;        // vector length
    bool isOwnData;             // data ownership flag
    Eigen::VectorXd* vecPtr;    // data array
};

// Forward declaration
typedef struct _N_VectorContent_Eigen* N_VectorContent_Eigen;


// --------------------------------------- //
// Declare functions used by NVector_Eigen //
// --------------------------------------- //

/**
 * Allocates memory for a new NVector and attaches the operations and content structure with empty data
 *
 * @param [in] length Length of the NVector
 * @return Returns an empty NVector
 */
SUNDIALS_EXPORT N_Vector N_VNewEmpty_Eigen(sunindextype length);

/**
 * Creates data in an empty VNVector using the Eigen copy constructor
 *
 * @param [in] length Length of the NVector
 * @return Returns an Nvector with data allocated
 */
SUNDIALS_EXPORT N_Vector N_VNew_Eigen(sunindextype length);

/**
 * unction to create an Eigen N_Vector with user data component
 *
 * @param [in] length Length of the NVector
 * @param [in] v_data The user data
 * @return Returns an Nvector with data allocated
 */
SUNDIALS_EXPORT N_Vector N_VMake_Eigen(sunindextype length, realtype *v_data);

/**
 * Get the vector ID of the NVector
 * @param [in] w Input vector
 * @return Returns the vector type identifier for the vector w
 */
SUNDIALS_EXPORT N_Vector_ID N_VGetVectorID_Eigen(N_Vector w);

/**
 * Returns a pointer to a double array from the N_Vector v.
 * Note that this assumes that the internal data in N Vector is a contiguous array of double.
 * This routine is only used in the solver-specific interfaces to the dense and banded (serial) linear solvers,
 * the sparse linear solvers (serial and threaded), and in the interfaces to the banded (serial)
 * and band-block- diagonal (parallel) preconditioner modules provided with sundials.
 * @param [in] v Input vector
 * @return Returns a pointer to a double array from the N_Vector v
 */
SUNDIALS_EXPORT realtype* N_VGetArrayPointer_Eigen(N_Vector v);

/**
 * Creates a new N_Vector of the same type as an existing vector w
 * and sets the ops field. It does not allocate storage for data.
 *
 * @param [in] w Input Vector
 * @return Returns a new N_Vector similar to w with no data
 */
SUNDIALS_EXPORT N_Vector N_VCloneEmpty_Eigen(N_Vector w);

/**
 * Creates a new N Vector of the same type as an existing vector w
 * and sets the ops field. It does not copy the vector,
 * but rather allocates storage for the new vector.
 *
 * @param [in] w Input vector
 * @return Returns a cloned vector with allocated memory
 */
SUNDIALS_EXPORT N_Vector N_VClone_Eigen(N_Vector w);

/**
 * Destroys the N_Vector v and frees memory allocated for its internal data.
 *
 * @param [in] v Input vector
 */
SUNDIALS_EXPORT void N_VDestroy_Eigen(N_Vector v);

/**
 * Returns storage requirements for one N Vector.
 * lrw contains the number of double words and liw contains the number
 * of integer words. This function is advisory only,
 * for use in determining a user’s total space requirements;
 * it could be a dummy function in a user-supplied nvector
 * module if that information is not of interest.
 *
 * @param [in] v Input vector
 * @param [out] lrw Number of double (double) words
 * @param [out] liw Number of integer words
 */
SUNDIALS_EXPORT void N_VSpace_Eigen(N_Vector v, sunindextype *lrw, sunindextype *liw);

/**
 * Performs the operation z = ax + by, where a and b are double scalars and x and y
 * are of type N Vector: @f$ zi = a*x_i + b*y_i, i=0,...,n−1 @f$.
 *
 * @param [in] a scalar
 * @param [in] x Input vector
 * @param [in] b scalar
 * @param [in] y Input vector
 * @param [out] z Output vector
 */
SUNDIALS_EXPORT void N_VLinearSum_Eigen(realtype a, N_Vector x, realtype b, N_Vector y, N_Vector z);

/**
 * Sets all components of the N_Vector z to double c: @f$ z_i = c, i = 0,...,n− 1@f$.
 *
 * @param [in] c double scalar
 * @param [out] z Output vector
 */
SUNDIALS_EXPORT void N_VConst_Eigen(realtype c, N_Vector z);

/**
 * Sets the N_Vector z to be the component-wise product of the N_Vector inputs x and y:
 * @f$ z_i = x_i*y_i, i=0,...,n−1.
 *
 * @param [in] x Input vector
 * @param [in] y Input vector
 * @param [out] z Output vector
 */
SUNDIALS_EXPORT void N_VProd_Eigen(N_Vector x, N_Vector y, N_Vector z);

/**
 * Sets the N Vector z to be the component-wise ratio of the N_Vector inputs x and y:
 * @f$ z_i = x_i/y_i, i=0,...,n−1@f$. The y_i may not be tested for 0 values.
 * It should only be called with a y that is guaranteed to have all non-zero components.
 *
 * @param [in] x Input vector
 * @param [in] y Input vector
 * @param [out] z Output vector
 */
SUNDIALS_EXPORT void N_VDiv_Eigen(N_Vector x, N_Vector y, N_Vector z);

/**
 * Sets the N Vector z to be the component-wise ratio of the N Vector inputs x and y:
 * @f$ z_i = x_i/y_i, i=0,...,n−1@f$.
 * The y_i may not be tested for 0 values.
 * It should only be called with a y that is guaranteed to have all nonzero components.
 *
 * @param [in] c Realtype scalar
 * @param [in] x Input vector
 * @param [out] z Output vector
 */
SUNDIALS_EXPORT void N_VScale_Eigen(realtype c, N_Vector x, N_Vector z);

/**
 * Sets the components of the N_Vector z to be the absolute values of the components
 * of the N_Vector x: y_i = |x_i|, i = 0,...,n−1.
 *
 * @param [in] x Input vector
 * @param [out] z Output vector
 */
SUNDIALS_EXPORT void N_VAbs_Eigen(N_Vector x, N_Vector z);

/**
 * Sets the components of the N_Vector z to be the inverses of the components of
 * the N_Vector x: z_i = 1.0/x_i, i = 0,...,n − 1.
 * This routine may not check for division by 0.
 * It should be called only with an x which is guaranteed to have all nonzero components.
 *
 * @param [in] x Input vector
 * @param [out] z Output vector
 */
SUNDIALS_EXPORT void N_VInv_Eigen(N_Vector x, N_Vector z);

/**
 * Adds the double scalar b to all components of x and returns the result
 * in the N_Vector z: z_i = x_i + b, i = 0, ... , n−1.
 *
 * @param [in] x Input vector
 * @param [in] b Realtype scalar
 * @param [in] z Output vector
 */
SUNDIALS_EXPORT void N_VAddConst_Eigen(N_Vector x, realtype b, N_Vector z);

/**
 * Returns the value of the ordinary dot product of x and y: @f$ d = 􏰌\sum_{i=0}^{n−1} x_i*y_i @f$.
 *
 * @param [in] x Input vector
 * @param [in] y Input vector
 * @return Returns the dot product of x and y
 */
SUNDIALS_EXPORT realtype N_VDotProd_Eigen(N_Vector x, N_Vector y);

/**
 * Returns the maximum norm of the N_Vector x: @f$ m = max_i |x_i| @f$.
 *
 * @param [in] x Input vector
 * @return Returns the maximum norm
 */
SUNDIALS_EXPORT realtype N_VMaxNorm_Eigen(N_Vector x);

/**
 * Compute the weighted root-mean-square norm of the N_Vector x with double weight vector w:
 * @f$ m = \sqrt{ \left( \sum_{i=0}^{n-1} (x_iw_i)^2 \right)/n } @f$
 *
 * @param [in] x Input vector
 * @param [in] w Input weight vector
 * @return Returns the weighted root-mean-square norm of the N_Vector x with􏰆 double weight vector w
 */
SUNDIALS_EXPORT realtype N_VWrmsNorm_Eigen(N_Vector x, N_Vector w);

/**
 * Computes the weighted root mean square norm of the N Vector x with double weight vector
 * w built using only the elements of x corresponding to positive elements of the N Vector id:
 * @f$ m = \sqrt{ \left( \sum_{i=0}^{n-1} (x_iw_ih(id_i))^2 \right)/n }, where h_i = 1 if x_i>0, otherwise h_i = 0. @f$
 *
 * @param [in] x Input vector
 * @param [in] w Input weight vector
 * @param [in] id Input ID vector corresponding to positive elements that should be included.
 * @return Returns the weighted root mean square norm of the N Vector x with double weight vector
 *      w built using only the elements of x corresponding to positive elements of the N Vector
 */
SUNDIALS_EXPORT realtype N_VWrmsNormMask_Eigen(N_Vector x, N_Vector w, N_Vector id);

/**
 * Returns the smallest element of the N_Vector x: m = min_i x_i.
 * @param [in] x Input vector
 * @return Returns the smallest element of x
 */
SUNDIALS_EXPORT realtype N_VMin_Eigen(N_Vector x);

/**
 * Returns the weighted Euclidean l2 norm of the N_Vector x with double weight vector w:
 * @f$ \sqrt{ \sum_{i=0}^{n-1} (x_i w_i)^2 } @f$
 *
 * @param [in] x Input vector
 * @param [in] w Input weight vector
 * @return Returns the weighted Euclidean l2 norm of x
 */
SUNDIALS_EXPORT realtype N_VWL2Norm_Eigen(N_Vector x, N_Vector w);

/**
 * Returns the l1 norm of the N Vector x: @f$ m = 􏰌\sum_{i=0}^{n-1} |x_i| @f$.
 *
 * @param [in] x Input vector
 * @return returns the l1 norm of the N Vector x
 */
SUNDIALS_EXPORT realtype N_VL1Norm_Eigen(N_Vector x);

/**
 * Compares the components of the N_Vector x to the double scalar c
 * and returns an NVector z such that:
 * @f$ z_i =1.0 if |x_i| ≥ cand z_i = 0.0 otherwise.
 *
 * @param [in] c double scalar
 * @param [in] x Input vector
 * @param [out] z Output vector.
 */
SUNDIALS_EXPORT void N_VCompare_Eigen(realtype c, N_Vector x, N_Vector z);

/**
 * Sets the components of the N Vector z to be the inverses of the components
 * of the N Vector x, with prior testing for zero values:
 * @f$ zi = 1.0/x_i, i = 0, . . . , n − 1 @f$.
 * This routine returns a boolean assigned to true if all components
 * of x are nonzero (successful inversion) and returns false otherwise.
 *
 * @param [in] x Input vector
 * @param [out] z Output vector
 * @return true if successfull inversion, otherwise false
 */
SUNDIALS_EXPORT booleantype N_VInvTest_Eigen(N_Vector x, N_Vector z);

/**
 * Performs the following constraint tests: xi > 0 if ci = 2,
 * xi ≥ 0 if c_i = 1, x_i ≤ 0 if c_i =−1, x_i < 0 if c_i = −2. There is no constraint on x_i if c_i = 0.
 * This routine returns a boolean assigned to false if any element failed the constraint test
 * and assigned to SUNTRUE if all passed.
 * It also sets a mask vector m, with elements equal to 1.0 where the constraint test failed,
 * and 0.0 where the test passed. This routine is used only for constraint checking.
 *
 * @param [in] c Input vector
 * @param [in] x Input vector
 * @param [in] m Input vector
 * @return Returns SUNTrue if all constrained test passed, otherwise SUNFalse
 */
SUNDIALS_EXPORT booleantype N_VConstrMask_Eigen(N_Vector c, N_Vector x, N_Vector m);

/**
 * This routine returns the minimum of the quotients obtained by
 * term-wise dividing num_i by denom_i. A zero element in denom will be skipped.
 * If no such quotients are found, then the large value BIG REAL
 * (defined in the header file sundials types.h) is returned.
 *
 * @param [in] num Input vector
 * @param [in] denom Input vector
 * @return The minimum quotient
 */
SUNDIALS_EXPORT realtype N_VMinQuotient_Eigen(N_Vector num, N_Vector denom);


//
// Imitate the macro functions used in SUNDIALS
//

SUNDIALS_EXPORT sunindextype N_VGetLength_Eigen(N_Vector v);

/**
 * Get the size of the N_Vector
 * @param [in] v Input vector
 * @return The length of v
 */
sunindextype NV_LENGTH_S(N_Vector v);

/**
 * Get a dereferenced VectorXd from an NVector
 *
 * @param [in] v Input vector
 * @return Returns a dereferenced VectorXd from an NVector
 */
SUNDIALS_EXPORT Eigen::VectorXd& NV_DATA_S(N_Vector v);

/**
 * Flag that returns true if the vector v owns its own data
 *
 * @param [in] v Input vector
 * @return Returns true if v owns the data, otherwise false
 */
SUNDIALS_EXPORT bool& NV_OWN_DATA_S(N_Vector v);


#endif //SUNEIGENMATRIX_NVECTOREIGEN_HPP

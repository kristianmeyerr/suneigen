
#include <iostream>

#include "Eigen/Dense"
#include "nvector_eigen.hpp"
#include <sundials/sundials_math.h>

N_Vector_ID N_VGetVectorID_Eigen(N_Vector w)
{
    (void)w;  // Use w for nothing
    return SUNDIALS_NVEC_CUSTOM;
}

N_Vector N_VNewEmpty_Eigen(sunindextype length)
{

    /* Create an empty vector object */
    N_Vector v = nullptr;
    v = N_VNewEmpty();
    if (v == nullptr) return(nullptr);

    /* Create content */
    N_VectorContent_Eigen content = nullptr;
    content = static_cast<N_VectorContent_Eigen>(malloc(sizeof *content));
    if (content == nullptr) { N_VDestroy(v); return(nullptr); }

    // Attach operations
    N_Vector_Ops options = nullptr;
    options = static_cast<N_Vector_Ops>(malloc(sizeof(struct _generic_N_Vector_Ops)));
    if (options == nullptr) {
        free(v);
        return(nullptr);
    }

    /* Attach operations */

    /* constructors, destructors, and utility operations */
    v->ops->nvgetvectorid     = N_VGetVectorID_Eigen;
    v->ops->nvclone           = N_VClone_Eigen;
    v->ops->nvcloneempty      = N_VCloneEmpty_Eigen;
    v->ops->nvdestroy         = N_VDestroy_Eigen;
    v->ops->nvspace           = N_VSpace_Eigen;
    v->ops->nvgetarraypointer = N_VGetArrayPointer_Eigen;
    v->ops->nvsetarraypointer = nullptr;
    v->ops->nvgetlength       = N_VGetLength_Eigen;

    /* standard vector operations */
    v->ops->nvlinearsum    = N_VLinearSum_Eigen;
    v->ops->nvconst        = N_VConst_Eigen;
    v->ops->nvprod         = N_VProd_Eigen;
    v->ops->nvdiv          = N_VDiv_Eigen;
    v->ops->nvscale        = N_VScale_Eigen;
    v->ops->nvabs          = N_VAbs_Eigen;
    v->ops->nvinv          = N_VInv_Eigen;
    v->ops->nvaddconst     = N_VAddConst_Eigen;
    v->ops->nvdotprod      = N_VDotProd_Eigen;
    v->ops->nvmaxnorm      = N_VMaxNorm_Eigen;
    v->ops->nvwrmsnormmask = N_VWrmsNormMask_Eigen;
    v->ops->nvwrmsnorm     = N_VWrmsNorm_Eigen;
    v->ops->nvmin          = N_VMin_Eigen;
    v->ops->nvwl2norm      = N_VWL2Norm_Eigen;
    v->ops->nvl1norm       = N_VL1Norm_Eigen;
    v->ops->nvcompare      = N_VCompare_Eigen;
    v->ops->nvinvtest      = N_VInvTest_Eigen;
    v->ops->nvconstrmask   = N_VConstrMask_Eigen;
    v->ops->nvminquotient  = N_VMinQuotient_Eigen;

    /* Attach content */
    v->content = content;

    /* Initialize content */
    content->length = length;
    content->isOwnData = SUNFALSE;
    content->vecPtr = nullptr;

    return v;
}

N_Vector N_VNew_Eigen(sunindextype length)
{

    N_Vector v = nullptr;
    v = N_VNewEmpty_Eigen(length);
    if (v == nullptr) return(nullptr);

    // Create data
    if (length > 0)
    {
        // Attach data
        NV_OWN_DATA_S(v) = SUNTRUE;
        reinterpret_cast<N_VectorContent_Eigen>(v->content)->vecPtr = new Eigen::VectorXd(length);
    }

    return v;
}

N_Vector N_VMake_Eigen(sunindextype length, realtype* v_data)
{
    N_Vector v = nullptr;

    v = N_VNewEmpty_Eigen(length);
    if (v == nullptr) return(nullptr);

    if (length > 0) {
        /* Attach data */
        NV_OWN_DATA_S(v) = SUNFALSE;
    }
    (void) v_data;
    return(v);
}

void N_VSpace_Eigen(N_Vector v, sunindextype* lrw, sunindextype* liw)
{
    *lrw = NV_LENGTH_S(v);
    *liw = 1;
}

void N_VLinearSum_Eigen(realtype a, N_Vector x, realtype b, N_Vector y, N_Vector z)
{
    NV_DATA_S(z).array() = a * NV_DATA_S(x).array() + b * NV_DATA_S(y).array();
}

void N_VConst_Eigen(realtype c, N_Vector z)
{
    NV_DATA_S(z).array().setConstant(c);
}

void N_VProd_Eigen(N_Vector x, N_Vector y, N_Vector z)
{
    NV_DATA_S(z) = NV_DATA_S(x).array() * NV_DATA_S(y).array();
}

void N_VDiv_Eigen(N_Vector x, N_Vector y, N_Vector z)
{
    NV_DATA_S(z) = NV_DATA_S(x).array()/NV_DATA_S(y).array();
}

void N_VScale_Eigen(double c, N_Vector x, N_Vector z)
{
    NV_DATA_S(z) = c * NV_DATA_S(x);
}

void N_VAbs_Eigen(N_Vector x, N_Vector z)
{
    NV_DATA_S(z) = NV_DATA_S(x).array().abs();
}

void N_VInv_Eigen(N_Vector x, N_Vector z)
{
    NV_DATA_S(z) = NV_DATA_S(x).array().inverse();
}

void N_VAddConst_Eigen(N_Vector x, realtype b, N_Vector z)
{
    NV_DATA_S(z) = NV_DATA_S(x) + Eigen::VectorXd::Constant(NV_LENGTH_S(x), b);
}

realtype N_VDotProd_Eigen(N_Vector x, N_Vector y)
{
    return NV_DATA_S(x).dot(NV_DATA_S(y));
}

realtype N_VMaxNorm_Eigen(N_Vector x)
{
    return NV_DATA_S(x).lpNorm<Eigen::Infinity>();
}

realtype N_VWrmsNorm_Eigen(N_Vector x, N_Vector w)
{
    return std::sqrt(((((NV_DATA_S(x).array()*NV_DATA_S(w).array()).pow(Eigen::ArrayXd::Constant(NV_DATA_S(x).size(), 2.0))).sum())/NV_DATA_S(x).size()));
}

realtype N_VWrmsNormMask_Eigen(N_Vector x, N_Vector w, N_Vector id)
{
    double prodi = 0;
    double sum = 0;
    for (int i = 0; i < NV_DATA_S(x).size(); i++) {
        if (NV_DATA_S(id)[i] > 0.0) {
            prodi = NV_DATA_S(x)[i] * NV_DATA_S(w)[i];
            sum += std::pow(prodi, 2.0);
        }
    }
    return std::sqrt(sum / NV_DATA_S(x).size());
}

realtype N_VWL2Norm_Eigen(N_Vector x, N_Vector w)
{
    return (NV_DATA_S(x).array() * NV_DATA_S(w).array()).matrix().norm();
}

realtype N_VL1Norm_Eigen(N_Vector x)
{
    return NV_DATA_S(x).lpNorm<1>();
}

void N_VCompare_Eigen(double c, N_Vector x, N_Vector z)
{
    Eigen::Array<bool, Eigen::Dynamic, 1> temp = NV_DATA_S(x).array().abs() >= c;
    NV_DATA_S(z) = temp.cast<double>();
}

booleantype N_VInvTest_Eigen(N_Vector x, N_Vector z)
{
    booleantype no_zero_found = true;
    for (int i = 0; i < NV_DATA_S(x).size(); i++) {
        if (NV_DATA_S(x)[i] == 0.0)
            no_zero_found = false;
        else
            NV_DATA_S(z)[i] = 1.0/NV_DATA_S(x)[i];
    }
    return no_zero_found;
}

booleantype N_VConstrMask_Eigen(N_Vector c, N_Vector x, N_Vector m)
{
    sunindextype N = NV_LENGTH_S(x);

    realtype* xd = N_VGetArrayPointer(x);
    realtype* cd = N_VGetArrayPointer(c);
    realtype* md = N_VGetArrayPointer(m);

    realtype temp = 0.0;
    booleantype test;
    for (int i = 0; i < N; i++) {
        md[i] = 0.0;

        /* Continue if no constraints were set for the variable */
        if (cd[i] == 0.0)
            continue;

        /* Check if a set constraint has been violated */
        test = (fabs(cd[i]) > 1.5 && xd[i]*cd[i] <= 0.0) ||
               (fabs(cd[i]) > 0.5 && xd[i]*cd[i] <  0.0);
        if (test) {
            temp = md[i] = 1.0;
        }
    }

    // Return false if any constraint was violated
    return temp != 1.0;
}

realtype N_VMinQuotient_Eigen(N_Vector num, N_Vector denom)
{
    sunindextype N  = NV_LENGTH_S(num);
    realtype* nd = N_VGetArrayPointer(num);
    realtype* dd = N_VGetArrayPointer(denom);

    booleantype notEvenOnce = true;
    realtype min = BIG_REAL;

    for (sunindextype i = 0; i < N; i++) {
        if (dd[i] == 0.0) continue;
        else {
            if (!notEvenOnce) min = SUNMIN(min, nd[i]/dd[i]);
            else {
                min = nd[i]/dd[i];
                notEvenOnce = false;
            }
        }
    }
    return min;
}

realtype N_VMin_Eigen(N_Vector x)
{
    return NV_DATA_S(x).minCoeff();
}

realtype* N_VGetArrayPointer_Eigen(N_Vector v)
{
    return NV_DATA_S(v).data();
}

N_Vector N_VCloneEmpty_Eigen(N_Vector w)
{

    if (w == nullptr) return(nullptr);

    // New generic vector
    N_Vector v = nullptr;
    v = N_VNewEmpty();
    if (v == nullptr) return nullptr;

    // Attach operations
    if (N_VCopyOps(w, v)) { N_VDestroy(v); return(nullptr); }

    // Create content
    N_VectorContent_Eigen content = nullptr;
    content = static_cast<N_VectorContent_Eigen>(malloc(sizeof *content));
    if (content == nullptr) {N_VDestroy(v); return(nullptr);}

    // Attach content
    v->content = content;

    // Fill content
    content->length = NV_LENGTH_S(w);
    content->isOwnData = SUNFALSE;
    content->vecPtr = nullptr;

    // Returned the empty cloned vector
    return v;
}

N_Vector N_VClone_Eigen(N_Vector w)
{
    // Clone w
    N_Vector v = nullptr;
    v = N_VCloneEmpty_Eigen(w);
    if (v == nullptr) return(nullptr);

    sunindextype length = NV_LENGTH_S(w);

    // Allocate memory
    if (length > 0) {
        reinterpret_cast<N_VectorContent_Eigen>(v->content)->vecPtr = new Eigen::VectorXd(length);
        NV_OWN_DATA_S(v) = true;
    }

    return(v);
}

void N_VDestroy_Eigen(N_Vector v)
{

    if (v == nullptr) return;

    /* free content */
    if (v->content != nullptr) {
        /* free data array if it's owned by the vector */
        if (NV_OWN_DATA_S(v)) {
            delete reinterpret_cast<N_VectorContent_Eigen>(v->content)->vecPtr;
            reinterpret_cast<N_VectorContent_Eigen>(v->content)->vecPtr = nullptr;
        }
        free(v->content);
        v->content = nullptr;
    }

    /* free ops and vector */
    if (v->ops != nullptr) { free(v->ops); v->ops = nullptr; }
    free(v); v = nullptr;
}

sunindextype N_VGetLength_Eigen(N_Vector v)
{
    return NV_LENGTH_S(v);
}

sunindextype NV_LENGTH_S(N_Vector v)
{
    return NV_DATA_S(v).size();
}

Eigen::VectorXd& NV_DATA_S(N_Vector v)
{
    return *(reinterpret_cast<N_VectorContent_Eigen>(v->content)->vecPtr);
}

bool& NV_OWN_DATA_S(N_Vector v)
{
    return reinterpret_cast<N_VectorContent_Eigen>(v->content)->isOwnData;
}

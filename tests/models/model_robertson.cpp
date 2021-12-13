#include "model_robertson.h"

using namespace suneigen;

namespace suneigen::model_robertson{

    void xdot_model_robertson(realtype* xdot, const realtype /*t*/, const realtype* x,
                              const realtype *p, const realtype */*k*/, const realtype */*h*/) {
        xdot[0] = -p[0]*x[0] + p[1]*x[1]*x[2];
        xdot[2] = p[2]*x[1]*x[1];
        xdot[1] = -xdot[0] - xdot[2];
    }

    void JSparse_model_robertson(SUNMatrixContent_Sparse JSparse, realtype /*t*/, const realtype *x,
                                 const realtype *p, const realtype */*k*/, const realtype */*h*/){

        JSparse->indexvals[0] = 0;
        JSparse->indexvals[1] = 1;
        JSparse->indexvals[2] = 2;
        JSparse->indexvals[3] = 0;
        JSparse->indexvals[4] = 1;
        JSparse->indexvals[5] = 2;
        JSparse->indexvals[6] = 0;
        JSparse->indexvals[7] = 1;
        JSparse->indexvals[8] = 2;

        JSparse->indexptrs[0] = 0;
        JSparse->indexptrs[1] = 3;
        JSparse->indexptrs[2] = 6;
        JSparse->indexptrs[3] = 9;

        JSparse->data[0] = -p[0];
        JSparse->data[1] = p[0];
        JSparse->data[2] = 0.0;
        JSparse->data[3] = p[1]*x[2];
        JSparse->data[4] = - p[1]*x[2] - 2*p[2]*x[1];
        JSparse->data[5] = 2*p[2]*x[1];
        JSparse->data[6] = p[1]*x[1];
        JSparse->data[7] = -p[1]*x[1];
        JSparse->data[8] = 0.0;
    }

    void fx0_model_robertson(realtype* x0, const realtype /*t*/, const realtype */*p*/, const realtype *k){
        x0[0] = k[0];
        x0[1] = 0.0;
        x0[2] = 0.0;
    }

    void root_model_robertson(realtype *root, const realtype /*t*/, const realtype *x,
                              const realtype */*p*/, const realtype */*k*/, const realtype */*h*/) {
        root[0] = x[2] - RCONST(0.01);
        root[1] = RCONST(0.0001) - x[0];
    }

    void fz_model_robertson(double *z, const unsigned int ie, const realtype t, const realtype */*x*/,
                            const realtype */*p*/, const realtype */*k*/, const realtype */*h*/){
        switch(ie){
            case 0: {
                z[0] = t;
                break;
            }

            case 1: {
                z[1] = t;
                break;
            }
        }
    }

} // namespace suneigen

namespace suneigen::generic_model {

    std::unique_ptr<suneigen::Model> getModel() {
        return std::unique_ptr<suneigen::Model>(
                new suneigen::model_robertson::Model_robertson());
    }

} // namespace suneigen

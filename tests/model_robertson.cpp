#include "model_robertson.h"

using namespace suneigen;

namespace suneigen::model_robertson{

    void xdot_model_robertson(realtype* xdot, const realtype t, const realtype* x) {
        (void) t;
        double y1 = x[0];
        double y2 = x[1];
        double y3 = x[2];
        xdot[0] = RCONST(-0.04)*y1 + RCONST(1.0e4)*y2*y3;
        xdot[2] = RCONST(3.0e7)*y2*y2;
        xdot[1] = -xdot[0] - xdot[2];

    }

    void JSparse_model_robertson(SUNMatrixContent_Sparse JSparse, realtype t, const realtype *x){
        (void) t;
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

        JSparse->data[0] = -0.04;
        JSparse->data[1] = 0.04;
        JSparse->data[2] = 0.0;
        JSparse->data[3] = 1.0e4*x[2];
        JSparse->data[4] = (RCONST(-1.0e4)*x[2]) - (RCONST(6.0e7)*x[1]);
        JSparse->data[5] = RCONST(6.0e7)*x[1];
        JSparse->data[6] = RCONST(1.0e4)*x[1];
        JSparse->data[7] = RCONST(-1.0e4)*x[1];
        JSparse->data[8] = 0.0;
    }

    void fx0_model_robertson(realtype* x0, const realtype t){
        (void) t;
        x0[0] = 1.0;
        x0[1] = 0;
        x0[2] = 0;
    }
} // namespace suneigen

namespace suneigen::generic_model {

    std::unique_ptr<suneigen::Model> getModel() {
        return std::unique_ptr<suneigen::Model>(
                new suneigen::model_robertson::Model_robertson());
    }

} // namespace suneigen

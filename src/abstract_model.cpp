#include "abstract_model.h"

namespace suneigen {

    void AbstractModel::fx0(realtype* /*x0*/,
                       const realtype /*t*/,
                       const realtype* /*p*/,
                       const realtype* /*k*/)
    {
        throw SunException("Requested functionality is not supported as %s is "
                           "not implemented for this model!",
                           __func__);
    }

    void AbstractModel::fsx0(realtype* /*sx0*/,
                        const realtype /*t*/,
                        const realtype* /*x0*/,
                        const realtype* /*p*/,
                        const unsigned int /*ip*/)
    {
        throw SunException("Requested functionality is not supported as %s is "
                           "not implemented for this model!",
                           __func__);
    }

    void AbstractModel::fstau(realtype* /*stau*/,
                         const realtype /*t*/,
                         const realtype* /*x*/,
                         const realtype* /*p*/,
                         const realtype* /*h*/,
                         const realtype* /*sx*/,
                         const unsigned int /*ip*/,
                         const unsigned int /*ie*/)
    {
        throw SunException("Requested functionality is not supported as %s is "
                           "not implemented for this model!",
                           __func__);
    }

    void AbstractModel::fdeltax(realtype* /*deltax*/,
                           const realtype /*t*/,
                           const realtype* /*x*/,
                           const realtype* /*p*/,
                           const realtype* /*h*/,
                           const unsigned int /*ie*/,
                           const realtype* /*xdot*/,
                           const realtype* /*xdot_old*/)
    {
        throw SunException("Requested functionality is not supported as %s is "
                           "not implemented for this model!",
                           __func__);
    }


    void AbstractModel::fz(realtype* /*z*/,
                           const unsigned int /*ie*/,
                           const realtype /*t*/,
                           const realtype* /*x*/,
                           const realtype* /*p*/,
                           const realtype* /*k*/,
                           const realtype* /*h*/)
    {
        throw SunException("Requested functionality is not supported as %s is "
                           "not implemented for this model!",
                           __func__);
    }

}

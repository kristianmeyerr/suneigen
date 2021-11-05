#include "vector.h"

#include <functional>

namespace suneigen {

    Vector &Vector::operator=(Vector const &other) {
        vec_ = other.vec_;
        synchroniseNVector();
        return *this;
    }

    realtype *Vector::data() { return vec_.data(); }

    const realtype *Vector::data() const { return vec_.data(); }

    N_Vector Vector::getNVector() { return nvec_; }

    const_N_Vector Vector::getNVector() const { return nvec_; }

    std::vector<realtype> const &Vector::getVector() const { return vec_; }

    int Vector::getLength() const { return static_cast<int>(vec_.size()); }

    void Vector::zero() { set(0.0); }

    void Vector::minus() {
        std::transform(vec_.begin(), vec_.end(),
                       vec_.begin(), std::negate<realtype>());
    }

    void Vector::set(realtype val) { std::fill(vec_.begin(), vec_.end(), val); }

    realtype &Vector::operator[](int pos) {
        return vec_.at(static_cast<decltype(vec_)::size_type>(pos));
    }

    realtype &Vector::at(int pos) {
        return vec_.at(static_cast<decltype(vec_)::size_type>(pos));
    }

    const realtype &Vector::at(int pos) const {
        return vec_.at(static_cast<decltype(vec_)::size_type>(pos));
    }

    void Vector::copy(const Vector &other) {
        if(getLength() != other.getLength())
            throw SunException("Dimension of Vector (%i) does not "
                               "match input dimension (%i)",
                               getLength(), other.getLength());
        std::copy(other.vec_.begin(), other.vec_.end(), vec_.begin());
        synchroniseNVector();
    }

    void Vector::synchroniseNVector() {
        if (nvec_)
            N_VDestroy_Serial(nvec_);
        nvec_ = N_VMake_Serial(static_cast<long int>(vec_.size()), vec_.data());
    }

    Vector::~Vector() {
        if (nvec_)
            N_VDestroy_Serial(nvec_);
    }

    VectorArray::VectorArray(long int length_inner, long int length_outer)
            : vec_array_(length_outer, Vector(length_inner)) {
        nvec_array_.resize(length_outer);
        for (int idx = 0; idx < length_outer; idx++) {
            nvec_array_.at(idx) = vec_array_.at(idx).getNVector();
        }
    }

    VectorArray &VectorArray::operator=(VectorArray const &other) {
        vec_array_ = other.vec_array_;
        nvec_array_.resize(other.getLength());
        for (int idx = 0; idx < other.getLength(); idx++) {
            nvec_array_.at(idx) = vec_array_.at(idx).getNVector();
        }
        return *this;
    }

    VectorArray::VectorArray(const VectorArray &vaold)
            : vec_array_(vaold.vec_array_) {
        nvec_array_.resize(static_cast<unsigned long>(vaold.getLength()));
        auto n = static_cast<unsigned long>(vaold.getLength());
        for (unsigned long idx = 0; idx < n; idx++) {
            nvec_array_.at(idx) = vec_array_.at(idx).getNVector();
        }
    }

    realtype *VectorArray::data(int pos) { return vec_array_.at(static_cast<unsigned long>(pos)).data(); }

    const realtype *VectorArray::data(int pos) const {
        return vec_array_.at(static_cast<unsigned long>(pos)).data();
    }

    realtype &VectorArray::at(int ipos, int jpos) {
        return vec_array_.at(static_cast<unsigned long>(jpos)).at(ipos);
    }

    const realtype &VectorArray::at(int ipos, size_t jpos) const {
        return vec_array_.at(jpos).at(ipos);
    }

    N_Vector *VectorArray::getNVectorArray() { return nvec_array_.data(); }

    N_Vector VectorArray::getNVector(int pos) { return nvec_array_.at(static_cast<unsigned long>(pos)); }

    const_N_Vector VectorArray::getNVector(int pos) const { return nvec_array_.at(static_cast<unsigned long>(pos)); }

    Vector &VectorArray::operator[](int pos) { return vec_array_.at(static_cast<unsigned long>(pos)); }

    const Vector &VectorArray::operator[](int pos) const {
        return vec_array_.at(static_cast<unsigned long>(pos));
    }

    int VectorArray::getLength() const {
        return static_cast<int>(vec_array_.size());
    }

    void VectorArray::zero() {
        for (auto &v : vec_array_)
            v.zero();
    }

    void VectorArray::flatten_to_vector(std::vector<realtype> &vec) const {

        size_t n_outer = vec_array_.size();
        if (n_outer == 0)
            return; // nothing to do ...

        auto n_inner = static_cast<size_t>(vec_array_.at(0).getLength());

        if (vec.size() != n_inner * n_outer) {
            throw SunException("Dimension of VectorArray (%ix%i) does not "
                               "match target vector dimension (%u)",
                               n_inner, n_outer, vec.size());
        }

        for (size_t outer = 0; outer < n_outer; ++outer) {
            for (size_t inner = 0; inner < n_inner; ++inner)
                vec.at(inner + outer * n_inner)  =static_cast<size_t>(this->at(inner, outer));
        }
    }

    void VectorArray::copy(const VectorArray &other) {
        if (getLength() != other.getLength())
            throw SunException("Dimension of VectorArray (%i) does not "
                               "match input dimension (%i)",
                               getLength(), other.getLength());

        for (int iv = 0; iv < getLength(); ++iv) {
            vec_array_.at(iv).copy(other.vec_array_.at(iv));
            nvec_array_[iv] = vec_array_.at(iv).getNVector();
        }
    }

} // namespace amici
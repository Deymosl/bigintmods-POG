//
// Created by User on 19.10.2018.
//

#include "vector.h"
#include <cstring>

fVector::fVector(size_t n, uint32_t q) {
    sz = n;
    mode = true;
    if (n > SZ) {
        becomeBig(ensureCapacity(sz), SZ);
    }
    std::fill(getPtr(), getPtr() + sz, q);
}

fVector::fVector(const fVector &other) {
    sz = other.sz;
    mode = other.mode;
    if (mode) {
        memcpy(data.arr, other.data.arr, sz * sizeof(uint32_t));
    } else {
        new(&data.vec) vector(other.data.vec);
    }
}

fVector& fVector::operator=(fVector other) {
    swap(other);
    return *this;
}

fVector::~fVector() {
    if (!mode) {
        data.vec.~vector();
    }
}

void fVector::push_back(uint32_t q) {
    check();
    if (sz == getCapacity()) {
        becomeBig(ensureCapacity(sz + 1), sz);
    }
    operator[](sz++) = q;
}

void fVector::pop_back() {
    check();
    sz--;
}

uint32_t &fVector::back() {
    return operator[](sz - 1);
}

const uint32_t &fVector::back() const {
    return operator[](sz - 1);
}

uint32_t &fVector::operator[](size_t ind) {
    check();
    return *(getPtr() + ind);
}

const uint32_t &fVector::operator[](size_t pos) const {
    return *(getPtr() + pos);
}

size_t fVector::size() const {
    return sz;
}

bool operator==(const fVector &q, const fVector &w) {
    if (q.size() != w.size())
        return false;
    return memcmp(q.getPtr(), w.getPtr(), sizeof(uint32_t) * q.size()) == 0;
}

void fVector::swap(fVector &other) {
    size_t temp = sz; // std::swap isn't working, hotfix incaming
    sz = other.sz;
    other.sz = temp;
    bool m = mode;
    mode = other.mode;
    other.mode = m;
    static char buffer[sizeof(Data)];
    memcpy(buffer, &data, sizeof(Data));
    memcpy(&data, &other.data, sizeof(Data));
    memcpy(&other.data, buffer, sizeof(Data));
}

size_t fVector::getCapacity() {
    return (mode ? SZ : data.vec.cap);
}

void fVector::check() {
    if (mode || data.vec.ptr.unique()) {
        return;
    }
    becomeBig(data.vec.cap, data.vec.cap);
}

uint32_t *fVector::getPtr() const {
    return (mode ? const_cast<uint32_t *>(data.arr) : data.vec.ptr.get());
}

size_t fVector::ensureCapacity(size_t n) {
    return (n < CAP ? CAP : n * 2);
}

void fVector::becomeBig(size_t ensured, size_t copy) {
    auto tmp = new uint32_t[ensured];
    memcpy(tmp, getPtr(), copy * sizeof(uint32_t));
    if (!mode)
        data.vec.~vector();
    else {
        data.~Data();
    }
    new(&data.vec) vector(ensured, std::shared_ptr<uint32_t>(tmp, std::default_delete<uint32_t[]>()));
    mode = false;
}

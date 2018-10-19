//
// Created by User on 19.10.2018.
//

#ifndef BIGINT_VECTOR_H
#define BIGINT_VECTOR_H

#include <cstdlib>
#include <memory>
#include <cstring>

const size_t CAP = 16;
const size_t SZ = 8;

class fVector {
public:

    fVector() = default;

    fVector(size_t n, uint32_t q);

    fVector(const fVector &source);

    fVector &operator=(fVector other);

    ~fVector();

    void push_back(uint32_t q);

    void pop_back();

    uint32_t &back();

    const uint32_t &back() const;

    uint32_t &operator[](size_t pos);

    const uint32_t &operator[](size_t pos) const;

    size_t size() const;

    friend bool operator==(fVector const &q, fVector const &w);

    void swap(fVector &dest);

    size_t getCapacity();

private:
    size_t sz = 0;
    bool mode = true;

    void check();

    uint32_t *getPtr() const;

    size_t ensureCapacity(size_t n);

    void becomeBig(size_t ensure, size_t copy);

    struct vector {
        size_t cap = CAP;
        std::shared_ptr<uint32_t> ptr = nullptr;
        vector(size_t capacity, std::shared_ptr<uint32_t> pointer) : cap(capacity), ptr(pointer) {}
        vector() = default;
    };

    union Data {
        uint32_t arr[SZ];
        vector vec;

        Data() {}

        ~Data() { }
    } data;
};

#endif //BIGINT_VECTOR_H

//
//

#ifndef BIG_INT_BIGINT_H
#define BIG_INT_BIGINT_H

#include <algorithm>
#include <string>
#include <vector>
#include <climits>
#include <cmath>
#include <cstdint>
#include "vector.h"

class BigInteger {
public:
    BigInteger();

    BigInteger(const BigInteger &source);

    BigInteger(int source);

    explicit BigInteger(std::string source);

    ~BigInteger() = default;

    BigInteger &operator=(const BigInteger &source) = default;

    BigInteger operator+() const;

    BigInteger operator-() const;

    BigInteger operator~() const;

    BigInteger &operator+=(const BigInteger &source);

    BigInteger &operator-=(const BigInteger &source);

    BigInteger &operator/=(const BigInteger &source);

    BigInteger &operator*=(const BigInteger &source);

    BigInteger &operator%=(const BigInteger &source);

    BigInteger &operator|=(const BigInteger &source);

    BigInteger &operator&=(const BigInteger &source);

    BigInteger &operator^=(const BigInteger &source);

    BigInteger &operator<<=(int rhs);

    BigInteger &operator>>=(int lhs);

    friend BigInteger operator+(BigInteger q, const BigInteger &w);

    friend BigInteger operator-(BigInteger q, const BigInteger &w);

    friend BigInteger operator*(BigInteger q, const BigInteger &w);

    friend BigInteger operator/(BigInteger q, const BigInteger &w);

    friend BigInteger operator%(BigInteger q, const BigInteger &w);


    friend BigInteger operator&(BigInteger q, const BigInteger &w);

    friend BigInteger operator|(BigInteger q, const BigInteger &w);

    friend BigInteger operator^(BigInteger q, const BigInteger &w);

    friend BigInteger operator<<(BigInteger, int w);

    friend BigInteger operator>>(BigInteger, int w);

    friend bool operator==(const BigInteger &q, const BigInteger &w);

    friend bool operator!=(const BigInteger &q, const BigInteger &w);

    friend bool operator>=(const BigInteger &q, const BigInteger &w);

    friend bool operator<=(const BigInteger &q, const BigInteger &w);

    friend bool operator>(const BigInteger &q, const BigInteger &w);

    friend bool operator<(const BigInteger &q, const BigInteger &w);

    friend std::string to_string(const BigInteger &q);

    friend std::ostream &operator<<(std::ostream &s, const BigInteger &q);

    size_t length() const;

    uint32_t operator[](size_t q) const;

private:
    bool sign;
    fVector data;

    BigInteger inverse() const;

    fVector quotient(uint32_t q) const;

    template<class Func>
    BigInteger &applyBitwiseOperation(const BigInteger &q, Func function);

    void shrink();
};

fVector &correct(fVector &q);

void divide(const fVector &q, const fVector &w, fVector &res);

bool isLess(const fVector &q, const fVector &w);

void add(const fVector &q, const fVector &w, fVector &res);

void sub(const fVector &q, const fVector &w, fVector &res);

void mul(const fVector &q, const fVector &w, fVector &res);

#endif //BIG_INT_BIGINT_H

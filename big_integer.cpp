#include <iostream>
#include <cstdint>
#include <stdint-gcc.h>
#include "big_integer.h"

typedef unsigned long long ull;

uint32_t get(const fVector &q, size_t w) {
    return (w < q.size() ? q[w] : 0);
}

void copy(int l, int r, int from, const fVector &source, fVector &dest) {
    for (size_t i = l; i < r; i++) {
        dest[from + i - l] = source[i];
    }
}

size_t BigInteger::length() const {
    return data.size();
}

void BigInteger::shrink() {
    correct(data);
    if (data.back() == 0)
        sign = false;
}

fVector &correct(fVector &q) {
    static int cnt = 0;
    while (true) {
        cnt++;
        if (q.size() <= 1 || q.back() != 0)
            break;
        q.pop_back();
    }
    return q;
}

BigInteger::BigInteger() {
    sign = false;
    data.push_back(0);
}

BigInteger::BigInteger(int source) {
    sign = source < 0;
    if (sign) {
        source = -(source + 1);
        data.push_back(static_cast<uint32_t>(source) + 1);
    } else data.push_back(source);
}

BigInteger::BigInteger(ull q) {
    sign = false;
    while (q) {
        data.push_back(q & UINT32_MAX);
        q >>= sizeof(uint32_t) * 8;
    }
}

BigInteger::BigInteger(unsigned int q) {
    sign = false;
    data.push_back(static_cast<uint32_t>(q));
}

BigInteger::BigInteger(const BigInteger &source) {
    sign = source.sign;
    data = source.data;
}

BigInteger::BigInteger(std::string source) {
    sign = false;
    BigInteger result(0);
    reverse(source.begin(), source.end());
    if (source.back() == '-') {
        source.pop_back();
        sign = true;
    }
    while (!source.empty()) {
        result *= 10;
        result += source.back() - '0';
        source.pop_back();
    }
    data = result.data;
    shrink();
}

BigInteger BigInteger::operator-() const {
    if (*this == 0)
        return *this;
    BigInteger result = BigInteger(*this);
    result.sign = !sign;
    return result;
}

BigInteger BigInteger::operator+() const {
    return *this;
}

BigInteger BigInteger::operator~() const {
    BigInteger res(*this);
    return -(res + 1);
}

BigInteger &BigInteger::operator+=(const BigInteger &q) {
    if (sign == q.sign) {
        add(data, q.data, data);
    } else {
        bool s = sign;
        const BigInteger &a1 = s ? q : *this;
        const BigInteger &a2 = !s ? q : *this;
        sign = isLess(a1.data, a2.data);
        sub(a1.data, a2.data, data);
    }
    return *this;
}

fVector BigInteger::quotient(uint32_t q) const {
    BigInteger res;
    res.data = fVector(data.size(), 0);
    uint64_t carry = 0;
    for (int i = data.size() - 1; i >= 0; i--) {
        uint64_t temp = data[i] + (carry << 32);
        res.data[i] = temp / q;
        carry = temp % q;
    }
    correct(res.data);
    return res.data;
}

bool isLess(const fVector &q, const fVector &w) {
    if (q.size() < w.size())
        return true;
    if (q.size() > w.size())
        return false;
    for (int i = q.size() - 1; i >= 0; i--) {
        if (q[i] < w[i])
            return true;
        else if (q[i] > w[i])
            return false;
    }
    return false;
}

void mul(const fVector &q, const fVector &w, fVector &dest) {
    int n = q.size();
    fVector res(q.size() + w.size(), 0);
    res.back() = 0;
    uint32_t carry = 0;
    int j = 0;
    uint64_t cur = 0;
    for (int i = n - 1; i >= 0; i--) {
        uint32_t t = q[i];
        for (j = 0, carry = 0; j < (int) w.size() || carry > 0; j++) {
            cur = (uint64_t) res[i + j] * (j != 0) + (uint64_t) t * (uint64_t) get(w, j);
            res[i + j] = cur + carry;
            carry = (cur + carry) >> 32;
        }
    }
    correct(res);
    dest = res;
}

void fmul(const fVector &q, uint32_t w, fVector &dest) {
    dest = q;
    uint32_t carry = 0;
    for (size_t i = 0 ; i < q.size(); i++) {
        uint64_t cur = carry + q[i] *1ull * w;
        dest[i] = (uint32_t) (cur);
        carry = (uint32_t) (cur >> 32);
    }
    if (carry)
        dest.push_back(carry);
}

void add(const fVector &q, const fVector &w, fVector &dest) {
    fVector res(std::max(q.size(), w.size()) + 1, 0);
    uint32_t carry = 0;
    for (unsigned int i = 0; i < std::max(q.size(), w.size()); i++) {
        uint32_t qi = get(q, i);
        uint32_t wi = get(w, i);
        res[i] = qi + wi + carry;
        carry = (UINT32_MAX - qi < wi + carry); //hm
    }
    if (carry)
        res.back() = 1;
    correct(res);
    dest = res;
}

void sub(const fVector &q, const fVector &w, fVector &dest) {
    bool less = isLess(q, w);
    const fVector &a = (less ? w : q);
    const fVector &b = (less ? q : w);
    fVector res(std::max(a.size(), b.size()), 0);
    for (int i = a.size() - 1; i >= 0; i--) {
        uint32_t ai = get(a, i);
        uint32_t bi = get(b, i);
        res[i] = ai - bi + (ai < bi) * UINT32_MAX + (ai < bi);
        if (ai < bi) {
            int j = i + 1;
            while (res[j] == 0) {
                res[j]--;
                j++;
            }
            res[j]--;
        }
    }
    correct(res);
    dest = res;
}

void divide(const fVector &q, const fVector &w, fVector &dest) {
    const unsigned int fac = (unsigned int) ((1ll << 32) / ((unsigned long long) (w.back() + 1)));
    fVector res(0, 0);
    fVector qFactorized;
    fmul(q, fac, qFactorized);
    fVector wFactorized;
    fmul(w, fac, wFactorized);
    int n = qFactorized.size();
    fVector temp;
    for (int i = n - 1; i >= 0; i--) {
        temp.push_back(qFactorized[i]);
        for (int s = temp.size() - 1; s > 0; s--) {
            temp[s] = temp[s - 1];
        }
        temp[0] = qFactorized[i];
        while (temp.size() > 1 && temp.back() == 0) {
            temp.pop_back();
        }

        if (isLess(temp, wFactorized)) {
            if (temp.back() == 0)
                res.push_back(0);
            continue;
        }
        uint64_t dividend = uint64_t(temp.back());
        if (temp.size() != wFactorized.size()) {
            dividend <<= 32;
            dividend += temp[temp.size() - 2];
        }
        uint64_t divider = uint64_t(wFactorized.back());
        uint64_t quot = dividend / divider;
        quot = std::min(quot, ((uint64_t) 1 << 32) - 1);
        fVector buffer;
        fmul(wFactorized, quot, buffer);
        if (isLess(temp, buffer)) {
            sub(buffer, wFactorized, buffer);
            quot--;
        }
        if (isLess(temp, buffer)) {
            uint64_t l = 0, r = quot;
            while (r - l > 1) {
                uint64_t mid = (r + l) / 2;
                fVector tt;
                fmul(wFactorized, mid, tt);
                if (isLess(temp, tt)) {
                    r = mid;
                } else {
                    l = mid;
                }
            }
            quot = l;
            mul(wFactorized, fVector(1, quot), buffer);
        }
        sub(temp, buffer, temp);
        res.push_back(quot);
    }
    for (size_t i = 0; i < res.size() / 2; i++)
    {
        uint32_t f = res[i];
        res[i] = res[res.size() - 1 - i];
        res[res.size() - 1 - i] = f;
    }
    dest = res;
}

BigInteger &BigInteger::operator/=(const BigInteger &q) {
    if (q.data.size() == 1) {
        data = quotient(q.data[0]);
        sign ^= q.sign;
    } else if (q != 0) {
        if (isLess(data, q.data)) {
            data = fVector(1, 0);
            sign = false;
        } else {
            sign ^= q.sign;
            fVector res;
            divide(data, q.data, res);
            data = res;
        }
    }
    shrink();
    return *this;
}

BigInteger &BigInteger::operator*=(const BigInteger &q) {
    sign ^= q.sign;
    mul(data, q.data, data);
    if (data.back() == 0)
        sign = false;
    return *this;
}

BigInteger &BigInteger::operator%=(const BigInteger &q) {
    BigInteger temp(*this);
    temp /= q;
    temp *= q;
    *this -= temp;
    return *this;
}

bool operator==(const BigInteger &q, const BigInteger &w) {
    if (q.sign != w.sign || q.length() != w.length())
        return false;
    for (unsigned int i = q.length(); i > 0; i--) {
        if (q[i - 1] != w[i - 1])
            return false;
    }
    return true;
}

bool operator!=(const BigInteger &q, const BigInteger &w) {
    return !(q == w);
}

bool operator<(const BigInteger &q, const BigInteger &w) {
    if (q.sign != w.sign)
        return (q.sign > w.sign);
    if (q.length() != w.length()) {
        return (q.length() < w.length()) ^ q.sign;
    }
    return isLess(q.sign ? w.data : q.data, q.sign ? q.data : w.data);
}

bool operator>(const BigInteger &q, const BigInteger &w) {
    return w < q;
}

bool operator>=(const BigInteger &q, const BigInteger &w) {
    return !(q < w);
}

bool operator<=(const BigInteger &q, const BigInteger &w) {
    return !(q > w);
}

template<class Func>
BigInteger &BigInteger::applyBitwiseOperation(const BigInteger &q, Func function) {
    fVector res = fVector(std::max(data.size(), q.data.size()), 0);
    if (sign)
        data = inverse().data;
    const BigInteger &right = (q.sign ? q.inverse() : q);
    for (size_t i = 0; i < length(); i++)
        data[i] = function(data[i], right[i]);
    correct(data);
    sign = function(sign, q.sign);
    if (sign)
        data = inverse().data;
    return *this;

}

BigInteger BigInteger::inverse() const {
    BigInteger res(*this);
    for (size_t i = 0; i < res.data.size(); i++)
        res.data[i] = UINT32_MAX - res.data[i];
    res.sign ^= true;
    return res + 1;
}

BigInteger &BigInteger::operator&=(const BigInteger &q) {
    return applyBitwiseOperation(q, std::bit_and<uint32_t>());
}

BigInteger operator&(BigInteger q, const BigInteger &w) {
    return q &= w;
}

BigInteger &BigInteger::operator|=(const BigInteger &q) {
    return applyBitwiseOperation(q, std::bit_or<uint32_t>());
}

BigInteger operator|(BigInteger q, const BigInteger &w) {
    return q |= w;
}

BigInteger &BigInteger::operator^=(const BigInteger &q) {
    return applyBitwiseOperation(q, std::bit_xor<uint32_t>());
}

BigInteger operator^(BigInteger q, const BigInteger &w) {
    return q ^= w;
}

int getDeg(uint32_t q) {
    int res = 0;
    while (q) {
        res++;
        q /= 2;
    }
    return res;
}

BigInteger &BigInteger::operator<<=(unsigned int q) {
    int maxDeg = getDeg(data.back());
    int rem = q % 32;
    int carry = (maxDeg + rem > 32);
    fVector res(data.size() + q / 32 + carry, 0);
    copy(0, data.size(), 0, data, res);
    for (int i = (int) res.size() - 1; i >= (int) res.size() - (int) data.size(); i--) {
        res[i] <<= rem;
        if (i != (int) res.size() - (int) data.size())
            res[i] += res[i - 1] >> (32 - rem);
    }
    data = res;
    return *this;
}

BigInteger operator<<(BigInteger q, int w) {
    return (q <<= w);
}

void shift_right(const fVector &q, int w, fVector &dest) {
    fVector res = fVector(std::max(0, (int) q.size() - w / 32), 0);
    copy(w / 32, q.size(), 0, q, res);
    if (res.size() == 0)
        res.push_back(0);
    else {
        uint32_t next = 0;
        for (int i = 0; i < w; i++) next *= 2, next += 1;
        for (int i = 0; i < res.size(); i++) {
            res[i] >>= w;
            if (i != res.size() - 1)
                res[i] += (res[i + 1] & next) << (sizeof(uint32_t) * 8 - w);
        }
    }
    correct(res);
    dest = res;
}

BigInteger &BigInteger::operator>>=(unsigned int q) {
    shift_right(data, q, data);
    if (sign)
        add(data, fVector(1, 1), data);
    shrink();
    return *this;
}

BigInteger operator>>(BigInteger q, int w) {
    BigInteger f = (q >>= w);
    return f;
}

uint32_t BigInteger::operator[](size_t pos) const {
    if (pos < length())
        return data[pos];
    return 0;
}

BigInteger &BigInteger::operator-=(const BigInteger &q) {
    return *this += -q;
}

BigInteger operator+(BigInteger q, const BigInteger &w) {
    return q += w;
}

BigInteger operator-(BigInteger q, const BigInteger &w) {
    return q -= w;
}

BigInteger operator*(BigInteger q, const BigInteger &w) {
    return q *= w;
}

BigInteger operator/(BigInteger q, const BigInteger &w) {
    return q /= w;
}

BigInteger operator%(BigInteger q, const BigInteger &w) {
    return q %= w;
}

std::string to_string(BigInteger const &q) {
    if (q == 0)
        return "0";
    std::string res = "";
    BigInteger w = q;
    while (w != 0) {
        res += to_string((w % 10).data.back());
        w /= 10;
    }
    if (q.sign)
        res.push_back('-');
    reverse(res.begin(), res.end());
    return res;
}

std::ostream &operator<<(std::ostream &os, const BigInteger &q) {
    std::string res = to_string(q);
    os << res;
    return os;
}
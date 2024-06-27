#pragma once

#include <Bow/Types.h>
#include <Bow/Utils/Logging.h>
#include <oneapi/tbb.h>

namespace Bow {

template <class T, class TM, class TV>
class ConjugateGradient {

    /** All notations adopted from Wikipedia,
         * q denotes A*p in general */
    TV r, p, q, temp;
    TV mr, s;

public:
    T tolerance;
    T relative_tolerance;
    int max_iterations;

    ConjugateGradient(const int max_it_input)
        : max_iterations(max_it_input)
    {
        setTolerance(std::is_same<T, float>::value ? (T)1e-6 : (T)1e-12);
    }

    ~ConjugateGradient() {}

    void setTolerance(const T tolerance_input = 16 * std::numeric_limits<T>::epsilon()) { tolerance = tolerance_input; }

    void setRelativeTolerance(const T tolerance_input = 1) { relative_tolerance = tolerance_input; }

    void reinitialize(const TV& b)
    {
        r = b; // r.resizeLike(b);
        p = b; // p.resizeLike(b);
        q = b; // q.resizeLike(b);
        temp = b; // temp.resizeLike(b);

        mr = b; // mr.resizeLike(b);
        s = b; // s.resizeLike(b);
    }

    T dotProduct(const TV& A, const TV& B)
    {
        return (A.array() * B.array()).sum();
    }

    int solve(const TM& A, TV& x, const TV& b, const bool verbose = false)
    {
        reinitialize(x);
        int cnt = 0;
        T alpha, beta, residual_preconditioned_norm, zTrk, zTrk_last;

        //NOTE: requires that the input x has been projected
        A.multiply(x, temp);
        r = b - temp;
        A.project(r);
        A.precondition(r, q); //NOTE: requires that preconditioning matrix is projected
        p = q;

        zTrk = std::abs(dotProduct(r, q));
        residual_preconditioned_norm = std::sqrt(zTrk);
        T local_tolerance = std::min(relative_tolerance * residual_preconditioned_norm, tolerance);
        for (cnt = 0; cnt < max_iterations; ++cnt) {
            if (residual_preconditioned_norm <= local_tolerance) {
                Bow::Logging::info("\tCG terminates at ", cnt, "; (preconditioned norm) residual = ", residual_preconditioned_norm);
                return cnt;
            }

            if (cnt % 50 == 0) {
                Bow::Logging::info("\tCG iter ", cnt, "; (preconditioned norm) residual = ", residual_preconditioned_norm);
            }

            A.multiply(p, temp);
            A.project(temp);
            alpha = zTrk / dotProduct(temp, p);

            x = x + alpha * p;
            r = r - alpha * temp;
            A.precondition(r, q); //NOTE: requires that preconditioning matrix is projected

            zTrk_last = zTrk;
            zTrk = dotProduct(q, r);
            beta = zTrk / zTrk_last;

            p = q + beta * p;

            residual_preconditioned_norm = std::sqrt(zTrk);
        }
        Bow::Logging::info("ConjugateGradient max iterations reached ", max_iterations);
        return max_iterations;
    }
};
} // namespace Bow

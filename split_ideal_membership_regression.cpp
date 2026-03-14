#include <cassert>
#include <iostream>

#include "fast_quaternions.hpp"
#include "quaternions.hpp"

int main()
{
    quatalg Bp {Integer(3), Integer(1)};
    Integer norm(5);
    quatlat O0 = Bp.maximal_order(false);
    Integer a = NTL::SqrRootMod(Integer(norm - Bp.q), norm);

    // Pick a split-case generator that still makes the old comma-expression misbranch:
    // prod[3] == 0 mod norm while prod[0] and prod[1] are non-zero mod norm.
    quat gen = {{Integer(1), norm - a, a, Integer(1), Integer(1)}, Bp};
    assert(O0.contains(gen));
    assert(gen.integral_norm() % norm == 0);
    assert(NTL::GCD(norm, (gen[2] * gen[2] + gen[3] * gen[3]) % norm) == norm);

    quat special_gen = {{a, Integer(1), Integer(0), Integer(0), Integer(1)}, Bp};
    quat prod = gen * special_gen;

    bool old_branch = ((prod[0] % norm == 0 && prod[1] % norm == 0 && prod[2] % norm == 0), prod[3] % norm == 0);
    bool new_branch = (prod[0] % norm == 0 && prod[1] % norm == 0 && prod[2] % norm == 0 && prod[3] % norm == 0);

    assert(old_branch);
    assert(!new_branch);

    quatlat slow = create_from_generator_O0(gen, norm);
    assert(slow.basis[1][0] == 2 * a);
    assert(slow.basis[3][1] == norm - a);
    assert(slow.basis[3][2] == norm - a);

    FastQuatAlg fast_Bp(Bp);
    FastQuat fast_gen(gen, fast_Bp);
    FastQuat fast_special_gen = {{convert(a), 1, 0, 0, 1}, fast_Bp};
    FastQuat fast_prod = fast_gen * fast_special_gen;

    bool fast_new_branch = (fast_prod[0] % convert(norm) == 0 &&
                            fast_prod[1] % convert(norm) == 0 &&
                            fast_prod[2] % convert(norm) == 0 &&
                            fast_prod[3] % convert(norm) == 0);
    assert(!fast_new_branch);

    FastQuatLat fast = create_from_generator_O0(fast_gen, convert(norm));
    assert(fast.basis[1][0] == 2 * convert(a));
    assert(fast.basis[3][1] == convert(norm - a));
    assert(fast.basis[3][2] == convert(norm - a));

    std::cout << "split ideal membership regression passed\n";
    return 0;
}

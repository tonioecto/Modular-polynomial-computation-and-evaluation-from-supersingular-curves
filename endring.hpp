//////////////////////////////////////////////////////////////////////////////////////////////////////
/////// This code includes functions to compute the endomorphism ring of an elliptic curve
//////////////////////////////////////////////////////////////////////////////////////////////////////


#pragma once

#include <NTL/lzz_pXFactoring.h>

#include <iomanip>
#include <queue>

#include "quaternions.hpp"
#include "isog.hpp"
#include "utils.hpp"
#include "Fp2k.hpp"
#include "id2iso.hpp"

static inline Fp_elem tofp(FpE_elem const &v)
{
    assert(NTL::power(v, Fp_elem::modulus()) == v);   // defined over Fp
    if (NTL::IsZero(v))
        return {};
    return NTL::rep(v)[0];
}

bool inline is_surface(ec const &E)
{
    // check curve is defined over Fp
    assert(NTL::power(E.a(), Fp_elem::modulus()) == E.a());
    assert(NTL::power(E.b(), Fp_elem::modulus()) == E.b());

    FpX_elem f;
    NTL::SetCoeff(f, 3);
    f[1] = tofp(E.a());
    f[0] = tofp(E.b());

    auto facs = NTL::NewDDF(f, NTL::PowerXMod(Fp_elem::modulus(), f));

    return facs.length() == 1 && facs[0].b == 1; // splits completely
}

bool inline is_isomorphic(ec const &E1, ec const &E2)
{
    assert(E1.j_invariant() == E2.j_invariant());  // caller must check

    // check curves are defined over Fp
    assert(NTL::power(E1.a(), Fp_elem::modulus()) == E1.a());
    assert(NTL::power(E1.b(), Fp_elem::modulus()) == E1.b());
    assert(NTL::power(E2.a(), Fp_elem::modulus()) == E2.a());
    assert(NTL::power(E2.b(), Fp_elem::modulus()) == E2.b());

    assert(NTL::IsZero(E1.a()) == NTL::IsZero(E2.a()));
    assert(NTL::IsZero(E1.b()) == NTL::IsZero(E2.b()));

    if (!NTL::IsZero(E1.a())) {
        FpX_elem f;
        NTL::SetCoeff(f, 4);
        // std::cerr << " \x1b[35m " << E1 << " " << E2 << "\x1b[0m" << std::endl;
        f[0] = -tofp(E2.a())/tofp(E1.a());
        auto facs = NTL::NewDDF(f, NTL::PowerXMod(Fp_elem::modulus(), f));
        bool good = false;
        for (auto const &[_,d]: facs)
            if (d == 1) {
                good = true;
                break;
            }
        if (!good)
            return false;
    }

    if (!NTL::IsZero(E1.b())) {
        FpX_elem f;
        NTL::SetCoeff(f, 6);
        f[0] = -tofp(E2.b())/tofp(E1.b());
        auto facs = NTL::NewDDF(f, NTL::PowerXMod(Fp_elem::modulus(), f));
        bool good = false;
        for (auto const &[_,d]: facs)
            if (d == 1) {
                good = true;
                break;
            }
        if (!good)
            return false;
    }

    return true;
}

std::pair<ec, quatlat> inline starting_curve(quatalg const &alg, bool surface=false)
{
    if ((alg.q > 3)) {
        std::cout << "!!!!!! CAREFUL! YOU PROBABLY NEED TO CALL THE OTHER FUNCTION WITH FIELD EXTS !!!!!!!" << std::endl;
    }
    ec E0;
    if (!(NTL::IsOne(alg.q))) {
        assert (alg.p % 4 != 3);
        // Everything is on surface in this case
        if (alg.q == 3) {
            FpE_elem a(0), b(1);
            E0 = ec(a, b);
        } else if (alg.q == 7) {
            FpE_elem a(35), b(98);
            a = -a;
            E0 = ec(a, b);
        } else if (alg.q == 11) {
            FpE_elem a(1056), b(13552);
            a = -a;
            E0 = ec(a, b);
        } else if (alg.q == 19) {
            FpE_elem a(152), b(722);
            a = -a;
            E0 = ec(a, b);
        } else {
            throw;
        }
    } else {
        FpE_elem a(surface ? -1 : +1), b;
        E0 = ec(a, b);
        // std::cerr << E0 << std::endl;
    }
    quatlat O0 = alg.maximal_order(surface);
    assert(O0.is_order());

    // computing the norm
    auto test= O0.norm();
    assert(test.first == test.second && test.first == 1);

    return {E0, O0};
}

std::pair<ec, quatlat> inline starting_curve(quatalg const &alg, std::map<unsigned,Fp2k> &Fexts, bool surface=false)
{
    ec E0;
    if (!(NTL::IsOne(alg.q))) {
        assert (alg.p % 4 != 3);
        // Everything is on surface in this case
        if (alg.q == 3) {
            FpE_elem a(0), b(1);
            E0 = ec(a, b);
        } else if (alg.q == 7) {
            FpE_elem a(35), b(98);
            a = -a;
            E0 = ec(a, b);
        } else if (alg.q == 11) {
            FpE_elem a(1056), b(13552);
            a = -a;
            E0 = ec(a, b);
        } else if (alg.q == 19) {
            FpE_elem a(152), b(722);
            a = -a;
            E0 = ec(a, b);
        } else {
            throw;
        }
    } else {
        FpE_elem a(surface ? -1 : +1), b;
        E0 = ec(a, b);
        // std::cerr << E0 << std::endl;
    }
    quat alpha_check(alg);
    quatlat O0 = alg.maximal_order_with_quat_for_check(&alpha_check, surface);
    if (alg.q > 3) {
        int q_int;
        NTL::conv(q_int, alg.q);
        unsigned deg = torsionToFieldDegree(alg.q);
        auto jt = Fexts.find(deg);
        assert(jt != Fexts.end());
        auto bas = E0.torsionBasis(jt->second, q_int, 1);

        if (!(evalEndo(alpha_check, bas.first, alg.q).is_identity()) || !(evalEndo(alpha_check, bas.second, alg.q).is_identity())) {
            quat quat_i({{NTL::ZZ(0), NTL::ZZ(1), NTL::ZZ(0), NTL::ZZ(0), NTL::ZZ(1)}, alg});
            quat quat_i_inv = quat_i;
            quat_i_inv.invert();
            O0 = quat_i*O0*quat_i_inv;
        }
    }

    assert(O0.is_order());

    // computing the norm
    O0.reset_norm();
    auto test= O0.norm();
    assert(test.first == test.second && test.first == 1);

    return {E0, O0};
}

quatlat inline transport_endring(ec const &E0, quatlat const &O0, ec const &E1, std::map<unsigned,Fp2k> &Fexts, uint64_t avoid_prime=0)
{
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    //// Sssumption: j in O0 corresponds to Frobenius; in particular j^2=-p
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////

    using expvec = std::array<int16_t, 9>;  //TODO: figure out how many primes to use

    auto const &p = O0.alg.p;

    // check curves are defined over Fp
    assert(NTL::power(E0.a(), p) == E0.a());
    assert(NTL::power(E0.b(), p) == E0.b());
    assert(NTL::power(E1.a(), p) == E1.a());
    assert(NTL::power(E1.b(), p) == E1.b());

    assert(is_surface(E0) == is_surface(E1));

    if (E0 == E1)
        return O0;

    ////////////////////////////////////////////////

    std::vector<std::pair<NTL::ZZ,NTL::ZZ>> evs;
    NTL::PrimeSeq primes;
    primes.next();  // 2
    int max_k = Fexts.size();
    while (evs.size() < expvec{}.size()) {
        NTL::ZZ l(primes.next());
        if (l == avoid_prime)
            continue;
        ssize_t deg = torsionToFieldDegree(l);
        assert(deg == -1 || deg >= 1);
        if (deg == -1) continue;
        if (deg > max_k) continue;
        if (NTL::Jacobi(-O0.alg.p % l, l) != 1)
            continue;
        auto ev = NTL::SqrRootMod(-O0.alg.p % l, l);
        evs.emplace_back(l, ev);
        // std::cerr << "p=" << p << " ℓ=" << evs.back().first << " deg=" << deg << " λ=" << evs.back().second << "\n";
    }
    assert(!evs.empty());

    ////////////////////////////////////////////////
    std::unordered_map<FpE_elem, std::pair<expvec,ec>> tab0, tab1;
    tab0[E0.j_invariant()] = std::make_pair(expvec{}, E0);
    tab1[E1.j_invariant()] = std::make_pair(expvec{}, E1);

    std::queue<std::tuple<int,expvec,ec>> qu;
    qu.push({0, {}, E0});
    qu.push({1, {}, E1});

    size_t count = 0;

    while (!qu.empty()) {
        auto const [which, path, E] = std::move(qu.front()); qu.pop();

        for (size_t idx = 0; idx < evs.size(); ++idx) {
            auto [l,ev] = evs[idx];
            unsigned deg = torsionToFieldDegree(l);
            auto jt = Fexts.find(deg);
            assert(jt != Fexts.end());
            auto const &ext = jt->second;
            assert(ext.k == deg);

            while (true) {
                ecp ker {std::make_shared<const ec>(E), ext};
                {
                    FpE_push push(ext.F);
                    ecp pt {std::make_shared<const ec>(E), ext};
                    while (true) {
                        auto x = random_FpE_elem();
                        auto maybe_pt = E.lift_x(ext, x);
                        if (maybe_pt) {
                            pt = *maybe_pt;
                            break;
                        }
                    }

                    auto ord = NTL::power(p, deg) - (deg%2 ? -1 : +1);
                    assert(NTL::IsZero(ord % l));
                    assert(!(ord*pt));
                    pt = ord/l * pt;
                    if (!pt)
                        continue;
                    assert(!(l*pt));

                    auto ptp = ecp(std::make_shared<const ec>(E), ext, NTL::power(pt.get_x(), p), NTL::power(pt.get_y(), p), NTL::power(pt.get_z(), p));
                    ker = ptp + ev*pt;
                    if (!ker)
                        continue;
                #ifndef NDEBUG
                    {
                        auto kerp = ecp(std::make_shared<const ec>(E), ext, NTL::power(ker.get_x(), p), NTL::power(ker.get_y(), p), NTL::power(ker.get_z(), p));
                        assert(!(kerp - ev*ker));
                    }
                #endif
                }

                auto isog = isogeny(ker, NTL::to_long(l));

                auto const &EE = isog.get_codomain();
                auto const j = EE.j_invariant();

                auto path1 = path;
                ++path1[idx];
                auto &lookup_tab = which ? tab0 : tab1;

                ++count;
                if (auto it = lookup_tab.find(j); it != lookup_tab.end()) {
                    auto const &[path2,EEE] = it->second;
                    if (is_isomorphic(EE, EEE)) {   //TODO is this needed? the l-isogeny neighbourhood of E and its twist should look the same
                        expvec conn;
                        for (size_t i = 0; i < conn.size(); ++i)
                            conn[i] = which ? path2[i] - path1[i] : path1[i] - path2[i];

                        NTL::ZZ the_ev(0), the_norm(1);
                        for (size_t i = 0; i < conn.size(); ++i) {
                            auto c = conn[i];
                            if (!c)
                                continue;

                            auto [l,v] = evs[i];
                            if (c < 0) {
                                c = -c;
                                v = l-v;
                            }
                            assert(c > 0);
                            assert(0 <= v && v < l);
                            // std::cerr << l << " " << v << " " << c << std::endl;

                            // Hensel for the powers of l
                            auto m = l, ll = NTL::power(l, c);
                            for (decltype(c) cc = 1; cc < c; cc *= 2) {
                                m *= m;
                                v = (v - (v*v + p) * NTL::InvMod(2*v, m)) % m;
                            }
                            // std::cerr << ll << " " << v << " " << c << std::endl;
                            assert(0 <= v && v < ll);

                            // now CRT it with the rest
                            CRT(the_ev, the_norm, v, ll);
                            // std::cerr << the_norm << " " << the_ev << std::endl;
                        }

                        quat elt {{-the_ev, NTL::ZZ(0), NTL::ZZ(1), NTL::ZZ(0), NTL::ZZ(1)}, O0.alg};
                        // std::cerr << the_norm << " " << elt << std::endl;
                        auto idl = O0*the_norm + O0*elt;
                        // std::cerr << idl << std::endl;
                        // std::cerr << idl.left_order() << std::endl;
                        // std::cerr << idl.right_order() << std::endl;
                        // std::cerr << std::endl;
                        idl.reset_norm();
                        assert(idl.norm().first == the_norm && NTL::IsOne(idl.norm().second));
                        //return idl.right_order();
                        return idl._compute_order(true); //The other one is not always correct
                    }
                }

                auto &insert_tab = which ? tab1 : tab0;

                if (insert_tab.find(j) == insert_tab.end()) {
                    insert_tab[j] = std::make_pair(path1, EE);
                    qu.push({which, path1, EE});
                }
                break;
            }
        }
    }

    assert(false);  // should never happen
    __builtin_unreachable();
}


quatlat inline compute_endring(ec const &E, quatalg const &alg, std::map<unsigned,Fp2k> &Fexts, uint64_t avoid_prime=0)
{
    assert(NTL::power(E.a(), alg.p) == E.a());
    assert(NTL::power(E.b(), alg.p) == E.b());
    bool surf = is_surface(E);
    auto const [E0,O0] = starting_curve(alg, Fexts, surf);
    return transport_endring(E0, O0, E, Fexts, avoid_prime);
}


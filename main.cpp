#include "modpoly.hpp"

#include <cassert>
#include <ctime>

int main(int argc, char **argv)
{
    std::vector<std::string> args;
    for (int i = 1; i < argc; ++i)
        args.push_back(argv[i]);

    if (args.size() < 5) {
        bad_args:
        std::cerr << "usage:\n";
        std::cerr << "  " << argv[0] << " bigl j [p] [l] [j]\n";
        std::cerr << "  " << argv[0] << " bigl w [p] [l] [w]\n";
        std::cerr << "  " << argv[0] << " bigc w [p] [l] [w]\n";
        return 1;
    }
    if (args[0] != "bigc" and args[0] != "bigl")
        goto bad_args;
    if (args[1] != "j" and args[1] != "w")
        goto bad_args;

    NTL::ZZ p, l;
    NTL::ZZ_p inv;
    NTL::ZZ_pX F;
    if (!(std::istringstream(args[2]) >> p))
        goto bad_args;

    NTL::ZZ_p::init(p);
    NTL::ZZ_pX f;
    BuildIrred(f, 2);
    NTL::ZZ_pE::init(f);


    if (!(std::istringstream(args[3]) >> l))
        goto bad_args;
    if (!(std::istringstream(args[4]) >> inv))
        goto bad_args;

    if (args[0] == "bigl" && args[1] == "j")
        F = ModEvalBigLevel(p, NTL::ZZ_pE(inv), NTL::conv<long>(l));
    else if (args[0] == "bigl" && args[1] == "w")
        F = ModEvalBigLevelWeber(p, NTL::ZZ_pE(inv), NTL::conv<long>(l));
    else if (args[0] == "bigc" && args[1] == "w") {
        F = ModEvalBigCharacteristicWeber(p, inv, l);
    }
    else
        goto bad_args;

    for(int i = 0; i <= NTL::deg(F); i++)
        std::cout << "Coefficient of x^" << i << " is: " << NTL::coeff(F,i) << "\n";
    std::cout << std::flush;

    return 0;
}
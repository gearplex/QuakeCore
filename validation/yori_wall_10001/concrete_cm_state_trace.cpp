#include "quake/concrete_cm.hpp"

#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

using quake::ConcreteCM;
using quake::ConcreteCMParameters;
using quake::ConcreteCMState;

static void dump(std::ostream& out, int i, const ConcreteCMState& s) {
    out << i << ',' << std::setprecision(17)
        << s.strain << ',' << s.stress << ',' << s.tangent << ','
        << static_cast<int>(s.rule) << ',' << s.increment << ','
        << (s.has_positive_to_negative_reversal ? 1 : 0) << ','
        << (s.has_second_negative_to_positive_reversal ? 1 : 0) << ','
        << s.unloading_strain << ',' << s.unloading_stress << ','
        << s.zero_stress_strain << ',' << s.zero_stress_tangent << ','
        << s.tension_zero_strain << ',' << s.tension_peak_strain << ','
        << s.tension_peak_stress << ',' << s.tension_new_stress << ','
        << s.tension_new_tangent << ',' << s.tension_rejoin_strain << ','
        << s.positive_reversal_strain << ',' << s.positive_reversal_stress << ','
        << s.positive_zero_stress_strain << ',' << s.positive_zero_stress_tangent << ','
        << s.compression_new_stress << ',' << s.compression_new_tangent << ','
        << s.compression_rejoin_strain << ','
        << s.nested_positive_origin_strain << ','
        << s.nested_positive_target_strain << ','
        << s.nested_negative_target_strain << '\n';
}

int main(int argc, char** argv) {
    if (argc != 3 && argc != 13) {
        std::cerr << "usage: concrete_cm_state_trace strains.txt trace.csv "
                     "[fc epsc Ec rc xcrn ft et rt xcrp gap_close]\n";
        return 2;
    }
    std::ifstream in(argv[1]);
    std::ofstream out(argv[2]);
    if (!in || !out) return 2;
    std::vector<double> strains;
    double e = 0.0;
    while (in >> e) strains.push_back(e);
    if (strains.empty()) return 2;

    ConcreteCMParameters p;
    if (argc == 13) {
        p.fc = std::stod(argv[3]);
        p.epsc = std::stod(argv[4]);
        p.Ec = std::stod(argv[5]);
        p.rc = std::stod(argv[6]);
        p.xcrn = std::stod(argv[7]);
        p.ft = std::stod(argv[8]);
        p.et = std::stod(argv[9]);
        p.rt = std::stod(argv[10]);
        p.xcrp = std::stod(argv[11]);
        p.gap_close = std::stoi(argv[12]) != 0;
    } else {
        // Legacy default: story-3 confined concrete used by the original
        // wall-5/fiber-0 tangent diagnostic.
        p.fc = -7.17649;
        p.epsc = -0.00304075;
        p.Ec = 4828.707488552191;
        p.rc = 7.0;
        p.xcrn = 1.015;
        p.ft = 0.0635356;
        p.et = 2.631577918133538e-05;
        p.rt = 1.2;
        p.xcrp = 10000.0;
        p.gap_close = true;
    }

    ConcreteCM material(p);
    auto state = material.initial_state();
    const double first = strains.front();
    for (int k = 1; k <= 100; ++k) {
        state = material.trial(first * static_cast<double>(k) / 100.0, state).state;
    }

    out << "index,strain,stress,tangent,rule,increment,has_p2n,has_second_n2p,"
           "eunn,funn,espln,Epln,Te0,Teunp,Tfunp,fnewp,Enewp,esrep,"
           "er0p,fr0p,esplp,Eplp,fnewn,Enewn,esren,origin12,Tea,Teb\n";
    for (std::size_t i = 0; i < strains.size(); ++i) {
        state = material.trial(strains[i], state).state;
        dump(out, static_cast<int>(i), state);
    }
    return 0;
}

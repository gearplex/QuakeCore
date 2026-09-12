#include "quake/concrete_cm.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

using quake::ConcreteCM;
using quake::ConcreteCMEnvelope;
using quake::ConcreteCMParameters;
using quake::ConcreteCMRule;
using quake::ConcreteCMState;

namespace {

bool near(double a, double b, double rel = 1.0e-11, double abs = 1.0e-13) {
    return std::abs(a - b) <= abs + rel * std::max(std::abs(a), std::abs(b));
}

void check(bool ok, const char* message) {
    if (!ok) throw std::runtime_error(message);
}

std::vector<double> compression_protocol_prefix() {
    std::vector<double> strains{0.0};
    const double targets[] = {-0.0005, -0.0020, -0.0040};
    for (double target : targets) {
        const double start = strains.back();
        for (int i = 1; i <= 20; ++i) {
            strains.push_back(start + (target - start) * static_cast<double>(i) / 20.0);
        }
    }
    return strains;
}

void verify_oracle_prefix(const ConcreteCMEnvelope& material,
                          const std::vector<double>& expected_stress,
                          const std::vector<double>& expected_tangent,
                          double expected_work) {
    const auto strains = compression_protocol_prefix();
    const int steps[] = {0, 1, 20, 40, 41, 60};
    double work = 0.0;
    auto previous = material.compression(strains.front());
    for (std::size_t i = 1; i < strains.size(); ++i) {
        const auto current = material.compression(strains[i]);
        work += 0.5 * (previous.stress + current.stress) * (strains[i] - strains[i - 1]);
        previous = current;
    }
    for (std::size_t i = 0; i < std::size(steps); ++i) {
        const auto response = material.compression(strains[static_cast<std::size_t>(steps[i])]);
        check(near(response.stress, expected_stress[i]), "ConcreteCM compression stress disagrees with frozen OpenSees oracle");
        check(near(response.tangent, expected_tangent[i]), "ConcreteCM compression tangent disagrees with frozen OpenSees oracle");
    }
    check(near(work, expected_work, 1.0e-10, 1.0e-13), "ConcreteCM compression work disagrees with frozen OpenSees oracle");
}

ConcreteCMParameters compression_parameters(double fc, double epsc, double Ec, double xcrn) {
    // Tension fields do not participate in these Gate 4 compression-only
    // regressions. They remain neutral placeholders until the exact story-1
    // tensile parameters are recovered and tension/cyclic rules are admitted.
    return {fc, epsc, Ec, 7.0, xcrn, 1.0, 1.0 / Ec, 1.2, 10000.0, true};
}

struct OraclePoint {
    double strain{};
    double stress{};
    double tangent{};
};

ConcreteCMState commit_compression_prefix(const ConcreteCM& material) {
    auto committed = material.initial_state();
    for (double strain : compression_protocol_prefix()) {
        committed = material.trial(strain, committed).state;
    }
    return committed;
}

void verify_first_compression_reversal(const ConcreteCM& material,
                                       const std::vector<OraclePoint>& points,
                                       double work_at_step60,
                                       double expected_final_work,
                                       double first_disallowed_strain) {
    auto committed = commit_compression_prefix(material);
    check(near(committed.strain, -0.004), "ConcreteCM cyclic prefix did not reach step 60");
    check(committed.rule == ConcreteCMRule::CompressionEnvelope,
          "ConcreteCM cyclic prefix did not remain on compression envelope");

    // A trial evaluation must not mutate the committed state. Repeating the
    // first reversal point from the same committed step must be identical.
    const auto probe_a = material.trial(points.front().strain, committed);
    const auto probe_b = material.trial(points.front().strain, committed);
    check(near(probe_a.response.stress, probe_b.response.stress) &&
          near(probe_a.response.tangent, probe_b.response.tangent),
          "ConcreteCM trial evaluation mutated committed state");
    check(near(committed.strain, -0.004),
          "ConcreteCM trial evaluation changed committed strain");

    double work = work_at_step60;
    double previous_strain = committed.strain;
    double previous_stress = committed.stress;
    for (const auto& point : points) {
        const auto trial = material.trial(point.strain, committed);
        check(trial.state.rule == ConcreteCMRule::CompressionUnloading,
              "ConcreteCM first reversal left rule 3 before oracle boundary");
        check(near(trial.response.stress, point.stress),
              "ConcreteCM rule-3 stress disagrees with frozen OpenSees oracle");
        check(near(trial.response.tangent, point.tangent),
              "ConcreteCM rule-3 tangent disagrees with frozen OpenSees oracle");
        work += 0.5 * (previous_stress + trial.response.stress) *
                (point.strain - previous_strain);
        previous_strain = point.strain;
        previous_stress = trial.response.stress;
        committed = trial.state;
    }
    check(near(work, expected_final_work, 1.0e-10, 1.0e-13),
          "ConcreteCM rule-3 signed work disagrees with frozen OpenSees oracle");

    bool rejected_next_rule = false;
    try {
        (void)material.trial(first_disallowed_strain, committed);
    } catch (const std::logic_error&) {
        rejected_next_rule = true;
    }
    check(rejected_next_rule,
          "ConcreteCM admitted an unvalidated cyclic rule beyond rule 3");
}

} // namespace

int main() try {
    const auto unconfined_parameters = compression_parameters(
        -6.5, -0.002, 4595.486916530173, 1.030);
    const ConcreteCMEnvelope unconfined_envelope(unconfined_parameters);
    verify_oracle_prefix(
        unconfined_envelope,
        {0.0, -0.11453308056925031, -2.163921164843064, -6.5,
         -6.465420005089847, -5.52346369744675},
        {4595.486916530173, 4567.2031827561395, 4075.5367970685274, 0.0,
         -495.76647770689294, -495.76647770689294},
        0.019733877567269662);

    const ConcreteCM unconfined(unconfined_parameters);
    verify_first_compression_reversal(
        unconfined,
        {
            {-0.00385, -4.868521484436172, 4175.211733628288},
            {-0.0037, -4.267361077384484, 3846.494219064618},
            {-0.00355, -3.7132300094049464, 3545.2910171029707},
            {-0.0034000000000000002, -3.2029527640932196, 3260.67086343076},
            {-0.0032500000000000003, -2.7344455326361974, 2987.7784731615448},
            {-0.0031, -2.3061709629764637, 2723.882092880909},
            {-0.00295, -1.916919085854849, 2467.2392618015333},
            {-0.0028, -1.5656971917906932, 2216.6467734407483},
            {-0.00265, -1.2516667738582887, 1971.2267364300578},
            {-0.0025, -0.9741041379206195, 1730.3120674648326},
            {-0.0023499999999999997, -0.7323742027125686, 1493.3797849861303},
            {-0.0021999999999999997, -0.5259122226935871, 1260.0095852839318},
            {-0.0020499999999999997, -0.3542105508682081, 1029.8567981564838},
            {-0.0018999999999999998, -0.21680875682339806, 802.6340033899764},
            {-0.0017499999999999998, -0.11328606254884122, 578.0981117898541},
            {-0.0015999999999999999, -0.04325542930827808, 356.0410294104313},
            {-0.00145, -0.00635885151364235, 136.28274884804978},
        },
        0.019733877567269662,
        0.015001103439301009,
        -0.0013);

    const auto confined_parameters = compression_parameters(
        -8.01435, -0.00432976, 5102.805419570689, 1.015);
    const ConcreteCMEnvelope confined_envelope(confined_parameters);
    verify_oracle_prefix(
        confined_envelope,
        {0.0, -0.12640952032925176, -2.155578462184863, -5.881318217785812,
         -6.046356343523506, -7.958377184018189},
        {5102.805419570689, 5010.378571312529, 3642.3236647637405,
         1687.047305323382, 1614.318243147867, 330.2002092395013},
        0.021221688703467297);

    const ConcreteCM confined(confined_parameters);
    verify_first_compression_reversal(
        confined,
        {
            {-0.00385, -7.196074037808563, 5048.586944230898},
            {-0.0037, -6.4465506765236045, 4937.462100768252},
            {-0.00355, -5.716873322188831, 4785.3750230120595},
            {-0.0034000000000000002, -5.012669999847136, 4598.578531941104},
            {-0.0032500000000000003, -4.338845058854822, 4380.840240288793},
            {-0.0031, -3.699833105149228, 4134.777548501832},
            {-0.00295, -3.099729060896795, 3862.3561062045037},
            {-0.0028, -2.542366163834129, 3565.127598355512},
            {-0.00265, -2.031367381406704, 3244.3603779518453},
            {-0.0025, -1.5701815404253887, 2901.118299928332},
            {-0.0023499999999999997, -1.1621099249108768, 2536.3116243024037},
            {-0.0021999999999999997, -0.8103265636671555, 2150.7316127604904},
            {-0.0020499999999999997, -0.5178941438555391, 1745.0750088974564},
            {-0.0018999999999999998, -0.28777677988356487, 1319.9619332352258},
            {-0.0017499999999999998, -0.12285045305610343, 875.9493207090854},
            {-0.0015999999999999999, -0.025911682694002458, 413.5412413501326},
        },
        0.021221688703467297,
        0.013939549806617615,
        -0.00145);

    std::cout << "ConcreteCM compression envelope and first reversal match frozen OpenSees 3.8.0 oracle.\n";
    return 0;
} catch (const std::exception& e) {
    std::cerr << e.what() << '\n';
    return 1;
}

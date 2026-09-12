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

std::vector<double> protocol_through_second_compression_peak() {
    std::vector<double> strains{0.0};
    const double targets[] = {
        -0.0005, -0.0020, -0.0040, -0.0010, 0.0002, 0.0010, 0.0025,
        0.0, -0.0060, -0.0120,
    };
    for (double target : targets) {
        const double start = strains.back();
        for (int i = 1; i <= 20; ++i) {
            strains.push_back(start + (target - start) * static_cast<double>(i) / 20.0);
        }
    }
    return strains;
}

std::vector<double> compression_protocol_prefix() {
    auto all = protocol_through_second_compression_peak();
    all.resize(61);
    return all;
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
        check(near(response.stress, expected_stress[i]),
              "ConcreteCM compression stress disagrees with frozen OpenSees oracle");
        check(near(response.tangent, expected_tangent[i]),
              "ConcreteCM compression tangent disagrees with frozen OpenSees oracle");
    }
    check(near(work, expected_work, 1.0e-10, 1.0e-13),
          "ConcreteCM compression work disagrees with frozen OpenSees oracle");
}

ConcreteCMParameters yori_parameters(double fc,
                                     double epsc,
                                     double Ec,
                                     double xcrn,
                                     double ft) {
    return {fc, epsc, Ec, 7.0, xcrn, ft, 2.0 * ft / Ec,
            1.2, 10000.0, true};
}

ConcreteCMState commit_through_step(const ConcreteCM& material, int final_step) {
    const auto strains = protocol_through_second_compression_peak();
    auto committed = material.initial_state();
    for (int step = 0; step <= final_step; ++step) {
        committed = material.trial(strains[static_cast<std::size_t>(step)], committed).state;
    }
    return committed;
}

struct Checkpoint {
    int step{};
    double stress{};
    double tangent{};
    ConcreteCMRule rule{};
};

void verify_first_full_excursion(const ConcreteCM& material,
                                 const std::vector<Checkpoint>& checkpoints,
                                 double expected_work_step200,
                                 double rule9_strain,
                                 double rule9_stress,
                                 double rule9_tangent) {
    const auto strains = protocol_through_second_compression_peak();
    auto committed = material.initial_state();
    double work = 0.0;
    double previous_strain = 0.0;
    double previous_stress = 0.0;

    for (int step = 0; step <= 200; ++step) {
        const double strain = strains[static_cast<std::size_t>(step)];
        const auto trial = material.trial(strain, committed);
        for (const auto& point : checkpoints) {
            if (point.step != step) continue;
            check(near(trial.response.stress, point.stress),
                  "ConcreteCM full excursion stress disagrees with frozen OpenSees oracle");
            check(near(trial.response.tangent, point.tangent),
                  "ConcreteCM full excursion tangent disagrees with frozen OpenSees oracle");
            check(trial.state.rule == point.rule,
                  "ConcreteCM full excursion rule classification disagrees with OpenSees path");
        }
        if (step > 0) {
            work += 0.5 * (previous_stress + trial.response.stress) *
                    (strain - previous_strain);
        }
        previous_strain = strain;
        previous_stress = trial.response.stress;
        committed = trial.state;
    }
    check(near(work, expected_work_step200, 1.0e-10, 1.0e-13),
          "ConcreteCM full excursion work disagrees with frozen OpenSees oracle");

    const auto step60 = commit_through_step(material, 60);
    const auto probe_a = material.trial(rule9_strain, step60);
    const auto probe_b = material.trial(rule9_strain, step60);
    check(probe_a.state.rule == ConcreteCMRule::CompressionToTension,
          "ConcreteCM targeted rule-9 checkpoint did not enter rule 9");
    check(near(probe_a.response.stress, rule9_stress) &&
          near(probe_a.response.tangent, rule9_tangent),
          "ConcreteCM targeted rule-9 checkpoint disagrees with OpenSees 3.8.0 formulation");
    check(near(probe_a.response.stress, probe_b.response.stress) &&
          near(probe_a.response.tangent, probe_b.response.tangent),
          "ConcreteCM targeted rule-9 trial mutated committed state");

    bool rejected_reversal = false;
    try {
        (void)material.trial(-0.0115, committed); // protocol step 201
    } catch (const std::logic_error&) {
        rejected_reversal = true;
    }
    check(rejected_reversal,
          "ConcreteCM admitted unvalidated second negative-to-positive reversal behavior");
}

} // namespace

int main() try {
    const auto unconfined_parameters = yori_parameters(
        -6.5, -0.002, 4595.486916530173, 1.030, 0.0604669);
    const ConcreteCMEnvelope unconfined_envelope(unconfined_parameters);
    verify_oracle_prefix(
        unconfined_envelope,
        {0.0, -0.11453308056925031, -2.163921164843064, -6.5,
         -6.465420005089847, -5.52346369744675},
        {4595.486916530173, 4567.2031827561395, 4075.5367970685274, 0.0,
         -495.76647770689294, -495.76647770689294},
        0.019733877567269662);

    const ConcreteCM unconfined(unconfined_parameters);
    verify_first_full_excursion(
        unconfined,
        {
            {77, -0.00635885151364235, 136.28274884804978,
             ConcreteCMRule::CompressionUnloading},
            {78, 0.05062261221935397, 305.9572489863933,
             ConcreteCMRule::TensionRejoining},
            {140, 0.012642738562086833, -0.926485027524689,
             ConcreteCMRule::TensionEnvelope},
            {141, 0.008690544844429719, 26.415420513127174,
             ConcreteCMRule::TensionUnloading},
            {144, 0.000442655318296601, 19.189810520831088,
             ConcreteCMRule::TensionUnloading},
            {145, -0.0023379185693622605, 28.172758549492265,
             ConcreteCMRule::TensionToCompression},
            {160, -0.4085274336946412, 449.9597538601236,
             ConcreteCMRule::TensionToCompression},
            {173, -4.639665697801751, 1788.214999460805,
             ConcreteCMRule::TensionToCompression},
            {174, -5.110196383181005, 1065.475810265215,
             ConcreteCMRule::CompressionRejoining},
            {175, -5.255450336808381, -99.28121074510159,
             ConcreteCMRule::CompressionRejoining},
            {176, -5.126850515281235, -495.76647770689294,
             ConcreteCMRule::CompressionEnvelope},
            {200, -1.5573318757916055, -495.76647770689294,
             ConcreteCMRule::CompressionEnvelope},
        },
        0.05211791880439542,
        -0.0013340890634169552,
        0.011560052531841973,
        819.4827307488272);

    const auto confined_parameters = yori_parameters(
        -8.01435, -0.00432976, 5102.805419570689, 1.015, 0.0671422);
    const ConcreteCMEnvelope confined_envelope(confined_parameters);
    verify_oracle_prefix(
        confined_envelope,
        {0.0, -0.12640952032925176, -2.155578462184863, -5.881318217785812,
         -6.046356343523506, -7.958377184018189},
        {5102.805419570689, 5010.378571312529, 3642.3236647637405,
         1687.047305323382, 1614.318243147867, 330.2002092395013},
        0.021221688703467297);

    const ConcreteCM confined(confined_parameters);
    verify_first_full_excursion(
        confined,
        {
            {76, -0.025911682694002458, 413.5412413501326,
             ConcreteCMRule::CompressionUnloading},
            {77, 0.06575190271163192, -275.9944477614152,
             ConcreteCMRule::TensionEnvelope},
            {140, 0.013911759078056883, -0.9847616179524155,
             ConcreteCMRule::TensionEnvelope},
            {141, 0.009642230354284938, 28.527765281919528,
             ConcreteCMRule::TensionUnloading},
            {144, 0.0007380029132003692, 20.710376657319102,
             ConcreteCMRule::TensionUnloading},
            {145, -0.0020993237393562374, 28.167253842147105,
             ConcreteCMRule::TensionToCompression},
            {160, -0.5345040951139467, 625.8758431284255,
             ConcreteCMRule::TensionToCompression},
            {173, -6.983062707613081, 2834.9727606129,
             ConcreteCMRule::TensionToCompression},
            {174, -7.584188657570592, 1153.1804831201682,
             ConcreteCMRule::CompressionRejoining},
            {176, -7.945980732674475, 198.9481591447111,
             ConcreteCMRule::CompressionRejoining},
            {177, -7.961575585545971, -71.5479505326459,
             ConcreteCMRule::CompressionEnvelope},
            {200, -7.467894726870714, -71.5479505326459,
             ConcreteCMRule::CompressionEnvelope},
        },
        0.08868877526304741,
        -0.0014778637051757941,
        0.01487206449348923,
        1622.5199509108309);

    std::cout << "ConcreteCM first complete cyclic excursion matches admitted OpenSees 3.8.0 references.\n";
    return 0;
} catch (const std::exception& e) {
    std::cerr << e.what() << '\n';
    return 1;
}

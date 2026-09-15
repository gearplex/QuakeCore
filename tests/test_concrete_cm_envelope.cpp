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

std::vector<double> full_protocol() {
    std::vector<double> strains{0.0};
    const double targets[] = {
        -0.0005, -0.0020, -0.0040, -0.0010, 0.0002, 0.0010, 0.0025,
        0.0, -0.0060, -0.0120, -0.0020, 0.0,
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
    auto all = full_protocol();
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
    const auto strains = full_protocol();
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

void verify_full_protocol(const ConcreteCM& material,
                          const std::vector<Checkpoint>& checkpoints,
                          double expected_work_step240,
                          double rule9_strain,
                          double rule9_stress,
                          double rule9_tangent) {
    const auto strains = full_protocol();
    check(strains.size() == 241, "ConcreteCM frozen protocol should contain 241 points");

    auto committed = material.initial_state();
    double work = 0.0;
    double previous_strain = 0.0;
    double previous_stress = 0.0;

    for (int step = 0; step <= 240; ++step) {
        const double strain = strains[static_cast<std::size_t>(step)];
        const auto trial = material.trial(strain, committed);
        for (const auto& point : checkpoints) {
            if (point.step != step) continue;
            check(near(trial.response.stress, point.stress),
                  "ConcreteCM full protocol stress disagrees with frozen OpenSees oracle");
            check(near(trial.response.tangent, point.tangent),
                  "ConcreteCM full protocol tangent disagrees with frozen OpenSees oracle");
            check(trial.state.rule == point.rule,
                  "ConcreteCM full protocol rule classification disagrees with OpenSees path");
        }
        if (step > 0) {
            work += 0.5 * (previous_stress + trial.response.stress) *
                    (strain - previous_strain);
        }
        previous_strain = strain;
        previous_stress = trial.response.stress;
        committed = trial.state;
    }
    check(near(work, expected_work_step240, 1.0e-10, 1.0e-13),
          "ConcreteCM full protocol work disagrees with frozen OpenSees oracle");

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

    const auto step200 = commit_through_step(material, 200);
    const auto rebound_a = material.trial(strains[201], step200);
    const auto rebound_b = material.trial(strains[201], step200);
    check(rebound_a.state.has_second_negative_to_positive_reversal,
          "ConcreteCM step-201 trial did not mark the admitted second rebound");
    check(near(rebound_a.response.stress, rebound_b.response.stress) &&
          near(rebound_a.response.tangent, rebound_b.response.tangent),
          "ConcreteCM step-201 pure trial mutated committed state");
    check(step200.rule == ConcreteCMRule::CompressionEnvelope &&
          !step200.has_second_negative_to_positive_reversal,
          "ConcreteCM step-201 pure trial mutated the supplied committed state");

    bool rejected_reversal = false;
    try {
        (void)material.trial(-0.0001, committed);
    } catch (const std::logic_error&) {
        rejected_reversal = true;
    }
    check(rejected_reversal,
          "ConcreteCM admitted an unvalidated reversal after the frozen protocol");
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
    verify_full_protocol(
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
            {145, -0.0023379185693622605, 28.172758549492265,
             ConcreteCMRule::TensionToCompression},
            {174, -5.110196383181005, 1065.475810265215,
             ConcreteCMRule::CompressionRejoining},
            {176, -5.126850515281235, -495.76647770689294,
             ConcreteCMRule::CompressionEnvelope},
            {200, -1.5573318757916055, -495.76647770689294,
             ConcreteCMRule::CompressionEnvelope},
            {201, -0.8834864310566235, 935.8020713450278,
             ConcreteCMRule::CompressionUnloading},
            {206, -1.1760681025485553e-05, 2.135992113229804,
             ConcreteCMRule::CompressionUnloading},
            {207, 0.0023479104685855184, 9.60034348183496,
             ConcreteCMRule::CompressionToTension},
            {208, 0.009601581445333164, 19.414387173501645,
             ConcreteCMRule::CompressionToTension},
            {209, 0.011970462484724314, 0.7817652245539541,
             ConcreteCMRule::TensionRejoining},
            {210, 0.011891238122979231, -0.6885299881381444,
             ConcreteCMRule::TensionEnvelope},
            {220, 0.009800016755910655, -0.2645286736917171,
             ConcreteCMRule::TensionEnvelope},
            {240, 0.009332851314804602, -0.2068250932057661,
             ConcreteCMRule::TensionEnvelope},
        },
        0.050919645670333524,
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
    verify_full_protocol(
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
            {145, -0.0020993237393562374, 28.167253842147105,
             ConcreteCMRule::TensionToCompression},
            {174, -7.584188657570592, 1153.1804831201682,
             ConcreteCMRule::CompressionRejoining},
            {177, -7.961575585545971, -71.5479505326459,
             ConcreteCMRule::CompressionEnvelope},
            {200, -7.467894726870714, -71.5479505326459,
             ConcreteCMRule::CompressionEnvelope},
            {201, -5.682890704010626, 2998.9047347411415,
             ConcreteCMRule::CompressionUnloading},
            {210, -0.027426053914848758, 141.29222857050627,
             ConcreteCMRule::CompressionUnloading},
            {211, 0.00034971837576002927, 4.154970951684991,
             ConcreteCMRule::CompressionToTension},
            {212, 0.004919158972324763, 14.200654886952332,
             ConcreteCMRule::CompressionToTension},
            {213, 0.012490426557548964, 3.773005534097152,
             ConcreteCMRule::TensionRejoining},
            {214, 0.013247029450263513, 0.11029334064593144,
             ConcreteCMRule::TensionRejoining},
            {215, 0.012982384714799415, -0.7039778689077104,
             ConcreteCMRule::TensionEnvelope},
            {220, 0.011644515255360165, -0.41197942990929065,
             ConcreteCMRule::TensionEnvelope},
            {240, 0.010940821678976793, -0.30179975078300614,
             ConcreteCMRule::TensionEnvelope},
        },
        0.07702171917559245,
        -0.0014778637051757941,
        0.01487206449348923,
        1622.5199509108309);

    std::cout << "ConcreteCM full frozen cyclic protocol matches admitted OpenSees 3.8.0 references.\n";
    return 0;
} catch (const std::exception& e) {
    std::cerr << e.what() << '\n';
    return 1;
}

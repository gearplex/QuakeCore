# QuakeCore Phase 9D Architecture Prototype

Status: **architecture prototype; not an externally validated release**

Baseline source archive: uploaded `quakecore-phase9a-alpha1.zip`.
Phase 9B/9C validation findings were used as design requirements, but their later
working-tree source was not present in the uploaded files. Merge this patch into
the current development branch rather than treating this archive as the new
canonical branch tip.

## Implemented

- Explicit `ASCE41BackboneShape` policy:
  - `ResearchExtendedCDE` (legacy/default; backward compatible)
  - `StraightCE` (NIST ASCE 41-17 benchmark topology)
- Optional explicit point-C strength (`posMc`, `negMc`) so capping strength is
  independent of zero-length penalty stiffness.
- New `RCColumnASCE41Provider` with edition/provenance metadata.
- NIST benchmark adapter uses straight C->E and E=F component-failure semantics.
- ASCE 41-23 / ACI 369.1-22 adapter requires caller-resolved authorized values;
  no proprietary table coefficients are guessed or embedded.
- Cyclic degradation parameters remain separate from code-backbone provenance.
- Berkeley A/B versus C/D longitudinal reinforcement assignment corrected in the
  available validation source.
- Regression tests added for straight C->E, explicit C strength, provider
  semantics, and no embedded-code-table claim.

## Verification

- CMake Release build: PASS
- Full `quake_tests`: PASS (`All quake-core tests passed.`)

## Next gate

- Merge onto latest Phase 9C working tree.
- Add section-level RC-column inputs + authorized code-rule lookup.
- Add converged maximum-compression/two-pass parameter update.
- Add moving-surface P-M return mapping as a distinct constitutive layer.
- Run legacy vs NIST-benchmark vs production-code three-way analyses without
  EDP tuning.

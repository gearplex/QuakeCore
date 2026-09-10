# QuakeCore Phase 9D.2 Architecture Notes

Implemented in this prototype:

- factual `RCColumnSectionInput`;
- `RCColumnDemandState` for compression/tension/shear envelopes;
- abstract `RCColumnRulesResolver` and callback adapter;
- auditable `RCColumnRuleResolution` with per-parameter provenance;
- provider overloads that preserve demand/audit metadata;
- iterative axial-demand-sensitive ASCE 41-23 / ACI 369 model regeneration;
- explicit two-pass convenience workflow;
- convergence checks on both demand and generated hinge parameters;
- optional under-relaxation;
- new regression tests;
- successful `quake_tests` run;
- successful build of `ucb_validation_phase8`.

No ASCE/ACI code-table coefficients were guessed or embedded.

Next recommended increment: whole-building simultaneous column regeneration
manager, followed by Berkeley legacy/NIST/production three-way comparison.

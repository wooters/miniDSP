## Context

The miniDSP test suite (`tests/test_minidsp.c`) currently has ~200 tests covering correctness via mathematical properties. Coverage of trivial/degenerate inputs is inconsistent: some functions (`MD_dot`, `MD_energy`, `MD_impulse`) are well-covered with zero/single-element/negative cases, while others lack obvious sanity checks. All tests live in a single file following the `static int test_xxx(void)` + `RUN_TEST()` convention.

## Goals / Non-Goals

**Goals:**
- Systematically add tests for trivial inputs (zero signal, constant signal, single element) where missing.
- Add tests for identity/degenerate parameters (zero amplitude, zero delay, width=1, same-rate resample, feedback=0) where missing.
- Ensure every public API function has at least one "obvious" test case where the expected output is self-evident.
- Follow existing test conventions exactly — no new infrastructure.

**Non-Goals:**
- Not adding fuzz testing or property-based testing infrastructure.
- Not changing library code — purely additive test changes.
- Not testing `assert()`-guarded preconditions (those are structural misuse, not "bad input"). For example, we won't pass NULL pointers or N=0 where asserts protect against it.
- Not adding performance/stress tests.
- Not reorganizing or refactoring existing tests.

## Decisions

**Test only valid-but-degenerate inputs, not assert-guarded misuse.**
The library uses `assert()` for structural preconditions (NULL pointers, N=0 in most functions). Testing those would require a signal handler or fork — too complex for the benefit. Instead, we test degenerate-but-valid inputs: zero-amplitude signals, identity parameters (delay=0, width=1), single-element arrays (where N=1 is allowed), same-in/out ranges.

*Alternative considered:* Fork-and-check-crash pattern for assert tests. Rejected — platform-dependent and the asserts are simple enough to verify by inspection.

**Add tests inline in existing sections, not in a new section.**
New tests go alongside existing tests for each function, registered in the same `main()` section. This keeps related tests together.

*Alternative considered:* A separate "sanity tests" section at the end. Rejected — scatters related tests and makes it harder to see coverage per function.

**Use stack VLAs for small buffers, malloc for larger.**
Consistent with existing convention. Sanity tests use small inputs, so most will be stack-allocated.

## Risks / Trade-offs

**Risk: Test file grows larger.** → The file is already ~5000 lines. Adding ~40-60 tests adds ~500-800 lines. Still manageable as a single file. No mitigation needed.

**Risk: Some identity tests may be tautological.** → For example, testing `MD_mix` with weight_a=0, weight_b=0 (all zeros). These are still worth having as regression guards — if the implementation changes, they catch unexpected behavior on trivial inputs.

## 1. Fix MD_adjust_dblevel

- [x] 1.1 Add zero-energy guard in `src/minidsp_core.c`: if `input_energy == 0.0`, copy `in` to `out` and return early
- [x] 1.2 Add `test_adjust_dblevel_silence` to `tests/test_minidsp.c`: all-zero input produces all-zero output
- [x] 1.3 Add `test_adjust_dblevel_near_zero` to `tests/test_minidsp.c`: extremely quiet input produces finite output (no NaN/Inf)

## 2. Fix BiQuad_new

- [x] 2.1 Add `sin(omega) == 0` guard in `src/biquad.c`: if `sn == 0.0`, free the struct and return NULL
- [x] 2.2 Add `test_biquad_freq_zero` to `tests/test_minidsp.c`: `BiQuad_new` with freq=0 returns NULL
- [x] 2.3 Add `test_biquad_freq_nyquist` to `tests/test_minidsp.c`: `BiQuad_new` with freq=srate/2 returns NULL
- [x] 2.4 Add `test_biquad_freq_near_boundaries` to `tests/test_minidsp.c`: freq=1.0 and freq=srate/2-1 return non-NULL

## 3. Build and Verify

- [x] 3.1 Build and run full test suite — all tests pass including new ones

# HemePure_coupled
Coupling HemePure to another instance of HemePure.
Build dependencies before attempting to build HemePure.

## DEPENDENCIES #
1) Create `dep/build/`.
2) In `dep/build/` run `ccmake -B. -H../' or 'ccmake ..`.
3) Configure using CMake.
4) Run `make` in `dep/build/`.

## SOURCE #
1) Create `src/build/`.
2) In `src/build/` run `ccmake -B. -H../` or `ccmake ..`.
3) Configure using CMake.
4) Run `make` in `src/build/`.

## CASES #
- Run with: `run.sh`.
- Case should run with the default build options.
- Ensure that you are running with enough MPI ranks to satisfy HEMELB_READING_GROUP_SIZE.

## NOTES #
- Compilation of some dependencies is disabled:
  - CppUnit: compilation of unit tests is disabled.

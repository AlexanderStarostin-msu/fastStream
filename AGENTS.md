# Project requirements

Before making project changes, read `doc/req/nfunc.req` and the relevant
functional requirements in `doc/req/`. Treat them as project-wide constraints.

- Keep application and utility code under `src/` and `utils/`, respectively.
- For C++ classes, provide matching `.H` and `.C` files.
- For Python scripts, include top-of-file `USAGE` and `REQUIREMENTS` sections.
  Write `REQUIREMENTS` so it can guide regeneration of the script.
- Preserve existing behavior and update directly related requirements or
  documentation when implementation changes affect them.

TTF total time of flight along the trajectory
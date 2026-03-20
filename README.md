# Simulation of Hydrogen Degassing in a Vacuum Arc Degasser
Multiphase CFD simulation of hydrogen degassing in a steelmaking Vacuum Arc Degasser (VAD) using Ansys Fluent. Argon is injected through a bottom plug to drive mixing and promote hydrogen transfer from liquid steel to the gas phase under vacuum conditions.

## Physics

- Eulerian multiphase model with discrete population balance for argon bubble size distribution
- Species transport for hydrogen between liquid steel and argon phases
- Sieverts law based mass transfer implemented via UDF
- Hydrostatic argon density correction via UDF
- Standard k-ε turbulence model

## UDFs
UDFs written in C: hydrogen_transport (Sieverts law mass transfer) and argon_phase_density (hydrostatic density), combined in one file

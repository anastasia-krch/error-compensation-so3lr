# Nearest-neighbour ion-pair lifetimes

For each EMIM cation \(i\), the nearest NTf2 anion was identified in every
saved frame using the minimum-image distance between the NTf2 central N atom
and the periodic midpoint of the two EMIM ring-N atoms. If \(k_i(t)\) is the
identity of that nearest anion, the pair indicator for an origin \(t_0\) is

\[
h_i(t_0,t)=\begin{cases}
1,& k_i(t)=k_i(t_0),\\
0,& \text{otherwise}.
\end{cases}
\]

The intermittent nearest-neighbour correlation function was evaluated as

\[
C(t)=\left\langle h_i(t_0,t_0+t)\right\rangle_{i,t_0}.
\]

Reduced coordinates were saved every 1 ps. At each force field and temperature,
the included trajectories were restricted to their longest common validated
window before averaging over cations and time origins. The exact windows and
seed counts are recorded in the seed-level CSV.

Pair re-formation is therefore permitted: after an exchange, \(h_i\) returns
to unity whenever the original anion again becomes the nearest neighbour. No
distance, RDF, or coordination-number cutoff was used. The correlation is
cation-centred; every cation has exactly one nearest anion at every frame, so
\(C(0)=1\) by construction.

The ion-pair lifetime \(\tau_{1/e}\) is the first time at which \(C(t)\)
reaches \(e^{-1}\). The crossing time is obtained by linear interpolation
between the two adjacent sampled times. No parametric tail fit is used for
this observable. If a trajectory does not reach \(e^{-1}\) over its available
analysis window, its result is reported as the lower bound
\(\tau_{1/e}>t_\mathrm{max}\), rather than as a finite extrapolated lifetime.

At each temperature, values are calculated independently for all available
seeds; error bars are the sample standard deviation over seeds for temperatures
where every included seed reaches the crossing. When one or more seeds are
right-censored, seed-level values and lower bounds are retained and no finite
temperature mean is reported.

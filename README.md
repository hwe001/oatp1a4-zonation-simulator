# Oatp1a4 Zonation Simulator

A minimal compartmental transport model of the rat hepatic lobule, built to test three
competing hypotheses about the spatial (portal vein → central vein) distribution of the
uptake transporter Oatp1a4.

**Companion code/tool for the Perspective:**
> *A Minimal Transport Model Reconciles Conflicting Reports of Oatp1a4 Zonation in Rat
> Liver* — submitted to *CPT: Pharmacometrics & Systems Pharmacology* (referenced there
> under Data and Code Availability).

## Why this exists

Two published studies disagree about whether Oatp1a4 is zonated in rat liver:

- **Akanuma et al. (2019)**, *Drug Metab Pharmacokinet* 34:239–246 — IHC and in vivo
  SR-101 imaging showing Oatp1a4 sharply enriched around the central vein.
- **van de Graaf et al. (2024)**, *Acta Physiol*, doi:[10.1111/apha.14239](https://doi.org/10.1111/apha.14239)
  — argues Oatp1a4 does *not* follow a zonated expression pattern.

Both can't be right about the same protein in the same species. Rather than guess which
paper to believe, this model asks a narrower question: given the kinetic parameters
Akanuma et al. themselves measured, what would *each* hypothesis actually predict for
the experiment they ran — and does either one, followed through, contradict data already
in hand?

## What's in this repo

- **`index.html`** — the interactive browser simulator. Open it directly, or serve it
  locally (`python3 -m http.server`) and visit `index.html`. Lets you vary SR-101 dose,
  digoxin concentration, and observation time, and compare all three zonation hypotheses
  live against a discriminability table.
- **`model.py`** — a dependency-free Python reimplementation of the identical model
  (pure standard library, no numpy required), used to generate the manuscript's Figure 1.
  Run it directly: `python3 model.py`.

## The model, briefly

24 well-mixed compartments in series (portal vein → central vein), blood flowing PV→CV
at a fixed rate, each compartment exchanging with a hepatocyte pool via saturable
Michaelis-Menten influx (`Vmax(x)·C/(Km+C)`, Km = 8.2 µM from Akanuma et al.'s isolated
rat hepatocyte data) plus first-order basolateral efflux. Three candidate `Vmax(x)`
profiles — sharp pericentral, spatially uniform, and a mild continuous gradient — are
each normalized to identical total Oatp1a4 capacity, so the comparison isolates spatial
*arrangement*, not overall expression level.

**Headline finding**: because blood is progressively depleted flowing PV→CV, a spatially
uniform ("non-zonated") transporter does not, in general, predict a spatially uniform
hepatocyte signal — under high-extraction conditions it predicts the *opposite* bias
(periportal) from what Akanuma et al. observed. Separately, a digoxin dose-response
confirms transport specificity but is mathematically incapable of discriminating between
the three zonation hypotheses (the fractional suppression is identical across all three).

## Status

This is a hypothesis-generating model, not a fitted PBPK model — see the manuscript's
Limitations section. Parameters are physiologically plausible, not fitted to any dataset
beyond Akanuma et al.'s reported Km.

## License

MIT — see [LICENSE](LICENSE).

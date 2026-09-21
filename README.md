# In-silico prioritisation of phytochemical effects on mitochondrial membrane potential and metabolic-regulator binding under claim–evidence gates

**Thesis #16** (series label R5). Computational research, set out in Nile University B.Sc. chapter order for handoff.

**Author:** Kelechi Emeka Ogbonna  
**Email:** kelechiogbonna300@gmail.com  
**GitHub:** https://github.com/cloudynirvana  
**Date:** 21 September 2026

Mitochondrial-reprogramming roadmaps prioritise simulating ΔΨm and AMPK/PI3K/GLUT1 binding, but docking scores and tip-shift anecdotes are easily smuggled into kinetic Θ. How can a frozen plant-ligand library be screened so scores remain gated evidence and never become identified parameters?

Twelve toy rows are hashed and scored with a frozen linear surrogate. AutoDock Vina, Glide, GOLD and RDKit were not run. The PAINS-flagged row wins three regulator columns and stays off the priority list. Three promotion calls, including a roadmap sentence about G_tip 0.238 to 0.245, return refused, and the SHA-256 of θ does not change. A separate kinetic schedule identifies the four rates without reading the scores. The tip shift is motivation only. It is not recomputed, and it is not a result of Thesis #7 either.

This is distinct from Thesis #7, which asks whether a frozen TNBC tip ODE is identifiable under known forcings, and from Thesis #8, which treats a papaya AgNP assay as an observation channel.

This is research only. It is not a medical device, not clinical decision support, not a dose, and not a cure. No document DOI is registered.

See [DISCLAIMER.md](DISCLAIMER.md). The manuscript is [THESIS.md](THESIS.md).

## Files

| Path | Role |
| --- | --- |
| `THESIS.md` | Manuscript (Chapters 1 to 5, Vancouver citations) |
| `THESIS.pdf` | PDF built from the Markdown |
| `build_pdf.py` | Regenerates `THESIS.pdf` |
| `CITATION.cff` | Citation metadata, no document DOI |
| `DISCLAIMER.md` | Research-only boundary |
| `sim/screen.py` | Seeded screen, gates, and kinetic profile (seed 20260921) |
| `sim/library.json` | Canonical frozen library |
| `sim/results.json` | Numbers cited in Chapter Four |
| `sim/figures/` | Scores, ranks, trajectories, spectrum, profile |

## Reproduce

```bash
python3 -m pip install -r sim/requirements.txt
python3 sim/screen.py
python3 build_pdf.py
```

NumPy, SciPy and Matplotlib are required for the screen. The PDF step also needs the `markdown` and `weasyprint` packages. Regenerating the script rewrites `sim/results.json`, `sim/library.json` and `sim/figures/`. The library SHA-256 is pinned inside `sim/screen.py`. A descriptor edit that does not update the pin raises.

## Cite

Ogbonna KE. In-silico prioritisation of phytochemical effects on mitochondrial membrane potential and metabolic-regulator binding under claim–evidence gates [Internet]. Thesis #16 computational research thesis. 21 September 2026 [cited YYYY Mon DD]. Available from: https://github.com/cloudynirvana/thesis-16-mitochondrial-dpsim-phytochemical-screen

Machine-readable fields are in `CITATION.cff`. Add a document DOI there only after one exists.

Hub index, for cataloguing only: [research-theses-hub](https://github.com/cloudynirvana/research-theses-hub).

## Licence

Text and sketch code are MIT, with attribution. Computational research only.

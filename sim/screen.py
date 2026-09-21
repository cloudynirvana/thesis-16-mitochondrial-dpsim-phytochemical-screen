#!/usr/bin/env python3
"""Frozen toy-ligand screen with claim–evidence gates.

The scores are a declared linear surrogate. AutoDock Vina, Glide, GOLD,
and RDKit were not run. A score is an evidence record. The promotion
function that would copy a score into kinetic theta has no success branch.

Seed 20260921. Research only. Not a dose, a device, or a measured affinity.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares
from scipy.stats import chi2

ROOT = Path(__file__).resolve().parent
FIG = ROOT / "figures"
SEED = 20260921
SIGMA = 0.015
N_REP = 6
CHI2_95 = float(chi2.ppf(0.95, 1))

# Pinned after the canonical bytes below were hashed. A descriptor edit
# must change this constant on purpose; the script will not update it.
PINNED_LIBRARY_SHA256 = "46d685fbeccf5b7b8facf06623c6d46f34f7aa1c91fc3bed2fec0cd0fcc68a16"

TARGETS = ("psi", "ampk", "pi3k", "glut1")
TARGET_LABEL = {
    "psi": "ΔΨm proxy",
    "ampk": "AMPK",
    "pi3k": "PI3K",
    "glut1": "GLUT1",
}

# Standardization constants are declared. They are not the sample mean
# of this library, so a new row would not refit the scale.
CENTER = {
    "mw": 350.0,
    "logp": 2.5,
    "tpsa": 90.0,
    "hbd": 3.0,
    "hba": 5.0,
    "rings": 3.0,
}
SCALE = {
    "mw": 80.0,
    "logp": 1.2,
    "tpsa": 40.0,
    "hbd": 2.0,
    "hba": 2.0,
    "rings": 1.0,
}
DESC_KEYS = ("mw", "logp", "tpsa", "hbd", "hba", "rings")

# Frozen surrogate weights. More negative mimics a docking score's
# direction. The numbers are not kcal/mol and were not fit to an assay.
WEIGHTS = {
    "psi": {
        "logp": -0.90,
        "tpsa": 0.35,
        "mw": 0.10,
        "hbd": -0.15,
        "hba": 0.05,
        "rings": -0.40,
        "bias": -0.20,
    },
    "ampk": {
        "logp": -0.25,
        "tpsa": -0.55,
        "mw": 0.15,
        "hbd": -0.70,
        "hba": -0.30,
        "rings": -0.10,
        "bias": 0.05,
    },
    "pi3k": {
        "logp": -0.40,
        "tpsa": -0.20,
        "mw": -0.35,
        "hbd": -0.10,
        "hba": -0.45,
        "rings": -0.50,
        "bias": 0.10,
    },
    "glut1": {
        "logp": 0.15,
        "tpsa": -0.60,
        "mw": -0.25,
        "hbd": -0.20,
        "hba": -0.35,
        "rings": 0.10,
        "bias": -0.05,
    },
}

# Mnemonic labels only. Not a statement that the natural product was
# docked, purchased, or assayed. pains=1 is a frozen toy flag, not the
# output of a substructure filter.
LIBRARY = [
    {"id": "L01", "label": "querol", "scaffold": "flavonol", "mw": 302.0, "logp": 1.5, "tpsa": 131.0, "hbd": 5, "hba": 7, "rings": 3, "pains": 0},
    {"id": "L02", "label": "kampol", "scaffold": "flavonol", "mw": 286.0, "logp": 1.9, "tpsa": 111.0, "hbd": 4, "hba": 6, "rings": 3, "pains": 0},
    {"id": "L03", "label": "curmin", "scaffold": "diarylheptanoid", "mw": 368.0, "logp": 3.2, "tpsa": 93.0, "hbd": 2, "hba": 6, "rings": 2, "pains": 0},
    {"id": "L04", "label": "resvol", "scaffold": "stilbene", "mw": 228.0, "logp": 3.1, "tpsa": 60.0, "hbd": 3, "hba": 3, "rings": 2, "pains": 0},
    {"id": "L05", "label": "egallate", "scaffold": "flavan", "mw": 458.0, "logp": 1.2, "tpsa": 197.0, "hbd": 8, "hba": 11, "rings": 4, "pains": 1},
    {"id": "L06", "label": "berbin", "scaffold": "alkaloid", "mw": 336.0, "logp": 2.6, "tpsa": 41.0, "hbd": 0, "hba": 5, "rings": 5, "pains": 0},
    {"id": "L07", "label": "apigen", "scaffold": "flavone", "mw": 270.0, "logp": 1.7, "tpsa": 91.0, "hbd": 3, "hba": 5, "rings": 3, "pains": 0},
    {"id": "L08", "label": "luteol", "scaffold": "flavone", "mw": 286.0, "logp": 2.5, "tpsa": 112.0, "hbd": 4, "hba": 6, "rings": 3, "pains": 0},
    {"id": "L09", "label": "genist", "scaffold": "isoflavone", "mw": 270.0, "logp": 2.7, "tpsa": 87.0, "hbd": 3, "hba": 5, "rings": 3, "pains": 0},
    {"id": "L10", "label": "ellag", "scaffold": "polyphenol", "mw": 302.0, "logp": 1.1, "tpsa": 134.0, "hbd": 4, "hba": 8, "rings": 4, "pains": 0},
    {"id": "L11", "label": "piperin", "scaffold": "amide", "mw": 285.0, "logp": 3.5, "tpsa": 39.0, "hbd": 0, "hba": 4, "rings": 3, "pains": 0},
    {"id": "L12", "label": "ursol", "scaffold": "triterpene", "mw": 457.0, "logp": 7.3, "tpsa": 58.0, "hbd": 2, "hba": 3, "rings": 5, "pains": 0},
]

# Kinetic theta. These four rates are the only estimands. Scores are not
# among them. k_feed and k_drain are known constants, not estimates.
THETA_NAMES = ("k_resp", "k_leak", "k_atp", "k_use")
THETA_TRUE = {
    "k_resp": 0.80,
    "k_leak": 1.10,
    "k_atp": 0.70,
    "k_use": 0.90,
}
K_FEED = 0.50
K_DRAIN = 0.35
Y0 = np.array([0.05, 0.02, 0.75])
T_OBS = np.linspace(0.0, 12.0, 31)

# Stored as motivation metadata. Not an input. Not recomputed.
ROADMAP_ANECDOTE = {
    "id": "roadmap_G_tip",
    "text": "An unpublished mitochondrial-reprogramming roadmap mentioned G_tip moving from 0.238 to 0.245.",
    "role": "motivation_only",
    "recomputed_here": False,
    "may_enter_theta": False,
}

CLAIM_TEXT = {
    0: "row_frozen",
    1: "surrogate_score_recorded",
    2: "ranked_inside_this_surrogate",
    3: "nominate_future_observation_channel",
    4: "identified_kinetic_parameter",
}


class FreezeError(RuntimeError):
    pass


class Ledger:
    """Evidence and theta live in different fields. Promotion only appends a refusal."""

    def __init__(self, theta: dict):
        self.theta = {k: float(theta[k]) for k in THETA_NAMES}
        self.evidence = []
        self.refusals = []
        self._theta_token = theta_token(self.theta)

    def record_evidence(self, record: dict) -> None:
        if record.get("may_enter_theta", False):
            raise FreezeError("evidence record tried to set may_enter_theta true")
        if int(record.get("claim", 0)) > 3:
            raise FreezeError("claim level 4 is not a writable claim")
        record = dict(record)
        record["may_enter_theta"] = False
        self.evidence.append(record)

    def promote_score_to_theta(self, score_record: dict, theta_name: str) -> str:
        if theta_name not in THETA_NAMES:
            reason = "name_is_not_in_theta"
        else:
            reason = "surrogate_score_is_evidence_not_a_parameter"
        self.refusals.append(
            {
                "gate": "G3",
                "action": "promote_score_to_theta",
                "ligand_id": score_record.get("id"),
                "label": score_record.get("label"),
                "target": score_record.get("target"),
                "score": score_record.get("score"),
                "theta_name": theta_name,
                "status": "refused",
                "reason": reason,
            }
        )
        if theta_token(self.theta) != self._theta_token:
            raise FreezeError("theta changed inside a refusal")
        return "refused"

    def promote_anecdote_to_theta(self, anecdote: dict, theta_name: str) -> str:
        self.refusals.append(
            {
                "gate": "G4",
                "action": "promote_anecdote_to_theta",
                "anecdote_id": anecdote.get("id"),
                "theta_name": theta_name,
                "status": "refused",
                "reason": "roadmap_anecdote_is_not_a_parameter",
            }
        )
        if theta_token(self.theta) != self._theta_token:
            raise FreezeError("theta changed inside an anecdote refusal")
        return "refused"


def theta_token(theta: dict) -> str:
    payload = {k: round(float(theta[k]), 12) for k in THETA_NAMES}
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(raw).hexdigest()


def canonical_library_bytes(rows: list[dict]) -> bytes:
    payload = []
    for row in rows:
        payload.append(
            {
                "id": row["id"],
                "label": row["label"],
                "scaffold": row["scaffold"],
                "mw": row["mw"],
                "logp": row["logp"],
                "tpsa": row["tpsa"],
                "hbd": row["hbd"],
                "hba": row["hba"],
                "rings": row["rings"],
                "pains": row["pains"],
            }
        )
    return json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()


def sha256_bytes(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def library_sha256(rows: list[dict]) -> str:
    return sha256_bytes(canonical_library_bytes(rows))


def z_features(row: dict) -> dict:
    return {k: (float(row[k]) - CENTER[k]) / SCALE[k] for k in DESC_KEYS}


def surrogate_score(row: dict, target: str) -> float:
    z = z_features(row)
    w = WEIGHTS[target]
    total = float(w["bias"])
    for k in DESC_KEYS:
        total += float(w[k]) * z[k]
    return float(total)


def score_library(rows: list[dict], expected_hash: str) -> dict:
    got = library_sha256(rows)
    if expected_hash != "UNSET" and got != expected_hash:
        raise FreezeError(f"library hash {got} != pinned {expected_hash}")
    scores = {}
    for row in rows:
        scores[row["id"]] = {t: surrogate_score(row, t) for t in TARGETS}
    return scores


def ranks_from_scores(scores: dict) -> dict:
    """Rank 1 is the most negative surrogate score on that target."""
    out = {lid: {} for lid in scores}
    for target in TARGETS:
        ordered = sorted(scores, key=lambda lid: (scores[lid][target], lid))
        for rank, lid in enumerate(ordered, start=1):
            out[lid][target] = rank
    return out


def utility(rank_row: dict) -> float:
    n = len(LIBRARY)
    return float(np.mean([(n - rank_row[t] + 1) / n for t in TARGETS]))


def rhs(_t, y, theta_vec):
    psi, a, r = y
    k_resp, k_leak, k_atp, k_use = theta_vec
    dpsi = k_resp * r - k_leak * psi
    da = k_atp * psi - k_use * a
    dr = K_FEED - (K_DRAIN + k_resp) * r
    return (dpsi, da, dr)


def theta_vec(theta: dict) -> np.ndarray:
    return np.array([theta[k] for k in THETA_NAMES], dtype=float)


def simulate(theta: dict, t: np.ndarray = T_OBS) -> np.ndarray:
    sol = solve_ivp(
        lambda tt, y: rhs(tt, y, theta_vec(theta)),
        (float(t[0]), float(t[-1])),
        Y0,
        t_eval=t,
        method="LSODA",
        rtol=1e-8,
        atol=1e-10,
    )
    if not sol.success:
        raise RuntimeError(sol.message)
    # Observed channels are the membrane proxy and the ATP proxy. r is hidden.
    return np.concatenate([sol.y[0], sol.y[1]])


def hidden_redox(theta: dict, t: np.ndarray = T_OBS) -> np.ndarray:
    sol = solve_ivp(
        lambda tt, y: rhs(tt, y, theta_vec(theta)),
        (float(t[0]), float(t[-1])),
        Y0,
        t_eval=t,
        method="LSODA",
        rtol=1e-8,
        atol=1e-10,
    )
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol.y[2]


def make_data(theta: dict, rng: np.random.Generator) -> np.ndarray:
    clean = simulate(theta)
    draws = clean + rng.normal(0.0, SIGMA, size=(N_REP, clean.size))
    return draws


def residuals(theta: dict, data: np.ndarray) -> np.ndarray:
    pred = simulate(theta)
    return ((data - pred) / SIGMA).ravel()


def rss(theta: dict, data: np.ndarray) -> float:
    r = residuals(theta, data)
    return float(np.dot(r, r))


def fit_theta(data: np.ndarray, starts: list[np.ndarray]) -> dict:
    def fun(x):
        theta = {k: float(v) for k, v in zip(THETA_NAMES, x)}
        return residuals(theta, data)

    best = None
    for x0 in starts:
        fit = least_squares(fun, x0, bounds=(0.05, 4.0), xtol=1e-12, ftol=1e-12, gtol=1e-12)
        if best is None or fit.cost < best.cost:
            best = fit
    theta = {k: float(v) for k, v in zip(THETA_NAMES, best.x)}
    jac = best.jac
    fim = jac.T @ jac
    return {"theta": theta, "rss": rss(theta, data), "success": bool(best.success), "fim": fim, "nfev": int(best.nfev)}


def fisher_spectrum(fim: np.ndarray) -> dict:
    eig = np.linalg.eigvalsh(fim)
    eig = np.sort(np.real(eig))[::-1]
    smax = float(eig[0])
    rank_num = int(np.sum(eig > 1e-8 * smax))
    rank_prac = int(np.sum(eig > 1e-3 * smax))
    smin = float(eig[-1])
    cond = float(smax / smin) if smin > 0 else float("inf")
    return {
        "eigenvalues": [float(v) for v in eig],
        "rank_numerical": rank_num,
        "rank_practical_1e-3": rank_prac,
        "condition_number": cond,
        "n_parameters": int(fim.shape[0]),
    }


def expected_fim(theta: dict) -> np.ndarray:
    """Gaussian information of the mean map, multiplied by the replicate count."""
    base = simulate(theta)
    jac = np.zeros((base.size, len(THETA_NAMES)))
    for j, name in enumerate(THETA_NAMES):
        step = 1e-5 * max(1.0, abs(theta[name]))
        plus = dict(theta)
        minus = dict(theta)
        plus[name] = theta[name] + step
        minus[name] = theta[name] - step
        jac[:, j] = (simulate(plus) - simulate(minus)) / (2.0 * step)
    return (N_REP / SIGMA**2) * (jac.T @ jac)


def profile_k_resp(data: np.ndarray, mle: dict, grid: np.ndarray) -> dict:
    name_i = THETA_NAMES.index("k_resp")
    others = [n for n in THETA_NAMES if n != "k_resp"]
    rows = []
    rss_hat = mle["rss"]
    for value in grid:
        x0 = np.array([mle[n] for n in others])

        def fun(x, value=value):
            theta = {n: float(v) for n, v in zip(others, x)}
            theta["k_resp"] = float(value)
            return residuals(theta, data)

        fit = least_squares(fun, x0, bounds=(0.05, 4.0), xtol=1e-12, ftol=1e-12, gtol=1e-12)
        theta = {n: float(v) for n, v in zip(others, fit.x)}
        theta["k_resp"] = float(value)
        this = rss(theta, data)
        rows.append(
            {
                "k_resp": float(value),
                "rss": this,
                "delta_chi2": this - rss_hat,
                "k_leak": theta["k_leak"],
                "k_atp": theta["k_atp"],
                "k_use": theta["k_use"],
            }
        )
    lo, hi = interpolate_profile_bounds(rows, CHI2_95)
    return {
        "grid": rows,
        "ci95_low": lo,
        "ci95_high": hi,
        "threshold": CHI2_95,
        "profiled": "k_resp",
        "bound": "linear_interpolation_on_the_scan",
    }


def interpolate_profile_bounds(rows: list[dict], threshold: float) -> tuple[float | None, float | None]:
    lo = hi = None
    for left, right in zip(rows, rows[1:]):
        y1, y2 = left["delta_chi2"], right["delta_chi2"]
        if y1 == y2:
            continue
        if (y1 - threshold) * (y2 - threshold) > 0:
            continue
        x = left["k_resp"] + (threshold - y1) * (right["k_resp"] - left["k_resp"]) / (y2 - y1)
        if y1 > threshold >= y2 or (y1 > threshold and y2 <= threshold):
            lo = float(x)
        if y1 <= threshold < y2 or (y1 <= threshold and y2 > threshold):
            hi = float(x)
    return lo, hi


def forbidden_smuggle_for_contrast(theta: dict, score: float) -> dict:
    """A write the ledger is not allowed to perform.

    The returned vector is a contrast. It is not a mitochondrial estimate
    and it is not stored as theta.
    """
    out = {k: float(theta[k]) for k in THETA_NAMES}
    # Arbitrary map from a negative surrogate score onto k_resp.
    # Declared here so the size of the forbidden move is inspectable.
    out["k_resp"] = 0.20 + 0.50 * (-float(score))
    return out


def nll_ignoring_score(theta: dict, data: np.ndarray, score: float) -> float:
    del score
    return 0.5 * rss(theta, data)


def style_plots() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.size": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 0.6,
            "figure.dpi": 140,
            "savefig.dpi": 160,
        }
    )


def plot_scores(rows, scores) -> None:
    style_plots()
    colors = {"psi": "#1b4f72", "ampk": "#0e6655", "pi3k": "#b9770e", "glut1": "#6c3483"}
    fig, ax = plt.subplots(figsize=(8.2, 5.4))
    y = np.arange(len(rows))
    height = 0.18
    for i, target in enumerate(TARGETS):
        vals = [scores[row["id"]][target] for row in rows]
        ax.barh(
            y + (i - 1.5) * height,
            vals,
            height=height,
            color=colors[target],
            label=TARGET_LABEL[target],
        )
    labels = []
    for row in rows:
        mark = "  · PAINS flag" if row["pains"] else ""
        labels.append(f"{row['id']} {row['label']}{mark}")
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.axvline(0.0, color="#444", lw=0.6)
    ax.set_xlabel("Surrogate score (arbitrary units; more negative mimics a docking convention)")
    ax.set_title("Frozen surrogate scores. Not a docking run. Not kcal/mol.")
    ax.legend(frameon=False, ncol=4, loc="lower right")
    ax.invert_yaxis()
    fig.tight_layout()
    fig.savefig(FIG / "surrogate_scores.png")
    plt.close(fig)


def plot_ranks(rows, rank, shortlist_ids) -> None:
    style_plots()
    mat = np.array([[rank[row["id"]][t] for t in TARGETS] for row in rows], dtype=float)
    fig, ax = plt.subplots(figsize=(6.4, 5.2))
    im = ax.imshow(mat, cmap="viridis_r", vmin=1, vmax=len(rows), aspect="auto")
    ax.set_xticks(range(len(TARGETS)))
    ax.set_xticklabels([TARGET_LABEL[t] for t in TARGETS])
    ylabels = []
    for row in rows:
        tag = " shortlist" if row["id"] in shortlist_ids else ""
        if row["pains"]:
            tag = " PAINS, ineligible"
        ylabels.append(f"{row['id']} {row['label']}{tag}")
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(ylabels)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            ax.text(j, i, f"{int(mat[i, j])}", ha="center", va="center", color="white", fontsize=8)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="Rank (1 = most negative score)")
    ax.set_title("Ranks inside this surrogate. A rank is claim C2, not θ.")
    fig.tight_layout()
    fig.savefig(FIG / "rank_matrix.png")
    plt.close(fig)


def plot_trajectories(data, mle, truth) -> None:
    style_plots()
    n_t = len(T_OBS)
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.6), sharex=True)
    pred = simulate(mle)
    true = simulate(truth)
    channels = [("Membrane proxy ψ", 0), ("ATP proxy a", 1)]
    for ax, (title, block) in zip(axes, channels):
        sl = slice(block * n_t, (block + 1) * n_t)
        for rep in data:
            ax.plot(T_OBS, rep[sl], color="#b0b0b0", lw=0.7)
        ax.plot(T_OBS, true[sl], color="#1a1a1a", lw=1.6, label="generator")
        ax.plot(T_OBS, pred[sl], color="#1b4f72", lw=1.4, ls="--", label="kinetic MLE")
        ax.set_title(title)
        ax.set_xlabel("Time (arbitrary)")
        ax.set_ylabel("Dimensionless state")
    axes[0].legend(frameon=False, loc="best")
    fig.suptitle("Kinetic observations used for θ. Ligand scores are not on these axes.", y=1.02)
    fig.tight_layout()
    fig.savefig(FIG / "kinetic_trajectories.png", bbox_inches="tight")
    plt.close(fig)


def plot_spectrum(spec: dict) -> None:
    style_plots()
    eig = np.array(spec["eigenvalues"])
    fig, ax = plt.subplots(figsize=(5.6, 3.6))
    ax.bar([f"λ{i+1}" for i in range(len(eig))], eig, color="#1b4f72")
    ax.set_yscale("log")
    ax.set_ylabel("Fisher eigenvalue at the generator")
    ax.set_title(
        f"Kinetic block rank {spec['rank_numerical']} of {spec['n_parameters']}"
        f"  ·  condition {spec['condition_number']:.3g}"
    )
    fig.tight_layout()
    fig.savefig(FIG / "fisher_spectrum.png")
    plt.close(fig)


def plot_profile(profile: dict, mle_k: float, true_k: float, refused_k: float) -> None:
    style_plots()
    xs = [r["k_resp"] for r in profile["grid"]]
    ys = [r["delta_chi2"] for r in profile["grid"]]
    fig, ax = plt.subplots(figsize=(6.4, 3.8))
    ax.plot(xs, ys, color="#1b4f72", lw=1.6, label="profile of k_resp")
    ax.axhline(profile["threshold"], color="#666", ls=":", lw=1.0, label="χ² 95% (1 df)")
    ax.axvline(true_k, color="#1a1a1a", lw=0.8, label="generator")
    ax.axvline(mle_k, color="#0e6655", ls="--", lw=1.0, label="kinetic MLE")
    ax.set_xlim(min(xs), max(xs))
    ax.text(
        0.98,
        0.95,
        f"Refused write: k_resp = {refused_k:.3f}\noutside this panel, not stored in θ",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=8,
        color="#922b21",
    )
    ax.set_xlabel("k_resp on the kinetic profile")
    ax.set_ylabel("Δχ² against the kinetic MLE")
    ax.set_title("Kinetic profile of k_resp. The score was not written into θ.")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "profile_kresp_refusal.png")
    plt.close(fig)


def assign_claims(rows, rank) -> tuple[dict, list]:
    """C3 requires a clear PAINS flag and top-quartile rank on at least two targets.

    Top quartile of 12 is rank <= 3. C4 is never assigned.
    """
    claims = {}
    shortlist = []
    for row in rows:
        n_top = sum(1 for t in TARGETS if rank[row["id"]][t] <= 3)
        if row["pains"]:
            level = 2
        elif n_top >= 2:
            level = 3
            shortlist.append(row["id"])
        else:
            level = 2
        claims[row["id"]] = {
            "claim": level,
            "claim_text": CLAIM_TEXT[level],
            "n_targets_in_top_quartile": int(n_top),
            "pains": int(row["pains"]),
            "eligible": bool(level == 3),
        }
    return claims, shortlist


def main() -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    lib_hash = library_sha256(LIBRARY)
    if PINNED_LIBRARY_SHA256 != "UNSET" and lib_hash != PINNED_LIBRARY_SHA256:
        raise FreezeError(
            f"pinned library hash mismatch: got {lib_hash}, pinned {PINNED_LIBRARY_SHA256}"
        )

    mutated = [dict(row) for row in LIBRARY]
    mutated[0]["logp"] = float(mutated[0]["logp"]) + 0.01
    if library_sha256(mutated) == lib_hash:
        raise FreezeError("hash failed to see a descriptor edit")
    mutation_refused = False
    try:
        score_library(mutated, lib_hash)
    except FreezeError:
        mutation_refused = True
    if not mutation_refused:
        raise FreezeError("mutated library was scored against the pin")

    definition = {
        "library_sha256": lib_hash,
        "weights": WEIGHTS,
        "center": CENTER,
        "scale": SCALE,
    }
    definition_sha = sha256_bytes(
        json.dumps(definition, sort_keys=True, separators=(",", ":")).encode()
    )

    scores = score_library(LIBRARY, PINNED_LIBRARY_SHA256)
    rank = ranks_from_scores(scores)
    claims, shortlist = assign_claims(LIBRARY, rank)
    by_id = {row["id"]: row for row in LIBRARY}
    if any(rec["claim"] > 3 for rec in claims.values()):
        raise FreezeError("claim ceiling breached")
    if any(int(by_id[lid]["pains"]) == 1 for lid in shortlist):
        raise FreezeError("a PAINS-flagged row reached claim C3")
    utilities = {lid: utility(rank[lid]) for lid in rank}
    column_winners = {}
    for target in TARGETS:
        winner = min(scores, key=lambda lid: (scores[lid][target], lid))
        column_winners[target] = {
            "id": winner,
            "label": by_id[winner]["label"],
            "score": scores[winner][target],
            "pains": int(by_id[winner]["pains"]),
            "rank": 1,
        }

    priority = sorted(
        [lid for lid, rec in claims.items() if rec["eligible"]],
        key=lambda lid: (-utilities[lid], lid),
    )
    ineligible_order = sorted(
        [lid for lid, rec in claims.items() if not rec["eligible"]],
        key=lambda lid: (-utilities[lid], lid),
    )

    # --- kinetic theta, independent of every score ---------------------------
    rng = np.random.default_rng(SEED)
    data = make_data(THETA_TRUE, rng)
    starts = [
        np.array([1.30, 0.60, 1.10, 0.50]),
        np.array([0.40, 1.80, 0.30, 1.40]),
        np.array([0.95, 0.85, 0.55, 1.20]),
    ]
    fitted = fit_theta(data, starts)
    mle = fitted["theta"]
    fim_gen = expected_fim(THETA_TRUE)
    spec = fisher_spectrum(fim_gen)
    spec_mle = fisher_spectrum(fitted["fim"])
    cov = np.linalg.inv(fim_gen)
    se = np.sqrt(np.diag(cov))
    cramer_rao = {
        name: {
            "se": float(se_j),
            "relative_se": float(se_j / THETA_TRUE[name]),
        }
        for name, se_j in zip(THETA_NAMES, se)
    }

    grid = np.linspace(0.68, 0.84, 65)
    profile = profile_k_resp(data, {"rss": fitted["rss"], **mle}, grid)

    # Ledger starts from the kinetic MLE. Refusals must leave it byte-stable.
    ledger = Ledger(mle)
    token_before = theta_token(ledger.theta)
    for row in LIBRARY:
        for target in TARGETS:
            ledger.record_evidence(
                {
                    "kind": "surrogate_docking_score",
                    "id": row["id"],
                    "label": row["label"],
                    "target": target,
                    "score": scores[row["id"]][target],
                    "unit": "arbitrary_surrogate",
                    "engine": None,
                    "claim": 1,
                    "may_enter_theta": False,
                }
            )
        ledger.record_evidence(
            {
                "kind": "rank_claim",
                "id": row["id"],
                "label": row["label"],
                "claim": claims[row["id"]]["claim"],
                "claim_text": claims[row["id"]]["claim_text"],
                "may_enter_theta": False,
            }
        )

    ampk_winner = column_winners["ampk"]
    ampk_record = {
        "id": ampk_winner["id"],
        "label": ampk_winner["label"],
        "target": "ampk",
        "score": ampk_winner["score"],
    }
    status_1 = ledger.promote_score_to_theta(ampk_record, "k_resp")

    # Also refuse the best eligible membrane-column score, if one exists.
    psi_candidates = [lid for lid in priority] or [
        lid for lid in scores if by_id[lid]["pains"] == 0
    ]
    psi_best = min(psi_candidates, key=lambda lid: scores[lid]["psi"])
    status_2 = ledger.promote_score_to_theta(
        {"id": psi_best, "label": by_id[psi_best]["label"], "target": "psi", "score": scores[psi_best]["psi"]},
        "k_leak",
    )
    status_3 = ledger.promote_anecdote_to_theta(ROADMAP_ANECDOTE, "g_tip")
    token_after = theta_token(ledger.theta)
    if (status_1, status_2, status_3) != ("refused", "refused", "refused"):
        raise FreezeError("a promotion returned something other than refused")
    if token_before != token_after or ledger.theta != mle:
        raise FreezeError("theta moved")
    if "g_tip" in ledger.theta or any(abs(ledger.theta[k] - mle[k]) > 0 for k in THETA_NAMES):
        raise FreezeError("anecdote or score entered theta")

    # Derivative of the kinetic objective w.r.t. a score that it does not read.
    probe = float(ampk_winner["score"])
    step = 1e-3
    d_nll = (
        nll_ignoring_score(mle, data, probe + step) - nll_ignoring_score(mle, data, probe - step)
    ) / (2.0 * step)
    if d_nll != 0.0:
        raise FreezeError(f"kinetic objective depended on a score: {d_nll}")

    contrast_theta = forbidden_smuggle_for_contrast(mle, probe)
    contrast_rss = rss(contrast_theta, data)
    contrast_delta = contrast_rss - fitted["rss"]
    ci_lo, ci_hi = profile["ci95_low"], profile["ci95_high"]
    contrast_inside = bool(ci_lo is not None and ci_lo <= contrast_theta["k_resp"] <= ci_hi)
    # The contrast value must not be written back.
    if theta_token(ledger.theta) != token_before:
        raise FreezeError("contrast write leaked into the ledger")

    plot_scores(LIBRARY, scores)
    plot_ranks(LIBRARY, rank, set(shortlist))
    plot_trajectories(data, mle, THETA_TRUE)
    plot_spectrum(spec)
    plot_profile(profile, mle["k_resp"], THETA_TRUE["k_resp"], contrast_theta["k_resp"])

    score_table = []
    for row in LIBRARY:
        item = {
            "id": row["id"],
            "label": row["label"],
            "scaffold": row["scaffold"],
            "pains": int(row["pains"]),
            "utility": utilities[row["id"]],
            "claim": claims[row["id"]]["claim"],
            "claim_text": claims[row["id"]]["claim_text"],
            "n_top_quartile": claims[row["id"]]["n_targets_in_top_quartile"],
            "scores": scores[row["id"]],
            "ranks": rank[row["id"]],
        }
        score_table.append(item)

    results = {
        "seed": SEED,
        "docking_engine": None,
        "scoring": "frozen_linear_surrogate",
        "surrogate_not_docking": True,
        "library_sha256": lib_hash,
        "library_pin_status": "unset" if PINNED_LIBRARY_SHA256 == "UNSET" else "pinned",
        "definition_sha256": definition_sha,
        "n_ligands": len(LIBRARY),
        "targets": list(TARGETS),
        "claim_ceiling": 3,
        "claim_4_reachable": False,
        "shortlist_rule": "pains==0 and rank<=3 on at least 2 of 4 targets",
        "column_winners": column_winners,
        "shortlist_ids": shortlist,
        "priority_ids": priority,
        "ineligible_ids_by_utility": ineligible_order,
        "rows": score_table,
        "refusals": ledger.refusals,
        "n_refusals": len(ledger.refusals),
        "theta_names": list(THETA_NAMES),
        "theta_generating": THETA_TRUE,
        "theta_mle": mle,
        "theta_after_refusals": ledger.theta,
        "theta_sha256_before_refusals": token_before,
        "theta_sha256_after_refusals": token_after,
        "theta_unchanged": token_before == token_after,
        "known_constants": {"k_feed": K_FEED, "k_drain": K_DRAIN},
        "observation": {
            "channels": ["psi", "a"],
            "hidden": ["r"],
            "sigma": SIGMA,
            "n_replicates": N_REP,
            "n_times": int(T_OBS.size),
            "t_end": float(T_OBS[-1]),
        },
        "fit": {"rss": fitted["rss"], "success": fitted["success"], "nfev": fitted["nfev"]},
        "fisher_at_generator": spec,
        "cramer_rao_at_generator": cramer_rao,
        "fisher_at_mle_from_jacobian": {
            "eigenvalues": spec_mle["eigenvalues"],
            "rank_numerical": spec_mle["rank_numerical"],
            "rank_practical_1e-3": spec_mle["rank_practical_1e-3"],
            "condition_number": spec_mle["condition_number"],
        },
        "profile_k_resp": {
            "threshold": profile["threshold"],
            "ci95_low": profile["ci95_low"],
            "ci95_high": profile["ci95_high"],
            "grid": profile["grid"],
        },
        "dnll_dscore_at_mle": d_nll,
        "anecdote": ROADMAP_ANECDOTE,
        "contrast_not_a_result": {
            "description": "If the refused AMPK-column score were written onto k_resp by the declared forbidden map, theta would move. That vector was not stored as theta.",
            "source_ligand": ampk_winner["id"],
            "source_score": ampk_winner["score"],
            "map": "k_resp = 0.20 + 0.50 * (-score)",
            "k_resp_contrast": contrast_theta["k_resp"],
            "rss_mle": fitted["rss"],
            "rss_contrast": contrast_rss,
            "delta_rss": contrast_delta,
            "inside_profile_95": contrast_inside,
            "written_to_ledger_theta": False,
        },
    }
    (ROOT / "library.json").write_text(
        json.dumps(json.loads(canonical_library_bytes(LIBRARY).decode()), indent=2) + "\n",
        encoding="utf-8",
    )
    (ROOT / "results.json").write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")

    print("library_sha256", lib_hash)
    print("definition_sha256", definition_sha)
    print("theta_unchanged", results["theta_unchanged"])
    print("theta_sha", token_after)
    print("mle", {k: round(mle[k], 6) for k in THETA_NAMES})
    print("fisher_rank", spec["rank_numerical"], "cond", round(spec["condition_number"], 4))
    print("mle_fisher_rank", spec_mle["rank_numerical"], "cond", round(spec_mle["condition_number"], 4))
    print("profile_ci", profile["ci95_low"], profile["ci95_high"], "threshold", round(CHI2_95, 4))
    print("shortlist", shortlist)
    print("priority", priority)
    print("winners", column_winners)
    print("refusals", len(ledger.refusals), [r["reason"] for r in ledger.refusals])
    print("contrast", results["contrast_not_a_result"]["k_resp_contrast"], "delta_rss", round(contrast_delta, 3), "inside", contrast_inside)
    print("dnll_dscore", d_nll)
    print("rss", round(fitted["rss"], 4), "success", fitted["success"])


if __name__ == "__main__":
    main()

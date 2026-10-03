"""CellShield physics – fast analytical + surrogate models for a PEMFC cathode U-bend.

All functions run in milliseconds so the UI reacts instantly to sliders.
The 2D flow field is a SURROGATE (reduced-order stand-in) for the 3D Navier–Stokes
solution that run_quanscient_allsolve() is meant to produce; the UI labels it as such.

Conventions: SI units internally, *_mm fields for geometry in the UI.
"""
from __future__ import annotations

import json
import math
import time
from dataclasses import asdict, dataclass

import numpy as np

# --------------------------------------------------------------------------------------------
# Constants
# --------------------------------------------------------------------------------------------
U_CRIT = 1.0                 # m/s, below this shear cannot detach water droplets (flooding)
T_COOLANT = 65.0             # °C
K_GDL_IP = 1.5               # W/(m·K), GDL in-plane conductivity (lateral path to ribs)
K_GDL_TP = 0.5               # W/(m·K), GDL through-plane conductivity
T_GDL = 0.2e-3               # m
R_CONTACT = 1.5e-4           # m²K/W, GDL–rib thermal contact resistance
DEPTH = 0.5e-3               # m, channel depth
LEG = 10.0e-3                # m, straight leg length
E_TN = 1.25                  # V, thermoneutral voltage (LHV); 1.23 V would omit entropic heat
ETA_COMP = 0.6               # compressor isentropic efficiency
AIR_RHO, AIR_MU = 1.04, 2.03e-5   # air at ~65 °C
PLATE_K = 16.0               # W/(m·K), stainless steel metallic bipolar plate
DE_SEP = 30.0                # Dean number around which inner-wall separation appears

T_WARN, T_CRIT = 80.0, 85.0          # °C
FLOOD_WARN, FLOOD_BAD = 15.0, 30.0   # % of channel area
DP_WARN, DP_BAD = 5e3, 10e3          # Pa, full serpentine channel

PRESETS = {
    "High Load (Critical Flooding & Heat)": dict(i_acm2=2.0, v_in=3.5, v_cell=0.60),
    "Nominal Load": dict(i_acm2=1.2, v_in=2.5, v_cell=0.65),
    "Low Load / Idle": dict(i_acm2=0.5, v_in=1.2, v_cell=0.78),
}


@dataclass(frozen=True)
class Design:
    w_mm: float = 1.0          # channel width
    r_in_mm: float = 0.6       # bend inner radius
    v_in: float = 2.5          # m/s, mean inlet velocity
    i_acm2: float = 1.2        # A/cm²
    v_cell: float = 0.65       # V
    n_bends: int = 40          # U-bends per serpentine channel
    n_channels: int = 20       # parallel channels per cell
    n_cells: int = 300         # cells in stack

    @property
    def rib_mm(self) -> float:
        """Geometric lock: the rib between the two legs ends in the bend's inner radius."""
        return 2.0 * self.r_in_mm


# --------------------------------------------------------------------------------------------
# Hydraulics
# --------------------------------------------------------------------------------------------
def shah_london_fre(alpha: float) -> float:
    """Darcy f·Re for laminar flow in a rectangular duct (Shah & London 1978 fit), alpha ≤ 1."""
    a = alpha
    return 96 * (1 - 1.3553 * a + 1.9467 * a**2 - 1.7012 * a**3 + 0.9564 * a**4 - 0.2537 * a**5)


def exact_rect_fre(alpha: float, n_terms: int = 60) -> float:
    """Darcy f·Re from the exact Fourier-series solution of Poiseuille flow in a rectangle.
    Independent of the Shah & London fit – used as the validation reference."""
    n = np.arange(1, 2 * n_terms, 2, dtype=float)
    series = np.sum(np.tanh(n * np.pi / (2 * alpha)) / n**5)
    fanning = 24.0 / ((1 + alpha) ** 2 * (1 - 192 * alpha / np.pi**5 * series))
    return 4.0 * fanning


def idelchik_k180(rc_over_dh: float) -> float:
    """180° smooth-bend loss coefficient, Idelchik form K = A1·B1 with A1(180°) = 1.4."""
    r = max(rc_over_dh, 0.5)
    b1 = 0.21 / math.sqrt(r) if r >= 1.0 else 0.21 / r**2.5
    return 1.4 * b1


def hydraulics(d: Design) -> dict:
    w, rin = d.w_mm * 1e-3, d.r_in_mm * 1e-3
    rc = rin + w / 2
    dh = 2 * w * DEPTH / (w + DEPTH)
    alpha = min(w, DEPTH) / max(w, DEPTH)
    re = AIR_RHO * d.v_in * dh / AIR_MU
    de = re * math.sqrt(dh / (2 * rc))
    fre = shah_london_fre(alpha)
    lc = 2 * LEG + math.pi * rc
    dyn = 0.5 * AIR_RHO * d.v_in**2
    dp_fric = fre / re * (lc / dh) * dyn
    k_bend = idelchik_k180(rc / dh)
    dp_bend = k_bend * dyn
    dp_unit = dp_fric + dp_bend
    dp_channel = dp_unit * d.n_bends
    q_stack = d.v_in * w * DEPTH * d.n_channels * d.n_cells
    return dict(Dh_mm=dh * 1e3, alpha=alpha, Re=re, De=de, laminar=re < 2000, fRe=fre,
                Lc_mm=lc * 1e3, dyn_pa=dyn, dp_fric=dp_fric, K_bend=k_bend, dp_bend=dp_bend,
                dp_unit=dp_unit, dp_channel=dp_channel, Q_stack=q_stack,
                P_comp=dp_channel * q_stack / ETA_COMP)


# --------------------------------------------------------------------------------------------
# Thermal (analytic)
# --------------------------------------------------------------------------------------------
def thermal(d: Design) -> dict:
    q = d.i_acm2 * 1e4 * (E_TN - d.v_cell)                    # W/m²
    w = d.w_mm * 1e-3
    dt_contact = q * R_CONTACT
    dt_through = q * T_GDL / K_GDL_TP
    dt_rib = dt_contact + dt_through
    dt_lat = q * w**2 / (8 * K_GDL_IP * T_GDL)
    return dict(q=q, dT_contact=dt_contact, dT_through=dt_through, dT_rib=dt_rib,
                dT_lat=dt_lat, T_rib=T_COOLANT + dt_rib, T_max=T_COOLANT + dt_rib + dt_lat)


# --------------------------------------------------------------------------------------------
# Geometry of the unit cell (channel + inner rib + half outer ribs)
# --------------------------------------------------------------------------------------------
def geometry(d: Design, nx: int = 220) -> dict:
    """Flow: upper leg (+x) -> bend around (L, 0) -> lower leg (-x). Coordinates in mm.
    eta = 0 at the inner wall, 1 at the outer wall."""
    w, R, L = d.w_mm, d.r_in_mm, LEG * 1e3
    ro, rc = R + w, R + w / 2
    edge = ro + R                                  # half outer rib = R (periodic pitch)
    x = np.linspace(0, L + edge, nx)
    y = np.linspace(-edge, edge, max(50, int(nx * 2 * edge / (L + edge))))
    X, Y = np.meshgrid(x, y)
    r = np.hypot(X - L, Y)
    leg = X <= L
    radial = np.where(leg, np.abs(Y), r)
    chan = (radial >= R) & (radial <= ro)
    theta = np.pi / 2 - np.arctan2(Y, X - L)                        # 0 top, π bottom
    S = np.where(leg, np.where(Y > 0, X, L + np.pi * rc + (L - X)), L + theta * rc)
    eta = (radial - R) / w
    return dict(x=x, y=y, X=X, Y=Y, chan=chan, rib=~chan, S=S, eta=eta,
                s_b0=L, s_b1=L + np.pi * rc, s_tot=2 * L + np.pi * rc)


# --------------------------------------------------------------------------------------------
# Surrogate flow field
# --------------------------------------------------------------------------------------------
def _profile(eta, s, d: Design, geo: dict, hyd: dict):
    """Unnormalised depth-averaged speed shape g(eta, s)."""
    w, R, L = d.w_mm, d.r_in_mm, LEG * 1e3
    de = hyd["De"]
    n = 2 * max(1.0, w / (DEPTH * 1e3))                         # flatter profile in shallow ducts
    g = 1 - np.abs(2 * eta - 1) ** n
    # Dean skew toward the outer wall inside the bend, decaying downstream over ~2w
    wb = np.where(s < geo["s_b0"], 0.0,
                  np.where(s <= geo["s_b1"], 1.0, np.exp(-(s - geo["s_b1"]) / (2 * w))))
    g = g * (1 + min(0.6, de / 200) * wb * (2 * eta - 1))
    # Inner-wall separation bubble downstream of the bend
    onset = 1 / (1 + math.exp(-(de - DE_SEP) / 8))
    lr = min(0.8 * L, w * 0.05 * de * math.sqrt(w / R))
    bw = min(0.7, 0.5 * (w / R) ** 0.3)
    bubble = (onset * np.exp(-((s - geo["s_b1"] - 0.45 * lr) / (0.55 * lr + 1e-9)) ** 2)
              * np.exp(-(eta / bw) ** 2))
    return np.clip(g * (1 - 0.97 * bubble), 0, None)


def velocity_field(d: Design, geo: dict, hyd: dict) -> np.ndarray:
    """|u| [m/s], NaN outside the channel. Each cross-section carries the inlet flow rate,
    so a separation bubble accelerates the core – continuity is enforced by construction."""
    s_q = np.linspace(0, geo["s_tot"], 400)
    eta_q = (np.arange(64) + 0.5) / 64
    gmean = _profile(eta_q[None, :], s_q[:, None], d, geo, hyd).mean(axis=1)
    eta = np.clip(geo["eta"], 0, 1)
    U = d.v_in * _profile(eta, geo["S"], d, geo, hyd) / np.interp(geo["S"], s_q, gmean)
    return np.where(geo["chan"], U, np.nan)


def flooding_fraction(d: Design, geo: dict, hyd: dict, ns: int = 600, ne: int = 120) -> float:
    """% of channel area with |u| < U_CRIT, integrated in channel coordinates (s, eta)
    so the result is free of pixel quantisation. Bend cells are weighted by r / r_c."""
    s = (np.arange(ns) + 0.5) / ns * geo["s_tot"]
    eta = (np.arange(ne) + 0.5) / ne
    gq = _profile(eta[None, :], s[:, None], d, geo, hyd)
    U = d.v_in * gq / gq.mean(axis=1, keepdims=True)
    rc = d.r_in_mm + d.w_mm / 2
    in_bend = (s >= geo["s_b0"]) & (s <= geo["s_b1"])
    wgt = np.where(in_bend[:, None], (d.r_in_mm + eta[None, :] * d.w_mm) / rc, 1.0)
    return 100.0 * float(np.sum(wgt * (U < U_CRIT)) / np.sum(wgt))


def continuity_error(U: np.ndarray, geo: dict, d: Design) -> float:
    """Max relative flux error over leg cross-sections of the gridded field (numerical check)."""
    errs, dy = [], geo["y"][1] - geo["y"][0]
    for j in np.linspace(0.15, 0.85, 6) * np.searchsorted(geo["x"], LEG * 1e3):
        col = U[:, int(j)]
        for half in (geo["y"] > 0, geo["y"] < 0):
            flux = np.nansum(col[half]) * dy
            errs.append(abs(flux / (d.v_in * d.w_mm) - 1))
    return float(max(errs))


# --------------------------------------------------------------------------------------------
# Membrane temperature fields
# --------------------------------------------------------------------------------------------
def membrane_map(d: Design, geo: dict, th: dict) -> np.ndarray:
    """In-plane membrane temperature over the unit cell: parabolic bump under the channel,
    rib level under the ribs."""
    eta = np.clip(geo["eta"], 0, 1)
    return th["T_rib"] + np.where(geo["chan"], th["dT_lat"] * 4 * eta * (1 - eta), 0.0)


def cross_section(d: Design, th: dict, nx: int = 160, nz: int = 24) -> dict:
    """GDL cross-section over one pitch (half rib | channel | half rib).
    z = 0 membrane face, z = t_GDL rib/channel face. Analytic superposition of the
    1D through-plane drop and the 1D lateral fin bump."""
    w, h = d.w_mm, d.r_in_mm                     # half rib = R_in
    xs = np.linspace(-h, w + h, nx)
    zs = np.linspace(0, T_GDL * 1e3, nz)
    xi = np.clip(xs / w, 0, 1)
    lat = np.where((xs > 0) & (xs < w), th["dT_lat"] * 4 * xi * (1 - xi), 0.0)
    through = th["q"] * (T_GDL - zs * 1e-3) / K_GDL_TP
    T = T_COOLANT + th["dT_contact"] + through[:, None] + lat[None, :]
    return dict(x=xs, z=zs, T=T, rib_half=h, w=w)


# --------------------------------------------------------------------------------------------
# Full evaluation
# --------------------------------------------------------------------------------------------
@dataclass
class Result:
    design: Design
    hyd: dict
    th: dict
    geo: dict
    U: np.ndarray
    T_mem: np.ndarray
    flooding_risk_area_pct: float
    continuity_err: float
    compute_ms: float


def evaluate(d: Design, nx: int = 220, fields: bool = True) -> Result:
    t0 = time.perf_counter()
    hyd, th = hydraulics(d), thermal(d)
    geo = geometry(d, nx if fields else 40)
    U = velocity_field(d, geo, hyd) if fields else None
    flood = flooding_fraction(d, geo, hyd)
    T_mem = membrane_map(d, geo, th) if fields else None
    cont = continuity_error(U, geo, d) if fields else float("nan")
    return Result(d, hyd, th, geo, U, T_mem, flood, cont, (time.perf_counter() - t0) * 1e3)


def pareto_sweep(d: Design, v_values=None, r_values=None, nx: int = 80) -> list[dict]:
    v_values = np.linspace(1.0, 5.0, 17) if v_values is None else v_values
    r_values = [0.2, 0.4, 0.6, 1.0, 1.5, 2.0] if r_values is None else r_values
    pts = []
    for r in r_values:
        for v in v_values:
            dd = Design(**{**asdict(d), "v_in": float(v), "r_in_mm": float(r)})
            res = evaluate(dd, nx=nx, fields=False)
            pts.append(dict(v_in=float(v), r_in_mm=float(r), flood=res.flooding_risk_area_pct,
                            dp_unit=res.hyd["dp_unit"], P_comp=res.hyd["P_comp"]))
    return pts


def pareto_front(pts: list[dict]) -> list[dict]:
    """Non-dominated set minimising (dp_unit, flood)."""
    front, best = [], float("inf")
    for p in sorted(pts, key=lambda p: (p["dp_unit"], p["flood"])):
        if p["flood"] < best - 1e-9:
            front.append(p)
            best = p["flood"]
    return front


# --------------------------------------------------------------------------------------------
# Engineering assessment
# --------------------------------------------------------------------------------------------
def assessment(res: Result) -> list[tuple[str, str]]:
    d, h, t = res.design, res.hyd, res.th
    out = []
    if t["T_max"] > T_CRIT:
        out.append(("error", f"CRITICAL PINHOLE RISK: membrane peaks at {t['T_max']:.1f} °C under the channel centre. "
                             f"Lateral GDL conduction adds {t['dT_lat']:.1f} K and scales with w² – narrow the channel. "
                             "Raising air velocity will not fix this."))
    elif t["T_max"] >= T_WARN:
        out.append(("warning", f"Membrane at {t['T_max']:.1f} °C – above the 80 °C durability limit. "
                               "Reduce channel width or improve GDL/rib contact."))
    if res.flooding_risk_area_pct >= FLOOD_WARN:
        lvl = "error" if res.flooding_risk_area_pct >= FLOOD_BAD else "warning"
        fix = (f"increase the bend radius (R/w = {d.r_in_mm / d.w_mm:.2f}) to shorten the separation bubble"
               if d.r_in_mm / d.w_mm < 1.0 else "increase inlet velocity, at the cost of pressure drop")
        out.append((lvl, f"Flooding risk: {res.flooding_risk_area_pct:.0f}% of the channel is below "
                         f"{U_CRIT} m/s, where droplets stay attached. Try: {fix}."))
    if h["dp_channel"] >= DP_BAD:
        out.append(("error", f"Channel pressure drop {h['dp_channel'] / 1e3:.1f} kPa – compressor parasitic "
                             f"load {h['P_comp']:.0f} W. Reduce velocity or widen the duct."))
    elif h["dp_channel"] >= DP_WARN:
        out.append(("warning", f"Channel pressure drop {h['dp_channel'] / 1e3:.1f} kPa is elevated."))
    if not h["laminar"]:
        out.append(("warning", f"Re = {h['Re']:.0f} > 2000 – laminar correlations are not valid."))
    if not out:
        out.append(("success", "Design within limits: membrane below 80 °C, small flooding-risk area, "
                               "moderate pressure drop."))
    return out


# --------------------------------------------------------------------------------------------
# Quanscient Allsolve hook
# --------------------------------------------------------------------------------------------
def allsolve_payload(d: Design) -> dict:
    """Solver-agnostic specification of the 3D conjugate model. Each block maps to one stage
    of the Allsolve SDK workflow (geometry -> materials -> physics -> BCs -> mesh -> solve)."""
    th = thermal(d)
    return {
        "project": {"name": f"CellShield w{d.w_mm:.2f} R{d.r_in_mm:.2f} v{d.v_in:.2f}",
                    "description": "PEMFC cathode U-bend, conjugate laminar flow + heat transfer"},
        "geometry": {"type": "u_bend_extrusion_3d", "units": "mm",
                     "channel_width": d.w_mm, "channel_depth": DEPTH * 1e3,
                     "inner_radius": d.r_in_mm, "rib_width": d.rib_mm,
                     "leg_length": LEG * 1e3, "gdl_thickness": T_GDL * 1e3,
                     "domains": ["air_channel", "gdl", "plate"]},
        "materials": {"air_channel": {"rho": AIR_RHO, "mu": AIR_MU, "cp": 1008, "k": 0.029},
                      "gdl": {"k_in_plane": K_GDL_IP, "k_through_plane": K_GDL_TP},
                      "plate": {"k": PLATE_K}},
        "physics": [{"type": "laminar_navier_stokes", "domains": ["air_channel"], "study": "stationary"},
                    {"type": "heat_transfer", "domains": ["air_channel", "gdl", "plate"],
                     "convection_from": "laminar_navier_stokes"}],
        "boundary_conditions": {
            "inlet": {"normal_velocity": d.v_in, "temperature": T_COOLANT},
            "outlet": {"pressure": 0.0},
            "channel_walls": "no_slip",
            "plate_coolant_face": {"temperature": T_COOLANT},
            "gdl_rib_interface": {"contact_resistance": R_CONTACT},
            "membrane_face": {"heat_flux": th["q"]},
            "lateral_faces": "symmetry"},
        "mesh": {"max_element_size_mm": d.w_mm / 6, "boundary_layers": 3},
        "outputs": ["velocity_magnitude@z=depth/2", "temperature@membrane_face",
                    "pressure@channel_centerline"],
    }


def run_quanscient_allsolve(params: Design | dict) -> dict:
    """Submit the 3D model to Allsolve. Creates the cloud project via the SDK; the
    geometry/physics/solve stages are mapped from allsolve_payload() block by block.
    Never raises – returns a status dict so the UI can fall back to the surrogate."""
    d = params if isinstance(params, Design) else Design(**params)
    payload = allsolve_payload(d)
    try:
        import allsolve  # pip install allsolve  (Python ≥ 3.10), credentials in .env
    except ImportError:
        return {"status": "sdk_not_installed", "payload": payload,
                "message": "pip install allsolve and add ALLSOLVE_ACCESS_KEY / SECRET_KEY / HOST to .env"}
    try:
        client = allsolve.Client(dotenv_file=".env")
        project = client.create_project(name=payload["project"]["name"],
                                        description=json.dumps(payload["geometry"]))
        return {"status": "project_created", "project_url": client.get_url(project), "payload": payload,
                "message": "Cloud project created. Next: geometry, physics, mesh and solve stages."}
    except Exception as e:  # noqa: BLE001
        return {"status": "error", "payload": payload, "message": f"{type(e).__name__}: {e}"}

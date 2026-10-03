"""CellShield – PEMFC cathode U-bend: flooding, membrane hotspot and pressure drop.

Run:  pip install -r requirements.txt
      streamlit run app.py
"""
from dataclasses import asdict

import numpy as np
import plotly.graph_objects as go
import streamlit as st

import physics as P

st.set_page_config(page_title="CellShield", page_icon="🛡️", layout="wide")

GREEN, YELLOW, RED, GREY = "#2ea86b", "#d9a400", "#d9473f", "#5b6b78"
PLOT = dict(template="plotly_dark", paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)",
            margin=dict(l=10, r=10, t=50, b=10))

st.markdown("""
<style>
.kpi {background:#17212a;border-radius:10px;padding:14px 16px;border-left:5px solid var(--c);}
.kpi .lbl {font-size:0.85rem;color:#9fb0bd;}
.kpi .val {font-size:1.9rem;font-weight:600;color:#e3e9ee;line-height:1.2;}
.kpi .sub {font-size:0.8rem;color:#9fb0bd;}
.badge {display:inline-block;padding:2px 9px;border-radius:999px;font-size:0.75rem;font-weight:600;
        color:#0e1419;background:var(--c);margin-top:6px;}
</style>""", unsafe_allow_html=True)


def kpi(col, label, value, sub, color, badge=None):
    b = f'<div class="badge">{badge}</div>' if badge else ""
    col.markdown(f'<div class="kpi" style="--c:{color}"><div class="lbl">{label}</div>'
                 f'<div class="val">{value}</div><div class="sub">{sub}</div>{b}</div>',
                 unsafe_allow_html=True)


def level(value, warn, bad):
    return GREEN if value < warn else (YELLOW if value < bad else RED)


# --- Sidebar ----------------------------------------------------------------------------------------
def apply_preset():
    st.session_state.update(P.PRESETS[st.session_state.preset])


if "i_acm2" not in st.session_state:
    st.session_state.preset = "High Load (Critical Flooding & Heat)"
    apply_preset()

with st.sidebar:
    st.header("🛡️ CellShield")
    st.selectbox("Operating preset", list(P.PRESETS), key="preset", on_change=apply_preset)
    st.subheader("Geometry")
    w = st.slider("Channel width w [mm]", 0.5, 2.0, 1.0, 0.05)
    r_in = st.slider("Bend inner radius R_in [mm]", 0.2, 2.0, 0.6, 0.05)
    st.info(f"Linked rib width: w_rib = 2·R_in = **{2 * r_in:.2f} mm**  \n"
            f"Channel depth {P.DEPTH * 1e3:.1f} mm · legs {P.LEG * 1e3:.0f} mm")
    with st.expander("Operating point (set by preset)"):
        st.slider("Current density i [A/cm²]", 0.1, 2.5, step=0.05, key="i_acm2")
        st.slider("Inlet velocity v_in [m/s]", 0.5, 6.0, step=0.1, key="v_in")
        st.slider("Cell voltage [V]", 0.50, 0.90, step=0.01, key="v_cell")
    with st.expander("Stack scale (for ΔP and compressor power)"):
        n_bends = st.number_input("U-bends per channel", 1, 300, 40)
        n_channels = st.number_input("Parallel channels per cell", 1, 200, 20)
        n_cells = st.number_input("Cells in stack", 1, 1000, 300)

design = P.Design(w_mm=w, r_in_mm=r_in, v_in=st.session_state.v_in, i_acm2=st.session_state.i_acm2,
                  v_cell=st.session_state.v_cell, n_bends=int(n_bends), n_channels=int(n_channels),
                  n_cells=int(n_cells))


@st.cache_data(show_spinner=False)
def cached_pareto(d_dict: dict):
    d = P.Design(**d_dict)  # v_in / r_in_mm defaults are overridden by the sweep
    pts = P.pareto_sweep(d)
    return pts, P.pareto_front(pts)


res = P.evaluate(design)
h, th = res.hyd, res.th

# --- Header + KPIs -----------------------------------------------------------------------------------
st.title("CellShield · cathode U-bend optimizer")
st.caption(f"Fast surrogate + analytic models · evaluated in {res.compute_ms:.0f} ms · "
           f"Re {h['Re']:.0f} · Dean {h['De']:.0f} · {'laminar' if h['laminar'] else 'NOT laminar'} · "
           "3D verification via Quanscient Allsolve (see bottom)")

tmax = th["T_max"]
t_badge = ("CRITICAL PINHOLE RISK" if tmax > P.T_CRIT else "ABOVE 80 °C LIMIT" if tmax >= P.T_WARN
           else "WITHIN LIMIT")
c1, c2, c3, c4 = st.columns(4)
kpi(c1, "Peak membrane temperature", f"{tmax:.1f} °C",
    f"rib {th['T_rib']:.1f} °C + lateral {th['dT_lat']:.1f} K", level(tmax, P.T_WARN, P.T_CRIT + 1e-9), t_badge)
fl = res.flooding_risk_area_pct
kpi(c2, "Flooding-risk zone", f"{fl:.1f} %", f"channel area with |u| < {P.U_CRIT:.1f} m/s",
    level(fl, P.FLOOD_WARN, P.FLOOD_BAD), "HIGH" if fl >= P.FLOOD_BAD else "ELEVATED" if fl >= P.FLOOD_WARN else "LOW")
kpi(c3, "Pressure drop per U-bend", f"{h['dp_unit']:.0f} Pa",
    f"full channel {h['dp_channel'] / 1e3:.2f} kPa", level(h["dp_channel"], P.DP_WARN, P.DP_BAD))
kpi(c4, "Parasitic compressor power", f"{h['P_comp']:.0f} W",
    f"flow-field share, η = {P.ETA_COMP:.0%}", GREY)

st.write("")

# --- Tabs -----------------------------------------------------------------------------------------------
tab1, tab2, tab3 = st.tabs(["Velocity & Flooding Risk", "MEA Temperature Profile", "Pareto Trade-off"])
g = res.geo

with tab1:
    fig = go.Figure()
    fig.add_trace(go.Heatmap(x=g["x"], y=g["y"], z=np.where(g["rib"], 1.0, np.nan), showscale=False,
                             colorscale=[[0, "#2b3640"], [1, "#2b3640"]], hoverinfo="skip"))
    fig.add_trace(go.Heatmap(x=g["x"], y=g["y"], z=res.U, colorscale="Viridis", hoverongaps=False,
                             colorbar=dict(title="|u| m/s"),
                             hovertemplate="x %{x:.2f} mm<br>y %{y:.2f} mm<br>|u| %{z:.2f} m/s<extra></extra>"))
    fig.add_trace(go.Contour(x=g["x"], y=g["y"], z=res.U, showscale=False, hoverinfo="skip",
                             contours=dict(start=P.U_CRIT, end=P.U_CRIT, size=1, coloring="lines"),
                             line=dict(color="#ff5a5f", width=2.5)))
    fig.update_layout(**PLOT, height=470,
                      title=f"Depth-averaged |u| – red isoline = {P.U_CRIT} m/s (stagnation pockets inside)",
                      xaxis_title="x [mm]", yaxis=dict(title="y [mm]", scaleanchor="x"))
    st.plotly_chart(fig, width="stretch")
    st.caption("Grey = ribs (inner rib width 2·R_in). Flow enters the upper leg, turns, and exits the lower leg; "
               "the separation bubble forms on the inner wall downstream of the bend.")

with tab2:
    left, right = st.columns([3, 2])
    cs = P.cross_section(design, th)
    with left:
        fc = go.Figure(go.Contour(x=cs["x"], y=cs["z"], z=cs["T"], colorscale="Inferno",
                                  contours=dict(showlabels=True, labelfont=dict(color="white")),
                                  colorbar=dict(title="°C"),
                                  hovertemplate="x %{x:.2f} mm<br>z %{y:.3f} mm<br>T %{z:.2f} °C<extra></extra>"))
        top = cs["z"][-1]
        hgt = 0.6 * top
        for x0, x1, col, txt in [(-cs["rib_half"], 0, "#8a96a0", "rib"), (cs["w"], cs["w"] + cs["rib_half"], "#8a96a0", "rib"),
                                 (0, cs["w"], "#1f5f8b", "air channel")]:
            fc.add_shape(type="rect", x0=x0, x1=x1, y0=top, y1=top + hgt, fillcolor=col, line_width=0)
            fc.add_annotation(x=(x0 + x1) / 2, y=top + hgt / 2, text=txt, showarrow=False, font=dict(color="white"))
        imax = np.unravel_index(np.argmax(cs["T"]), cs["T"].shape)
        fc.add_trace(go.Scatter(x=[cs["x"][imax[1]]], y=[cs["z"][imax[0]]], mode="markers+text", text=["T max"],
                                textposition="top center", marker=dict(symbol="x", size=13, color="#3fb6c9"),
                                showlegend=False))
        fc.update_layout(**PLOT, height=430, title="GDL cross-section (membrane at z = 0, thickness exaggerated)",
                         xaxis_title="x across one pitch [mm]", yaxis_title="z [mm]")
        st.plotly_chart(fc, width="stretch")
    with right:
        fm = go.Figure(go.Heatmap(x=g["x"], y=g["y"], z=res.T_mem, colorscale="Inferno",
                                  colorbar=dict(title="°C"),
                                  hovertemplate="x %{x:.2f} mm<br>y %{y:.2f} mm<br>T %{z:.2f} °C<extra></extra>"))
        fm.update_layout(**PLOT, height=430, title="Membrane temperature, in-plane",
                         xaxis_title="x [mm]", yaxis=dict(title="y [mm]", scaleanchor="x"))
        st.plotly_chart(fm, width="stretch")
    st.caption(f"Heat flux q = i·(E_tn − V_cell) = {th['q'] / 1e4:.2f} W/cm². Peak sits mid-channel because heat must "
               f"travel {w / 2:.2f} mm sideways through the GDL to the ribs: ΔT_lat = q·w²/(8·k·t) = {th['dT_lat']:.1f} K. "
               "Air velocity does not appear in this budget.")

with tab3:
    sweep_key = {k: v for k, v in asdict(design).items() if k not in ("v_in", "r_in_mm")}  # sweep varies these
    pts, front = cached_pareto(sweep_key)
    fp = go.Figure()
    fp.add_trace(go.Scatter(
        x=[p["dp_unit"] for p in pts], y=[p["flood"] for p in pts], mode="markers", name="Sweep (v_in × R_in)",
        marker=dict(size=8, color=[p["r_in_mm"] for p in pts], colorscale="Tealgrn",
                    colorbar=dict(title="R_in mm"), line=dict(width=0)),
        customdata=[[p["v_in"], p["r_in_mm"], p["P_comp"]] for p in pts],
        hovertemplate="v_in %{customdata[0]:.2f} m/s<br>R_in %{customdata[1]:.1f} mm<br>ΔP %{x:.0f} Pa"
                      "<br>flooding %{y:.1f}%<br>compressor %{customdata[2]:.0f} W<extra></extra>"))
    fp.add_trace(go.Scatter(x=[p["dp_unit"] for p in front], y=[p["flood"] for p in front], mode="lines",
                            name="Pareto front", line=dict(color="#e3e9ee", width=2, dash="dash")))
    fp.add_trace(go.Scatter(x=[h["dp_unit"]], y=[fl], mode="markers", name="Current design",
                            marker=dict(symbol="star", size=20, color="#ff5a5f", line=dict(color="white", width=1))))
    fp.add_hrect(y0=P.FLOOD_BAD, y1=max(p["flood"] for p in pts) + 2, fillcolor=RED, opacity=0.08, line_width=0)
    fp.update_layout(**PLOT, height=470, title=f"Flooding risk vs pressure drop (w = {w:.2f} mm fixed)",
                     xaxis_title="ΔP per U-bend [Pa]", yaxis_title="Flooding-risk area [%]",
                     legend=dict(orientation="h", y=-0.18))
    st.plotly_chart(fp, width="stretch")
    nearest = min(front, key=lambda p: abs(p["dp_unit"] - h["dp_unit"]))
    gap = fl - nearest["flood"]
    if gap > 0.5:
        st.info(f"At the same pressure budget, the front reaches {nearest['flood']:.1f}% flooding "
                f"(v_in {nearest['v_in']:.2f} m/s, R_in {nearest['r_in_mm']:.1f} mm) – {gap:.1f} points better.")
    else:
        st.success("Current design sits on the Pareto front.")

# --- Assessment -----------------------------------------------------------------------------------------
st.subheader("Engineering assessment")
for lvl, text in P.assessment(res):
    getattr(st, lvl)(text)

# --- Validation -----------------------------------------------------------------------------------------
with st.expander("Validation – laminar friction and conservation checks"):
    exact = P.exact_rect_fre(h["alpha"])
    st.markdown(f"""
| Check | Model | Reference | Deviation |
|---|---|---|---|
| f·Re, aspect {h['alpha']:.2f} | Shah & London fit **{h['fRe']:.2f}** | exact Fourier series **{exact:.2f}** | {100 * (h['fRe'] / exact - 1):+.2f} % |
| Laminar regime | Re = **{h['Re']:.0f}** | Re < 2000 | {'✅' if h['laminar'] else '❌'} |
| Continuity of the 2D field | max flux error **{100 * res.continuity_err:.2f} %** | 0 % | leg cross-sections |
| Bend loss (Idelchik 180°) | K = **{h['K_bend']:.3f}** | R_c/D_h = {(design.r_in_mm + w / 2) / h['Dh_mm']:.2f} | {100 * h['dp_bend'] / h['dp_unit']:.1f} % of ΔP |
| Thermal budget | rib {th['dT_rib']:.2f} K + lateral {th['dT_lat']:.2f} K | T_coolant {P.T_COOLANT:.0f} °C | T_max {th['T_max']:.1f} °C |
""")
    st.caption(f"ΔP per U-bend = friction {h['dp_fric']:.1f} Pa + bend {h['dp_bend']:.1f} Pa over "
               f"centreline length {h['Lc_mm']:.1f} mm, D_h = {h['Dh_mm']:.3f} mm.")

# --- Allsolve hook -------------------------------------------------------------------------------------
with st.expander("☁️ Quanscient Allsolve – 3D conjugate verification"):
    st.write("The surrogate is fast enough for live exploration; the selected design is verified with a full 3D "
             "laminar Navier–Stokes + heat transfer model in the cloud.")
    if st.button("Submit current design to Allsolve"):
        with st.spinner("Contacting Allsolve…"):
            out = P.run_quanscient_allsolve(design)
        {"project_created": st.success, "sdk_not_installed": st.warning}.get(out["status"], st.error)(
            f"{out['status']}: {out['message']}")
        if out.get("project_url"):
            st.markdown(f"[Open project in Allsolve]({out['project_url']})")
    st.json(P.allsolve_payload(design), expanded=False)

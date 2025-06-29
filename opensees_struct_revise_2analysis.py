import openseespy.opensees as ops
import opsvis as opsv
import math
import matplotlib.pyplot as plt

import matplotlib.axes as maxes

# patch set_zlabel if missing (you already have this)
if not hasattr(maxes.Axes, "set_zlabel"):
    maxes.Axes.set_zlabel = lambda self, label: None

# patch view_init for 2D Axes
if not hasattr(maxes.Axes, "view_init"):
    maxes.Axes.view_init = lambda self, azim=None, elev=None: None

# now patch quiver to accept 6 args by dropping Z & W
_orig_quiver = maxes.Axes.quiver
def _quiver_3d_compat(self, *args, **kwargs):
    # opsvis will pass (X, Y, Z, U, V, W, **)
    if len(args) >= 6:
        X, Y, Z, U, V, W = args[:6]
        return _orig_quiver(self, X, Y, U, V, **kwargs)
    else:
        return _orig_quiver(self, *args, **kwargs)

maxes.Axes.quiver = _quiver_3d_compat


# -----------------------------------------------------------------------------
# 0) BASIC LOAD MAGNITUDES
# ---------------------------------------------------------------------------
dead_psf, live_psf = 60.0, 100.0
roof_dead_psf, roof_live_psf = 20, 20
wind_x_windward, wind_x_leeward = 50.0, 35.0
wind_y_windward, wind_y_leeward = 60.0, 40.0
wind_z_flat, wind_z_windward, wind_z_leeward = -50.0, -55.0, -45.0

# -----------------------------------------------------------------------------
# 1) USER INPUT – bay sizes & roof type/pitch
# -----------------------------------------------------------------------------
nBays_x = int(input("Number of bays in X : "))
nBays_y = int(input("Number of bays in Y : "))

print(f"Enter {nBays_x} X-bay widths [m], comma-separated:")
w_bays_x = [float(v) for v in input("  > ").split(",")]
assert len(w_bays_x) == nBays_x

print(f"Enter {nBays_y} Y-bay widths [m], comma-separated:")
w_bays_y = [float(v) for v in input("  > ").split(",")]
assert len(w_bays_y) == nBays_y

roof_type = input("Roof type [flat | monoslope | gable] : ").strip().lower()
assert roof_type in ("flat","monoslope","gable")
roof_slope_deg = 0.0 if roof_type=="flat" else float(input("Roof pitch (deg): "))
roof_slope_rad = math.radians(roof_slope_deg)

# build cumulative plan coordinates
x_coords = [0.0]
for w in w_bays_x: x_coords.append(x_coords[-1] + w)
y_coords = [0.0]
for w in w_bays_y: y_coords.append(y_coords[-1] + w)
Lx_total = x_coords[-1]

# tributary widths
def trib_y(j):
    if j == 0:       return 0.5*w_bays_y[0]
    if j == len(w_bays_y): return 0.5*w_bays_y[-1]
    return 0.5*(w_bays_y[j-1] + w_bays_y[j])

def trib_x(i):
    if i == 0:       return 0.5*w_bays_x[0]
    if i == len(w_bays_x): return 0.5*w_bays_x[-1]
    return 0.5*(w_bays_x[i-1] + w_bays_x[i])

# -----------------------------------------------------------------------------
# 2) LOAD-CASE DICTIONARY
# -----------------------------------------------------------------------------
if roof_type=="flat":
    load_cases = {
        "Dead":dead_psf,"Live":live_psf,
        "Wind_X_windward":wind_x_windward,"Wind_X_leeward":wind_x_leeward,
        "Wind_Y_windward":wind_y_windward,"Wind_Y_leeward":wind_y_leeward,
        "Wind_Z":wind_z_flat
    }
else:
    load_cases = {
        "Dead":dead_psf,"Live":live_psf,
        "Wind_X_windward":wind_x_windward,"Wind_X_leeward":wind_x_leeward,
        "Wind_Y_windward":wind_y_windward,"Wind_Y_leeward":wind_y_leeward,
        "Wind_Z_windward":wind_z_windward,"Wind_Z_leeward":wind_z_leeward
    }

# -----------------------------------------------------------------------------
# 3) FRAME GEOMETRY
# -----------------------------------------------------------------------------
E,A,Iz,Iy,J = 200e9,0.0144,553.6e-6,63.2e-6,2*553.6e-6
G = E/(2*(1+0.3))
h_story = [3,3,3,3]
z_levels = [0.0]
for h in h_story: z_levels.append(z_levels[-1]+h)
nLevels = len(z_levels)

def roof_rise(x):
    if roof_type=="flat": return 0.0
    slope = math.tan(roof_slope_rad)
    if roof_type=="monoslope": return slope*x
    xm = Lx_total/2
    return slope*x if x<=xm else slope*(Lx_total-x)

# -----------------------------------------------------------------------------
# 4) INITIALIZE MODEL & GEOM TRANSF
# -----------------------------------------------------------------------------
ops.wipe()
ops.model("basic","-ndm",3,"-ndf",6)
IDCol,IDBeamX,IDBeamY = 1,2,3
ops.geomTransf("Linear", IDCol   , 0, -1, 0)
ops.geomTransf("Linear", IDBeamX , 0, -1, 0)
ops.geomTransf("Linear", IDBeamY , -1, 0,  0)

# -----------------------------------------------------------------------------
# 5) BUILD NODES
# -----------------------------------------------------------------------------
nodes = {}
tag_id = 1
for k, z0 in enumerate(z_levels):
    for j, y in enumerate(y_coords):
        for i, x in enumerate(x_coords):
            z = z0 + (roof_rise(x) if k==nLevels-1 else 0.0)
            ops.node(tag_id, x, y, z)
            nodes[(i,j,k)] = tag_id
            tag_id += 1

nX, nY = len(x_coords), len(y_coords)

# -----------------------------------------------------------------------------
# 6) COLUMNS
# -----------------------------------------------------------------------------
eleTag = 1
col_elems = {}
for i in range(nX):
    for j in range(nY):
        for k in range(nLevels-1):
            n1 = nodes[(i,j,k)]
            n2 = nodes[(i,j,k+1)]
            ops.element("elasticBeamColumn", eleTag, n1, n2,
                        A, E, G, J, Iy, Iz, IDCol)
            col_elems[(i,j,k)] = eleTag
            eleTag += 1

# -----------------------------------------------------------------------------
# 7) X-BEAMS (subdivided)
# -----------------------------------------------------------------------------
nSeg = 6
xBeamNodes = {}
for k in range(1,nLevels):
    for j in range(nY):
        for i in range(nX-1):
            Lx = w_bays_x[i]
            beam_nodes = []
            for t in range(nSeg+1):
                if t==0:
                    nid = nodes[(i,j,k)]
                elif t==nSeg:
                    nid = nodes[(i+1,j,k)]
                else:
                    x_new = x_coords[i] + t*(Lx/nSeg)
                    y_new = y_coords[j]
                    z_new = z_levels[k] + (roof_rise(x_new) if k==nLevels-1 else 0.0)
                    nid = tag_id
                    ops.node(nid, x_new, y_new, z_new)
                    tag_id += 1
                beam_nodes.append(nid)
            xBeamNodes[(i,j,k)] = beam_nodes

            for m in range(len(beam_nodes)-1):
                ops.element("elasticBeamColumn", eleTag,
                            beam_nodes[m], beam_nodes[m+1],
                            A, E, G, J, Iy, Iz, IDBeamX)
                eleTag += 1

# -----------------------------------------------------------------------------
# 8) Y-BEAMS + DEAD/LIVE storage
# -----------------------------------------------------------------------------
load_case_loads = {"Dead":[], "Live":[]}
beam_y_ids = []

for k in range(1,nLevels-1):
    for i in range(nX-1):
        for j in range(nY-1):
            row = xBeamNodes[(i,j,k)]
            rowp= xBeamNodes[(i,j+1,k)]
            seg_len = w_bays_x[i]/nSeg
            for t in range(nSeg+1):
                n1,n2 = row[t], rowp[t]
                c1,c2 = ops.nodeCoord(n1), ops.nodeCoord(n2)
                mid  = tag_id
                ops.node(mid,
                         0.5*(c1[0]+c2[0]),
                         0.5*(c1[1]+c2[1]),
                         0.5*(c1[2]+c2[2]))
                tag_id += 1

                wx = seg_len * (0.5 if t in (0,nSeg) else 1.0)
                # first half
                for case in ("Dead","Live"):
                    load_case_loads[case].append((eleTag, load_cases[case]*wx))
                ops.element("elasticBeamColumn", eleTag, n1, mid,
                            A, E, G, J, Iy, Iz, IDBeamY)
                beam_y_ids.append(eleTag)
                eleTag += 1
                # second half
                for case in ("Dead","Live"):
                    load_case_loads[case].append((eleTag, load_cases[case]*wx))
                ops.element("elasticBeamColumn", eleTag, mid, n2,
                            A, E, G, J, Iy, Iz, IDBeamY)
                beam_y_ids.append(eleTag)
                eleTag += 1

# ---------------------------------------------------------------------------
# 9)  WIND-Z ON ROOF BEAMS
# ---------------------------------------------------------------------------
roof_beam_wind_loads = {}
k_roof = nLevels - 1
if roof_type == "flat":
    roof_beam_wind_loads["Wind_Z"] = []
else:
    roof_beam_wind_loads["Wind_Z_windward"] = []
    roof_beam_wind_loads["Wind_Z_leeward"]  = []

for i in range(nX-1):
    for j in range(nY-1):
        row  = xBeamNodes[(i, j, k_roof)]
        rowp = xBeamNodes[(i, j+1, k_roof)]
        Lx   = w_bays_x[i]
        seg  = Lx / nSeg

        for t in range(nSeg+1):
            n1, n2 = row[t], rowp[t]
            c1, c2 = ops.nodeCoord(n1), ops.nodeCoord(n2)
            midx = 0.5*(c1[0] + c2[0])
            midy = 0.5*(c1[1] + c2[1])
            midz = 0.5*(c1[2] + c2[2])

            # create the mid-node
            mid = tag_id
            ops.node(mid, midx, midy, midz)
            tag_id += 1

            # tributary width in X for this strip
            wx = seg * (0.5 if t in (0, nSeg) else 1.0)

            # determine bucket and q
            if roof_type == "flat":
                key = "Wind_Z"
            else:
                key = "Wind_Z_windward" if midx >= Lx_total/2.0 else "Wind_Z_leeward"
            q = load_cases[key] * wx

            # split into two half-beams *and* record loads on both
            tag1 = eleTag
            ops.element("elasticBeamColumn", tag1, n1, mid,
                        A, E, G, J, Iy, Iz, IDBeamY)
            roof_beam_wind_loads[key].append((tag1, q))
            eleTag += 1

            tag2 = eleTag
            ops.element("elasticBeamColumn", tag2, mid, n2,
                        A, E, G, J, Iy, Iz, IDBeamY)
            roof_beam_wind_loads[key].append((tag2, q))
            eleTag += 1
# ─────────────────────────────────────────────────────────────
# 9a) ROOF DEAD & ROOF LIVE on the top‐level Y-beams only
# ─────────────────────────────────────────────────────────────
# Prepare storage
roof_dead_loads = []
roof_live_loads = []

# Add to the load_cases dict so we can refer to them by name later
load_cases["RoofDead"] = roof_dead_psf
load_cases["RoofLive"] = roof_live_psf

# k_roof already defined as nLevels-1
for i in range(nX-1):
    for j in range(nY-1):
        # grab the node chains along X at the roof level
        chain1 = xBeamNodes[(i,   j,   k_roof)]
        chain2 = xBeamNodes[(i, j+1, k_roof)]
        Lx    = w_bays_x[i]
        seg   = Lx / nSeg

        for t in range(nSeg+1):
            n1 = chain1[t]
            n2 = chain2[t]
            c1, c2 = ops.nodeCoord(n1), ops.nodeCoord(n2)

            # create mid‐node
            mid = tag_id
            ops.node(mid,
                     0.5*(c1[0]+c2[0]),
                     0.5*(c1[1]+c2[1]),
                     0.5*(c1[2]+c2[2]))
            tag_id += 1

            # tributary width in X
            wx = seg * (0.5 if t in (0, nSeg) else 1.0)

            # first half‐beam
            tag1 = eleTag
            ops.element("elasticBeamColumn", tag1, n1, mid,
                        A, E, G, J, Iy, Iz, IDBeamY)
            # store roof dead + roof live
            roof_dead_loads.append((tag1, roof_dead_psf * wx))
            roof_live_loads.append((tag1, roof_live_psf * wx))
            eleTag += 1

            # second half‐beam
            tag2 = eleTag
            ops.element("elasticBeamColumn", tag2, mid, n2,
                        A, E, G, J, Iy, Iz, IDBeamY)
            roof_dead_loads.append((tag2, roof_dead_psf * wx))
            roof_live_loads.append((tag2, roof_live_psf * wx))
            eleTag += 1



# -----------------------------------------------------------------------------
# 10) FIX BASE – pin ALL 6 DOF!
# -----------------------------------------------------------------------------
for i in range(nX):
    for j in range(nY):
        ops.fix(nodes[(i,j,0)], 1,1,1, 0,0,0)   # ← the critical fix

# -----------------------------------------------------------------------------
# 11) QUICK VISUAL (no tags)
# -----------------------------------------------------------------------------
def plot_no_tags(title):
    opsv.plot_model(fig_wi_he=(40,22))
    for txt in plt.gca().texts:
        txt.set_visible(False)
    plt.title(title); plt.xlabel("X"); plt.ylabel("Y")
    if hasattr(plt.gca(),'set_zlabel'):
        plt.gca().set_zlabel("Z")
    plt.show()

plot_no_tags("3-D Frame (no tags)")

# -----------------------------------------------------------------------------
# 12) LOAD PATTERN FUNCTIONS
# -----------------------------------------------------------------------------
def apply_load_pattern(loads, ts, pat):
    try: ops.remove("timeSeries",ts)
    except: pass
    try: ops.remove("pattern",pat)
    except: pass
    ops.timeSeries("Constant",ts); ops.pattern("Plain",pat,ts)
    for e,q in loads:
        ops.eleLoad("-ele",e,"-type","-beamUniform",q,0,0)

def apply_load_pattern_windY(loads, ts, pat):
    try: ops.remove("timeSeries",ts)
    except: pass
    try: ops.remove("pattern",pat)
    except: pass
    ops.timeSeries("Constant",ts); ops.pattern("Plain",pat,ts)
    for e,q in loads:
        ops.eleLoad("-ele",e,"-type","-beamUniform",0,q,0)

# -----------------------------------------------------------------------------
# 13) INDIVIDUAL LOAD VISUALIZATIONS
# -----------------------------------------------------------------------------
# Dead / Live
for idx, case in enumerate(("Dead","Live")):
    ts,pat = 10+idx,10+idx
    apply_load_pattern(load_case_loads[case],ts,pat)
    plt.figure()
    opsv.plot_loads_3d(
      nep=8, fig_lbrt=(0,0,1,1),
      fmt_model_loads={"color":"k","linewidth":2,"markersize":6},
      node_supports=True, truss_node_offset=0.001,
      ax=False, sfac=3.0, fig_wi_he=(40,22), local_axes=False
    )
    plt.title(f"{case} Loads"); plt.show()
    try: ops.remove("pattern",pat); ops.remove("timeSeries",ts)
    except: pass

# ─────────────────────────────────────────────────────────────
# 13.a) ROOF DEAD & ROOF LIVE on roof beams only
# ─────────────────────────────────────────────────────────────
for idx, (case, loads) in enumerate((("RoofDead", roof_dead_loads),
                                     ("RoofLive", roof_live_loads))):
    ts = 20 + idx    # pick unused tags (e.g. 20,21)
    pat = ts
    apply_load_pattern(loads, ts, pat)
    plt.figure()
    opsv.plot_loads_3d(
        nep=8, fig_lbrt=(0,0,1,1),
        fmt_model_loads={"color":"c","linewidth":2,"markersize":6},
        node_supports=True, truss_node_offset=0.001,
        ax=False, sfac=3.0, fig_wi_he=(40,22), local_axes=False
    )
    plt.title(f"{case} Loads on Roof Beams")
    plt.show()
    try:
        ops.remove("pattern", pat)
        ops.remove("timeSeries", ts)
    except:
        pass    

# Wind X on columns
windX_leew = wind_columns_loads = {"Wind_X_leeward":[], "Wind_X_windward":[]}
for (i,j,k), e in col_elems.items():
    eff_y = trib_y(j)
    if i==0:
        wind_columns_loads["Wind_X_leeward"].append((e, load_cases["Wind_X_leeward"]*eff_y))
    if i==nX-1:
        wind_columns_loads["Wind_X_windward"].append((e, load_cases["Wind_X_windward"]*eff_y))

ts,pat = 20,20
apply_load_pattern(wind_columns_loads["Wind_X_leeward"]+wind_columns_loads["Wind_X_windward"], ts,pat)
plt.figure()
opsv.plot_loads_3d(
  nep=8, fig_lbrt=(0,0,1,1),
  fmt_model_loads={"color":"r","linewidth":2,"markersize":6},
  node_supports=True, truss_node_offset=0.001,
  ax=False, sfac=3.0, fig_wi_he=(40,22), local_axes=False
)
plt.title("Wind X on Columns"); plt.show()
try: ops.remove("pattern",pat); ops.remove("timeSeries",ts)
except: pass

# Wind Y on columns
ts,pat = 30,30
windY_loads=[]
for (i,j,k), e in col_elems.items():
    eff_x= trib_x(i)
    if j==0:
        windY_loads.append((e, load_cases["Wind_Y_windward"]*eff_x))
    if j==nY-1:
        windY_loads.append((e, load_cases["Wind_Y_leeward"]*eff_x))

apply_load_pattern_windY(windY_loads, ts,pat)
plt.figure()
opsv.plot_loads_3d(
  nep=8, fig_lbrt=(0,0,1,1),
  fmt_model_loads={"color":"m","linewidth":2,"markersize":6},
  node_supports=True, truss_node_offset=0.001,
  ax=False, sfac=3.0, fig_wi_he=(40,22), local_axes=False
)
plt.title("Wind Y on Columns"); plt.show()
try: ops.remove("pattern",pat); ops.remove("timeSeries",ts)
except: pass

# Wind Z on roof beams
ts,pat = 40,40
if roof_type=="flat":
    wz = roof_beam_wind_loads["Wind_Z"]
else:
    wz = roof_beam_wind_loads["Wind_Z_leeward"] + roof_beam_wind_loads["Wind_Z_windward"]

apply_load_pattern(wz, ts,pat)
plt.figure()
opsv.plot_loads_3d(
  nep=8, fig_lbrt=(0,0,1,1),
  fmt_model_loads={"color":"b","linewidth":2,"markersize":6},
  node_supports=True, truss_node_offset=0.001,
  ax=False, sfac=3.0, fig_wi_he=(40,22), local_axes=False
)
plt.title("Wind Z on Roof Beams"); plt.show()
try: ops.remove("pattern",pat); ops.remove("timeSeries",ts)
except: pass


# ─────────────────────────────────────────────────────────────
# 13.5) BUILD GLOBAL_LOADS DICTIONARY
# ─────────────────────────────────────────────────────────────
global_loads = {
    "Dead":    load_case_loads["Dead"],
    "Live":    load_case_loads["Live"],
    "Wind_X_leeward": wind_columns_loads["Wind_X_leeward"],
    "Wind_X_windward":wind_columns_loads["Wind_X_windward"],
    "Wind_Y_leeward": windY_loads,
    "Wind_Y_windward":windY_loads,
}

# Finally, add these two new lists into your global_loads dict:
global_loads["RoofDead"] = roof_dead_loads
global_loads["RoofLive"] = roof_live_loads

if roof_type == "flat":
    global_loads["Wind_Z"] = roof_beam_wind_loads["Wind_Z"]
else:
    global_loads["Wind_Z_leeward"]  = roof_beam_wind_loads["Wind_Z_leeward"]
    global_loads["Wind_Z_windward"] = roof_beam_wind_loads["Wind_Z_windward"]

# -----------------------------------------------------------------------------
# 14) SET UP ANALYSIS HANDLER ONCE
# -----------------------------------------------------------------------------
# ─────────────────────────────────────────────────────────────
# 14) LOAD COMBINATIONS DEFINITION
# ─────────────────────────────────────────────────────────────
if roof_type == "flat":
    load_combinations = {
        "Combo1": {"Dead":1.0, "Live":1.0, "Wind_X_leeward":0.0, "Wind_X_windward":0.0},
        "Combo2": {"Dead":1.0, "Live":0.0, "Wind_X_leeward":0.6, "Wind_X_windward":0.6, "Wind_Z":0.6},
        "Combo3": {"Dead":1.0, "Live":0.0, "Wind_Y_leeward":0.6, "Wind_Y_windward":0.6, "Wind_Z":0.6},
        "Combo4": {"Dead":0.0, "Live":0.0, "Wind_Z":1.0},
    }
else:
    load_combinations = {
        "Combo1": {"Dead":1.0, "Live":0.0, "Wind_X_leeward":1.0, "Wind_X_windward":1.0},
        "Combo2": {"Dead":1.0, "Live":0.0, "Wind_Y_leeward":1.0, "Wind_Y_windward":1.0},
        "Combo3": {"Dead":0.0, "Live":0.0, "Wind_Z_leeward":1.0, "Wind_Z_windward":1.0},
    }


ops.system("SparseGeneral")
ops.numberer("RCM")
ops.constraints("Transformation")
ops.integrator("LoadControl",1.0)
ops.algorithm("Linear")
ops.analysis("Static")

# -----------------------------------------------------------------------------
# 15) LOAD COMBINATIONS & ANALYSIS
# -----------------------------------------------------------------------------
load_combo_loads = {}        # <<-- Added
analysis_results = {}

for idx,(combo,factors) in enumerate(load_combinations.items()):
    combined = {}
    for case, f in factors.items():
        for e,q in global_loads.get(case, []):
            combined[e] = combined.get(e,0) + f*q
    loads = list(combined.items())
    load_combo_loads[combo] = loads    # <<-- Added

    ts,pat = 100+idx,100+idx
    apply_load_pattern(loads, ts,pat)
    if ops.analyze(1)==0:
        analysis_results[combo] = {n: ops.nodeDisp(n) for n in nodes.values()}
    else:
        print(f"⚠️ Analysis failed for {combo}")
    try: ops.remove("pattern",pat); ops.remove("timeSeries",ts)
    except: pass

print("Available Load Cases:", list(load_cases.keys()))
print("Available Load Combos:", list(load_combinations.keys()))

# -----------------------------------------------------------------------------
# 16) LOOP FOR DEFORMED SHAPE DISPLAY
# -----------------------------------------------------------------------------

while True:
    choice = input("Enter case/combo to plot (or q to quit): ").strip()
    if choice.lower() == "q":
        break

    ts, pat = 200, 200
    # apply the loads
    if choice in global_loads:
        if choice.startswith("Wind_Y"):
            apply_load_pattern_windY(global_loads[choice], ts, pat)
        else:
            apply_load_pattern(global_loads[choice], ts, pat)
        scale = 10.0       # small magnification for gravity‐only cases
    elif choice in load_combo_loads:
        apply_load_pattern(load_combo_loads[choice], ts, pat)
        scale = 250.0      # large magnification to see wind effects
    else:
        print("Invalid choice")
        continue

    # run the analysis
    if ops.analyze(1) != 0:
        print("⚠️ Convergence failed")

    # 1) Deformed shape
    plt.figure(figsize=(30,20))
    opsv.plot_defo(scale, 9, az_el=(-68,39), endDispFlag=0, fig_wi_he=(30,20))
    plt.title(f"Deformed Shape – {choice} (scale×{scale:g})")
    plt.show()

    # 2) Bending Moment Diagram
    mFac = 5e-6
    plt.figure(figsize=(50,20))
    opsv.section_force_diagram_2d(
        'M', mFac, fig_wi_he=(50,20),
        fmt_secforce1={'color':'green'},
        fmt_secforce2={'color':'green'}
    )
    plt.title(f"Bending Moment Diagram – {choice}")
    plt.xlabel("Distance (m)")
    plt.ylabel("Moment (N·m)")
    plt.grid(True)
    plt.show()

    # 3) Shear Force Diagram
    vFac = 15e-6
    plt.figure(figsize=(50,20))
    opsv.section_force_diagram_2d(
        'V', vFac, fig_wi_he=(50,20),
        fmt_secforce1={'color':'red'},
        fmt_secforce2={'color':'red'}
    )
    plt.title(f"Shear Force Diagram – {choice}")
    plt.xlabel("Distance (m)")
    plt.ylabel("Shear (N)")
    plt.grid(True)
    plt.show()

    # 4) Extract & print the overall maxima
    maxM, maxV = 0.0, 0.0
    maxM_ele, maxV_ele = None, None
    for e in range(1, eleTag):
        try:
            f = ops.eleForce(e)
            Mval = max(abs(f[4]), abs(f[5]))
            Vval = max(abs(f[1]), abs(f[2]))
            if Mval > maxM:
                maxM, maxM_ele = Mval, e
            if Vval > maxV:
                maxV, maxV_ele = Vval, e
        except:
            continue

    print(f"🔹 {choice} ⇒ max bending M = {maxM:.3f} at element #{maxM_ele}")
    print(f"🔸 {choice} ⇒ max shear  V = {maxV:.3f} at element #{maxV_ele}")
    print("---------------------------------------------------------\n")

    # 4a) Prompt for a specific element
    while True:
        sel = input("Enter an element tag to query its M, V and deflection (or 'n' to skip): ").strip()
        if sel.lower()=='n':
            break
        try:
            e_sel = int(sel)
        except ValueError:
            print("  → please enter a valid integer element tag or 'n'")
            continue

        # fetch section‐forces for that element
        try:
            fsel = ops.eleForce(e_sel)
        except Exception:
            print(f"  → element {e_sel} not found.")
            continue

        M_sel = max(abs(fsel[4]), abs(fsel[5]))
        V_sel = max(abs(fsel[1]), abs(fsel[2]))
        # fetch its end‐node tags
        n1, n2 = ops.eleNodes(e_sel)
        # vertical deflections of those nodes
        d1z = abs(ops.nodeDisp(n1)[2])
        d2z = abs(ops.nodeDisp(n2)[2])
        defl_sel = max(d1z, d2z)

        print(f"  • Element {e_sel}:   M = {M_sel:.3f},   V = {V_sel:.3f},   max Δz = {defl_sel:.6f}")
        break

    # clean up for next choice
    try:
        ops.remove("pattern", pat)
        ops.remove("timeSeries", ts)
    except:
        pass

print("✅ Done.")


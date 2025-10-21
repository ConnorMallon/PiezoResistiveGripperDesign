import gmsh
import math

# 1) Initialize Gmsh
gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 1)

# 2) Domain dimensions
Lx, Ly = 0.15, 0.07


factor = 1

# 3) Mesh sizes
lc_outer = 0.003 * factor    # coarse background
lc_ref   = 0.003 * factor # refined near features
lc_extra_ref = lc_outer #0.0005 * factor # extra refinement around specific points

# 4) Define rectangle corners
p1 = gmsh.model.occ.addPoint(0.0,  0.0, 0.0, lc_outer)
p2 = gmsh.model.occ.addPoint(Lx,   0.0, 0.0, lc_outer)
p3 = gmsh.model.occ.addPoint(Lx,   Ly,  0.0, lc_outer)
p4 = gmsh.model.occ.addPoint(0.0,  Ly,  0.0, lc_outer)

# 5) Bottom boundary splits at x = [0, .01, .035, .045, .075, .085, Lx]
bx = [0.0, 0.01, 0.035, 0.045, 0.075, 0.085, Lx]
bottom_pts = []
for x in bx:
    lc = lc_ref if abs(x - 0.075) < 1e-8 else lc_outer
    bottom_pts.append(gmsh.model.occ.addPoint(x, 0.0, 0.0, lc))
bottom_lines = [
    gmsh.model.occ.addLine(bottom_pts[i], bottom_pts[i+1])
    for i in range(len(bottom_pts)-1)
]

# 6) Top boundary splits at x = [0, .12, .125, .145, Lx]
tx = [0.0, 0.12, 0.125, 0.14, Lx]
top_pts = []
for x in tx:
    lc = lc_ref if abs(x - 0.14) < 1e-8 else lc_outer
    top_pts.append(gmsh.model.occ.addPoint(x, Ly, 0.0, lc))
top_lines = [
    gmsh.model.occ.addLine(top_pts[i], top_pts[i+1])
    for i in range(len(top_pts)-1)
]

# 7) Left and right boundary lines
right_line = gmsh.model.occ.addLine(bottom_pts[-1], top_pts[-1])
left_line  = gmsh.model.occ.addLine(top_pts[0],    bottom_pts[0])

# 8) Create the main surface
boundary = [*bottom_lines, right_line, *reversed(top_lines), left_line]
loop     = gmsh.model.occ.addCurveLoop(boundary)
surf     = gmsh.model.occ.addPlaneSurface([loop])

# 9) Horizontal splits for D1 and D2 zones
y_D1 = 0.0046666667
d1A = gmsh.model.occ.addPoint(0.0,         y_D1, 0.0, lc_outer)
d1B = gmsh.model.occ.addPoint(Lx,          y_D1, 0.0, lc_outer)
d1_line = gmsh.model.occ.addLine(d1A, d1B)

d2A = gmsh.model.occ.addPoint(0.0,         0.06, 0.0, lc_outer)
d2B = gmsh.model.occ.addPoint(Lx,          0.06, 0.0, lc_outer)
d2_line = gmsh.model.occ.addLine(d2A, d2B)

# 10) Vertical splits inside those bands
vD1A = gmsh.model.occ.addPoint(0.085, 0.0,  0.0, lc_outer)
vD1B = gmsh.model.occ.addPoint(0.085, y_D1, 0.0, lc_outer)
vD2A = gmsh.model.occ.addPoint(0.12,  0.06,0.0, lc_outer)
vD2B = gmsh.model.occ.addPoint(0.12,  Ly,   0.0, lc_outer)
vD1_line = gmsh.model.occ.addLine(vD1A, vD1B)
vD2_line = gmsh.model.occ.addLine(vD2A, vD2B)



# 11) Interior feature-lines e3…e6 at y = y_D1
interior = [
    ((0.02, y_D1), (0.03, y_D1)),  # e3
    ((0.06, y_D1), (0.07, y_D1)),  # e4
    ((0.00, y_D1), (0.01, y_D1)),  # e5
    ((0.04, y_D1), (0.05, y_D1)),  # e6
]
e_lines = []
for (x1, y1), (x2, y2) in interior:
    pa = gmsh.model.occ.addPoint(x1, y1, 0.0, lc_outer)
    pb = gmsh.model.occ.addPoint(x2, y2, 0.0, lc_outer)
    e_lines.append(gmsh.model.occ.addLine(pa, pb))

# # 10-2) Extra vertical lines for D1 
# vD1x1 = gmsh.model.occ.addPoint(0.065, y_D1*2, 0.0, lc_outer)
# vD1x2 = gmsh.model.occ.addPoint(0.085, y_D1*2, 0.0, lc_outer)
# vD1_line2 = gmsh.model.occ.addLine(vD1B, vD1x2)
# vD1_line3 = gmsh.model.occ.addLine(vD1t, vD1x1)

# 12) Fragment the surface with all split lines
all_lines = (
    bottom_lines + top_lines +
    [right_line, left_line,
     d1_line, d2_line,
     vD1_line, vD2_line] +
    e_lines
)
gmsh.model.occ.fragment([(2, surf)], [(1, l) for l in all_lines])
gmsh.model.occ.synchronize()

# pg1 = gmsh.model.addPhysicalGroup(0, [p1])          # 0 = point‐dimensional group
# gmsh.model.setPhysicalName(0, pg1, "Gamma_p1")      # give it a name

# 13) Physical curve groups
groups = {
    "Gamma_d1":  [], "Gamma_d2":  [], "Gamma_e1":  [],
    "Gamma_e2":  [], "Gamma_f1":  [], "Gamma_e3":  [],
    "Gamma_e4":  [], "Gamma_e5":  [], "Gamma_e6":  [],
    "Gamma_D1s": [], "Gamma_D2s": [],
}
tol = 1e-6

# ———   ⇩ Fixed here ⇩   ———
line_entities = gmsh.model.occ.getEntities(1)
for dim, tag in line_entities:
    x, y, _ = gmsh.model.occ.getCenterOfMass(dim, tag)
    # bottom edge
    if abs(y) < tol:
        if   x <  0.01 + tol:            groups["Gamma_d1"].append(tag)
        elif 0.035 < x < 0.045 + tol:    groups["Gamma_d2"].append(tag)
        elif 0.075 < x < 0.085 + tol:    groups["Gamma_e1"].append(tag)
        if x <= 0.085 + tol:             groups["Gamma_D1s"].append(tag)
    # top edge
    elif abs(y - Ly) < tol:
        if   0.12 < x < 0.125 + tol:     groups["Gamma_e2"].append(tag)
        elif x >  0.14 - tol:           groups["Gamma_f1"].append(tag)
        if x >= 0.12 - tol:              groups["Gamma_D2s"].append(tag)
    # interior e-lines
    elif abs(y - y_D1) < tol:
        if   0.02 <= x <= 0.025:         groups["Gamma_e3"].append(tag)
        elif 0.06 <= x <= 0.065:         groups["Gamma_e4"].append(tag)
        elif 0.00 <= x <= 0.005:         groups["Gamma_e5"].append(tag)
        elif 0.04 <= x <= 0.045:         groups["Gamma_e6"].append(tag)

for name, tags in groups.items():
    if tags:
        pg = gmsh.model.addPhysicalGroup(1, tags)
        gmsh.model.setPhysicalName(1, pg, name)

# 14) Tag the top-right corner point
point_entities = gmsh.model.occ.getEntities(0)
for dim, tag in point_entities:
    x, y, _ = gmsh.model.occ.getCenterOfMass(dim, tag)
    if abs(x - Lx) < tol and abs(y - Ly) < tol:
        pg = gmsh.model.addPhysicalGroup(0, [tag])
        gmsh.model.setPhysicalName(0, pg, "Gamma_topright")
        break

# 15) Physical surface groups
surf_groups = {"Gamma_D1": [], "Gamma_D2": []}
surf_entities = gmsh.model.occ.getEntities(2)
for dim, tag in surf_entities:
    x, y, _ = gmsh.model.occ.getCenterOfMass(dim, tag)
    if 0.0  <= x <= 0.085 + tol and 0.0   <= y <= y_D1 + tol:
        surf_groups["Gamma_D1"].append(tag)
    if 0.12 <= x <= Lx   + tol and 0.065 <= y <= Ly    + tol:
        surf_groups["Gamma_D2"].append(tag)

for name, tags in surf_groups.items():
    if tags:
        pg = gmsh.model.addPhysicalGroup(2, tags)
        gmsh.model.setPhysicalName(2, pg, name)

# Omega = all surfaces
all_surfs = [tag for (_, tag) in surf_entities]
pg = gmsh.model.addPhysicalGroup(2, all_surfs)
gmsh.model.setPhysicalName(2, pg, "Omega")
# ——————————————————————————————————————————————————————————————————
# Build only the D1- and D2-gap segments as physical curve groups
tol = 1e-6

# 1) Cache all 1D entities’ centers
lines = { tag: gmsh.model.occ.getCenterOfMass(1, tag)[:2]
          for (_, tag) in gmsh.model.occ.getEntities(1) }

# 2) D1 gaps along y = y_D1
yD1    = y_D1        # your existing D1 level
e5_rng = (0.00, 0.005)
e3_rng = (0.02, 0.025)
e6_rng = (0.04, 0.045)
e4_rng = (0.06, 0.065)
xD1max = 0.085

D1_gaps = []
for tag, (x,y) in lines.items():
    if abs(y - yD1) < tol:
        # gap between e5 and e3
        if   e5_rng[1] + tol < x < e3_rng[0] - tol:   D1_gaps.append(tag)
        # gap between e3 and e6
        elif e3_rng[1] + tol < x < e6_rng[0] - tol:   D1_gaps.append(tag)
        # gap between e6 and e4
        elif e6_rng[1] + tol < x < e4_rng[0] - tol:   D1_gaps.append(tag)
        # gap between e4 and right border
        elif e4_rng[1] + tol < x < xD1max - tol:      D1_gaps.append(tag)

# 3) D2 gaps along y = 0.065
yD2    = 0.065
e2_rng = (0.12,  0.125)
f1_rng = (0.14, Lx)

D2_gaps = []
for tag, (x,y) in lines.items():
    if abs(y - yD2) < tol:
        # gap between e2 and f1
        if e2_rng[1] + tol < x < f1_rng[0] - tol:
            D2_gaps.append(tag)

# 4) Deduplicate & create the Physical Groups
D1_gaps = sorted(set(D1_gaps))
D2_gaps = sorted(set(D2_gaps))

if D1_gaps:
    pg = gmsh.model.addPhysicalGroup(1, D1_gaps)
    gmsh.model.setPhysicalName(1, pg, "Gamma_D1_gaps")

if D2_gaps:
    pg = gmsh.model.addPhysicalGroup(1, D2_gaps)
    gmsh.model.setPhysicalName(1, pg, "Gamma_D2_gaps")
# ——————————————————————————————————————————————————————————————————
# ——————————————————————————————————————————————————————————————————
# 5) Now pick out the two vertical D1/D2 border lines

# we already have:
#   yD1 = y_D1
#   lines = { tag: (x,y) for (_, tag), (x,y) in lines.items() }

D1_vert, D2_vert = [], []
xD1b = 0.085
xD2b = 0.12

for tag, (x, y) in lines.items():
    # D1 vertical from y=0 up to yD1
    if abs(x - xD1b) < tol and 0.0 - tol <= y <= yD1 + tol:
        D1_vert.append(tag)
    # D2 vertical from y=yD2 up to Ly
    if abs(x - xD2b) < tol and yD2 - tol <= y <= Ly + tol:
        D2_vert.append(tag)

# make them unique
D1_vert = sorted(set(D1_vert))
D2_vert = sorted(set(D2_vert))

# add the Physical Groups
if D1_vert:
    pg = gmsh.model.addPhysicalGroup(1, D1_vert)
    gmsh.model.setPhysicalName(1, pg, "Gamma_D1_vertical")

if D2_vert:
    pg = gmsh.model.addPhysicalGroup(1, D2_vert)
    gmsh.model.setPhysicalName(1, pg, "Gamma_D2_vertical")
# ——————————————————————————————————————————————————————————————————


# ————————————————————————————————————————————————————————————————
# Physical groups for the vertical sides of D1 (bottom band) and D2 (top band)

tol = 1e-6
Lx, Ly = 0.15, 0.07
y_D1 = 0.0046666667
y_D2 = 0.065

D1_vert_boundary = []
D2_vert_boundary = []

for dim, tag in gmsh.model.occ.getEntities(1):
    x, y, _ = gmsh.model.occ.getCenterOfMass(dim, tag)
    # left or right boundary at x=0 or x=Lx
    on_left_D1  = abs(x - 0.0) < tol
    on_right_D1 = abs(x - 0.085   ) < tol

    on_left_D2  = abs(x - 0.12) < tol
    on_right_D2 = abs(x - 0.15   ) < tol

    # segment of that boundary belonging to D1 (0 ≤ y ≤ y_D1)
    if (on_left_D1 or on_right_D1) and 0.0 - tol <= y <= y_D1 + tol:
        D1_vert_boundary.append(tag)

    # segment of that boundary belonging to D2 (y_D2 ≤ y ≤ Ly)
    if (on_left_D2 or on_right_D2) and y_D2 - tol <= y <= Ly + tol:
        D2_vert_boundary.append(tag)

# Make sure each group is unique
D1_vert_boundary = sorted(set(D1_vert_boundary))
D2_vert_boundary = sorted(set(D2_vert_boundary))

# Now create the physical groups
if D1_vert_boundary:
    pg1 = gmsh.model.addPhysicalGroup(1, D1_vert_boundary)
    gmsh.model.setPhysicalName(1, pg1, "Gamma_D1_boundary_vertical")

if D2_vert_boundary:
    pg2 = gmsh.model.addPhysicalGroup(1, D2_vert_boundary)
    gmsh.model.setPhysicalName(1, pg2, "Gamma_D2_boundary_vertical")
# ————————————————————————————————————————————————————————————————

# Add “Gamma_remainder” on the *outer* boundary only

tol = 1e-6

# 1) Collect all 1D curves whose centers lie on one of the 4 domain edges
boundary_lines = []
for dim, tag in gmsh.model.occ.getEntities(1):
    x, y, _ = gmsh.model.occ.getCenterOfMass(dim, tag)
    if abs(y) < tol        or \
       abs(y - Ly) < tol   or \
       abs(x) < tol        or \
       abs(x - Lx) < tol:
        boundary_lines.append(tag)

# 2) Subtract out all the ones you've already put in your named Gamma_* groups
used = set(sum(groups.values(), []))
remainder = sorted(set(boundary_lines) - used)

# 3) Make the Physical Group for those “leftover” outer‐boundary edges
if remainder:
    pg = gmsh.model.addPhysicalGroup(1, remainder)
    gmsh.model.setPhysicalName(1, pg, "Gamma_remainder")



# 16) Mesh-size fields for smooth refinement
f_dist = gmsh.model.mesh.field.add("Distance")
gmsh.model.mesh.field.setNumbers(f_dist, "EdgesList", sum(groups.values(), []))

f_thresh = gmsh.model.mesh.field.add("Threshold")
gmsh.model.mesh.field.setNumber(f_thresh, "IField",  f_dist)
gmsh.model.mesh.field.setNumber(f_thresh, "LcMin",   lc_ref)
gmsh.model.mesh.field.setNumber(f_thresh, "LcMax",   lc_outer)
gmsh.model.mesh.field.setNumber(f_thresh, "DistMin", 0.0)
gmsh.model.mesh.field.setNumber(f_thresh, "DistMax", 0.02)

f_box = gmsh.model.mesh.field.add("Box")
gmsh.model.mesh.field.setNumber(f_box, "XMin", 0.12)
gmsh.model.mesh.field.setNumber(f_box, "XMax", Lx)
gmsh.model.mesh.field.setNumber(f_box, "YMin", 0.065)
gmsh.model.mesh.field.setNumber(f_box, "YMax", Ly)
gmsh.model.mesh.field.setNumber(f_box, "VIn",  lc_outer)
gmsh.model.mesh.field.setNumber(f_box, "VOut", lc_outer)

# --- point‐based refinements around (0.14, 0.07) and (0.075, 0.0) ---

# Ball around top split at (0.14, 0.07)
f_ball1 = gmsh.model.mesh.field.add("Ball")
gmsh.model.mesh.field.setNumber(f_ball1, "XCenter", 0.14)
gmsh.model.mesh.field.setNumber(f_ball1, "YCenter", 0.07)
gmsh.model.mesh.field.setNumber(f_ball1, "ZCenter", 0.0)
gmsh.model.mesh.field.setNumber(f_ball1, "Radius",  0.005)   # adjust radius as needed
gmsh.model.mesh.field.setNumber(f_ball1, "VIn",     lc_extra_ref)
gmsh.model.mesh.field.setNumber(f_ball1, "VOut",    lc_outer)

# Ball around bottom split at (0.075, 0.0)
f_ball2 = gmsh.model.mesh.field.add("Ball")
gmsh.model.mesh.field.setNumber(f_ball2, "XCenter", 0.075)
gmsh.model.mesh.field.setNumber(f_ball2, "YCenter", 0.0)
gmsh.model.mesh.field.setNumber(f_ball2, "ZCenter", 0.0)
gmsh.model.mesh.field.setNumber(f_ball2, "Radius",  0.005*(3/4))
gmsh.model.mesh.field.setNumber(f_ball2, "VIn",     lc_extra_ref)
gmsh.model.mesh.field.setNumber(f_ball2, "VOut",    lc_outer)

# f_ball3 = gmsh.model.mesh.field.add("Ball")
# gmsh.model.mesh.field.setNumber(f_ball3, "XCenter", 0.12)
# gmsh.model.mesh.field.setNumber(f_ball3, "YCenter", 0.065)
# gmsh.model.mesh.field.setNumber(f_ball3, "ZCenter", 0.0)
# gmsh.model.mesh.field.setNumber(f_ball3, "Radius",  0.005)
# gmsh.model.mesh.field.setNumber(f_ball3, "VIn",     lc_extra_ref)
# gmsh.model.mesh.field.setNumber(f_ball3, "VOut",    lc_outer)

# now combine: distance‐threshold + box + both balls
f_min = gmsh.model.mesh.field.add("Min")
gmsh.model.mesh.field.setNumbers(f_min, "FieldsList",
    [ f_thresh,    # your distance‐based refinement
      f_box,       # your D2 box refinement
      f_ball1,     # top point refinement
      f_ball2      # bottom point refinement
    ]
)
gmsh.model.mesh.field.setAsBackgroundMesh(f_min)


# f_min = gmsh.model.mesh.field.add("Min")
# gmsh.model.mesh.field.setNumbers(f_min, "FieldsList", [f_thresh, f_box])
# gmsh.model.mesh.field.setAsBackgroundMesh(f_min)

# 17) Generate and write
gmsh.model.mesh.generate(2)
gmsh.write("domain_with_e_coarse.msh")
print("Wrote domain_with_e_coarse.msh")

# 18) Finalize
gmsh.finalize()

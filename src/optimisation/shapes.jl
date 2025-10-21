# 1) Circle symmetric about y = 0.035
function circle(R, x, δ)
    # base center on y = 0.035
    c0x, c0y = 0.025, 0.025#0.035
    cx, cy   = c0x + δ, c0y

    dx, dy = x[1] - cx, x[2] - cy
    # positive inside, zero on boundary, negative outside
    return -(sqrt(dx^2 + dy^2) - R)
end

# 1) Axis-aligned square, symmetric about y = 0.035
function square(R, x, δ)

    # compute center
    cx = 0.025 + δ
    cy = 0.025#0.035

    # half-width = R
    dx = abs(x[1] - cx) - R
    dy = abs(x[2] - cy) - R

    # level-set: negative outside, zero on edge, positive inside
    val = -max(dx, dy)
    return val
end

# 2) Diamond (45°-rotated square), symmetric about y = 0.035
function diamond(R, x, δ)
    # compute center
    cx = 0.025 + δ
    cy = 0.025#0.035

    # diamond uses L1 norm
    dx = abs(x[1] - cx)
    dy = abs(x[2] - cy)

    # level-set: negative outside, zero on edge, positive inside
    val = -(dx + dy - R)
    return val
end



# 4) Equilateral Triangle symmetric about y = 0.035
#    (apex to the right, base vertical)
function triangle(R, x, δ)
    # 1) Base center on y = 0.035, shifted by δ
    c0x, c0y = 0.025, 0.025
    cx, cy   = c0x + δ, c0y

    # 2) Local (u,v) coords = (x-cx, y-cy)
    u = x[1] - cx
    v = x[2] - cy

    # 3) Flipped triangle vertices in (u,v):
    #    Apex at (-R, 0), and the other two points mirrored accordingly
    h  = sqrt(3)/2 * R
    d  = R/2
    V1 = (-R,   0.0)   # apex now pointing left
    V2 = ( d,   h )
    V3 = ( d,  -h )

    # 4) Signed‐distance to line from a→b
    sd(a, b, p) = ((b[2] - a[2])*(p[1] - a[1]) - (b[1] - a[1])*(p[2] - a[2])) /
                  hypot(b[1] - a[1], b[2] - a[2])

    # 5) Determine “inside” side at (0,0)
    s1 = sign(sd(V1, V2, (0.0, 0.0)))
    s2 = sign(sd(V2, V3, (0.0, 0.0)))
    s3 = sign(sd(V3, V1, (0.0, 0.0)))

    # 6) Oriented half‐space distances at (u,v)
    d1 = -s1 * sd(V1, V2, (u, v))
    d2 = -s2 * sd(V2, V3, (u, v))
    d3 = -s3 * sd(V3, V1, (u, v))

    # 7) Combined level‐set: positive inside, zero on edges
    return -max(d1, d2, d3)
end


# 5-pointed star
function star5(R, x, δ)
    # center on y = 0.035, shifted horizontally by δ
    c0x, c0y = 0.025, 0.025#0.035
    cx, cy   = c0x + δ, c0y

    # inner radius ratio
    r_inner = 0.4 * R
    rot     = 0.0      # first outer tip to the right

    # build 10 vertices (5 outer, 5 inner)
    verts = Tuple{Float64,Float64}[]
    for i in 0:4
        θo = 2π * i / 5 + rot
        θi = θo + π/5
        push!(verts, (cx + R      * cos(θo), cy + R      * sin(θo)))
        push!(verts, (cx + r_inner * cos(θi), cy + r_inner * sin(θi)))
    end

    # point→segment distance
    function pt_seg_dist(p, a, b)
        vx, vy = b[1]-a[1], b[2]-a[2]
        wx, wy = p[1]-a[1], p[2]-a[2]
        c1     = vx*wx + vy*wy
        if c1 ≤ 0 return hypot(wx, wy) end
        c2 = vx^2 + vy^2
        if c2 ≤ c1 return hypot(p[1]-b[1], p[2]-b[2]) end
        t    = c1/c2
        proj = (a[1] + t*vx, a[2] + t*vy)
        return hypot(p[1]-proj[1], p[2]-proj[2])
    end

    # inside-test via ray-crossing
    function point_in_poly(p, vs)
        x0, y0 = p
        inside = false
        j      = length(vs)
        for i in 1:length(vs)
            xi, yi = vs[i]
            xj, yj = vs[j]
            if (yi>y0) != (yj>y0) && x0 < (xj-xi)*(y0-yi)/(yj-yi) + xi
                inside = !inside
            end
            j = i
        end
        return inside
    end

    # compute min distance to any edge
    dmin = Inf
    for i in 1:length(verts)
        j    = i == length(verts) ? 1 : i+1
        dmin = min(dmin, pt_seg_dist(x, verts[i], verts[j]))
    end

    # signed level-set: positive inside, negative outside
    return point_in_poly(x, verts) ? dmin : -dmin
end

# 6-pointed star
function star6(R, x, δ)
    # center on y = 0.035, shifted horizontally by δ
    c0x, c0y = 0.025, 0.025#0.035
    cx, cy   = c0x + δ, c0y

    # inner radius ratio
    r_inner = 0.5 * R
    rot     = 0.0      # first outer tip to the right

    # build 12 vertices (6 outer, 6 inner)
    verts = Tuple{Float64,Float64}[]
    for i in 0:5
        θo = 2π * i / 6 + rot
        θi = θo + π/6
        push!(verts, (cx + R      * cos(θo), cy + R      * sin(θo)))
        push!(verts, (cx + r_inner * cos(θi), cy + r_inner * sin(θi)))
    end

    # point→segment distance
    function pt_seg_dist(p, a, b)
        vx, vy = b[1]-a[1], b[2]-a[2]
        wx, wy = p[1]-a[1], p[2]-a[2]
        c1     = vx*wx + vy*wy
        if c1 ≤ 0 return hypot(wx, wy) end
        c2 = vx^2 + vy^2
        if c2 ≤ c1 return hypot(p[1]-b[1], p[2]-b[2]) end
        t    = c1/c2
        proj = (a[1] + t*vx, a[2] + t*vy)
        return hypot(p[1]-proj[1], p[2]-proj[2])
    end

    # inside-test via ray-crossing
    function point_in_poly(p, vs)
        x0, y0 = p
        inside = false
        j      = length(vs)
        for i in 1:length(vs)
            xi, yi = vs[i]
            xj, yj = vs[j]
            if (yi>y0) != (yj>y0) && x0 < (xj-xi)*(y0-yi)/(yj-yi) + xi
                inside = !inside
            end
            j = i
        end
        return inside
    end

    # compute min distance to any edge
    dmin = Inf
    for i in 1:length(verts)
        j    = i == length(verts) ? 1 : i+1
        dmin = min(dmin, pt_seg_dist(x, verts[i], verts[j]))
    end

    # signed level-set: positive inside, negative outside
    return point_in_poly(x, verts) ? dmin : -dmin
end

# 8-pointed star (octagram), centered at y = 0.025 (shifted horizontally by δ)
function star8(R, x, δ)
    # 1) Base center on y = 0.025
    c0x, c0y = 0.025, 0.025
    cx, cy   = c0x + δ, c0y

    # 2) Choose inner‐radius ratio (e.g. 0.4 of outer)
    r_inner = 0.7 * R
    rot     = 0.0  # first outer tip to the right

    # 3) Build 16 alternating vertices (8 outer, 8 inner)
    verts = Tuple{Float64,Float64}[]
    for i in 0:7
        θo = 2π * i / 8 + rot
        θi = θo + π/8
        push!(verts, (cx + R      * cos(θo), cy + R      * sin(θo)))
        push!(verts, (cx + r_inner * cos(θi), cy + r_inner * sin(θi)))
    end

    # 4) Distance from point p to segment a→b
    function pt_seg_dist(p, a, b)
        vx, vy = b[1]-a[1], b[2]-a[2]
        wx, wy = p[1]-a[1], p[2]-a[2]
        c1     = vx*wx + vy*wy
        if c1 ≤ 0
            return hypot(wx, wy)
        end
        c2 = vx^2 + vy^2
        if c2 ≤ c1
            return hypot(p[1]-b[1], p[2]-b[2])
        end
        t    = c1/c2
        proj = (a[1] + t*vx, a[2] + t*vy)
        return hypot(p[1]-proj[1], p[2]-proj[2])
    end

    # 5) Ray‐crossing point‐in‐polygon test
    function point_in_poly(p, vs)
        x0, y0 = p
        inside = false
        j      = length(vs)
        for i in 1:length(vs)
            xi, yi = vs[i]
            xj, yj = vs[j]
            if (yi>y0) != (yj>y0) && x0 < (xj-xi)*(y0-yi)/(yj-yi) + xi
                inside = !inside
            end
            j = i
        end
        return inside
    end

    # 6) Compute minimal distance to any edge
    dmin = Inf
    n    = length(verts)
    for i in 1:n
        j    = (i == n) ? 1 : i+1
        dmin = min(dmin, pt_seg_dist(x, verts[i], verts[j]))
    end

    # 7) Signed level‐set: positive inside, negative outside
    return point_in_poly(x, verts) ? dmin : -dmin
end

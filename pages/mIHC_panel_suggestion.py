import itertools
import pandas as pd
import streamlit as st

# =========================
# CHANNEL CONFIG
# =========================
CHANNELS = [480, 520, 540, 570, 620, 650, 690, 780]

PREFERRED = [520, 570, 620, 690]
INTERMEDIATE = [540, 650]
LATE = [780]

LOCATION_OPTIONS = ["nucleus", "cytoplasm", "membrane"]
STRENGTH_OPTIONS = ["weak", "medium", "strong"]


# =========================
# RISK SYSTEM
# =========================
def pair_risk(c1, c2, same_location=False, weak=False, strong=False, morph_diff=False):
    if c1 == c2:
        return 999

    if 780 in [c1, c2]:
        base = 0
    elif 480 in [c1, c2]:
        base = 1
    else:
        dist = abs(CHANNELS.index(c1) - CHANNELS.index(c2))
        mapping = {1: 5, 2: 4, 3: 3, 4: 2, 5: 1, 6: 0}
        base = mapping.get(dist, 0)

    if same_location and not morph_diff:
        base += 1
    if weak:
        base += 1
    if strong:
        base -= 1

    return max(0, base)


def spacing_rule_ok(channels):
    ordered = sorted(channels, key=lambda x: CHANNELS.index(x))
    for i in range(len(ordered) - 1):
        if abs(CHANNELS.index(ordered[i + 1]) - CHANNELS.index(ordered[i])) < 2:
            return False
    return True


def spread_penalty(channels):
    ordered = sorted(channels, key=lambda x: CHANNELS.index(x))
    penalty = 0
    for i in range(len(ordered) - 1):
        dist = abs(CHANNELS.index(ordered[i + 1]) - CHANNELS.index(ordered[i]))
        if dist == 2:
            penalty += 4
        elif dist == 3:
            penalty += 1
    return penalty


def late_channel_penalty(channels, n_markers):
    penalty = 0
    if 780 in channels:
        if n_markers <= 5:
            penalty += 12
        elif n_markers == 6:
            penalty += 6
    return penalty


def norm_pair(a, b):
    return tuple(sorted([a, b]))


def checkpoint_rule_ok(assign, df, checkpoint_pairs):
    check = set(norm_pair(a, b) for a, b in checkpoint_pairs if a and b and a != b)

    for i in range(len(df)):
        for j in range(i + 1, len(df)):
            m1 = df.iloc[i]["marker"]
            m2 = df.iloc[j]["marker"]

            if norm_pair(m1, m2) in check:
                idx1 = CHANNELS.index(assign[i])
                idx2 = CHANNELS.index(assign[j])

                if abs(idx1 - idx2) < 2:
                    return False
    return True


def total_risk(df, assign, checkpoint_pairs, morph_pairs):
    score = 0
    morph = set(norm_pair(a, b) for a, b in morph_pairs if a and b and a != b)
    check = set(norm_pair(a, b) for a, b in checkpoint_pairs if a and b and a != b)

    for i in range(len(df)):
        for j in range(i + 1, len(df)):
            a = df.iloc[i]
            b = df.iloc[j]

            same = a["location"] == b["location"]
            weak = "weak" in [a["strength"], b["strength"]]
            strong = "strong" in [a["strength"], b["strength"]]

            morph_diff = norm_pair(a["marker"], b["marker"]) in morph

            score += pair_risk(assign[i], assign[j], same, weak, strong, morph_diff)

    for i in range(len(df)):
        for j in range(i + 1, len(df)):
            a = df.iloc[i]["marker"]
            b = df.iloc[j]["marker"]

            if norm_pair(a, b) in check:
                dist = abs(CHANNELS.index(assign[i]) - CHANNELS.index(assign[j]))
                score += max(0, 4 - dist) * 2

    return score


def suggest_panel(df, seg_marker=None, fixed=None, checkpoint_pairs=None, morph_pairs=None):
    fixed = fixed or {}
    checkpoint_pairs = checkpoint_pairs or []
    morph_pairs = morph_pairs or []

    reserved = []
    reserved_channels = []

    if seg_marker and seg_marker in df["marker"].values:
        row = df[df["marker"] == seg_marker].copy()
        row["channel"] = 480
        reserved.append(row)
        reserved_channels.append(480)
        df = df[df["marker"] != seg_marker]

    for m, c in fixed.items():
        if m in df["marker"].values:
            row = df[df["marker"] == m].copy()
            row["channel"] = c
            reserved.append(row)
            reserved_channels.append(c)
            df = df[df["marker"] != m]

    n = len(df)

    base = [c for c in PREFERRED if c not in reserved_channels]
    mid = [c for c in INTERMEDIATE if c not in reserved_channels]
    late = [c for c in LATE if c not in reserved_channels]
    pool = base + mid + late

    if n == 0:
        return pd.concat(reserved, ignore_index=True), (0, 0, 0)

    if len(pool) < n:
        return None, None

    best = None
    best_tuple = None

    for combo in itertools.combinations(pool, n):
        if not spacing_rule_ok(combo):
            continue

        for perm in itertools.permutations(combo, n):
            if not checkpoint_rule_ok(perm, df, checkpoint_pairs):
                continue

            risk = total_risk(df, perm, checkpoint_pairs, morph_pairs)
            spread = spread_penalty(perm)
            late_pen = late_channel_penalty(perm, n)

            score_tuple = (late_pen, risk, spread)

            if best_tuple is None or score_tuple < best_tuple:
                best_tuple = score_tuple
                best = perm

    if best is None:
        for combo in itertools.combinations(pool, n):
            for perm in itertools.permutations(combo, n):
                if not checkpoint_rule_ok(perm, df, checkpoint_pairs):
                    continue

                risk = total_risk(df, perm, checkpoint_pairs, morph_pairs)
                spread = spread_penalty(perm)
                late_pen = late_channel_penalty(perm, n)

                score_tuple = (late_pen, risk, spread)

                if best_tuple is None or score_tuple < best_tuple:
                    best_tuple = score_tuple
                    best = perm

    if best is None:
        return None, None

    df = df.copy()
    df["channel"] = best
    result = pd.concat(reserved + [df], ignore_index=True) if reserved else df

    return result, best_tuple


# =========================
# STREAMLIT UI
# =========================
st.title("mIHC Panel Planner")

marker_count = st.number_input("Number of markers", min_value=1, max_value=8, value=4, step=1)

markers = []
st.subheader("Marker Information")
for i in range(marker_count):
    c1, c2, c3 = st.columns(3)
    with c1:
        name = st.text_input(f"Marker {i+1} name", key=f"name_{i}")
    with c2:
        location = st.selectbox(f"Location {i+1}", LOCATION_OPTIONS, key=f"loc_{i}")
    with c3:
        strength = st.selectbox(f"Strength {i+1}", STRENGTH_OPTIONS, key=f"str_{i}")

    markers.append({
        "marker": name.strip() if name.strip() else f"Marker{i+1}",
        "location": location,
        "strength": strength
    })

df = pd.DataFrame(markers)
names = df["marker"].tolist()

st.subheader("Segmentation Marker")
use_seg = st.radio("Use segmentation marker at 480?", ["No", "Yes"], horizontal=True)
seg_marker = None
if use_seg == "Yes":
    seg_marker = st.selectbox("Segmentation marker", names)

st.subheader("Fixed Channels")
fixed = {}
fixed_n = st.number_input("Number of fixed markers", min_value=0, max_value=8, value=0, step=1)
for i in range(fixed_n):
    c1, c2 = st.columns(2)
    with c1:
        m = st.selectbox(f"Fixed marker {i+1}", names, key=f"fixed_m_{i}")
    with c2:
        ch = st.selectbox(f"Fixed channel {i+1}", CHANNELS, key=f"fixed_c_{i}")
    fixed[m] = ch

st.subheader("Checkpoint Pairs")
checkpoint_pairs = []
cp_n = st.number_input("Number of checkpoint pairs", min_value=0, max_value=10, value=0, step=1)
for i in range(cp_n):
    c1, c2 = st.columns(2)
    with c1:
        a = st.selectbox(f"Checkpoint A{i+1}", names, key=f"cp_a_{i}")
    with c2:
        b = st.selectbox(f"Checkpoint B{i+1}", names, key=f"cp_b_{i}")
    checkpoint_pairs.append((a, b))

st.subheader("Morphology-Different Pairs")
morph_pairs = []
mp_n = st.number_input("Number of morphology-different pairs", min_value=0, max_value=10, value=0, step=1)
for i in range(mp_n):
    c1, c2 = st.columns(2)
    with c1:
        a = st.selectbox(f"Morph A{i+1}", names, key=f"mp_a_{i}")
    with c2:
        b = st.selectbox(f"Morph B{i+1}", names, key=f"mp_b_{i}")
    morph_pairs.append((a, b))

if st.button("Suggest Panel"):
    result, score = suggest_panel(
        df=df,
        seg_marker=seg_marker,
        fixed=fixed,
        checkpoint_pairs=checkpoint_pairs,
        morph_pairs=morph_pairs
    )

    if result is None:
        st.error("No valid panel found under current constraints.")
    else:
        st.success("Best panel found.")
        st.dataframe(result, use_container_width=True)
        st.write("Score tuple (late penalty, total risk, spread penalty):", score)

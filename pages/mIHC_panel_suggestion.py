import itertools
import pandas as pd
import streamlit as st

CHANNELS = [480, 520, 540, 570, 620, 650, 690, 780]

# Brightness groups
BRIGHT = [520, 650]
MEDIUM = [540, 570, 620]
DIM = [690]
LATE = [780]

LOCATION_OPTIONS = ["nucleus", "cytoplasm", "membrane"]
STRENGTH_OPTIONS = ["weak", "medium", "strong"]


def pair_risk(c1, c2, same_location=False, weak=False, strong=False, morph_diff=False):
    if c1 == c2:
        return 999

    if 780 in [c1, c2]:
        base = 0
    elif 480 in [c1, c2]:
        base = 1
    else:
        i1 = CHANNELS.index(c1)
        i2 = CHANNELS.index(c2)
        dist = abs(i1 - i2)

        mapping = {
            1: 5,
            2: 4,
            3: 3,
            4: 2,
            5: 1,
            6: 0,
        }
        base = mapping.get(dist, 0)

    if same_location and not morph_diff:
        base += 1
    if weak:
        base += 1
    if strong:
        base -= 1

    return max(0, base)


def brightness_penalty(marker_strength, channel):
    """
    Lower score = better match between marker intensity and Opal brightness.
    Main idea:
    - weak markers should prefer bright channels
    - medium markers should prefer medium channels
    - strong markers can tolerate dimmer channels better
    - 780 is discouraged unless needed
    - 480 is treated as a special segmentation slot, not a general preferred channel
    """
    if marker_strength == "weak":
        if channel in BRIGHT:
            return 0
        if channel in MEDIUM:
            return 3
        if channel in DIM:
            return 8
        if channel == 780:
            return 12
        if channel == 480:
            return 4

    elif marker_strength == "medium":
        if channel in MEDIUM:
            return 0
        if channel in BRIGHT:
            return 1
        if channel in DIM:
            return 3
        if channel == 780:
            return 7
        if channel == 480:
            return 3

    elif marker_strength == "strong":
        if channel in DIM:
            return 0
        if channel in MEDIUM:
            return 1
        if channel in BRIGHT:
            return 4
        if channel == 780:
            return 2
        if channel == 480:
            return 2

    return 0


def spacing_rule_ok(channels):
    ordered = sorted(channels, key=lambda x: CHANNELS.index(x))

    for i in range(len(ordered) - 1):
        c1 = ordered[i]
        c2 = ordered[i + 1]
        dist = abs(CHANNELS.index(c1) - CHANNELS.index(c2))

        if dist < 2:
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


def total_risk(df, assign, morph_pairs_set, checkpoint_pairs_set):
    score = 0

    # 1) Brightness matching for each marker-channel assignment
    for i in range(len(df)):
        score += brightness_penalty(df.iloc[i]["strength"], assign[i])

    # 2) Pairwise interaction / spectral-spacing logic
    for i in range(len(df)):
        for j in range(i + 1, len(df)):
            a = df.iloc[i]
            b = df.iloc[j]

            same = a["location"] == b["location"]
            weak = "weak" in [a["strength"], b["strength"]]
            strong = "strong" in [a["strength"], b["strength"]]

            morph_diff = norm_pair(a["marker"], b["marker"]) in morph_pairs_set

            score += pair_risk(assign[i], assign[j], same, weak, strong, morph_diff)

    # 3) Checkpoint pairs should be further apart
    idx = {df.iloc[i]["marker"]: i for i in range(len(df))}
    for m1, m2 in checkpoint_pairs_set:
        if m1 in idx and m2 in idx:
            i = idx[m1]
            j = idx[m2]
            dist = abs(CHANNELS.index(assign[i]) - CHANNELS.index(assign[j]))
            score += max(0, 4 - dist) * 2

    return score


st.title("mIHC panel suggestion")
st.write("Suggest Opal channel assignments for multiplex IHC panels.")

marker_count = st.number_input("Number of markers", min_value=1, max_value=8, value=4, step=1)

st.subheader("Marker information")

markers = []
for i in range(marker_count):
    col1, col2, col3 = st.columns(3)

    with col1:
        marker_name = st.text_input(f"Marker {i+1} name", key=f"marker_name_{i}")

    with col2:
        location = st.selectbox(
            f"Location {i+1}",
            LOCATION_OPTIONS,
            key=f"location_{i}",
        )

    with col3:
        strength = st.selectbox(
            f"Strength {i+1}",
            STRENGTH_OPTIONS,
            key=f"strength_{i}",
        )

    if marker_name.strip() == "":
        marker_name = f"Marker{i+1}"

    markers.append({
        "marker": marker_name,
        "location": location,
        "strength": strength,
    })

df = pd.DataFrame(markers)
names = df["marker"].tolist()

st.subheader("Logic settings")

seg_yes = st.radio("Use a segmentation marker for Opal 480?", ["No", "Yes"], horizontal=True)
seg_marker = None
if seg_yes == "Yes":
    seg_marker = st.selectbox("Segmentation marker", names, key="seg_marker")

fixed_yes = st.radio("Use fixed channels?", ["No", "Yes"], horizontal=True)
fixed_map = {}
if fixed_yes == "Yes":
    fixed_n = st.number_input(
        "Number of fixed markers",
        min_value=0,
        max_value=min(8, len(names)),
        value=0,
        step=1,
    )
    st.write("Define fixed marker-channel pairs")
    for i in range(fixed_n):
        c1, c2 = st.columns(2)
        with c1:
            m = st.selectbox(f"Fixed marker {i+1}", names, key=f"fixed_marker_{i}")
        with c2:
            c = st.selectbox(f"Fixed channel {i+1}", CHANNELS, key=f"fixed_channel_{i}")
        fixed_map[m] = c

checkpoint_pairs_set = set()
checkpoint_n = st.number_input("Number of checkpoint pairs", min_value=0, max_value=10, value=0, step=1)
if checkpoint_n > 0:
    st.write("Checkpoint pairs should be placed further apart")
for i in range(checkpoint_n):
    c1, c2 = st.columns(2)
    with c1:
        a = st.selectbox(f"Checkpoint pair {i+1} - marker A", names, key=f"checkpoint_a_{i}")
    with c2:
        b = st.selectbox(f"Checkpoint pair {i+1} - marker B", names, key=f"checkpoint_b_{i}")
    if a != b:
        checkpoint_pairs_set.add(norm_pair(a, b))

morph_pairs_set = set()
morph_n = st.number_input("Number of morphology-different pairs", min_value=0, max_value=10, value=0, step=1)
if morph_n > 0:
    st.write("These pairs can tolerate being closer because morphology is distinguishable")
for i in range(morph_n):
    c1, c2 = st.columns(2)
    with c1:
        a = st.selectbox(f"Morph pair {i+1} - marker A", names, key=f"morph_a_{i}")
    with c2:
        b = st.selectbox(f"Morph pair {i+1} - marker B", names, key=f"morph_b_{i}")
    if a != b:
        morph_pairs_set.add(norm_pair(a, b))

if st.button("Suggest panel"):
    work_df = df.copy()

    reserved = []
    reserved_channels = []

    if seg_marker is not None and seg_marker in work_df["marker"].values:
        row = work_df[work_df["marker"] == seg_marker].copy()
        row["channel"] = 480
        reserved.append(row)
        reserved_channels.append(480)
        work_df = work_df[work_df["marker"] != seg_marker]

    for m, c in fixed_map.items():
        if m in work_df["marker"].values:
            row = work_df[work_df["marker"] == m].copy()
            row["channel"] = c
            reserved.append(row)
            reserved_channels.append(c)
            work_df = work_df[work_df["marker"] != m]

    if len(set(reserved_channels)) != len(reserved_channels):
        st.error("Two reserved markers are using the same channel. Please revise segmentation/fixed settings.")
    else:
        n = len(work_df)

        base = [c for c in BRIGHT if c not in reserved_channels]
        mid = [c for c in MEDIUM if c not in reserved_channels]
        dim = [c for c in DIM if c not in reserved_channels]
        late = [c for c in LATE if c not in reserved_channels]

        pool = base + mid + dim + late

        if n == 0:
            result = pd.concat(reserved, ignore_index=True) if reserved else pd.DataFrame()
            st.success("All markers are already assigned.")
            st.dataframe(result, use_container_width=True)

        elif len(pool) < n:
            st.error("Not enough available channels for the number of unfixed markers.")

        else:
            best = None
            best_tuple = None

            for combo in itertools.combinations(pool, n):
                if not spacing_rule_ok(combo):
                    continue

                for perm in itertools.permutations(combo, n):
                    risk = total_risk(work_df, perm, morph_pairs_set, checkpoint_pairs_set)
                    spread = spread_penalty(perm)
                    late_pen = late_channel_penalty(perm, n)

                    # Priority order:
                    # 1. Avoid 780 if possible
                    # 2. Minimize overall risk (now includes brightness matching)
                    # 3. Maximize spread
                    score_tuple = (late_pen, risk, spread)

                    if best_tuple is None or score_tuple < best_tuple:
                        best_tuple = score_tuple
                        best = perm

            if best is None:
                for combo in itertools.combinations(pool, n):
                    for perm in itertools.permutations(combo, n):
                        risk = total_risk(work_df, perm, morph_pairs_set, checkpoint_pairs_set)
                        spread = spread_penalty(perm)
                        late_pen = late_channel_penalty(perm, n)

                        score_tuple = (late_pen, risk, spread)

                        if best_tuple is None or score_tuple < best_tuple:
                            best_tuple = score_tuple
                            best = perm

            result_df = work_df.copy()
            result_df["channel"] = best

            result = pd.concat(reserved + [result_df], ignore_index=True) if reserved else result_df
            result = result[["marker", "location", "strength", "channel"]]

            st.success("Suggested panel generated.")
            st.dataframe(result, use_container_width=True)

            if best_tuple is not None:
                st.caption(
                    f"Scoring summary — late channel penalty: {best_tuple[0]}, "
                    f"risk score: {best_tuple[1]}, spread penalty: {best_tuple[2]}"
                )

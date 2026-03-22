import itertools
import pandas as pd
import ipywidgets as widgets
from IPython.display import display, clear_output

CHANNELS = [480, 520, 540, 570, 620, 650, 690, 780]

PREFERRED = [520, 570, 620, 690]
INTERMEDIATE = [540, 650]
LATE = [780]

LOCATION_OPTIONS = ["nucleus", "cytoplasm", "membrane"]
STRENGTH_OPTIONS = ["weak", "medium", "strong"]

marker_widgets = []
checkpoint_pairs = []
morph_pairs = []
fixed_rows = []

logic_state = {}

title = widgets.HTML("<h3>mIHC Panel Planner</h3>")
marker_count = widgets.BoundedIntText(value=4, min=1, max=8, description="Markers")
generate_button = widgets.Button(description="Generate markers", button_style="info")
logic_button = widgets.Button(description="Build logic", button_style="warning")
suggest_button = widgets.Button(description="Suggest panel", button_style="success")

form_out = widgets.Output()
logic_out = widgets.Output()
result_out = widgets.Output()


# =========================
# RISK
# =========================
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
            6: 0
        }
        base = mapping.get(dist, 0)

    if same_location and not morph_diff:
        base += 1
    if weak:
        base += 1
    if strong:
        base -= 1

    return max(0, base)


def spacing_rule_ok(channels):
    """
    Global spacing rule:
    adjacent channels are not allowed anywhere in the panel.
    Example:
    520-540 -> invalid
    540-570 -> invalid
    570-620 -> invalid
    540-620 -> valid
    """
    ordered = sorted(channels, key=lambda x: CHANNELS.index(x))

    for i in range(len(ordered) - 1):
        c1 = ordered[i]
        c2 = ordered[i + 1]
        dist = abs(CHANNELS.index(c1) - CHANNELS.index(c2))

        if dist < 2:
            return False

    return True


def spread_penalty(channels):
    """
    Prefer more spread-out panels.
    """
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
    """
    Discourage using 780 too early.
    """
    penalty = 0
    if 780 in channels:
        if n_markers <= 5:
            penalty += 12
        elif n_markers == 6:
            penalty += 6
    return penalty


# =========================
# DATA
# =========================
def collect_df():
    rows = []
    for i, w in enumerate(marker_widgets):
        name = w["name"].value.strip()
        if name == "":
            name = f"Marker{i+1}"

        rows.append({
            "marker": name,
            "location": w["location"].value,
            "strength": w["strength"].value
        })

    return pd.DataFrame(rows)


def norm_pair(a, b):
    return tuple(sorted([a, b]))


def morph_set():
    s = set()
    for a, b in morph_pairs:
        if a.value and b.value and a.value != b.value:
            s.add(norm_pair(a.value, b.value))
    return s


def checkpoint_set():
    s = set()
    for a, b in checkpoint_pairs:
        if a.value and b.value and a.value != b.value:
            s.add(norm_pair(a.value, b.value))
    return s


# =========================
# NEW HARD RULE FOR CHECKPOINT PAIRS
# =========================
def checkpoint_rule_ok(assign, df):
    """
    HARD CONSTRAINT:
    checkpoint pairs must NOT be adjacent.
    They must be at least distance >= 2 in CHANNELS.

    Example:
    540 with checkpoint partner must be 620+ or 520-
    620 cannot be paired with 540 or 650 if they are checkpoint pairs
    """
    check = checkpoint_set()

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


# =========================
# SEGMENTATION
# =========================
def get_seg_marker():
    if "seg_yes" in logic_state and logic_state["seg_yes"].value == "Yes":
        m = logic_state["seg_text"].value.strip()
        return m if m else None
    return None


# =========================
# FIXED CHANNELS
# =========================
def get_fixed():
    fixed_map = {}
    if "fixed_yes" not in logic_state:
        return fixed_map

    if logic_state["fixed_yes"].value != "Yes":
        return fixed_map

    for m, c in fixed_rows:
        if m.value:
            fixed_map[m.value] = c.value

    return fixed_map


# =========================
# TOTAL RISK
# =========================
def total_risk(df, assign):
    score = 0
    morph = morph_set()
    check = checkpoint_set()

    for i in range(len(df)):
        for j in range(i + 1, len(df)):
            a = df.iloc[i]
            b = df.iloc[j]

            same = a["location"] == b["location"]
            weak = "weak" in [a["strength"], b["strength"]]
            strong = "strong" in [a["strength"], b["strength"]]

            morph_diff = norm_pair(a["marker"], b["marker"]) in morph

            score += pair_risk(assign[i], assign[j], same, weak, strong, morph_diff)

    # soft checkpoint penalty still kept
    for i in range(len(df)):
        for j in range(i + 1, len(df)):
            a = df.iloc[i]["marker"]
            b = df.iloc[j]["marker"]

            if norm_pair(a, b) in check:
                dist = abs(CHANNELS.index(assign[i]) - CHANNELS.index(assign[j]))
                score += max(0, 4 - dist) * 2

    return score


# =========================
# FORM
# =========================
def build_form(_):
    global marker_widgets
    marker_widgets = []

    with form_out:
        clear_output()

        header = widgets.HBox([
            widgets.HTML("<b style='width:160px'>Marker</b>"),
            widgets.HTML("<b style='width:160px'>Location</b>"),
            widgets.HTML("<b style='width:160px'>Strength</b>")
        ])

        rows = [header]

        for i in range(marker_count.value):
            name = widgets.Text(layout=widgets.Layout(width="160px"))
            loc = widgets.Dropdown(options=LOCATION_OPTIONS, layout=widgets.Layout(width="160px"))
            strength = widgets.Dropdown(options=STRENGTH_OPTIONS, layout=widgets.Layout(width="160px"))

            marker_widgets.append({
                "name": name,
                "location": loc,
                "strength": strength
            })

            rows.append(widgets.HBox([name, loc, strength]))

        display(widgets.VBox(rows))


# =========================
# LOGIC BOXES
# =========================
def build_logic(_):
    global checkpoint_pairs, morph_pairs, fixed_rows
    checkpoint_pairs = []
    morph_pairs = []
    fixed_rows = []

    df = collect_df()
    names = df["marker"].tolist()

    with logic_out:
        clear_output()

        display(widgets.HTML("<b>Segmentation marker (480)</b>"))
        seg_yes = widgets.RadioButtons(options=["Yes", "No"], value="No")
        seg_text = widgets.Text(description="Marker")
        seg_box = widgets.Output()

        display(seg_yes)
        display(seg_box)

        def seg_update(*args):
            with seg_box:
                clear_output()
                if seg_yes.value == "Yes":
                    display(seg_text)

        seg_yes.observe(seg_update, names="value")
        seg_update()

        logic_state["seg_yes"] = seg_yes
        logic_state["seg_text"] = seg_text

        display(widgets.HTML("<hr><b>Fixed channels</b>"))
        fixed_yes = widgets.RadioButtons(options=["Yes", "No"], value="No")
        fixed_n = widgets.BoundedIntText(value=0, min=0, max=8)
        fixed_box = widgets.Output()

        display(fixed_yes)
        display(fixed_n)
        display(fixed_box)

        def fixed_update(*args):
            global fixed_rows
            fixed_rows = []

            with fixed_box:
                clear_output()
                if fixed_yes.value == "Yes":
                    for i in range(fixed_n.value):
                        m = widgets.Dropdown(options=names)
                        c = widgets.Dropdown(options=CHANNELS)
                        fixed_rows.append((m, c))
                        display(widgets.HBox([m, c]))

        fixed_yes.observe(fixed_update, names="value")
        fixed_n.observe(fixed_update, names="value")
        fixed_update()

        logic_state["fixed_yes"] = fixed_yes

        display(widgets.HTML("<hr><b>Checkpoint pairs</b>"))
        pair_n = widgets.BoundedIntText(value=0, min=0, max=10)
        pair_box = widgets.Output()

        display(pair_n)
        display(pair_box)

        def pair_update(*args):
            global checkpoint_pairs
            checkpoint_pairs = []

            with pair_box:
                clear_output()
                for i in range(pair_n.value):
                    a = widgets.Dropdown(options=names)
                    b = widgets.Dropdown(options=names)
                    checkpoint_pairs.append((a, b))
                    display(widgets.HBox([a, b]))

        pair_n.observe(pair_update, names="value")
        pair_update()

        display(widgets.HTML("<hr><b>Morphology different pairs</b>"))
        morph_n = widgets.BoundedIntText(value=0, min=0, max=10)
        morph_box = widgets.Output()

        display(morph_n)
        display(morph_box)

        def morph_update(*args):
            global morph_pairs
            morph_pairs = []

            with morph_box:
                clear_output()
                for i in range(morph_n.value):
                    a = widgets.Dropdown(options=names)
                    b = widgets.Dropdown(options=names)
                    morph_pairs.append((a, b))
                    display(widgets.HBox([a, b]))

        morph_n.observe(morph_update, names="value")
        morph_update()


# =========================
# SUGGEST
# =========================
def suggest_panel(_):
    df = collect_df()

    seg = get_seg_marker()
    fixed = get_fixed()

    reserved = []
    reserved_channels = []

    # segmentation marker fixed to 480
    if seg in df["marker"].values:
        row = df[df["marker"] == seg].copy()
        row["channel"] = 480
        reserved.append(row)
        reserved_channels.append(480)
        df = df[df["marker"] != seg]

    # fixed channels
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
        result = pd.concat(reserved, ignore_index=True) if reserved else pd.DataFrame()
        with result_out:
            clear_output()
            display(result)
        return

    if len(pool) < n:
        with result_out:
            clear_output()
            print("Not enough available channels for the number of unfixed markers.")
        return

    best = None
    best_tuple = None

    for combo in itertools.combinations(pool, n):
        if not spacing_rule_ok(combo):
            continue

        for perm in itertools.permutations(combo, n):
            # NEW: hard checkpoint filter
            if not checkpoint_rule_ok(perm, df):
                continue

            risk = total_risk(df, perm)
            spread = spread_penalty(perm)
            late_pen = late_channel_penalty(perm, n)

            score_tuple = (late_pen, risk, spread)

            if best_tuple is None or score_tuple < best_tuple:
                best_tuple = score_tuple
                best = perm

    # fallback: relax spacing rule, but keep checkpoint hard rule
    if best is None:
        for combo in itertools.combinations(pool, n):
            for perm in itertools.permutations(combo, n):
                if not checkpoint_rule_ok(perm, df):
                    continue

                risk = total_risk(df, perm)
                spread = spread_penalty(perm)
                late_pen = late_channel_penalty(perm, n)

                score_tuple = (late_pen, risk, spread)

                if best_tuple is None or score_tuple < best_tuple:
                    best_tuple = score_tuple
                    best = perm

    if best is None:
        with result_out:
            clear_output()
            print("No valid panel found under current constraints.")
        return

    df["channel"] = best
    result = pd.concat(reserved + [df], ignore_index=True)

    with result_out:
        clear_output()
        display(result)


generate_button.on_click(build_form)
logic_button.on_click(build_logic)
suggest_button.on_click(suggest_panel)

display(
    widgets.VBox([
        title,
        marker_count,
        generate_button,
        form_out,
        logic_button,
        logic_out,
        suggest_button,
        result_out
    ])
)

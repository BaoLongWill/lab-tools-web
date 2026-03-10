import pandas as pd
import streamlit as st

st.title("Cell-to-Gene Marker Finder")

st.write("Find marker genes for a selected cell type by species.")

MARKER_FILE = "/data/Cell_marker_All.xlsx"


@st.cache_data
def load_marker_data():
    df = pd.read_excel(MARKER_FILE)

    for col in ["species", "cell_name", "cell_type", "Symbol"]:
        df[col] = df[col].fillna("").astype(str).str.strip()

    df = df[
        (df["species"] != "") &
        (df["Symbol"] != "")
    ].copy()

    df["cell_label"] = df["cell_name"]
    df.loc[df["cell_label"] == "", "cell_label"] = df.loc[df["cell_label"] == "", "cell_type"]

    df = df[df["cell_label"] != ""].copy()

    cell_options = sorted(df["cell_label"].unique().tolist())

    return df, cell_options


def find_genes(df, selected_cell, selected_species):
    if selected_cell == "No match found":
        return []

    result = df.copy()

    result = result[result["species"].str.lower() == selected_species.lower()]

    exact_name = result[result["cell_name"].str.lower() == selected_cell.lower()]

    if len(exact_name) > 0:
        result = exact_name
    else:
        result = result[result["cell_type"].str.lower() == selected_cell.lower()]

    genes = sorted(result["Symbol"].dropna().astype(str).str.strip().unique().tolist())

    return [g for g in genes if g != ""]


try:
    df, cell_options = load_marker_data()

    st.success(
        f"Loaded {len(df)} rows across species: "
        + ", ".join(sorted(df['species'].dropna().astype(str).unique().tolist()))
    )

    search_text = st.text_input(
        "Search cell type",
        placeholder="Type here if you cannot find the cell in the dropdown"
    )

    if search_text.strip() == "":
        filtered_options = cell_options
    else:
        keyword = search_text.strip().lower()
        filtered_options = [x for x in cell_options if keyword in x.lower()]

    if len(filtered_options) == 0:
        filtered_options = ["No match found"]

    selected_cell = st.selectbox("Choose cell", filtered_options)
    selected_species = st.selectbox("Species", ["Human", "Mouse"], index=0)

    if st.button("Find genes"):
        genes = find_genes(df, selected_cell, selected_species)

        st.write(f"**Selected cell:** {selected_cell}")
        st.write(f"**Species:** {selected_species}")

        if len(genes) == 0:
            st.info("No genes found.")
        else:
            st.subheader("Marker genes")
            st.text("\n".join(genes))

            result_df = pd.DataFrame({"Gene": genes})
            st.dataframe(result_df, use_container_width=True)

except Exception as e:
    st.error("Failed to load marker database.")
    st.code(str(e))

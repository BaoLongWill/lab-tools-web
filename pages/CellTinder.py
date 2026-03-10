import pandas as pd
import streamlit as st

st.title("CellTinder")

st.write("Find possible signaling pathways between two human cell types using curated marker and interaction databases.")

# -----------------------------
# File paths
# -----------------------------
MARKER_FILE = "/data/Cell_marker_All.xlsx"
INTERACTION_FILE = "/data/interaction_input.csv"
GENE_FILE = "/data/gene_input.csv"
COMPLEX_FILE = "/data/complex_input.csv"


@st.cache_data
def load_data():
    marker_df = pd.read_excel(MARKER_FILE)
    interaction_df = pd.read_csv(INTERACTION_FILE)
    gene_df = pd.read_csv(GENE_FILE)
    complex_df = pd.read_csv(COMPLEX_FILE)

    for col in ["species", "cell_name", "cell_type", "Symbol"]:
        marker_df[col] = marker_df[col].fillna("").astype(str).str.strip()

    marker_df = marker_df[
        (marker_df["species"] != "") &
        (marker_df["Symbol"] != "")
    ].copy()

    marker_df["cell_label"] = marker_df["cell_name"]
    marker_df.loc[
        marker_df["cell_label"] == "",
        "cell_label"
    ] = marker_df.loc[marker_df["cell_label"] == "", "cell_type"]

    marker_df = marker_df[marker_df["cell_label"] != ""].copy()
    marker_df = marker_df[marker_df["species"].str.lower() == "human"].copy()

    for col in ["uniprot", "hgnc_symbol"]:
        gene_df[col] = gene_df[col].fillna("").astype(str).str.strip()

    gene_df = gene_df[
        (gene_df["uniprot"] != "") &
        (gene_df["hgnc_symbol"] != "")
    ].copy()

    uniprot_to_gene = dict(zip(gene_df["uniprot"], gene_df["hgnc_symbol"]))

    complex_cols = ["uniprot_1", "uniprot_2", "uniprot_3", "uniprot_4", "uniprot_5"]

    for col in ["complex_name"] + complex_cols:
        complex_df[col] = complex_df[col].fillna("").astype(str).str.strip()

    complex_to_genes = {}

    for _, row in complex_df.iterrows():
        complex_name = row["complex_name"]
        genes = []

        for c in complex_cols:
            u = row[c]
            if u and u in uniprot_to_gene:
                genes.append(uniprot_to_gene[u])

        genes = sorted(list(set([g for g in genes if g != ""])))
        complex_to_genes[complex_name] = genes

    def partner_to_genes(partner):
        partner = str(partner).strip()

        if partner == "" or partner.lower() == "nan":
            return []

        if partner in uniprot_to_gene:
            return [uniprot_to_gene[partner]]

        if partner in complex_to_genes:
            return complex_to_genes[partner]

        return []

    interaction_df["partner_a"] = interaction_df["partner_a"].fillna("").astype(str).str.strip()
    interaction_df["partner_b"] = interaction_df["partner_b"].fillna("").astype(str).str.strip()
    interaction_df["classification"] = interaction_df["classification"].fillna("").astype(str).str.strip()
    interaction_df["interactors"] = interaction_df["interactors"].fillna("").astype(str).str.strip()

    interaction_df["genes_a"] = interaction_df["partner_a"].apply(partner_to_genes)
    interaction_df["genes_b"] = interaction_df["partner_b"].apply(partner_to_genes)

    interaction_clean = interaction_df[
        (interaction_df["genes_a"].apply(len) > 0) &
        (interaction_df["genes_b"].apply(len) > 0)
    ].copy()

    cell_options = sorted(marker_df["cell_label"].unique().tolist())

    return marker_df, interaction_clean, cell_options


def get_cell_markers(marker_df, cell_name):
    sub = marker_df[
        marker_df["cell_label"].str.lower() == str(cell_name).lower()
    ].copy()

    genes = sorted(
        sub["Symbol"].dropna().astype(str).str.strip().unique().tolist()
    )

    return set([g for g in genes if g != ""])


def find_possible_pathways(marker_df, interaction_clean, cell_a, cell_b):
    genes_a = get_cell_markers(marker_df, cell_a)
    genes_b = get_cell_markers(marker_df, cell_b)

    hits = []

    for _, row in interaction_clean.iterrows():
        row_genes_a = set(row["genes_a"])
        row_genes_b = set(row["genes_b"])

        forward_match = (
            len(genes_a.intersection(row_genes_a)) > 0 and
            len(genes_b.intersection(row_genes_b)) > 0
        )

        reverse_match = (
            len(genes_a.intersection(row_genes_b)) > 0 and
            len(genes_b.intersection(row_genes_a)) > 0
        )

        if forward_match or reverse_match:
            axis = f"{'/'.join(row['genes_a'])} ↔ {'/'.join(row['genes_b'])}"
            pathway = row["classification"] if row["classification"] != "" else "Unclassified"
            interactors = row["interactors"] if row["interactors"] != "" else axis

            hits.append({
                "axis": axis,
                "pathway": pathway,
                "interactors": interactors
            })

    hit_df = pd.DataFrame(hits).drop_duplicates()

    if len(hit_df) == 0:
        return hit_df, []

    pathway_list = sorted(hit_df["pathway"].dropna().unique().tolist())
    return hit_df, pathway_list


# -----------------------------
# Main app
# -----------------------------
try:
    marker_df, interaction_clean, cell_options = load_data()

    st.success(f"Loaded {len(marker_df)} human marker rows and {len(interaction_clean)} mapped interactions.")

    search_a = st.text_input("Search Cell A")
    filtered_a = [x for x in cell_options if search_a.lower() in x.lower()] if search_a else cell_options
    if not filtered_a:
        filtered_a = ["No match found"]
    cell_a = st.selectbox("Cell A", filtered_a)

    search_b = st.text_input("Search Cell B")
    filtered_b = [x for x in cell_options if search_b.lower() in x.lower()] if search_b else cell_options
    if not filtered_b:
        filtered_b = ["No match found"]
    cell_b = st.selectbox("Cell B", filtered_b)

    if st.button("Find possible pathways"):
        if cell_a == "No match found" or cell_b == "No match found":
            st.warning("Please select valid cell types.")
        else:
            hit_df, pathway_list = find_possible_pathways(
                marker_df,
                interaction_clean,
                cell_a,
                cell_b
            )

            st.write(f"**Cell A:** {cell_a}")
            st.write(f"**Cell B:** {cell_b}")

            if len(hit_df) == 0:
                st.info("No possible pathways found.")
            else:
                st.subheader("Possible pathways")
                for p in pathway_list:
                    st.write(f"- {p}")

                st.subheader("Example signaling axes")
                for axis in hit_df["axis"].head(20).tolist():
                    st.write(f"- {axis}")

                st.subheader("Detailed results")
                st.dataframe(hit_df, use_container_width=True)

except Exception as e:
    st.error("Failed to load CellTinder data.")
    st.code(str(e))

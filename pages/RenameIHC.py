import os
from datetime import datetime
import streamlit as st

st.title("EzRename IHC")

st.write("Upload IHC images and rename them using standardized format. Useful for a long day of experimenting, you might forget some details.")

save_folder = os.path.expanduser("~/Desktop/RenameIHC_Output")
os.makedirs(save_folder, exist_ok=True)

st.info(f"Renamed images will be saved to: {save_folder}")

uploaded_files = st.file_uploader(
    "Upload images",
    type=["tif", "tiff", "png", "jpg", "jpeg"],
    accept_multiple_files=True
)

marker = st.text_input("Marker")

dilution = st.selectbox(
    "Dilution",
    ["1:50", "1:100", "1:200", "1:300", "1:500", "1:800", "1:1000", "1:2000", "Other"]
)

dilution_other = ""
if dilution == "Other":
    dilution_other = st.text_input("Other dilution")

buffer = st.selectbox(
    "Buffer",
    ["Citrate", "HTTR", "EDTA", "Other"]
)

buffer_other = ""
if buffer == "Other":
    buffer_other = st.text_input("Other buffer")

magnification = st.selectbox(
    "Magnification",
    ["4X", "10X", "20X", "40X"]
)

user = st.text_input("User")

if st.button("Rename Images"):
    if not uploaded_files:
        st.warning("Upload images first.")
    else:
        m = marker.strip()

        d = dilution
        if d == "Other":
            d = dilution_other.strip()

        buf = buffer
        if buf == "Other":
            buf = buffer_other.strip()

        mag = magnification
        u = user.strip().upper()

        date_str = datetime.now().strftime("%Y%m%d")
        saved_names = []

        for i, uploaded_file in enumerate(uploaded_files, start=1):
            original_name = uploaded_file.name
            ext = os.path.splitext(original_name)[1]
            d_clean = d.replace(":", "-")

            new_name = f"{m}_{d_clean}_{buf}_{mag}_{u}_{date_str}_ROI{i:03d}{ext}"
            save_path = os.path.join(save_folder, new_name)

            with open(save_path, "wb") as f:
                f.write(uploaded_file.getbuffer())

            saved_names.append(new_name)

        st.success("Images renamed and saved successfully.")
        st.write("Saved files:")
        for name in saved_names:
            st.write(name)

        st.write(f"Saved folder: {save_folder}")

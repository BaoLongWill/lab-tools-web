import streamlit as st

st.set_page_config(page_title="Will Tools", layout="wide")

st.title("Will Lazy Tools")

st.write("A tiny toolbox for tired researchers doing spatial biology and marker-related.
It handles the boring parts of experiments and other repetitive tasks so that you can focus on the fun part.")

st.markdown("""
### Available tools

- EzRename IHC  
- mIHC panel suggestion  
- Cell-to-Gene Marker Finder
- CellTinder

Use the sidebar to select a tool.
""")

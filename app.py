import streamlit as st

st.set_page_config(page_title="Will Tools", layout="wide")

st.title("Will Lazy Tools")

st.write("LazyMe is a tiny toolbox for tired researchers doing spatial biology and IHC.
It handles the boring parts of experiments — organizing images, planning mIHC panels, and other repetitive tasks — so you can focus on the fun part: discovering biology.")

st.markdown("""
### Available tools

- EzRename IHC  
- mIHC panel suggestion  
- Cell-to-Gene Marker Finder
- CellTinder

Use the sidebar to select a tool.
""")

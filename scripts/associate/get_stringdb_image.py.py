# ---
# jupyter:
#   jupytext:
#     cell_metadata_json: true
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.5
#   kernelspec:
#     display_name: Python [conda env:anaconda-florian4]
#     language: python
#     name: conda-env-anaconda-florian4-py
# ---

# %%
from IPython.display import display, SVG, Image

# %%
import os
import numpy as np
import pandas as pd

import json
import yaml

import polars as pl


# %%
import stringdb
import requests

# %%
import plotnine as pn
# import seaborn as sns
# import networkx as nx

import matplotlib
import matplotlib.pyplot as plt

# %%
from rep.notebook_init import setup_plot_style
setup_plot_style()

# %%
# %matplotlib inline
# %config InlineBackend.figure_format='retina'

# %%
# import os
# # os.environ["RAY_ADDRESS"] = os.environ.get("RAY_ADDRESS", 'ray://192.168.16.30:10001')
# os.environ["RAY_ADDRESS"] = 'ray://192.168.16.28:10001'
# os.environ["RAY_ADDRESS"]

# %%
snakefile_path = os.getcwd() + "/../../Snakefile"
snakefile_path

# %%
#del snakemake

# %%
try:
    snakemake
except NameError:
    from snakemk_util import load_rule_args
    
    snakemake = load_rule_args(
        snakefile = snakefile_path,
        rule_name = 'associate__get_stringdb_image',
        default_wildcards={
            # "phenotype_col": "Asthma",
            # "phenotype_col": "severe_LDL",
            # "phenotype_col": "triglycerides_f30870_0_0",
            # "phenotype_col": "Aspartate_aminotransferase",
            "phenotype_col": "HDL_cholesterol",
            # "feature_set": "LOFTEE_pLoF",
            "feature_set": "AbExp_all_tissues",
            # "feature_set": "max_AbExp",
            # "covariates": "sex_age_genPC",
            # "covariates": "sex_age_genPC_CLMP_PRS",
            "covariates": "sex_age_genPC_BMI_smoking_CLMP_PRS",
        }
    )

# %%
from snakemk_util import pretty_print_snakemake
print(pretty_print_snakemake(snakemake))

# %%
if "plot_dpi" in snakemake.params:
    DPI = snakemake.params["plot_dpi"]
else:
    DPI=450

# %%
pval_cutoff = snakemake.params["pval_cutoff"]
pval_cutoff

# %% [markdown]
# # read input data

# %%
significant_genes = pl.read_parquet(snakemake.input["significant_genes_pq"]).filter(pl.col("padj") < snakemake.params["pval_cutoff"])["gene"].to_list()
significant_genes

# %%
mapped_string_ids = stringdb.get_string_ids(significant_genes, species=9606)
mapped_string_ids

# %%
network = stringdb.get_network(mapped_string_ids["stringId"])
network

# %% [raw]
# graph = nx.MultiGraph()
# for node in mapped_string_ids["preferredName"]:
#     graph.add_node(node)
#
# for idx, edge in network.iterrows():
#     graph.add_edge(edge["preferredName_A"], edge["preferredName_B"], weight=edge["score"])
# graph

# %% [raw]
# fig, ax = plt.subplots(figsize=(10,6))
# pos = nx.nx_agraph.graphviz_layout(graph)
# nx.draw(graph, pos=pos, ax=ax)
# nx.draw_networkx_labels(graph, pos=pos, ax=ax)
# nx.draw_networkx_nodes(graph, pos=pos, ax=ax)
# fig.show()

# %%
enrichment_df = (
    pl.DataFrame(stringdb.get_enrichment(mapped_string_ids["stringId"]))
    .with_columns([
        pl.col("inputGenes").str.split(","),
        pl.col("preferredNames").str.split(","),
    ])
)
enrichment_df

# %%
enrichment_df.write_parquet(snakemake.output["stringdb_enrichment_pq"])


# %% [markdown]
# # download StringDB network image

# %%
def download_stringdb_svg(string_ids, add_color_nodes=None, add_white_nodes=None):
    import requests

    # Download the SVG file using requests
    url = "https://string-db.org/api/svg/network"
    params = {
        "identifiers":  "%0d".join(mapped_string_ids["stringId"]),
        # "add_color_nodes": "10",
        "network_type": "physical",
        "species": "9606",
    }
    if add_color_nodes is not None:
        params["add_color_nodes"] = add_color_nodes
    if add_white_nodes is not None:
        params["add_white_nodes"] = add_white_nodes

    response = requests.get(url, params=params)
    assert response.ok, "Request failed: "
    svg_content = response.content
    
    return svg_content.decode("UTF-8")


# %%
svg_content = download_stringdb_svg(mapped_string_ids["stringId"].unique(), add_white_nodes=len(mapped_string_ids["stringId"].unique()))

# %%
from IPython.display import SVG
display(SVG(svg_content))


# %%
def get_node_colors(svg):
    import xml.etree.ElementTree as ET
    
    retval = {}
    
    # Parse the SVG file
    root = ET.fromstring(svg)

    for node in root.findall(".//{http://www.w3.org/2000/svg}g"):
        if node.get("class") == "nwnodecontainer":
            node_label = node.get("data-safe_div_label")
            for circle in node.findall(".//{http://www.w3.org/2000/svg}circle"):
                if circle.get("class") == "nwbubblecoloredcircle":
                    retval[node_label] = circle.get("fill")
    
    return retval


# %%
get_node_colors(svg_content)

# %%
import math

def get_point_on_circle(x, y, r, fraction=0.0):
    t = 2*math.pi * fraction
    adj_x = r*math.sin(t) + x
    adj_y = -r*math.cos(t) + y
    return (adj_x, adj_y)

def set_color(node, node_colors, default='rgb(255,255,255)'):
    import xml.etree.ElementTree as ET
    print(f"setting {node} to color(s) {node_colors}. Len: {len(node_colors)}")
    
    for circle in node.findall(".//{http://www.w3.org/2000/svg}circle"):
        if circle.get("class") != "nwbubblecoloredcircle":
            continue
        
        if len(node_colors) == 0:
            circle.set("fill", default)
        elif len(node_colors) == 1:
            # set the node color
            circle.set("fill", node_colors[0])
        else:
            circle.set("display", "none")
            r = float(node.get("data-radius"))
            x = float(node.get("data-x_pos"))
            y = float(node.get("data-y_pos"))
            
            parts = len(node_colors)
            for i in range(parts):
                fraction1 = 1/parts * i
                fraction2 = 1/parts * (i + 1)

                c_x1, c_y1 = get_point_on_circle(x, y, r, fraction1)
                c_x2, c_y2 = get_point_on_circle(x, y, r, fraction2)
                
                color = node_colors[i]
                if color is None:
                    color = default

                svg_path = ET.fromstring(
                    f'''<path xmlns="http://www.w3.org/2000/svg" d="
                        M {c_x1} {c_y1}
                        A {r} {r} 0 0 1 {c_x2} {c_y2} 
                        L {x} {y}
                    " fill="{color}" opacity="0.4"/>'''
                )
                
                node.append(svg_path)
    return node



# %%
def set_node_colors(svg, node_colors, default='rgb(255,255,255)'):
    import xml.etree.ElementTree as ET
    
    # Parse the SVG file
    root = ET.fromstring(svg)

    for node in root.findall(".//{http://www.w3.org/2000/svg}g"):
        if node.get("class") == "nwnodecontainer":
            node_label = node.get("data-safe_div_label")
            for circle in node.findall(".//{http://www.w3.org/2000/svg}circle"):
                if circle.get("class") == "nwbubblecoloredcircle":
                    # set the node color
                    if node_label in node_colors:
                        set_color(node, node_colors[node_label], default=default)
                    else:
                        circle.set("fill", default)
    
    edited_svg_content = ET.tostring(root, encoding='unicode')
    return edited_svg_content


# %%
def color_nodes(mapping, by=None, node_col="preferredName", palette=['rgb(255,0,0)', 'rgb(0,0,255)', "green"]):
    from collections import defaultdict

    node_colors = defaultdict(list)
    if by is not None and len(by) > 0:
        for idx, (key, df) in enumerate(mapped_string_ids.groupby(by)):
            color = palette[idx]
            nodes = df[node_col]
            for node in nodes:
                node_colors[node].append(color)
    else:
        for node in mapping[node_col]:
            node_colors[node].append(palette[0])
    
    return node_colors


# %%
node_colors = {
    "ABCA1": ['rgb(255,0,0)', 'rgb(0,0,255)', None, "green", None],
    "NUP100": [],
    "NUP159": [],
}
node_colors = color_nodes(mapped_string_ids, by=None)
node_colors

# %%
edited_svg_content = set_node_colors(svg_content, node_colors)
# print(edited_svg_content)

# %%
from IPython.display import SVG
display(SVG(
    edited_svg_content
))

# %%
# Save the image to file
path = snakemake.output["string_svg"]
with open(path, 'w') as fd:
    fd.write(edited_svg_content)


# %%
dict(snakemake.output)

# %%
import wand.image

with wand.image.Image(blob=edited_svg_content.encode(), format='svg', resolution=450) as img:
    # to save output to a file:
    with img.convert('png') as output_img:
        # output_img.save(filename=output_png_path)
        # or, to get the output data in a variable:
        png_data = img.make_blob(format='png')

# %%
from IPython.display import Image
display(Image(png_data))

# %%
# Save the image to file
path = snakemake.output["string_png"]
with open(path, 'wb') as fd:
    fd.write(png_data)

# %%

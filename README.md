App4PubMedinR
A web-based bibliometric analysis tool for PubMed data using FLCA–MajorSampling, AAC, and multi-domain visualization

Overview
App4PubMedinR is an R/Shiny application for automated bibliometric analysis using PubMed data. It enables users to retrieve records via PubMed queries or upload MEDLINE .txt files, and then generate structured author profiles, domain-level summaries, and visual network analytics.
The system integrates:


PubMed data retrieval and parsing


Metadata normalization


Network construction (node–edge)


FLCA–MajorSampling for Top-20 structure extraction


AAC (Absolute Advantage Coefficient) for dominance analysis


Multi-panel visualization and reporting


This tool is designed to provide a rapid yet reproducible overview of publication patterns across multiple domains.

Key Features


🔍 PubMed query-based analysis


Supports standard tags: [Author], [Affiliation], [AD], [TA], [MH], [Title/Abstract], [PDAT]




📂 Upload MEDLINE .txt files


Enables reproducible offline analysis




🧠 FLCA–MajorSampling algorithm


Extracts Top-20 nodes and one-link edges


Reduces network complexity while preserving structure




📊 AAC-based dominance analysis


Quantifies leadership strength across entities




🌐 Multi-domain analysis


Author


Journal / Year


Country


State / Province


Institute


Department


MeSH




📈 Visualization suite


Interactive network graphs (visNetwork)


SSplot (cluster structure + dominance)


Kano plots (value vs value2, ss vs a*)


Sankey diagrams (leader–follower flows)


Geographic maps (global + regional)




📦 Export & reporting


Download nodes / edges (CSV)


Export PNG figures


Generate HTML reports


ZIP package for reproducibility





System Architecture
The application is organized into five modules:


Input module


PubMed query or MEDLINE upload




Parsing & normalization


Extracts and standardizes metadata




Network construction


Builds node–edge structures across domains




Analytical engine


FLCA–MajorSampling


Silhouette (SS)


AAC computation




Visualization & export


Interactive plots and downloadable reports





Workflow
Step 1 — Input


Enter a PubMed query OR


Upload a MEDLINE .txt file


Step 2 — Fetch / Parse


Retrieve PMIDs and metadata


Parse authors, affiliations, MeSH, etc.


Step 3 — Run domain analysis


Author / Country / Institute / MeSH / etc.


Step 4 — Visualization


Network, SSplot, Kano, Sankey, Maps


Step 5 — Export


CSV, PNG, HTML report, ZIP package



Installation
Requirements


R (≥ 4.2 recommended)


RStudio (optional)


Required packages
install.packages(c(  "shiny", "rentrez", "DT", "visNetwork",  "igraph", "ggplot2", "maps", "htmlwidgets"))

Run the App
shiny::runApp("path/to/App4PubMedinR")
Or inside the project directory:
shiny::runApp()

Example Query
(Tsair-Wei Chien[Author]) AND (Taiwan[Affiliation])
Other examples:


"Artificial Intelligence"[MH] AND Taiwan[AD]


"machine learning"[Title/Abstract] AND 2020:2025[PDAT]



Core Methods
FLCA–MajorSampling


Identifies leader–follower structures


Extracts Top-20 representative nodes


Uses one-link edges to simplify visualization



AAC (Absolute Advantage Coefficient)
r=v1/v2v2/v3r = \frac{v_1 / v_2}{v_2 / v_3}r=v2​/v3​v1​/v2​​
AAC=r1+rAAC = \frac{r}{1 + r}AAC=1+rr​


v1,v2,v3v_1, v_2, v_3v1​,v2​,v3​: top 3 ranked values


Higher AAC → stronger dominance



SSplot


Visualizes:


Cluster structure


Silhouette scores


Dominance patterns





Output


Interactive dashboards (Shiny UI)


Tables:


Nodes (Top-20)


Edges (leader–follower)




Figures:


Network


SSplot


Kano


Sankey


Maps




Downloads:


CSV files


PNG images


HTML report


ZIP archive





Deployment
Local
Run via RStudio or R console
Cloud (shinyapps.io)


Ensure:


No browseURL() calls


Use tempdir() for writable paths




Deploy via:


rsconnect::deployApp()

Repository Structure
App4PubMedinR/│── app.R│── renderSSplot.R│── kano.R│── sankey.R│── pubmed_parse_biblio.R│── appAAC.R│── flca_ms_sil_module.R│── helper_ss_patch.R│── ipmodule.R│── README.md

Use Cases


Bibliometric research


Author collaboration analysis


Institutional performance analysis


MeSH-based topic exploration


Rapid manuscript background review



Citation
If you use this software, please cite:

App4PubMedinR: A web-based bibliometric analysis tool for PubMed using FLCA–MajorSampling and AAC.


Contact


📧 rasch.smile@gmail.com


📧 codingpaperabc@gmail.com



License
Specify your license here (e.g., MIT, GPL-3.0)

Acknowledgments


PubMed / NCBI API


R Shiny ecosystem


FLCA–MajorSampling framework


AAC dominance analysis



If you want, I can next:


tailor this README to SoftwareX submission style


or generate a downloadable README.md file


or add badges + screenshots + demo GIFs for GitHub polish


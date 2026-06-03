# App4PubMedinR

**App4PubMedinR** is an R/Shiny web application for PubMed-oriented bibliometric analysis, author collaboration profiling, and dominance visualization using the **Absolute Advantage Coefficient (AAC)** and **FLCA–MajorSampling**.

The app supports three practical workflows:

1. **Two-step PubMed fetching** through a PubMed query or uploaded MEDLINE `.txt` file.  
2. **Two-step AMA/PubMed, Google Scholar, or NCKU Pure/Scopus extraction** through upload or pasted normalized text.  
3. **Real-time visualization through a PubMed URL link**, allowing users to generate visual outputs directly from a PubMed search URL.

---

## 1. Two-Step Fetch PubMed Workflow

This workflow is used when users want to retrieve PubMed records directly and generate bibliometric dashboards.

### Step 1. Enter a PubMed query or upload a MEDLINE `.txt` file

Users can either enter a PubMed query, for example:

```text
(Tsair-Wei Chien[Author]) AND (Taiwan[Affiliation])
```

or upload a PubMed MEDLINE `.txt` file downloaded from PubMed.

### Step 2. Click **Fetch PubMed**

After clicking **Fetch PubMed**, the app retrieves or parses PubMed records and generates structured outputs, including:

- Author collaboration networks
- First/last-author analysis
- Journal, country, institution, department, year, and MeSH summaries
- AAC-based dominance indicators
- Network, SSplot, Kano, Sankey, and map visualizations
- Downloadable CSV, PNG, HTML, and ZIP reports

---

## 2. Two-Step AMA / Google Scholar / NCKU Pure-Scopus Extraction Workflow

This workflow is used when users already have formatted reference records or author lists.

### Step 1. Upload a reference file or paste records into the textarea box

Users can paste or upload records from:

- AMA/PubMed reference lists
- Google Scholar copied records
- NCKU Pure / Scopus normalized research-output records
- Semicolon-separated author lists

Example NCKU Pure / Scopus normalized record:

```text
Peer-assisted learning in critical care: a simulation-based approach for postgraduate medical training
Chiu, P. W., Chu, S. C., Yang, C. H., Lee, H. F., Hung, H. M. & Hsu, H. C., 2025, 於: Medical Education Online. 30, 1, 2497333.
```

Example semicolon-separated author list:

```text
Zeng J.; Ning X.; Lan L.; Peng W.; Zheng J.; Yang K.; Luo H.
```

### Step 2. Select the parser and run extraction

Select the appropriate parser checkbox, such as:

- **Parse AMA/PubMed reference text**
- **Parse Google Scholar normalized rows**
- **Parse NCKU Pure / Scopus research output normalized records**

Then click one of the extraction buttons:

- **Run AMA journal extraction**
- **Run AMA author/journal extraction**
- **Run AMA author AAC (1st+Last only)**

The app produces author/journal frequency tables, first/last-author networks, AAC summaries, h-index tables when citation counts are available, and downloadable results.

---

## 3. Real-Time Visuals via PubMed URL Link

App4PubMedinR can generate real-time visual analytics directly from a PubMed search URL.

### Example PubMed URL

```text
https://pubmed.ncbi.nlm.nih.gov/?term=%28Tsair-Wei+Chien%5BAuthor%5D%29+AND+%28Taiwan%5BAffiliation%5D%29
```

### Workflow

1. Copy a PubMed search URL.
2. Paste the URL into the URL-link input box in App4PubMedinR.
3. Click **Fetch PubMed**.
4. The app immediately retrieves PubMed records and generates interactive visual outputs.

This URL-based workflow is useful for rapid author-profile analysis, collaboration mapping, and dominance visualization without manually downloading files.

---

## Key Features

- PubMed query-based retrieval
- PubMed URL-link retrieval
- MEDLINE `.txt` upload
- AMA/PubMed reference extraction
- Google Scholar normalized-record extraction
- NCKU Pure / Scopus normalized-record extraction
- Semicolon-separated author-list parsing
- First/last-author collaboration analysis
- FLCA–MajorSampling Top-20 node extraction
- AAC-based dominance analysis
- Article-level and author-level citation metrics when citation counts are available
- Multi-domain visualizations for authors, journals, countries, institutions, departments, publication years, MeSH terms, U.S. states, and Chinese provinces
- Downloadable CSV, PNG, HTML, and ZIP outputs

---

## Core Methods

### FLCA–MajorSampling

FLCA–MajorSampling identifies leader–follower structures and extracts representative Top-20 nodes with one-link network structures.

### Absolute Advantage Coefficient

AAC quantifies dominance strength among ranked entities.

```text
r = (v1/v2) / (v2/v3)
AAC = r / (1 + r)
```

AAC helps classify whether a leading author, journal, country, institution, department, or MeSH term shows low, moderate, high, or very high dominance.

---

## Installation

Install the required R packages:

```r
install.packages(c(
  "shiny", "rentrez", "DT", "visNetwork", "igraph", "ggplot2",
  "maps", "htmlwidgets", "dplyr", "tibble", "tidyr", "readr",
  "stringr", "cluster", "boot", "htmltools", "plotly",
  "networkD3", "zip", "clue", "gridExtra", "ggraph",
  "uwot", "reticulate", "chorddiag", "circlize"
))
```

Optional packages may be required for some geographic or advanced visualization functions, depending on the deployment environment.

---

## Run the App Locally

Place all required source files in the same folder, then run:

```r
shiny::runApp()
```

---

## Deployment

The app can be deployed to Shinyapps.io:

```r
rsconnect::deployApp()
```

A PubMed API key can optionally be saved in `.Renviron`:

```text
PUBMED_API_KEY=your_ncbi_api_key
```

---

## Repository Structure

```text
app.R
appAAC.R
pubmed_parse_biblio.R
renderSSplot.R
kano.R
sankey.R
flca_ms_sil_module.R
flca_ma_sil_module.R
helper_ss_patch.R
ipmodule.R
README.md
```

---

## Outputs

App4PubMedinR produces:

- Interactive networks
- SSplots
- Kano plots
- Sankey diagrams
- Geographic maps
- Slope and trend plots
- Author, journal, country, institution, department, year, and MeSH tables
- First/last-author AAC reports
- Citation and h-index summaries when citation counts are included
- Downloadable CSV, PNG, HTML, and ZIP reports

---

## Contact

- rasch.smile@gmail.com
- codingpaperabc@gmail.com

---

## License

Please specify the repository license, such as MIT, GPL-3.0, or another license appropriate for your project.

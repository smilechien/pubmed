# App4PubMedinR

A web-based bibliometric analysis tool for PubMed data using FLCA–MajorSampling, AAC, and multi-domain visualization

## Overview

App4PubMedinR is an R/Shiny application for automated bibliometric analysis using PubMed data. It enables users to retrieve records via PubMed queries or upload MEDLINE `.txt` files, and then generate structured author profiles, domain-level summaries, and visual network analytics.

The system integrates:
- PubMed data retrieval and parsing  
- Metadata normalization  
- Network construction (node–edge)  
- FLCA–MajorSampling for Top-20 structure extraction  
- AAC (Absolute Advantage Coefficient) for dominance analysis  
- Multi-panel visualization and reporting  

## Key Features

- PubMed query-based analysis  
- Upload MEDLINE `.txt` files  
- FLCA–MajorSampling algorithm (Top-20 nodes & one-link edges)  
- AAC-based dominance analysis  
- Multi-domain analysis (Author, Journal, Country, Institute, Department, MeSH)  
- Visualization: Network, SSplot, Kano, Sankey, Maps  
- Export: CSV, PNG, HTML report, ZIP  

## Installation

```r
install.packages(c(
  "shiny", "rentrez", "DT", "visNetwork",
  "igraph", "ggplot2", "maps", "htmlwidgets"
))
```

## Run the App

```r
shiny::runApp()
```

## Example Query

```
(Tsair-Wei Chien[Author]) AND (Taiwan[Affiliation])
```

## Core Methods

### FLCA–MajorSampling
Identifies leader–follower structures and extracts Top-20 representative nodes.

### AAC (Absolute Advantage Coefficient)

r = (v1/v2)/(v2/v3)  
AAC = r / (1 + r)

## Output

- Interactive dashboards  
- Nodes and edges tables  
- Network, SSplot, Kano, Sankey, Maps  
- Downloadable CSV, PNG, HTML, ZIP  

## Deployment

Use shinyapps.io:
```r
rsconnect::deployApp()
```

## Repository Structure

- app.R  
- renderSSplot.R  
- kano.R  
- sankey.R  
- pubmed_parse_biblio.R  
- appAAC.R  
- flca_ms_sil_module.R  
- helper_ss_patch.R  
- ipmodule.R  

## Contact

- rasch.smile@gmail.com  
- codingpaperabc@gmail.com  

## License

Specify your license (e.g., MIT)



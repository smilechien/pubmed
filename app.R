# PATCH_2026_06_02_SEPARATE_JOURNAL_ONLY_AND_MIXED_AUTHOR_JOURNAL
# VERIFIED_NO_STANDALONE_BURST_TABS_2026_05_28
# PATCHED_FROM_USER_UPLOADED_VERIFIED2_WITH_SLOPE_3D_DASHBOARD_HEATMAP3_2026_05_28
# PATCHED_AUTHOR_BYLINE_COUNTS_ANY_BYLINE_AND_FIRST_LAST_2026_05_30
if (file.exists(".Renviron")) {
  readRenviron(".Renviron")
}

# ---- IP module (shared gate) ----
APP_DIR <- tryCatch(normalizePath(getwd(), winslash = "/", mustWork = FALSE), error = function(e) getwd())

# ---- SAFE: source plotting modules if present (stable app keeps homepage unchanged) ----
.safe_source_if_exists <- function(paths) {
  for (p in paths) {
    if (file.exists(p)) {
      try(source(p, local = FALSE, encoding = "UTF-8"), silent = TRUE)
      return(invisible(TRUE))
    }
  }
  invisible(FALSE)
}
.safe_source_if_exists(c(file.path(APP_DIR, "renderSSplot(81).R"), file.path(APP_DIR, "renderSSplot(80).R"), file.path(APP_DIR, "renderSSplot(79).R"), file.path(APP_DIR, "renderSSplot.R"), file.path(APP_DIR, "renderSSplot(78).R")))
.safe_source_if_exists(c(file.path(APP_DIR, "kano(63).R"), file.path(APP_DIR, "kano.R"), file.path(APP_DIR, "kano(62).R")))
.safe_source_if_exists(c(file.path(APP_DIR, "flca_ms_sil_module.R"), file.path(APP_DIR, "flca_ma_sil_module.R")))
if (exists("run_flca_ms_sil_runner", mode = "function")) run_flca_ms_sil <- run_flca_ms_sil_runner


# ---- Load SSplot/Kano modules (support both stable and uploaded filenames) ----
for (.ss_file in c("renderSSplot(81).R", "renderSSplot(80).R", "renderSSplot(79).R", "renderSSplot.R", "renderSSplot(78).R")) {
  .ss_path <- file.path(APP_DIR, .ss_file)
  if (file.exists(.ss_path)) try(source(.ss_path, local = FALSE, encoding = "UTF-8"), silent = TRUE)
}
for (.kano_file in c("kano(63).R", "kano.R", "kano(62).R")) {
  .kano_path <- file.path(APP_DIR, .kano_file)
  if (file.exists(.kano_path)) try(source(.kano_path, local = FALSE, encoding = "UTF-8"), silent = TRUE)
}
# ---- Load FLCA-MA-SIL module if present ----
for (.flca_file in c("flca_ms_sil_module.R", "flca_ma_sil_module.R")) {
  .flca_path <- file.path(APP_DIR, .flca_file)
  if (file.exists(.flca_path)) try(source(.flca_path, local = FALSE, encoding = "UTF-8"), silent = TRUE)
}
if (exists("run_flca_ms_sil_runner", mode = "function")) run_flca_ms_sil <- run_flca_ms_sil_runner
if (exists("run_flca_ma_sil_runner", mode = "function")) run_flca_ms_sil <- run_flca_ma_sil_runner

if (file.exists(file.path(APP_DIR, "ipmodule.R"))) {
  source(file.path(APP_DIR, "ipmodule.R"), local = FALSE)
}
PERM_PUBMED_XML <- file.path(tempdir(), "pubmed.xml") 
PERM_PUBMED_MEDLINE <- file.path(tempdir(), "pubmed_medline.txt") 
PERM_PUBMED_TXT <- file.path(tempdir(), "uploaded_pubmed.txt")
options(shiny.error = function() {
  message("SHINY ERROR: ", geterrmessage())
})
# ----------------------------
# Helpers: AAC (Top1 / value2_strength)
# ----------------------------
.calc_aac_from_top3 <- function(x){
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x) & !is.na(x)]
  x <- sort(x, decreasing = TRUE)
  if (length(x) < 3) return(NA_real_)
  v1 <- x[1]; v2 <- x[2]; v3 <- x[3]
  if (v2 <= 0 || v3 <= 0) return(NA_real_)
  r <- (v1 * v3) / (v2 * v2)  # (v1/v2)/(v2/v3)
  if (!is.finite(r)) return(NA_real_)
  r / (1 + r)
}

source(file.path(APP_DIR, "appAAC.R"), local = FALSE)
# ==== PubMed MEDLINE TXT parser (added) ====
parse_pmids_from_medline_txt <- function(path){
  if (is.null(path) || !nzchar(path) || !file.exists(path)) return(character())
  x <- readLines(path, warn = FALSE, encoding = "UTF-8")
  if (!length(x)) return(character())
  pm <- x[grepl("^PMID-\\s*", x)]
  pm <- sub("^PMID-\\s*", "", pm)
  pm <- trimws(pm)
  pm <- pm[nzchar(pm)]
  unique(pm)
}

PUBMED_API_KEY <- Sys.getenv("PUBMED_API_KEY")

# Permanent path for uploaded MEDLINE txt (prevents Shiny temp expiry)
PERM_PUBMED_TXT <- file.path(getwd(), "uploaded_pubmed.txt")

# --- FIX: wrap zip logic into a function (prevents unmatched braces) ---
.webzip_zip_dir <- function(download_dir, zipfile) {
  old <- getwd()
  on.exit(setwd(old), add = TRUE)

  setwd(download_dir)
  utils::zip(
    zipfile = zipfile,
    files   = list.files(".", recursive = TRUE, all.files = TRUE, no.. = TRUE)
  )
  invisible(zipfile)
}

# ================================
# WebZIP helpers (offline package)

# ================================
# WebZIP figure export (AUTO)
# ================================
.webzip_dl_dir <- function() .get_dl_dir()

.webzip_save_png <- function(path, plot_fun, width=1600, height=1000, res=150){
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  grDevices::png(path, width=width, height=height, res=res)
  on.exit(grDevices::dev.off(), add=TRUE)
  try(plot_fun(), silent=TRUE)
  invisible(path)
}

.webzip_export_domain_figs <- function(domain, nodes, edges, rv, dl_dir){
  if (is.null(nodes) || !is.data.frame(nodes) || nrow(nodes)==0) return(invisible(FALSE))
  dom_safe <- gsub("[^A-Za-z0-9_]+","_", domain)

  # 1) Kano: value vs value2 (prefer pretty kano.R functions)
  if (exists("kano_plot", mode="function")) {
    .webzip_save_png(file.path(dl_dir, sprintf("kano_%s_value_vs_value2.png", dom_safe)), function(){
      kano_plot(nodes, edges, title_txt = sprintf("Kano (%s): value vs value2", domain))
    })
  } else {
    .webzip_save_png(file.path(dl_dir, sprintf("kano_%s_value_vs_value2.png", dom_safe)), function(){
      plot(nodes$value2, nodes$value, xlab="value2", ylab="value", main=sprintf("Kano (%s): value vs value2", domain))
      if ("name" %in% names(nodes)) text(nodes$value2, nodes$value, labels=nodes$name, pos=3, cex=0.7)
      abline(h=median(nodes$value, na.rm=TRUE), v=median(nodes$value2, na.rm=TRUE), lty=2)
    })
  }

  # 2) Kano: ss vs a* (if columns exist)
  if (all(c("ssi","a_star1") %in% names(nodes))) {
    if (exists("kano_plot_ss_astar", mode="function")) {
      .webzip_save_png(file.path(dl_dir, sprintf("kano_%s_ss_vs_astar.png", dom_safe)), function(){
        kano_plot_ss_astar(nodes, edges, title_txt = sprintf("Kano (%s): ss vs a*", domain))
      })
    } else {
      .webzip_save_png(file.path(dl_dir, sprintf("kano_%s_ss_vs_astar.png", dom_safe)), function(){
        plot(nodes$a_star1, nodes$ssi, xlab="a*", ylab="ss", main=sprintf("Kano (%s): ss vs a*", domain))
        if ("name" %in% names(nodes)) text(nodes$a_star1, nodes$ssi, labels=nodes$name, pos=3, cex=0.7)
        abline(h=median(nodes$ssi, na.rm=TRUE), v=median(nodes$a_star1, na.rm=TRUE), lty=2)
      })
    }
  }

  # 3) SSplot (simple scatter)
  .webzip_save_png(file.path(dl_dir, sprintf("ssplot_%s.png", dom_safe)), function(){
    x <- nodes$value2; y <- nodes$value
    plot(x, y, xlab="value2", ylab="value", main=sprintf("SSplot (%s): value vs value2", domain))
    if ("name" %in% names(nodes)) text(x, y, labels=nodes$name, pos=3, cex=0.7)
  })

  # 4) Slopegraph (Top20, recent 10y) - simple line plot if year-count helpers exist
  if (exists("compute_combo_year_counts", mode="function")) {
    try({
      yc <- compute_combo_year_counts(rv, domains = c(domain), recent_n_years = 10)
      if (!is.null(yc) && nrow(yc)) {
        top_items <- nodes$name
        yc <- yc[yc$domain == domain & yc$item %in% top_items, , drop=FALSE]
        if (nrow(yc)) {
          yrs <- sort(unique(as.integer(as.character(yc$Year))))
          items <- unique(yc$item)
          .webzip_save_png(file.path(dl_dir, sprintf("slope_%s_top20.png", dom_safe)), function(){
            plot(range(yrs), c(0, max(yc$Count, na.rm=TRUE)), type="n",
                 xlab="Year", ylab="Count", main=sprintf("Slopegraph (%s): Top20 recent 10y", domain))
            for (it in items) {
              d <- yc[yc$item==it, , drop=FALSE]
              d <- d[order(as.integer(as.character(d$Year))), , drop=FALSE]
              lines(as.integer(as.character(d$Year)), d$Count)
            }
            legend("topleft", legend=utils::head(items, 10), cex=0.6, bty="n")
          })
        }
      }
    }, silent=TRUE)
  }

  # Sankey: placeholder (your sankey.R is placeholder)
  .webzip_save_png(file.path(dl_dir, sprintf("sankey_%s.png", dom_safe)), function(){
    plot.new()
    text(0.5,0.5, sprintf("Sankey (%s): sankey.R is placeholder\\nUpload real sankey.R to export Sankey.", domain))
  })

  invisible(TRUE)
}

# - Serves download/ as static path
# - Generates download/index.html for offline viewing
# - Zips download/ into a .zip
# ================================
.get_dl_dir <- function() file.path(getwd(), "download")

.list_download_files <- function(dl_dir) {
  if (!dir.exists(dl_dir)) return(character(0))
  list.files(dl_dir, recursive = TRUE, full.names = FALSE, all.files = TRUE, no.. = TRUE)
}

# ===== FINAL WALK-IN FIX: NEVER-NA AAC in Report =====
safe_AAC <- function(x){
  if (is.null(x)) return(0)
  x <- x[is.finite(x)]
  if (!length(x)) return(0)
  .compute_AAC(x)
}


# ============================================================
# WebZIP helpers: offline index.html + zip download/
# ============================================================
.write_offline_index <- function(download_dir = "download", rv = NULL, title="Author profile analysis (offline)", ...){
  if (!dir.exists(download_dir)) dir.create(download_dir, recursive = TRUE, showWarnings = FALSE)

  term <- ""
  hits <- ""
  if (!is.null(rv)) {
    if (!is.null(rv$term)) term <- as.character(rv$term)
    if (!is.null(rv$hits)) hits <- as.character(rv$hits)
  }

  files <- list.files(download_dir, recursive = TRUE, full.names = FALSE)
  imgs  <- files[grepl("\\.(png|jpg|jpeg|gif|svg)$", files, ignore.case = TRUE)]
  tabs  <- files[grepl("\\.(csv|tsv|txt|html)$", files, ignore.case = TRUE) & !grepl("index\\.html$", files, ignore.case=TRUE)]

  esc <- function(x) htmltools::htmlEscape(x)

  lines <- c(
    "<!doctype html>",
    "<html><head><meta charset='utf-8'/>",
    "<meta name='viewport' content='width=device-width, initial-scale=1'/>",
    "<title>%s</title>",
    "<style>",
    "body{font-family:Arial,Helvetica,sans-serif;margin:18px;line-height:1.35;}",
    ".box{padding:10px 12px;border:1px solid #ddd;border-radius:10px;margin:10px 0;}",
    "img{max-width:100%;height:auto;border:1px solid #eee;border-radius:8px;margin:10px 0;}",
    "code{background:#f7f7f7;padding:2px 5px;border-radius:6px;}",
    "</style></head><body>",
    "<h2>Author profile analysis (offline)</h2>",
    sprintf("<div class='box'><b>Generated:</b> %s<br><b>Query:</b> %s<br><b>Hits:</b> %s</div>",
            esc(format(Sys.time(), "%Y-%m-%d %H:%M:%S")), esc(term), esc(hits))
  )

  if (length(imgs)) {
    lines <- c(lines, "<h3>Figures</h3>")
    for (f in imgs) {
      lines <- c(lines, sprintf("<div class='box'><div><code>%s</code></div><img src='%s'/></div>", esc(f), esc(f)))
    }
  } else {
    lines <- c(lines, "<div class='box'>No figures found in <code>download/</code> yet.</div>")
  }

  if (length(tabs)) {
    lines <- c(lines, "<h3>Tables & files</h3><ul>")
    for (f in tabs) lines <- c(lines, sprintf("<li><a href='%s'>%s</a></li>", esc(f), esc(f)))
    lines <- c(lines, "</ul>")
  }

  lines <- c(lines,
             "<hr><div class='box'>Tip: keep this folder structure. Open <code>index.html</code> locally in a browser.</div>",
             "</body></html>")

  writeLines(lines, file.path(download_dir, "index.html"), useBytes = TRUE)
  invisible(file.path(download_dir, "index.html"))
}

.zip_download_dir <- function(zipfile, download_dir="download"){
  if (!dir.exists(download_dir)) dir.create(download_dir, recursive = TRUE, showWarnings = FALSE)
  fpaths <- list.files(download_dir, recursive = TRUE, full.names = TRUE)
  if (!length(fpaths)) {
    readme <- file.path(download_dir, "README.txt")
    writeLines("No files in download/ yet. Run analysis first to generate figures/tables.", readme, useBytes = TRUE)
    fpaths <- readme
  }
  if (requireNamespace("zip", quietly = TRUE)) {
    # root=download_dir keeps paths relative to download/
    zip::zipr(zipfile = zipfile, files = fpaths, root = download_dir)
  } else {
    owd <- getwd(); on.exit(setwd(owd), add=TRUE)
    setwd(download_dir)
    utils::zip(
    zipfile = zipfile,
    files   = list.files(".", recursive = TRUE, all.files = TRUE, no.. = TRUE)
  )
  }
  invisible(zipfile)
}

# ------------------------------------------------------------
# Use the external module runner (DO NOT modify the module)
# - The module defines: run_flca_ms_sil_runner(), as_LFW(), etc.
# - We alias run_flca_ms_sil() to the module runner, so the rest of
#   the app can call run_flca_ms_sil(...) exactly as you specified.
# ------------------------------------------------------------
try({
  if (file.exists("flca_ms_sil_module.R")) {
    source("flca_ms_sil_module.R", local = FALSE, encoding = "UTF-8")
source("helper_ss_patch.R")
  }
}, silent = TRUE)
if (exists("run_flca_ms_sil_runner")) run_flca_ms_sil <- run_flca_ms_sil_runner



# ============================================================
# Unified runner: FLCA + MajorSampling + Silhouette(SS) + a*
# (Mirrors the proven NON-Shiny script; does NOT modify the module.)
# Returns: list(nodes=nodes20, data=edges20, nodes_full=nodes_full, edges_full=edges0)
# ============================================================
run_flca_ms_sil_internal <- function(nodes, edges0, cfg, verbose = FALSE){
  log <- function(...) if (isTRUE(verbose)) message(...)
  stopifnot(is.data.frame(nodes), is.data.frame(edges0))
  cfg <- modifyList(list(top_clusters=5, base_per_cluster=4, target_n=20,
                         intra_delta=2, inter_delta=5, eps=1e-9), cfg %||% list())

  # normalize + enforce schema
  nodes0 <- normalize_nodes(nodes)
  nodes0 <- nodes0[, intersect(c("name","value"), names(nodes0)), drop=FALSE]
  nodes0$name  <- as.character(nodes0$name)
  nodes0$value <- suppressWarnings(as.numeric(nodes0$value))

  edges0 <- normalize_edges(edges0)
  # allow Source/Target -> Leader/Follower
  if (all(c("Source","Target") %in% names(edges0)) && !all(c("Leader","Follower") %in% names(edges0))) {
    names(edges0)[match(c("Source","Target"), names(edges0))] <- c("Leader","Follower")
  }
  # allow follower -> Follower
  if ("follower" %in% names(edges0) && !("Follower" %in% names(edges0))) {
    names(edges0)[names(edges0)=="follower"] <- "Follower"
  }
  stopifnot(all(c("Leader","Follower","WCD") %in% names(edges0)))
  edges0 <- edges0[, c("Leader","Follower","WCD"), drop=FALSE]
  edges0$Leader   <- as.character(edges0$Leader)
  edges0$Follower <- as.character(edges0$Follower)
  edges0$WCD      <- suppressWarnings(as.numeric(edges0$WCD))
  edges0 <- edges0[is.finite(edges0$WCD) & edges0$WCD > 0, , drop=FALSE]

  log("[run_flca_ms_sil] nodes=", nrow(nodes0), " edges=", nrow(edges0))

  if (nrow(nodes0) == 0 || nrow(edges0) == 0) {
    # keep empty but consistent columns
    out_nodes <- nodes0
    out_nodes$carac <- if (!"carac" %in% names(out_nodes)) 1 else out_nodes$carac
    out_nodes$value2 <- if (!"value2" %in% names(out_nodes)) out_nodes$value else out_nodes$value2
    out_nodes$ssi <- 0; out_nodes$a_i <- 0; out_nodes$b_i <- 0; out_nodes$a_star1 <- 0
    out_nodes$SSi <- out_nodes$ssi; out_nodes$a_star <- out_nodes$a_star1
    return(list(nodes=out_nodes, data=edges0, nodes_full=out_nodes, edges_full=edges0))
  }

  # symmetric edges for clustering (FLCA wants directed)
  edges_sym <- rbind(
    edges0,
    data.frame(Leader=edges0$Follower, Follower=edges0$Leader, WCD=edges0$WCD, stringsAsFactors = FALSE)
  )
  edges_sym$follower <- edges_sym$Follower
  edges_sym <- edges_sym[, c("Leader","follower","WCD"), drop=FALSE]

  fl <- run_flca_ms_sil(list(nodes = nodes0, data = edges_sym))
  nodes_full <- fl$nodes
  # value2 strength from ORIGINAL undirected edges0
  nodes_full <- compute_value2_strength(nodes_full, edges0)

  ms <- major_sample_flca_top20(
    nodes_full, edges0,
    top_clusters = cfg$top_clusters,
    base_per_cluster = cfg$base_per_cluster,
    target_n = cfg$target_n
  )
  nodes20 <- ms$nodes
  edges20 <- ms$data

  if ("follower" %in% names(edges20) && !("Follower" %in% names(edges20))) {
    names(edges20)[names(edges20)=="follower"] <- "Follower"
  }
  # Some pipelines output Source/Target
  if (all(c("Source","Target") %in% names(edges20)) && !all(c("Leader","Follower") %in% names(edges20))) {
    names(edges20)[match(c("Source","Target"), names(edges20))] <- c("Leader","Follower")
  }
  if (nrow(edges20)) {
    edges20 <- edges20[, c("Leader","Follower","WCD"), drop=FALSE]
  } else {
    edges20 <- data.frame(Leader=character(), Follower=character(), WCD=numeric(), stringsAsFactors = FALSE)
  }

  # recompute within-sample strength
  nodes20 <- compute_value2_strength(nodes20, edges20)

  # silhouette + a*
  sil_df <- compute_silhouette_df(nodes20, edges20,
                                 intra_delta = cfg$intra_delta,
                                 inter_delta = cfg$inter_delta,
                                 eps = cfg$eps)
  nodes20 <- merge(nodes20, sil_df, by="name", all.x=TRUE)
  nodes20$ssi[is.na(nodes20$ssi) | !is.finite(nodes20$ssi)] <- 0
  nodes20$a_i[is.na(nodes20$a_i) | !is.finite(nodes20$a_i)] <- 0
  nodes20$b_i[is.na(nodes20$b_i) | !is.finite(nodes20$b_i)] <- 0
  nodes20$a_star1[is.na(nodes20$a_star1) | !is.finite(nodes20$a_star1)] <- 0
  # back-compat aliases
  if (!("SSi" %in% names(nodes20))) nodes20$SSi <- nodes20$ssi
  if (!("a_star" %in% names(nodes20))) nodes20$a_star <- nodes20$a_star1

  list(nodes=nodes20, data=edges20, nodes_full=nodes_full, edges_full=edges0)
}

# ---- HS / OneLink silhouette module (prefer user's silhouette.R; no igraph required) ----
try({
  if (FALSE && file.exists("silhouette.R")) {
    source("silhouette.R", local = TRUE, encoding = "UTF-8")
    message("[BOOT] sourced: silhouette.R (HS/OneLink)")
  } else if (FALSE && file.exists(file.path(getwd(), "silhouette.R"))) {
    source(file.path(getwd(), "silhouette.R"), local = TRUE, encoding = "UTF-8")
    message("[BOOT] sourced: silhouette.R (HS/OneLink)")
  } else {
    message("[BOOT] silhouette.R not found; will use fallback silhouette (may require igraph).")
  }
}, silent = TRUE)



# ==== REPORT DIAGNOSTICS HELPERS (SAFE; pure functions) ====
# NOTE: Do NOT touch `output$` at top-level. Bind inside server().

bind_yearcount_top20_dt <- function(output, input, session, rv) {
  output$tbl_yearcount_top20 <- DT::renderDT({
    df <- if (!is.null(rv$pubmeta)) rv$pubmeta else NULL
    shiny::validate(shiny::need(is.data.frame(df) && "Year" %in% names(df), "Year not available."))

    domain_col <- .guess_domain_col(df)
    shiny::validate(shiny::need(!is.null(domain_col), "找不到 metadata 欄位。"))

    # Top20 terms: prefer rv$nodes20$name, else take first 20 terms in metadata
    if (!is.null(rv$nodes20) && is.data.frame(rv$nodes20) && "name" %in% names(rv$nodes20)) {
      top20 <- as.character(rv$nodes20$name)
    } else {
      terms_all <- unique(unlist(lapply(df[[domain_col]], .split_terms)))
      top20 <- utils::head(terms_all, 20)
    }
    top20 <- top20[!is.na(top20) & nzchar(top20)]
    shiny::validate(shiny::need(length(top20) > 0, "Top20 名稱不存在。"))

    yc <- compute_top20_year_counts(df, top20, domain_col)
    shiny::validate(shiny::need(!is.null(yc), "未滿 10 年，不產生年度次數。"))

    st <- summarize_prepost_ttest(yc)
    if (is.null(st) || !nrow(st)) {
      st <- data.frame(item=unique(yc$item), mean_pre=NA_real_, mean_post=NA_real_, pval=NA_real_, inc=NA, stringsAsFactors=FALSE)
    }

    # Wide counts (base R reshape; avoids extra deps)
    wide <- reshape(yc, idvar="item", timevar="Year", direction="wide")
    names(wide) <- sub("^Count\\.", "", names(wide))
    wide <- merge(st, wide, by="item", all.x=TRUE)
    if ("pval" %in% names(wide)) wide <- wide[order(is.na(wide$pval), wide$pval, -wide$mean_post), , drop=FALSE]

    DT::datatable(wide, options=list(pageLength=10, scrollX=TRUE))
  })
}

# ==== END REPORT DIAGNOSTICS HELPERS ====
.safe_nrow <- function(x){
  if (is.null(x)) return(0L)
  if (!is.data.frame(x)) return(0L)
  nrow(x)
}
as_scalar_or_na <- function(x, na=NA){ if (is.null(x) || length(x)==0 || all(is.na(x))) return(na); x[[1]] }
.compute_AAC <- function(x){
  # Absolute Advantage Coefficient (AAC)
  # gamma = (v1/v2)/(v2/v3) = v1*v3 / v2^2
  # AAC = gamma/(1+gamma)
  if (is.null(x) || length(x)==0 || all(is.na(x))) return(NA_real_)
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x) & !is.na(x)]
  if (length(x) < 3) return(NA_real_)
  x <- sort(abs(x), decreasing = TRUE)
  v1 <- x[1]; v2 <- x[2]; v3 <- x[3]
  eps <- 1e-9
  gamma <- (v1 * (v3 + eps)) / ((v2 + eps)^2)
  if (!is.finite(gamma) || gamma < 0) return(NA_real_)
  gamma / (1 + gamma)
}

.get_obj_any <- function(candidates, env = parent.frame()){
  for (cand in candidates){
    if (is.character(cand) && grepl("^[A-Za-z.][A-Za-z0-9._]*$", cand)){
      obj <- tryCatch(get0(cand, envir=env, inherits=TRUE), error=function(e) NULL)
      if (!is.null(obj)) return(obj)
    }
    obj <- tryCatch(eval(parse(text=cand), envir=env), error=function(e) NULL)
    if (!is.null(obj)) return(obj)
  }
  NULL
}
.get_col <- function(df, col){
  if (is.null(df) || !is.data.frame(df)) return(NULL)
  if (!col %in% names(df)) return(NULL)
  df[[col]]
}
.aac_if_edges <- function(nodes, edges, col){
  if (.safe_nrow(nodes) == 0L) return(NA_real_)
  if (.safe_nrow(edges) == 0L) return(NA_real_)
  .compute_AAC(.get_col(nodes, col))
}
.get_first_nonnull_col <- function(df, cols){
  for (cc in cols){
    v <- .get_col(df, cc)
    if (!is.null(v)) return(v)
  }
  NULL
}
.aac_a_star_if_edges <- function(nodes, edges){
  if (.safe_nrow(nodes) == 0L) return(NA_real_)
  if (.safe_nrow(edges) == 0L) return(NA_real_)
  .compute_AAC(.get_first_nonnull_col(nodes, c('a_star1','a_star')))
}

# ---- Silhouette fallback (no extra deps) ----
# Returns data.frame(name, sil_width, a_i, b_i, a_star1)
compute_silhouette_df_fallback <- function(nodes, data,
                                  intra_delta = 2,
                                  inter_delta = 5,
                                  eps = 1e-9) {

  stopifnot(is.data.frame(nodes), nrow(nodes) > 0)

  # --- choose a stable key column for names (prefer orig_name if present) ---
  key_col <- if ("orig_name" %in% names(nodes)) "orig_name" else "name"
  if (!key_col %in% names(nodes)) stop("nodes must have a name column ('name' or 'orig_name').", call. = FALSE)
  if (!"carac" %in% names(nodes)) stop("nodes must have 'carac' (cluster) column.", call. = FALSE)

  nd <- tibble::as_tibble(nodes) %>%
    transmute(
      name      = trimws(as.character(.data[[key_col]])),
      carac_raw = suppressWarnings(as.integer(as.character(.data[["carac"]])))
    ) %>%
    filter(nzchar(name)) %>%
    distinct(name, .keep_all = TRUE)

  if (anyNA(nd$carac_raw)) stop("nodes$carac contains NA after coercion.", call. = FALSE)

  # If <2 clusters, silhouette is not meaningful; return zeros (and NA a/b)
  if (dplyr::n_distinct(nd$carac_raw) < 2) {
    return(nd %>%
      transmute(
        name = name,
        ssi  = 0,
        a_i  = NA_real_,
        b_i  = NA_real_,
        a_star1 = 0
      ))
  }

  stopifnot(is.data.frame(data), nrow(data) > 0)
  if (!all(c("Leader","follower","WCD") %in% names(data))) {
    stop("data must have columns: Leader, follower, WCD", call. = FALSE)
  }

  edges <- tibble::as_tibble(data) %>%
    transmute(
      Leader   = trimws(as.character(Leader)),
      follower = trimws(as.character(follower)),
      WCD      = suppressWarnings(as.numeric(as.character(WCD)))
    ) %>%
    filter(nzchar(Leader), nzchar(follower), is.finite(WCD), WCD > 0)

  nms <- nd$name
  edges <- edges %>% filter(Leader %in% nms, follower %in% nms)

  # --- undirected: merge multi-edges by SUM weight on unordered pair ---
  edges_und <- edges %>%
    transmute(a = pmin(Leader, follower),
              b = pmax(Leader, follower),
              WCD = WCD) %>%
    group_by(a, b) %>%
    summarise(WCD = sum(WCD, na.rm = TRUE), .groups = "drop")

  # --- full W matrix (0 = missing) ---
  full_pairs <- expand.grid(a = nms, b = nms, stringsAsFactors = FALSE)
  W_df <- full_pairs %>%
    left_join(edges_und, by = c("a","b")) %>%
    mutate(WCD = tidyr::replace_na(WCD, 0))

  W <- W_df %>%
    pivot_wider(names_from = b, values_from = WCD) %>%
    column_to_rownames("a") %>%
    as.matrix()

  # symmetry + zero diagonal
  W <- pmax(W, t(W))
  diag(W) <- 0

  # --- cost / distance with penalties on EMPTY edges ---
  cost <- matrix(NA_real_, nrow(W), ncol(W), dimnames = dimnames(W))
  has_edge <- (W > 0)
  cost[has_edge] <- 1 / (W[has_edge] + eps)
  max_cost <- if (any(has_edge)) max(cost[has_edge], na.rm = TRUE) else 1

  cl <- setNames(nd$carac_raw, nd$name)[nms]
  same_cluster <- outer(cl, cl, FUN = "==")

  D <- matrix(NA_real_, nrow(W), ncol(W), dimnames = dimnames(W))
  D[has_edge] <- cost[has_edge]
  D[!has_edge &  same_cluster] <- max_cost * (1 + intra_delta)  # EMPTY intra
  D[!has_edge & !same_cluster] <- max_cost * (1 + inter_delta)  # EMPTY inter
  diag(D) <- 0

  D_sym <- pmin(D, t(D))
  diag(D_sym) <- 0
  d_obj <- as.dist(D_sym)

  # silhouette requires cluster ids 1..K
  cl_int <- as.integer(factor(cl, levels = sort(unique(cl))))
  sil <- cluster::silhouette(cl_int, d_obj)

  ssi <- as.numeric(sil[, "sil_width"])
  names(ssi) <- rownames(D_sym)

  # a(i), b(i)
  diss <- as.matrix(D_sym)
  a_i <- numeric(length(nms))
  b_i <- numeric(length(nms))

  for (i in seq_along(nms)) {
    same_idx <- which(cl_int == cl_int[i])
    same_idx <- setdiff(same_idx, i)
    a_i[i] <- if (length(same_idx) > 0) mean(diss[i, same_idx]) else NA_real_

    other <- setdiff(unique(cl_int), cl_int[i])
    if (length(other) > 0) {
      d_to <- vapply(other, function(k) mean(diss[i, which(cl_int == k)]), numeric(1))
      b_i[i] <- min(d_to)
    } else {
      b_i[i] <- NA_real_
    }
  }

  # a* (your requested transform)
  a_star1 <- 1 / (1 + a_i + eps)

  tibble(
    name    = nms,
    ssi     = as.numeric(ssi[nms]),
    a_i     = a_i,
    b_i     = b_i,
    a_star1 = as.numeric(a_star1)
  )
}

# Pick silhouette function if renderSSplot.R didn't provide one
if (!exists('compute_silhouette_df', mode='function')) {
  compute_silhouette_df <- compute_silhouette_df_fallback
}

# ---- AAC updater (compute whenever a domain run completes) ----
.update_domain_aac <- function(rv, prefix, nodes, edges){
  # Always compute AAC for value/value2 when nodes exist
  if (is.null(nodes) || !is.data.frame(nodes) || nrow(nodes) == 0) {
    rv[[paste0('AAC_value_',prefix)]]  <- 0
    rv[[paste0('AAC_value2_',prefix)]] <- 0
    rv[[paste0('AAC_ss_',prefix)]]     <- 0
    rv[[paste0('AAC_a_star_',prefix)]] <- 0
    return(invisible(NULL))
  }

  if ('value' %in% names(nodes))  rv[[paste0('AAC_value_',prefix)]]  <- .compute_AAC(nodes$value)
  if ('value2' %in% names(nodes)) rv[[paste0('AAC_value2_',prefix)]] <- .compute_AAC(nodes$value2)

  # Only compute ss/a* AAC if these columns exist AND have some non-NA values
  if ('ssi' %in% names(nodes) && any(is.finite(nodes$ssi))) {
    rv[[paste0('AAC_ss_',prefix)]] <- .compute_AAC(nodes$ssi)
  } else {
    rv[[paste0('AAC_ss_',prefix)]] <- 0
  }
  if ('a_star1' %in% names(nodes) && any(!is.na(nodes$a_star1))) {
    rv[[paste0('AAC_a_star_',prefix)]] <- .compute_AAC(nodes$a_star1)
  } else if ('a_star' %in% names(nodes) && any(!is.na(nodes$a_star))) {
    rv[[paste0('AAC_a_star_',prefix)]] <- .compute_AAC(nodes$a_star)
  } else {
    rv[[paste0('AAC_a_star_',prefix)]] <- 0
  }
  invisible(NULL)
}


# app.R — Interactive bibliometrics dashboard (shinyapps.io ready)
# PubMed (term) → Author (First–Last) FLCA+MajorSampling → interactive network
# plus Country / Institute / MeSH dashboards ALSO via FLCA+MajorSampling (one-link structure)
# and a one-click self-contained HTML report.

suppressPackageStartupMessages({
  library(shiny)
  library(rentrez)
  library(DT)
  library(visNetwork)
  library(igraph)
  library(ggplot2)
  library(maps)
  library(htmlwidgets)
  if (requireNamespace("reticulate", quietly = TRUE)) library(reticulate)
  if (requireNamespace("uwot", quietly = TRUE)) library(uwot)
  if (requireNamespace("ggraph", quietly = TRUE)) library(ggraph)

# =========================================================
# Shinyapps-safe writable paths
# =========================================================
is_cloud <- function() {
  nzchar(Sys.getenv("RSCONNECT_SERVER")) ||
    nzchar(Sys.getenv("RSCONNECT_ACCOUNT")) ||
    nzchar(Sys.getenv("SHINY_PORT")) ||
    nzchar(Sys.getenv("SHINY_SERVER_VERSION"))
}

APP_DIR <- tryCatch(normalizePath(getwd(), winslash="/", mustWork=FALSE), error=function(e) getwd())

PERM_PUBMED_TXT <- if (is_cloud()) file.path(tempdir(), "uploaded_pubmed.txt") else PERM_PUBMED_TXT
PERM_PUBMED_XML <- if (is_cloud()) file.path(tempdir(), "pubmed.xml") else PERM_PUBMED_XML
PERM_PUBMED_MEDLINE <- if (is_cloud()) file.path(tempdir(), "pubmed_medline.txt") else PERM_PUBMED_MEDLINE

RUNTIME_DIR <- if (is_cloud()) file.path(tempdir(), "runtime") else file.path(APP_DIR, "runtime")
dir.create(RUNTIME_DIR, showWarnings = FALSE, recursive = TRUE)

source("ipmodule.R")

})


# ----------------------------
# AAC inline helper (used for value_strength)
# AAC = r/(1+r), r=(v1/v2)/(v2/v3) with v sorted desc (top3)
# ----------------------------
.AAC_INLINE <- function(v){
  v <- suppressWarnings(as.numeric(v))
  v <- v[is.finite(v) & !is.na(v)]
  v <- sort(v, decreasing=TRUE)
  if (length(v) < 3) return(NA_real_)
  v1 <- v[1]; v2 <- v[2]; v3 <- v[3]
  if (v2 <= 0 || v3 <= 0) return(NA_real_)
  r <- (v1/v2)/(v2/v3)
  r/(1+r)
}




# ---- BCa/EIS helpers: conditional follower bootstrap for AAC/EIS elbow identification ----
# IMPORTANT METHOD CHANGE:
# For each candidate rank i, x[i] is fixed as the tested value.
# Only the lower-ranked followers x[(i+1):n] are bootstrapped.
# This prevents the candidate/top value from being duplicated or removed during bootstrap.
# BCa/EIS is computed only when there are at least 3 followers after the candidate.
aac3 <- function(x1, x2, x3) {
  x <- suppressWarnings(as.numeric(c(x1, x2, x3)))
  if (any(!is.finite(x)) || any(x <= 0)) return(NA_real_)
  (x[1] / x[2]) / (x[2] / x[3])
}

.scan_one_conditional_bca <- function(self_value, followers, R = 500, cutoff = 2.0, conf = 0.95) {
  followers <- suppressWarnings(as.numeric(followers))
  followers <- sort(followers[is.finite(followers) & followers > 0], decreasing = TRUE)

  if (!is.finite(self_value) || self_value <= 0) {
    return(list(
      AAC = NA_real_, BCa_lower = NA_real_, BCa_upper = NA_real_,
      decision = "No BCa/EIS", EIS_point = NA_real_,
      reason = "Candidate value is not a positive finite number."
    ))
  }

  if (length(followers) < 3) {
    return(list(
      AAC = NA_real_, BCa_lower = NA_real_, BCa_upper = NA_real_,
      decision = "No BCa/EIS", EIS_point = NA_real_,
      reason = paste0("Fewer than 3 followers after this rank (n_followers=", length(followers), ").")
    ))
  }

  obs_r <- aac3(self_value, followers[1], followers[2])

  stat_fun <- function(data, idx) {
    d <- sort(data[idx], decreasing = TRUE)
    d <- d[is.finite(d) & d > 0]
    if (length(d) < 2) return(NA_real_)
    aac3(self_value, d[1], d[2])
  }

  bt <- boot::boot(followers, stat_fun, R = R)
  boot_vals <- suppressWarnings(as.numeric(bt$t[, 1]))
  boot_vals <- boot_vals[is.finite(boot_vals)]

  ci <- tryCatch(
    boot::boot.ci(bt, type = "bca", conf = conf),
    error = function(e) NULL
  )

  if (!is.null(ci) && !is.null(ci$bca)) {
    lower <- suppressWarnings(as.numeric(ci$bca[4]))
    upper <- suppressWarnings(as.numeric(ci$bca[5]))
    reason <- "BCa interval computed by fixing the candidate and bootstrapping only lower-ranked followers."
  } else if (length(boot_vals) >= 10) {
    lower <- as.numeric(stats::quantile(boot_vals, probs = (1 - conf) / 2, na.rm = TRUE, names = FALSE))
    upper <- as.numeric(stats::quantile(boot_vals, probs = 1 - (1 - conf) / 2, na.rm = TRUE, names = FALSE))
    reason <- "BCa unavailable; percentile bootstrap fallback used by fixing the candidate and resampling followers."
  } else {
    lower <- NA_real_
    upper <- NA_real_
    reason <- "Bootstrap interval unavailable because too few finite bootstrap AAC values were produced."
  }

  decision <- ifelse(
    is.finite(lower) && lower > cutoff,
    "Significant elbow",
    ifelse(is.finite(upper) && upper < cutoff, "Not elbow", "Uncertain")
  )

  list(
    AAC = obs_r,
    BCa_lower = lower,
    BCa_upper = upper,
    decision = decision,
    EIS_point = lower / cutoff,
    reason = reason
  )
}

scan_aac_bca <- function(x, R = 500, cutoff = 2.0, conf = 0.95, progress = NULL) {
  if (!requireNamespace("boot", quietly = TRUE)) {
    stop("Package 'boot' is required for BCa/EIS. Please install.packages('boot').", call. = FALSE)
  }

  x <- suppressWarnings(as.numeric(x))
  x <- sort(x[is.finite(x) & x > 0], decreasing = TRUE)
  n <- length(x)
  if (n < 5) stop("BCa/EIS requires at least n >= 5 positive values; n >= 20 is recommended for stable intervals.", call. = FALSE)

  # Scan every ranked candidate. BCa/EIS is reported only where >=3 followers are available.
  total_points <- n
  out <- lapply(seq_len(n), function(i) {
    if (is.function(progress)) progress(i, total_points)

    self_value <- x[i]
    followers <- if (i < n) x[(i + 1):n] else numeric(0)
    f1 <- if (length(followers) >= 1) followers[1] else NA_real_
    f2 <- if (length(followers) >= 2) followers[2] else NA_real_
    f3 <- if (length(followers) >= 3) followers[3] else NA_real_

    z <- .scan_one_conditional_bca(
      self_value = self_value,
      followers = followers,
      R = R,
      cutoff = cutoff,
      conf = conf
    )

    data.frame(
      point_index = i,
      point_value = self_value,
      x_self = self_value,
      follower1 = f1,
      follower2 = f2,
      follower3 = f3,
      n_followers = length(followers),
      AAC = z$AAC,
      BCa_lower = z$BCa_lower,
      BCa_upper = z$BCa_upper,
      cutoff = cutoff,
      decision = z$decision,
      EIS_point = z$EIS_point,
      BCa_reason = z$reason,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, out)
}

# ---- BCa/EIS manual-value parser ----
.parse_bca_values <- function(txt) {
  if (is.null(txt) || !nzchar(txt)) return(numeric(0))
  txt <- enc2utf8(as.character(txt))
  txt <- gsub("[，；、]+", " ", txt)
  txt <- gsub("[,;\\t\\r\\n]+", " ", txt)
  z <- unlist(strsplit(txt, "[[:space:]]+", perl = TRUE))
  z <- trimws(z)
  z <- z[nzchar(z)]
  x <- suppressWarnings(as.numeric(z))
  x <- x[is.finite(x) & x > 0]
  sort(x, decreasing = TRUE)
}


# ---- BCa/EIS summary helper ----
# EIS = max(BCa_lower / cutoff). This converts the BCa evidence for an elbow
# into an interpretable evidence score: EIS > 1 means BCa_lower exceeds cutoff,
# therefore the elbow is statistically supported.
summarize_bca_eis <- function(df, cutoff = 2.0) {
  if (is.null(df) || !is.data.frame(df) || !nrow(df)) {
    return(data.frame(
      selected_point = NA_integer_, selected_value = NA_real_, AAC = NA_real_,
      BCa_lower = NA_real_, BCa_upper = NA_real_, EIS = NA_real_,
      classification = "No result", rule = "No BCa/EIS result available.",
      stringsAsFactors = FALSE
    ))
  }

  df$point_index <- suppressWarnings(as.numeric(df$point_index))
  df$point_value <- suppressWarnings(as.numeric(df$point_value))
  df$AAC <- suppressWarnings(as.numeric(df$AAC))
  df$BCa_lower <- suppressWarnings(as.numeric(df$BCa_lower))
  df$BCa_upper <- suppressWarnings(as.numeric(df$BCa_upper))
  df$cutoff <- cutoff
  df$EIS <- df$BCa_lower / cutoff

  sig <- which(is.finite(df$BCa_lower) & df$BCa_lower > cutoff)
  if (length(sig) > 0) {
    k <- sig[1]
    class <- "Significant elbow"
    rule <- paste0("Selected the first point where BCa lower bound > cutoff (",
                   round(df$BCa_lower[k], 3), " > ", cutoff, ").")
  } else {
    cand <- which(is.finite(df$EIS))
    if (!length(cand)) cand <- which(is.finite(df$AAC))
    if (length(cand)) {
      k <- cand[which.max(df$EIS[cand])]
      class <- "No significant elbow"
      rule <- "No point had BCa lower bound > cutoff. The displayed point is the strongest non-significant candidate based on max EIS = BCa_lower/cutoff."
    } else {
      k <- 1
      class <- "No significant elbow"
      rule <- "No finite BCa/EIS candidate was available."
    }
  }

  data.frame(
    selected_point = df$point_index[k],
    selected_value = df$point_value[k],
    AAC = df$AAC[k],
    BCa_lower = df$BCa_lower[k],
    BCa_upper = df$BCa_upper[k],
    EIS = df$EIS[k],
    classification = class,
    rule = rule,
    stringsAsFactors = FALSE
  )
}



# ---- Top-1 conditional BCa/EIS helper ----
# Purpose: test whether the top-ranked value is a stable dominance elbow.
# The top-1 value is fixed, and only the remaining lower-ranked values are bootstrapped.
# This avoids the unstable 3-point bootstrap case where the top value may be duplicated or removed.
.aac_coef_from_or <- function(r) {
  r <- suppressWarnings(as.numeric(r))
  r / (1 + r)
}

scan_top1_bca_eis <- function(x, R = 500, cutoff_coef = 0.70, conf = 0.95) {
  if (!requireNamespace("boot", quietly = TRUE)) {
    stop("Package 'boot' is required for Top-1 BCa/EIS. Please install.packages('boot').", call. = FALSE)
  }

  x <- suppressWarnings(as.numeric(x))
  x <- sort(x[is.finite(x) & x > 0], decreasing = TRUE)
  n <- length(x)
  if (n < 5) stop("Top-1 BCa/EIS requires at least 5 positive values.", call. = FALSE)

  top1 <- x[1]
  pool <- x[-1]
  pool <- pool[is.finite(pool) & pool > 0]
  if (length(pool) < 3) stop("Top-1 BCa/EIS requires at least 3 followers below top 1.", call. = FALSE)

  obs_or <- aac3(top1, pool[1], pool[2])
  obs_coef <- .aac_coef_from_or(obs_or)

  stat_fun <- function(data, idx) {
    d <- sort(data[idx], decreasing = TRUE)
    if (length(d) < 2 || any(!is.finite(d[1:2])) || any(d[1:2] <= 0)) return(NA_real_)
    .aac_coef_from_or(aac3(top1, d[1], d[2]))
  }

  bt <- boot::boot(pool, stat_fun, R = R)
  boot_vals <- suppressWarnings(as.numeric(bt$t[, 1]))
  boot_vals <- boot_vals[is.finite(boot_vals)]

  ci <- tryCatch(
    boot::boot.ci(bt, type = "bca", conf = conf),
    error = function(e) NULL
  )

  if (!is.null(ci) && !is.null(ci$bca)) {
    lower <- suppressWarnings(as.numeric(ci$bca[4]))
    upper <- suppressWarnings(as.numeric(ci$bca[5]))
    ci_method <- "BCa"
  } else if (length(boot_vals) >= 10) {
    lower <- as.numeric(stats::quantile(boot_vals, probs = (1 - conf) / 2, na.rm = TRUE, names = FALSE))
    upper <- as.numeric(stats::quantile(boot_vals, probs = 1 - (1 - conf) / 2, na.rm = TRUE, names = FALSE))
    ci_method <- "Percentile fallback"
  } else {
    lower <- NA_real_; upper <- NA_real_; ci_method <- "CI unavailable"
  }

  eis <- lower / cutoff_coef
  decision <- ifelse(
    is.finite(lower) && lower >= cutoff_coef,
    "BCa-confirmed top-1 elbow",
    ifelse(is.finite(obs_coef) && obs_coef >= cutoff_coef,
           "Strong observed top-1 elbow; not BCa-confirmed",
           "Not top-1 elbow")
  )

  data.frame(
    top1 = top1,
    x2 = pool[1],
    x3 = pool[2],
    AAC_OR = obs_or,
    AAC_coef = obs_coef,
    CI_lower = lower,
    CI_upper = upper,
    cutoff_coef = cutoff_coef,
    EIS = eis,
    decision = decision,
    ci_method = ci_method,
    R = R,
    conf = conf,
    stringsAsFactors = FALSE
  )
}



# ---- General Top-k conditional BCa/EIS helper (k = 1, 2, or 3) ----
# The candidate at rank k is fixed. Only lower-ranked followers after k are resampled.
# AAC_OR = (candidate/follower1)/(follower1/follower2). EIS_OR = CI_lower_OR / AAC_OR cutoff.
scan_topk_bca_eis <- function(x, k = 1, R = 500, cutoff_or = 2.0, conf = 0.95) {
  if (!requireNamespace("boot", quietly = TRUE)) {
    stop("Package 'boot' is required for Top-k BCa/EIS. Please install.packages('boot').", call. = FALSE)
  }
  x <- suppressWarnings(as.numeric(x))
  x <- sort(x[is.finite(x) & x > 0], decreasing = TRUE)
  n <- length(x)
  k <- suppressWarnings(as.integer(k[1]))
  if (!is.finite(k) || !(k %in% c(1L, 2L, 3L))) k <- 1L
  if (n < 5) stop("Top-k BCa/EIS requires at least 5 positive values.", call. = FALSE)
  if (n - k < 3) stop(paste0("Top-", k, " BCa/EIS requires at least 3 followers below the candidate."), call. = FALSE)

  candidate <- x[k]
  followers <- x[(k + 1):n]
  z <- .scan_one_conditional_bca(
    self_value = candidate,
    followers = followers,
    R = R,
    cutoff = cutoff_or,
    conf = conf
  )
  coef <- .aac_coef_from_or(z$AAC)
  lower_coef <- .aac_coef_from_or(z$BCa_lower)
  upper_coef <- .aac_coef_from_or(z$BCa_upper)

  data.frame(
    selected_rank = k,
    candidate = candidate,
    follower1 = followers[1],
    follower2 = followers[2],
    follower3 = followers[3],
    n_followers = length(followers),
    AAC_OR = z$AAC,
    AAC_coef = coef,
    CI_lower_OR = z$BCa_lower,
    CI_upper_OR = z$BCa_upper,
    CI_lower_coef = lower_coef,
    CI_upper_coef = upper_coef,
    AAC_OR_cutoff = cutoff_or,
    EIS_OR = z$BCa_lower / cutoff_or,
    decision = z$decision,
    BCa_reason = z$reason,
    stringsAsFactors = FALSE
  )
}

.bca_caption_text <- function(z, source_label = "Editable values") {
  if (is.null(z) || !is.data.frame(z) || !nrow(z)) return("No Top-k BCa/EIS result available.")
  paste0(
    "Top-", z$selected_rank[1], " conditional BCa/EIS result (source: ", source_label, "). ",
    "The candidate value was fixed at ", round(z$candidate[1], 3),
    ", and only lower-ranked followers were resampled. ",
    "AAC_OR = ", round(z$AAC_OR[1], 3),
    ", BCa CI_OR = [", round(z$CI_lower_OR[1], 3), ", ", round(z$CI_upper_OR[1], 3), "]",
    ", AAC_OR cutoff = ", round(z$AAC_OR_cutoff[1], 3),
    ", and EIS_OR = BCa_lower_OR / cutoff = ", round(z$EIS_OR[1], 3), ". ",
    "The decision was: ", z$decision[1], ". ",
    "Criterion: the elbow is statistically supported when BCa lower bound of AAC_OR exceeds the AAC_OR cutoff, equivalently EIS_OR > 1."
  )
}

.plot_topk_dominance <- function(x, z, main_title = "Top-k Dominance vs Rest: Conditional BCa/EIS") {
  x <- suppressWarnings(as.numeric(x))
  x <- sort(x[is.finite(x) & x > 0], decreasing = TRUE)
  if (length(x) < 5) stop("Top-k dominance plot requires at least 5 positive values.", call. = FALSE)
  if (is.null(z) || !is.data.frame(z) || !nrow(z)) stop("Top-k result is unavailable.", call. = FALSE)
  k <- suppressWarnings(as.integer(z$selected_rank[1]))
  if (!is.finite(k) || k < 1 || k > length(x)) k <- 1L

  ranks <- seq_along(x)
  y_max <- max(x, na.rm = TRUE)
  y_min <- min(x[x > 0], na.rm = TRUE)
  ylim <- c(max(0.1, y_min * 0.75), y_max * 1.35)

  plot(
    ranks, x,
    type = "b", pch = 19, lwd = 2, log = "y",
    xlab = "Rank position",
    ylab = "Entity count / score (log scale)",
    main = main_title,
    ylim = ylim
  )

  if (length(x) >= 4) {
    rest_r <- ranks[(k + 1):length(x)]
    rest_x <- x[(k + 1):length(x)]
    if (length(rest_x) >= 4) {
      sm <- stats::lowess(rest_r, rest_x, f = 2/3, iter = 1)
      lines(sm$x, sm$y, lwd = 3, lty = 1)
    }
  }

  points(k, x[k], pch = 19, cex = 2.4, col = "red")
  text(k, x[k], labels = paste0("Top ", k, "\n", round(x[k], 1)), pos = 4, col = "red", cex = 0.95)

  comp_idx <- (k + 1):min(length(x), k + 3)
  if (length(comp_idx)) {
    points(comp_idx, x[comp_idx], pch = 19, cex = 1.7, col = "orange")
    text(comp_idx, x[comp_idx], labels = paste0("f", seq_along(comp_idx), "\n", round(x[comp_idx], 1)), pos = 4, col = "orange", cex = 0.8)
    segments(k, x[k], comp_idx[1], x[comp_idx[1]], lwd = 2, col = "red")
    if (length(comp_idx) >= 2) segments(comp_idx[1], x[comp_idx[1]], comp_idx[2], x[comp_idx[2]], lwd = 2, col = "orange")
  }

  info <- paste0(
    "Top-", k, " candidate\n",
    "AAC_OR = ", round(z$AAC_OR[1], 2), "\n",
    "BCa CI_OR = [", round(z$CI_lower_OR[1], 2), ", ", round(z$CI_upper_OR[1], 2), "]\n",
    "EIS_OR = ", round(z$EIS_OR[1], 2), "\n",
    z$decision[1]
  )
  legend("topright", legend = info, bty = "n", cex = 0.9)
  legend("bottomleft",
         legend = c("Ranked values (auto-sorted)", "LOWESS of followers", "Fixed candidate", "Follower comparators"),
         lty = c(1, 1, NA, NA), pch = c(19, NA, 19, 19),
         lwd = c(2, 3, NA, NA), col = c("black", "black", "red", "orange"), bty = "n")
}

# ---- country name harmonization for world map (adaptive to map_data names) ----
.normalize_country_for_map <- function(x, world_regions){
  x <- trimws(as.character(x))
  x_upper <- toupper(x)
  out <- x

  # Decide canonical labels used by the current map dataset
  us_label <- if ("USA" %in% world_regions) "USA" else if ("United States" %in% world_regions) "United States" else "USA"
  uk_label <- if ("UK" %in% world_regions) "UK" else if ("United Kingdom" %in% world_regions) "United Kingdom" else "UK"

  out[x_upper %in% c("US","USA","U.S.","U.S","UNITED STATES","UNITED STATES OF AMERICA")] <- us_label
  out[x_upper %in% c("UK","U.K.","U.K","UNITED KINGDOM","UNITED KINGDOM OF GREAT BRITAIN AND NORTHERN IRELAND",
                     "GREAT BRITAIN","ENGLAND","SCOTLAND","WALES","NORTHERN IRELAND")] <- uk_label
  out
}
# ---------------------------------------------------------------------------

# =========================
# MAP (DESCRIPTIVE, PRE-FLCA) - GLOBAL SETUP
# =========================
.world_map_data <- local({
  w <- map_data("world")
  subset(w, region != "Antarctica" & region != "French Southern and Antarctic Lands")
})


# ✅ world map region names used for harmonization
world_regions <- sort(unique(.world_map_data$region))
.country_fix_map <- c(
  "USA" = "United States",
  "UK"  = "United Kingdom",
  "South Korea" = "Korea"
)
source("pubmed_parse_biblio.R")
tryCatch({ if (file.exists("kano(63).R")) source("kano(63).R", local = FALSE, encoding = "UTF-8") else if (file.exists("kano.R")) source("kano.R", local = FALSE, encoding = "UTF-8") else if (file.exists("kano(62).R")) source("kano(62).R", local = FALSE, encoding = "UTF-8"); cat("[BOOT] sourced: kano.R
") }, error=function(e){ cat("[BOOT][ERR] kano.R:", e$message, "
") })
try({ if (file.exists("sankey(59).R")) source("sankey(59).R", local = FALSE, encoding = "UTF-8") else source("sankey.R", local = FALSE, encoding = "UTF-8"); cat("[BOOT] sourced: sankey.R\n") }, silent=TRUE)

# ---- Sankey color series (cluster number -> fixed color) ----
specified_colors <- c(
  "#FF0000", "#0000FF", "#998000", "#008000", "#800080",
  "#FFC0CB", "#000000", "#ADD8E6", "#FF4500", "#A52A2A",
  "#8B4513", "#FF8C00", "#32CD32", "#4682B4", "#9400D3",
  "#FFD700", "#C0C0C0", "#DC143C", "#1E90FF"
)


# ---- Sankey renderer (override) : align node colors with nodes$carac ----
render_author_sankey <- function(edges, nodes = NULL) {
  if (!requireNamespace("htmltools", quietly = TRUE)) return(NULL)

  if (!requireNamespace("networkD3", quietly = TRUE)) {
    return(htmltools::tags$div(class="small-note",
                              "networkD3 not installed; Sankey preview disabled. Use the SankeyMATIC code below."))
  }

  if (is.null(edges) || !is.data.frame(edges) || nrow(edges) == 0) {
    return(htmltools::tags$div(class="small-note", "No Sankey edges available."))
  }

  # normalize edges -> source/target/value
  e <- edges
  if (all(c("Leader","follower","WCD") %in% names(e))) {
    e2 <- dplyr::transmute(e, source = as.character(Leader), target = as.character(follower),
                          value = suppressWarnings(as.numeric(WCD)))
  } else if (all(c("Leader","Follower","WCD") %in% names(e))) {
    e2 <- dplyr::transmute(e, source = as.character(Leader), target = as.character(Follower),
                          value = suppressWarnings(as.numeric(WCD)))
  } else if (all(c("source","target","value") %in% names(e))) {
    e2 <- dplyr::transmute(e, source = as.character(source), target = as.character(target),
                          value = suppressWarnings(as.numeric(value)))
  } else {
    e2 <- e[,1:3, drop=FALSE]
    colnames(e2) <- c("source","target","value")
    e2$source <- as.character(e2$source); e2$target <- as.character(e2$target)
    e2$value <- suppressWarnings(as.numeric(e2$value))
  }
  e2 <- e2 %>% dplyr::filter(nzchar(source), nzchar(target), !is.na(value), value > 0)

  # nodes: use provided nodes (Top-20) if possible, else derive from edges
  if (is.null(nodes) || !is.data.frame(nodes) || nrow(nodes) == 0) {
    nd <- data.frame(name = sort(unique(c(e2$source, e2$target))), stringsAsFactors = FALSE)
    nd$carac <- "1"
  } else {
    nd <- nodes
    if (!("name" %in% names(nd))) nd$name <- as.character(nd[[1]])
    nd$name <- as.character(nd$name)

    if (!("carac" %in% names(nd))) {
      if ("cluster" %in% names(nd)) nd$carac <- nd$cluster
      else if ("group" %in% names(nd)) nd$carac <- nd$group
      else nd$carac <- "1"
    }
    nd$carac <- as.character(nd$carac)

    miss <- setdiff(unique(c(e2$source, e2$target)), nd$name)
    if (length(miss)) {
      nd <- rbind(nd, data.frame(name = miss, carac = "1", stringsAsFactors = FALSE))
    }
  }

  nd <- nd[!is.na(nd$name) & nzchar(nd$name), , drop=FALSE]
  nd <- nd[!duplicated(nd$name), , drop=FALSE]

  # index mapping
  idx <- seq_len(nrow(nd)) - 1L
  map <- setNames(idx, nd$name)
  links <- data.frame(
    source = unname(map[e2$source]),
    target = unname(map[e2$target]),
    value  = e2$value,
    stringsAsFactors = FALSE
  )
  links <- links[stats::complete.cases(links), , drop=FALSE]

  # colour scale aligned to carac (dark -> light blues)
  lv <- unique(nd$carac)
  lv_num <- suppressWarnings(as.integer(lv))
  if (all(!is.na(lv_num))) lv <- as.character(sort(lv_num)) else lv <- sort(lv)

  pal <- specified_colors
  if (all(!is.na(lv_num))) {
    cols <- vapply(as.integer(lv), function(k) pal[(k - 1L) %% length(pal) + 1L], character(1))
  } else {
    cols <- rep(pal, length.out = length(lv))
  }
  colourScale <- paste0(
    "d3.scaleOrdinal().domain([",
    paste0(sprintf('"%s"', lv), collapse = ","),
    "]).range([",
    paste0(sprintf('"%s"', cols), collapse = ","),
    "])"
  )

  networkD3::sankeyNetwork(
    Links = links,
    Nodes = nd,
    Source = "source",
    Target = "target",
    Value  = "value",
    NodeID = "name",
    NodeGroup = "carac",
    fontSize = 12,
    nodeWidth = 28,
    sinksRight = TRUE,
    colourScale = colourScale,
    width = NULL, height = 520
  )
}

# # ---- Prefer renderSSplot.R for SS panel drawing (REAL SSplot) ----
# Prefer newest uploaded renderer, then stable fallbacks.
# IMPORTANT: renderSSplot(81).R contains the one-AAC-per-cluster fix.
if (file.exists("renderSSplot(81).R")) {
  source("renderSSplot(81).R", local = FALSE, encoding = "UTF-8")
  cat("[BOOT] sourced: renderSSplot(81).R\\n")
} else if (file.exists("renderSSplot(80).R")) {
  source("renderSSplot(80).R", local = FALSE, encoding = "UTF-8")
  cat("[BOOT] sourced: renderSSplot(80).R\\n")
} else if (file.exists("renderSSplot(79).R")) {
  source("renderSSplot(79).R", local = FALSE, encoding = "UTF-8")
  cat("[BOOT] sourced: renderSSplot(79).R\\n")
} else if (file.exists("renderSSplot.R")) {
  source("renderSSplot.R", local = FALSE, encoding = "UTF-8")
  cat("[BOOT] sourced: renderSSplot.R\\n")
} else if (file.exists("renderSSplot(78).R")) {
  source("renderSSplot(78).R", local = FALSE, encoding = "UTF-8")
  cat("[BOOT] sourced: renderSSplot(78).R\\n")
} else {
  stop("renderSSplot.R not found in app folder")
}

# ---- lock SSplot renderer from renderSSplot.R (avoid any local overrides) ----
.render_panel_ss <- render_panel
if (!("font_scale" %in% names(formals(.render_panel_ss)))) {
  stop("renderSSplot.R render_panel() is not the expected version (font_scale missing).")
}

default_term <- "(Tsair-Wei Chien[Author]) AND (Taiwan[Affiliation])"

examples <- list(
  "Author + Country (Affiliation)" = '(Jane Doe[Author]) AND (Taiwan[Affiliation])',
  "Journal (TA = journal abbreviation)" = 'BMC Med Inform Decis Mak[TA] AND (Taiwan[Affiliation])',
  "Country in Affiliation (AD)" = 'Taiwan[AD] AND (machine learning[Title/Abstract])',
  "MeSH term (MH)" = '"Artificial Intelligence"[MH] AND Taiwan[AD]',
  "Title/Abstract keywords" = '"author collaboration"[Title/Abstract] OR "co-authorship"[Title/Abstract]',
  "Date range (PDAT)" = '("2018/01/01"[PDAT] : "2025/12/31"[PDAT]) AND Taiwan[AD]'
)


ui <- fluidPage(
  tags$head(tags$style(HTML("
    .small-note { font-size: 12px; color: #555; }
    .mono { font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, 'Liberation Mono', 'Courier New', monospace; }
    .card { border: 1px solid #e5e5e5; border-radius: 10px; padding: 14px; margin-bottom: 12px; background: #fff; }
    .nav-tabs > li > a[data-value='AAC'],
    .nav-tabs > li > a[data-value='BCa/EIS'],
    .nav-tabs > li > a[data-value='TAAA2']{ color:#d00000 !important; font-weight:700 !important; }
    .nav-tabs > li.active > a[data-value='AAC'],
    .nav-tabs > li.active > a[data-value='BCa/EIS'],
    .nav-tabs > li.active > a[data-value='TAAA2'],
    .nav-tabs > li.active > a[data-value='AAC']:focus,
    .nav-tabs > li.active > a[data-value='BCa/EIS']:focus,
    .nav-tabs > li.active > a[data-value='TAAA2']:focus{ background:#fff0f0 !important; border-top:3px solid #d00000 !important; color:#b00000 !important; }
    .bca-progress-note{margin-top:8px; padding:8px 10px; border-left:4px solid #f0ad4e; background:#fff8e8;}

	    /* Floating contact button (always visible, does not alter homepage layout) */
    .contact-fab{
      position: fixed;
      right: 18px;
      bottom: 18px;
      z-index: 9999;
      border-radius: 999px;
      padding: 10px 14px;
      background: #2d7ef7;
      color: #fff;
      border: 0;
      box-shadow: 0 6px 18px rgba(0,0,0,0.18);
    }
    .contact-fab:hover{ filter: brightness(0.95); }
    #cmc { background:#e6f4ff; border:2px solid #7bb6ff; }
  "))),

    tags$script(HTML("
      // ---- CMC: remember last value in localStorage (client-side only) ----
      (function(){
        const KEY = 'acsaas_cmc_last';
        function save(){
          const el = document.getElementById('cmc');
          if(!el) return;
          try{ localStorage.setItem(KEY, el.value || ''); }catch(e){}
        }
        function restore(){
          const el = document.getElementById('cmc');
          if(!el) return;
          let v = '';
          try{ v = localStorage.getItem(KEY) || ''; }catch(e){}
          if(v && !el.value){
            el.value = v;
            // also notify Shiny input binding
            if(window.Shiny && Shiny.setInputValue){
              Shiny.setInputValue('cmc', v, {priority:'event'});
            }
          }
        }
        document.addEventListener('DOMContentLoaded', function(){
          restore();
        // ---- helper: programmatically click buttons (used for Demo / autorun) ----
        if (window.Shiny) {
          Shiny.addCustomMessageHandler('jsClick', function(id){
            try { var el = document.getElementById(id); if (el) el.click(); } catch(e) {}
          });
        }
          const el = document.getElementById('cmc');
          if(el){
            el.addEventListener('input', save);
            el.addEventListener('change', save);
          }
        });
      })()
      // ---- TERM: remember last search term in localStorage ----
      (function(){
        const KEY = 'acsaas_term_last';
        function save(){
          const el = document.getElementById('term');
          if(!el) return;
          try{ localStorage.setItem(KEY, el.value || ''); }catch(e){}
        }
        function restore(){
          const el = document.getElementById('term');
          if(!el) return;
          let v = '';
          try{ v = localStorage.getItem(KEY) || ''; }catch(e){}
          if(v && !el.value){
            el.value = v;
            if(window.Shiny && Shiny.setInputValue){
              Shiny.setInputValue('term', v, {priority:'event'});
            }
          }
        }
        document.addEventListener('DOMContentLoaded', function(){
          restore();
          const el = document.getElementById('term');
          if(el){
            el.addEventListener('input', save);
            el.addEventListener('change', save);
          }
        });
      })();

      // ---- Contact authors: open a Shiny modal (no popup blockers) ----
      function openContactWindow(){
        if(window.Shiny && Shiny.setInputValue){
          Shiny.setInputValue('contact_btn', Date.now(), {priority:'event'});
        }
        return false;
      }
        const html = `
<!doctype html>
<html><head><meta charset='utf-8'><title>Contact authors / Request CMC</title>
<style>
  body{font-family:Arial, Helvetica, sans-serif; margin:18px; line-height:1.4;}
  h2{margin:0 0 10px 0;}
  .muted{color:#555; font-size:13px;}
  .grid{display:grid; grid-template-columns:repeat(auto-fit, minmax(240px, 1fr)); gap:12px; margin-top:12px;}
  .card{border:1px solid #ddd; border-radius:10px; padding:12px; background:#fff;}
  code{background:#f5f5f5; padding:2px 6px; border-radius:6px;}
  .cta{display:inline-block; padding:10px 14px; border-radius:10px; background:#0b6efd; color:#fff; text-decoration:none; font-weight:700;}
  .qr img{width:180px; height:180px; border:1px solid #eee; border-radius:8px;}
  .topbar{display:flex; justify-content:space-between; align-items:center;}
  .close{border:0; background:#eee; padding:6px 10px; border-radius:10px; cursor:pointer;}
</style>
</head><body>
<div class='topbar'>
  <h2>Contact authors / Request CMC</h2>
  <button class='close' onclick='window.close()'>Close</button>
</div>
<p>Please choose one of the following ways to contact the authors.</p>
<div class='grid'>
  <div class='card'>
    <h3>1. LINE account</h3>
    <p>LINE Official Account ID:</p>
    <p><code>@onq5657t</code></p>
    <a class='qr' href='https://line.me/R/ti/p/%40onq5657t' target='_blank' rel='noopener noreferrer'>
      <img src='https://api.qrserver.com/v1/create-qr-code/?size=200x200&data=https%3A%2F%2Fline.me%2FR%2Fti%2Fp%2F%2540onq5657t' alt='LINE QR code for @onq5657t'>
    </a>
    <p class='muted'>In LINE, choose “Add friends → QR code” and scan this image, or tap the QR code to open the add-friend page.</p>
  </div>

  <div class='card'>
    <h3>2. Email</h3>
    <p>You can send an email to either of the following addresses:</p>
    <ul class='muted'>
      <li><code>rasch.smile@gmail.com</code></li>
      <li><code>codingpaperabc@gmail.com</code></li>
    </ul>
    <p class='muted'>Copy and paste the address into your preferred email client or webmail service.</p>
  </div>

  <div class='card'>
    <h3>3. Join the ChatGPT group chatroom</h3>
    <p>We host a shared ChatGPT room where you can discuss TAAA / Basket / FLCA results together with the authors.</p>
    <ol class='muted'>
      <li>Copy this email address:<br><code>rasch.smile@gmail.com</code></li>
      <li>Send an email with the subject:<br><code>[Join TAAA FLCA ChatGPT Room]</code></li>
      <li>You will receive a reply containing the group-chat link.</li>
    </ol>
    <p class='muted'>Open the invite link in your browser and sign in with your ChatGPT account to enter the chatroom.</p>
  </div>

  <div class='card'>
    <p class='muted'>The authors personally cover the Google App Engine usage fees required to keep this app online. If this tool is helpful, voluntary donations are appreciated to support maintenance and further development.</p>
    <a class='cta' href='https://payment.ecpay.com.tw/QuickCollect/PayData?D50jzB3Lqk68BhoyYtTffB90EJpM0b4XmYFUwS4pAMI%3d' target='_blank' rel='noopener noreferrer'>🚀 Donate</a>
  </div>
</div>
</body></html>`;
        w.document.open();
        w.document.write(html);
        w.document.close();
        try{ w.focus(); }catch(e){}
      }
      window.openContactWindow = openContactWindow;
    ")),
  
    tags$script(HTML("
      // ---- URL params: cmc / title / term (fill only; NEVER autorun) ----
      (function(){
        function decodePlus(s){
          try { return decodeURIComponent((s||'').split('+').join(' ')); }
          catch(e){ return (s||'').split('+').join(' '); }
        }
        function setInput(id, val){
          if(val === null || typeof val === 'undefined') return;
          var v = decodePlus(val);
          var el = document.getElementById(id);
          if(el){
            el.value = v;
            el.dispatchEvent(new Event('input', {bubbles:true}));
            el.dispatchEvent(new Event('change', {bubbles:true}));
          }
          if(window.Shiny && Shiny.setInputValue){
            Shiny.setInputValue(id, v, {priority:'event'});
          }
        }
        var p = new URLSearchParams(window.location.search || '');
        if(!p) return;
        var cmc  = p.get('cmc');
        var term = p.get('term');
        var tit  = p.get('title');

        function tryUncheckUploaded(){
          // Prefer query over uploaded MEDLINE when URL provides `term`
          var cb = document.getElementById('use_uploaded_txt');
          if(cb){
            cb.checked = false;
            cb.dispatchEvent(new Event('change', {bubbles:true}));
          }
          if(window.Shiny && Shiny.setInputValue){
            Shiny.setInputValue('use_uploaded_txt', false, {priority:'event'});
          }
          return !!cb;
        }

        function applyParams(){
          // Fill only; never autorun
          if(cmc)  setInput('cmc', cmc);
          if(term) setInput('term', term);
          if(tit)  setInput('title', tit);

          // If term is provided, uncheck upload toggle (retry until the checkbox exists)
          if(term){
            var tries = 0;
            var iv = setInterval(function(){
              tries++;
              if(tryUncheckUploaded() || tries >= 50){ // ~5 seconds
                clearInterval(iv);
              }
            }, 100);
          }
        }

        if (document.readyState === 'loading') {
          document.addEventListener('DOMContentLoaded', applyParams);
        } else {
          applyParams();
        }
      })();")),
titlePanel("All-in-one Bibliometrics Dashboard (PubMed → Author/Country/Institute/MeSH)"),
  sidebarLayout(
    sidebarPanel(

      div(class="card",
          h4("Upload PubMed MEDLINE .txt"),
          fileInput("pubmed_txt", "Upload PubMed .txt (MEDLINE export)", accept = ".txt"),
          checkboxInput('use_uploaded_txt', "Use uploaded MEDLINE (ignore query)", value = TRUE),
          actionButton("btn_clear_uploaded", "Clear uploaded file"),
          tags$hr(),
          h4("AMA/PubMed/Google Scholar/NCKU Scopus journal extraction"),
          numericInput("ama_top_n", "Top N journals", value = 20, min = 1, max = 100, step = 1),
          actionButton("run_ama_journals", "Run AMA journal extraction", class = "btn-info", width = "100%"),
          tags$div(style="margin-top:8px;"),
          actionButton("run_ama_author_journal", "Run AMA author/journal extraction", class = "btn-warning", width = "100%"),
          tags$div(style="margin-top:8px;"),
          actionButton("run_ama_author_aac", "Run AMA author AAC (1st+Last only)", class = "btn-danger", width = "100%"),
          tags$div(class="small-note", "This extra run uses only the first and last authors from each reference/article; journals are excluded."),
          tags$div(style="margin-top:8px;"),
          checkboxInput("use_ama_textarea", "Run journal/authors from textarea below", value = FALSE),
          tags$div(class="small-note", "Select which pasted-reference parsers are allowed. Use one source at a time for the cleanest extraction; mixed mode is allowed but only normalized records are kept."),
          checkboxInput("ref_src_ama", "Parse AMA / PubMed MEDLINE / Vancouver references", value = TRUE),
          checkboxInput("ref_src_google", "Parse Google Scholar normalized 4-line records", value = TRUE),
          checkboxInput("ref_src_scopus", "Parse NCKU Pure / Scopus research output normalized records", value = TRUE),
          textAreaInput(
            "ama_refs_text",
            "Paste normalized AMA/PubMed, Google Scholar, or NCKU Pure/Scopus records for journal/author extraction",
            value = "",
            rows = 8,
            width = "100%",
            placeholder = "Paste normalized records only. Supported: (1) AMA/Vancouver references; (2) PubMed MEDLINE/NBIB records; (3) Google Scholar 4-line rows: Title / Authors / Journal-source / Cited-by+Year; (4) NCKU Pure/Scopus rows: Title / Authors+Year+於:Journal / 研究成果 / citation count after each reference. Non-normalized page text is ignored. When checked, this textarea is used instead of uploaded TXT."
          ),
          actionButton("btn_clear_ama_refs_text", "Clear pasted-reference textarea", class = "btn-default", width = "100%"),
          tags$div(class="small-note", "Clears only the pasted-reference textarea; uploaded TXT and PubMed query are unchanged."),
          tags$div(style="margin-top:8px;"),
          verbatimTextOutput("ama_processing_status"),
          downloadButton("dl_ama_journal_csv", "Download journal CSV", width = "100%"),
          downloadButton("dl_ama_journal_pie", "Download pie PNG", width = "100%"),
          tags$div(class="small-note","After upload, the file is saved to ./uploaded_pubmed.txt automatically. Then click Run AMA journal extraction to analyze journals from this text file."),
          tags$hr()
      ),
      div(class="card",
          h4("1) PubMed query"),
          textAreaInput("term", "Search term", value = default_term, rows = 3, width = "100%"),
          checkboxInput("batch_author_aac", "Series author AAC mode: treat each pasted name as an [Author] query", value = FALSE),
          tags$div(class="small-note", "When checked, paste one author per line in the search box. If the AMA/PubMed textarea option above is checked, the app uses that textarea instead and can extract author names from full references. AAC is computed only when PubMed hits are fully fetched (pubmed_hits <= fetched_pmids)."),
          downloadButton("dl_batch_author_aac", "Download series author AAC CSV", width = "100%"),
          tags$div(class="small-note",
                   "Common tags: ",
                   tags$span(class="mono","[Author]"), ", ",
                   tags$span(class="mono","[Affiliation]"), ", ",
                   tags$span(class="mono","[AD]"), ", ",
                   tags$span(class="mono","[TA]"), ", ",
                   tags$span(class="mono","[MH]"), ", ",
                   tags$span(class="mono","[Title/Abstract]"), ", ",
                   tags$span(class="mono","[PDAT]"), "."
          ),
          selectInput("example_pick", "Insert an example", choices = c("— choose —" = "", names(examples)), selected = ""),
          actionButton("load_example", "Insert", width = "100%")
      ),
      div(class="card",
          h4("2) Run"),
          numericInput("retmax", "retmax (PMIDs to fetch)", value = 1500, min = 10, max = 100000, step = 50),
	          passwordInput("cmc", "CMC or test (optional)", value = "", placeholder = "type test for 7-day trial, or enter 10-digit CMC"),
          textInput("title", "Title / Author name", value = ""),
          checkboxInput("save_xml", "Save pubmed.xml + pubmed_medline.txt", value = TRUE),
          numericInput("seed", "Layout seed", value = 1, min = 1, max = 99999, step = 1),
          sliderInput("label_size", "Label size (author dashboard)", min = 12, max = 40, value = 24),
          actionButton("run", "Fetch PubMed (Author + metadata)", class = "btn-primary", width = "100%"),
          # Demo button removed: demo is now auto via PubMed query when CMC invalid
          tags$div(class="small-note","Run domains one by one (faster + less freezing)."),
          actionButton("run_country",   "Run Country",        width="100%"),
          actionButton("run_stateprov", "Run State/Province", width="100%"),
          actionButton("run_inst",      "Run Institute",      width="100%"),
          actionButton("run_dept",      "Run Department",     width="100%"),
          actionButton("run_mesh",      "Run MeSH",           width="100%"),
          tags$hr(),
          downloadButton("dl_report", "Download HTML report", width = "100%"),
          downloadButton("dl_nodes", "Download nodes.csv", width = "100%"),
          downloadButton("dl_edges", "Download data_edges.csv", width = "100%")
      ),
      div(class="card",
          h4("3) Notes"),
          tags$ul(
            tags$li(tags$b("Author:"), " first author → last author (directed), FLCA + major sampling Top20."),
            tags$li(tags$b("Country/Institute/MeSH:"), " parsed from XML; each tab also uses FLCA + major sampling (one-link)."),
            tags$li(tags$b("MeSH filtering:"), " excludes demographics/general terms (Adult, Male, Humans, etc.) to keep medical terms focused."),
            tags$li(tags$b("Country:"), " dictionary-based; if not a country term, it is dropped (no city labels).")
          )
      )
    ),
    mainPanel(
      div(class="card",
          h4("Status"),
          verbatimTextOutput("log", placeholder = TRUE)
      ),
      tabsetPanel(id="main_tabs",
        tabPanel("Summary",

tags$hr(),
h3("Performance report PNG (8 domains)"),
fluidRow(
  column(4,
         actionButton("perf_use_auto", "Use fetched PubMed (auto)", class="btn-success", width="100%"),
         br(),
         fileInput("perf_file", "Upload Long Format CSV (Domain, Element, Value)"),
         downloadButton("perf_download_png", "Download performance PNG"),
         tags$p(class="small-note", "Tip: CSV needs columns Domain, Element, and Value (or value).")
  ),
  column(8,
         uiOutput("perf_plot_ui")
  )
),
tags$hr(),
                 h4("Summary report (Top 5 + AAC)"),
                 downloadButton("dl_summary_png", "Download summary PNG"),
                 br(), br(),
                 uiOutput("summary_btns"),
                 tags$hr(),
                 uiOutput("summary_status"),
                 DTOutput("tbl_summary_top10")

,
tags$hr(),
h4("Summary HTML report (preview)"),
downloadButton("dl_summary_html", "Download summary HTML"),
br(), br(),
div(style="border:1px solid #e6e6e6; border-radius:12px; padding:10px; background:#fff;",
    uiOutput("summary_html_preview")
)
        ),
        


        tabPanel("AMA Journals",
                 h4("Top 20 journals from uploaded AMA/PubMed or Google Scholar text"),
                 tags$div(class = "small-note",
                          "Upload a PubMed MEDLINE/AMA .txt file in the sidebar, then click Run AMA journal extraction. ",
                          "The parser uses MEDLINE JT/TA fields when available; otherwise it extracts the journal name from AMA-style references."),
                 tags$hr(),
                 h4("Summary"),
                 verbatimTextOutput("ama_journal_summary"),
                 tags$hr(),
                 h4("Google Scholar / NCKU Scopus citation metrics"),
                 tags$div(class = "small-note",
                          "Shown when pasted Google Scholar or NCKU Pure/Scopus rows include article citation counts and publication years. ",
                          "The Since column is computed from article publication year in the pasted rows."),
                 DTOutput("tbl_ama_gs_metrics"),
                 tags$hr(),
                 h4("Year bar chart from normalized Google Scholar / NCKU Scopus references"),
                 tags$div(class = "small-note",
                          "This chart uses only normalized parsed reference rows with a valid publication year."),
                 plotOutput("plt_ama_gs_year_bar", height = "420px"),
                 DTOutput("tbl_ama_gs_year_bar"),
                 tags$hr(),
                 h4("All-reference h-index from parsed Google Scholar / NCKU Scopus rows"),
                 tags$div(class = "small-note",
                          "This table uses all parsed article/reference citation counts once, not only first/last authors."),
                 DTOutput("tbl_ama_gs_all_reference_hindex"),
                 tags$hr(),
                 h4("Author h-index from Google Scholar / NCKU Scopus references"),
                 tags$div(class = "small-note",
                          "Two author-level h-index tables are reported separately: ",
                          "(1) first/last authors only, and (2) all byline authors. ",
                          "Both use the article citation counts parsed from the pasted Google Scholar rows."),
                 fluidRow(
                   column(6,
                          h4("1st/last author-based h-index"),
                          DTOutput("tbl_ama_gs_author_hindex_first_last")),
                   column(6,
                          h4("All-author based h-index"),
                          DTOutput("tbl_ama_gs_author_hindex_all"))
                 ),
                 tags$hr(),
                 h4("Top Journal Frequencies"),
                 DTOutput("tbl_ama_top_journals"),
                 tags$hr(),
                 h4("Pie Plot"),
                 plotOutput("plt_ama_journal_pie", height = "650px"),
                 tags$hr(),
                 tags$div(class = "small-note",
                          "This tab is journal-only. For the mixed Journal + Author co-word network, SSplot, and cluster AAC, use the AMA Author/Journal tab."),
                 tags$hr(),
                 h4("Extracted journal by reference"),
                 DTOutput("tbl_ama_journal_refs"),
                 tags$hr(),
                 h4("Article citations parsed from Google Scholar / NCKU Scopus"),
                 DTOutput("tbl_ama_gs_articles")
        ),

        tabPanel("AMA Author/Journal",
                 h4("Top 20 author/journal elements from uploaded AMA/PubMed or Google Scholar records"),
                 tags$div(class = "small-note",
                          "Upload a PubMed MEDLINE/AMA .txt file, then click Run AMA author/journal extraction. ",
                          "The app extracts journals and authors, builds Journal + Author co-occurrence edges, applies Top20 selection, and draws the same Top20 using the FLCA-process plotting logic. Cluster AAC is shown once per cluster in the SSplot."),
                 tags$hr(),
                 verbatimTextOutput("ama_ref_status"),
                 tags$hr(),
                 h4("All-reference h-index"),
                 tags$div(class = "small-note",
                          "This reports profile-reported Google Scholar metrics when the copied profile summary is present, plus parsed normalized reference metrics once per article. It is independent of byline position."),
                 DTOutput("tbl_ama_ref_all_reference_hindex"),
                 tags$hr(),
                 h4("Year bar chart from normalized Google Scholar / NCKU Scopus references"),
                 tags$div(class = "small-note",
                          "This chart uses only normalized parsed reference rows with a valid publication year."),
                 plotOutput("plt_ama_ref_year_bar", height = "420px"),
                 DTOutput("tbl_ama_ref_year_bar"),
                 tags$hr(),
                 h4("Reference-level author h-index"),
                 tags$div(class = "small-note",
                          "When Google Scholar or NCKU Pure/Scopus rows are pasted, this reports author-level h-index from the same article citation counts: first/last authors only vs all byline authors."),
                 fluidRow(
                   column(6,
                          h4("1st/last author-based h-index"),
                          DTOutput("tbl_ama_ref_author_hindex_first_last")),
                   column(6,
                          h4("All-author based h-index"),
                          DTOutput("tbl_ama_ref_author_hindex_all"))
                 ),
                 tags$hr(),
                 h4("Top 20 check list"),
                 tags$div(class = "small-note", "This table is the FLCA-MA-SIL Top20 used by Network, SSplot, and Kano. Numeric values are displayed to 2 decimals."),
                 verbatimTextOutput("tbl_ama_ref_top20_check"),
                 tags$hr(),
                 h4("Interactive network"),
                 tags$div(style = "margin-bottom:8px;",
                          downloadButton("dl_ama_ref_network_png", "Download Network PNG"),
                          downloadButton("dl_ama_ref_network_html", "Download Interactive Network HTML")),
                 visNetworkOutput("vn_ama_ref", height = "760px"),
                 tags$hr(),
                 h4("Clustered network from FLCA-MA-SIL"),
                 tags$div(class = "small-note",
                          "This second network fixes node positions by FLCA cluster (carac), so the clusters used in SSplot can be checked directly in the network view."),
                 tags$div(style = "margin-bottom:8px;",
                          downloadButton("dl_ama_ref_cluster_network_png", "Download Cluster Network PNG"),
                          downloadButton("dl_ama_ref_cluster_network_html", "Download Cluster Network HTML")),
                 visNetworkOutput("vn_ama_ref_cluster", height = "760px"),
                 tags$hr(),
                 h4("SSplot from FLCA process / renderSSplot.R"),
                 tags$div(style = "margin-bottom:8px;",
                          downloadButton("dl_ama_ref_ss_png", "Download SSplot PNG")),
                 plotOutput("plt_ama_ref_ss", height = "900px"),
                 tags$hr(),
                 h4("Kano plot from FLCA process / kano.R"),
                 tags$div(style = "margin-bottom:8px;",
                          downloadButton("dl_ama_ref_kano_png", "Download Kano PNG")),
                 plotOutput("plt_ama_ref_kano", height = "1250px"),
                 tags$hr(),
                 h4("Top 20 nodes"),
                 verbatimTextOutput("tbl_ama_ref_nodes"),
                 tags$hr(),
                 h4("Edges"),
                 verbatimTextOutput("tbl_ama_ref_edges"),
                 tags$hr(),
                 h4("Parsed Journal/Author table"),
                 verbatimTextOutput("tbl_ama_ref_wide"),
                 tags$hr(),
                 downloadButton("dl_ama_ref_nodes", "Download AMA author/journal nodes CSV"),
                 downloadButton("dl_ama_ref_edges", "Download AMA author/journal edges CSV"),
                 downloadButton("dl_ama_ref_wide", "Download parsed journal/authors CSV")
        ),

        tabPanel("AMA Author AAC",
                 h4("AMA Author AAC: 1st and last authors only"),
                 tags$div(class = "small-note",
                          "Click the red sidebar button: Run AMA author AAC (1st+Last only). ",
                          "This tab excludes journal nodes and keeps only first/last author terms per article. ",
                          "The same FLCA-MA-SIL Top20 logic is then used for network, clustered network, SSplot, Kano, and AAC."),
                 tags$hr(),
                 verbatimTextOutput("ama_author_aac_status"),
                 tags$hr(),
                 h4("Highlighted h-index and AAC"),
                 tags$div(class = "small-note",
                          "The upper rows summarize profile-reported, all-reference, and first/last author-based citation metrics. ",
                          "AAC is computed from the Top20 first/last-author node values."),
                 DTOutput("tbl_ama_author_aac_highlights"),
                 tags$hr(),
                 h4("All-reference Google Scholar / NCKU h-index"),
                 tags$div(class = "small-note",
                          "This additional table uses all parsed references and citation counts, including NCKU Pure/Scopus citation counts found after each reference; it is not limited to first/last-author rows."),
                 DTOutput("tbl_ama_author_aac_all_reference_hindex"),
                 tags$hr(),
                 h4("1st/last author-based Google Scholar / NCKU h-index"),
                 tags$div(class = "small-note",
                          "This table uses the same parsed reference citation counts, but assigns each reference only to its first and last authors."),
                 DTOutput("tbl_ama_author_aac_first_last_reference_hindex"),
                 tags$hr(),
                 h4("Top 20 first/last-author check list"),
                 verbatimTextOutput("tbl_ama_author_aac_top20_check"),
                 tags$hr(),
                 h4("Interactive network"),
                 tags$div(style = "margin-bottom:8px;",
                          downloadButton("dl_ama_author_aac_network_png", "Download Network PNG"),
                          downloadButton("dl_ama_author_aac_network_html", "Download Interactive Network HTML")),
                 visNetworkOutput("vn_ama_author_aac", height = "760px"),
                 tags$hr(),
                 h4("Clustered network from FLCA-MA-SIL"),
                 tags$div(class = "small-note", "Node positions are fixed by carac; each follower has one edge to its own cluster leader."),
                 tags$div(style = "margin-bottom:8px;",
                          downloadButton("dl_ama_author_aac_cluster_network_png", "Download Cluster Network PNG"),
                          downloadButton("dl_ama_author_aac_cluster_network_html", "Download Cluster Network HTML")),
                 visNetworkOutput("vn_ama_author_aac_cluster", height = "760px"),
                 tags$hr(),
                 h4("SSplot"),
                 tags$div(style = "margin-bottom:8px;",
                          downloadButton("dl_ama_author_aac_ss_png", "Download SSplot PNG")),
                 plotOutput("plt_ama_author_aac_ss", height = "900px"),
                 tags$hr(),
                 h4("Kano plot"),
                 tags$div(style = "margin-bottom:8px;",
                          downloadButton("dl_ama_author_aac_kano_png", "Download Kano PNG")),
                 plotOutput("plt_ama_author_aac_kano", height = "1250px"),
                 tags$hr(),
                 h4("Top20 first/last-author nodes"),
                 verbatimTextOutput("tbl_ama_author_aac_nodes"),
                 tags$hr(),
                 h4("First-last author edges"),
                 verbatimTextOutput("tbl_ama_author_aac_edges"),
                 tags$hr(),
                 h4("Parsed first/last-author table"),
                 verbatimTextOutput("tbl_ama_author_aac_wide"),
                 tags$hr(),
                 downloadButton("dl_ama_author_aac_nodes", "Download first/last-author nodes CSV"),
                 downloadButton("dl_ama_author_aac_edges", "Download first/last-author edges CSV"),
                 downloadButton("dl_ama_author_aac_wide", "Download parsed first/last-author CSV")
        ),

  tabPanel("Author",
                 h4("Interactive network (Author, Top20)"),
                 visNetworkOutput("vn_author", height = "650px"),
                 tags$hr(),
                 h4("Sampled nodes / edges"),
                 DTOutput("tbl_author_nodes"),        DTOutput("tbl_author_edges")
        ),

        tabPanel("Journal/Year",
                 h4("Journal / Year (Top20)"),
                 visNetworkOutput("vn_journal_year", height = "520px"),
                 DTOutput("tbl_journal_year")
        ),

        tabPanel("Year/Articles",
                 h4("Publications by Year"),
                 plotOutput("plt_year_articles", height = "520px"),
                 DTOutput("tbl_year_articles")
        )

        , tabPanel("Slope",

      tags$hr(),
      h4("Combo entity Top10 slope graph over years"),
      tags$div(class = "small-note", "Top 10 selected entity elements by yearly counts."),
      selectInput("slope_simple_entity_domain", "Entity/domain", 
                  choices = c("Author","Journal","Country","State/Province","Institute","Department","MeSH"),
                  selected = "Author"),
      sliderInput("slope_simple_recent_years", "Recent years", min = 5, max = 20, value = 10, step = 1),
      verbatimTextOutput("slope_combo_entity_note"),
      plotOutput("slope_combo_entity_plot", height = "760px"),
      DTOutput("slope_combo_entity_table"),

         h4("Top10 各名稱年度次數 / Slopegraph (10年門檻 + 前5/後5 t-test；single-term only)"),
         uiOutput("yearcount_ui"),
         tags$hr(),
         h4("Combo entity slopegraph from app(708).R style"),
         tags$p(class = "small-note",
                "Choose one entity domain. The graph shows a real Tufte-style slopegraph for the Top 10 elements of that entity over recent years."),
         fluidRow(
           column(4,
                  selectInput("slope_combo_entity_domain",
                              "Entity domain",
                              choices = c("Author","Journal","Country","State/Province","Institute","Department","MeSH"),
                              selected = "Author")
           ),
           column(4,
                  sliderInput("slope_combo_recent_years",
                              "Recent years",
                              min = 5, max = 20, value = 10, step = 1)
           ),
           column(4,
                  tags$div(class = "small-note", style = "margin-top:25px;",
                           "This selector alone controls the real slopegraph; no Combo-tab sync is used.")
           )
         ),
         plotOutput("plot_slope_combo_entity_top10", height = "760px"),
         DT::DTOutput("tbl_slope_combo_entity_top10"),
         tags$hr(),
         h4("Metadata combo slopegraph"),
         plotOutput("plt_slope_top2_country", height = "520px"),
         DT::DTOutput("tbl_yearcount_top20"),
         tags$hr(),
         h4("Author metadata slopegraph"),
         plotOutput("plt_slope_top2_author", height = "520px"),
         DT::DTOutput("tbl_yearcount_top20_author"),
         tags$hr(),
         h4("MeSH metadata slopegraph"),
         plotOutput("plt_slope_top2_mesh", height = "520px"),
         DT::DTOutput("tbl_yearcount_top20_mesh"),

         tags$hr(),
         h4("Year-count matrix inside Slope tab"),
         tags$p(class = "small-note",
                "Uses the same selected entity/domain and recent-year window as the Combo entity slopegraph above. Highlight point = Count greater than that element's own mean count."),
         plotOutput("slope_burst_heatmap_plot", height = "650px"),

         tags$hr(),
         h4("Year-count spot plot inside Slope tab"),
         tags$p(class = "small-note",
                "Top10 long table is generated first, then plotted as a year-by-element spot timeline. Red spots indicate Count greater than the element mean."),
         plotOutput("slope_burst_timeline_plot", height = "650px"),

         tags$hr(),
         h4("Strategic map inside Slope tab"),
         tags$p(class = "small-note",
                "Top10 elements from the same selected entity/domain. x=maturity(value), y=influence(value2), point size/label reflects recent trend from recent timepoints."),
         plotly::plotlyOutput("slope_mesh_3d_plot", height = "720px"),

         tags$hr(),
         h4("Top10 long table used by Slope-tab extra plots"),
         verbatimTextOutput("slope_burst_long_table")
),
        tabPanel("Term/Year",
                 h4("1st+Last Authors: Terms / Year counts (slopegraph)"),
                 uiOutput("yearcount_ui_fl"),
                 plotOutput("plt_slope_top2_author_fl", height = "520px"),
                 DT::DTOutput("tbl_yearcount_top20_fl")
        ),

        tabPanel("Sankey",
                 h4("Sankey (Author)"),
                 tags$p("Shows an in-app Sankey when available; always provides a Sankeymatic link + code for copy/paste."),
                 fluidRow(
                   column(4,
                          selectInput("sankey_domain", "Domain (Top-20 nodes)",
                                      choices = c("Author","Journal/Year","Country","Institute","Department","MeSH"),
                                      selected = "Author")
                   ),
                   column(8, tags$div(style="padding-top:28px;", ""))
                 ),
                 uiOutput("ui_sankey_author"),
                 tags$p(tags$strong("Sankeymatic URL")),
                 uiOutput("ui_sankey_link"),
                 tags$p(tags$strong("Sankeymatic text (copy/paste)")),
                 tags$pre(style="white-space:pre-wrap;", textOutput("txt_sankey_code"))
        ),

        tabPanel("SSplot",
                 h4("SSplot (Silhouette panel)"),
                 fluidRow(
                   column(4,
                          selectInput("ss_domain", "Domain (Top-20 nodes)",
                                      choices = c("Author","Journal/Year","Country","Institute","Department","MeSH"),
                                      selected = "Author")
                   ),
                   column(4,
                          sliderInput("ss_font_scale", "Font scale (SSplot panel)",
                                      min = 0.6, max = 2.2, value = 1.3, step = 0.1)
                   ),
                   column(4, tags$div(style="padding-top:28px;", ""))
                 ),
                 plotOutput("ssplot_panel", height = "820px")
        ),



        tabPanel("SS kano",
                 h4("SS Kano (Silhouette SS vs a*)"),
                 fluidRow(
                   column(4,
                          selectInput("ss_kano_domain", "Domain (Top-20 nodes)",
                                      choices = c("Author","Journal/Year","Country","Institute","Department","MeSH"),
                                      selected = "Author")
                   ),
                   column(4,
                          sliderInput("ss_kano_label_size", "Label size (SS Kano)",
                                      min = 2, max = 10, value = 4, step = 0.5)
                   ),
                   column(4, tags$div(style="padding-top:28px;", "Kano: x=a*, y=SS; bubble size=value; color=cluster."))
                 ),
                 plotOutput("ss_kano_plot", height = "1050px")
        ),

        tabPanel("TAAA",
          h4("TAAA: row-level representative topic (MeSH domain)"),
          tags$p(class = "small-note",
                 "Each article row is mapped to the dominant FLCA cluster among its MeSH terms. The representative theme name is the within-cluster leader term (highest value). If a Profile column is available in future uploads, agreement tables can also be shown."),
          tableOutput("taaa_table"),
          tags$hr(),
          h4("Frequency table of article themes"),
          tableOutput("taaa_theme_freq_table"),
          tags$hr(),
          h4("Profile vs article theme mapping"),
          tableOutput("taaa_profile_map_table"),
          tags$hr(),
          h4("Kappa summary after Hungarian mapping"),
          tableOutput("taaa_kappa_table"),
          tags$hr(),
          h4("K x K table after Hungarian mapping"),
          tableOutput("taaa_confusion_table")
        ),
        tabPanel("TAAA2",
          h4("TAAA2: semantic article clustering using PubMedBERT / SPECTER2"),
          tags$p(class = "small-note",
                 "TAAA2 extends the TAAA idea from rule-based article-theme mapping to embedding-based semantic clustering. Titles and abstracts are converted into PubMedBERT or SPECTER2 embeddings, clustered by cosine similarity, labeled by TF-IDF terms, and visualized as UMAP and semantic network plots."),
          fluidRow(
            column(4,
                   selectInput(
                     "taaa2_model",
                     "Embedding model",
                     choices = c(
                       "PubMedBERT / BiomedBERT" = "microsoft/BiomedNLP-BiomedBERT-base-uncased-abstract-fulltext",
                       "SPECTER2" = "allenai/specter2"
                     ),
                     selected = "microsoft/BiomedNLP-BiomedBERT-base-uncased-abstract-fulltext"
                   )
            ),
            column(3,
                   numericInput("taaa2_k", "Number of clusters", value = 3, min = 2, max = 8, step = 1)
            ),
            column(3,
                   numericInput("taaa2_top_edges", "Top links per article", value = 2, min = 1, max = 5, step = 1)
            ),
            column(2,
                   tags$div(style = "padding-top:25px;",
                            actionButton("run_taaa2", "Run TAAA2", class = "btn-danger"))
            )
          ),
          tags$hr(),
          h4("TAAA2 UMAP semantic cluster map"),
          plotOutput("taaa2_umap_plot", height = "560px"),
          tags$hr(),
          h4("TAAA2 semantic similarity network"),
          plotOutput("taaa2_network_plot", height = "660px"),
          tags$hr(),
          h4("Article-level TAAA2 results"),
          tableOutput("taaa2_article_table"),
          tags$hr(),
          h4("Cluster labels and element terms"),
          tableOutput("taaa2_cluster_label_table")
        ),
        tabPanel("Lotka",
                 h4("Lotka's law for author publication distribution"),
                 tags$p(class = "small-note",
                        "Based on author publication counts from the current PubMed run. Left: log-log Lotka plot. Right: observed vs expected frequencies. The tables report the chi-square goodness-of-fit test."),
                 plotOutput("lotka_plot", height = "460px"),
                 tags$hr(),
                 h4("Lotka chi-square summary"),
                 tableOutput("lotka_test_table"),
                 tags$hr(),
                 h4("Observed vs expected table"),
                 tableOutput("lotka_observed_expected_table"),
                 tags$hr(),
                 h4("Merged tail table used for chi-square test"),
                 tableOutput("lotka_chisq_table")
        ),

        
        tabPanel("Summary",
                 fluidRow(
                   column(12,
                          h4("Summary report (Top 5 + AAC)"),
                          downloadButton("download_summary_png", "Download summary PNG"),
                          tags$hr()
                   )
                 ),
                 fluidRow(
                   column(6,
                          div(class="card",
                              h4("Country"),
                              tags$div(style="font-weight:bold;", textOutput("aac_country")),
                              DTOutput("sum_country")
                          ),
                          div(class="card",
                              h4("Institute"),
                              tags$div(style="font-weight:bold;", textOutput("aac_inst")),
                              DTOutput("sum_inst")
                          ),
                          div(class="card",
                              h4("Department"),
                              tags$div(style="font-weight:bold;", textOutput("aac_dept")),
                              DTOutput("sum_dept")
                          ),
                          div(class="card",
                              h4("Author"),
                              tags$div(style="font-weight:bold;", textOutput("aac_author")),
                              DTOutput("sum_author")
                          )
                   ),
                   column(6,
                          div(class="card",
                              h4("Journal"),
                              tags$div(style="font-weight:bold;", textOutput("aac_journal")),
                              DTOutput("sum_journal")
                          ),
                          div(class="card",
                              h4("Year"),
                              tags$div(style="font-weight:bold;", textOutput("aac_year")),
                              DTOutput("sum_year")
                          ),
                          div(class="card",
                              h4("State/Province"),
                              tags$div(style="font-weight:bold;", textOutput("aac_stateprov")),
                              DTOutput("sum_stateprov")
                          ),
                          div(class="card",
                              h4("MeSH Term"),
                              tags$div(style="font-weight:bold;", textOutput("aac_mesh")),
                              DTOutput("sum_mesh")
                          )
                   )
                 )
        ),
  tabPanel("AAC",
    aac_ui("aac")   # <<<<<< 這行就是插入點
  ),
  tabPanel("BCa/EIS",
    h4("BCa/EIS: Bootstrap BCa confidence interval for elbow identification"),
    tags$p(class = "small-note",
           "BCa/EIS scans ranked author publication counts using conditional follower bootstrap. For each candidate, the candidate value is fixed and only lower-ranked followers are resampled. AAC_OR = (candidate/follower1)/(follower1/follower2); BCa/EIS judges whether the BCa lower bound exceeds the AAC_OR cutoff."),
    fluidRow(
      column(3, numericInput("bca_R", "Bootstrap replications R", value = 500, min = 100, max = 5000, step = 100)),
      column(3, numericInput("bca_cutoff", "AAC OR cutoff", value = 2.0, min = 1, max = 10, step = 0.1)),
      column(3, numericInput("bca_conf", "Confidence level", value = 0.95, min = 0.80, max = 0.99, step = 0.01)),
      column(3, selectInput("bca_top_k", "Dominance target", choices = c("Top 1" = 1, "Top 2" = 2, "Top 3" = 3), selected = 1))
    ),
    fluidRow(
      column(4, tags$div(style = "padding-top:10px;", downloadButton("dl_bca_eis", "Download BCa/EIS CSV", width = "100%"))),
      column(4, tags$div(style = "padding-top:10px;", downloadButton("dl_bca_topk_png", "Download Top-k PNG", width = "100%"))),
      column(4, tags$div(style = "padding-top:10px;", downloadButton("dl_bca_topk_figure_zip", "Download PNG + Caption ZIP", width = "100%")))
    ),
    tags$hr(),
    h4("Editable values for BCa/EIS demonstration or re-analysis"),
    tags$p(class = "small-note",
           "Paste any entity counts or scores below (Author, Country, Institute, Journal, MeSH, Department, etc.). Values can be separated by new lines, commas, spaces, or semicolons. They are auto-sorted decreasingly before BCa/EIS; no manual sorting is required."),
    textAreaInput(
      "bca_values_text",
      "Editable values",
      value = "1000
100
50
49
48
47
46
45
44
43
42
41
40
39
38
37
36
35
34
33
32
31
30
29
28
27
26
25
24
23",
      rows = 7,
      width = "100%"
    ),
    actionButton("bca_recompute_manual", "Recompute BCa/EIS from editable values", class = "btn-warning"),
    br(), br(),
    verbatimTextOutput("txt_bca_run_status"),
    tags$div(class = "bca-progress-note small-note",
             "After pressing Recompute, please wait until the progress bar finishes. Large R or many values can be slow because BCa is re-estimated at every internal point."),
    tags$p(class = "small-note",
           "Note: BCa/EIS now uses conditional follower bootstrap: each candidate value is fixed, and only its lower-ranked followers are resampled. Rows with fewer than 3 followers are marked with the reason and are not assigned BCa/EIS."),
    tags$hr(),
    h4("Automatic elbow decision and EIS"),
    verbatimTextOutput("txt_bca_eis_decision"),
    DT::DTOutput("tbl_bca_eis_summary"),
    tags$p(class = "small-note",
           "EIS is defined as BCa_lower / AAC_OR cutoff for the automatically selected elbow. EIS > 1 means the BCa lower bound exceeds the AAC_OR cutoff and the elbow is statistically supported."),
    tags$hr(),
    plotOutput("plt_bca_eis", height = "560px"),
    tags$hr(),
    DT::DTOutput("tbl_bca_eis"),
    tags$hr(),
    h4("Top-k dominance vs rest"),
    tags$p(class = "small-note",
           "Choose Top 1, Top 2, or Top 3 above. The selected candidate rank is fixed, and only its lower-ranked followers are resampled. This allows AAC/BCa/EIS testing for different entity elements, not only authors."),
    verbatimTextOutput("txt_top1_bca_eis"),
    plotOutput("plt_top1_dominance", height = "620px"),
    DT::DTOutput("tbl_top1_bca_eis")
  ),
tabPanel("Download List",
                 h4("Download nodes & edges"),
                 downloadButton("dl_nodes_edges_zip", "Download ZIP (nodes.csv + edges.csv)"),
                 tags$p(class="small-note", "Includes Top20 nodes/edges after FLCA + Major Sampling.")
        ),

        
        
        tabPanel("WebZIP",
                 h4("WebZIP (offline package)"),
                 tags$p("Wrap the current download/ folder into a single ZIP for offline viewing (index.html + figures + tables)."),
                 actionButton("btn_refresh_webzip", "Refresh index.html + file list"),
                 br(), br(),
                 downloadButton("dl_webzip", "Download WebZIP (download/)"),
                 br(), br(),
                 fluidRow(
                   column(4,
                          h4("Files in download/"),
                          uiOutput("webzip_file_list")
                   ),
                   column(8,
                          h4("Offline index.html preview"),
                          uiOutput("webzip_index_preview"),
                          br(),
                          h4("ZIP manifest (what will be included)"),
                          tableOutput("webzip_manifest")
                   )
                 )
        ),

tabPanel("IP",
                 h4("IP pass status"),
                 textOutput("ip_pass_status"),
                 uiOutput("ip_access_mode_ui"),
                 verbatimTextOutput("ip_status"),
                 tags$div(class="small-note",
                          "Access is allowed if your IP is in IPlist.txt, or if CMC/test validation passes."),
                 DTOutput("ip_allowlist")
        ),
tabPanel("Report",
      tableOutput("tbl_report_domains"),
                 h4("HTML report"),
                 downloadButton("dl_report_html", "Download HTML report"),
                 tags$p(class="small-note", "Self-contained HTML with key tables and plots.")
        ),


        tabPanel("Most cited",
                 h4("NIH iCite (Most cited / RCR analysis)"),
                 tags$p(class="small-note",
                        "Opens NIH iCite analysis using the PMIDs fetched in this run."),
                 uiOutput("ui_icite_link"),
                 tags$hr(),
                 h5("PMIDs used"),
                 downloadButton("dl_pmids_txt", "Download PMID list (.txt)"),
                 downloadButton("dl_pmids_csv", "Download PMID list (.csv)"),
                 tags$p(class="small-note", "PMIDs are taken directly from the latest PubMed query result or uploaded MEDLINE file used in this run."),
                 DTOutput("tbl_pmids")
        ),

        tabPanel("ReadMe",
                 h4("ReadMe"),
                 htmlOutput("readme_html")
        ),

        tabPanel("Kano",
                 h4("Kano plots"),
                 fluidRow(
                   column(4,
                          selectInput("kano_domain", "Domain (Top-20 nodes)",
                                      choices = c("Author","Journal/Year","Country","Institute","Department","MeSH"),
                                      selected = "Author")
                   ),
                   column(4,
                          sliderInput("kano_label_size", "Label size (Kano plots)",
                                      min = 2, max = 10, value = 4, step = 0.5)
                   ),
                   column(4,
                          tags$div(style="padding-top:28px; font-weight:600;",
                                   textOutput("kano_aac_header"))
                   )
                 ),
                 tags$p(class="small-note",
                        "Kano plots are drawn from the selected domain's Top-20 nodes after FLCA+MajorSampling, using value/value2/SS/a* for that domain."),
                 h5("Kano: value vs value2"),
                 plotOutput("kano_vv2_plot", height = "1050px"),
                 tags$hr(),
                 h5("Kano: SS vs a*"),
                 plotOutput("kano_ss_astar_plot", height = "1050px")
        ),

        tabPanel("Country",
                 h4("Interactive network (Country, Top20)"),
                 visNetworkOutput("vn_country", height = "650px"),
                 DTOutput("tbl_country_nodes"),
                 DTOutput("tbl_country_edges")
        ),
        
tabPanel("Map",
         h4("World map (Country nodes, pre-FLCA)"),
         tags$p(class="small-note",
                "This map uses Country nodes BEFORE FLCA (distribution only; edges not required)."),
         plotOutput("country_map", height = "650px")
),
tabPanel("State/Province",
                 h4("Interactive network (State/Province, Top20)"),
                 tags$p(class="small-note",
                        "For US: states; for China: provinces; all other nodes keep the country name (US/UK abbreviated in Country tab)."),
                 visNetworkOutput("vn_stateprov", height = "650px"),
                 DTOutput("tbl_stateprov_nodes"),
                 DTOutput("tbl_stateprov_edges")
        ),
        
        tabPanel("China",
                 h4("China (Provinces / Municipalities / SARs)"),
                 tags$p(class="small-note",
                        "Built from State/Province term-list before FLCA. Includes Mainland China provinces, municipalities, Taiwan, Hong Kong, Macau, and South China Sea Islands."),
                 # Use highcharter output directly (avoids htmlwidgets::renderWidget)
                 uiOutput("china_map"),
                 DTOutput("tbl_china_counts")
        ),

        tabPanel("USA",
                 h4("USA (States)"),
                 tags$p(class="small-note",
                        "Built from State/Province term-list before FLCA (USA only)."),
                 if (requireNamespace("plotly", quietly = TRUE)) plotly::plotlyOutput("usa_map", height = "650px") else uiOutput("usa_map_placeholder"),
                 DTOutput("tbl_usa_counts")
        ),

tabPanel("Institute",
                 h4("Interactive network (Institute, Top20)"),
                 visNetworkOutput("vn_inst", height = "650px"),
                 DTOutput("tbl_inst_nodes"),
                 DTOutput("tbl_inst_edges")
        ),

        tabPanel("Department",
                 h4("Interactive network (Department, Top20)"),
                 visNetworkOutput("vn_dept", height = "650px"),
                 DTOutput("tbl_dept_nodes"),
                 DTOutput("tbl_dept_edges")
        ),

        tabPanel("MeSH",
                 h4("Interactive network (MeSH, Top20; demographics removed)"),
                 visNetworkOutput("vn_mesh", height = "650px"),
                 DTOutput("tbl_mesh_nodes"),
                 DTOutput("tbl_mesh_edges")
        )
        ,

        tabPanel("Combo",
          fluidRow(
            column(3,
              wellPanel(
                selectInput("combo_domain", "Choose metadata domain", 
                            choices = c("Author","Journal/Year","Country","State/Province","Institute","Department","MeSH"),
                            selected = "Author"),
                helpText("Uses the existing Top20 nodes/edges from each domain tab; this view just lets you switch domains in one place.")
              )
            ),
            column(9,
              visNetworkOutput("vn_combo", height = "650px"),
              tags$hr(),
              h4("Nodes (Top20)"),
              DTOutput("combo_nodes"),
              tags$hr(),
              h4("Edges (Top20)"),
              DTOutput("combo_edges")
            )
          )
        ),

        tabPanel("Chord",
          h4("Chord diagram (from the selected domain Top20 edges)"),
          fluidRow(
            column(4,
              selectInput("chord_domain", "Choose metadata domain (Chord)",
                          choices = c("Author","Journal/Year","Country","State/Province","Institute","Department","MeSH"),
                          selected = "Author")
            ),
            column(4,
              checkboxInput("chord_follow_combo", "Sync with Combo selection", value = TRUE)
            )
          ),
          tags$p(class="small-note",
                 "If 'chorddiag' is installed, an interactive chord is shown; otherwise a static 'circlize' chord is used. Colors follow node clusters (carac/cluster)."),
          uiOutput("chord_ui"),
          tags$hr(),
          verbatimTextOutput("chord_debug")
        ),


        tabPanel("IP",
          tags$div(class="card",
            h4("Access status"),
            tableOutput("ip_status"),
            tags$p(class="small-note", "Access types: IP pass / Trial / CMC pass.")
          ),
          tags$div(class="card",
            h4("Allowlist (iplist.txt)"),
            DTOutput("ip_allowlist"),
            tags$p(class="small-note", "If iplist.txt exists and contains IPs, only allowlisted IPs can use the app.")
          ),
          tags$div(class="card",
            h4("IP log (ip.txt)"),
            DTOutput("ip_log"),
            tags$p(class="small-note", "Tracks IP, recent_date, total_count, and last CMC used.")
          )
        ),

    ),
    actionButton("contact_btn", "Contact authors / Request CMC", class = "contact-fab")
  )
)
)



  .term_stats <- function(term_list){
    # term_list: list(article -> character vector of terms)
    if (is.null(term_list) || !length(term_list)) return(list(max_items=0L, median_items=0, n_ge2=0L))
    lens <- vapply(term_list, function(x){
      x <- as.character(x); x <- trimws(x); x <- x[nzchar(x)]
      length(unique(x))
    }, integer(1))
    if (!length(lens)) return(list(max_items=0L, median_items=0, n_ge2=0L))
    list(
      max_items = if (length(lens)) max(lens) else 0L,
      median_items = if (length(lens)) as.numeric(stats::median(lens)) else 0,
      n_ge2 = sum(lens >= 2, na.rm = TRUE)
    )
  }

# ---- Dual-circle Kano plot helpers (base R; no extra deps) ----
.plot_kano_dual_base <- function(df, xcol, ycol, sizecol=NULL, groupcol=NULL,
                                 title_txt="", xlab=NULL, ylab=NULL,
                                 labelcol="name") {
  x <- suppressWarnings(as.numeric(df[[xcol]]))
  y <- suppressWarnings(as.numeric(df[[ycol]]))
  ok <- is.finite(x) & is.finite(y)
  if (!any(ok)) { plot.new(); text(0.5,0.5,"Kano: no finite data"); return(invisible()) }

  x <- x[ok]; y <- y[ok]
  lab <- if (labelcol %in% names(df)) as.character(df[[labelcol]][ok]) else rep("", length(x))

  # bubble size
  if (!is.null(sizecol) && sizecol %in% names(df)) {
    sz <- suppressWarnings(as.numeric(df[[sizecol]][ok]))
    sz[!is.finite(sz) | sz <= 0] <- 1
  } else {
    sz <- rep(1, length(x))
  }
  cex <- 0.8 + 3 * sqrt(sz / max(sz, na.rm=TRUE))

  # color by group
  if (!is.null(groupcol) && groupcol %in% names(df)) {
    grp <- as.factor(df[[groupcol]][ok])
  } else {
    grp <- factor(rep("1", length(x)))
  }
  colv <- as.integer(grp)

  # center for circles: medians (robust)
  cx <- stats::median(x, na.rm=TRUE)
  cy <- stats::median(y, na.rm=TRUE)

  # radii based on data spread
  rx <- diff(range(x, na.rm=TRUE))
  ry <- diff(range(y, na.rm=TRUE))
  r1 <- 0.18 * sqrt(rx^2 + ry^2)
  r2 <- 0.36 * sqrt(rx^2 + ry^2)
  if (!is.finite(r1) || r1 <= 0) r1 <- 1
  if (!is.finite(r2) || r2 <= 0) r2 <- 2

  op <- par(no.readonly=TRUE); on.exit(par(op), add=TRUE)
  par(mar=c(4,4,2,1)); par(xpd=NA)

  plot(x, y, pch=19, cex=cex, col=colv,
       xlab=if (is.null(xlab)) xcol else xlab,
       ylab=if (is.null(ylab)) ycol else ylab,
       main=title_txt)

  # quadrant lines (Kano-style)
  abline(v=cx, h=cy, lty=2)

  # "wings": gentle diagonals through the center (visual guide)
  abline(a=cy - 0.5*(cx), b=0.5, lty=3)
  abline(a=cy + 0.5*(cx), b=-0.5, lty=3)

  # dual circles
  th <- seq(0, 2*pi, length.out=361)
  lines(cx + r1*cos(th), cy + r1*sin(th), lty=2)
  lines(cx + r2*cos(th), cy + r2*sin(th), lty=2)

  # labels
  text(x, y, labels=lab, pos=4, cex=0.7)

  legend("topleft", legend=levels(grp), col=seq_along(levels(grp)), pch=19, bty="n",
         title=if (!is.null(groupcol)) groupcol else NULL)

  invisible(TRUE)
}

plot_kano_ss_astar <- function(nodes_df, title_txt="Kano: SS vs a*") {
  df <- as.data.frame(nodes_df, stringsAsFactors = FALSE)
  if (!("name" %in% names(df)) && ("id" %in% names(df))) df$name <- df$id
  # normalize columns
  if (!("ss" %in% names(df))) {
    if ("ssi" %in% names(df)) df$ss <- df$ssi
    if ("SSi" %in% names(df)) df$ss <- df$SSi
    if ("sil_width" %in% names(df)) df$ss <- df$sil_width
  }
  if (!("a_star" %in% names(df))) {
    if ("a_star1" %in% names(df)) df$a_star <- df$a_star1
  }

  # Prefer the "real" 2-wings + 2-circles ggplot from kano.R if available
  if (exists("kano_plot_ss_astar", mode="function")) {
    return(kano_plot_ss_astar(df, title_txt = title_txt))
  }
  if (exists("plot_kano_real_xy", mode="function")) {
    # kano.R provides this too
    return(plot_kano_real_xy(df, edges=NULL, xcol="ss", ycol="a_star", sizecol="value",
                             title_txt=title_txt, xlab="SS (ssi)", ylab="a*"))
  }

  # Fallback: simple base version
  .plot_kano_dual_base(df, xcol="ss", ycol="a_star",
                       sizecol="ss", groupcol=if ("carac" %in% names(df)) "carac" else NULL,
                       title_txt=title_txt, xlab="SS (ssi)", ylab="a*")
  invisible(NULL)
}


plot_kano_value_value2 <- function(nodes_df, title_txt="Kano: value vs value2") {
  df <- as.data.frame(nodes_df, stringsAsFactors = FALSE)
  if (!("name" %in% names(df)) && ("id" %in% names(df))) df$name <- df$id

  # Prefer the "real" 2-wings + 2-circles ggplot from kano.R if available
  if (exists("kano_plot", mode="function")) {
    return(kano_plot(df, edges = NULL, title_txt = title_txt, xlab = "value2", ylab = "value"))
  }
  if (exists("plot_kano_real_xy", mode="function")) {
    return(plot_kano_real_xy(df, edges=NULL, xcol="value2", ycol="value", sizecol="value",
                             title_txt=title_txt, xlab="value2", ylab="value"))
  }

  # Fallback: simple base version
  .plot_kano_dual_base(df, xcol="value2", ycol="value",
                       sizecol="value", groupcol=if ("carac" %in% names(df)) "carac" else NULL,
                       title_txt=title_txt, xlab="value2", ylab="value")
  invisible(NULL)
}

# =========================================================
# Lotka helpers (author publication distribution)
# =========================================================
.build_lotka_input_from_author_lists <- function(author_lists_all) {
  if (is.null(author_lists_all) || !length(author_lists_all)) {
    return(data.frame(papers = integer(), authors = integer(), stringsAsFactors = FALSE))
  }

  per_article <- lapply(author_lists_all, function(a) {
    a <- trimws(as.character(a))
    a <- a[!is.na(a) & nzchar(a)]
    unique(a)
  })

  all_authors <- unlist(per_article, use.names = FALSE)
  all_authors <- trimws(as.character(all_authors))
  all_authors <- all_authors[!is.na(all_authors) & nzchar(all_authors)]
  if (!length(all_authors)) {
    return(data.frame(papers = integer(), authors = integer(), stringsAsFactors = FALSE))
  }

  author_pub <- sort(table(all_authors), decreasing = TRUE)
  tb <- as.data.frame(table(as.integer(author_pub)), stringsAsFactors = FALSE)
  names(tb) <- c("papers", "authors")
  tb$papers <- as.integer(as.character(tb$papers))
  tb$authors <- as.integer(tb$authors)
  tb <- tb[order(tb$papers), , drop = FALSE]
  rownames(tb) <- NULL
  tb
}

test_lotka <- function(df,
                       merge_tail   = TRUE,
                       min_expected = 5,
                       make_plot    = TRUE,
                       new_device   = FALSE,
                       device_width = 12,
                       device_height= 5,
                       quiet        = FALSE) {
  stopifnot(is.data.frame(df))
  stopifnot(all(c("papers", "authors") %in% names(df)))

  df <- df[order(df$papers), c("papers", "authors")]
  df <- df[is.finite(df$papers) & is.finite(df$authors), ]
  df <- df[df$papers > 0 & df$authors > 0, ]

  if (nrow(df) < 3) {
    stop("Need at least 3 non-zero categories to test Lotka's law.")
  }

  fit <- lm(log(authors) ~ log(papers), data = df)
  a <- unname(coef(fit)[1])
  b <- unname(coef(fit)[2])
  c_est  <- -b
  A1_est <- exp(a)

  df$expected_raw <- A1_est / (df$papers ^ c_est)
  df$expected <- df$expected_raw / sum(df$expected_raw) * sum(df$authors)
  df$residual <- df$authors - df$expected

  df_test <- df[, c("papers", "authors", "expected")]
  if (merge_tail) {
    while (nrow(df_test) > 3 && any(df_test$expected < min_expected)) {
      i <- nrow(df_test)
      df_test$papers[i - 1]   <- paste0(df_test$papers[i - 1], "+")
      df_test$authors[i - 1]  <- df_test$authors[i - 1] + df_test$authors[i]
      df_test$expected[i - 1] <- df_test$expected[i - 1] + df_test$expected[i]
      df_test <- df_test[-i, ]
    }
  }

  k <- nrow(df_test)
  chisq_stat <- sum((df_test$authors - df_test$expected)^2 / df_test$expected)
  dfree <- k - 3
  p_value <- if (dfree > 0) pchisq(chisq_stat, df = dfree, lower.tail = FALSE) else NA_real_

  if (isTRUE(make_plot)) {
    if (isTRUE(new_device)) {
      try(dev.new(width = device_width, height = device_height), silent = TRUE)
    }
    oldpar <- par(c("mfrow", "mar", "oma", "mgp", "las", "xpd"))
    on.exit(par(oldpar), add = TRUE)
    layout(matrix(c(1, 2), nrow = 1, byrow = TRUE))
    par(mar = c(4.2, 4.2, 2.5, 1.2), oma = c(0, 0, 0.5, 0), mgp = c(2.3, 0.8, 0), las = 1, xpd = NA)

    plot(df$papers, df$authors,
         log  = "xy", pch = 16,
         xlab = "Number of papers (log scale)",
         ylab = "Number of authors (log scale)",
         main = "Lotka log-log plot")
    lines(df$papers, df$expected, lwd = 2, lty = 1)
    legend("topright", legend = c("Observed", "Fitted"), pch = c(16, NA), lty = c(NA, 1), lwd = c(NA, 2), bty = "n")

    ymax <- max(c(df$authors, df$expected), na.rm = TRUE) * 1.08
    plot(df$papers, df$authors,
         type = "b", pch = 16, ylim = c(0, ymax),
         xlab = "Number of papers", ylab = "Number of authors",
         main = "Observed vs expected")
    lines(df$papers, df$expected, type = "b", pch = 1, lty = 2, lwd = 2)
    legend("topright", legend = c("Observed", "Expected"), pch = c(16, 1), lty = c(1, 2), lwd = c(1, 2), bty = "n")
  }

  if (!isTRUE(quiet)) {
    cat("\n=============================\n")
    cat("Lotka's law test results\n")
    cat("=============================\n")
    cat(sprintf("Estimated exponent c  = %.4f\n", c_est))
    cat(sprintf("Estimated A1          = %.4f\n", A1_est))
    cat(sprintf("R-squared             = %.4f\n", summary(fit)$r.squared))
    cat(sprintf("Chi-square statistic  = %.4f\n", chisq_stat))
    cat(sprintf("Degrees of freedom    = %s\n", ifelse(is.na(dfree), "NA", as.character(dfree))))
    cat(sprintf("P-value               = %s\n", ifelse(is.na(p_value), "NA", format(p_value, digits = 6))))
    if (!is.na(p_value)) {
      if (p_value > 0.05) {
        cat("Conclusion            = Data do not significantly differ from Lotka's law.\n")
      } else {
        cat("Conclusion            = Data significantly differ from Lotka's law.\n")
      }
    } else {
      cat("Conclusion            = Too few grouped categories for chi-square test.\n")
    }
  }

  invisible(list(
    model = fit,
    exponent_c = c_est,
    A1_est = A1_est,
    r_squared = summary(fit)$r.squared,
    chi_square = chisq_stat,
    df = dfree,
    p_value = p_value,
    observed_expected_table = df,
    test_table = df_test
  ))
}

.plot_lotka_result <- function(res) {
  if (is.null(res) || is.null(res$observed_expected_table)) {
    plot.new(); text(0.5, 0.5, "No Lotka result available")
    return(invisible(NULL))
  }
  df <- as.data.frame(res$observed_expected_table, stringsAsFactors = FALSE)
  oldpar <- par(c("mfrow", "mar", "oma", "mgp", "las", "xpd"))
  on.exit(par(oldpar), add = TRUE)
  layout(matrix(c(1, 2), nrow = 1, byrow = TRUE))
  par(mar = c(4.2, 4.2, 2.5, 1.2), oma = c(0, 0, 1.4, 0), mgp = c(2.3, 0.8, 0), las = 1, xpd = NA)

  plot(df$papers, df$authors,
       log = "xy", pch = 16,
       xlab = "Number of papers (log scale)",
       ylab = "Number of authors (log scale)",
       main = "Lotka log-log plot")
  lines(df$papers, df$expected, lwd = 2)
  legend("topright", legend = c("Observed", "Fitted"), pch = c(16, NA), lty = c(NA, 1), lwd = c(NA, 2), bty = "n")

  ymax <- max(c(df$authors, df$expected), na.rm = TRUE) * 1.08
  plot(df$papers, df$authors,
       type = "b", pch = 16, ylim = c(0, ymax),
       xlab = "Number of papers", ylab = "Number of authors",
       main = "Observed vs expected")
  lines(df$papers, df$expected, type = "b", pch = 1, lty = 2, lwd = 2)
  legend("topright", legend = c("Observed", "Expected"), pch = c(16, 1), lty = c(1, 2), lwd = c(1, 2), bty = "n")

  chi_txt <- paste0(
    "Chi-square = ", format(round(res$chi_square, 4), nsmall = 4),
    "   |   df = ", ifelse(is.na(res$df), "NA", as.character(res$df)),
    "   |   p = ", ifelse(is.na(res$p_value), "NA", format(res$p_value, digits = 6)),
    "   |   c = ", format(round(res$exponent_c, 4), nsmall = 4)
  )
  mtext(chi_txt, outer = TRUE, cex = 0.95, font = 2)
  invisible(NULL)
}

.lotka_summary_table <- function(res) {
  if (is.null(res)) return(NULL)
  conclusion <- if (is.na(res$p_value)) {
    "Too few grouped categories for chi-square test"
  } else if (res$p_value > 0.05) {
    "Data do not significantly differ from Lotka's law"
  } else {
    "Data significantly differ from Lotka's law"
  }
  data.frame(
    Metric = c("Estimated exponent c", "Estimated A1", "R-squared", "Chi-square statistic", "Degrees of freedom", "P-value", "Conclusion"),
    Value = c(
      sprintf("%.4f", res$exponent_c),
      sprintf("%.4f", res$A1_est),
      sprintf("%.4f", res$r_squared),
      sprintf("%.4f", res$chi_square),
      ifelse(is.na(res$df), "NA", as.character(res$df)),
      ifelse(is.na(res$p_value), "NA", format(res$p_value, digits = 6)),
      conclusion
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

.augment_author_nodes <- function(nodes, edges){
  # Safe Author-node augmentation.
  # This function must never stop the Shiny run. In some Top20 fallbacks,
  # the edge table can have zero valid Leader/Follower matches; merge/aggregate
  # then creates no FA/LA column and `$<-` may crash with:
  # "replacement has 0 rows, data has 20". Use explicit vectors instead.
  tryCatch({
    if (is.null(nodes) || !is.data.frame(nodes) || nrow(nodes) == 0) return(nodes)
    nd <- nodes

    if (!("name" %in% names(nd))) {
      if ("label" %in% names(nd)) nd$name <- as.character(nd$label)
      else if ("id" %in% names(nd)) nd$name <- as.character(nd$id)
      else nd$name <- as.character(seq_len(nrow(nd)))
    }
    nd$name <- trimws(as.character(nd$name))

    # defaults that are safe for zero-edge / fallback Top20 cases.
    # If full FA/LA metrics were already merged into the nodes, preserve them;
    # do not overwrite them using the reduced Top20 one-link display edges.
    .has_existing_fala <- all(c("FA", "LA") %in% names(nd)) &&
      any(is.finite(suppressWarnings(as.numeric(nd$FA)) + suppressWarnings(as.numeric(nd$LA))) &
            (suppressWarnings(as.numeric(nd$FA)) + suppressWarnings(as.numeric(nd$LA))) > 0)
    if (!("FA" %in% names(nd))) nd$FA <- 0
    if (!("LA" %in% names(nd))) nd$LA <- 0
    if (!("single_author" %in% names(nd))) {
      if (all(c("value", "value2") %in% names(nd))) {
        nd$single_author <- pmax(0, suppressWarnings(as.numeric(nd$value) - as.numeric(nd$value2)))
      } else {
        nd$single_author <- 0
      }
    }
    nd$single_author <- suppressWarnings(as.numeric(nd$single_author))
    nd$single_author[!is.finite(nd$single_author)] <- 0
    if (!("n_pubmed" %in% names(nd))) nd$n_pubmed <- NA_real_

    # Normalize edge schema only if usable.
    ed <- edges
    if (!.has_existing_fala && !is.null(ed) && is.data.frame(ed) && nrow(ed) > 0) {
      if (!("Leader" %in% names(ed)) && "leader" %in% names(ed)) ed$Leader <- ed$leader
      if (!("Follower" %in% names(ed)) && "follower" %in% names(ed)) ed$Follower <- ed$follower
      if (!("Leader" %in% names(ed)) && "from" %in% names(ed)) ed$Leader <- ed$from
      if (!("Follower" %in% names(ed)) && "to" %in% names(ed)) ed$Follower <- ed$to
      wcol <- if ("WCD" %in% names(ed)) "WCD" else if ("wcd" %in% names(ed)) "wcd" else if ("weight" %in% names(ed)) "weight" else NULL
      if (is.null(wcol)) { ed$WCD <- 1; wcol <- "WCD" }

      if (all(c("Leader", "Follower") %in% names(ed))) {
        ed$Leader <- trimws(as.character(ed$Leader))
        ed$Follower <- trimws(as.character(ed$Follower))
        ed[[wcol]] <- suppressWarnings(as.numeric(ed[[wcol]]))
        ed <- ed[nzchar(ed$Leader) & nzchar(ed$Follower) & is.finite(ed[[wcol]]) & ed[[wcol]] > 0, , drop = FALSE]

        if (nrow(ed) > 0) {
          fa <- tapply(ed[[wcol]], ed$Leader, sum, na.rm = TRUE)
          la <- tapply(ed[[wcol]], ed$Follower, sum, na.rm = TRUE)
          fa_v <- suppressWarnings(as.numeric(fa[match(nd$name, names(fa))]))
          la_v <- suppressWarnings(as.numeric(la[match(nd$name, names(la))]))
          fa_v[!is.finite(fa_v)] <- 0
          la_v[!is.finite(la_v)] <- 0
          nd$FA <- fa_v
          nd$LA <- la_v
        }
      }
    }

    # value2 = first + last author edge strength; single-author self loop excluded.
    nd$value2 <- suppressWarnings(as.numeric(nd$FA)) + suppressWarnings(as.numeric(nd$LA))
    nd$value2[!is.finite(nd$value2)] <- 0

    # value = first/last author appearance including single-author papers once.
    nd$value <- nd$value2 + suppressWarnings(as.numeric(nd$single_author))
    nd$value[!is.finite(nd$value)] <- 0

    if (!("n_byline" %in% names(nd))) {
      nd$n_byline <- nd$n_pubmed
    } else {
      nd$n_byline <- ifelse(is.na(nd$n_byline) & !is.na(nd$n_pubmed), nd$n_pubmed, nd$n_byline)
    }
    if (!("n_byline_non_single" %in% names(nd))) {
      nd$n_byline_non_single <- ifelse(is.na(nd$n_pubmed), NA_real_, pmax(0, suppressWarnings(as.numeric(nd$n_pubmed)) - suppressWarnings(as.numeric(nd$single_author))))
    }

    nd$value2_strength <- .AAC_INLINE(nd$value2)
    nd
  }, error = function(e) {
    warning("[AUTHOR] .augment_author_nodes skipped safely: ", conditionMessage(e), call. = FALSE)
    nodes
  })
}


# ---- AUTHOR edge guard: exactly one leader edge per follower ----
# Input edges must come from first-author/last-author pairs only.
# This guard re-orients candidate pairs by node value (higher value = leader),
# aggregates duplicate leader-follower pairs, and retains only the strongest
# single incoming edge for every follower. It is used after FLCA-MA-SIL and
# also for deterministic fallback, so Sankey/network displays remain one-link.
.author_enforce_single_edge_per_follower <- function(nodes, edges) {
  if (is.null(edges) || !is.data.frame(edges) || nrow(edges) == 0) {
    return(data.frame(Leader=character(), Follower=character(), WCD=numeric(), stringsAsFactors=FALSE))
  }
  ed <- as.data.frame(edges, stringsAsFactors = FALSE)
  if (!"Leader" %in% names(ed) && "leader" %in% names(ed)) ed$Leader <- ed$leader
  if (!"Follower" %in% names(ed) && "follower" %in% names(ed)) ed$Follower <- ed$follower
  if (!"Leader" %in% names(ed) && "Source" %in% names(ed)) ed$Leader <- ed$Source
  if (!"Follower" %in% names(ed) && "Target" %in% names(ed)) ed$Follower <- ed$Target
  if (!"Leader" %in% names(ed) && "from" %in% names(ed)) ed$Leader <- ed$from
  if (!"Follower" %in% names(ed) && "to" %in% names(ed)) ed$Follower <- ed$to
  if (!"WCD" %in% names(ed)) {
    if ("wcd" %in% names(ed)) ed$WCD <- ed$wcd
    else if ("weight" %in% names(ed)) ed$WCD <- ed$weight
    else if ("value" %in% names(ed)) ed$WCD <- ed$value
    else ed$WCD <- 1
  }
  if (!all(c("Leader","Follower","WCD") %in% names(ed))) {
    return(data.frame(Leader=character(), Follower=character(), WCD=numeric(), stringsAsFactors=FALSE))
  }
  ed$Leader <- trimws(as.character(ed$Leader))
  ed$Follower <- trimws(as.character(ed$Follower))
  ed$WCD <- suppressWarnings(as.numeric(ed$WCD))
  ed <- ed[nzchar(ed$Leader) & nzchar(ed$Follower) & ed$Leader != ed$Follower & is.finite(ed$WCD) & ed$WCD > 0,
           c("Leader","Follower","WCD"), drop=FALSE]
  if (!nrow(ed)) return(data.frame(Leader=character(), Follower=character(), WCD=numeric(), stringsAsFactors=FALSE))

  # Keep only displayed/selected nodes when supplied.
  nd <- nodes
  if (!is.null(nd) && is.data.frame(nd) && nrow(nd) > 0) {
    if (!"name" %in% names(nd)) {
      if ("label" %in% names(nd)) nd$name <- nd$label
      else if ("id" %in% names(nd)) nd$name <- nd$id
      else nd$name <- as.character(seq_len(nrow(nd)))
    }
    nd$name <- trimws(as.character(nd$name))
    nd <- nd[!is.na(nd$name) & nzchar(nd$name), , drop=FALSE]
    nd <- nd[!duplicated(nd$name), , drop=FALSE]
    keep <- nd$name
    ed <- ed[ed$Leader %in% keep & ed$Follower %in% keep, , drop=FALSE]
    if (!nrow(ed)) return(data.frame(Leader=character(), Follower=character(), WCD=numeric(), stringsAsFactors=FALSE))

    # FLCA leader direction: higher node value becomes Leader; ties are stable by name.
    val_col <- if ("value" %in% names(nd)) "value" else if ("value2" %in% names(nd)) "value2" else NULL
    if (!is.null(val_col)) {
      vv <- suppressWarnings(as.numeric(nd[[val_col]])); vv[!is.finite(vv)] <- 0
      names(vv) <- nd$name
      vL <- vv[ed$Leader]; vF <- vv[ed$Follower]
      vL[!is.finite(vL)] <- 0; vF[!is.finite(vF)] <- 0
      swap <- (vF > vL) | (vF == vL & ed$Follower < ed$Leader)
      if (any(swap, na.rm=TRUE)) {
        tmp <- ed$Leader[swap]
        ed$Leader[swap] <- ed$Follower[swap]
        ed$Follower[swap] <- tmp
      }
      ed$.leader_value <- vv[ed$Leader]
      ed$.leader_value[!is.finite(ed$.leader_value)] <- 0
    } else {
      ed$.leader_value <- 0
    }
  } else {
    ed$.leader_value <- 0
  }

  # Aggregate duplicate Leader-Follower pairs first.
  agg <- stats::aggregate(WCD ~ Leader + Follower, data = ed, FUN = function(x) sum(x, na.rm=TRUE))
  lv <- stats::aggregate(.leader_value ~ Leader + Follower, data = ed, FUN = function(x) max(x, na.rm=TRUE))
  agg <- merge(agg, lv, by=c("Leader","Follower"), all.x=TRUE)
  agg$.leader_value[!is.finite(agg$.leader_value)] <- 0

  # Exactly one incoming edge per follower: strongest WCD, then stronger leader, then name.
  agg <- agg[order(agg$Follower, -agg$WCD, -agg$.leader_value, agg$Leader), , drop=FALSE]
  out <- agg[!duplicated(agg$Follower), c("Leader","Follower","WCD"), drop=FALSE]
  out <- out[order(-out$WCD, out$Leader, out$Follower), , drop=FALSE]
  rownames(out) <- NULL
  out
}


# ---- AUTHOR Top20 guard: force at least two clusters for SS calculation ----
# Rationale: a network may visibly have separate components, but silhouette needs
# at least two *cluster labels* in nodes$carac. FLCA/MA can still return one
# dominant cluster when the Top20 authors are all linked to one leader. This
# guard preserves the FLCA Top20 nodes/edges, but if only one cluster label is
# present it derives a second label from Top20 graph components first, then from
# the two strongest leader seeds as a deterministic fallback.
.force_author_min2_clusters_and_ss <- function(nodes20, edges20, data_edges = NULL, cfg = list()) {
  if (is.null(nodes20) || !is.data.frame(nodes20) || nrow(nodes20) == 0) {
    return(list(nodes = nodes20, edges = edges20, changed = FALSE, reason = "no nodes"))
  }

  nd <- nodes20
  if (!"name" %in% names(nd)) {
    if ("label" %in% names(nd)) nd$name <- as.character(nd$label)
    else if ("id" %in% names(nd)) nd$name <- as.character(nd$id)
    else nd$name <- as.character(seq_len(nrow(nd)))
  }
  nd$name <- trimws(as.character(nd$name))
  nd <- nd[!is.na(nd$name) & nzchar(nd$name), , drop = FALSE]
  nd <- nd[!duplicated(nd$name), , drop = FALSE]
  if (!nrow(nd)) return(list(nodes = nd, edges = edges20, changed = FALSE, reason = "empty names"))

  ed <- edges20
  if (is.null(ed) || !is.data.frame(ed)) {
    ed <- data.frame(Leader = character(), Follower = character(), WCD = numeric(), stringsAsFactors = FALSE)
  }
  if (nrow(ed)) {
    if (!"Leader" %in% names(ed) && "leader" %in% names(ed)) ed$Leader <- ed$leader
    if (!"Follower" %in% names(ed) && "follower" %in% names(ed)) ed$Follower <- ed$follower
    if (!"WCD" %in% names(ed)) {
      if ("wcd" %in% names(ed)) ed$WCD <- ed$wcd
      else if ("weight" %in% names(ed)) ed$WCD <- ed$weight
      else if ("value" %in% names(ed)) ed$WCD <- ed$value
      else ed$WCD <- 1
    }
    ed$Leader <- trimws(as.character(ed$Leader))
    ed$Follower <- trimws(as.character(ed$Follower))
    ed$WCD <- suppressWarnings(as.numeric(ed$WCD))
    ed$WCD[!is.finite(ed$WCD) | is.na(ed$WCD)] <- 1
    keep_names <- nd$name
    ed <- ed[ed$Leader %in% keep_names & ed$Follower %in% keep_names & nzchar(ed$Leader) & nzchar(ed$Follower), c("Leader","Follower","WCD"), drop = FALSE]
  } else {
    ed <- data.frame(Leader = character(), Follower = character(), WCD = numeric(), stringsAsFactors = FALSE)
  }

  if (!"carac" %in% names(nd)) nd$carac <- 1L
  nd$carac <- as.character(nd$carac)
  nd$carac[is.na(nd$carac) | !nzchar(nd$carac)] <- "1"
  k0 <- length(unique(nd$carac))
  changed <- FALSE
  reason <- "FLCA labels already had >=2 clusters"

  if (nrow(nd) >= 2 && k0 < 2) {
    # 1) Prefer true graph components of the displayed Top20 network.
    comp_lab <- rep(NA_integer_, nrow(nd))
    if (nrow(ed) > 0 && requireNamespace("igraph", quietly = TRUE)) {
      gg <- tryCatch({
        igraph::graph_from_data_frame(ed[, c("Leader", "Follower"), drop = FALSE], directed = FALSE,
                                      vertices = data.frame(name = nd$name, stringsAsFactors = FALSE))
      }, error = function(e) NULL)
      if (!is.null(gg)) {
        cc <- igraph::components(gg)$membership
        comp_lab <- as.integer(cc[match(nd$name, names(cc))])
      }
    }
    if (length(unique(comp_lab[is.finite(comp_lab)])) >= 2) {
      nd$carac <- as.integer(factor(comp_lab, levels = sort(unique(comp_lab))))
      changed <- TRUE
      reason <- "derived from Top20 graph components"
    } else {
      # 2) If the Top20 graph is one connected component, split by two strongest seeds.
      if (!"value" %in% names(nd)) nd$value <- 1
      if (!"value2" %in% names(nd)) nd$value2 <- nd$value
      val <- suppressWarnings(as.numeric(nd$value)); val[!is.finite(val)] <- 0
      val2 <- suppressWarnings(as.numeric(nd$value2)); val2[!is.finite(val2)] <- 0
      ord <- order(-val, -val2, nd$name)
      leaders <- nd$name[ord[seq_len(min(2L, length(ord)))]]
      nd$carac <- NA_integer_
      nd$carac[match(leaders, nd$name)] <- seq_along(leaders)

      # use all author edges for seed assignment when available, then displayed Top20 edges
      ed_assign <- data_edges
      if (is.null(ed_assign) || !is.data.frame(ed_assign) || !nrow(ed_assign)) ed_assign <- ed
      if (is.data.frame(ed_assign) && nrow(ed_assign)) {
        if (!"Leader" %in% names(ed_assign) && "leader" %in% names(ed_assign)) ed_assign$Leader <- ed_assign$leader
        if (!"Follower" %in% names(ed_assign) && "follower" %in% names(ed_assign)) ed_assign$Follower <- ed_assign$follower
        if (!"WCD" %in% names(ed_assign)) {
          if ("wcd" %in% names(ed_assign)) ed_assign$WCD <- ed_assign$wcd else ed_assign$WCD <- 1
        }
        ed_assign$Leader <- trimws(as.character(ed_assign$Leader))
        ed_assign$Follower <- trimws(as.character(ed_assign$Follower))
        ed_assign$WCD <- suppressWarnings(as.numeric(ed_assign$WCD))
        ed_assign$WCD[!is.finite(ed_assign$WCD) | is.na(ed_assign$WCD)] <- 1
      } else {
        ed_assign <- data.frame(Leader = character(), Follower = character(), WCD = numeric(), stringsAsFactors = FALSE)
      }

      for (ii in which(is.na(nd$carac))) {
        nm <- nd$name[ii]
        hit <- ed_assign[(ed_assign$Leader == nm & ed_assign$Follower %in% leaders) |
                           (ed_assign$Follower == nm & ed_assign$Leader %in% leaders), , drop = FALSE]
        if (nrow(hit) > 0) {
          hit$seed <- ifelse(hit$Leader %in% leaders, hit$Leader, hit$Follower)
          hit <- hit[order(-hit$WCD), , drop = FALSE]
          nd$carac[ii] <- match(hit$seed[1], leaders)
        }
      }
      # any node not linked to a seed is alternated to keep >=2 labels without randomness
      miss <- which(is.na(nd$carac))
      if (length(miss)) nd$carac[miss] <- ((seq_along(miss) - 1L) %% max(1L, length(leaders))) + 1L
      if (length(unique(nd$carac)) < 2 && nrow(nd) >= 2) nd$carac[seq(2, nrow(nd), by = 2)] <- 2L
      changed <- TRUE
      reason <- "derived from two strongest author seeds"
    }
  }

  nd$carac <- as.integer(factor(as.character(nd$carac), levels = unique(as.character(nd$carac))))

  # Recompute silhouette/a* after the cluster-label repair.
  for (cc in c("ssi", "a_i", "b_i", "a_star1")) if (!cc %in% names(nd)) nd[[cc]] <- 0
  if (nrow(ed) > 0 && length(unique(nd$carac)) >= 2 && exists("compute_silhouette_df", mode = "function")) {
    sil <- tryCatch({
      ed2 <- ed
      if (!"follower" %in% names(ed2)) ed2$follower <- ed2$Follower
      compute_silhouette_df(nd, ed2,
                            intra_delta = cfg$intra_delta %||% 2,
                            inter_delta = cfg$inter_delta %||% 5,
                            eps = cfg$eps %||% 1e-9)
    }, error = function(e) {
      message("[AUTHOR min2 cluster] silhouette recompute failed: ", conditionMessage(e))
      NULL
    })
    if (is.data.frame(sil) && nrow(sil) && "name" %in% names(sil)) {
      drop_cols <- intersect(c("ssi", "a_i", "b_i", "a_star1", "SSi", "a_star"), names(nd))
      nd[drop_cols] <- NULL
      nd <- merge(nd, sil, by = "name", all.x = TRUE, sort = FALSE)
      for (cc in c("ssi", "a_i", "b_i", "a_star1")) {
        if (!cc %in% names(nd)) nd[[cc]] <- 0
        nd[[cc]][!is.finite(suppressWarnings(as.numeric(nd[[cc]]))) | is.na(nd[[cc]])] <- 0
      }
      nd$SSi <- nd$ssi
      nd$a_star <- nd$a_star1
    }
  }
  list(nodes = nd, edges = ed, changed = changed, reason = reason)
}

bypass_record <- FALSE

.profile_col_index <- function(df) {
  if (is.null(df) || !is.data.frame(df) || !ncol(df)) return(NA_integer_)
  nms <- trimws(tolower(as.character(names(df))))
  idx <- which(nms %in% c("profile", "profiles"))
  if (length(idx)) return(idx[1])
  NA_integer_
}

.solve_assignment_max <- function(tab) {
  m <- as.matrix(tab)
  storage.mode(m) <- "numeric"
  m[!is.finite(m)] <- 0
  nr <- nrow(m); nc <- ncol(m)
  if (nr < 1L || nc < 1L) return(integer())
  n <- max(nr, nc)
  mm <- matrix(0, nrow = n, ncol = n)
  mm[seq_len(nr), seq_len(nc)] <- m
  if (requireNamespace("clue", quietly = TRUE)) {
    perm <- clue::solve_LSAP(mm, maximum = TRUE)
    return(as.integer(perm)[seq_len(nr)])
  }
  out <- integer(nr)
  used <- rep(FALSE, n)
  for (i in seq_len(nr)) {
    cand <- which(!used)
    if (!length(cand)) cand <- seq_len(n)
    j <- cand[which.max(mm[i, cand])][1]
    if (!length(j) || is.na(j)) j <- cand[1]
    used[j] <- TRUE
    out[i] <- j
  }
  out
}

.kappa_from_table <- function(tab) {
  m <- as.matrix(tab)
  storage.mode(m) <- "numeric"
  m[!is.finite(m)] <- 0
  n <- sum(m)
  if (!is.finite(n) || n <= 0) return(NA_real_)
  po <- sum(diag(m)) / n
  pe <- sum(rowSums(m) * colSums(m)) / (n * n)
  if (!is.finite(pe) || abs(1 - pe) < .Machine$double.eps) return(NA_real_)
  as.numeric((po - pe) / (1 - pe))
}

.build_taaa_profile_agreement <- function(row_df) {
  if (is.null(row_df) || !is.data.frame(row_df) || !nrow(row_df)) return(NULL)
  if (!all(c("Profile", "Cluster") %in% names(row_df))) return(NULL)
  prof <- suppressWarnings(as.integer(as.character(row_df$Profile)))
  cl <- trimws(as.character(row_df$Cluster))
  keep <- is.finite(prof) & prof > 0L & nzchar(cl)
  if (!any(keep)) return(NULL)
  prof <- prof[keep]
  cl <- cl[keep]
  k <- suppressWarnings(as.integer(max(prof, na.rm = TRUE)))
  if (!is.finite(k) || k < 1L) return(NULL)
  clu_levels <- sort(unique(cl))
  if (length(clu_levels) != k) return(NULL)
  tab0 <- table(factor(prof, levels = seq_len(k)), factor(cl, levels = clu_levels))
  perm <- .solve_assignment_max(tab0)
  if (!length(perm)) return(NULL)
  map_df <- data.frame(
    Profile = seq_len(k),
    Cluster = clu_levels[perm[seq_len(min(k, length(perm)))]],
    stringsAsFactors = FALSE
  )
  map_vec <- stats::setNames(map_df$Profile, map_df$Cluster)
  mapped <- unname(map_vec[cl])
  tab1 <- table(factor(prof, levels = seq_len(k)), factor(mapped, levels = seq_len(k)))
  conf_df <- as.data.frame.matrix(tab1, stringsAsFactors = FALSE)
  conf_df <- cbind(Profile = rownames(conf_df), conf_df, stringsAsFactors = FALSE)
  rownames(conf_df) <- NULL
  kap <- .kappa_from_table(tab1)
  n <- sum(tab1)
  po <- if (n > 0) sum(diag(tab1)) / n else NA_real_
  pe <- if (n > 0) sum(rowSums(tab1) * colSums(tab1)) / (n * n) else NA_real_
  kappa_df <- data.frame(K = k, N = as.integer(n), Agreement = round(po, 4), Expected = round(pe, 4), Kappa = round(kap, 4), stringsAsFactors = FALSE, check.names = FALSE)
  list(mapping_table = map_df, kappa_table = kappa_df, confusion_table = conf_df)
}

.build_taaa_from_term_list <- function(term_list, nodes_df, profile_vec = NULL) {
  if (is.null(term_list) || !length(term_list)) return(NULL)
  if (is.null(nodes_df) || !is.data.frame(nodes_df) || !nrow(nodes_df)) return(NULL)
  if (!all(c("name", "carac") %in% names(nodes_df))) return(NULL)

  nd <- as.data.frame(nodes_df, stringsAsFactors = FALSE)
  nd$name <- trimws(as.character(nd$name))
  nd$carac <- trimws(as.character(nd$carac))
  val_col <- if ("value" %in% names(nd)) "value" else names(nd)[1]
  nd$value <- suppressWarnings(as.numeric(nd[[val_col]]))
  nd$value[!is.finite(nd$value)] <- 0
  nd <- nd[nzchar(nd$name) & nzchar(nd$carac), , drop = FALSE]
  if (!nrow(nd)) return(NULL)

  rep_df <- do.call(rbind, lapply(split(nd, nd$carac), function(z) {
    z <- z[order(-z$value, z$name), , drop = FALSE]
    data.frame(carac = as.character(z$carac[1]), theme_name = as.character(z$name[1]), cluster_n = nrow(z), leader_value = as.numeric(z$value[1]), stringsAsFactors = FALSE)
  }))
  rep_df$cluster_n[!is.finite(rep_df$cluster_n)] <- 0
  rep_df$leader_value[!is.finite(rep_df$leader_value)] <- 0

  term_to_cluster <- stats::setNames(as.character(nd$carac), nd$name)
  if (is.null(profile_vec)) profile_vec <- rep(NA_character_, length(term_list))

  out <- vector("list", length(term_list))
  kk <- 0L
  for (i in seq_along(term_list)) {
    vv <- unique(trimws(as.character(term_list[[i]] %||% character(0))))
    vv <- vv[nzchar(vv)]
    if (!length(vv)) next
    cl <- unname(term_to_cluster[vv])
    cl <- cl[!is.na(cl) & nzchar(cl)]
    if (!length(cl)) next
    tb <- as.data.frame(table(cl), stringsAsFactors = FALSE)
    names(tb) <- c("carac", "count")
    tb$count <- suppressWarnings(as.numeric(tb$count))
    tb$count[!is.finite(tb$count)] <- 0
    tb$theme_name <- rep_df$theme_name[match(tb$carac, rep_df$carac)]
    tb$cluster_n  <- rep_df$cluster_n[match(tb$carac, rep_df$carac)]
    tb$leader_value <- rep_df$leader_value[match(tb$carac, rep_df$carac)]
    tb$theme_name[is.na(tb$theme_name)] <- tb$carac[is.na(tb$theme_name)]
    tb$cluster_n[!is.finite(tb$cluster_n)] <- 0
    tb$leader_value[!is.finite(tb$leader_value)] <- 0
    carac_num <- suppressWarnings(as.numeric(tb$carac))
    carac_num[!is.finite(carac_num)] <- Inf
    tb <- tb[order(-tb$count, carac_num, -tb$leader_value, tb$theme_name), , drop = FALSE]
    chosen <- as.character(tb$carac[1])
    theme <- as.character(tb$theme_name[1])
    theme_count_str <- paste0(tb$theme_name, "(", tb$count, ")", collapse = ", ")
    kk <- kk + 1L
    out[[kk]] <- data.frame(Row = i, Profile = as.character(profile_vec[i] %||% ""), 主題名稱 = theme %||% "", 主題詞頻 = theme_count_str, Cluster = chosen %||% "", stringsAsFactors = FALSE, check.names = FALSE)
  }
  if (kk == 0L) return(NULL)
  row_df <- do.call(rbind, out[seq_len(kk)])
  rownames(row_df) <- NULL
  freq_df <- as.data.frame(table(row_df$主題名稱), stringsAsFactors = FALSE)
  names(freq_df) <- c("主題名稱", "Frequency")
  freq_df$Frequency <- suppressWarnings(as.numeric(freq_df$Frequency))
  freq_df$Frequency[!is.finite(freq_df$Frequency)] <- 0
  freq_df$Cluster <- rep_df$carac[match(freq_df$主題名稱, rep_df$theme_name)]
  freq_df$群大小 <- rep_df$cluster_n[match(freq_df$主題名稱, rep_df$theme_name)]
  freq_df$群首值 <- rep_df$leader_value[match(freq_df$主題名稱, rep_df$theme_name)]
  freq_df$群大小[!is.finite(freq_df$群大小)] <- 0
  freq_df$群首值[!is.finite(freq_df$群首值)] <- 0
  freq_df <- freq_df[order(-freq_df$Frequency, -freq_df$群大小, -freq_df$群首值, freq_df$主題名稱), , drop = FALSE]
  rownames(freq_df) <- NULL
  agree <- .build_taaa_profile_agreement(row_df)
  list(row_table = row_df, freq_table = freq_df, profile_map_table = agree$mapping_table %||% NULL, kappa_table = agree$kappa_table %||% NULL, confusion_table = agree$confusion_table %||% NULL)
}



# ============================================================
# AMA/PubMed reference-list journal extractor (Top 20 + pie)
# Added for uploaded PubMed/AMA .txt files
# ============================================================
.read_text_any_encoding <- function(path){
  if (is.null(path) || !nzchar(path) || !file.exists(path)) return("")
  encs <- c("UTF-8", "UTF-8-BOM", "CP950", "Big5", "latin1")
  for (enc in encs) {
    txt <- tryCatch(readr::read_file(path, locale = readr::locale(encoding = enc)),
                    error = function(e) NULL)
    if (!is.null(txt) && nchar(txt) > 0) return(enc2utf8(txt))
  }
  paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
}

.split_ama_refs <- function(txt){
  txt <- paste(txt, collapse = "\n")
  txt <- gsub("\r\n|\r", "\n", txt)
  txt <- trimws(txt)
  if (!nzchar(txt)) return(character(0))

  # Numbered AMA/PubMed summary style: 1: ... or 1. ...
  refs <- unlist(strsplit(txt, "\n(?=\\s*\\d+\\s*[:\\.])", perl = TRUE))
  refs <- trimws(gsub("\\s+", " ", refs))
  refs <- refs[nzchar(refs)]

  # Fallback: blank-line records, useful for PubMed "Summary/Text" exports
  if (length(refs) <= 1) {
    refs2 <- unlist(strsplit(txt, "\n\\s*\n+", perl = TRUE))
    refs2 <- trimws(gsub("\\s+", " ", refs2))
    refs2 <- refs2[nzchar(refs2)]
    if (length(refs2) > length(refs)) refs <- refs2
  }
  refs
}


.clean_journal_candidate <- function(x){
  x <- trimws(as.character(x))
  x <- gsub("\\s+", " ", x)
  x <- gsub("^[[:punct:] ]+|[[:punct:] ]+$", "", x)
  x <- trimws(x)

  bad <- is.na(x) | x == "" |
    grepl("^\\d{4}$", x) |
    grepl("^\\d{4}[;:. ]", x) |
    grepl("^(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)\\b", x, ignore.case = TRUE) |
    grepl("^(doi|pmid|pmcid|abstract|free article)\\b", x, ignore.case = TRUE) |
    grepl("\\b(accessed|retrieved|available from|available at)\\b", x, ignore.case = TRUE, perl = TRUE) |
    grepl("(取自|檢索|检索|瀏覽|浏览|取得|查詢)", x, ignore.case = TRUE, perl = TRUE) |
    grepl("https?://|www[.]", x, ignore.case = TRUE, perl = TRUE) |
    grepl("^[0-9]+$", x) |
    nchar(x) < 3
  x[bad] <- NA_character_
  x
}

.extract_medline_field <- function(txt, tag){
  txt <- paste(txt, collapse = "\n")
  txt <- gsub("\r\n|\r", "\n", txt)
  # MEDLINE continuation lines begin with six spaces; join them to the field.
  pat <- paste0("(?ms)^", tag, "\\s*-\\s*(.*?)(?=\\n[A-Z0-9]{2,4}\\s*-|\\nPMID-\\s*|\\z)")
  m <- gregexpr(pat, txt, perl = TRUE)
  z <- regmatches(txt, m)[[1]]
  if (!length(z)) return(character(0))
  z <- sub(paste0("(?s)^", tag, "\\s*-\\s*"), "", z, perl = TRUE)
  z <- gsub("\n\\s+", " ", z, perl = TRUE)
  z <- trimws(gsub("\\s+", " ", z))
  z[nzchar(z)]
}

.extract_journals_from_medline_txt <- function(txt){
  # Correct source for PubMed MEDLINE exports:
  # JT = full journal title; TA = ISO/abbreviated journal title; SO = source fallback.
  jt <- .extract_medline_field(txt, "JT")
  ta <- .extract_medline_field(txt, "TA")
  so <- .extract_medline_field(txt, "SO")

  n <- max(length(jt), length(ta), length(so), 0)
  if (n == 0) return(data.frame())

  pad <- function(x, n) c(x, rep("", max(0, n - length(x))))[seq_len(n)]
  jt <- pad(jt, n); ta <- pad(ta, n); so <- pad(so, n)

  journal <- ifelse(nzchar(jt), jt, ifelse(nzchar(ta), ta, ""))
  # SO usually starts like: Journal title. 2024 Jan;...
  so_j <- sub("\\.\\s*\\d{4}.*$", "", so)
  journal <- ifelse(!nzchar(journal) & nzchar(so_j), so_j, journal)
  journal <- .clean_journal_candidate(journal)

  out <- data.frame(
    reference_no = seq_len(n),
    reference = so,
    journal = journal,
    stringsAsFactors = FALSE
  )
  out[!is.na(out$journal) & out$journal != "", , drop = FALSE]
}

.extract_journal_from_pubmed_summary_ref <- function(ref){
  ref0 <- gsub("\r\n|\r", "\n", ref)
  lines <- trimws(unlist(strsplit(ref0, "\n+", perl = TRUE)))
  lines <- lines[nzchar(lines)]
  lines <- gsub("^\\s*\\d+\\s*[:.]\\s*", "", lines)

  # PubMed Summary/Text export often has a line like:
  # Medicine (Baltimore). 2022 Jan 14;101(2):e285xx.
  cand <- character(0)
  for (ln in lines) {
    z <- stringr::str_match(ln, "^(.+?)\\.\\s*(?:Epub\\s+)?\\d{4}(?:\\s|;|:|\\.|$)")[, 2]
    if (!is.na(z) && nzchar(z)) cand <- c(cand, z)
  }
  cand <- .clean_journal_candidate(cand)
  cand <- cand[!is.na(cand)]
  if (length(cand)) return(cand[1])

  # Single-line AMA reference: Authors. Title. Journal. 2026;...
  one <- trimws(gsub("\\s+", " ", ref0))
  one <- sub("^\\s*\\d+\\s*[:.]\\s*", "", one)
  parts <- trimws(unlist(strsplit(one, "\\.\\s+", perl = TRUE)))
  parts <- parts[nzchar(parts)]
  yr_pos <- which(grepl("^\\d{4}(\\s|;|:|\\.|$)", parts))
  if (length(yr_pos) > 0 && yr_pos[1] > 1) {
    j <- .clean_journal_candidate(parts[yr_pos[1] - 1])
    if (!is.na(j) && nzchar(j)) return(j)
  }

  # Last fallback for compact AMA: . Journal. 2026;
  j <- stringr::str_match(one, "\\.\\s*([^.;:]+?)\\.\\s*\\d{4}(?:\\s|;|:|\\.|$)")[, 2]
  j <- .clean_journal_candidate(j)
  ifelse(is.na(j), "", j)
}

.split_ama_refs <- function(txt){
  txt <- paste(txt, collapse = "\n")
  txt <- gsub("\r\n|\r", "\n", txt)
  txt <- trimws(txt)
  if (!nzchar(txt)) return(character(0))

  refs <- unlist(strsplit(txt, "\n(?=\\s*\\d+\\s*[:.])", perl = TRUE))
  refs <- trimws(refs)
  refs <- refs[nzchar(refs)]

  if (length(refs) <= 1) {
    refs2 <- unlist(strsplit(txt, "\n\\s*\n+", perl = TRUE))
    refs2 <- trimws(refs2)
    refs2 <- refs2[nzchar(refs2)]
    if (length(refs2) > length(refs)) refs <- refs2
  }
  refs
}

.extract_ama_pubmed_journal_results <- function(txt, top_n = 20){
  top_n <- max(1, suppressWarnings(as.integer(top_n))[1] %||% 20L)

  # Always try MEDLINE JT/TA/SO first. This avoids the old error where years
  # such as 2022/2023/2024 were counted as journal names.
  medline_df <- if (grepl("(?m)^(JT|TA|SO)\\s*-\\s*", txt, perl = TRUE)) {
    .extract_journals_from_medline_txt(txt)
  } else data.frame()

  if (is.data.frame(medline_df) && nrow(medline_df) > 0) {
    journal_df <- medline_df
  } else {
    refs <- .split_ama_refs(txt)
    journal_df <- data.frame(
      reference_no = seq_along(refs),
      reference = refs,
      journal = vapply(refs, .extract_journal_from_pubmed_summary_ref, character(1)),
      stringsAsFactors = FALSE
    )
    journal_df$journal <- .clean_journal_candidate(journal_df$journal)
    journal_df <- journal_df[!is.na(journal_df$journal) & journal_df$journal != "", , drop = FALSE]
  }

  if (!nrow(journal_df)) {
    return(list(journals = journal_df,
                top = data.frame(journal = character(), frequency = integer()),
                total_refs = 0L))
  }

  top_df <- journal_df |>
    dplyr::count(journal, name = "frequency", sort = TRUE) |>
    dplyr::slice_head(n = top_n)

  list(journals = journal_df, top = top_df, total_refs = nrow(journal_df))
}

.plot_ama_journal_pie <- function(top_df){
  if (is.null(top_df) || !is.data.frame(top_df) || !nrow(top_df)) {
    plot.new(); text(0.5, 0.5, "No journal data available."); return(invisible(NULL))
  }
  df <- top_df |>
    dplyr::mutate(
      journal_short = ifelse(nchar(journal) > 30, paste0(substr(journal, 1, 30), "..."), journal),
      pct = frequency / sum(frequency),
      label = paste0(journal_short, "\n", frequency, " (", round(pct * 100, 1), "%)")
    )

  ggplot2::ggplot(df, ggplot2::aes(x = "", y = frequency, fill = journal_short)) +
    ggplot2::geom_col(width = 1, color = "white") +
    ggplot2::coord_polar(theta = "y") +
    ggplot2::geom_text(ggplot2::aes(label = label),
                       position = ggplot2::position_stack(vjust = 0.5),
                       size = 3) +
    ggplot2::labs(title = "Top Journal Frequency Pie Plot", fill = "Journal") +
    ggplot2::theme_void(base_size = 13) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
                   legend.position = "right")
}



# ============================================================
# AMA/PubMed author + journal extraction for FLCA-MA-SIL
# ============================================================
.ref_is_doi_fragment <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- trimws(tolower(x))
  out <- x == "" |
    grepl("^10\\.", x) |
    grepl("^\\d{3,}/", x) |
    grepl("journal\\.", x) |
    grepl("/", x)
  out[is.na(out)] <- TRUE
  out
}

.ref_is_access_fragment <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  z <- trimws(tolower(enc2utf8(x)))
  out <- z == "" |
    grepl("\\b(accessed|retrieved|cited|available from|available at|last accessed|access date)\\b", z, perl = TRUE) |
    grepl("\\b(accessed|retrieved)\\s+(jan|feb|mar|apr|may|jun|jul|aug|sep|sept|oct|nov|dec|january|february|march|april|june|july|august|september|october|november|december|[0-9])", z, perl = TRUE) |
    grepl("\\bpublished\\s+((online|ahead of print)\\s+)?(19|20)[0-9]{2}\\b", z, perl = TRUE) |
    grepl("(取自|檢索|检索|瀏覽|浏览|取得|查詢|查询)", z, perl = TRUE) |
    grepl("https?://|www[.]", z, perl = TRUE)
  out[is.na(out)] <- TRUE
  out
}


# Reject article-title fragments that are often misread as journal names
# when pasted AMA/PubMed references contain a title followed by a year.
.ref_title_like_journal_fragment <- function(x) {
  x0 <- enc2utf8(as.character(x %||% ""))
  x0[is.na(x0)] <- ""
  z <- tolower(x0)
  z <- gsub("[.]+", " ", z, perl = TRUE)
  z <- gsub("[[:space:]]+", " ", z, perl = TRUE)
  z <- trimws(z)

  # Title fragments frequently end with a preposition immediately before the year,
  # e.g., "A bibliometric analysis from 2024" or "impact of 2024".
  # Real journal abbreviations normally end with a capitalized journal token or a
  # balanced parenthetical qualifier, e.g., "Medicine (Baltimore)".
  ends_with_prep <- grepl("\\b(from|of|for|with|without|among|between|in|on|to|by|about|against|during|after|before|under|over|into|through|via)$", z, perl = TRUE)

  # Lowercase running-sentence fragments are rejected unless they look like an
  # accepted journal title pattern. This keeps title text out while still allowing
  # official names such as "Journal of Advanced Nursing".
  words0 <- unlist(strsplit(trimws(x0), "[[:space:]]+", perl = TRUE), use.names = FALSE)
  words0 <- words0[nzchar(words0)]
  first_lower <- length(words0) > 0 && grepl("^[a-z]", words0[1], perl = TRUE)
  lower_words <- sum(grepl("^[a-z]", words0, perl = TRUE))
  allowed_lower <- c("of", "and", "the", "in", "for")
  bad_lower_sentence <- first_lower || any(grepl("^[a-z]", words0, perl = TRUE) & !(tolower(words0) %in% allowed_lower))

  out <- !nzchar(z) |
    ends_with_prep |
    bad_lower_sentence |
    grepl("\\bbibliometric analysis\\b", z, perl = TRUE) |
    grepl("\\b(systematic|scoping|narrative|umbrella) review\\b", z, perl = TRUE) |
    grepl("\\bmeta[ -]?analysis\\b", z, perl = TRUE) |
    grepl("\\bcase report\\b|\\bresearch protocol\\b", z, perl = TRUE) |
    grepl("\\b(from|between)\\s+(19|20)[0-9]{2}\\b.*\\b(to|and|through)\\b.*\\b(19|20)[0-9]{2}\\b", z, perl = TRUE) |
    grepl("^(a|an|the)\\s+.+\\b(analysis|study|review|report|survey|assessment|evaluation)\\b", z, perl = TRUE) |
    grepl("\\b(title|abstract|background|objective|methods?|results?|conclusions?)\\s*[:：]", z, perl = TRUE)
  # Unbalanced punctuation such as "A bibliometric analysis (" is not a journal title.
  out <- out | (grepl("\\(", x0, perl = TRUE) != grepl("\\)", x0, perl = TRUE)) |
    (grepl("\\[", x0, perl = TRUE) != grepl("\\]", x0, perl = TRUE)) |
    grepl("[\\(\\[]\\s*$", x0, perl = TRUE)
  out[is.na(out)] <- TRUE
  out
}

# Final journal-display normalization: remove period symbols inside journal names
# (e.g., "J. Biol. Chem." -> "J Biol Chem") while keeping meaningful balanced
# parenthetical qualifiers such as "Medicine (Baltimore)".
.ref_clean_journal_symbols <- function(x) {
  x <- enc2utf8(as.character(x %||% ""))
  x[is.na(x)] <- ""
  x <- gsub("[./／]+", " ", x, perl = TRUE)
  x <- gsub("[,;:；：，]+", " ", x, perl = TRUE)
  x <- gsub("[[:space:]]+", " ", x, perl = TRUE)
  trimws(x)
}

.ref_split_refs_smart <- function(x) {
  if (length(x) > 1L) x <- paste(x, collapse = "\n")
  x <- enc2utf8(x)
  x <- gsub("\r\n|\r", "\n", x)
  x <- trimws(x)
  if (!nzchar(x)) return(character(0))

  if (grepl("(?m)^PMID\\s*-", x, perl = TRUE)) {
    refs <- unlist(strsplit(x, "(?m)(?=^PMID\\s*-)", perl = TRUE), use.names = FALSE)
    refs <- trimws(refs)
    refs <- refs[nzchar(refs)]
    if (length(refs)) return(refs)
  }

  chunks <- unlist(strsplit(x, "\\n\\s*\\n", perl = TRUE), use.names = FALSE)
  chunks <- trimws(chunks)
  chunks <- chunks[nzchar(chunks)]
  if (length(chunks) > 1L) {
    chunks <- gsub("\\s*\\n\\s*", " ", chunks, perl = TRUE)
    return(trimws(gsub("\\s+", " ", chunks)))
  }

  extract_numbered <- function(txt) {
    m <- gregexpr(
      "(?ms)^\\s*(?:\\[?\\d+\\]?\\s*[:;.\\)]|\\d+\\s+)\\s*.*?(?=(?:\\n\\s*(?:\\[?\\d+\\]?\\s*[:;.\\)]|\\d+\\s+))|\\Z)",
      txt, perl = TRUE
    )
    if (m[[1]][1] == -1) return(character(0))
    refs <- regmatches(txt, m)[[1]]
    refs <- gsub("\\s*\\n\\s*", " ", refs, perl = TRUE)
    trimws(gsub("\\s+", " ", refs[nzchar(trimws(refs))]))
  }

  refs <- extract_numbered(x)
  if (length(refs) > 1L) return(refs)

  x_fix <- gsub("(?<!\\n)\\s(?=(?:\\[?\\d+\\]?\\s*[:;.\\)]))", "\n", x, perl = TRUE)
  refs2 <- extract_numbered(x_fix)
  if (length(refs2) > 1L) return(refs2)

  lines <- trimws(unlist(strsplit(x, "\n", fixed = TRUE), use.names = FALSE))
  lines <- lines[nzchar(lines)]
  trimws(gsub("\\s+", " ", lines))
}

.ref_medline_field <- function(block, tag) {
  pat <- paste0("(?m)^", tag, "\\s*-\\s*(.*)$")
  m <- gregexpr(pat, block, perl = TRUE)
  if (m[[1]][1] == -1) return(character(0))
  vals <- regmatches(block, m)[[1]]
  vals <- sub(paste0("^", tag, "\\s*-\\s*"), "", vals, perl = TRUE)
  vals <- trimws(gsub("\\s+", " ", vals))
  vals[nzchar(vals)]
}

.ref_parse_medline <- function(block) {
  au <- .ref_medline_field(block, "AU")
  if (!length(au)) au <- .ref_medline_field(block, "FAU")
  jt <- .ref_medline_field(block, "JT")
  ta <- .ref_medline_field(block, "TA")
  journal <- if (length(jt)) jt[1] else if (length(ta)) ta[1] else NA_character_
  list(journal = journal, authors = au)
}

.ref_strip_tails <- function(s) {
  sub("(published[[:space:]]+((online|ahead[[:space:]]+of[[:space:]]+print)[[:space:]]+)?(19|20)[0-9]{2}|doi[[:space:]]*:|pmid[[:space:]]*:|pmcid[[:space:]]*:|epub|accessed|retrieved|available[[:space:]]+(from|at)[[:space:]]*:|https?://|www[.]|取自|檢索|检索|瀏覽|浏览|取得|查詢|查询).*$",
      "", s, perl = TRUE, ignore.case = TRUE)
}

.ref_norm_journal <- function(j) {
  j0 <- trimws(gsub("\\s+", " ", as.character(j %||% "")))
  j0 <- gsub("\\[Internet\\]", "", j0, ignore.case = TRUE)
  j0 <- sub("(published[[:space:]]+((online|ahead[[:space:]]+of[[:space:]]+print)[[:space:]]+)?(19|20)[0-9]{2}|accessed|retrieved|available[[:space:]]+(from|at)|取自|檢索|检索|瀏覽|浏览|取得|查詢|查询).*$", "", j0, perl = TRUE, ignore.case = TRUE)
  j0 <- trimws(sub("[\\.,;:/／；：，]+\\s*$", "", j0))
  if (grepl("[:;]", j0)) {
    parts <- trimws(unlist(strsplit(j0, "[:;]", perl = TRUE)))
    parts <- parts[nzchar(parts) & !.ref_is_access_fragment(parts)]
    if (length(parts)) j0 <- parts[length(parts)]
  }
  j <- .ref_clean_journal_symbols(j0)
  if (!nzchar(j) || .ref_is_doi_fragment(j0) || .ref_is_access_fragment(j0) || .ref_title_like_journal_fragment(j0)) NA_character_ else j
}


.ref_journal_like <- function(x) {
  if (is.null(x) || length(x) == 0 || is.na(x[1])) return(FALSE)
  x <- trimws(as.character(x[1]))
  if (!nzchar(x) || .ref_is_doi_fragment(x) || .ref_is_access_fragment(x) || .ref_title_like_journal_fragment(x)) return(FALSE)
  # Authors usually contain commas; journal titles normally do not.
  if (grepl(",", x)) return(FALSE)
  if (grepl("^[a-z0-9\\-]+$", x)) return(FALSE)
  if (grepl("\\b(from|of|for|with|without|among|between|in|on|to)$", x, ignore.case = TRUE, perl = TRUE)) return(FALSE)
  # Allow English journal names, short abbreviations such as JAMA, and Chinese journal names.
  has_letter_or_han <- grepl("[A-Za-z]", x) || grepl("[\u4e00-\u9fff]", x)
  if (!has_letter_or_han) return(FALSE)
  grepl("\\s", x) || grepl("\\(", x) || nchar(x) >= 4 || grepl("^[A-Z]{2,10}$", x) || grepl("[\u4e00-\u9fff]", x)
}



.ref_extract_authors_ama <- function(ref) {
  clean <- trimws(gsub("\\s+", " ", gsub("^\\s*(?:\\[\\d+\\]|\\(?\\d+\\)?|\\d+)\\s*[:;.\\)]?\\s*", "", ref, perl = TRUE)))
  # Author segment ends at English period or Chinese full stop.
  authors_str <- sub("^([^\\.。]+?)[\\.。].*$", "\\1", clean, perl = TRUE)
  authors_str <- gsub("\\band\\b", ",", authors_str, ignore.case = TRUE)
  # Split English AMA authors and Chinese author delimiters.
  parts <- trimws(unlist(strsplit(authors_str, "\\s*,\\s*|、|；|;", perl = TRUE)))
  parts <- parts[nzchar(parts)]
  parts <- sub("[\\.。]+$", "", parts)
  parts <- parts[!grepl("^et\\s*al\\.?$", parts, ignore.case = TRUE)]
  # Keep English author/group names and Chinese personal/organization names.
  keep <- grepl("[A-Za-z]", parts) | grepl("[\u4e00-\u9fff]", parts)
  unique(parts[keep])
}


.ref_extract_journal_ama <- function(ref) {
  s <- trimws(gsub("\\s+", " ", gsub("^\\s*(?:\\[\\d+\\]|\\(?\\d+\\)?|\\d+)\\s*[:;.\\)]?\\s*", "", ref, perl = TRUE)))
  s <- .ref_strip_tails(s)
  s <- gsub("[；：]", ";", s, perl = TRUE)
  s <- gsub("([!?])+(\\s*)", ".\\2", s, perl = TRUE)

  # Strict AMA/Vancouver rule:
  #   Authors. Title. Journal Name. 2024;...
  #   Authors. Title. Journal Name 2024;...
  #   Authors. Title. Journal Name. 2024 Nov 3;...
  #   Authors. Title. Journal Name 2024 Nov 3;...
  #   Authors. Title. Journal Name. 2024/... or Journal Name 2024/...
  # Also accept a space/end after year/date (Name. 2024 or Name 2024), but reject
  # title fragments ending in prepositions such as from/of before the year.
  .try_journal_candidate <- function(cand) {
    cand0 <- .ref_norm_spaces2(cand)
    # Reject title fragments ending with a preposition immediately before the
    # publication year/date, e.g. "A bibliometric analysis from 2025" or
    # "impact of 2025". Do not trim the preposition and retry, because that
    # could turn a title fragment into a false journal name.
    if (grepl("\\b(from|of|for|with|without|among|between|in|on|to|by|about|against|during|after|before|under|over|into|through|via)\\s*$", cand0, ignore.case = TRUE, perl = TRUE)) {
      return(NA_character_)
    }
    cand <- .ref_norm_journal(cand0)
    if (.ref_journal_like(cand)) return(cand)
    NA_character_
  }

  token <- "(?:[A-Z(][A-Za-z&'()\\-/]*\\.?|[A-Z]{1,10}\\.?|[\u4e00-\u9fff]{2,})(?:[[:space:]]+|$)"
  month_pat <- "(?:Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Sept|Oct|Nov|Dec)[A-Za-z.]*"
  # Valid date/year markers. Examples: 2025; 2025/ 2025 Nov 3; 2025 Oct 24; 2025 Nov.
  # The trailing look-ahead accepts semicolon/slash/comma/colon, whitespace, or end of string.
  year_mark <- paste0("(?:19|20)[0-9]{2}(?:[[:space:]]+", month_pat, ")?(?:[[:space:]]+[0-3]?[0-9])?[[:space:]]*(?=$|[;,:/／；：，]|[[:space:]])")

  pat_dot_year <- paste0("(?:^|[.!?][[:space:]]+)((?:", token, "){1,12})\\.[[:space:]]*", year_mark)
  m <- gregexpr(pat_dot_year, s, perl = TRUE, ignore.case = FALSE)
  hit <- regmatches(s, m)[[1]]
  if (length(hit) && hit[1] != "") {
    for (mm in rev(hit)) {
      cand <- sub("^[.!?][[:space:]]+", "", mm, perl = TRUE)
      cand <- sub(paste0("\\.[[:space:]]*", year_mark, ".*$"), "", cand, perl = TRUE)
      j <- .try_journal_candidate(cand)
      if (!is.na(j) && nzchar(j)) return(j)
    }
  }

  pat_space_year <- paste0("(?:^|[.!?][[:space:]]+)((?:", token, "){1,12})", year_mark)
  m <- gregexpr(pat_space_year, s, perl = TRUE, ignore.case = FALSE)
  hit <- regmatches(s, m)[[1]]
  if (length(hit) && hit[1] != "") {
    for (mm in rev(hit)) {
      cand <- sub("^[.!?][[:space:]]+", "", mm, perl = TRUE)
      cand <- sub(paste0(year_mark, ".*$"), "", cand, perl = TRUE)
      j <- .try_journal_candidate(cand)
      if (!is.na(j) && nzchar(j)) return(j)
    }
  }

  # Chinese or mixed references: journal segment often appears immediately before a year,
  # e.g. 醫療資訊雜誌 2025；34(1)：42-52。
  segs_cn <- trimws(unlist(strsplit(s, "[。]\\s*", perl = TRUE)))
  segs_cn <- segs_cn[nzchar(segs_cn)]
  if (length(segs_cn)) {
    for (seg in rev(segs_cn)) {
      if (grepl(paste0(year_mark), seg, perl = TRUE, ignore.case = TRUE)) {
        cand <- sub(paste0("\\s*", year_mark, ".*$"), "", seg, perl = TRUE)
        cand <- trimws(gsub("[；;:：,，/／]+$", "", cand, perl = TRUE))
        j <- .try_journal_candidate(cand)
        if (!is.na(j) && nzchar(j)) return(j)
      }
    }
  }

  # Conservative fallback: inspect only the text immediately before a valid
  # publication year/date followed by a volume/date separator or space/end.
  valid_year_pat <- year_mark
  m <- regexpr(valid_year_pat, s, perl = TRUE, ignore.case = TRUE)
  if (m[1] > 0) {
    pre <- substr(s, 1, m[1] - 1)
    pre <- trimws(sub("[\\.[:space:]]+$", "", pre, perl = TRUE))
    # Try suffixes after sentence boundaries while preserving journal abbreviations.
    pieces <- trimws(unlist(strsplit(pre, "\\.[[:space:]]+", perl = TRUE)))
    pieces <- pieces[nzchar(pieces)]
    if (length(pieces)) {
      for (k in seq_along(pieces)) {
        cand <- paste(pieces[k:length(pieces)], collapse = " ")
        j <- .try_journal_candidate(cand)
        if (!is.na(j) && nzchar(j)) return(j)
      }
    }
  }
  NA_character_
}



.ref_parse_one <- function(ref) {
  if (grepl("(?m)^PMID\\s*-", ref, perl = TRUE)) return(.ref_parse_medline(ref))
  list(journal = .ref_extract_journal_ama(ref), authors = .ref_extract_authors_ama(ref))
}

.ref_to_wide <- function(refs, fallback_author_to_journal = TRUE) {
  # IMPORTANT: preserve line breaks in MEDLINE/NBIB blocks.
  # Collapsing whitespace here makes AU/JT/TA lines disappear and causes DOI/SO fragments to become nodes.
  refs <- enc2utf8(as.character(refs))
  refs <- trimws(refs)
  refs <- refs[nzchar(refs)]
  if (!length(refs)) stop("No usable reference entries detected.", call. = FALSE)

  parsed <- lapply(refs, .ref_parse_one)
  journal <- vapply(parsed, function(z) z$journal %||% NA_character_, character(1))
  authors_list <- lapply(parsed, function(z) unique(trimws(as.character(z$authors %||% character(0)))))
  authors_list <- lapply(authors_list, function(v) v[nzchar(v) & !is.na(v)])

  if (isTRUE(fallback_author_to_journal)) {
    for (i in seq_along(journal)) {
      if (is.na(journal[i]) || !nzchar(journal[i]) || .ref_is_doi_fragment(journal[i])) {
        first_author <- if (length(authors_list[[i]]) >= 1) authors_list[[i]][1] else NA_character_
        if (!is.na(first_author) && nzchar(first_author)) {
          journal[i] <- first_author
          authors_list[[i]] <- authors_list[[i]][-1]
        }
      }
    }
  }

  max_authors <- max(1L, max(lengths(authors_list), na.rm = TRUE))
  out <- data.frame(Journal = journal, stringsAsFactors = FALSE, check.names = FALSE)
  for (k in seq_len(max_authors)) {
    out[[paste0("Author_", k)]] <- vapply(authors_list, function(v) if (length(v) >= k) v[k] else NA_character_, character(1))
  }
  out$Journal <- trimws(as.character(out$Journal))
  out$Journal[is.na(out$Journal) | .ref_is_doi_fragment(out$Journal)] <- ""
  if ("Author_1" %in% names(out)) {
    same <- !is.na(out$Author_1) & nzchar(out$Author_1) & out$Author_1 == out$Journal
    out$Author_1[same] <- NA_character_
  }
  out[is.na(out)] <- ""
  out
}

.ref_valid_node <- function(x) {
  x <- trimws(gsub("^doi\\s*:?\\s*", "", as.character(x %||% ""), ignore.case = TRUE, perl = TRUE))
  if (!nzchar(x) || .ref_is_doi_fragment(x)) return(FALSE)
  if (grepl("^PMID|^PMCID", x, ignore.case = TRUE)) return(FALSE)
  grepl("[A-Za-z]", x)
}

.ref_build_nodes_edges <- function(wide, top_n = Inf) {
  rows_terms <- lapply(seq_len(nrow(wide)), function(i) {
    # Keep Journal and Author columns explicitly.  This is important for
    # Google Scholar profile rows because the source line is the journal/source
    # term and should participate in the co-word network with authors.
    j <- character(0)
    if ("Journal" %in% names(wide)) {
      j <- trimws(as.character(wide[i, "Journal", drop = TRUE]))
      j <- j[!is.na(j) & nzchar(j)]
      # Journal names are already cleaned by .ref_gs_clean_source_to_journal()
      # or .ref_clean_journal_strict(); do not reject them with title-fragment
      # rules that can remove lowercase real journal names.
      j <- j[!vapply(j, .ref_is_doi_fragment, logical(1)) & !vapply(j, .ref_is_access_fragment, logical(1))]
    }
    author_cols <- grep("^Author_", names(wide), value = TRUE)
    a <- character(0)
    if (length(author_cols)) {
      a <- unlist(wide[i, author_cols, drop = FALSE], use.names = FALSE)
      a <- trimws(as.character(a))
      a <- a[!is.na(a) & nzchar(a)]
      a <- a[vapply(a, .ref_valid_node, logical(1))]
    }
    v <- unique(c(j, a))
    v
  })
  rows_terms <- rows_terms[lengths(rows_terms) > 0]
  if (!length(rows_terms)) stop("No usable author/journal terms after cleaning.", call. = FALSE)

  # Keep a journal dictionary so real journal/source names are not removed by
  # author-oriented node filters.  This is essential for Google Scholar profile
  # rows, where the third line is the journal/source term.
  journal_set <- character(0)
  if ("Journal" %in% names(wide)) {
    journal_set <- trimws(as.character(wide$Journal))
    journal_set <- journal_set[!is.na(journal_set) & nzchar(journal_set)]
    journal_set <- unique(journal_set)
  }

  tb <- sort(table(unlist(rows_terms, use.names = FALSE)), decreasing = TRUE)
  nodes <- data.frame(name = names(tb), value = as.numeric(tb), stringsAsFactors = FALSE, check.names = FALSE)
  nodes$term_type <- ifelse(nodes$name %in% journal_set, "Journal", "Author")
  ok_node <- nodes$term_type == "Journal" | vapply(nodes$name, .ref_valid_node, logical(1))
  nodes <- nodes[ok_node, , drop = FALSE]
  nodes <- nodes[order(-nodes$value, nodes$term_type, nodes$name), , drop = FALSE]
  if (is.finite(top_n) && nrow(nodes) > top_n) nodes <- nodes[seq_len(top_n), , drop = FALSE]

  selected_all <- nodes$name
  edge_list <- list(); kk <- 0L
  for (terms in rows_terms) {
    terms <- sort(unique(terms[terms %in% selected_all]))
    if (length(terms) < 2) next
    for (pair in utils::combn(terms, 2, simplify = FALSE)) {
      kk <- kk + 1L
      edge_list[[kk]] <- data.frame(Leader = pair[1], Follower = pair[2], WCD = 1, stringsAsFactors = FALSE)
    }
  }
  if (length(edge_list)) {
    edge_df <- dplyr::bind_rows(edge_list)
    edges <- stats::aggregate(WCD ~ Leader + Follower, data = edge_df, FUN = sum)
    edges <- edges[order(-edges$WCD, edges$Leader, edges$Follower), , drop = FALSE]
  } else {
    edges <- data.frame(Leader = character(0), Follower = character(0), WCD = numeric(0), stringsAsFactors = FALSE)
  }
  nodes$carac <- 1L
  nodes$value2 <- nodes$value
  rownames(nodes) <- NULL
  rownames(edges) <- NULL
  list(nodes = nodes, edges = edges, rows_terms = rows_terms)
}

.ref_run_flca_top20 <- function(nodes, edges) {
  if (!is.data.frame(nodes) || !nrow(nodes)) stop("No nodes available for FLCA.", call. = FALSE)
  if (!is.data.frame(edges)) edges <- data.frame(Leader=character(0), Follower=character(0), WCD=numeric(0))

  nodes <- as.data.frame(nodes, stringsAsFactors = FALSE, check.names = FALSE)
  if (!"name" %in% names(nodes)) stop("nodes must contain name.", call. = FALSE)
  if (!"value" %in% names(nodes)) nodes$value <- 1
  nodes$name <- trimws(as.character(nodes$name))
  nodes$value <- suppressWarnings(as.numeric(nodes$value))
  nodes$value[!is.finite(nodes$value)] <- 1
  keep_cols <- intersect(c("name", "value", "term_type"), names(nodes))
  nodes <- nodes[nzchar(nodes$name), keep_cols, drop = FALSE]
  if (!"term_type" %in% names(nodes)) nodes$term_type <- "Author"
  nodes$term_type[is.na(nodes$term_type) | !nzchar(nodes$term_type)] <- "Author"
  nodes <- nodes[!duplicated(nodes$name), , drop = FALSE]

  edges <- as.data.frame(edges, stringsAsFactors = FALSE, check.names = FALSE)
  if ("follower" %in% names(edges) && !("Follower" %in% names(edges))) names(edges)[names(edges) == "follower"] <- "Follower"
  if (!("Leader" %in% names(edges)) && "Source" %in% names(edges)) names(edges)[names(edges) == "Source"] <- "Leader"
  if (!("Follower" %in% names(edges)) && "Target" %in% names(edges)) names(edges)[names(edges) == "Target"] <- "Follower"
  if (!("WCD" %in% names(edges))) edges$WCD <- 1
  if (!all(c("Leader","Follower","WCD") %in% names(edges))) edges <- data.frame(Leader=character(0), Follower=character(0), WCD=numeric(0))
  edges$Leader <- trimws(as.character(edges$Leader))
  edges$Follower <- trimws(as.character(edges$Follower))
  edges$WCD <- suppressWarnings(as.numeric(edges$WCD))
  edges$WCD[!is.finite(edges$WCD)] <- 1
  edges <- edges[nzchar(edges$Leader) & nzchar(edges$Follower) & edges$Leader != edges$Follower, c("Leader","Follower","WCD"), drop = FALSE]
  edges <- edges[edges$Leader %in% nodes$name & edges$Follower %in% nodes$name, , drop = FALSE]

  # Google Scholar / pasted author-journal mode:
  # use the raw author+journal frequency Top20 as the source of truth so high-frequency
  # journals (e.g., Medicine) cannot be dropped by major-sampling of author clusters.
  # We still compute value2, clusters, and silhouette-like columns for the same
  # Network / SSplot / Kano downstream panels.
  if ("term_type" %in% names(nodes) && any(nodes$term_type == "Journal", na.rm = TRUE)) {
    if (nrow(edges) > 0) {
      strength_df <- rbind(
        data.frame(name = edges$Leader,   value2 = edges$WCD, stringsAsFactors = FALSE),
        data.frame(name = edges$Follower, value2 = edges$WCD, stringsAsFactors = FALSE)
      )
      strength_df <- stats::aggregate(value2 ~ name, strength_df, sum)
      nodes$value2 <- strength_df$value2[match(nodes$name, strength_df$name)]
      nodes$value2[!is.finite(nodes$value2) | is.na(nodes$value2)] <- nodes$value[!is.finite(nodes$value2) | is.na(nodes$value2)]
    } else {
      nodes$value2 <- nodes$value
    }
    # Frequency-first Top20; journals and authors are allowed together.
    nodes <- nodes[order(-nodes$value, -nodes$value2, nodes$term_type, nodes$name), , drop = FALSE]
    nd <- utils::head(nodes, 20)
    selected <- as.character(nd$name)
    ed <- edges[edges$Leader %in% selected & edges$Follower %in% selected, , drop = FALSE]

    k <- min(5L, nrow(nd))
    leaders <- nd$name[seq_len(k)]
    nd$carac <- match(nd$name, leaders)
    nd$carac[is.na(nd$carac)] <- NA_integer_
    if (nrow(ed) > 0 && length(leaders) > 0) {
      for (ii in which(is.na(nd$carac))) {
        nm <- nd$name[ii]
        hit <- ed[(ed$Leader == nm & ed$Follower %in% leaders) | (ed$Follower == nm & ed$Leader %in% leaders), , drop = FALSE]
        if (nrow(hit) > 0) {
          hit$leader_hit <- ifelse(hit$Leader %in% leaders, hit$Leader, hit$Follower)
          hit <- hit[order(-hit$WCD), , drop = FALSE]
          nd$carac[ii] <- match(hit$leader_hit[1], leaders)
        }
      }
    }
    nd$carac[is.na(nd$carac)] <- ((seq_len(nrow(nd))[is.na(nd$carac)] - 1L) %% max(1L, k)) + 1L
    nd$ssi <- 0; nd$a_i <- 0; nd$b_i <- 0; nd$a_star1 <- 0
    if (exists("compute_silhouette_df", mode = "function") && nrow(ed) > 0 && length(unique(nd$carac)) >= 2) {
      sil_df <- tryCatch({
        ed2 <- ed
        if ("Follower" %in% names(ed2) && !("follower" %in% names(ed2))) ed2$follower <- ed2$Follower
        compute_silhouette_df(nd, ed2)
      }, error = function(e) NULL)
      if (is.data.frame(sil_df) && nrow(sil_df) && "name" %in% names(sil_df)) {
        sx <- match(nd$name, as.character(sil_df$name))
        for (cc in c("ssi", "sil_width", "a_i", "b_i", "a_star1", "a_star")) {
          if (cc %in% names(sil_df)) {
            val <- suppressWarnings(as.numeric(sil_df[[cc]][sx]))
            val[!is.finite(val)] <- 0
            if (cc == "sil_width") nd$ssi <- val
            else if (cc == "a_star") nd$a_star1 <- val
            else nd[[cc]] <- val
          }
        }
      }
    }
    nd$SSi <- nd$ssi
    nd$a_star <- nd$a_star1
    nd <- nd[order(-nd$value, -nd$value2, nd$term_type, nd$name), , drop = FALSE]
    rownames(nd) <- NULL
    rownames(ed) <- NULL
    return(list(nodes = nd, edges = ed, engine = "Journal/Author Frequency Top20 + FLCA-style clustering"))
  }

  # Prefer the existing app wrapper, which calls flca_ms_sil_module.R / flca_ma_sil_module.R.
  if (exists("run_flca_ms_sil_internal", mode = "function") && nrow(edges) > 0 && nrow(nodes) >= 3) {
    out <- tryCatch(
      run_flca_ms_sil_internal(nodes, edges, cfg = list(top_clusters = 5, base_per_cluster = 4, target_n = 20), verbose = FALSE),
      error = function(e) {
        message("[AMA FLCA-MA-SIL] run_flca_ms_sil_internal failed: ", conditionMessage(e))
        NULL
      }
    )
    if (!is.null(out) && is.data.frame(out$nodes) && nrow(out$nodes)) {
      nd <- out$nodes
      ed <- out$data
      if ("follower" %in% names(ed) && !("Follower" %in% names(ed))) names(ed)[names(ed) == "follower"] <- "Follower"
      return(list(nodes = nd, edges = ed, engine = "FLCA-MA-SIL"))
    }
  }

  # Direct module runner fallback if the internal wrapper was not available.
  if (exists("run_flca_ms_sil", mode = "function") && nrow(edges) > 0 && nrow(nodes) >= 3) {
    edges_sym <- rbind(
      data.frame(Leader = edges$Leader, follower = edges$Follower, WCD = edges$WCD, stringsAsFactors = FALSE),
      data.frame(Leader = edges$Follower, follower = edges$Leader, WCD = edges$WCD, stringsAsFactors = FALSE)
    )
    out <- tryCatch(run_flca_ms_sil(list(nodes = nodes, data = edges_sym)), error = function(e) {
      message("[AMA FLCA-MA-SIL] direct module failed: ", conditionMessage(e)); NULL
    })
    if (!is.null(out) && is.data.frame(out$nodes) && nrow(out$nodes)) {
      nd <- out$nodes
      if (!"value2" %in% names(nd)) nd$value2 <- nd$value
      if (!"ssi" %in% names(nd)) nd$ssi <- 0
      if (!"a_i" %in% names(nd)) nd$a_i <- 0
      if (!"b_i" %in% names(nd)) nd$b_i <- 0
      if (!"a_star1" %in% names(nd)) nd$a_star1 <- 0
      nd <- nd[order(-suppressWarnings(as.numeric(nd$value)), as.character(nd$name)), , drop = FALSE]
      nd <- utils::head(nd, 20)
      selected <- as.character(nd$name)
      ed <- edges[edges$Leader %in% selected & edges$Follower %in% selected, , drop = FALSE]
      return(list(nodes = nd, edges = ed, engine = "FLCA-MA-SIL"))
    }
  }

  # ------------------------------------------------------------
  # Built-in fallback: FLCA-MA-SIL-like Top20 if external module is absent.
  # This prevents the Shiny app from shutting down when flca_ms_sil_module.R
  # is not present in the app folder. It still forces the same Top20 source
  # for Network, SSplot, Kano, node table, and edge table.
  # ------------------------------------------------------------
  message("[AMA FLCA-MA-SIL] External module unavailable; using built-in safe FLCA-MA-SIL fallback.")

  # value2 = incident edge strength, then major-sampling Top20 by value/value2
  if (nrow(edges) > 0) {
    strength_df <- rbind(
      data.frame(name = edges$Leader,   value2 = edges$WCD, stringsAsFactors = FALSE),
      data.frame(name = edges$Follower, value2 = edges$WCD, stringsAsFactors = FALSE)
    )
    strength_df <- stats::aggregate(value2 ~ name, strength_df, sum)
    nodes$value2 <- strength_df$value2[match(nodes$name, strength_df$name)]
    nodes$value2[!is.finite(nodes$value2) | is.na(nodes$value2)] <- nodes$value[!is.finite(nodes$value2) | is.na(nodes$value2)]
  } else {
    nodes$value2 <- nodes$value
  }

  nodes <- nodes[order(-nodes$value, -nodes$value2, nodes$name), , drop = FALSE]
  nd <- utils::head(nodes, 20)
  selected <- as.character(nd$name)
  ed <- edges[edges$Leader %in% selected & edges$Follower %in% selected, , drop = FALSE]

  # Choose up to 5 leaders and assign followers by strongest connection to leaders.
  k <- min(5L, nrow(nd))
  leaders <- nd$name[seq_len(k)]
  nd$carac <- match(nd$name, leaders)
  nd$carac[is.na(nd$carac)] <- NA_integer_

  if (nrow(ed) > 0 && length(leaders) > 0) {
    for (ii in which(is.na(nd$carac))) {
      nm <- nd$name[ii]
      hit <- ed[(ed$Leader == nm & ed$Follower %in% leaders) | (ed$Follower == nm & ed$Leader %in% leaders), , drop = FALSE]
      if (nrow(hit) > 0) {
        hit$leader_hit <- ifelse(hit$Leader %in% leaders, hit$Leader, hit$Follower)
        hit <- hit[order(-hit$WCD), , drop = FALSE]
        nd$carac[ii] <- match(hit$leader_hit[1], leaders)
      }
    }
  }
  nd$carac[is.na(nd$carac)] <- ((seq_len(nrow(nd))[is.na(nd$carac)] - 1L) %% max(1L, k)) + 1L

  # Silhouette-like defaults. If compute_silhouette_df exists and edges are enough, use it safely.
  nd$ssi <- 0; nd$a_i <- 0; nd$b_i <- 0; nd$a_star1 <- 0
  if (exists("compute_silhouette_df", mode = "function") && nrow(ed) > 0 && length(unique(nd$carac)) >= 2) {
    sil_df <- tryCatch({
      ed2 <- ed
      if ("Follower" %in% names(ed2) && !("follower" %in% names(ed2))) ed2$follower <- ed2$Follower
      compute_silhouette_df(nd, ed2)
    }, error = function(e) {
      message("[AMA FLCA-MA-SIL fallback] silhouette failed: ", conditionMessage(e)); NULL
    })
    if (is.data.frame(sil_df) && nrow(sil_df) && "name" %in% names(sil_df)) {
      sx <- match(nd$name, as.character(sil_df$name))
      for (cc in c("ssi", "sil_width", "a_i", "b_i", "a_star1", "a_star")) {
        if (cc %in% names(sil_df)) {
          val <- suppressWarnings(as.numeric(sil_df[[cc]][sx]))
          val[!is.finite(val)] <- 0
          if (cc == "sil_width") nd$ssi <- val
          else if (cc == "a_star") nd$a_star1 <- val
          else nd[[cc]] <- val
        }
      }
    }
  }
  nd$SSi <- nd$ssi
  nd$a_star <- nd$a_star1
  nd <- nd[order(-nd$value, -nd$value2, nd$name), , drop = FALSE]
  rownames(nd) <- NULL
  rownames(ed) <- NULL
  return(list(nodes = nd, edges = ed, engine = "Internal FLCA-MA-SIL fallback"))
}


# ---- Safe coercion for AMA author/journal plots ----
# Prevents htmltools/grid/text functions from receiving list/data-frame/factor labels.
.ref_safe_text <- function(x) {
  x <- as.character(unlist(x, use.names = FALSE))
  x[is.na(x)] <- ""
  x <- enc2utf8(x)
  x <- gsub("[\r\n\t]+", " ", x, perl = TRUE)
  x <- gsub("\\s+", " ", x, perl = TRUE)
  trimws(x)
}

.ref_safe_num <- function(x, default = 0) {
  y <- suppressWarnings(as.numeric(x))
  y[!is.finite(y) | is.na(y)] <- default
  y
}

.ref_sanitize_plot_payload <- function(nodes, edges) {
  nodes <- as.data.frame(nodes, stringsAsFactors = FALSE, check.names = FALSE)
  if (!"name" %in% names(nodes)) nodes$name <- seq_len(nrow(nodes))
  nodes$name <- .ref_safe_text(nodes$name)
  nodes <- nodes[nzchar(nodes$name), , drop = FALSE]
  nodes <- nodes[!duplicated(nodes$name), , drop = FALSE]

  if (!"value" %in% names(nodes)) nodes$value <- 1
  nodes$value <- .ref_safe_num(nodes$value, default = 1)
  if (!"value2" %in% names(nodes)) nodes$value2 <- nodes$value
  nodes$value2 <- .ref_safe_num(nodes$value2, default = nodes$value)
  if (!"carac" %in% names(nodes)) nodes$carac <- 1
  nodes$carac <- as.integer(.ref_safe_num(nodes$carac, default = 1))
  nodes$carac[is.na(nodes$carac)] <- 1L
  # compact cluster labels for plotting/display; keeps FLCA grouping but avoids non-consecutive IDs
  nodes$carac <- as.integer(factor(as.character(nodes$carac), levels = unique(as.character(nodes$carac))))
  if (!"ssi" %in% names(nodes) && "SSi" %in% names(nodes)) nodes$ssi <- nodes$SSi
  if (!"ssi" %in% names(nodes)) nodes$ssi <- 0
  nodes$ssi <- .ref_safe_num(nodes$ssi, default = 0)
  if (!"a_star1" %in% names(nodes) && "a_star" %in% names(nodes)) nodes$a_star1 <- nodes$a_star
  if (!"a_star1" %in% names(nodes)) nodes$a_star1 <- 0
  nodes$a_star1 <- .ref_safe_num(nodes$a_star1, default = 0)

  nodes <- nodes[order(-nodes$value, nodes$name), , drop = FALSE]
  rownames(nodes) <- NULL

  if (is.null(edges) || !is.data.frame(edges) || !nrow(edges)) {
    edges <- data.frame(Leader = character(0), Follower = character(0), WCD = numeric(0), stringsAsFactors = FALSE)
  } else {
    edges <- as.data.frame(edges, stringsAsFactors = FALSE, check.names = FALSE)
    if ("follower" %in% names(edges) && !"Follower" %in% names(edges)) names(edges)[names(edges) == "follower"] <- "Follower"
    if (!"Leader" %in% names(edges) && "from" %in% names(edges)) names(edges)[names(edges) == "from"] <- "Leader"
    if (!"Follower" %in% names(edges) && "to" %in% names(edges)) names(edges)[names(edges) == "to"] <- "Follower"
    if (!"WCD" %in% names(edges)) edges$WCD <- 1
    edges$Leader <- .ref_safe_text(edges$Leader)
    edges$Follower <- .ref_safe_text(edges$Follower)
    edges$WCD <- .ref_safe_num(edges$WCD, default = 1)
    edges <- edges[nzchar(edges$Leader) & nzchar(edges$Follower) & edges$Leader %in% nodes$name & edges$Follower %in% nodes$name, , drop = FALSE]
    edges <- edges[edges$Leader != edges$Follower, , drop = FALSE]
    rownames(edges) <- NULL
  }

  list(nodes = nodes, edges = edges)
}

.ref_format_numeric_2 <- function(df) {
  if (is.null(df) || !is.data.frame(df) || !nrow(df)) return(df)
  int_cols <- intersect(c("rank", "carac", "cluster", "membership"), names(df))
  num_cols <- names(df)[vapply(df, is.numeric, logical(1))]
  for (cc in num_cols) {
    if (cc %in% int_cols) {
      df[[cc]] <- as.integer(round(df[[cc]], 0))
    } else {
      df[[cc]] <- round(df[[cc]], 2)
    }
  }
  df
}

.ref_top20_check_df <- function(r) {
  if (is.null(r) || is.null(r$nodes) || !is.data.frame(r$nodes)) {
    return(data.frame(Message = "Click Run AMA author/journal extraction first.", stringsAsFactors = FALSE))
  }
  nd <- .ref_sanitize_plot_payload(r$nodes, r$edges)$nodes
  if (!is.data.frame(nd) || !nrow(nd)) {
    return(data.frame(Message = "No Top20 nodes available.", stringsAsFactors = FALSE))
  }
  # carac is a cluster label from FLCA-MA-SIL. Recode only the labels to compact 1..K
  # for display and plotting; this avoids confusing large/non-consecutive cluster IDs.
  if ("carac" %in% names(nd)) {
    nd$carac <- as.integer(factor(as.character(nd$carac), levels = unique(as.character(nd$carac))))
  }
  keep <- intersect(c("name", "value", "value2", "carac", "ssi", "a_i", "b_i", "a_star1"), names(nd))
  out <- nd[, keep, drop = FALSE]
  out$rank <- seq_len(nrow(out))
  out <- out[, c("rank", setdiff(names(out), "rank")), drop = FALSE]
  .ref_format_numeric_2(out)
}



# ============================================================
# STRICT PubMed/AMA author+journal parser override
# Purpose: ensure Top20 contains ONLY journal titles and author names.
# It rejects title words, DOI fragments, SO/PST/AID lines, and one-letter nodes.
# ============================================================
.ref_norm_spaces2 <- function(x) {
  x <- enc2utf8(as.character(x %||% ""))
  x[is.na(x)] <- ""
  x <- gsub("[‘’ʼ]", "'", x, perl = TRUE)
  x <- gsub('[“”]', '"', x, perl = TRUE)
  x <- gsub("[–—−]", "-", x, perl = TRUE)
  x <- gsub("\\s+", " ", x, perl = TRUE)

  trimws(x)
}

# ---- Google Scholar profile pasted-text parser ----
# Supports rows copied from a Google Scholar author profile table:
#   Title
#   AU1, AU2, AU3
#   Journal/source volume(issue), pages
#   cited_by<TAB>year    OR    year
# It is intentionally local/offline: paste or browser-export text into the textarea.
.ref_gs_is_year_line <- function(x) {
  z <- .ref_norm_spaces2(x)
  grepl("^(?:[0-9,]+[[:space:]]+)?(?:19|20)\\d{2}$", z, perl = TRUE) ||
    grepl("^[0-9,]+[[:space:]]*\\t[[:space:]]*(?:19|20)\\d{2}$", z, perl = TRUE)
}

.ref_gs_author_line_like <- function(x) {
  z <- .ref_norm_spaces2(x)
  if (!nzchar(z)) return(FALSE)
  if (.ref_gs_is_year_line(z)) return(FALSE)
  if (grepl("^(Journals|Articles|Title|Cited by|Year)$", z, ignore.case = TRUE, perl = TRUE)) return(FALSE)
  # Google Scholar author rows are usually initials + surname separated by commas;
  # allow a single-author row such as "TW Chien" as well.
  has_comma <- grepl(",", z, fixed = TRUE)
  m <- gregexpr("\\b[A-Z]{1,5}[[:space:]]+[A-Z][A-Za-z'\\-]+\\b", z, perl = TRUE)[[1]]
  n_initial_names <- if (length(m) == 1L && m[1] == -1L) 0L else length(m)
  has_initial_name <- n_initial_names >= 1L
  has_initial_name && (has_comma || n_initial_names == 1L)
}

.ref_gs_clean_source_to_journal <- function(x) {
  z0 <- .ref_norm_spaces2(x)
  if (!nzchar(z0)) return(NA_character_)
  z0 <- gsub("[[:space:]]*…[[:space:]]*$", "", z0, perl = TRUE)
  z0 <- gsub("[[:space:]]*\\.\\.\\.[[:space:]]*$", "", z0, perl = TRUE)
  # Remove common Google Scholar source suffixes: volume(issue), pages / volume, pages / volume, article id.
  j <- sub("[[:space:]]+[0-9]+[A-Za-z]?[[:space:]]*(\\([^)]*\\))?[[:space:]]*,.*$", "", z0, perl = TRUE)
  j <- sub("[[:space:]]+[0-9]+[A-Za-z]?[[:space:]]*(\\([^)]*\\))?[[:space:]]*$", "", j, perl = TRUE)
  j <- sub("[[:space:]]+e[0-9]{3,}.*$", "", j, perl = TRUE, ignore.case = TRUE)
  j <- .ref_norm_spaces2(j)
  if (!nzchar(j)) j <- z0
  j <- .ref_clean_journal_symbols(j)
  j <- trimws(gsub("[.,;:/／；：，]+$", "", j, perl = TRUE))
  if (!nzchar(j) || .ref_is_access_fragment(j) || .ref_is_doi_fragment(j)) return(NA_character_)
  j
}

.ref_gs_split_authors <- function(x) {
  z <- .ref_norm_spaces2(x)
  z <- gsub("\\s+(and|&)\\s+", ", ", z, ignore.case = TRUE, perl = TRUE)
  z <- gsub("\\bet\\s+al\\.?", "", z, ignore.case = TRUE, perl = TRUE)
  z <- gsub("\\.{2,}|…", "", z, perl = TRUE)
  z <- gsub("[。；;、]", ",", z, perl = TRUE)
  parts <- trimws(unlist(strsplit(z, "\\s*,\\s*", perl = TRUE), use.names = FALSE))
  parts <- parts[nzchar(parts)]
  parts <- parts[!grepl("^\\.$", parts)]
  unique(na.omit(.ref_clean_author_strict(parts)))
}

.ref_parse_google_scholar_block <- function(block) {
  lines <- unlist(strsplit(gsub("\\r\\n|\\r", "\\n", as.character(block %||% "")), "\\n", fixed = TRUE), use.names = FALSE)
  lines <- .ref_norm_spaces2(lines)
  lines <- lines[nzchar(lines)]
  if (length(lines) < 4L) return(NULL)
  # A Scholar block is title / authors / source / citation-year.
  if (!.ref_gs_author_line_like(lines[2]) || !.ref_gs_is_year_line(lines[4])) return(NULL)
  authors <- .ref_gs_split_authors(lines[2])
  journal <- .ref_gs_clean_source_to_journal(lines[3])
  if ((!length(authors)) && (is.na(journal) || !nzchar(journal))) return(NULL)
  list(journal = ifelse(length(journal) && !is.na(journal), journal, ""), authors = authors)
}

.ref_split_google_scholar_profile_blocks <- function(txt) {
  txt <- enc2utf8(as.character(txt %||% ""))
  txt <- gsub("\r\n|\r", "\n", txt, perl = TRUE)
  lines <- trimws(unlist(strsplit(txt, "\n", fixed = TRUE), use.names = FALSE))
  lines <- .ref_norm_spaces2(lines)
  lines <- lines[nzchar(lines)]
  if (length(lines) < 4L) return(character(0))
  # Drop common copied table headers/buttons.
  bad <- grepl("^(Title|Authors|Publication|Journal|Cited by|Year|All|Since|Sort by|Show more|顯示更多|匯出|Export)$", lines, ignore.case = TRUE, perl = TRUE)
  lines <- lines[!bad]
  out <- character(0)
  i <- 1L
  while (i <= length(lines) - 3L) {
    if (.ref_gs_author_line_like(lines[i + 1L]) && .ref_gs_is_year_line(lines[i + 3L])) {
      out <- c(out, paste(lines[i:(i + 3L)], collapse = "\n"))
      i <- i + 4L
    } else {
      i <- i + 1L
    }
  }
  out <- out[nzchar(out)]
  # Require at least two records to switch the whole textarea to Scholar mode.
  if (length(out) >= 2L) out else character(0)
}

.ref_split_pubmed_blocks_strict <- function(txt) {
  txt <- enc2utf8(as.character(txt %||% ""))
  txt <- gsub("\r\n|\r", "\n", txt, perl = TRUE)
  txt <- trimws(txt)
  if (!nzchar(txt)) return(character(0))

  # Google Scholar author-profile pasted rows: title / authors / source / cited-year.
  gs_blocks <- .ref_split_google_scholar_profile_blocks(txt)
  if (length(gs_blocks) >= 2L) return(gs_blocks)

  # Best case: PubMed MEDLINE/NBIB blocks with PMID-
  if (grepl("(?m)^PMID\\s*-", txt, perl = TRUE)) {
    z <- unlist(strsplit(txt, "(?m)(?=^PMID\\s*-)", perl = TRUE), use.names = FALSE)
    z <- trimws(z)
    return(z[nzchar(z)])
  }

  # PubMed-like field blocks without PMID: collect until SO - line.
  if (grepl("(?m)^(AU|FAU|JT|TA|SO|AID|PST)\\s*-", txt, perl = TRUE)) {
    lines <- unlist(strsplit(txt, "\n", fixed = TRUE), use.names = FALSE)
    blocks <- character(0); cur <- character(0); seen_field <- FALSE
    for (ln in lines) {
      if (grepl("^(AU|FAU|JT|TA|SO|AID|PST|TI|DP)\\s*-", ln, perl = TRUE)) seen_field <- TRUE
      if (seen_field) cur <- c(cur, ln)
      if (grepl("^SO\\s*-", ln, perl = TRUE) && length(cur)) {
        blocks <- c(blocks, paste(cur, collapse = "\n"))
        cur <- character(0); seen_field <- FALSE
      }
    }
    if (length(cur)) blocks <- c(blocks, paste(cur, collapse = "\n"))
    blocks <- trimws(blocks)
    blocks <- blocks[nzchar(blocks)]
    if (length(blocks)) return(blocks)
  }

  # AMA/Vancouver numbered references. Supports:
  # [1] text, 1. text, 1: text, 1) text, and mixed Chinese/English references.
  marker <- "(?:\\[\\d+\\]|\\(?\\d+\\)?|\\d+)\\s*[:;.\\)]?\\s+"
  m <- gregexpr(paste0("(?ms)^\\s*", marker, ".*?(?=(?:\n\\s*", marker, ")|\\Z)"), txt, perl = TRUE)
  if (m[[1]][1] != -1) {
    refs <- regmatches(txt, m)[[1]]
    refs <- trimws(gsub("\\s*\n\\s*", " ", refs, perl = TRUE))
    refs <- refs[nzchar(refs)]
    if (length(refs)) return(refs)
  }

  # Blank-line records, useful for pasted references without numbering.
  chunks <- unlist(strsplit(txt, "\n\\s*\n", perl = TRUE), use.names = FALSE)
  chunks <- trimws(chunks); chunks <- chunks[nzchar(chunks)]
  if (length(chunks) > 1L) return(chunks)

  # Last fallback: one reference per non-empty line.
  lines <- trimws(unlist(strsplit(txt, "\n", fixed = TRUE), use.names = FALSE))
  lines <- lines[nzchar(lines)]
  if (length(lines) > 1L) return(lines)

  character(0)
}


.ref_get_fields_strict <- function(block, tag) {
  lines <- unlist(strsplit(gsub("\r\n|\r", "\n", as.character(block)), "\n", fixed = TRUE), use.names = FALSE)
  vals <- character(0); cur <- NULL
  for (ln in lines) {
    if (grepl(paste0("^", tag, "\\s*-"), ln, perl = TRUE)) {
      if (!is.null(cur)) vals <- c(vals, cur)
      cur <- sub(paste0("^", tag, "\\s*-\\s*"), "", ln, perl = TRUE)
    } else if (!is.null(cur) && grepl("^\\s{2,}\\S", ln, perl = TRUE) && !grepl("^[A-Z0-9]{2,4}\\s*-", ln, perl = TRUE)) {
      cur <- paste(cur, trimws(ln))
    } else if (grepl("^[A-Z0-9]{2,4}\\s*-", ln, perl = TRUE)) {
      if (!is.null(cur)) vals <- c(vals, cur)
      cur <- NULL
    }
  }
  if (!is.null(cur)) vals <- c(vals, cur)
  vals <- .ref_norm_spaces2(vals)
  vals[nzchar(vals)]
}

.ref_journal_from_so_strict <- function(so) {
  so <- .ref_norm_spaces2(so)
  if (!length(so) || !nzchar(so[1])) return(NA_character_)
  x <- so[1]
  # SO - Medicine (Baltimore). 2025 Oct 24;104(43):e...
  # Also accept SO strings without the period before year, or with slash/space after year:
  #   Journal Name 2025;...   Journal Name. 2025/...   Journal Name. 2025 Nov 3;...
  month_pat <- "(?:Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Sept|Oct|Nov|Dec)[A-Za-z.]*"
  date_pat <- paste0("(?:19|20)\\d{2}(?:[[:space:]]+", month_pat, ")?(?:[[:space:]]+[0-3]?\\d)?")
  j <- sub(paste0("(?:\\.[[:space:]]*|[[:space:]]+)", date_pat, "(?:[[:space:]]*(?:[;,:/／；：，]|$)|[[:space:]]+).*$"), "", x, perl = TRUE, ignore.case = TRUE)
  j <- .ref_norm_spaces2(j)
  if (!nzchar(j)) NA_character_ else j
}

.ref_clean_author_strict <- function(x) {
  x <- .ref_norm_spaces2(x)
  x <- gsub("^[[:punct:] ]+|[[:punct:] ]+$", "", x, perl = TRUE)
  x <- trimws(gsub("\\s+", " ", x, perl = TRUE))

  # English AMA/Google Scholar style:
  #   Aiken LH, TW Chien, HF Lee, W Chou, R Core Team, Posit PBC.
  # PATCH 2026-06-02: allow one-letter initial + surname (e.g., "W Chou").
  # The old pattern required at least two leading letters, so "W Chou" was
  # dropped; then the previous valid middle author could be mistaken as "last".
  x <- gsub("\\b([A-Z])\\.[[:space:]]*", "\\1 ", x, perl = TRUE)
  x <- trimws(gsub("\\s+", " ", x, perl = TRUE))
  ok_en <- grepl("^[A-Z][A-Za-z' -]{1,50}(?:\\s+[A-Z][A-Za-z' -]{0,20}|\\s+[A-Z]{1,8})+$", x, perl = TRUE)
  ok_initial_surname <- grepl("^[A-Z]{1,8}\\s+[A-Z][A-Za-z'\\-]{1,50}(?:\\s+[A-Z][A-Za-z'\\-]{1,50})?$", x, perl = TRUE)
  ok_en <- ok_en | ok_initial_surname
  # Chinese personal names and organization/group authors.
  ok_zh <- grepl("^[\u4e00-\u9fff]{2,40}$", x, perl = TRUE)

  bad <- is.na(x) | !nzchar(x) |
    grepl("^(URL|http|https|www|doi|pmid|pmcid)\\b", x, ignore.case = TRUE, perl = TRUE) |
    .ref_is_access_fragment(x) |
    grepl("(?:研究|分析|比較|標準公告|公式|資訊公開|評鑑基準|資料開放平台).{4,}", x, perl = TRUE)

  x[bad | !(ok_en | ok_zh)] <- NA_character_
  x
}


.ref_clean_journal_strict <- function(x) {
  x0 <- .ref_norm_spaces2(x)
  x0 <- gsub("\\[Internet\\]", "", x0, ignore.case = TRUE, perl = TRUE)
  x <- .ref_clean_journal_symbols(x0)
  x <- trimws(gsub("[.,;:/／；：，]+$", "", x, perl = TRUE))
  bad <- is.na(x) | !nzchar(x) |
    .ref_is_doi_fragment(x0) | .ref_is_doi_fragment(x) |
    .ref_is_access_fragment(x0) | .ref_is_access_fragment(x) |
    .ref_title_like_journal_fragment(x0) | .ref_title_like_journal_fragment(x) |
    grepl("\\[doi\\]|\\b(PST|SO|AID|PMID|PMCID)\\s*-", x0, ignore.case = TRUE, perl = TRUE) |
    grepl("^(P|S|N)$", x, perl = TRUE) |
    grepl("^[0-9]", x, perl = TRUE) |
    grepl("case fatality|observational study|usability study|meta-analysis", x, ignore.case = TRUE, perl = TRUE)
  x[bad] <- NA_character_
  x
}



# ---- Exact journal-list / journal-frequency-table fallback ----
# Some uploaded text already contains journal names (one per line) or a
# journal-frequency table, rather than full AMA/PubMed references. The strict
# AMA parser intentionally requires a year marker, so those exact journal lines
# were previously dropped. This fallback accepts already-extracted journal names
# only when the input clearly looks like a journal list/table, not a reference list.
.ref_extract_existing_journal_records <- function(txt) {
  empty <- data.frame(reference_no = integer(0), journal = character(0), stringsAsFactors = FALSE)
  x <- enc2utf8(as.character(txt %||% ""))
  x <- gsub("\r\n|\r", "\n", x, perl = TRUE)
  if (!nzchar(trimws(x))) return(empty)

  # Do not use this fallback for true MEDLINE/NBIB blocks or normal references.
  medline_like <- grepl("(?m)^(PMID|AU|FAU|JT|TA|SO|TI|DP|AID|PST)\\s*-", x, perl = TRUE)

  lines <- trimws(unlist(strsplit(x, "\n", fixed = TRUE), use.names = FALSE))
  lines <- lines[nzchar(lines)]
  if (!length(lines)) return(empty)

  # Remove common table headers/footers.
  header_like <- grepl("^\\s*(journal|journal[ _-]*name|source[ _-]*title|frequency|freq|count|n|%)\\b", lines, ignore.case = TRUE, perl = TRUE) |
    grepl("^\\s*total\\b", lines, ignore.case = TRUE, perl = TRUE)
  work <- lines[!header_like]
  work <- trimws(work)
  work <- work[nzchar(work)]
  if (!length(work)) return(empty)

  has_table_header <- any(grepl("\\bjournal\\b", lines, ignore.case = TRUE, perl = TRUE) &
                            grepl("\\b(frequency|freq|count|%)\\b", lines, ignore.case = TRUE, perl = TRUE))

  parse_table_row <- function(ln) {
    z <- trimws(ln)
    z <- gsub("[|]", "\t", z, perl = TRUE)
    z <- gsub(",", "\t", z, perl = TRUE)
    # Prefer separated columns: journal <tab/2+ spaces> frequency [percent]
    m <- regexec("^(.+?)(?:\\t+| {2,})\\s*([0-9]+)(?:\\s+[0-9]+(?:\\.[0-9]+)?%?)?\\s*$", z, perl = TRUE)
    r <- regmatches(z, m)[[1]]
    if (length(r) >= 3) return(list(journal = r[2], freq = suppressWarnings(as.integer(r[3]))))

    # Fallback for copied Word/Excel rows collapsed to single spaces.
    # Keep it table-only to avoid confusing full references with page numbers.
    m <- regexec("^(.+?)\\s+([0-9]+)(?:\\s+[0-9]+(?:\\.[0-9]+)?%?)?\\s*$", z, perl = TRUE)
    r <- regmatches(z, m)[[1]]
    if (length(r) >= 3 && nchar(r[2]) <= 90 && !grepl("(19|20)\\d{2}\\s*[;:/／]", r[2], perl = TRUE)) {
      return(list(journal = r[2], freq = suppressWarnings(as.integer(r[3]))))
    }
    NULL
  }

  parsed <- lapply(work, parse_table_row)
  good <- vapply(parsed, function(z) !is.null(z) && is.finite(z$freq) && z$freq > 0, logical(1))
  table_mode <- !isTRUE(medline_like) && (isTRUE(has_table_header) || sum(good) >= 2L)

  if (isTRUE(table_mode)) {
    parsed <- parsed[good]
    journals <- vapply(parsed, function(z) z$journal, character(1))
    freqs <- vapply(parsed, function(z) z$freq, integer(1))
    journals <- .ref_clean_journal_strict(journals)
    ok <- !is.na(journals) & nzchar(journals) & vapply(journals, .ref_journal_like, logical(1))
    journals <- journals[ok]
    freqs <- freqs[ok]
    if (!length(journals)) return(empty)
    # Replicate by frequency so the existing count() code stays unchanged.
    out <- rep(journals, times = pmax(1L, freqs))
    return(data.frame(reference_no = seq_along(out), journal = out, stringsAsFactors = FALSE))
  }

  # Journal-list mode: exact journal names, one per line, no frequencies.
  # This accepts lines such as "Medicine (Baltimore)" or
  # "Int J Environ Res Public Health" but rejects full article references.
  if (!isTRUE(medline_like)) {
    cand_raw <- work
    too_reference_like <- grepl("(19|20)\\d{2}\\s*[;:/／]", cand_raw, perl = TRUE) |
      grepl("\\bdoi\\b|https?://|PMID|PMCID", cand_raw, ignore.case = TRUE, perl = TRUE) |
      nchar(cand_raw) > 90
    cand <- .ref_clean_journal_strict(cand_raw)
    ok <- !too_reference_like & !is.na(cand) & nzchar(cand) & vapply(cand, .ref_journal_like, logical(1))
    if (sum(ok) >= 1L && mean(ok) >= 0.50) {
      out <- cand[ok]
      return(data.frame(reference_no = seq_along(out), journal = out, stringsAsFactors = FALSE))
    }
  }

  empty
}

.ref_parse_block_strict <- function(block) {
  # Google Scholar author-profile rows pasted from browser/table export.
  gs <- .ref_parse_google_scholar_block(block)
  if (!is.null(gs)) return(gs)

  # MEDLINE/NBIB tags first.
  au <- .ref_get_fields_strict(block, "AU")
  if (!length(au)) au <- .ref_get_fields_strict(block, "FAU")
  au <- unique(na.omit(.ref_clean_author_strict(au)))

  jt <- .ref_get_fields_strict(block, "JT")
  ta <- .ref_get_fields_strict(block, "TA")
  so <- .ref_get_fields_strict(block, "SO")
  journal <- if (length(jt)) jt[1] else if (length(ta)) ta[1] else .ref_journal_from_so_strict(so)
  journal <- .ref_clean_journal_strict(journal)

  # Fallback for AMA/Vancouver single-line reference.
  if ((!length(au) || is.na(journal) || !nzchar(journal)) && !grepl("(?m)^(AU|FAU|JT|TA|SO)\\s*-", block, perl = TRUE)) {
    j2 <- .ref_extract_journal_ama(block)
    journal <- .ref_clean_journal_strict(j2)
    au2 <- .ref_extract_authors_ama(block)
    au <- unique(na.omit(.ref_clean_author_strict(au2)))
  }

  list(journal = ifelse(length(journal) && !is.na(journal), journal, ""), authors = au)
}


# ---- Robust Google Scholar profile pasted-row parser ----
# Handles rows copied from Google Scholar author profile as:
#   title / authors / journal-source / cited-by + year
# This is deliberately independent of the stricter AMA parser so journal/source
# rows such as "Medicine 102 (25), e34063" are kept as journal nodes.
.ref_to_wide_google_scholar_force <- function(txt) {
  txt <- enc2utf8(as.character(txt %||% ""))
  txt <- gsub("\r\n|\r", "\n", txt, perl = TRUE)
  lines <- trimws(unlist(strsplit(txt, "\n", fixed = TRUE), use.names = FALSE))
  lines <- .ref_norm_spaces2(lines)
  lines <- lines[nzchar(lines)]
  if (length(lines) < 4L) return(data.frame())

  bad <- grepl("^(Title|Authors|Publication|Journal|Cited by|Year|All|Since|Sort by|Show more|顯示更多|匯出|Export)$", lines,
               ignore.case = TRUE, perl = TRUE)
  lines <- lines[!bad]
  if (length(lines) < 4L) return(data.frame())

  is_year_line <- function(x) grepl("^(?:[0-9]+[[:space:]]+)?(?:19|20)[0-9]{2}$", trimws(x), perl = TRUE)
  is_author_line <- function(x) {
    z <- trimws(x)
    if (!nzchar(z)) return(FALSE)
    if (grepl("[0-9]{4}|http|doi|PMID|PMCID", z, ignore.case = TRUE, perl = TRUE)) return(FALSE)
    # Initials + surname, Chinese delimiters, or comma-separated byline.
    m <- gregexpr("\\b[A-Z]{1,5}[[:space:]]+[A-Z][A-Za-z'\\-]+\\b", z, perl = TRUE)[[1]]
    n0 <- if (length(m) == 1L && m[1] == -1L) 0L else length(m)
    n0 >= 1L && (grepl(",|；|、", z, perl = TRUE) || n0 == 1L)
  }
  clean_journal <- function(x) {
    z <- .ref_norm_spaces2(x)
    z <- gsub("[[:space:]]*…[[:space:]]*$", "", z, perl = TRUE)
    z <- gsub("[[:space:]]*\\.\\.\\.[[:space:]]*$", "", z, perl = TRUE)
    # Remove volume/issue/pages/article-id after the source name.
    z <- sub("[[:space:]]+[0-9]+[A-Za-z]?[[:space:]]*(\\([^)]*\\))?[[:space:]]*,.*$", "", z, perl = TRUE)
    z <- sub("[[:space:]]+[0-9]+[A-Za-z]?[[:space:]]*(\\([^)]*\\))?[[:space:]]*$", "", z, perl = TRUE)
    z <- sub("[[:space:]]+e[0-9]{3,}.*$", "", z, perl = TRUE, ignore.case = TRUE)
    z <- trimws(gsub("[.,;:/／；：，]+$", "", z, perl = TRUE))
    if (!nzchar(z) || .ref_is_doi_fragment(z) || .ref_is_access_fragment(z)) return(NA_character_)
    z
  }
  split_authors <- function(x) {
    z <- .ref_norm_spaces2(x)
    z <- gsub("\\s+(and|&)\\s+", ", ", z, ignore.case = TRUE, perl = TRUE)
    z <- gsub("\\bet\\s+al\\.?", "", z, ignore.case = TRUE, perl = TRUE)
    z <- gsub("\\.{2,}|…", "", z, perl = TRUE)
    z <- gsub("[。；;、]", ",", z, perl = TRUE)
    parts <- trimws(unlist(strsplit(z, "\\s*,\\s*", perl = TRUE), use.names = FALSE))
    parts <- parts[nzchar(parts)]
    parts <- parts[grepl("[A-Za-z]|[\u4e00-\u9fff]", parts, perl = TRUE)]
    parts <- parts[!grepl("^(and|et al)$", parts, ignore.case = TRUE, perl = TRUE)]
    unique(parts)
  }

  rows <- list(); i <- 1L
  while (i <= length(lines) - 3L) {
    if (is_author_line(lines[i + 1L]) && is_year_line(lines[i + 3L])) {
      journal <- clean_journal(lines[i + 2L])
      authors <- split_authors(lines[i + 1L])
      if (length(authors) || (!is.na(journal) && nzchar(journal))) {
        rows[[length(rows) + 1L]] <- list(journal = ifelse(is.na(journal), "", journal), authors = authors)
      }
      i <- i + 4L
    } else {
      i <- i + 1L
    }
  }
  if (length(rows) < 2L) return(data.frame())
  max_authors <- max(1L, max(vapply(rows, function(z) length(z$authors), integer(1))))
  out <- data.frame(Journal = vapply(rows, function(z) z$journal, character(1)), stringsAsFactors = FALSE, check.names = FALSE)
  for (k in seq_len(max_authors)) {
    out[[paste0("Author_", k)]] <- vapply(rows, function(z) if (length(z$authors) >= k) z$authors[k] else "", character(1))
  }

  out[is.na(out)] <- ""
  out
}

# ---- Google Scholar article-citation parser + h-index/i10 summary ----
# This fixes journal-only extraction when pasted Google Scholar rows are used.
# It keeps the source/journal row as Journal and uses the citation/year row
# to compute a Google-Scholar-like citation table.
.ref_gs_parse_cite_year <- function(x) {
  z <- .ref_norm_spaces2(x)
  z <- gsub("\t+", " ", z, perl = TRUE)
  z <- gsub("[,，]", "", z, perl = TRUE)
  z <- gsub("\\s+", " ", z, perl = TRUE)
  z <- trimws(z)
  m <- regexec("^(?:([0-9]+)\\s+)?((?:19|20)[0-9]{2})$", z, perl = TRUE)
  r <- regmatches(z, m)[[1]]
  if (length(r) >= 3) {
    cites <- ifelse(nzchar(r[2]), suppressWarnings(as.integer(r[2])), 0L)
    year  <- suppressWarnings(as.integer(r[3]))
    if (!is.finite(cites)) cites <- 0L
    if (!is.finite(year)) year <- NA_integer_
    return(list(citations = cites, year = year))
  }
  NULL
}

.ref_gs_is_author_line_relaxed <- function(x) {
  z <- .ref_norm_spaces2(x)
  if (!nzchar(z)) return(FALSE)
  if (!is.null(.ref_gs_parse_cite_year(z))) return(FALSE)
  if (grepl("^(Title|Authors|Publication|Journal|Source|Cited by|Year|All|Since|Sort by|Show more|顯示更多|匯出|Export)$", z, ignore.case = TRUE, perl = TRUE)) return(FALSE)
  if (grepl("https?://|doi|PMID|PMCID", z, ignore.case = TRUE, perl = TRUE)) return(FALSE)
  latin <- gregexpr("\\b[A-Z]{1,6}[[:space:]]+[A-Z][A-Za-z'\\-]+\\b", z, perl = TRUE)[[1]]
  n_latin <- if (length(latin) == 1L && latin[1] == -1L) 0L else length(latin)
  # Google Scholar bylines are usually comma-delimited; Chinese bylines may use 、/，/；.
  has_delim <- grepl(",|，|、|；|;| and | & ", z, ignore.case = TRUE, perl = TRUE)
  has_zh_name <- grepl("[\u4e00-\u9fff]{2,4}(?:[,，、；;]|$)", z, perl = TRUE)
  n_latin >= 1L || has_zh_name || has_delim
}

.ref_google_scholar_articles <- function(txt) {
  empty <- data.frame(
    reference_no = integer(0), title = character(0), authors = character(0),
    journal = character(0), citations = integer(0), year = integer(0),
    source = character(0), normalized_record = logical(0), record_style = character(0),
    stringsAsFactors = FALSE
  )
  txt <- enc2utf8(as.character(txt %||% ""))
  txt <- gsub("\r\n|\r", "\n", txt, perl = TRUE)
  lines <- trimws(unlist(strsplit(txt, "\n", fixed = TRUE), use.names = FALSE))
  lines <- .ref_norm_spaces2(lines)
  lines <- lines[nzchar(lines)]
  if (length(lines) < 4L) return(empty)

  bad <- grepl("^(Title|Authors|Publication|Journal|Source|Cited by|Year|All|Since|Sort by|Show more|顯示更多|匯出|Export|文章|引用次數|共同作者)$", lines,
               ignore.case = TRUE, perl = TRUE)
  lines <- lines[!bad]
  if (length(lines) < 4L) return(empty)

  rows <- list(); i <- 1L
  while (i <= length(lines) - 3L) {
    cy <- .ref_gs_parse_cite_year(lines[i + 3L])
    if (.ref_gs_is_author_line_relaxed(lines[i + 1L]) && !is.null(cy)) {
      title <- lines[i]
      authors <- lines[i + 1L]
      journal_line <- lines[i + 2L]
      journal <- .ref_gs_clean_source_to_journal(journal_line)
      bad_title <- .ref_ncku_pure_bad_title(title) || grepl("^(查看|新增|正在追蹤|個別顯示|讓同事|文章|引用次數|共同作者|標題|隱私權|服務條款|說明)$", title, perl = TRUE)
      bad_source <- grepl("^(查看|新增|正在追蹤|文章|引用次數|共同作者|標題|隱私權|服務條款|說明)$", journal_line, perl = TRUE)
      if (!bad_title && !bad_source && !is.na(journal) && nzchar(journal) && !.ref_title_like_journal_fragment(journal)) {
        rows[[length(rows) + 1L]] <- data.frame(
          reference_no = length(rows) + 1L,
          title = title,
          authors = authors,
          journal = journal,
          citations = as.integer(cy$citations),
          year = as.integer(cy$year),
          source = "Google Scholar",
          normalized_record = TRUE,
          record_style = "Google Scholar 4-line normalized record",
          stringsAsFactors = FALSE
        )
      }
      i <- i + 4L
    } else {
      i <- i + 1L
    }
  }
  if (!length(rows)) return(empty)
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

.ref_google_scholar_articles_to_wide <- function(articles) {
  if (is.null(articles) || !is.data.frame(articles) || !nrow(articles)) return(data.frame())
  rows <- lapply(seq_len(nrow(articles)), function(i) {
    authors <- .ref_gs_split_authors(articles$authors[i])
    list(journal = as.character(articles$journal[i]), authors = authors)
  })
  max_authors <- max(1L, max(vapply(rows, function(z) length(z$authors), integer(1))))
  out <- data.frame(Journal = vapply(rows, function(z) z$journal, character(1)), stringsAsFactors = FALSE, check.names = FALSE)
  for (k in seq_len(max_authors)) {
    out[[paste0("Author_", k)]] <- vapply(rows, function(z) if (length(z$authors) >= k) z$authors[k] else "", character(1))
  }
  out[is.na(out)] <- ""
  out
}

.ref_gs_h_index <- function(citations) {
  x <- suppressWarnings(as.integer(citations))
  x <- sort(x[is.finite(x) & !is.na(x)], decreasing = TRUE)
  if (!length(x)) return(0L)
  ok <- which(x >= seq_along(x))
  if (!length(ok)) 0L else max(ok)
}


# ---- Google Scholar profile-summary parser ----
# Supports copied profile summary lines such as:
#   全部   自 2021 年
#   引文   3259   2396
#   H 指數 30     25
#   i10 指數 103  87
# or English UI labels: All / Since 2021 / Citations / h-index / i10-index.
.ref_google_scholar_profile_summary <- function(txt) {
  empty <- data.frame(
    Metric = character(0), All = integer(0), Since = integer(0),
    Since_Year = integer(0), Source = character(0),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  txt <- enc2utf8(as.character(txt %||% ""))
  txt <- gsub("\r\n|\r", "\n", txt, perl = TRUE)
  lines <- trimws(unlist(strsplit(txt, "\n", fixed = TRUE), use.names = FALSE))
  lines <- .ref_norm_spaces2(lines)
  lines <- lines[nzchar(lines)]
  if (!length(lines)) return(empty)

  # Since year from header, if available.
  sy <- NA_integer_
  hdr <- lines[grepl("自[[:space:]]*(?:19|20)[0-9]{2}|Since[[:space:]]*(?:19|20)[0-9]{2}", lines, ignore.case = TRUE, perl = TRUE)]
  if (length(hdr)) {
    sy <- suppressWarnings(as.integer(sub("^.*?((?:19|20)[0-9]{2}).*$", "\\1", hdr[1], perl = TRUE)))
    if (!is.finite(sy)) sy <- NA_integer_
  }

  pick_nums <- function(pattern) {
    hit <- lines[grepl(pattern, lines, ignore.case = TRUE, perl = TRUE)]
    if (!length(hit)) return(NULL)
    # Prefer a line with at least one number and exclude article rows by using labels.
    z <- hit[1]
    nums <- regmatches(z, gregexpr("[0-9][0-9,]*", z, perl = TRUE))[[1]]
    nums <- suppressWarnings(as.integer(gsub(",", "", nums, perl = TRUE)))
    nums <- nums[is.finite(nums)]
    if (!length(nums)) return(NULL)
    # If the metric line itself contains a since-year, remove that year from metric values.
    if (length(nums) >= 3 && any(nums %in% sy)) nums <- nums[!(nums %in% sy)]
    nums
  }

  cit <- pick_nums("^(引文|Citations)")
  hix <- pick_nums("^(H[[:space:]]*指數|h[ -]?index)")
  i10 <- pick_nums("^(i10[[:space:]]*指數|i10[ -]?index)")
  if (is.null(cit) && is.null(hix) && is.null(i10)) return(empty)

  metric_row <- function(label, vals) {
    if (is.null(vals) || !length(vals)) vals <- c(NA_integer_, NA_integer_)
    data.frame(
      Metric = label,
      All = as.integer(vals[1] %||% NA_integer_),
      Since = as.integer(vals[2] %||% NA_integer_),
      Since_Year = as.integer(ifelse(is.finite(sy), sy, NA_integer_)),
      Source = "Google Scholar profile summary",
      stringsAsFactors = FALSE, check.names = FALSE
    )
  }
  out <- rbind(
    metric_row("Google Scholar profile citations", cit),
    metric_row("Google Scholar profile h-index", hix),
    metric_row("Google Scholar profile i10-index", i10)
  )
  out
}

# Build a clear combined h-index table: profile-reported metrics first (when present),
# then metrics computed from normalized parsed reference rows.
.ref_reference_hindex_table <- function(articles, profile_summary = NULL) {
  rows <- list()
  if (is.data.frame(profile_summary) && nrow(profile_summary)) {
    ps <- profile_summary
    # One row per profile metric, plus Since column when known.
    for (i in seq_len(nrow(ps))) {
      m <- as.character(ps$Metric[i])
      rows[[length(rows)+1L]] <- data.frame(
        Metric = m,
        Value = as.character(ps$All[i]),
        Source = "profile-reported",
        stringsAsFactors = FALSE, check.names = FALSE
      )
      if ("Since" %in% names(ps) && is.finite(suppressWarnings(as.numeric(ps$Since[i])))) {
        sy <- suppressWarnings(as.integer(ps$Since_Year[i])); sy_txt <- ifelse(is.finite(sy), paste0(" since ", sy), " since period")
        rows[[length(rows)+1L]] <- data.frame(
          Metric = paste0(m, sy_txt),
          Value = as.character(ps$Since[i]),
          Source = "profile-reported",
          stringsAsFactors = FALSE, check.names = FALSE
        )
      }
    }
  }

  if (is.null(articles) || !is.data.frame(articles) || !nrow(articles)) {
    if (!length(rows)) {
      return(data.frame(Metric = "Google Scholar/NCKU normalized article-citation rows", Value = "Not detected", Source = "parsed normalized references", stringsAsFactors = FALSE, check.names = FALSE))
    }
  } else {
    cites <- suppressWarnings(as.integer(articles$citations))
    cites[!is.finite(cites) | is.na(cites)] <- 0L
    rows[[length(rows)+1L]] <- data.frame(Metric="Parsed normalized references", Value=as.character(length(cites)), Source="parsed normalized references", stringsAsFactors=FALSE, check.names=FALSE)
    rows[[length(rows)+1L]] <- data.frame(Metric="Parsed normalized-reference citations", Value=as.character(sum(cites, na.rm=TRUE)), Source="parsed normalized references", stringsAsFactors=FALSE, check.names=FALSE)
    rows[[length(rows)+1L]] <- data.frame(Metric="Parsed normalized-reference h-index", Value=as.character(.ref_gs_h_index(cites)), Source="parsed normalized references", stringsAsFactors=FALSE, check.names=FALSE)
    rows[[length(rows)+1L]] <- data.frame(Metric="Parsed normalized-reference i10-index", Value=as.character(sum(cites >= 10L, na.rm=TRUE)), Source="parsed normalized references", stringsAsFactors=FALSE, check.names=FALSE)
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

.ref_reference_year_counts <- function(articles) {
  if (is.null(articles) || !is.data.frame(articles) || !nrow(articles) || !("year" %in% names(articles))) {
    return(data.frame(Year=integer(0), Articles=integer(0), Citations=integer(0), stringsAsFactors=FALSE))
  }
  yy <- suppressWarnings(as.integer(articles$year))
  cc <- suppressWarnings(as.integer(articles$citations)); cc[!is.finite(cc) | is.na(cc)] <- 0L
  ok <- is.finite(yy) & !is.na(yy)
  if (!any(ok)) return(data.frame(Year=integer(0), Articles=integer(0), Citations=integer(0), stringsAsFactors=FALSE))
  d <- data.frame(Year=yy[ok], Citations=cc[ok], stringsAsFactors=FALSE)
  out <- aggregate(list(Articles=rep(1L, nrow(d)), Citations=d$Citations), by=list(Year=d$Year), FUN=sum)
  out <- out[order(out$Year), , drop=FALSE]
  rownames(out) <- NULL
  out
}

.ref_plot_reference_year_bar <- function(articles, main = "Parsed references by publication year") {
  yc <- .ref_reference_year_counts(articles)
  if (!is.data.frame(yc) || !nrow(yc)) {
    plot.new(); text(0.5, 0.5, "No normalized parsed reference years available."); return(invisible(NULL))
  }
  op <- par(mar = c(5, 5, 4, 2)); on.exit(par(op), add = TRUE)
  barplot(height = yc$Articles, names.arg = yc$Year, las = 2,
          xlab = "Publication year", ylab = "Number of normalized references",
          main = main)
  invisible(yc)
}

.ref_google_scholar_metrics <- function(articles, since_year = NULL) {
  if (is.null(since_year) || !is.finite(since_year)) since_year <- as.integer(format(Sys.Date(), "%Y")) - 5L
  if (is.null(articles) || !is.data.frame(articles) || !nrow(articles)) {
    return(data.frame(metric = "Google Scholar/NCKU article-citation rows", All = "Not detected", check.names = FALSE))
  }
  cites <- suppressWarnings(as.integer(articles$citations)); cites[!is.finite(cites) | is.na(cites)] <- 0L
  yrs <- suppressWarnings(as.integer(articles$year))
  keep_since <- is.finite(yrs) & !is.na(yrs) & yrs >= since_year
  cites_since <- cites[keep_since]
  out <- data.frame(
    metric = c("Citations", "h-index", "i10-index"),
    All = c(sum(cites), .ref_gs_h_index(cites), sum(cites >= 10L)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  out[[paste0("Since ", since_year, "*")]] <- c(sum(cites_since), .ref_gs_h_index(cites_since), sum(cites_since >= 10L))
  out
}


# ---- All-reference Google Scholar/NCKU h-index summary ----
# This is article/reference-level only: each parsed reference contributes its
# own citation count once. It is independent of the first/last-author and
# all-author h-index tables.
.ref_google_scholar_all_reference_hindex <- function(articles, profile_summary = NULL) {
  .ref_reference_hindex_table(articles, profile_summary = profile_summary)
}


# ---- Google Scholar author-based h-index tables ----
# basis = "first_last": each article contributes only its first and last byline authors.
# basis = "all": each article contributes every cleaned byline author.
.ref_google_scholar_author_hindex <- function(articles, basis = c("first_last", "all"), since_year = NULL) {
  basis <- match.arg(basis)
  if (is.null(since_year) || !is.finite(since_year)) since_year <- as.integer(format(Sys.Date(), "%Y")) - 5L
  empty <- data.frame(
    author = character(0), papers = integer(0), citations = integer(0),
    h_index = integer(0), i10_index = integer(0),
    first_author_papers = integer(0), last_author_papers = integer(0),
    since_year = integer(0), papers_since = integer(0), citations_since = integer(0),
    h_index_since = integer(0), i10_index_since = integer(0),
    basis = character(0), stringsAsFactors = FALSE, check.names = FALSE
  )
  if (is.null(articles) || !is.data.frame(articles) || !nrow(articles)) return(empty)

  rows <- list()
  for (i in seq_len(nrow(articles))) {
    aa <- .ref_gs_split_authors(articles$authors[i] %||% "")
    aa <- aa[!is.na(aa) & nzchar(aa)]
    aa <- unique(aa)
    if (!length(aa)) next

    first_a <- aa[1]
    last_a  <- aa[length(aa)]
    use_a <- if (basis == "first_last") unique(c(first_a, last_a)) else aa
    cites <- suppressWarnings(as.integer(articles$citations[i])); if (!is.finite(cites) || is.na(cites)) cites <- 0L
    yr <- suppressWarnings(as.integer(articles$year[i])); if (!is.finite(yr) || is.na(yr)) yr <- NA_integer_
    refno <- suppressWarnings(as.integer(articles$reference_no[i])); if (!is.finite(refno) || is.na(refno)) refno <- i

    for (au in use_a) {
      rows[[length(rows) + 1L]] <- data.frame(
        author = au, reference_no = refno, citations = cites, year = yr,
        first_author = identical(au, first_a), last_author = identical(au, last_a),
        stringsAsFactors = FALSE, check.names = FALSE
      )
    }
  }
  if (!length(rows)) return(empty)
  long <- do.call(rbind, rows)
  long <- long[!duplicated(long[, c("author", "reference_no")]), , drop = FALSE]

  split_by_author <- split(long, long$author)
  out <- lapply(split_by_author, function(d) {
    cites <- suppressWarnings(as.integer(d$citations)); cites[!is.finite(cites) | is.na(cites)] <- 0L
    yrs <- suppressWarnings(as.integer(d$year))
    keep_since <- is.finite(yrs) & !is.na(yrs) & yrs >= since_year
    cites_since <- cites[keep_since]
    data.frame(
      author = as.character(d$author[1]),
      papers = length(unique(d$reference_no)),
      citations = sum(cites, na.rm = TRUE),
      h_index = .ref_gs_h_index(cites),
      i10_index = sum(cites >= 10L, na.rm = TRUE),
      first_author_papers = sum(isTRUE(d$first_author) | d$first_author, na.rm = TRUE),
      last_author_papers = sum(isTRUE(d$last_author) | d$last_author, na.rm = TRUE),
      since_year = since_year,
      papers_since = sum(keep_since, na.rm = TRUE),
      citations_since = sum(cites_since, na.rm = TRUE),
      h_index_since = .ref_gs_h_index(cites_since),
      i10_index_since = sum(cites_since >= 10L, na.rm = TRUE),
      basis = ifelse(basis == "first_last", "1st/last authors only", "all byline authors"),
      stringsAsFactors = FALSE, check.names = FALSE
    )
  })
  out <- do.call(rbind, out)
  out <- out[order(-out$h_index, -out$citations, -out$papers, out$author), , drop = FALSE]
  rownames(out) <- NULL
  out
}



# ---- NCKU Pure / researchoutput.ncku.edu.tw Scopus pasted-row parser ----
# Supports rows copied from the NCKU Research Output / Pure page, for example:
#   Title
#   Lee, H. F., Chiang, H. Y. & Kuo, H. T., 2019 1月, 於: Journal of Nursing Management. 27, 1, p. 52-65
#   研究成果: Article › 同行評審
#   ...
#   121
#     連結會在新分頁中打開
#   引文
#   斯高帕斯（Scopus）
# A copied normalized row may also provide only a numeric citation count after
# each record; that count is used for h-index computation when unambiguous.
# The output intentionally matches .ref_google_scholar_articles(), so the same
# h-index, i10-index, journal-frequency, first/last-author, and all-author
# reports can be reused without changing downstream plotting code.
.ref_ncku_pure_empty_articles <- function() {
  data.frame(
    reference_no = integer(0), title = character(0), authors = character(0),
    journal = character(0), citations = integer(0), year = integer(0),
    source = character(0), normalized_record = logical(0), record_style = character(0),
    stringsAsFactors = FALSE, check.names = FALSE
  )
}

.ref_ncku_pure_initial_token <- function(x) {
  z <- .ref_norm_spaces2(x)
  if (!nzchar(z)) return(FALSE)
  zz <- gsub("[^A-Za-z]", "", z, perl = TRUE)
  nzchar(zz) && nchar(zz) <= 8L && grepl("^[A-Z]+$", zz, perl = TRUE)
}

.ref_ncku_pure_author_one <- function(surname, initials = "") {
  s <- .ref_norm_spaces2(surname)
  s <- gsub("\\([^)]*(Editor|編輯)[^)]*\\)", "", s, ignore.case = TRUE, perl = TRUE)
  # Prefer English transliteration inside parentheses when available.
  m <- regexec("\\(([A-Za-z][A-Za-z .'-]+)\\)", s, perl = TRUE)
  r <- regmatches(s, m)[[1]]
  if (length(r) >= 2L && nzchar(r[2])) s <- r[2]
  s <- gsub("^[[:punct:] ]+|[[:punct:] ]+$", "", s, perl = TRUE)
  s <- .ref_norm_spaces2(s)
  if (!nzchar(s)) return(NA_character_)

  ii <- toupper(gsub("[^A-Za-z]", "", .ref_norm_spaces2(initials), perl = TRUE))
  if (nzchar(ii) && grepl("^[A-Za-z][A-Za-z' -]+$", s, perl = TRUE)) {
    return(.ref_clean_author_strict(paste(ii, s)))
  }

  # Already-initialized forms such as "Lee H. F." or "HF Lee".
  m2 <- regexec("^([A-Z][A-Za-z'\\-]+)\\s+((?:[A-Z]\\.?[- ]*){1,6})$", s, perl = TRUE)
  r2 <- regmatches(s, m2)[[1]]
  if (length(r2) >= 3L) {
    ii2 <- toupper(gsub("[^A-Za-z]", "", r2[3], perl = TRUE))
    return(.ref_clean_author_strict(paste(ii2, r2[2])))
  }

  .ref_clean_author_strict(s)
}

.ref_ncku_pure_split_authors <- function(x) {
  z <- .ref_norm_spaces2(x)
  if (!nzchar(z)) return(character(0))
  z <- gsub("\\bet\\s+al\\.?", "", z, ignore.case = TRUE, perl = TRUE)
  z <- gsub("\\s+(&|and)\\s+", ", ", z, ignore.case = TRUE, perl = TRUE)
  z <- gsub("[；;、]", ",", z, perl = TRUE)
  toks <- trimws(unlist(strsplit(z, "\\s*,\\s*", perl = TRUE), use.names = FALSE))
  toks <- toks[nzchar(toks)]
  if (!length(toks)) return(character(0))

  out <- character(0)
  i <- 1L
  while (i <= length(toks)) {
    if (i < length(toks) && .ref_ncku_pure_initial_token(toks[i + 1L])) {
      out <- c(out, .ref_ncku_pure_author_one(toks[i], toks[i + 1L]))
      i <- i + 2L
    } else {
      out <- c(out, .ref_ncku_pure_author_one(toks[i], ""))
      i <- i + 1L
    }
  }
  out <- unique(out[!is.na(out) & nzchar(out)])
  out
}

.ref_ncku_pure_journal_from_source <- function(src) {
  z <- .ref_norm_spaces2(src)
  if (!grepl("於\\s*:", z, perl = TRUE)) return(NA_character_)
  j <- sub("^.*?於\\s*:\\s*", "", z, perl = TRUE)
  # Remove volume/issue/page/article-number suffix after the journal title.
  j <- sub("\\.\\s*(?:[0-9]+|p\\.|e[0-9]|[0-9]{4}|[A-Za-z]*[[:space:]]*[0-9]|\\(Accepted|卷|頁).*$", "", j, ignore.case = TRUE, perl = TRUE)
  j <- sub("[[:space:]]+(?:[0-9]+|p\\.|e[0-9]).*$", "", j, ignore.case = TRUE, perl = TRUE)
  j <- .ref_clean_journal_strict(j)
  if (length(j) && !is.na(j) && nzchar(j)) as.character(j) else NA_character_
}

.ref_ncku_pure_parse_source_line <- function(src) {
  z <- .ref_norm_spaces2(src)
  y <- suppressWarnings(as.integer(sub("^.*?((?:19|20)[0-9]{2}).*$", "\\1", z, perl = TRUE)))
  if (!is.finite(y)) y <- NA_integer_
  byline <- sub("\\s*,?\\s*(?:19|20)[0-9]{2}.*$", "", z, perl = TRUE)
  journal <- .ref_ncku_pure_journal_from_source(z)
  list(authors = .ref_ncku_pure_split_authors(byline), journal = journal, year = y)
}

.ref_ncku_pure_bad_title <- function(x) {
  z <- .ref_norm_spaces2(x)
  !nzchar(z) ||
    grepl("^(搜尋結果|國立成功大學|在 國立成功大學|首頁|概要|研究單位|研究成果|專案|學生論文|設備|獎項|活動|更多|概覽|指紋|網路|開啟存取|引文|斯高帕斯|技術支援|網站使用|登入 Pure|關於無障礙|通報網頁|聯繫我們)$", z, perl = TRUE) ||
    grepl("^(Article|Conference contribution|Chapter|Review article|Book)$", z, ignore.case = TRUE, perl = TRUE) ||
    grepl("^(?:19|20)[0-9]{2}$", z, perl = TRUE) ||
    grepl("^[0-9]+%$", z, perl = TRUE) ||
    grepl("^Article has an altmetric score", z, ignore.case = TRUE, perl = TRUE) ||
    grepl("連結會在新分頁中打開", z, perl = TRUE) ||
    grepl("^研究成果:", z, perl = TRUE) ||
    grepl("^貢獻的翻譯標題", z, perl = TRUE)
}

.ref_ncku_pure_extract_citations <- function(lines, scan_rng) {
  # NCKU/Pure pages may place the Scopus citation count after each record as:
  #   121 / 引文 / 斯高帕斯（Scopus）
  # or as a compact line such as "Cited by 121" / "121 citations".
  # Use explicit labels first; if no label is present, use a single numeric-only
  # line after the source row as the record-level citation count.  This allows
  # the parsed reference citations to support all-reference and first/last-author
  # h-index computations.
  if (!length(scan_rng)) return(0L)
  ln <- .ref_norm_spaces2(lines[scan_rng])
  ln[is.na(ln)] <- ""
  local_num <- function(x) {
    x <- gsub("[,，]", "", .ref_norm_spaces2(x), perl = TRUE)
    z <- regmatches(x, regexpr("[0-9][0-9]*", x, perl = TRUE))
    if (!length(z) || is.na(z) || !nzchar(z)) return(NA_integer_)
    suppressWarnings(as.integer(z))
  }

  # Same-line labelled forms.
  lab_pat <- "(引文|斯高帕斯|Scopus|Cited by|Citations|Citation count|被引用|引用)"
  lab_idx <- which(grepl(lab_pat, ln, ignore.case = TRUE, perl = TRUE) & grepl("[0-9]", ln, perl = TRUE))
  if (length(lab_idx)) {
    val <- local_num(ln[lab_idx[1]])
    if (is.finite(val) && !is.na(val)) return(as.integer(val))
  }

  # Numeric-only line near labelled lines, e.g. 121 / link-open text / 引文 / Scopus.
  num_idx <- which(grepl("^[0-9][0-9,]*$", ln, perl = TRUE))
  if (length(num_idx)) {
    for (ii in num_idx) {
      lo <- max(1L, ii - 3L)
      hi <- min(length(ln), ii + 8L)
      look <- ln[seq.int(lo, hi)]
      if (any(grepl(lab_pat, look, ignore.case = TRUE, perl = TRUE))) {
        val <- local_num(ln[ii])
        if (is.finite(val) && !is.na(val)) return(as.integer(val))
      }
    }
  }

  # Fallback: copied normalized records sometimes have exactly one number after
  # each reference but omit the visible label.  Treat that number as citations.
  # Avoid years, percentages, and page-like long fragments.
  if (length(num_idx) == 1L) {
    val <- local_num(ln[num_idx[1]])
    if (is.finite(val) && !is.na(val) && val >= 0L) return(as.integer(val))
  }
  0L
}

.ref_ncku_pure_articles <- function(txt) {
  empty <- .ref_ncku_pure_empty_articles()
  txt <- enc2utf8(as.character(txt %||% ""))
  txt <- gsub("\r\n|\r", "\n", txt, perl = TRUE)
  if (!grepl("研究成果:|斯高帕斯|Scopus|Pure|researchoutput\\.ncku|於\\s*:", txt, ignore.case = TRUE, perl = TRUE)) return(empty)

  lines <- trimws(unlist(strsplit(txt, "\n", fixed = TRUE), use.names = FALSE))
  lines <- .ref_norm_spaces2(lines)
  lines <- lines[nzchar(lines)]
  if (length(lines) < 3L) return(empty)

  # Flexible normalized-record detector.  NCKU/Pure copied rows may keep
  # "authors + year + 於: journal" on one line, or may split the author/year
  # part and the "於: journal" part across adjacent lines.  Keep the first
  # line index so the title can still be found immediately above it.
  src_rows <- list()
  for (ii in seq_along(lines)) {
    cand <- c(lines[ii])
    if (ii + 1L <= length(lines)) cand <- c(cand, paste(lines[ii], lines[ii + 1L]))
    if (ii + 2L <= length(lines)) cand <- c(cand, paste(lines[ii], lines[ii + 1L], lines[ii + 2L]))
    hit <- cand[grepl("(?:19|20)[0-9]{2}.*於\\s*:", cand, perl = TRUE)]
    if (length(hit)) {
      src_rows[[length(src_rows) + 1L]] <- data.frame(si = ii, src = hit[1], stringsAsFactors = FALSE)
    }
  }
  if (!length(src_rows)) return(empty)
  src_df <- do.call(rbind, src_rows)
  src_df <- src_df[!duplicated(paste(src_df$si, src_df$src, sep = "\r")), , drop = FALSE]

  rows <- list()
  for (jj in seq_len(nrow(src_df))) {
    si <- src_df$si[jj]
    src <- src_df$src[jj]
    meta <- .ref_ncku_pure_parse_source_line(src)
    if (is.na(meta$journal) || !nzchar(meta$journal) || !length(meta$authors)) next

    # Prefer the closest previous valid line as the title.  This is more robust
    # than requiring si - 1, because Pure pages sometimes insert translated-title
    # or open-access/status rows between title and source metadata.
    title <- NA_character_
    prev_rng <- if (si > 1L) rev(seq.int(max(1L, si - 6L), si - 1L)) else integer(0)
    for (ti in prev_rng) {
      if (!.ref_ncku_pure_bad_title(lines[ti]) &&
          !grepl("(?:19|20)[0-9]{2}.*於\\s*:", lines[ti], perl = TRUE)) {
        title <- lines[ti]
        break
      }
    }
    if (is.na(title) || !nzchar(title)) next

    # A research-output type line strengthens confidence, but should not be
    # mandatory for normalized copied Scopus/Pure rows.  Citation extraction is
    # bounded by the next source line to avoid leaking counts from the next item.
    next_si <- if (jj < nrow(src_df)) src_df$si[jj + 1L] else (length(lines) + 1L)
    scan_to <- max(si + 1L, next_si - 1L)
    scan_rng <- if (si + 1L <= scan_to) seq.int(si + 1L, scan_to) else integer(0)
    cites <- .ref_ncku_pure_extract_citations(lines, scan_rng)

    rows[[length(rows) + 1L]] <- data.frame(
      reference_no = length(rows) + 1L,
      title = title,
      authors = paste(meta$authors, collapse = ", "),
      journal = as.character(meta$journal),
      citations = as.integer(cites),
      year = as.integer(meta$year),
      source = "NCKU Pure / Scopus",
      normalized_record = TRUE,
      record_style = "NCKU Pure/Scopus normalized record",
      stringsAsFactors = FALSE, check.names = FALSE
    )
  }
  if (!length(rows)) return(empty)
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

.ref_reference_citation_articles <- function(txt, sources = c("ama", "google", "scopus")) {
  txt0 <- enc2utf8(as.character(txt %||% ""))
  sources <- unique(tolower(as.character(sources %||% c("ama", "google", "scopus"))))
  if (!length(sources)) sources <- c("ama", "google", "scopus")
  use_google <- any(sources %in% c("google", "gs", "scholar", "google_scholar"))
  use_scopus <- any(sources %in% c("scopus", "ncku", "pure", "ncku_scopus"))

  rows <- list()
  if (isTRUE(use_scopus)) {
    nk <- tryCatch(.ref_ncku_pure_articles(txt0), error = function(e) data.frame())
    if (is.data.frame(nk) && nrow(nk) > 0) rows[[length(rows) + 1L]] <- nk
  }
  if (isTRUE(use_google)) {
    gs <- tryCatch(.ref_google_scholar_articles(txt0), error = function(e) data.frame())
    if (is.data.frame(gs) && nrow(gs) > 0) {
      if (!"source" %in% names(gs)) gs$source <- "Google Scholar"
      rows[[length(rows) + 1L]] <- gs
    }
  }
  if (!length(rows)) return(.ref_ncku_pure_empty_articles())

  out <- do.call(rbind, rows)
  # De-duplicate only exact normalized records. This prevents copied page headers
  # or repeated pages from inflating citations/h-index.
  key <- paste(tolower(.ref_norm_spaces2(out$title)),
               tolower(.ref_norm_spaces2(out$authors)),
               suppressWarnings(as.integer(out$year)),
               sep = "\r")
  out <- out[!duplicated(key), , drop = FALSE]
  out$reference_no <- seq_len(nrow(out))
  rownames(out) <- NULL
  out
}

.ref_empty_author_hindex_message <- function() {
  data.frame(Message = "No Google Scholar/NCKU author-citation rows detected. Paste Google Scholar rows or NCKU Pure/Scopus rows with authors, journal/source, citations, and year.", stringsAsFactors = FALSE, check.names = FALSE)
}

.ref_to_wide_from_text_strict <- function(txt) {
  blocks <- .ref_split_pubmed_blocks_strict(txt)
  if (!length(blocks)) stop("No usable PubMed/AMA/Google Scholar/NCKU Scopus entries were detected.", call. = FALSE)
  parsed <- lapply(blocks, .ref_parse_block_strict)
  journals <- vapply(parsed, function(z) as.character(z$journal %||% ""), character(1))
  authors_list <- lapply(parsed, function(z) unique(na.omit(as.character(z$authors %||% character(0)))))
  keep <- nzchar(journals) | lengths(authors_list) > 0
  journals <- journals[keep]
  authors_list <- authors_list[keep]
  if (!length(journals)) stop("No usable journal/author terms after parsing.", call. = FALSE)

  max_authors <- max(1L, max(lengths(authors_list), na.rm = TRUE))
  out <- data.frame(Journal = journals, stringsAsFactors = FALSE, check.names = FALSE)
  for (k in seq_len(max_authors)) {
    out[[paste0("Author_", k)]] <- vapply(authors_list, function(v) if (length(v) >= k) v[k] else "", character(1))
  }
  out[is.na(out)] <- ""
  out
}

.ref_valid_node <- function(x) {
  x <- .ref_norm_spaces2(x)
  if (!nzchar(x)) return(FALSE)
  if (.ref_is_doi_fragment(x) || .ref_is_access_fragment(x) || .ref_title_like_journal_fragment(x)) return(FALSE)
  if (grepl("\\[doi\\]|\\b(PST|SO|AID|PMID|PMCID)\\s*-", x, ignore.case = TRUE, perl = TRUE)) return(FALSE)
  if (grepl("^(P|S|N)$", x, perl = TRUE)) return(FALSE)
  if (grepl("^[0-9]", x, perl = TRUE)) return(FALSE)
  if (grepl("case fatality|observational study|usability study|meta-analysis", x, ignore.case = TRUE, perl = TRUE)) return(FALSE)
  TRUE
}



# ============================================================
# SAFE Slope Combo Entity helpers
# Avoid htmltools::HTML()/renderUI() to prevent:
# "不是所有的 is.character(txt) 都是 TRUE"
# ============================================================

.slope_safe_chr <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- iconv(x, from = "", to = "UTF-8", sub = "")
  x[is.na(x)] <- ""
  trimws(x)
}

.slope_pick_year_col <- function(df) {
  nms <- names(df)
  cand <- c("Year", "year", "PY", "Publication Year", "Publication.Year", "publication_year")
  hit <- cand[cand %in% nms]
  if (length(hit)) return(hit[1])
  hit <- nms[grepl("year|publication.*year", tolower(nms))]
  if (length(hit)) hit[1] else NA_character_
}

.slope_pick_entity_col <- function(df, domain) {
  nms <- names(df)
  domain <- as.character(domain %||% "")
  map <- list(
    Country = c("Country","country","FAcountry","CAcountry"),
    Institute = c("Institute","institute","FAinstitute","CAinstitute"),
    Department = c("Department","department","FADept","CAdept"),
    Author = c("Author","author","FAauthor","CAauthor"),
    Journal = c("Journal","journal","SO","Source Title","Source.Title"),
    Keyword = c("Keyword","keyword","KeywordsPlus","Keywords Plus","DE","ID"),
    Region = c("Region","region","FAregion","CAregion"),
    `State/Province` = c("Region","region","FAregion","CAregion"),
    DocumentType = c("DocumentType","Document Type","DT","document_type"),
    Year = c("Year","year")
  )
  cand <- map[[domain]]
  if (is.null(cand)) cand <- c(domain, tolower(domain))
  hit <- cand[cand %in% nms]
  if (length(hit)) return(hit[1])
  hit <- nms[tolower(nms) %in% tolower(cand)]
  if (length(hit)) hit[1] else NA_character_
}

.slope_split_items <- function(x) {
  x <- .slope_safe_chr(x)
  x <- unlist(strsplit(x, "\\s*;\\s*|\\s*\\|\\s*", perl = TRUE), use.names = FALSE)
  x <- .slope_safe_chr(x)
  x[nzchar(x)]
}

.slope_get_meta_df <- function(rv = NULL) {
  # Try common reactiveValues names used in app_709/app_708 without failing.
  cand <- list()
  if (!is.null(rv)) {
    for (nm in c("metatable","meta","meta_df","pubmeta","woslong32","long32","woswide32","wide32")) {
      val <- tryCatch(rv[[nm]], error=function(e) NULL)
      if (is.data.frame(val) && nrow(val)) cand[[length(cand)+1]] <- val
    }
  }
  if (length(cand)) return(cand[[1]])
  NULL
}

.slope_build_combo_entity_year <- function(df, domain = "Author", top_n = 10, recent_n_years = 10) {
  if (is.null(df) || !is.data.frame(df) || !nrow(df)) return(data.frame())
  ycol <- .slope_pick_year_col(df)
  ecol <- .slope_pick_entity_col(df, domain)
  if (is.na(ycol) || is.na(ecol) || !(ycol %in% names(df)) || !(ecol %in% names(df))) return(data.frame())

  yy <- suppressWarnings(as.integer(as.character(df[[ycol]])))
  ok <- is.finite(yy)
  df <- df[ok, , drop = FALSE]
  yy <- yy[ok]
  if (!length(yy)) return(data.frame())

  ymax <- max(yy, na.rm = TRUE)
  ymin <- ymax - as.integer(recent_n_years) + 1L
  keep <- yy >= ymin & yy <= ymax
  df <- df[keep, , drop = FALSE]
  yy <- yy[keep]
  if (!nrow(df)) return(data.frame())

  rows <- lapply(seq_len(nrow(df)), function(i) {
    items <- .slope_split_items(df[[ecol]][i])
    if (!length(items)) return(NULL)
    data.frame(Year = yy[i], Entity = items, stringsAsFactors = FALSE)
  })
  rows <- do.call(rbind, rows)
  if (is.null(rows) || !nrow(rows)) return(data.frame())

  rows <- unique(rows)
  top <- as.data.frame(sort(table(rows$Entity), decreasing = TRUE), stringsAsFactors = FALSE)
  names(top) <- c("Entity","Total")
  top <- top[seq_len(min(top_n, nrow(top))), , drop = FALSE]
  rows <- rows[rows$Entity %in% top$Entity, , drop = FALSE]

  out <- as.data.frame(table(rows$Year, rows$Entity), stringsAsFactors = FALSE)
  names(out) <- c("Year","Entity","Count")
  out$Year <- suppressWarnings(as.integer(as.character(out$Year)))
  out$Count <- suppressWarnings(as.numeric(out$Count))
  out <- out[out$Entity %in% top$Entity, , drop = FALSE]
  out <- merge(out, top, by = "Entity", all.x = TRUE)
  out <- out[order(out$Entity, out$Year), , drop = FALSE]
  rownames(out) <- NULL
  out
}

.slope_plot_combo_entity_top10 <- function(dat, domain = "Entity") {
  if (is.null(dat) || !is.data.frame(dat) || !nrow(dat)) {
    plot.new()
    text(0.5, 0.5, "No slope data available.\nRun analysis first or choose another entity/domain.", cex = 1.1)
    return(invisible(NULL))
  }
  yrs <- sort(unique(dat$Year))
  ents <- unique(dat$Entity[order(-dat$Total, dat$Entity)])
  ymax <- max(dat$Count, na.rm = TRUE)
  plot(range(yrs), c(0, ymax * 1.15), type = "n",
       xlab = "Year", ylab = "Number of papers",
       main = paste0("Top 10 ", domain, " slope graph over years"))
  for (e in ents) {
    d <- dat[dat$Entity == e, , drop = FALSE]
    d <- d[order(d$Year), , drop = FALSE]
    lines(d$Year, d$Count, lwd = 2)
    points(d$Year, d$Count, pch = 19)
    if (nrow(d)) {
      text(d$Year[nrow(d)], d$Count[nrow(d)], labels = e, pos = 4, cex = 0.75, xpd = NA)
    }
  }
}

server <- function(input, output, session) {

  output$slope_combo_entity_note <- renderText({
    "Shows the yearly slope graph for the Top 10 selected entity elements. Run analysis first, then choose a domain/entity."
  })

  output$slope_combo_entity_plot <- renderPlot({
    domain <- input$slope_combo_entity_domain %||% input$combo_domain %||% input$domain %||% "Author"
    recent_n <- input$slope_combo_recent_years %||% 10
    df <- .slope_get_meta_df(rv)
    dat <- .slope_build_combo_entity_year(df, domain = domain, top_n = 10, recent_n_years = recent_n)
    .slope_plot_combo_entity_top10(dat, domain = domain)
  })

  output$slope_combo_entity_table <- DT::renderDT({
    domain <- input$slope_combo_entity_domain %||% input$combo_domain %||% input$domain %||% "Author"
    recent_n <- input$slope_combo_recent_years %||% 10
    df <- .slope_get_meta_df(rv)
    dat <- .slope_build_combo_entity_year(df, domain = domain, top_n = 10, recent_n_years = recent_n)
    DT::datatable(dat, rownames = FALSE, options = list(pageLength = 20, scrollX = TRUE))
  })


  aac_server("aac")
  free_demo_run <- reactiveVal(TRUE)
  demo_term <- "asthma[Title/Abstract]"   # Demo 預設查詢（可改）

  ama_processing_status <- reactiveVal("Ready: upload a PubMed/AMA TXT file or paste normalized Google Scholar/NCKU Scopus rows, then click Run journal extraction.")

  output$ama_processing_status <- renderText({
    as.character(ama_processing_status() %||% "")
  })

  observeEvent(input$btn_clear_ama_refs_text, {
    updateTextAreaInput(session, "ama_refs_text", value = "")
    ama_processing_status("Textarea cleared. Uploaded TXT and PubMed query were not changed.")
    if (exists("ama_ref_status", inherits = FALSE) && is.function(ama_ref_status)) {
      ama_ref_status("Textarea cleared. Uploaded TXT and PubMed query were not changed.")
    }
    if (exists("ama_author_aac_status", inherits = FALSE) && is.function(ama_author_aac_status)) {
      ama_author_aac_status("Textarea cleared. Paste normalized records, then click Run AMA author AAC (1st+Last only).")
    }
    showNotification("Pasted-reference textarea cleared.", type = "message")
  }, ignoreInit = TRUE)

  # Store AMA results immediately when the button is clicked.
  # Do NOT use eventReactive here, because hidden/suspended outputs may prevent
  # the extraction from running until the user opens a result tab.
  ama_journal_results_val <- reactiveVal(NULL)

  # ---- Console logger for AMA extraction buttons ----
  .ama_console_log <- function(step, ..., .button = "AMA") {
    msg <- paste0(..., collapse = "")
    message(sprintf("[%s] [%s] %s%s",
                    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                    .button,
                    step,
                    if (nzchar(msg)) paste0(": ", msg) else ""))
    invisible(NULL)
  }

  .get_ama_input_text <- function(status_fun = ama_processing_status) {
    # If checked, use pasted textarea references instead of uploaded TXT.
    if (isTRUE(input$use_ama_textarea)) {
      txt_area <- paste(input$ama_refs_text %||% "", collapse = "
")
      txt_area <- enc2utf8(as.character(txt_area))
      txt_area <- gsub("\r\n|\r", "\n", txt_area, perl = TRUE)
      if (!nzchar(trimws(txt_area))) {
        status_fun("Error: textarea option is checked, but no AMA/PubMed references were pasted.")
        showNotification("Textarea option is checked, but the textarea is empty.", type = "error")
        return(NULL)
      }
      return(txt_area)
    }

    path <- NULL
    if (!is.null(input$pubmed_txt) && !is.null(input$pubmed_txt$datapath) && file.exists(input$pubmed_txt$datapath)) {
      path <- input$pubmed_txt$datapath
    } else if (!is.null(rv$uploaded_perm_path) && file.exists(rv$uploaded_perm_path)) {
      path <- rv$uploaded_perm_path
    } else if (exists("PERM_PUBMED_TXT", inherits = TRUE) && file.exists(PERM_PUBMED_TXT)) {
      path <- PERM_PUBMED_TXT
    }

    if (is.null(path) || !file.exists(path)) {
      status_fun("Error: please upload a PubMed/AMA .txt file first, or check the textarea option and paste references.")
      showNotification("Please upload a TXT file or use the textarea option.", type = "error")
      return(NULL)
    }
    .read_text_any_encoding(path)
  }

  .ama_selected_ref_sources <- function() {
    src <- character(0)
    if (isTRUE(input$ref_src_ama))    src <- c(src, "ama")
    if (isTRUE(input$ref_src_google)) src <- c(src, "google")
    if (isTRUE(input$ref_src_scopus)) src <- c(src, "scopus")
    if (!length(src)) {
      showNotification("No reference source parser was selected; all three parsers will be used.", type = "warning")
      src <- c("ama", "google", "scopus")
    }
    unique(src)
  }

  .ama_split_refs_selected <- function(txt, sources) {
    art <- tryCatch(.ref_reference_citation_articles(txt, sources = sources), error = function(e) data.frame())
    if (is.data.frame(art) && nrow(art) > 0) {
      return(paste(art$title, art$authors, art$journal, art$year, sep = "\n"))
    }
    if ("ama" %in% sources) return(.ref_split_pubmed_blocks_strict(txt))
    character(0)
  }

  .ama_parse_articles_selected <- function(txt, sources) {
    tryCatch(.ref_reference_citation_articles(txt, sources = sources), error = function(e) data.frame())
  }

  observeEvent(input$run_ama_journals, {
    .ama_console_log("Button pressed", "Run AMA journal extraction", .button = "AMA Journal")
    ama_processing_status("Processing: checking uploaded TXT file...")
    ama_journal_results_val(NULL)
    .ama_console_log("Status", "cleared previous journal results", .button = "AMA Journal")

    tryCatch({
      shiny::withProgress(message = "Running AMA journal extraction", value = 0, {
        shiny::incProgress(0.10, detail = "Checking input source")
        .ama_console_log("Input", if (isTRUE(input$use_ama_textarea)) "using pasted textarea" else "using uploaded/permanent TXT file", .button = "AMA Journal")
        ama_processing_status(if (isTRUE(input$use_ama_textarea)) "Processing: reading pasted textarea references..." else "Processing: reading uploaded TXT file...")

        txt <- .get_ama_input_text(ama_processing_status)
        if (is.null(txt) || !nzchar(trimws(txt))) stop("No uploaded or pasted AMA/PubMed/Google Scholar/NCKU Scopus text was found.", call. = FALSE)
        .ama_console_log("Input loaded", sprintf("%d characters; %d lines", nchar(txt), length(strsplit(txt, "\n", fixed = TRUE)[[1]])), .button = "AMA Journal")

        shiny::incProgress(0.15, detail = "Splitting references")
        selected_sources <- .ama_selected_ref_sources()
        .ama_console_log("Parsers", paste(selected_sources, collapse = ", "), .button = "AMA Journal")
        refs <- .ama_split_refs_selected(txt, selected_sources)
        .ama_console_log("References", sprintf("detected normalized/selected entries = %d", length(refs)), .button = "AMA Journal")
        if (!length(refs)) stop("No usable PubMed/AMA/Google Scholar/NCKU Scopus entries were detected.", call. = FALSE)

        shiny::incProgress(0.20, detail = "Extracting journal names")
        ama_processing_status("Processing: extracting journal names from AMA/PubMed/Google Scholar/NCKU Scopus text...")
        .ama_console_log("Parser", "running strict/mixed PubMed-AMA/Google Scholar/NCKU Scopus parser", .button = "AMA Journal")

        # Prefer Google Scholar/NCKU row parsing before the journal-list fallback.
        # Otherwise article titles copied from Google Scholar/NCKU pages can be mistaken for journal names.
        gs_profile_summary <- if ("google" %in% selected_sources) tryCatch(.ref_google_scholar_profile_summary(txt), error = function(e) data.frame()) else data.frame()
        gs_articles <- .ama_parse_articles_selected(txt, selected_sources)
        if (is.data.frame(gs_articles) && nrow(gs_articles) > 0) {
          .ama_console_log("Parser", sprintf("Google Scholar/NCKU article rows detected: rows = %d; unique journals = %d", nrow(gs_articles), dplyr::n_distinct(gs_articles$journal)), .button = "AMA Journal")
          wide_strict <- .ref_google_scholar_articles_to_wide(gs_articles)
          journal_df2 <- data.frame(
            reference_no = seq_len(nrow(gs_articles)),
            journal = gs_articles$journal,
            citations = gs_articles$citations,
            year = gs_articles$year,
            title = gs_articles$title,
            stringsAsFactors = FALSE,
            check.names = FALSE
          )
          journal_df2$journal <- .ref_clean_journal_strict(journal_df2$journal)
          journal_df2 <- journal_df2[!is.na(journal_df2$journal) & nzchar(journal_df2$journal), , drop = FALSE]
          .ama_console_log("Parser done", sprintf("Google Scholar/NCKU wide rows = %d; columns = %s", nrow(wide_strict), paste(names(wide_strict), collapse = ", ")), .button = "AMA Journal")
        } else {
          exact_journal_df <- if ("ama" %in% selected_sources) .ref_extract_existing_journal_records(txt) else data.frame()
          if (is.data.frame(exact_journal_df) && nrow(exact_journal_df) > 0) {
            .ama_console_log("Parser", sprintf("detected existing journal list/table; exact journal records = %d", nrow(exact_journal_df)), .button = "AMA Journal")
            journal_df2 <- exact_journal_df
            wide_strict <- data.frame(Journal = journal_df2$journal, Author_1 = "", stringsAsFactors = FALSE, check.names = FALSE)
            .ama_console_log("Parser done", sprintf("exact journal rows = %d; columns = %s", nrow(wide_strict), paste(names(wide_strict), collapse = ", ")), .button = "AMA Journal")
          } else {
            if (!("ama" %in% selected_sources)) stop("No normalized Google Scholar/NCKU rows were extracted, and AMA/PubMed parser is unchecked.", call. = FALSE)
            wide_strict <- .ref_to_wide_from_text_strict(txt)
            .ama_console_log("Parser done", sprintf("wide rows = %d; columns = %s", nrow(wide_strict), paste(names(wide_strict), collapse = ", ")), .button = "AMA Journal")

            journal_vec <- .ref_clean_journal_strict(wide_strict$Journal)
            journal_vec <- journal_vec[!is.na(journal_vec) & nzchar(journal_vec)]
            journal_df2 <- data.frame(reference_no = seq_along(journal_vec), journal = journal_vec, stringsAsFactors = FALSE)
          }
        }
        .ama_console_log("Journals", sprintf("extracted records = %d; unique journals = %d", nrow(journal_df2), dplyr::n_distinct(journal_df2$journal)), .button = "AMA Journal")

        if (!is.data.frame(journal_df2) || nrow(journal_df2) == 0) {
          stop("No journal names were extracted. Please check MEDLINE, AMA/PubMed, Google Scholar, or NCKU Pure/Scopus pasted format.", call. = FALSE)
        }

        top_df2 <- journal_df2 |>
          dplyr::count(journal, name = "frequency", sort = TRUE) |>
          dplyr::slice_head(n = input$ama_top_n %||% 20)

        gs_since_year <- as.integer(format(Sys.Date(), "%Y")) - 5L
        gs_metrics <- if (exists("gs_articles", inherits = FALSE) && is.data.frame(gs_articles) && nrow(gs_articles) > 0) .ref_google_scholar_metrics(gs_articles, since_year = gs_since_year) else data.frame(metric = "Google Scholar/NCKU article-citation rows", All = "Not detected", stringsAsFactors = FALSE, check.names = FALSE)
        gs_author_hindex_first_last <- if (exists("gs_articles", inherits = FALSE) && is.data.frame(gs_articles) && nrow(gs_articles) > 0) .ref_google_scholar_author_hindex(gs_articles, basis = "first_last", since_year = gs_since_year) else data.frame()
        gs_author_hindex_all <- if (exists("gs_articles", inherits = FALSE) && is.data.frame(gs_articles) && nrow(gs_articles) > 0) .ref_google_scholar_author_hindex(gs_articles, basis = "all", since_year = gs_since_year) else data.frame()
        res <- list(journals = journal_df2, top = top_df2, total_refs = nrow(journal_df2), gs_articles = if (exists("gs_articles", inherits = FALSE)) gs_articles else data.frame(), gs_profile_summary = if (exists("gs_profile_summary", inherits = FALSE)) gs_profile_summary else data.frame(), gs_metrics = gs_metrics, gs_author_hindex_first_last = gs_author_hindex_first_last, gs_author_hindex_all = gs_author_hindex_all, gs_since_year = gs_since_year)

        # Journal-only mode: do not build or overwrite the mixed Author/Journal network here.
        # The mixed Journal + Author co-word network, Network plot, SSplot, Kano, and cluster AAC
        # are produced only by the separate Run AMA author/journal extraction button.
        shiny::incProgress(0.10, detail = "Saving journal-only results")
        ama_journal_results_val(res)
        ama_processing_status(paste0("Completed: extracted ", nrow(res$journals),
                                     " journal records; unique journals = ",
                                     dplyr::n_distinct(res$journals$journal), "."))
        .ama_console_log("Completed", sprintf("journal records = %d; unique journals = %d; Top-N = %d", nrow(res$journals), dplyr::n_distinct(res$journals$journal), nrow(res$top)), .button = "AMA Journal")
        showNotification("AMA journal extraction completed.", type = "message")
        updateTabsetPanel(session, "main_tabs", selected = "AMA Journals")
        shiny::incProgress(0.05, detail = "Done")
      })
    }, error = function(e) {
      .ama_console_log("Error", conditionMessage(e), .button = "AMA Journal")
      ama_processing_status(paste("Error:", conditionMessage(e)))
      showNotification(paste("AMA journal extraction failed:", conditionMessage(e)), type = "error", duration = 10)
      ama_journal_results_val(NULL)
      # Keep author/journal results untouched if they already exist; do not shut the app down.
    })
  }, ignoreInit = TRUE)

  output$ama_journal_summary <- renderText({
    r <- ama_journal_results_val()
    paste0(
      "Input type: uploaded PubMed MEDLINE or AMA/PubMed reference text\n",
      "References/records with journal extracted: ", nrow(r$journals), "\n",
      "Unique journals: ", dplyr::n_distinct(r$journals$journal), "\n",
      "Top N shown: ", input$ama_top_n, "\n\n",
      "Most frequent journal: ",
      ifelse(nrow(r$top) > 0, paste0(r$top$journal[1], " (n=", r$top$frequency[1], ")"), "None")
    )
  })

  output$tbl_ama_gs_metrics <- DT::renderDT({
    r <- ama_journal_results_val()
    if (is.null(r) || is.null(r$gs_metrics) || !is.data.frame(r$gs_metrics) || !nrow(r$gs_metrics)) {
      return(DT::datatable(data.frame(metric = "Google Scholar/NCKU article-citation rows", All = "Not detected", stringsAsFactors = FALSE, check.names = FALSE), rownames = FALSE,
                           options = list(dom = "t", scrollX = TRUE)))
    }
    DT::datatable(r$gs_metrics, rownames = FALSE,
                  options = list(dom = "t", scrollX = TRUE))
  })


  output$tbl_ama_gs_all_reference_hindex <- DT::renderDT({
    r <- ama_journal_results_val()
    art <- if (!is.null(r)) r$gs_articles %||% data.frame() else data.frame()
    df <- .ref_google_scholar_all_reference_hindex(art, profile_summary = (if (!is.null(r)) r$gs_profile_summary %||% data.frame() else data.frame()))
    DT::datatable(df, rownames = FALSE,
                  options = list(dom = "t", scrollX = TRUE))
  })

  output$plt_ama_gs_year_bar <- renderPlot({
    r <- ama_journal_results_val()
    art <- if (!is.null(r)) r$gs_articles %||% data.frame() else data.frame()
    .ref_plot_reference_year_bar(art, main = "Normalized Google Scholar / NCKU Scopus references by year")
  })

  output$tbl_ama_gs_year_bar <- DT::renderDT({
    r <- ama_journal_results_val()
    art <- if (!is.null(r)) r$gs_articles %||% data.frame() else data.frame()
    df <- .ref_reference_year_counts(art)
    if (!nrow(df)) df <- data.frame(Message = "No normalized parsed reference years available.", stringsAsFactors = FALSE)
    DT::datatable(df, rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE))
  })

  output$tbl_ama_gs_author_hindex_first_last <- DT::renderDT({
    r <- ama_journal_results_val()
    df <- if (!is.null(r)) r$gs_author_hindex_first_last %||% data.frame() else data.frame()
    if (is.null(r) || !is.data.frame(df) || !nrow(df)) {
      return(DT::datatable(.ref_empty_author_hindex_message(), rownames = FALSE,
                           options = list(dom = "t", scrollX = TRUE)))
    }
    show_cols <- intersect(c("author", "papers", "citations", "h_index", "i10_index", "first_author_papers", "last_author_papers", "papers_since", "citations_since", "h_index_since", "i10_index_since", "basis"), names(df))
    DT::datatable(df[, show_cols, drop = FALSE], rownames = FALSE,
                  options = list(pageLength = 10, scrollX = TRUE))
  })

  output$tbl_ama_gs_author_hindex_all <- DT::renderDT({
    r <- ama_journal_results_val()
    df <- if (!is.null(r)) r$gs_author_hindex_all %||% data.frame() else data.frame()
    if (is.null(r) || !is.data.frame(df) || !nrow(df)) {
      return(DT::datatable(.ref_empty_author_hindex_message(), rownames = FALSE,
                           options = list(dom = "t", scrollX = TRUE)))
    }
    show_cols <- intersect(c("author", "papers", "citations", "h_index", "i10_index", "first_author_papers", "last_author_papers", "papers_since", "citations_since", "h_index_since", "i10_index_since", "basis"), names(df))
    DT::datatable(df[, show_cols, drop = FALSE], rownames = FALSE,
                  options = list(pageLength = 10, scrollX = TRUE))
  })


  output$tbl_ama_ref_all_reference_hindex <- DT::renderDT({
    r <- ama_ref_results_val()
    art <- if (!is.null(r)) r$gs_articles %||% data.frame() else data.frame()
    df <- .ref_google_scholar_all_reference_hindex(art, profile_summary = (if (!is.null(r)) r$gs_profile_summary %||% data.frame() else data.frame()))
    DT::datatable(df, rownames = FALSE,
                  options = list(dom = "t", scrollX = TRUE))
  })

  output$plt_ama_ref_year_bar <- renderPlot({
    r <- ama_ref_results_val()
    art <- if (!is.null(r)) r$gs_articles %||% data.frame() else data.frame()
    .ref_plot_reference_year_bar(art, main = "Normalized Google Scholar / NCKU Scopus references by year")
  })

  output$tbl_ama_ref_year_bar <- DT::renderDT({
    r <- ama_ref_results_val()
    art <- if (!is.null(r)) r$gs_articles %||% data.frame() else data.frame()
    df <- .ref_reference_year_counts(art)
    if (!nrow(df)) df <- data.frame(Message = "No normalized parsed reference years available.", stringsAsFactors = FALSE)
    DT::datatable(df, rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE))
  })

  output$tbl_ama_ref_author_hindex_first_last <- DT::renderDT({
    r <- ama_ref_results_val()
    art <- if (!is.null(r)) r$gs_articles %||% data.frame() else data.frame()
    sy <- if (!is.null(r)) r$gs_since_year %||% (as.integer(format(Sys.Date(), "%Y")) - 5L) else (as.integer(format(Sys.Date(), "%Y")) - 5L)
    df <- if (is.data.frame(art) && nrow(art)) .ref_google_scholar_author_hindex(art, basis = "first_last", since_year = sy) else data.frame()
    if (!is.data.frame(df) || !nrow(df)) {
      return(DT::datatable(.ref_empty_author_hindex_message(), rownames = FALSE,
                           options = list(dom = "t", scrollX = TRUE)))
    }
    show_cols <- intersect(c("author", "papers", "citations", "h_index", "i10_index", "first_author_papers", "last_author_papers", "papers_since", "citations_since", "h_index_since", "i10_index_since", "basis"), names(df))
    DT::datatable(df[, show_cols, drop = FALSE], rownames = FALSE,
                  options = list(pageLength = 10, scrollX = TRUE))
  })

  output$tbl_ama_ref_author_hindex_all <- DT::renderDT({
    r <- ama_ref_results_val()
    art <- if (!is.null(r)) r$gs_articles %||% data.frame() else data.frame()
    sy <- if (!is.null(r)) r$gs_since_year %||% (as.integer(format(Sys.Date(), "%Y")) - 5L) else (as.integer(format(Sys.Date(), "%Y")) - 5L)
    df <- if (is.data.frame(art) && nrow(art)) .ref_google_scholar_author_hindex(art, basis = "all", since_year = sy) else data.frame()
    if (!is.data.frame(df) || !nrow(df)) {
      return(DT::datatable(.ref_empty_author_hindex_message(), rownames = FALSE,
                           options = list(dom = "t", scrollX = TRUE)))
    }
    show_cols <- intersect(c("author", "papers", "citations", "h_index", "i10_index", "first_author_papers", "last_author_papers", "papers_since", "citations_since", "h_index_since", "i10_index_since", "basis"), names(df))
    DT::datatable(df[, show_cols, drop = FALSE], rownames = FALSE,
                  options = list(pageLength = 10, scrollX = TRUE))
  })

  output$tbl_ama_top_journals <- DT::renderDT({
    DT::datatable(ama_journal_results_val()$top, rownames = FALSE,
                  options = list(pageLength = 20, scrollX = TRUE))
  })

  output$tbl_ama_journal_author_top20_check <- renderPrint({
    # IMPORTANT: do not use renderText() here.
    # renderText()/cat() cannot print a data.frame/list safely and causes:
    # "'cat' cannot handle argument type 'list'".
    tryCatch({
      cat("Top 20 author/journal elements for checking\n")
      cat("This uses the strict PubMed/AMA parser and should contain only journal titles and author names.\n")
      r <- ama_ref_results_val()
      df <- .ref_top20_check_df(r)
      print(.ama_plain_table(df, max_rows = 20), row.names = FALSE, right = FALSE)
    }, error = function(e) {
      cat("Top 20 author/journal elements for checking\n")
      cat("Run AMA author/journal extraction first. ", conditionMessage(e), "\n", sep = "")
    })
  })

  output$tbl_ama_journal_refs <- DT::renderDT({
    DT::datatable(ama_journal_results_val()$journals, rownames = FALSE,
                  options = list(pageLength = 10, scrollX = TRUE))
  })

  output$tbl_ama_gs_articles <- DT::renderDT({
    r <- ama_journal_results_val()
    if (is.null(r) || is.null(r$gs_articles) || !is.data.frame(r$gs_articles) || !nrow(r$gs_articles)) {
      return(DT::datatable(data.frame(Message = "No Google Scholar/NCKU article-citation rows detected.", stringsAsFactors = FALSE), rownames = FALSE,
                           options = list(dom = "t", scrollX = TRUE)))
    }
    show_cols <- intersect(c("reference_no", "source", "record_style", "normalized_record", "title", "authors", "journal", "citations", "year"), names(r$gs_articles))
    DT::datatable(r$gs_articles[, show_cols, drop = FALSE], rownames = FALSE,
                  options = list(pageLength = 10, scrollX = TRUE))
  })

  output$plt_ama_journal_pie <- renderPlot({
    .plot_ama_journal_pie(ama_journal_results_val()$top)
  })

  output$plt_ama_journal_ss <- renderPlot({
    tryCatch({
      if (is.null(ama_journal_results_val())) stop("Click Run AMA journal extraction first.", call. = FALSE)
      if (is.null(ama_ref_results_val()) || is.null(ama_ref_results_val()$nodes) || !is.data.frame(ama_ref_results_val()$nodes)) {
        stop("No AMA Journal network nodes/edges are available yet. Run AMA journal extraction again.", call. = FALSE)
      }
      .ama705_draw_ssplot(font_scale = 1.55,
                          footer_label = "AMA Journal extraction: Top20 author/journal nodes and edges")
    }, error = function(e) {
      plot.new()
      text(0.5, 0.5, paste("SSplot error:", conditionMessage(e)), col = "red", cex = 1.1, font = 2)
    })
  }, height = 900)

  output$dl_ama_journal_csv <- downloadHandler(
    filename = function(){ paste0("top20_ama_pubmed_journals_", Sys.Date(), ".csv") },
    content = function(file){ readr::write_csv(ama_journal_results_val()$top, file) }
  )

  output$dl_ama_journal_pie <- downloadHandler(
    filename = function(){ paste0("top20_ama_pubmed_journal_pie_", Sys.Date(), ".png") },
    content = function(file){
      ggplot2::ggsave(file, plot = .plot_ama_journal_pie(ama_journal_results_val()$top),
                      width = 10, height = 8, dpi = 300)
    }
  )

  output$dl_ama_journal_ss_png <- downloadHandler(
    filename = function(){ paste0("ama_journal_top20_ssplot_", Sys.Date(), ".png") },
    content = function(file){
      grDevices::png(file, width = 3400, height = 2600, res = 260)
      on.exit(grDevices::dev.off(), add = TRUE)
      if (is.null(ama_journal_results_val())) stop("Click Run AMA journal extraction first.", call. = FALSE)
      .ama705_draw_ssplot(font_scale = 1.65,
                          footer_label = "AMA Journal extraction: Top20 author/journal nodes and edges")
    }
  )


  # ---- AMA/PubMed author + journal extraction: FLCA-MA-SIL Top20 ----
  ama_ref_status <- reactiveVal("Ready: upload a PubMed/AMA TXT file or paste normalized Google Scholar/NCKU Scopus rows, then click Run author/journal extraction.")
  ama_ref_results_val <- reactiveVal(NULL)

  output$ama_ref_status <- renderText({
    as.character(ama_ref_status() %||% "")
  })

  observeEvent(input$run_ama_author_journal, {
    .ama_console_log("Button pressed", "Run AMA author/journal extraction", .button = "AMA Author/Journal")
    ama_ref_status("Processing: checking input and running FLCA-MA-SIL...")
    ama_ref_results_val(NULL)
    .ama_console_log("Status", "cleared previous author/journal results", .button = "AMA Author/Journal")

    tryCatch({
      shiny::withProgress(message = "Running AMA author/journal extraction", value = 0, {
        shiny::incProgress(0.10, detail = "Reading input")
        .ama_console_log("Input", "reading AMA/PubMed text", .button = "AMA Author/Journal")
        txt <- .get_ama_input_text(ama_ref_status)
        if (is.null(txt) || !nzchar(trimws(txt))) stop("No uploaded or pasted AMA/PubMed/Google Scholar/NCKU Scopus text was found.", call. = FALSE)
        .ama_console_log("Input loaded", sprintf("%d characters; %d lines", nchar(txt), length(strsplit(txt, "\n", fixed = TRUE)[[1]])), .button = "AMA Author/Journal")

        shiny::incProgress(0.20, detail = "Splitting references")
        selected_sources_ref <- .ama_selected_ref_sources()
        .ama_console_log("Parsers", paste(selected_sources_ref, collapse = ", "), .button = "AMA Author/Journal")
        refs <- .ama_split_refs_selected(txt, selected_sources_ref)
        .ama_console_log("References", sprintf("detected normalized/selected entries = %d", length(refs)), .button = "AMA Author/Journal")
        if (!length(refs)) stop("No usable PubMed/AMA/Google Scholar/NCKU Scopus entries were detected.", call. = FALSE)

        shiny::incProgress(0.20, detail = "Extracting journals and authors")
        gs_profile_summary_ref <- if ("google" %in% selected_sources_ref) tryCatch(.ref_google_scholar_profile_summary(txt), error = function(e) data.frame()) else data.frame()
        gs_articles_ref <- .ama_parse_articles_selected(txt, selected_sources_ref)
        if (is.data.frame(gs_articles_ref) && nrow(gs_articles_ref) > 0) {
          wide <- .ref_google_scholar_articles_to_wide(gs_articles_ref)
          .ama_console_log("Parser override", sprintf("Google Scholar/NCKU article parser used: rows = %d; unique journals = %d", nrow(wide), length(unique(wide$Journal[nzchar(wide$Journal)]))), .button = "AMA Author/Journal")
        } else {
          if (!("ama" %in% selected_sources_ref)) stop("No normalized Google Scholar/NCKU rows were extracted, and AMA/PubMed parser is unchecked.", call. = FALSE)
          wide <- .ref_to_wide_from_text_strict(txt)
          wide_gs_force <- if ("google" %in% selected_sources_ref) tryCatch(.ref_to_wide_google_scholar_force(txt), error = function(e) data.frame()) else data.frame()
          if (is.data.frame(wide_gs_force) && nrow(wide_gs_force) >= 2L) {
            # Prefer the force parser for Google Scholar profile pasted rows.
            # This keeps the journal/source line as an actual Journal node.
            wide <- wide_gs_force
            .ama_console_log("Parser override", sprintf("Google Scholar force parser used: rows = %d; unique journals = %d", nrow(wide), length(unique(wide$Journal[nzchar(wide$Journal)]))), .button = "AMA Author/Journal")
          }
        }
        .ama_console_log("Parser done", sprintf("wide rows = %d; columns = %s", nrow(wide), paste(names(wide), collapse = ", ")), .button = "AMA Author/Journal")

        shiny::incProgress(0.15, detail = "Building nodes and edges")
        ne <- .ref_build_nodes_edges(wide, top_n = Inf)
        if ("Journal" %in% names(wide) && "name" %in% names(ne$nodes)) {
          jset <- unique(trimws(as.character(wide$Journal)))
          jset <- jset[!is.na(jset) & nzchar(jset)]
          if (!"term_type" %in% names(ne$nodes)) ne$nodes$term_type <- "Author"
          ne$nodes$term_type[ne$nodes$name %in% jset] <- "Journal"
        }
        .ama_console_log("Network data done", sprintf("raw nodes = %d; raw edges = %d; journal nodes = %d", nrow(ne$nodes), nrow(ne$edges), if ("term_type" %in% names(ne$nodes)) sum(ne$nodes$term_type == "Journal", na.rm = TRUE) else 0L), .button = "AMA Author/Journal")

        shiny::incProgress(0.25, detail = "Running FLCA-MA-SIL Top20")
        .ama_console_log("FLCA-MA-SIL", "running Top20 selection", .button = "AMA Author/Journal")
        fl <- .ref_run_flca_top20(ne$nodes, ne$edges)
        .ama_console_log("FLCA-MA-SIL done", sprintf("engine = %s", fl$engine %||% "unknown"), .button = "AMA Author/Journal")
        # Accept both the original FLCA-MA-SIL engine and the Google Scholar
        # Journal/Author Frequency Top20 engine.  The latter deliberately
        # preserves high-frequency journal nodes (e.g., Medicine) before plotting.
        payload <- .ref_sanitize_plot_payload(fl$nodes, fl$edges)
        nd <- payload$nodes
        ed <- payload$edges
        if (!is.data.frame(nd) || nrow(nd) == 0) {
          stop("FLCA-MA-SIL did not return a valid Top20 result.", call. = FALSE)
        }
        .ama_console_log("Payload", sprintf("Top20 nodes = %d; edges = %d", nrow(nd), nrow(ed)), .button = "AMA Author/Journal")
        if (!is.data.frame(nd) || !nrow(nd)) stop("FLCA-MA-SIL returned no Top20 nodes.", call. = FALSE)

        gs_since_year_ref <- as.integer(format(Sys.Date(), "%Y")) - 5L
        ama_ref_results_val(list(
          refs = refs,
          wide = wide,
          nodes = nd,
          edges = ed,
          engine = fl$engine,
          raw_nodes = ne$nodes,
          raw_edges = ne$edges,
          gs_articles = if (exists("gs_articles_ref", inherits = FALSE)) gs_articles_ref else data.frame(),
          gs_profile_summary = if (exists("gs_profile_summary_ref", inherits = FALSE)) gs_profile_summary_ref else data.frame(),
          gs_since_year = gs_since_year_ref
        ))

        ama_ref_status(paste0("Completed via ", fl$engine, ": parsed ", length(refs),
                              " references; ", nrow(wide), " journal/author rows; ",
                              nrow(nd), " Top20 nodes; ", nrow(ed), " edges."))
        .ama_console_log("Completed", sprintf("references = %d; wide rows = %d; Top20 nodes = %d; edges = %d", length(refs), nrow(wide), nrow(nd), nrow(ed)), .button = "AMA Author/Journal")
        showNotification("AMA author/journal extraction completed via FLCA-MA-SIL.", type = "message")
        updateTabsetPanel(session, "main_tabs", selected = "AMA Author/Journal")
        shiny::incProgress(0.10, detail = "Done")
      })
    }, error = function(e) {
      .ama_console_log("Error", conditionMessage(e), .button = "AMA Author/Journal")
      ama_ref_status(paste("Error:", conditionMessage(e)))
      showNotification(paste("AMA author/journal extraction failed:", conditionMessage(e)), type = "error", duration = 10)
      # keep app alive; show a minimal empty payload so lower outputs do not crash
      ama_ref_results_val(NULL)
    })
  }, ignoreInit = TRUE)

  # ---- AMA Author/Journal outputs via FLCA-process plotting logic ----
  # These bindings deliberately use tableOutput/renderTable and plotOutput/renderPlot.
  # No DT, no htmlwidgets, no imageOutput, and no htmltools::HTML are used here,
  # which avoids the repeated "not all is.character(txt)" display error.

  .ama_display_df <- function(df, max_rows = 100) {
    if (is.null(df) || !is.data.frame(df) || !nrow(df)) {
      return(data.frame(Message = "No records available.", stringsAsFactors = FALSE))
    }
    df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
    df <- .ref_format_numeric_2(df)
    for (j in seq_along(df)) {
      if (is.list(df[[j]]) || is.data.frame(df[[j]])) {
        df[[j]] <- vapply(df[[j]], function(z) paste(as.character(unlist(z)), collapse = "; "), character(1))
      }
      if (is.numeric(df[[j]])) {
        if (names(df)[j] %in% c("rank", "carac", "cluster", "membership")) {
          df[[j]] <- as.integer(round(df[[j]], 0))
        } else {
          df[[j]] <- round(df[[j]], 2)
        }
      } else {
        df[[j]] <- as.character(df[[j]])
        df[[j]][is.na(df[[j]])] <- ""
      }
    }
    if (nrow(df) > max_rows) df <- df[seq_len(max_rows), , drop = FALSE]
    rownames(df) <- NULL
    df
  }

  .ama_payload <- function() {
    r <- ama_ref_results_val()
    validate(need(!is.null(r), "Click Run AMA author/journal extraction first."))
    p <- .ref_sanitize_plot_payload(r$nodes, r$edges)
    list(nodes = p$nodes, edges = p$edges, wide = r$wide)
  }

  # ---- AMA Author/Journal plain table outputs ----
  # Use renderPrint instead of renderText/DT/HTML. This avoids htmltools::HTML(txt)
  # and the repeated "not all is.character(txt)" error in the lower tables.

  .ama_plain_table <- function(df, max_rows = 100) {
    if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) {
      return(data.frame(Message = "No records available.", stringsAsFactors = FALSE))
    }
    df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
    if (nrow(df) > max_rows) df <- df[seq_len(max_rows), , drop = FALSE]

    # Convert list/data-frame columns to character safely.
    for (j in seq_along(df)) {
      if (is.list(df[[j]]) || is.data.frame(df[[j]])) {
        df[[j]] <- vapply(df[[j]], function(z) paste(as.character(unlist(z, use.names = FALSE)), collapse = "; "), character(1))
      }
    }

    # Numeric formatting for display only.
    int_cols <- intersect(c("rank", "carac", "cluster", "membership"), names(df))
    for (nm in names(df)) {
      if (is.numeric(df[[nm]])) {
        if (nm %in% int_cols) {
          df[[nm]] <- as.integer(round(df[[nm]], 0))
        } else {
          df[[nm]] <- sprintf("%.2f", df[[nm]])
        }
      } else {
        df[[nm]] <- as.character(df[[nm]])
        df[[nm]][is.na(df[[nm]])] <- ""
      }
    }
    rownames(df) <- NULL
    df
  }

  .ama_top20_plain <- function() {
    r <- ama_ref_results_val()
    if (is.null(r) || is.null(r$nodes) || !is.data.frame(r$nodes)) {
      return(data.frame(Message = "Click Run AMA author/journal extraction first.", stringsAsFactors = FALSE))
    }
    nd <- as.data.frame(r$nodes, stringsAsFactors = FALSE, check.names = FALSE)
    if (!"name" %in% names(nd)) nd$name <- seq_len(nrow(nd))
    nd$name <- as.character(nd$name)
    nd$name[is.na(nd$name)] <- ""
    if (!"value" %in% names(nd)) nd$value <- 0
    if (!"value2" %in% names(nd)) nd$value2 <- nd$value
    if (!"carac" %in% names(nd)) nd$carac <- 1L
    if (!"ssi" %in% names(nd) && "SSi" %in% names(nd)) nd$ssi <- nd$SSi
    if (!"ssi" %in% names(nd)) nd$ssi <- 0
    if (!"a_i" %in% names(nd)) nd$a_i <- 0
    if (!"b_i" %in% names(nd)) nd$b_i <- 0
    if (!"a_star1" %in% names(nd) && "a_star" %in% names(nd)) nd$a_star1 <- nd$a_star
    if (!"a_star1" %in% names(nd)) nd$a_star1 <- 0

    # Keep module-returned carac; do not recode it here.
    for (cc in c("value", "value2", "carac", "ssi", "a_i", "b_i", "a_star1")) {
      nd[[cc]] <- suppressWarnings(as.numeric(nd[[cc]]))
      nd[[cc]][!is.finite(nd[[cc]]) | is.na(nd[[cc]])] <- 0
    }
    keep <- intersect(c("name", "value", "value2", "carac", "ssi", "a_i", "b_i", "a_star1"), names(nd))
    out <- nd[, keep, drop = FALSE]
    out$rank <- seq_len(nrow(out))
    out <- out[, c("rank", setdiff(names(out), "rank")), drop = FALSE]
    .ama_plain_table(out, max_rows = 20)
  }

  output$tbl_ama_ref_top20_check <- renderPrint({
    print(.ama_top20_plain(), row.names = FALSE, right = FALSE)
  })

  output$tbl_ama_ref_nodes <- renderPrint({
    print(.ama_top20_plain(), row.names = FALSE, right = FALSE)
  })

  output$tbl_ama_ref_edges <- renderPrint({
    r <- ama_ref_results_val()
    if (is.null(r) || is.null(r$edges) || !is.data.frame(r$edges)) {
      print(data.frame(Message = "No edges available.", stringsAsFactors = FALSE), row.names = FALSE, right = FALSE)
    } else {
      print(.ama_plain_table(r$edges, max_rows = 100), row.names = FALSE, right = FALSE)
    }
  })

  output$tbl_ama_ref_wide <- renderPrint({
    r <- ama_ref_results_val()
    if (is.null(r) || is.null(r$wide) || !is.data.frame(r$wide)) {
      print(data.frame(Message = "No parsed journal/author table available.", stringsAsFactors = FALSE), row.names = FALSE, right = FALSE)
    } else {
      print(.ama_plain_table(r$wide, max_rows = 100), row.names = FALSE, right = FALSE)
    }
  })

  .ama_flca_plot_network <- function(nodes, edges) {
    nodes <- as.data.frame(nodes, stringsAsFactors = FALSE)
    edges <- as.data.frame(edges, stringsAsFactors = FALSE)
    nodes$name <- as.character(nodes$name)
    nodes$value <- suppressWarnings(as.numeric(nodes$value)); nodes$value[!is.finite(nodes$value)] <- 1
    if (!"carac" %in% names(nodes)) nodes$carac <- 1L
    nodes$carac <- suppressWarnings(as.integer(nodes$carac)); nodes$carac[!is.finite(nodes$carac)] <- 1L
    n <- nrow(nodes)
    if (n < 1) { plot.new(); text(0.5, 0.5, "No nodes"); return(invisible()) }
    theta <- seq(0, 2*pi, length.out = n + 1)[seq_len(n)]
    xy <- data.frame(name = nodes$name, x = cos(theta), y = sin(theta), stringsAsFactors = FALSE)
    rownames(xy) <- xy$name
    plot(xy$x, xy$y, type = "n", axes = FALSE, xlab = "", ylab = "",
         main = "Network from FLCA-MA-SIL Top20; numbers match Top20 rank")
    if (nrow(edges)) {
      if ("follower" %in% names(edges) && !"Follower" %in% names(edges)) names(edges)[names(edges)=="follower"] <- "Follower"
      edges$Leader <- as.character(edges$Leader); edges$Follower <- as.character(edges$Follower)
      edges <- edges[edges$Leader %in% nodes$name & edges$Follower %in% nodes$name, , drop = FALSE]
      if (nrow(edges)) {
        for (i in seq_len(nrow(edges))) {
          segments(xy[edges$Leader[i], "x"], xy[edges$Leader[i], "y"],
                   xy[edges$Follower[i], "x"], xy[edges$Follower[i], "y"], col = "grey75")
        }
      }
    }
    cols <- grDevices::hcl.colors(max(3, length(unique(nodes$carac))), "Dark 3")
    cfac <- as.integer(factor(nodes$carac))
    cexv <- 1.2 + 2.8 * sqrt(nodes$value) / max(sqrt(nodes$value), na.rm = TRUE)
    points(xy$x, xy$y, pch = 21, bg = cols[cfac], cex = cexv)
    text(xy$x, xy$y, labels = as.character(seq_len(n)), cex = 0.75)
    legend("topright", legend = paste("Cluster", sort(unique(nodes$carac))),
           pt.bg = cols[seq_along(sort(unique(nodes$carac)))], pch = 21, bty = "n", cex = 0.8)
  }

  output$vn_ama_ref <- renderPlot({
    tryCatch({
      payload <- .ama_payload()
      .ama_flca_plot_network(payload$nodes, payload$edges)
    }, error = function(e) { plot.new(); text(0.5, 0.5, paste("Network error:", conditionMessage(e)), col = "red") })
  })

  output$plt_ama_ref_ss <- renderPlot({
    tryCatch({
      payload <- .ama_payload()
      nodes <- payload$nodes; edges <- payload$edges
      nodes$name <- as.character(nodes$name)
      nodes$carac <- suppressWarnings(as.integer(nodes$carac)); nodes$carac[!is.finite(nodes$carac)] <- 1L
      if (!"sil_width" %in% names(nodes)) {
        nodes$sil_width <- if ("ssi" %in% names(nodes)) suppressWarnings(as.numeric(nodes$ssi)) else 0
      }
      nodes$sil_width[!is.finite(nodes$sil_width)] <- 0
      ok <- FALSE
      if (exists("render_panel", mode = "function")) {
        ok <- tryCatch({
          draw_expr <- quote(render_panel(sil_df = nodes, nodes0 = nodes, nodes = nodes, top_n = 20,
                                          footer_label = "AMA Journals/Authors via FLCA-MA-SIL"))
          if (exists(".with_clean_ssplot_aac_labels", mode = "function")) {
            .with_clean_ssplot_aac_labels(eval(draw_expr))
          } else {
            eval(draw_expr)
          }
          TRUE
        }, error = function(e) FALSE)
      }
      if (!isTRUE(ok)) {
        ss <- suppressWarnings(as.numeric(nodes$sil_width)); ss[!is.finite(ss)] <- 0
        plot(ss, seq_along(ss), xlim = c(-1, 1), pch = 21, bg = "lightblue",
             xlab = "Silhouette score", ylab = "Top20 rank",
             main = "SSplot from FLCA-MA-SIL Top20")
        abline(v = 0, lty = 2, col = "grey50")
        text(ss, seq_along(ss), labels = seq_along(ss), cex = 0.75, pos = 3)
      }
    }, error = function(e) { plot.new(); text(0.5, 0.5, paste("SSplot error:", conditionMessage(e)), col = "red") })
  })

  output$plt_ama_ref_kano <- renderPlot({
    tryCatch({
      payload <- .ama_payload()
      nodes <- payload$nodes; edges <- payload$edges
      nodes$name <- as.character(nodes$name)
      nodes$value <- suppressWarnings(as.numeric(nodes$value)); nodes$value[!is.finite(nodes$value)] <- 0
      nodes$value2 <- suppressWarnings(as.numeric(nodes$value2)); nodes$value2[!is.finite(nodes$value2)] <- nodes$value[!is.finite(nodes$value2)]
      nodes$carac <- suppressWarnings(as.integer(nodes$carac)); nodes$carac[!is.finite(nodes$carac)] <- 1L
      ok <- FALSE
      if (exists("plot_kano_real", mode = "function")) {
        ok <- tryCatch({
          print(plot_kano_real(nodes = nodes, edges = edges, title_txt = "Kano plot from FLCA-MA-SIL Top20"))
          TRUE
        }, error = function(e) FALSE)
      } else if (exists("kano_plot", mode = "function")) {
        ok <- tryCatch({
          print(kano_plot(nodes, edges, title_txt = "Kano plot from FLCA-MA-SIL Top20"))
          TRUE
        }, error = function(e) FALSE)
      }
      if (!isTRUE(ok)) {
        plot(nodes$value2, nodes$value, pch = 21, bg = "lightblue",
             xlab = "value2", ylab = "value", main = "Kano plot from FLCA-MA-SIL Top20")
        abline(v = median(nodes$value2, na.rm = TRUE), h = median(nodes$value, na.rm = TRUE), lty = 2, col = "grey50")
        text(nodes$value2, nodes$value, labels = seq_len(nrow(nodes)), cex = 0.75, pos = 3)
      }
    }, error = function(e) { plot.new(); text(0.5, 0.5, paste("Kano error:", conditionMessage(e)), col = "red") })
  })



  # ============================================================
  # FINAL FIX: AMA Author/Journal plots as native inline SVG
  # ------------------------------------------------------------
  # The old plotOutput/renderPlot path can leave only <img alt="Plot object">
  # in the AMA tab.  These output IDs are now uiOutput/renderUI, so the plots
  # are real SVG elements in the page, not Shiny image files.
  # ============================================================

  .ama_svg_escape <- function(x) {
    x <- as.character(x)
    x[is.na(x)] <- ""
    htmltools::htmlEscape(x)
  }

  .ama_svg_num <- function(x, default = 0) {
    y <- suppressWarnings(as.numeric(x))
    y[!is.finite(y) | is.na(y)] <- default
    y
  }


  # ---- SSplot label cleanup for AMA Author/Journal ----
  # renderSSplot.R may print both a standalone AAC number and a cluster label
  # beginning with "AAC=...".  This wrapper keeps the standalone number but
  # removes the duplicated "AAC=" prefix in text/mtext labels.
  .ama_clean_aac_label <- function(labels) {
    z <- as.character(labels)
    z[is.na(z)] <- ""

    # Keep the small standalone AAC number printed beside each cluster,
    # but remove the duplicated AAC value embedded in the large red cluster label.
    # Examples converted:
    #   "AAC=0.89 + SS=0.32 ..." -> "SS=0.32 ..."
    #   "0.89 + SS=0.32 ..."     -> "SS=0.32 ..."
    #   "AAC=0.89 | SS=0.32 ..." -> "SS=0.32 ..."
    z <- gsub("^\\s*AAC\\s*=\\s*[-+]?\\d+(?:\\.\\d+)?\\s*(?:\\+|\\||,|;)\\s*", "", z, perl = TRUE)
    z <- gsub("^\\s*[-+]?\\d+(?:\\.\\d+)?\\s*(?:\\+|\\||,|;)\\s*(?=SS\\s*=)", "", z, perl = TRUE)

    # For a label that is only "AAC=0.xx", keep the numeric value.
    only_aac <- grepl("^\\s*AAC\\s*=\\s*([-+]?\\d+(?:\\.\\d+)?)\\s*$", z, perl = TRUE)
    z[only_aac] <- sub("^\\s*AAC\\s*=\\s*([-+]?\\d+(?:\\.\\d+)?)\\s*$", "\\1", z[only_aac], perl = TRUE)

    # Remove any remaining inline AAC=... fragment, leaving SS/a*/n text intact.
    z <- gsub("\\bAAC\\s*=\\s*[-+]?\\d+(?:\\.\\d+)?\\s*", "", z, perl = TRUE)
    z <- gsub("^\\s*\\|\\s*", "", z, perl = TRUE)
    z <- gsub("\\s{2,}", " ", z, perl = TRUE)
    trimws(z)
  }


  .with_clean_ssplot_aac_labels <- function(expr) {
    old_text_exists <- exists("text", envir = .GlobalEnv, inherits = FALSE)
    old_mtext_exists <- exists("mtext", envir = .GlobalEnv, inherits = FALSE)
    old_text <- if (old_text_exists) get("text", envir = .GlobalEnv, inherits = FALSE) else NULL
    old_mtext <- if (old_mtext_exists) get("mtext", envir = .GlobalEnv, inherits = FALSE) else NULL

    # Some renderSSplot.R versions call graphics::text()/graphics::mtext() explicitly.
    # Therefore, clean both the global free-variable lookup and the graphics namespace
    # temporarily, then restore them immediately after drawing.
    graphics_ns <- asNamespace("graphics")
    old_ns_text <- get("text", envir = graphics_ns)
    old_ns_mtext <- get("mtext", envir = graphics_ns)

    text_wrapper <- function(x, y = NULL, labels, ...) {
      if (missing(labels) || is.null(labels) || length(labels) == 0L) {
        old_ns_text(x = x, y = y, ...)
      } else {
        old_ns_text(x = x, y = y, labels = .ama_clean_aac_label(labels), ...)
      }
    }
    mtext_wrapper <- function(text, ...) {
      if (missing(text) || is.null(text) || length(text) == 0L) {
        old_ns_mtext(...)
      } else {
        old_ns_mtext(text = .ama_clean_aac_label(text), ...)
      }
    }

    assign("text", text_wrapper, envir = .GlobalEnv)
    assign("mtext", mtext_wrapper, envir = .GlobalEnv)

    ns_text_unlocked <- FALSE
    ns_mtext_unlocked <- FALSE
    try({ unlockBinding("text", graphics_ns); ns_text_unlocked <- TRUE; assign("text", text_wrapper, envir = graphics_ns); lockBinding("text", graphics_ns); ns_text_unlocked <- FALSE }, silent = TRUE)
    try({ unlockBinding("mtext", graphics_ns); ns_mtext_unlocked <- TRUE; assign("mtext", mtext_wrapper, envir = graphics_ns); lockBinding("mtext", graphics_ns); ns_mtext_unlocked <- FALSE }, silent = TRUE)

    on.exit({
      if (old_text_exists) assign("text", old_text, envir = .GlobalEnv) else if (exists("text", envir = .GlobalEnv, inherits = FALSE)) rm("text", envir = .GlobalEnv)
      if (old_mtext_exists) assign("mtext", old_mtext, envir = .GlobalEnv) else if (exists("mtext", envir = .GlobalEnv, inherits = FALSE)) rm("mtext", envir = .GlobalEnv)
      try({ if (bindingIsLocked("text", graphics_ns)) unlockBinding("text", graphics_ns); assign("text", old_ns_text, envir = graphics_ns); lockBinding("text", graphics_ns) }, silent = TRUE)
      try({ if (bindingIsLocked("mtext", graphics_ns)) unlockBinding("mtext", graphics_ns); assign("mtext", old_ns_mtext, envir = graphics_ns); lockBinding("mtext", graphics_ns) }, silent = TRUE)
    }, add = TRUE)

    force(expr)
  }


  # ---- Safe scalar conversion helpers for AMA Author/Journal outputs ----
  # Required by .ama705_payload(), .ama705_network_widget(), and downloads.
  .ama705_chr <- function(x) {
    x <- as.character(x)
    x[is.na(x)] <- ""
    x
  }

  .ama705_num <- function(x, default = 0) {
    z <- suppressWarnings(as.numeric(x))
    if (length(default) > 1L) {
      d <- suppressWarnings(as.numeric(default))
      if (length(d) != length(z)) d <- rep_len(d, length(z))
    } else {
      d <- rep_len(suppressWarnings(as.numeric(default)), length(z))
    }
    d[!is.finite(d) | is.na(d)] <- 0
    z[!is.finite(z) | is.na(z)] <- d[!is.finite(z) | is.na(z)]
    z
  }

  .ama705_escape <- function(x) {
    x <- .ama705_chr(x)
    x <- gsub("&", "&amp;", x, fixed = TRUE)
    x <- gsub("<", "&lt;", x, fixed = TRUE)
    x <- gsub(">", "&gt;", x, fixed = TRUE)
    x <- gsub('"', "&quot;", x, fixed = TRUE)
    x
  }

  .ama705_payload <- function() {
    r <- ama_ref_results_val()
    if (is.null(r) || is.null(r$nodes) || !is.data.frame(r$nodes)) {
      stop("Run AMA author/journal extraction first.", call. = FALSE)
    }
    nd <- as.data.frame(r$nodes, stringsAsFactors = FALSE, check.names = FALSE)
    if (!"name" %in% names(nd)) nd$name <- .ama705_chr(nd[[1]])
    nd$name <- trimws(.ama705_chr(nd$name))
    nd <- nd[nzchar(nd$name), , drop = FALSE]
    nd <- nd[!duplicated(nd$name), , drop = FALSE]
    if (!nrow(nd)) stop("No Top20 nodes available.", call. = FALSE)

    if (!"value" %in% names(nd)) nd$value <- 1
    if (!"value2" %in% names(nd)) nd$value2 <- nd$value
    if (!"carac" %in% names(nd)) nd$carac <- 1
    if (!"ssi" %in% names(nd) && "SSi" %in% names(nd)) nd$ssi <- nd$SSi
    if (!"ssi" %in% names(nd) && "sil_width" %in% names(nd)) nd$ssi <- nd$sil_width
    if (!"ssi" %in% names(nd)) nd$ssi <- 0
    if (!"a_i" %in% names(nd)) nd$a_i <- 0
    if (!"b_i" %in% names(nd)) nd$b_i <- 0
    if (!"a_star1" %in% names(nd) && "a_star" %in% names(nd)) nd$a_star1 <- nd$a_star
    if (!"a_star1" %in% names(nd)) nd$a_star1 <- 0
    nd$value <- .ama705_num(nd$value, 1)
    nd$value2 <- .ama705_num(nd$value2, nd$value)
    nd$carac <- suppressWarnings(as.integer(gsub("[^0-9-]", "", .ama705_chr(nd$carac))))
    nd$carac[!is.finite(nd$carac) | is.na(nd$carac)] <- 1L
    nd$ssi <- .ama705_num(nd$ssi, 0)
    nd$a_i <- .ama705_num(nd$a_i, 0)
    nd$b_i <- .ama705_num(nd$b_i, 0)
    nd$a_star1 <- .ama705_num(nd$a_star1, 0)
    if (!"rank" %in% names(nd)) nd$rank <- seq_len(nrow(nd))
    nd$rank <- suppressWarnings(as.integer(nd$rank))
    nd$rank[!is.finite(nd$rank) | is.na(nd$rank)] <- seq_len(nrow(nd))[!is.finite(nd$rank) | is.na(nd$rank)]
    nd <- nd[order(nd$rank, -nd$value, -nd$value2), , drop = FALSE]
    nd$rank <- seq_len(nrow(nd))
    nd <- utils::head(nd, 20)

    ed <- if (!is.null(r$edges) && is.data.frame(r$edges)) as.data.frame(r$edges, stringsAsFactors = FALSE, check.names = FALSE) else data.frame()
    if (nrow(ed)) {
      if (all(c("Source", "Target") %in% names(ed)) && !all(c("Leader", "Follower") %in% names(ed))) {
        names(ed)[match(c("Source", "Target"), names(ed))] <- c("Leader", "Follower")
      }
      if ("follower" %in% names(ed) && !"Follower" %in% names(ed)) names(ed)[names(ed) == "follower"] <- "Follower"
      if (!"WCD" %in% names(ed)) {
        wname <- intersect(c("weight", "value", "Count", "count", "n"), names(ed))[1]
        ed$WCD <- if (length(wname)) ed[[wname]] else 1
      }
      if (all(c("Leader", "Follower", "WCD") %in% names(ed))) {
        ed <- ed[, c("Leader", "Follower", "WCD"), drop = FALSE]
        ed$Leader <- trimws(.ama705_chr(ed$Leader))
        ed$Follower <- trimws(.ama705_chr(ed$Follower))
        ed$WCD <- .ama705_num(ed$WCD, 1)
        ed <- ed[nzchar(ed$Leader) & nzchar(ed$Follower) & ed$Leader != ed$Follower & ed$WCD > 0, , drop = FALSE]
        ed <- ed[ed$Leader %in% nd$name & ed$Follower %in% nd$name, , drop = FALSE]
      } else {
        ed <- data.frame(Leader = character(), Follower = character(), WCD = numeric(), stringsAsFactors = FALSE)
      }
    } else {
      ed <- data.frame(Leader = character(), Follower = character(), WCD = numeric(), stringsAsFactors = FALSE)
    }
    list(nodes = nd, edges = ed, raw = r)
  }

  .ama705_fmt_table <- function(df, max_rows = Inf) {
    df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
    if (is.finite(max_rows)) df <- utils::head(df, max_rows)
    for (cc in names(df)) {
      if (is.numeric(df[[cc]])) {
        if (cc == "rank" || cc == "carac") df[[cc]] <- as.integer(round(df[[cc]]))
        else df[[cc]] <- sprintf("%.2f", df[[cc]])
      }
    }
    df
  }

  .ama705_one_link_edges <- function(nodes, edges) {
    nd <- as.data.frame(nodes, stringsAsFactors = FALSE)
    ed <- as.data.frame(edges, stringsAsFactors = FALSE)
    if (!nrow(ed)) return(ed)
    ed$Leader <- .ama705_chr(ed$Leader)
    ed$Follower <- .ama705_chr(ed$Follower)
    ed$WCD <- .ama705_num(ed$WCD, 1)
    ed <- ed[ed$Leader %in% nd$name & ed$Follower %in% nd$name & ed$Leader != ed$Follower & ed$WCD > 0, , drop = FALSE]
    if (!nrow(ed)) return(ed)
    # Treat the stronger/higher-value endpoint as leader when direction is ambiguous.
    val <- setNames(.ama705_num(nd$value2, nd$value), nd$name)
    swap <- val[ed$Follower] > val[ed$Leader]
    if (any(swap, na.rm = TRUE)) {
      tmp <- ed$Leader[swap]
      ed$Leader[swap] <- ed$Follower[swap]
      ed$Follower[swap] <- tmp
    }
    ed <- ed[order(ed$Follower, -ed$WCD, -val[ed$Leader], ed$Leader), , drop = FALSE]
    ed <- ed[!duplicated(ed$Follower), , drop = FALSE]
    ed <- ed[order(match(ed$Follower, nd$name)), , drop = FALSE]
    rownames(ed) <- NULL
    ed
  }

  .ama705_network_widget <- function() {
    z <- .ama705_payload()
    nd <- z$nodes
    ed <- .ama705_one_link_edges(nd, z$edges)
    mx <- max(nd$value, na.rm = TRUE); if (!is.finite(mx) || mx <= 0) mx <- 1
    groups <- paste0("Cluster ", nd$carac)
    vn_nodes <- data.frame(
      id = .ama705_chr(nd$name),
      label = .ama705_chr(nd$name),
      group = .ama705_chr(groups),
      value = as.numeric(pmax(12, 20 + 55 * sqrt(pmax(nd$value, 0)) / sqrt(mx))),
      title = .ama705_escape(paste0("Rank: ", nd$rank,
                                    "<br>Name: ", nd$name,
                                    "<br>Cluster: ", nd$carac,
                                    "<br>value: ", sprintf("%.2f", nd$value),
                                    "<br>value2: ", sprintf("%.2f", nd$value2),
                                    "<br>SS: ", sprintf("%.2f", nd$ssi))),
      stringsAsFactors = FALSE
    )
    vn_edges <- data.frame(
      from = .ama705_chr(ed$Leader),
      to = .ama705_chr(ed$Follower),
      value = as.numeric(pmax(1, ed$WCD)),
      label = .ama705_chr(sprintf("%.0f", ed$WCD)),
      title = .ama705_escape(paste0(ed$Leader, " -> ", ed$Follower, "<br>WCD=", sprintf("%.2f", ed$WCD))),
      arrows = rep("to", nrow(ed)),
      stringsAsFactors = FALSE
    )
    visNetwork::visNetwork(vn_nodes, vn_edges, width = "100%", height = "760px",
                           main = "Interactive network: FLCA-MA-SIL Top20; one primary edge per follower") %>%
      visNetwork::visNodes(shape = "dot",
                            font = list(size = 26, face = "bold", color = "#111111", strokeWidth = 3)) %>%
      visNetwork::visEdges(smooth = list(enabled = TRUE, type = "dynamic"),
                            font = list(size = 18, face = "bold", align = "middle"),
                            color = list(color = "#888888", highlight = "#d62728")) %>%
      visNetwork::visOptions(highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE),
                              nodesIdSelection = TRUE,
                              selectedBy = list(variable = "group", multiple = TRUE)) %>%
      visNetwork::visPhysics(stabilization = TRUE, solver = "forceAtlas2Based") %>%
      visNetwork::visLayout(randomSeed = 705)
  }




  .ama705_modularity_results <- function(nd, edges) {
    nd <- as.data.frame(nd, stringsAsFactors = FALSE, check.names = FALSE)
    nd$name <- .ama705_chr(nd$name)
    nd$carac <- suppressWarnings(as.integer(nd$carac))
    nd$carac[!is.finite(nd$carac) | is.na(nd$carac)] <- 1L
    nd$sil_width <- .ama705_num(nd$sil_width, 0)
    clv <- sort(unique(stats::na.omit(nd$carac)))
    if (!length(clv)) clv <- 1L

    base_results <- do.call(rbind, lapply(clv, function(cc) {
      sub <- nd[nd$carac == cc, , drop = FALSE]
      data.frame(
        Cluster = paste0("C", cc),
        SS = mean(sub$sil_width, na.rm = TRUE),
        Qw = 0, Qu = 0, Q = 0,
        D_GiniSimpson = NA_real_, Q_over_D = 0,
        OneMinus_1_over_k = NA_real_, Q_over_Dmax_eff = 0,
        n = nrow(sub), stringsAsFactors = FALSE
      )
    }))

    overall <- data.frame(
      Cluster = "OVERALL",
      SS = mean(nd$sil_width, na.rm = TRUE),
      Qw = 0, Qu = 0, Q = 0,
      D_GiniSimpson = NA_real_, Q_over_D = 0,
      OneMinus_1_over_k = NA_real_, Q_over_Dmax_eff = 0,
      n = nrow(nd), stringsAsFactors = FALSE
    )

    # Modularity is computed from FLCA cluster membership (nodes$carac) and
    # the full Top20 edge set returned by the FLCA-MA-SIL process.  The one-link
    # edge set is used for the interactive dashboard only; Q should not be
    # forced to zero when the FLCA module has valid edges.
    ed <- as.data.frame(edges, stringsAsFactors = FALSE, check.names = FALSE)
    if (nrow(ed) && all(c("Leader", "Follower", "WCD") %in% names(ed))) {
      ed$Leader <- .ama705_chr(ed$Leader)
      ed$Follower <- .ama705_chr(ed$Follower)
      ed$WCD <- .ama705_num(ed$WCD, 1)
      ed <- ed[ed$Leader %in% nd$name & ed$Follower %in% nd$name &
                 ed$Leader != ed$Follower & ed$WCD > 0, , drop = FALSE]
      if (nrow(ed)) {
        ed$a <- pmin(ed$Leader, ed$Follower)
        ed$b <- pmax(ed$Leader, ed$Follower)
        edw <- stats::aggregate(WCD ~ a + b, data = ed, FUN = sum)
        ed1 <- edw
        ed1$WCD <- 1

        .q_contrib <- function(e2, weighted = TRUE) {
          w <- if (weighted) .ama705_num(e2$WCD, 1) else rep(1, nrow(e2))
          m <- sum(w, na.rm = TRUE)
          out <- setNames(rep(0, length(clv)), paste0("C", clv))
          if (!is.finite(m) || m <= 0) return(list(overall = 0, per = out))
          cl <- setNames(nd$carac, nd$name)
          deg <- setNames(rep(0, nrow(nd)), nd$name)
          for (i in seq_len(nrow(e2))) {
            deg[e2$a[i]] <- deg[e2$a[i]] + w[i]
            deg[e2$b[i]] <- deg[e2$b[i]] + w[i]
          }
          for (cc in clv) {
            members <- nd$name[nd$carac == cc]
            internal <- e2$a %in% members & e2$b %in% members
            l_c <- sum(w[internal], na.rm = TRUE)
            d_c <- sum(deg[members], na.rm = TRUE)
            q_c <- (l_c / m) - (d_c / (2 * m))^2
            if (!is.finite(q_c)) q_c <- 0
            out[paste0("C", cc)] <- q_c
          }
          list(overall = sum(out, na.rm = TRUE), per = out)
        }

        qw <- .q_contrib(edw, weighted = TRUE)
        qu <- .q_contrib(ed1, weighted = FALSE)
        base_results$Qw <- as.numeric(qw$per[paste0("C", clv)])
        base_results$Qu <- as.numeric(qu$per[paste0("C", clv)])
        base_results$Q <- base_results$Qw
        overall$Qw <- qw$overall
        overall$Qu <- qu$overall
        overall$Q <- qw$overall
      }
    }

    cc <- nd$carac[!is.na(nd$carac)]
    if (length(cc)) {
      pk <- as.numeric(table(cc)) / length(cc)
      D <- 1 - sum(pk^2)
      k <- length(unique(cc))
      Dmax <- if (k > 0) 1 - 1/k else NA_real_
      overall$D_GiniSimpson <- D
      overall$OneMinus_1_over_k <- Dmax
      overall$Q_over_D <- if (is.finite(D) && D > 0) overall$Qw / D else 0
      overall$Q_over_Dmax_eff <- if (is.finite(Dmax) && Dmax > 0) overall$Qw / Dmax else 0
    }
    base_results$SS[!is.finite(base_results$SS)] <- 0
    overall$SS[!is.finite(overall$SS)] <- 0
    rbind(overall, base_results)
  }

  .ama705_ss_inputs <- function() {
    z <- .ama705_payload()
    nd <- z$nodes
    ed1 <- .ama705_one_link_edges(nd, z$edges)
    nd$sil_width <- nd$ssi
    nd$wsel <- NA_real_
    nd$role <- NA_character_
    nd$neighbor_name <- nd$name
    nd$neighborC <- nd$carac
    if (nrow(ed1)) {
      for (i in seq_len(nrow(ed1))) {
        j <- match(ed1$Follower[i], nd$name)
        k <- match(ed1$Leader[i], nd$name)
        if (!is.na(j)) {
          nd$wsel[j] <- ed1$WCD[i]
          nd$role[j] <- "follower"
          nd$neighbor_name[j] <- ed1$Leader[i]
          nd$neighborC[j] <- if (!is.na(k)) nd$carac[k] else nd$carac[j]
        }
      }
    }
    leaders <- tapply(seq_len(nrow(nd)), nd$carac, function(ii) ii[order(-nd$value2[ii], -nd$value[ii], nd$rank[ii])][1])
    nd$role[as.integer(leaders)] <- "leader"
    results <- .ama705_modularity_results(nd, z$edges)
    list(sil_df = nd, nodes0 = nd, nodes = nd, results = results)
  }

  .ama705_draw_ssplot <- function(font_scale = 1.55,
                                  footer_label = "AMA Journals/Authors via FLCA-MA-SIL") {
    inp <- .ama705_ss_inputs()
    if (!exists("render_panel", mode = "function")) {
      stop("render_panel() not loaded from renderSSplot.R or renderSSplot(79).R.", call. = FALSE)
    }
    .with_clean_ssplot_aac_labels({
      render_panel(sil_df = inp$sil_df, nodes0 = inp$nodes0, results = inp$results,
                   nodes = inp$nodes, top_n = nrow(inp$sil_df), font_scale = font_scale,
                   aac_side = "left", neighbor_side = "right", neighbor_on_bar = TRUE,
                   footer_label = footer_label)
    })
  }


  .ama705_print_kano_object <- function(obj) {
    if (inherits(obj, "ggplot")) {
      obj <- obj +
        ggplot2::theme(
          aspect.ratio = 1.35,
          text = ggplot2::element_text(face = "bold"),
          plot.title = ggplot2::element_text(face = "bold", size = 28),
          axis.title = ggplot2::element_text(face = "bold", size = 22),
          axis.text = ggplot2::element_text(face = "bold", size = 16),
          legend.title = ggplot2::element_text(face = "bold", size = 18),
          legend.text = ggplot2::element_text(face = "bold", size = 15),
          plot.margin = ggplot2::margin(18, 24, 18, 24)
        )
      print(obj)
    } else if (inherits(obj, c("grob", "recordedplot"))) {
      print(obj)
    } else if (!is.null(obj)) {
      print(obj)
    }
    invisible(TRUE)
  }

  .ama705_draw_kano <- function() {
    z <- .ama705_payload()
    nd <- z$nodes
    ed <- .ama705_one_link_edges(nd, z$edges)
    ok <- FALSE
    if (exists("plot_kano_real", mode = "function")) {
      ok <- tryCatch({
        obj <- plot_kano_real(nodes = nd, edges = ed, title_txt = "Kano plot from FLCA-MA-SIL Top20")
        .ama705_print_kano_object(obj)
        TRUE
      }, error = function(e) FALSE)
    }
    if (!ok && exists("kano_plot", mode = "function")) {
      ok <- tryCatch({
        obj <- kano_plot(nd, ed, title_txt = "Kano plot from FLCA-MA-SIL Top20")
        .ama705_print_kano_object(obj)
        TRUE
      }, error = function(e) FALSE)
    }
    if (ok) return(invisible(TRUE))
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      nd$cluster <- factor(nd$carac)
      medx <- stats::median(nd$value2, na.rm = TRUE)
      medy <- stats::median(nd$value, na.rm = TRUE)
      p <- ggplot2::ggplot(nd, ggplot2::aes(x = value2, y = value, size = value, color = cluster, label = name)) +
        ggplot2::geom_vline(xintercept = medx, linetype = "dashed", linewidth = 1.2, color = "red") +
        ggplot2::geom_hline(yintercept = medy, linetype = "dashed", linewidth = 1.2, color = "red") +
        ggplot2::geom_point(alpha = 0.85) +
        ggplot2::geom_text(fontface = "bold", size = 5.5, vjust = -1.05, show.legend = FALSE) +
        ggplot2::scale_size_continuous(range = c(5, 18)) +
        ggplot2::labs(title = "Kano plot from FLCA-MA-SIL Top20",
                      subtitle = "x = value2 (sum of WCD); y = value; color = FLCA cluster",
                      x = "value2", y = "value", color = "Cluster", size = "value") +
        ggplot2::theme_bw(base_size = 21) +
        ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 28),
                       plot.subtitle = ggplot2::element_text(face = "bold", size = 18),
                       axis.title = ggplot2::element_text(face = "bold", size = 23),
                       axis.text = ggplot2::element_text(face = "bold", size = 17),
                       legend.title = ggplot2::element_text(face = "bold", size = 19),
                       legend.text = ggplot2::element_text(face = "bold", size = 16),
                       aspect.ratio = 1.35,
                       plot.margin = ggplot2::margin(18, 24, 18, 24))
      print(p)
    } else {
      cols <- grDevices::hcl.colors(max(3, length(unique(nd$carac))), "Dark 3")
      bg <- cols[as.integer(factor(nd$carac))]
      par(mar = c(5, 5, 5, 3))
      plot(nd$value2, nd$value, pch = 21, bg = bg, cex = 2.6,
           xlab = "value2", ylab = "value", main = "Kano plot from FLCA-MA-SIL Top20",
           cex.main = 1.9, cex.lab = 1.6, cex.axis = 1.35, font.main = 2, font.lab = 2)
      abline(v = stats::median(nd$value2, na.rm = TRUE), h = stats::median(nd$value, na.rm = TRUE), lty = 2, col = "red", lwd = 2)
      text(nd$value2, nd$value, labels = nd$name, cex = 1.05, font = 2, pos = 3)
      legend("topright", legend = paste("Cluster", sort(unique(nd$carac))), pt.bg = cols[seq_along(sort(unique(nd$carac)))], pch = 21, bty = "n", cex = 1.1)
    }
    invisible(TRUE)
  }



  .ama705_cluster_positions <- function(nd) {
    nd <- as.data.frame(nd, stringsAsFactors = FALSE, check.names = FALSE)
    nd$name <- .ama705_chr(nd$name)
    nd$carac <- suppressWarnings(as.integer(gsub("[^0-9-]", "", .ama705_chr(nd$carac))))
    nd$carac[!is.finite(nd$carac) | is.na(nd$carac)] <- 1L
    nd$value <- .ama705_num(nd$value, 1)
    if (!"value2" %in% names(nd)) nd$value2 <- nd$value
    nd$value2 <- .ama705_num(nd$value2, nd$value)
    if (!"ssi" %in% names(nd)) nd$ssi <- 0
    nd$ssi <- .ama705_num(nd$ssi, 0)
    nd$rank <- if ("rank" %in% names(nd)) suppressWarnings(as.integer(nd$rank)) else seq_len(nrow(nd))
    bad_rank <- !is.finite(nd$rank) | is.na(nd$rank)
    if (any(bad_rank)) nd$rank[bad_rank] <- seq_len(nrow(nd))[bad_rank]

    nd$cluster_leader <- FALSE
    clv <- sort(unique(nd$carac))
    k <- length(clv)
    if (!k) k <- 1L

    center_angle <- seq(pi/2, pi/2 - 2*pi + 2*pi/k, length.out = k)
    centers <- data.frame(carac = clv,
                          cx = 560 * cos(center_angle),
                          cy = 560 * sin(center_angle),
                          stringsAsFactors = FALSE)

    out <- nd
    out$x <- 0
    out$y <- 0
    for (cc in clv) {
      ii <- which(out$carac == cc)
      ii <- ii[order(-out$value2[ii], -out$value[ii], out$rank[ii], out$name[ii])]
      m <- length(ii)
      cen <- centers[centers$carac == cc, , drop = FALSE]
      leader <- ii[1]
      out$cluster_leader[leader] <- TRUE
      out$x[leader] <- cen$cx
      out$y[leader] <- cen$cy
      if (m > 1L) {
        followers <- ii[-1]
        local_r <- 115 + 22 * m
        ang <- seq(pi/2, pi/2 - 2*pi + 2*pi/length(followers), length.out = length(followers))
        out$x[followers] <- cen$cx + local_r * cos(ang)
        out$y[followers] <- cen$cy + local_r * sin(ang)
      }
    }
    out
  }

  .ama705_cluster_leader_edges <- function(nodes, edges) {
    # Build the clustered network as a strict star per FLCA-MA-SIL cluster:
    # the highest value2/value node in each carac is the cluster leader, and
    # every other node in that carac gets exactly one edge to that leader.
    nd <- as.data.frame(nodes, stringsAsFactors = FALSE, check.names = FALSE)
    if (!nrow(nd)) {
      return(data.frame(Leader = character(), Follower = character(), WCD = numeric(),
                        edge_type = character(), stringsAsFactors = FALSE))
    }
    nd$name <- .ama705_chr(nd$name)
    nd$carac <- suppressWarnings(as.integer(gsub("[^0-9-]", "", .ama705_chr(nd$carac))))
    nd$carac[!is.finite(nd$carac) | is.na(nd$carac)] <- 1L
    nd$value <- .ama705_num(nd$value, 1)
    if (!"value2" %in% names(nd)) nd$value2 <- nd$value
    nd$value2 <- .ama705_num(nd$value2, nd$value)
    nd$rank <- if ("rank" %in% names(nd)) suppressWarnings(as.integer(nd$rank)) else seq_len(nrow(nd))
    bad_rank <- !is.finite(nd$rank) | is.na(nd$rank)
    if (any(bad_rank)) nd$rank[bad_rank] <- seq_len(nrow(nd))[bad_rank]

    ed <- as.data.frame(edges, stringsAsFactors = FALSE, check.names = FALSE)
    if (nrow(ed)) {
      if (all(c("Source", "Target") %in% names(ed)) && !all(c("Leader", "Follower") %in% names(ed))) {
        names(ed)[match(c("Source", "Target"), names(ed))] <- c("Leader", "Follower")
      }
      if ("follower" %in% names(ed) && !"Follower" %in% names(ed)) names(ed)[names(ed) == "follower"] <- "Follower"
      if (!"WCD" %in% names(ed)) {
        wname <- intersect(c("weight", "value", "Count", "count", "n"), names(ed))[1]
        ed$WCD <- if (length(wname)) ed[[wname]] else 1
      }
      if (all(c("Leader", "Follower", "WCD") %in% names(ed))) {
        ed <- ed[, c("Leader", "Follower", "WCD"), drop = FALSE]
        ed$Leader <- .ama705_chr(ed$Leader)
        ed$Follower <- .ama705_chr(ed$Follower)
        ed$WCD <- .ama705_num(ed$WCD, 1)
        ed <- ed[ed$Leader %in% nd$name & ed$Follower %in% nd$name & ed$Leader != ed$Follower & ed$WCD > 0, , drop = FALSE]
      } else {
        ed <- data.frame(Leader = character(), Follower = character(), WCD = numeric(), stringsAsFactors = FALSE)
      }
    } else {
      ed <- data.frame(Leader = character(), Follower = character(), WCD = numeric(), stringsAsFactors = FALSE)
    }

    get_pair_w <- function(a, b) {
      if (!nrow(ed)) return(NA_real_)
      jj <- which((ed$Leader == a & ed$Follower == b) | (ed$Leader == b & ed$Follower == a))
      if (!length(jj)) return(NA_real_)
      max(ed$WCD[jj], na.rm = TRUE)
    }

    out <- list()
    q <- 0L
    for (cc in sort(unique(nd$carac))) {
      ii <- which(nd$carac == cc)
      ii <- ii[order(-nd$value2[ii], -nd$value[ii], nd$rank[ii], nd$name[ii])]
      if (length(ii) <= 1L) next
      leader <- nd$name[ii[1]]
      followers <- nd$name[ii[-1]]
      for (ff in followers) {
        w <- get_pair_w(leader, ff)
        et <- "observed leader-follower edge"
        if (!is.finite(w) || is.na(w) || w <= 0) {
          w <- 1
          et <- "cluster membership edge"
        }
        q <- q + 1L
        out[[q]] <- data.frame(Leader = leader, Follower = ff, WCD = w,
                               carac = cc, edge_type = et, stringsAsFactors = FALSE)
      }
    }
    if (!length(out)) {
      return(data.frame(Leader = character(), Follower = character(), WCD = numeric(),
                        carac = integer(), edge_type = character(), stringsAsFactors = FALSE))
    }
    ans <- do.call(rbind, out)
    ans <- ans[order(ans$carac, match(ans$Leader, nd$name), match(ans$Follower, nd$name)), , drop = FALSE]
    rownames(ans) <- NULL
    ans
  }

  .ama705_cluster_network_widget <- function() {
    z <- .ama705_payload()
    nd <- .ama705_cluster_positions(z$nodes)
    ed <- .ama705_cluster_leader_edges(nd, z$edges)
    mx <- max(nd$value, na.rm = TRUE); if (!is.finite(mx) || mx <= 0) mx <- 1
    groups <- paste0("Cluster ", nd$carac)
    vn_nodes <- data.frame(
      id = .ama705_chr(nd$name),
      label = paste0(nd$rank, ". ", .ama705_chr(nd$name)),
      group = .ama705_chr(groups),
      value = as.numeric(pmax(18, 26 + 66 * sqrt(pmax(nd$value, 0)) / sqrt(mx) + ifelse(nd$cluster_leader, 14, 0))),
      x = as.numeric(nd$x),
      y = as.numeric(nd$y),
      fixed = TRUE,
      borderWidth = ifelse(nd$cluster_leader, 5, 2),
      title = .ama705_escape(paste0("Rank: ", nd$rank,
                                    "<br>Name: ", nd$name,
                                    "<br>Cluster: ", nd$carac,
                                    "<br>Role: ", ifelse(nd$cluster_leader, "cluster leader", "follower"),
                                    "<br>value: ", sprintf("%.2f", nd$value),
                                    "<br>value2: ", sprintf("%.2f", nd$value2),
                                    "<br>SS: ", sprintf("%.2f", nd$ssi))),
      stringsAsFactors = FALSE
    )
    vn_edges <- data.frame(
      # Arrow is follower -> leader, so each follower points to its respective cluster leader.
      from = .ama705_chr(ed$Follower),
      to = .ama705_chr(ed$Leader),
      value = as.numeric(pmax(1, ed$WCD)),
      label = .ama705_chr(sprintf("%.0f", ed$WCD)),
      title = .ama705_escape(paste0("Follower: ", ed$Follower,
                                    "<br>Leader: ", ed$Leader,
                                    "<br>Cluster: ", ed$carac,
                                    "<br>WCD=", sprintf("%.2f", ed$WCD),
                                    "<br>", ed$edge_type)),
      arrows = rep("to", nrow(ed)),
      stringsAsFactors = FALSE
    )
    visNetwork::visNetwork(vn_nodes, vn_edges, width = "100%", height = "760px",
                           main = "Clustered network: FLCA-MA-SIL Top20; each follower has one edge to its cluster leader") %>%
      visNetwork::visNodes(shape = "dot",
                            font = list(size = 24, face = "bold", color = "#111111", strokeWidth = 3)) %>%
      visNetwork::visEdges(smooth = list(enabled = TRUE, type = "curvedCW", roundness = 0.18),
                            font = list(size = 18, face = "bold", align = "middle"),
                            color = list(color = "#777777", highlight = "#d62728")) %>%
      visNetwork::visOptions(highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE),
                              nodesIdSelection = TRUE,
                              selectedBy = list(variable = "group", multiple = TRUE)) %>%
      visNetwork::visPhysics(enabled = FALSE) %>%
      visNetwork::visLayout(randomSeed = 706)
  }

  .ama705_draw_cluster_network_png <- function() {
    z <- .ama705_payload()
    nd <- .ama705_cluster_positions(z$nodes)
    ed <- .ama705_cluster_leader_edges(nd, z$edges)
    cols <- grDevices::hcl.colors(max(3, length(unique(nd$carac))), "Dark 3")
    clv <- sort(unique(nd$carac))
    bg <- cols[as.integer(factor(nd$carac, levels = clv))]
    xlim <- range(nd$x, na.rm = TRUE) + c(-240, 240)
    ylim <- range(nd$y, na.rm = TRUE) + c(-200, 200)
    par(mar = c(1, 1, 4, 1))
    plot(nd$x, nd$y, type = "n", axes = FALSE, xlab = "", ylab = "",
         main = "Clustered network PNG: each follower points to its FLCA-MA-SIL cluster leader",
         cex.main = 1.65, font.main = 2, xlim = xlim, ylim = ylim, asp = 1)
    for (cc in clv) {
      sub <- nd[nd$carac == cc, , drop = FALSE]
      if (nrow(sub)) {
        cx <- sub$x[which(sub$cluster_leader)[1]]
        cy <- sub$y[which(sub$cluster_leader)[1]]
        if (!is.finite(cx) || !is.finite(cy)) { cx <- mean(sub$x, na.rm = TRUE); cy <- mean(sub$y, na.rm = TRUE) }
        r <- max(sqrt((sub$x - cx)^2 + (sub$y - cy)^2), 80, na.rm = TRUE) + 75
        symbols(cx, cy, circles = r, inches = FALSE, add = TRUE,
                bg = grDevices::adjustcolor(cols[match(cc, clv)], alpha.f = 0.10),
                fg = grDevices::adjustcolor(cols[match(cc, clv)], alpha.f = 0.70), lwd = 2)
        text(cx, cy + r + 35, labels = paste("Cluster", cc), cex = 1.15, font = 2,
             col = cols[match(cc, clv)])
      }
    }
    if (nrow(ed)) {
      mxw <- max(ed$WCD, na.rm = TRUE); if (!is.finite(mxw) || mxw <= 0) mxw <- 1
      for (i in seq_len(nrow(ed))) {
        a <- match(ed$Leader[i], nd$name)
        b <- match(ed$Follower[i], nd$name)
        if (!is.na(a) && !is.na(b)) {
          # Arrow is follower -> leader.
          arrows(nd$x[b], nd$y[b], nd$x[a], nd$y[a], length = 0.08,
                 lwd = 1.5 + 2 * ed$WCD[i] / mxw, col = "grey45")
          text((nd$x[a] + nd$x[b]) / 2, (nd$y[a] + nd$y[b]) / 2,
               labels = sprintf("%.0f", ed$WCD[i]), cex = 0.75, font = 2, col = "grey30")
        }
      }
    }
    cexv <- 1.5 + 3.0 * sqrt(pmax(nd$value, 0)) / max(sqrt(pmax(nd$value, 0)), na.rm = TRUE)
    cexv[nd$cluster_leader] <- cexv[nd$cluster_leader] + 0.8
    points(nd$x, nd$y, pch = 21, bg = bg, cex = cexv, lwd = ifelse(nd$cluster_leader, 4, 2))
    text(nd$x, nd$y, labels = nd$rank, cex = 0.95, font = 2)
    text(nd$x + 28, nd$y - 28, labels = nd$name, cex = 0.88, font = 2, adj = c(0, 1))
    legend("bottomleft", legend = paste("Cluster", clv),
           pt.bg = cols[seq_along(clv)], pch = 21, bty = "n", cex = 1.0, text.font = 2)
    legend("topright", legend = c("larger border = cluster leader", "arrow: follower -> leader"),
           bty = "n", cex = 1.0, text.font = 2)
  }

  .ama705_draw_network_png <- function() {
    z <- .ama705_payload()
    nd <- z$nodes
    ed <- .ama705_one_link_edges(nd, z$edges)
    n <- nrow(nd)
    theta <- seq(pi/2, pi/2 - 2*pi + 2*pi/n, length.out = n)
    xy <- data.frame(name = nd$name, x = cos(theta), y = sin(theta), stringsAsFactors = FALSE)
    rownames(xy) <- xy$name
    cols <- grDevices::hcl.colors(max(3, length(unique(nd$carac))), "Dark 3")
    bg <- cols[as.integer(factor(nd$carac))]
    par(mar = c(1, 1, 4, 1))
    plot(xy$x, xy$y, type = "n", axes = FALSE, xlab = "", ylab = "",
         main = "Interactive network PNG: FLCA-MA-SIL Top20; one primary edge per follower",
         cex.main = 1.7, font.main = 2, xlim = c(-1.4, 1.4), ylim = c(-1.3, 1.3))
    if (nrow(ed)) {
      for (i in seq_len(nrow(ed))) {
        arrows(xy[ed$Leader[i], "x"], xy[ed$Leader[i], "y"], xy[ed$Follower[i], "x"], xy[ed$Follower[i], "y"],
               length = 0.08, lwd = 1.5 + 2 * ed$WCD[i] / max(ed$WCD, na.rm = TRUE), col = "grey55")
      }
    }
    cexv <- 1.5 + 3.0 * sqrt(pmax(nd$value, 0)) / max(sqrt(pmax(nd$value, 0)), na.rm = TRUE)
    points(xy$x, xy$y, pch = 21, bg = bg, cex = cexv, lwd = 2)
    text(xy$x, xy$y, labels = nd$rank, cex = 0.95, font = 2)
    text(1.14 * xy$x, 1.14 * xy$y, labels = nd$name, cex = 0.9, font = 2)
    legend("bottomleft", legend = paste("Cluster", sort(unique(nd$carac))),
           pt.bg = cols[seq_along(sort(unique(nd$carac)))], pch = 21, bty = "n", cex = 1.0, text.font = 2)
  }

  output$tbl_ama_ref_top20_check <- renderPrint({
    tryCatch({
      cat("Top 20 check list\n")
      cat("This table is the FLCA-MA-SIL Top20 used by Network, SSplot, and Kano. Numeric values are displayed to 2 decimals.\n")
      z <- .ama705_payload()
      keep <- intersect(c("rank", "name", "value", "value2", "carac", "ssi", "a_i", "b_i", "a_star1"), names(z$nodes))
      print(.ama705_fmt_table(z$nodes[, keep, drop = FALSE], max_rows = 20), row.names = FALSE, right = FALSE)
    }, error = function(e) {
      cat("Top 20 check list\n")
      cat("Run AMA author/journal extraction first. ", conditionMessage(e), "\n", sep = "")
    })
  })

  output$tbl_ama_ref_nodes <- renderPrint({
    tryCatch({
      z <- .ama705_payload()
      keep <- intersect(c("rank", "name", "value", "value2", "carac", "ssi", "a_i", "b_i", "a_star1"), names(z$nodes))
      print(.ama705_fmt_table(z$nodes[, keep, drop = FALSE], max_rows = 20), row.names = FALSE, right = FALSE)
    }, error = function(e) print(data.frame(Message = conditionMessage(e)), row.names = FALSE, right = FALSE))
  })

  output$tbl_ama_ref_edges <- renderPrint({
    tryCatch({
      z <- .ama705_payload()
      print(.ama705_fmt_table(z$edges, max_rows = 120), row.names = FALSE, right = FALSE)
    }, error = function(e) print(data.frame(Message = conditionMessage(e)), row.names = FALSE, right = FALSE))
  })

  output$vn_ama_ref <- visNetwork::renderVisNetwork({
    .ama705_network_widget()
  })

  output$vn_ama_ref_cluster <- visNetwork::renderVisNetwork({
    .ama705_cluster_network_widget()
  })

  output$plt_ama_ref_ss <- renderPlot({
    tryCatch(.ama705_draw_ssplot(font_scale = 1.55),
             error = function(e) { plot.new(); text(0.5, 0.5, paste("SSplot error:", conditionMessage(e)), col = "red", cex = 1.25, font = 2) })
  }, height = 900)

  output$plt_ama_ref_kano <- renderPlot({
    tryCatch(.ama705_draw_kano(),
             error = function(e) { plot.new(); text(0.5, 0.5, paste("Kano plot error:", conditionMessage(e)), col = "red", cex = 1.25, font = 2) })
  }, height = 1250)

  output$dl_ama_ref_network_png <- downloadHandler(
    filename = function() paste0("ama_author_journal_network_", Sys.Date(), ".png"),
    content = function(file) {
      grDevices::png(file, width = 2600, height = 2200, res = 240)
      on.exit(grDevices::dev.off(), add = TRUE)
      .ama705_draw_network_png()
    }
  )

  output$dl_ama_ref_network_html <- downloadHandler(
    filename = function() paste0("ama_author_journal_interactive_network_", Sys.Date(), ".html"),
    content = function(file) {
      widget <- .ama705_network_widget()
      htmlwidgets::saveWidget(widget, file = file, selfcontained = TRUE)
    }
  )

  output$dl_ama_ref_cluster_network_png <- downloadHandler(
    filename = function() paste0("ama_author_journal_cluster_network_", Sys.Date(), ".png"),
    content = function(file) {
      grDevices::png(file, width = 2800, height = 2300, res = 240)
      on.exit(grDevices::dev.off(), add = TRUE)
      .ama705_draw_cluster_network_png()
    }
  )

  output$dl_ama_ref_cluster_network_html <- downloadHandler(
    filename = function() paste0("ama_author_journal_cluster_network_", Sys.Date(), ".html"),
    content = function(file) {
      widget <- .ama705_cluster_network_widget()
      htmlwidgets::saveWidget(widget, file = file, selfcontained = TRUE)
    }
  )

  output$dl_ama_ref_ss_png <- downloadHandler(
    filename = function() paste0("ama_author_journal_ssplot_", Sys.Date(), ".png"),
    content = function(file) {
      grDevices::png(file, width = 3400, height = 2600, res = 260)
      on.exit(grDevices::dev.off(), add = TRUE)
      .ama705_draw_ssplot(font_scale = 1.65)
    }
  )

  output$dl_ama_ref_kano_png <- downloadHandler(
    filename = function() paste0("ama_author_journal_kano_", Sys.Date(), ".png"),
    content = function(file) {
      grDevices::png(file, width = 3200, height = 4200, res = 260)
      on.exit(grDevices::dev.off(), add = TRUE)
      .ama705_draw_kano()
    }
  )



  # ============================================================
  # AMA Author AAC tab: first/last author extraction only
  # ============================================================
  ama_author_aac_status <- reactiveVal("Ready: paste/upload AMA/PubMed/Google Scholar/NCKU Scopus rows, then click Run AMA author AAC (1st+Last only).")
  ama_author_aac_results_val <- reactiveVal(NULL)

  output$ama_author_aac_status <- renderText({
    as.character(ama_author_aac_status() %||% "")
  })

  .ama_author_fl_clean_author <- function(x) {
    x <- trimws(.ama705_chr(x))
    x <- x[!is.na(x) & nzchar(x)]
    x <- x[vapply(x, .ref_valid_node, logical(1))]
    unique(x)
  }


  .ama_author_force_two_clusters_if_one <- function(nd) {
    # Keep AMA Author AAC SSplot/Kano interpretable when FLCA returns only one cluster.
    # This mirrors flca_ms_sil_module.R split_cluster_if_needed(): if only one
    # cluster exists, split the ranked nodes into two value-based halves.
    if (is.null(nd) || !is.data.frame(nd) || nrow(nd) < 2L || !("carac" %in% names(nd))) return(nd)
    car <- trimws(as.character(nd$carac))
    car <- car[!is.na(car) & nzchar(car) & toupper(car) != "NA"]
    if (length(unique(car)) <= 1L) {
      score <- if ("value2" %in% names(nd)) nd$value2 else if ("value" %in% names(nd)) nd$value else rep(0, nrow(nd))
      score <- suppressWarnings(as.numeric(score)); score[!is.finite(score)] <- 0
      v1 <- if ("value" %in% names(nd)) suppressWarnings(as.numeric(nd$value)) else rep(0, nrow(nd))
      v1[!is.finite(v1)] <- 0
      nm <- if ("name" %in% names(nd)) as.character(nd$name) else as.character(seq_len(nrow(nd)))
      ord <- order(-score, -v1, nm)
      split_index <- ceiling(length(ord) / 2)
      nd$carac <- 2L
      nd$carac[ord[seq_len(split_index)]] <- 1L
    }
    nd$carac <- as.integer(factor(as.character(nd$carac), levels = unique(as.character(nd$carac))))
    nd
  }

  .ama_author_first_last_wide <- function(wide) {
    if (is.null(wide) || !is.data.frame(wide) || !nrow(wide)) return(data.frame())
    author_cols <- grep("^Author_", names(wide), value = TRUE)
    if (!length(author_cols)) return(data.frame())
    rows <- lapply(seq_len(nrow(wide)), function(i) {
      a <- unlist(wide[i, author_cols, drop = FALSE], use.names = FALSE)
      a <- .ama_author_fl_clean_author(a)
      first <- if (length(a) >= 1L) a[1] else ""
      last  <- if (length(a) >= 2L) a[length(a)] else ""
      data.frame(reference_no = i, FirstAuthor = first, LastAuthor = last, stringsAsFactors = FALSE, check.names = FALSE)
    })
    out <- do.call(rbind, rows)
    out <- out[nzchar(out$FirstAuthor) | nzchar(out$LastAuthor), , drop = FALSE]
    rownames(out) <- NULL
    out
  }

  .ama_author_first_last_nodes_edges <- function(fl_wide) {
    if (is.null(fl_wide) || !is.data.frame(fl_wide) || !nrow(fl_wide)) {
      stop("No first/last authors were detected.", call. = FALSE)
    }
    row_terms <- lapply(seq_len(nrow(fl_wide)), function(i) {
      unique(.ama_author_fl_clean_author(c(fl_wide$FirstAuthor[i], fl_wide$LastAuthor[i])))
    })
    row_terms <- row_terms[lengths(row_terms) > 0]
    if (!length(row_terms)) stop("No valid first/last-author terms after cleaning.", call. = FALSE)
    tb <- sort(table(unlist(row_terms, use.names = FALSE)), decreasing = TRUE)
    nodes <- data.frame(name = names(tb), value = as.numeric(tb), term_type = "Author",
                        stringsAsFactors = FALSE, check.names = FALSE)
    nodes$value2 <- nodes$value
    nodes$carac <- 1L

    edge_list <- list(); kk <- 0L
    for (i in seq_len(nrow(fl_wide))) {
      fa <- .ama_author_fl_clean_author(fl_wide$FirstAuthor[i])
      la <- .ama_author_fl_clean_author(fl_wide$LastAuthor[i])
      if (length(fa) && length(la) && !identical(fa[1], la[1])) {
        kk <- kk + 1L
        edge_list[[kk]] <- data.frame(Leader = fa[1], Follower = la[1], WCD = 1,
                                      stringsAsFactors = FALSE)
      }
    }
    if (length(edge_list)) {
      ed0 <- dplyr::bind_rows(edge_list)
      edges <- stats::aggregate(WCD ~ Leader + Follower, data = ed0, FUN = sum)
      edges <- edges[order(-edges$WCD, edges$Leader, edges$Follower), , drop = FALSE]
    } else {
      edges <- data.frame(Leader = character(0), Follower = character(0), WCD = numeric(0), stringsAsFactors = FALSE)
    }
    list(nodes = nodes, edges = edges, rows_terms = row_terms)
  }

  .ama_author_aac_run_frequency_top20 <- function(nodes, edges) {
    # PATCH 2026-06-02:
    # For the AMA Author AAC tab, "Top20 by first/last authors" should mean
    # frequency-first Top20 from the cleaned first/last bylines.  FLCA-MA-SIL
    # is still used for cluster labels when available, but the high-frequency
    # authors are no longer dropped by major-sampling.  This is important for
    # one-initial names such as "W Chou" that can be legitimate last authors.
    if (!is.data.frame(nodes) || !nrow(nodes)) stop("No first/last-author nodes available.", call. = FALSE)
    if (!is.data.frame(edges)) edges <- data.frame(Leader = character(0), Follower = character(0), WCD = numeric(0))

    nd_all <- as.data.frame(nodes, stringsAsFactors = FALSE, check.names = FALSE)
    if (!"name" %in% names(nd_all)) stop("nodes must contain name.", call. = FALSE)
    if (!"value" %in% names(nd_all)) nd_all$value <- 1
    nd_all$name <- trimws(as.character(nd_all$name))
    nd_all$value <- suppressWarnings(as.numeric(nd_all$value))
    nd_all$value[!is.finite(nd_all$value)] <- 1
    nd_all <- nd_all[nzchar(nd_all$name), , drop = FALSE]
    nd_all <- nd_all[!duplicated(nd_all$name), , drop = FALSE]
    if (!"term_type" %in% names(nd_all)) nd_all$term_type <- "Author"

    ed_all <- as.data.frame(edges, stringsAsFactors = FALSE, check.names = FALSE)
    if ("follower" %in% names(ed_all) && !("Follower" %in% names(ed_all))) names(ed_all)[names(ed_all) == "follower"] <- "Follower"
    if (!("Leader" %in% names(ed_all)) && "Source" %in% names(ed_all)) names(ed_all)[names(ed_all) == "Source"] <- "Leader"
    if (!("Follower" %in% names(ed_all)) && "Target" %in% names(ed_all)) names(ed_all)[names(ed_all) == "Target"] <- "Follower"
    if (!("WCD" %in% names(ed_all))) ed_all$WCD <- 1
    if (!all(c("Leader", "Follower", "WCD") %in% names(ed_all))) {
      ed_all <- data.frame(Leader = character(0), Follower = character(0), WCD = numeric(0), stringsAsFactors = FALSE)
    }
    ed_all$Leader <- trimws(as.character(ed_all$Leader))
    ed_all$Follower <- trimws(as.character(ed_all$Follower))
    ed_all$WCD <- suppressWarnings(as.numeric(ed_all$WCD))
    ed_all$WCD[!is.finite(ed_all$WCD)] <- 1
    ed_all <- ed_all[nzchar(ed_all$Leader) & nzchar(ed_all$Follower) & ed_all$Leader != ed_all$Follower,
                     c("Leader", "Follower", "WCD"), drop = FALSE]
    ed_all <- ed_all[ed_all$Leader %in% nd_all$name & ed_all$Follower %in% nd_all$name, , drop = FALSE]

    # Incident edge strength as value2.
    if (nrow(ed_all) > 0) {
      strength_df <- rbind(
        data.frame(name = ed_all$Leader, value2 = ed_all$WCD, stringsAsFactors = FALSE),
        data.frame(name = ed_all$Follower, value2 = ed_all$WCD, stringsAsFactors = FALSE)
      )
      strength_df <- stats::aggregate(value2 ~ name, strength_df, sum)
      nd_all$value2 <- strength_df$value2[match(nd_all$name, strength_df$name)]
      nd_all$value2[!is.finite(nd_all$value2) | is.na(nd_all$value2)] <- nd_all$value[!is.finite(nd_all$value2) | is.na(nd_all$value2)]
    } else {
      nd_all$value2 <- nd_all$value
    }

    # FLCA-MA-SIL cluster labels on the full first/last network, if available.
    cluster_map <- NULL
    if (exists("run_flca_ms_sil_internal", mode = "function") && nrow(ed_all) > 0 && nrow(nd_all) >= 3) {
      full_flca <- tryCatch(
        run_flca_ms_sil_internal(nd_all, ed_all,
                                 cfg = list(top_clusters = 5, base_per_cluster = 4, target_n = 20),
                                 verbose = FALSE),
        error = function(e) NULL
      )
      if (!is.null(full_flca) && is.data.frame(full_flca$nodes_full) &&
          all(c("name", "carac") %in% names(full_flca$nodes_full))) {
        cluster_map <- full_flca$nodes_full[, c("name", "carac"), drop = FALSE]
      } else if (!is.null(full_flca) && is.data.frame(full_flca$nodes) &&
                 all(c("name", "carac") %in% names(full_flca$nodes))) {
        cluster_map <- full_flca$nodes[, c("name", "carac"), drop = FALSE]
      }
    }
    if (!is.null(cluster_map) && nrow(cluster_map)) {
      nd_all$carac <- cluster_map$carac[match(nd_all$name, cluster_map$name)]
    }
    if (!"carac" %in% names(nd_all) || all(is.na(nd_all$carac))) nd_all$carac <- NA_integer_

    # Fallback cluster labels by the highest-count leaders.
    if (any(is.na(nd_all$carac))) {
      seeds <- nd_all$name[order(-nd_all$value, -nd_all$value2, nd_all$name)]
      seeds <- utils::head(seeds, min(5L, length(seeds)))
      for (ii in which(is.na(nd_all$carac))) {
        nm <- nd_all$name[ii]
        hit <- ed_all[(ed_all$Leader == nm & ed_all$Follower %in% seeds) |
                        (ed_all$Follower == nm & ed_all$Leader %in% seeds), , drop = FALSE]
        if (nrow(hit) > 0) {
          hit$seed <- ifelse(hit$Leader %in% seeds, hit$Leader, hit$Follower)
          hit <- hit[order(-hit$WCD), , drop = FALSE]
          nd_all$carac[ii] <- match(hit$seed[1], seeds)
        }
      }
      still <- which(is.na(nd_all$carac))
      if (length(still)) nd_all$carac[still] <- ((seq_along(still) - 1L) %% max(1L, length(seeds))) + 1L
    }

    nd_all$carac <- as.integer(factor(as.character(nd_all$carac), levels = unique(as.character(nd_all$carac))))
    nd_all <- .ama_author_force_two_clusters_if_one(nd_all)
    nd_all <- nd_all[order(-nd_all$value, -nd_all$value2, nd_all$name), , drop = FALSE]
    nd <- utils::head(nd_all, 20)
    nd <- .ama_author_force_two_clusters_if_one(nd)
    selected <- as.character(nd$name)
    ed <- ed_all[ed_all$Leader %in% selected & ed_all$Follower %in% selected, , drop = FALSE]

    nd$ssi <- 0; nd$a_i <- 0; nd$b_i <- 0; nd$a_star1 <- 0
    if (exists("compute_silhouette_df", mode = "function") && nrow(ed) > 0 && length(unique(nd$carac)) >= 2) {
      sil_df <- tryCatch({
        ed2 <- ed
        if ("Follower" %in% names(ed2) && !("follower" %in% names(ed2))) ed2$follower <- ed2$Follower
        compute_silhouette_df(nd, ed2)
      }, error = function(e) NULL)
      if (is.data.frame(sil_df) && nrow(sil_df) && "name" %in% names(sil_df)) {
        sx <- match(nd$name, as.character(sil_df$name))
        for (cc in c("ssi", "sil_width", "a_i", "b_i", "a_star1", "a_star")) {
          if (cc %in% names(sil_df)) {
            val <- suppressWarnings(as.numeric(sil_df[[cc]][sx]))
            val[!is.finite(val)] <- 0
            if (cc == "sil_width") nd$ssi <- val
            else if (cc == "a_star") nd$a_star1 <- val
            else nd[[cc]] <- val
          }
        }
      }
    }
    nd$SSi <- nd$ssi
    nd$a_star <- nd$a_star1
    nd$rank <- seq_len(nrow(nd))
    rownames(nd) <- NULL
    rownames(ed) <- NULL
    list(nodes = nd, edges = ed, engine = "FLCA-MA-SIL cluster labels + first/last frequency Top20; one-cluster result forced to 2 clusters")
  }

  .ama_author_aac_payload <- function() {
    r <- ama_author_aac_results_val()
    if (is.null(r) || is.null(r$nodes) || !is.data.frame(r$nodes)) {
      stop("No AMA Author AAC result yet. Paste/upload normalized records, then click Run AMA author AAC (1st+Last only).", call. = FALSE)
    }
    p <- .ref_sanitize_plot_payload(r$nodes, r$edges)
    nd <- p$nodes
    ed <- p$edges
    if (!"rank" %in% names(nd)) nd$rank <- seq_len(nrow(nd))
    if (!"value2" %in% names(nd)) nd$value2 <- nd$value
    if (!"ssi" %in% names(nd)) nd$ssi <- 0
    if (!"a_i" %in% names(nd)) nd$a_i <- 0
    if (!"b_i" %in% names(nd)) nd$b_i <- 0
    if (!"a_star1" %in% names(nd)) nd$a_star1 <- 0
    nd$value <- .ama705_num(nd$value, 1)
    nd$value2 <- .ama705_num(nd$value2, nd$value)
    nd$carac <- as.integer(.ama705_num(nd$carac, 1))
    nd$ssi <- .ama705_num(nd$ssi, 0)
    nd$a_i <- .ama705_num(nd$a_i, 0)
    nd$b_i <- .ama705_num(nd$b_i, 0)
    nd$a_star1 <- .ama705_num(nd$a_star1, 0)
    nd <- nd[order(-nd$value, -nd$value2, nd$name), , drop = FALSE]
    nd <- utils::head(nd, 20)
    nd$rank <- seq_len(nrow(nd))
    ed <- ed[ed$Leader %in% nd$name & ed$Follower %in% nd$name, , drop = FALSE]
    list(nodes = nd, edges = ed, raw = r)
  }

  .ama_author_aac_network_widget <- function() {
    z <- .ama_author_aac_payload()
    nd <- z$nodes
    ed <- .ama705_one_link_edges(nd, z$edges)
    mx <- max(nd$value, na.rm = TRUE); if (!is.finite(mx) || mx <= 0) mx <- 1
    groups <- paste0("Cluster ", nd$carac)
    vn_nodes <- data.frame(
      id = .ama705_chr(nd$name), label = .ama705_chr(nd$name), group = .ama705_chr(groups),
      value = as.numeric(pmax(12, 20 + 55 * sqrt(pmax(nd$value, 0)) / sqrt(mx))),
      title = .ama705_escape(paste0("Rank: ", nd$rank, "<br>Name: ", nd$name,
                                    "<br>Cluster: ", nd$carac,
                                    "<br>First/last count: ", sprintf("%.2f", nd$value),
                                    "<br>value2: ", sprintf("%.2f", nd$value2),
                                    "<br>SS: ", sprintf("%.2f", nd$ssi))),
      stringsAsFactors = FALSE
    )
    vn_edges <- data.frame(
      from = .ama705_chr(ed$Leader), to = .ama705_chr(ed$Follower),
      value = as.numeric(pmax(1, ed$WCD)), label = .ama705_chr(sprintf("%.0f", ed$WCD)),
      title = .ama705_escape(paste0(ed$Leader, " -> ", ed$Follower, "<br>WCD=", sprintf("%.2f", ed$WCD))),
      arrows = rep("to", nrow(ed)), stringsAsFactors = FALSE
    )
    visNetwork::visNetwork(vn_nodes, vn_edges, width = "100%", height = "760px",
                           main = "AMA Author AAC: FLCA-MA-SIL Top20 from 1st/last authors only") %>%
      visNetwork::visNodes(shape = "dot", font = list(size = 26, face = "bold", color = "#111111", strokeWidth = 3)) %>%
      visNetwork::visEdges(smooth = list(enabled = TRUE, type = "dynamic"),
                            font = list(size = 18, face = "bold", align = "middle"),
                            color = list(color = "#888888", highlight = "#d62728")) %>%
      visNetwork::visOptions(highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE),
                              nodesIdSelection = TRUE,
                              selectedBy = list(variable = "group", multiple = TRUE)) %>%
      visNetwork::visPhysics(stabilization = TRUE, solver = "forceAtlas2Based") %>%
      visNetwork::visLayout(randomSeed = 807)
  }

  .ama_author_aac_cluster_network_widget <- function() {
    z <- .ama_author_aac_payload()
    nd <- .ama705_cluster_positions(z$nodes)
    ed <- .ama705_cluster_leader_edges(nd, z$edges)
    mx <- max(nd$value, na.rm = TRUE); if (!is.finite(mx) || mx <= 0) mx <- 1
    groups <- paste0("Cluster ", nd$carac)
    vn_nodes <- data.frame(
      id = .ama705_chr(nd$name), label = paste0(nd$rank, ". ", .ama705_chr(nd$name)), group = .ama705_chr(groups),
      value = as.numeric(pmax(18, 26 + 66 * sqrt(pmax(nd$value, 0)) / sqrt(mx) + ifelse(nd$cluster_leader, 14, 0))),
      x = as.numeric(nd$x), y = as.numeric(nd$y), fixed = TRUE,
      borderWidth = ifelse(nd$cluster_leader, 5, 2),
      title = .ama705_escape(paste0("Rank: ", nd$rank, "<br>Name: ", nd$name,
                                    "<br>Cluster: ", nd$carac,
                                    "<br>Role: ", ifelse(nd$cluster_leader, "cluster leader", "follower"),
                                    "<br>First/last count: ", sprintf("%.2f", nd$value),
                                    "<br>value2: ", sprintf("%.2f", nd$value2),
                                    "<br>SS: ", sprintf("%.2f", nd$ssi))),
      stringsAsFactors = FALSE
    )
    vn_edges <- data.frame(
      from = .ama705_chr(ed$Follower), to = .ama705_chr(ed$Leader),
      value = as.numeric(pmax(1, ed$WCD)), label = .ama705_chr(sprintf("%.0f", ed$WCD)),
      title = .ama705_escape(paste0("Follower: ", ed$Follower, "<br>Leader: ", ed$Leader,
                                    "<br>Cluster: ", ed$carac, "<br>WCD=", sprintf("%.2f", ed$WCD),
                                    "<br>", ed$edge_type)),
      arrows = rep("to", nrow(ed)), stringsAsFactors = FALSE
    )
    visNetwork::visNetwork(vn_nodes, vn_edges, width = "100%", height = "760px",
                           main = "Clustered network: 1st/last-author Top20; each follower points to its cluster leader") %>%
      visNetwork::visNodes(shape = "dot", font = list(size = 24, face = "bold", color = "#111111", strokeWidth = 3)) %>%
      visNetwork::visEdges(smooth = list(enabled = TRUE, type = "curvedCW", roundness = 0.18),
                            font = list(size = 18, face = "bold", align = "middle"),
                            color = list(color = "#777777", highlight = "#d62728")) %>%
      visNetwork::visOptions(highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE),
                              nodesIdSelection = TRUE,
                              selectedBy = list(variable = "group", multiple = TRUE)) %>%
      visNetwork::visPhysics(enabled = FALSE) %>%
      visNetwork::visLayout(randomSeed = 808)
  }

  .ama_author_aac_ss_inputs <- function() {
    z <- .ama_author_aac_payload()
    nd <- z$nodes
    ed1 <- .ama705_one_link_edges(nd, z$edges)
    nd$sil_width <- nd$ssi
    nd$wsel <- NA_real_
    nd$role <- NA_character_
    nd$neighbor_name <- nd$name
    nd$neighborC <- nd$carac
    if (nrow(ed1)) {
      for (i in seq_len(nrow(ed1))) {
        j <- match(ed1$Follower[i], nd$name)
        k <- match(ed1$Leader[i], nd$name)
        if (!is.na(j)) {
          nd$wsel[j] <- ed1$WCD[i]
          nd$role[j] <- "follower"
          nd$neighbor_name[j] <- ed1$Leader[i]
          nd$neighborC[j] <- if (!is.na(k)) nd$carac[k] else nd$carac[j]
        }
      }
    }
    leaders <- tapply(seq_len(nrow(nd)), nd$carac, function(ii) ii[order(-nd$value2[ii], -nd$value[ii], nd$rank[ii])][1])
    nd$role[as.integer(leaders)] <- "leader"
    results <- .ama705_modularity_results(nd, z$edges)
    list(sil_df = nd, nodes0 = nd, nodes = nd, results = results)
  }

  .ama_author_aac_draw_ssplot <- function(font_scale = 1.55) {
    inp <- .ama_author_aac_ss_inputs()
    if (!exists("render_panel", mode = "function")) stop("render_panel() not loaded from renderSSplot.R.", call. = FALSE)
    .with_clean_ssplot_aac_labels({
      render_panel(sil_df = inp$sil_df, nodes0 = inp$nodes0, results = inp$results,
                   nodes = inp$nodes, top_n = nrow(inp$sil_df), font_scale = font_scale,
                   aac_side = "left", neighbor_side = "right", neighbor_on_bar = TRUE,
                   footer_label = "AMA Author AAC: first/last authors only")
    })
  }

  .ama_author_aac_draw_kano <- function() {
    z <- .ama_author_aac_payload()
    nd <- z$nodes
    ed <- .ama705_one_link_edges(nd, z$edges)
    ok <- FALSE
    if (exists("plot_kano_real", mode = "function")) {
      ok <- tryCatch({
        obj <- plot_kano_real(nodes = nd, edges = ed, title_txt = "Kano plot: 1st/last-author AAC Top20")
        .ama705_print_kano_object(obj); TRUE
      }, error = function(e) FALSE)
    }
    if (!ok && exists("kano_plot", mode = "function")) {
      ok <- tryCatch({
        obj <- kano_plot(nd, ed, title_txt = "Kano plot: 1st/last-author AAC Top20")
        .ama705_print_kano_object(obj); TRUE
      }, error = function(e) FALSE)
    }
    if (ok) return(invisible(TRUE))
    nd$cluster <- factor(nd$carac)
    p <- ggplot2::ggplot(nd, ggplot2::aes(x = value2, y = value, size = value, color = cluster, label = name)) +
      ggplot2::geom_vline(xintercept = stats::median(nd$value2, na.rm = TRUE), linetype = "dashed", linewidth = 1.2, color = "red") +
      ggplot2::geom_hline(yintercept = stats::median(nd$value, na.rm = TRUE), linetype = "dashed", linewidth = 1.2, color = "red") +
      ggplot2::geom_point(alpha = 0.85) +
      ggplot2::geom_text(fontface = "bold", size = 5.5, vjust = -1.05, show.legend = FALSE) +
      ggplot2::scale_size_continuous(range = c(5, 18)) +
      ggplot2::labs(title = "Kano plot: 1st/last-author AAC Top20",
                    subtitle = "x = value2 (edge strength); y = first/last-author count; color = FLCA cluster",
                    x = "value2", y = "1st/last-author count", color = "Cluster", size = "Count") +
      ggplot2::theme_bw(base_size = 21) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 28),
                     plot.subtitle = ggplot2::element_text(face = "bold", size = 18),
                     axis.title = ggplot2::element_text(face = "bold", size = 23),
                     axis.text = ggplot2::element_text(face = "bold", size = 17),
                     legend.title = ggplot2::element_text(face = "bold", size = 19),
                     legend.text = ggplot2::element_text(face = "bold", size = 16),
                     aspect.ratio = 1.35,
                     plot.margin = ggplot2::margin(18, 24, 18, 24))
    print(p)
    invisible(TRUE)
  }

  .ama_author_aac_draw_network_png <- function() {
    z <- .ama_author_aac_payload(); nd <- z$nodes; ed <- .ama705_one_link_edges(nd, z$edges)
    n <- nrow(nd); theta <- seq(pi/2, pi/2 - 2*pi + 2*pi/n, length.out = n)
    xy <- data.frame(name = nd$name, x = cos(theta), y = sin(theta), stringsAsFactors = FALSE); rownames(xy) <- xy$name
    cols <- grDevices::hcl.colors(max(3, length(unique(nd$carac))), "Dark 3")
    bg <- cols[as.integer(factor(nd$carac))]
    par(mar = c(1, 1, 4, 1))
    plot(xy$x, xy$y, type = "n", axes = FALSE, xlab = "", ylab = "",
         main = "AMA Author AAC network PNG: 1st/last-author Top20", cex.main = 1.7, font.main = 2,
         xlim = c(-1.4, 1.4), ylim = c(-1.3, 1.3))
    if (nrow(ed)) for (i in seq_len(nrow(ed))) arrows(xy[ed$Leader[i], "x"], xy[ed$Leader[i], "y"], xy[ed$Follower[i], "x"], xy[ed$Follower[i], "y"], length = 0.08, lwd = 1.5 + 2 * ed$WCD[i] / max(ed$WCD, na.rm = TRUE), col = "grey55")
    cexv <- 1.5 + 3.0 * sqrt(pmax(nd$value, 0)) / max(sqrt(pmax(nd$value, 0)), na.rm = TRUE)
    points(xy$x, xy$y, pch = 21, bg = bg, cex = cexv, lwd = 2)
    text(xy$x, xy$y, labels = nd$rank, cex = 0.95, font = 2)
    text(1.14 * xy$x, 1.14 * xy$y, labels = nd$name, cex = 0.9, font = 2)
    legend("bottomleft", legend = paste("Cluster", sort(unique(nd$carac))), pt.bg = cols[seq_along(sort(unique(nd$carac)))], pch = 21, bty = "n", cex = 1.0, text.font = 2)
  }

  .ama_author_aac_draw_cluster_network_png <- function() {
    z <- .ama_author_aac_payload(); nd <- .ama705_cluster_positions(z$nodes); ed <- .ama705_cluster_leader_edges(nd, z$edges)
    cols <- grDevices::hcl.colors(max(3, length(unique(nd$carac))), "Dark 3")
    clv <- sort(unique(nd$carac)); bg <- cols[as.integer(factor(nd$carac, levels = clv))]
    xlim <- range(nd$x, na.rm = TRUE) + c(-240, 240); ylim <- range(nd$y, na.rm = TRUE) + c(-200, 200)
    par(mar = c(1, 1, 4, 1))
    plot(nd$x, nd$y, type = "n", axes = FALSE, xlab = "", ylab = "",
         main = "Clustered AMA Author AAC network: follower -> cluster leader", cex.main = 1.65, font.main = 2,
         xlim = xlim, ylim = ylim, asp = 1)
    for (cc in clv) {
      sub <- nd[nd$carac == cc, , drop = FALSE]
      cx <- sub$x[which(sub$cluster_leader)[1]]; cy <- sub$y[which(sub$cluster_leader)[1]]
      if (!is.finite(cx) || !is.finite(cy)) { cx <- mean(sub$x, na.rm = TRUE); cy <- mean(sub$y, na.rm = TRUE) }
      r <- max(sqrt((sub$x - cx)^2 + (sub$y - cy)^2), 80, na.rm = TRUE) + 75
      symbols(cx, cy, circles = r, inches = FALSE, add = TRUE,
              bg = grDevices::adjustcolor(cols[match(cc, clv)], alpha.f = 0.10),
              fg = grDevices::adjustcolor(cols[match(cc, clv)], alpha.f = 0.70), lwd = 2)
      text(cx, cy + r + 35, labels = paste("Cluster", cc), cex = 1.15, font = 2, col = cols[match(cc, clv)])
    }
    if (nrow(ed)) {
      mxw <- max(ed$WCD, na.rm = TRUE); if (!is.finite(mxw) || mxw <= 0) mxw <- 1
      for (i in seq_len(nrow(ed))) {
        a <- match(ed$Leader[i], nd$name); b <- match(ed$Follower[i], nd$name)
        if (!is.na(a) && !is.na(b)) arrows(nd$x[b], nd$y[b], nd$x[a], nd$y[a], length = 0.08, lwd = 1.5 + 2 * ed$WCD[i] / mxw, col = "grey45")
      }
    }
    cexv <- 1.5 + 3.0 * sqrt(pmax(nd$value, 0)) / max(sqrt(pmax(nd$value, 0)), na.rm = TRUE)
    cexv[nd$cluster_leader] <- cexv[nd$cluster_leader] + 0.8
    points(nd$x, nd$y, pch = 21, bg = bg, cex = cexv, lwd = ifelse(nd$cluster_leader, 4, 2))
    text(nd$x, nd$y, labels = nd$rank, cex = 0.95, font = 2)
    text(nd$x + 28, nd$y - 28, labels = nd$name, cex = 0.88, font = 2, adj = c(0, 1))
    legend("bottomleft", legend = paste("Cluster", clv), pt.bg = cols[seq_along(clv)], pch = 21, bty = "n", cex = 1.0, text.font = 2)
    legend("topright", legend = c("larger border = cluster leader", "arrow: follower -> leader"), bty = "n", cex = 1.0, text.font = 2)
  }

  .ama_author_aac_highlights_df <- function() {
    z <- .ama_author_aac_payload(); nd <- z$nodes; r <- z$raw
    aac_value <- .compute_AAC(nd$value); if (!is.finite(aac_value)) aac_value <- 0
    aac_value2 <- .compute_AAC(nd$value2); if (!is.finite(aac_value2)) aac_value2 <- 0
    aac_ss <- .compute_AAC(nd$ssi); if (!is.finite(aac_ss)) aac_ss <- 0
    aac_astar <- .compute_AAC(nd$a_star1); if (!is.finite(aac_astar)) aac_astar <- 0
    art <- r$gs_articles %||% data.frame()
    ps <- r$gs_profile_summary %||% data.frame()
    if (is.data.frame(ps) && nrow(ps)) {
      cit_all <- ps$All[match("Google Scholar profile citations", ps$Metric)]
      h_all <- ps$All[match("Google Scholar profile h-index", ps$Metric)]
      i10_all <- ps$All[match("Google Scholar profile i10-index", ps$Metric)]
      metric_source <- "Google Scholar profile-reported"
    } else if (is.data.frame(art) && nrow(art)) {
      gm <- .ref_google_scholar_metrics(art)
      h_all <- gm$All[match("h-index", gm$metric)]
      i10_all <- gm$All[match("i10-index", gm$metric)]
      cit_all <- gm$All[match("Citations", gm$metric)]
      metric_source <- "Parsed normalized references"
    } else {
      h_all <- NA; i10_all <- NA; cit_all <- NA; metric_source <- "Not detected"
    }
    parsed_cit <- parsed_h <- parsed_i10 <- NA
    fl_top_author <- fl_top_h <- fl_top_cit <- NA
    if (is.data.frame(art) && nrow(art)) {
      parsed_cites <- suppressWarnings(as.integer(art$citations)); parsed_cites[!is.finite(parsed_cites) | is.na(parsed_cites)] <- 0L
      parsed_cit <- sum(parsed_cites, na.rm = TRUE)
      parsed_h <- .ref_gs_h_index(parsed_cites)
      parsed_i10 <- sum(parsed_cites >= 10L, na.rm = TRUE)
      fl_hx <- tryCatch(.ref_google_scholar_author_hindex(art, basis = "first_last", since_year = r$gs_since_year %||% NULL), error = function(e) data.frame())
      if (is.data.frame(fl_hx) && nrow(fl_hx)) {
        fl_top_author <- as.character(fl_hx$author[1])
        fl_top_h <- as.character(fl_hx$h_index[1])
        fl_top_cit <- as.character(fl_hx$citations[1])
      }
    }
    data.frame(
      Metric = c("Citation metric source", "Profile-reported citations", "Profile-reported h-index", "Profile-reported i10-index",
                 "Parsed normalized-reference citations", "Parsed normalized-reference h-index", "Parsed normalized-reference i10-index",
                 "Top 1st/last-author reference h-index", "Top 1st/last-author reference author", "Top 1st/last-author reference citations",
                 "Author AAC by first/last count", "Author AAC by edge strength (value2)",
                 "AAC by SS", "AAC by a*", "Mean SS", "Clusters", "Top20 nodes", "First-last edges", "Engine"),
      Value = c(as.character(metric_source), as.character(cit_all), as.character(h_all), as.character(i10_all),
                as.character(parsed_cit), as.character(parsed_h), as.character(parsed_i10),
                as.character(fl_top_h), as.character(fl_top_author), as.character(fl_top_cit),
                sprintf("%.3f", aac_value), sprintf("%.3f", aac_value2),
                sprintf("%.3f", aac_ss), sprintf("%.3f", aac_astar),
                sprintf("%.3f", mean(nd$ssi, na.rm = TRUE)),
                as.character(length(unique(nd$carac))), as.character(nrow(nd)), as.character(nrow(z$edges)),
                as.character(r$engine %||% "")),
      stringsAsFactors = FALSE, check.names = FALSE
    )
  }

  observeEvent(input$run_ama_author_aac, {
    .ama_console_log("Button pressed", "Run AMA author AAC (1st+Last only)", .button = "AMA Author AAC")
    ama_author_aac_status("Processing: extracting first/last authors and running FLCA-MA-SIL...")
    ama_author_aac_results_val(NULL)
    tryCatch({
      shiny::withProgress(message = "Running AMA Author AAC", value = 0, {
        shiny::incProgress(0.10, detail = "Reading input")
        txt <- .get_ama_input_text(ama_author_aac_status)
        if (is.null(txt) || !nzchar(trimws(txt))) stop("No uploaded or pasted AMA/PubMed/Google Scholar/NCKU Scopus text was found.", call. = FALSE)
        selected_sources_aac <- .ama_selected_ref_sources()
        .ama_console_log("Parsers", paste(selected_sources_aac, collapse = ", "), .button = "AMA Author AAC")
        refs <- .ama_split_refs_selected(txt, selected_sources_aac)
        if (!length(refs)) stop("No usable selected-source normalized records were detected.", call. = FALSE)
        shiny::incProgress(0.20, detail = "Parsing authors")
        gs_profile_summary <- if ("google" %in% selected_sources_aac) tryCatch(.ref_google_scholar_profile_summary(txt), error = function(e) data.frame()) else data.frame()
        gs_articles <- .ama_parse_articles_selected(txt, selected_sources_aac)
        if (is.data.frame(gs_articles) && nrow(gs_articles) > 0) {
          wide0 <- .ref_google_scholar_articles_to_wide(gs_articles)
        } else {
          if (!("ama" %in% selected_sources_aac)) stop("No normalized Google Scholar/NCKU rows were extracted, and AMA/PubMed parser is unchecked.", call. = FALSE)
          wide0 <- .ref_to_wide_from_text_strict(txt)
          wide_gs_force <- if ("google" %in% selected_sources_aac) tryCatch(.ref_to_wide_google_scholar_force(txt), error = function(e) data.frame()) else data.frame()
          if (is.data.frame(wide_gs_force) && nrow(wide_gs_force) >= 2L) wide0 <- wide_gs_force
        }
        fl_wide <- .ama_author_first_last_wide(wide0)
        if (!is.data.frame(fl_wide) || !nrow(fl_wide)) stop("No first/last authors were detected.", call. = FALSE)
        shiny::incProgress(0.20, detail = "Building first-last author network")
        ne <- .ama_author_first_last_nodes_edges(fl_wide)
        shiny::incProgress(0.30, detail = "Running FLCA-MA-SIL Top20")
        fl <- .ama_author_aac_run_frequency_top20(ne$nodes, ne$edges)
        payload <- .ref_sanitize_plot_payload(fl$nodes, fl$edges)
        nd <- payload$nodes; ed <- payload$edges
        if (!is.data.frame(nd) || !nrow(nd)) stop("FLCA-MA-SIL returned no first/last-author Top20 nodes.", call. = FALSE)
        ama_author_aac_results_val(list(
          refs = refs, wide = fl_wide, nodes = nd, edges = ed, engine = fl$engine,
          raw_nodes = ne$nodes, raw_edges = ne$edges, gs_articles = gs_articles,
          gs_profile_summary = if (exists("gs_profile_summary", inherits = FALSE)) gs_profile_summary else data.frame(),
          gs_since_year = as.integer(format(Sys.Date(), "%Y")) - 5L
        ))
        ama_author_aac_status(paste0("Completed via ", fl$engine, ": parsed ", length(refs),
                                     " references; ", nrow(fl_wide), " first/last-author rows; ",
                                     nrow(nd), " Top20 nodes; ", nrow(ed), " edges."))
        showNotification("AMA Author AAC completed: first/last authors only.", type = "message")
        updateTabsetPanel(session, "main_tabs", selected = "AMA Author AAC")
        shiny::incProgress(0.20, detail = "Done")
      })
    }, error = function(e) {
      .ama_console_log("Error", conditionMessage(e), .button = "AMA Author AAC")
      ama_author_aac_status(paste("Error:", conditionMessage(e)))
      showNotification(paste("AMA Author AAC failed:", conditionMessage(e)), type = "error", duration = 10)
      ama_author_aac_results_val(NULL)
    })
  }, ignoreInit = TRUE)

  output$tbl_ama_author_aac_highlights <- DT::renderDT({
    df <- tryCatch(.ama_author_aac_highlights_df(), error = function(e) {
      data.frame(Metric = "Status", Value = "Paste/upload normalized records and click Run AMA author AAC (1st+Last only).", stringsAsFactors = FALSE, check.names = FALSE)
    })
    DT::datatable(df, rownames = FALSE, options = list(pageLength = 12, dom = "t", scrollX = TRUE))
  })


  output$tbl_ama_author_aac_all_reference_hindex <- DT::renderDT({
    r <- ama_author_aac_results_val()
    art <- if (!is.null(r)) r$gs_articles %||% data.frame() else data.frame()
    df <- .ref_google_scholar_all_reference_hindex(art, profile_summary = (if (!is.null(r)) r$gs_profile_summary %||% data.frame() else data.frame()))
    DT::datatable(df, rownames = FALSE,
                  options = list(dom = "t", scrollX = TRUE))
  })

  output$tbl_ama_author_aac_first_last_reference_hindex <- DT::renderDT({
    r <- ama_author_aac_results_val()
    art <- if (!is.null(r)) r$gs_articles %||% data.frame() else data.frame()
    sy <- if (!is.null(r)) r$gs_since_year %||% (as.integer(format(Sys.Date(), "%Y")) - 5L) else (as.integer(format(Sys.Date(), "%Y")) - 5L)
    df <- if (is.data.frame(art) && nrow(art)) {
      .ref_google_scholar_author_hindex(art, basis = "first_last", since_year = sy)
    } else {
      .ref_empty_author_hindex_message()
    }
    DT::datatable(df, rownames = FALSE, options = list(pageLength = 15, scrollX = TRUE))
  })

  output$tbl_ama_author_aac_top20_check <- renderPrint({
    tryCatch({
      z <- .ama_author_aac_payload()
      keep <- intersect(c("rank", "name", "value", "value2", "carac", "ssi", "a_i", "b_i", "a_star1"), names(z$nodes))
      cat("Top20 first/last-author check list\n")
      cat("value = first/last-author occurrence count; value2 = incident edge strength; carac = FLCA cluster.\n")
      print(.ama705_fmt_table(z$nodes[, keep, drop = FALSE], max_rows = 20), row.names = FALSE, right = FALSE)
    }, error = function(e) cat("No AMA Author AAC result yet. ", conditionMessage(e), "\n", sep = ""))
  })

  output$tbl_ama_author_aac_nodes <- renderPrint({
    tryCatch({ z <- .ama_author_aac_payload(); keep <- intersect(c("rank", "name", "value", "value2", "carac", "ssi", "a_i", "b_i", "a_star1"), names(z$nodes)); print(.ama705_fmt_table(z$nodes[, keep, drop = FALSE], max_rows = 20), row.names = FALSE, right = FALSE) },
             error = function(e) print(data.frame(Message = conditionMessage(e)), row.names = FALSE, right = FALSE))
  })
  output$tbl_ama_author_aac_edges <- renderPrint({
    tryCatch({ z <- .ama_author_aac_payload(); print(.ama705_fmt_table(z$edges, max_rows = 120), row.names = FALSE, right = FALSE) },
             error = function(e) print(data.frame(Message = conditionMessage(e)), row.names = FALSE, right = FALSE))
  })
  output$tbl_ama_author_aac_wide <- renderPrint({
    r <- ama_author_aac_results_val()
    if (is.null(r) || !is.data.frame(r$wide)) {
      print(data.frame(Message = "No AMA Author AAC result yet. Paste/upload normalized records, then click Run AMA author AAC (1st+Last only)."), row.names = FALSE, right = FALSE)
    } else {
      print(.ama705_fmt_table(r$wide, max_rows = 120), row.names = FALSE, right = FALSE)
    }
  })

  output$vn_ama_author_aac <- visNetwork::renderVisNetwork({ .ama_author_aac_network_widget() })
  output$vn_ama_author_aac_cluster <- visNetwork::renderVisNetwork({ .ama_author_aac_cluster_network_widget() })
  output$plt_ama_author_aac_ss <- renderPlot({
    tryCatch(.ama_author_aac_draw_ssplot(font_scale = 1.55),
             error = function(e) { plot.new(); text(0.5, 0.5, paste("SSplot error:", conditionMessage(e)), col = "red", cex = 1.25, font = 2) })
  }, height = 900)
  output$plt_ama_author_aac_kano <- renderPlot({
    tryCatch(.ama_author_aac_draw_kano(),
             error = function(e) { plot.new(); text(0.5, 0.5, paste("Kano plot error:", conditionMessage(e)), col = "red", cex = 1.25, font = 2) })
  }, height = 1250)

  output$dl_ama_author_aac_network_png <- downloadHandler(
    filename = function() paste0("ama_author_aac_network_", Sys.Date(), ".png"),
    content = function(file) { grDevices::png(file, width = 2600, height = 2200, res = 240); on.exit(grDevices::dev.off(), add = TRUE); .ama_author_aac_draw_network_png() }
  )
  output$dl_ama_author_aac_network_html <- downloadHandler(
    filename = function() paste0("ama_author_aac_interactive_network_", Sys.Date(), ".html"),
    content = function(file) { htmlwidgets::saveWidget(.ama_author_aac_network_widget(), file = file, selfcontained = TRUE) }
  )
  output$dl_ama_author_aac_cluster_network_png <- downloadHandler(
    filename = function() paste0("ama_author_aac_cluster_network_", Sys.Date(), ".png"),
    content = function(file) { grDevices::png(file, width = 2800, height = 2300, res = 240); on.exit(grDevices::dev.off(), add = TRUE); .ama_author_aac_draw_cluster_network_png() }
  )
  output$dl_ama_author_aac_cluster_network_html <- downloadHandler(
    filename = function() paste0("ama_author_aac_cluster_network_", Sys.Date(), ".html"),
    content = function(file) { htmlwidgets::saveWidget(.ama_author_aac_cluster_network_widget(), file = file, selfcontained = TRUE) }
  )
  output$dl_ama_author_aac_ss_png <- downloadHandler(
    filename = function() paste0("ama_author_aac_ssplot_", Sys.Date(), ".png"),
    content = function(file) { grDevices::png(file, width = 3400, height = 2600, res = 260); on.exit(grDevices::dev.off(), add = TRUE); .ama_author_aac_draw_ssplot(font_scale = 1.65) }
  )
  output$dl_ama_author_aac_kano_png <- downloadHandler(
    filename = function() paste0("ama_author_aac_kano_", Sys.Date(), ".png"),
    content = function(file) { grDevices::png(file, width = 3200, height = 4200, res = 260); on.exit(grDevices::dev.off(), add = TRUE); .ama_author_aac_draw_kano() }
  )
  output$dl_ama_author_aac_nodes <- downloadHandler(
    filename = function() paste0("ama_author_aac_nodes_", Sys.Date(), ".csv"),
    content = function(file) { z <- .ama_author_aac_payload(); utils::write.csv(z$nodes, file, row.names = FALSE, fileEncoding = "UTF-8") }
  )
  output$dl_ama_author_aac_edges <- downloadHandler(
    filename = function() paste0("ama_author_aac_edges_", Sys.Date(), ".csv"),
    content = function(file) { z <- .ama_author_aac_payload(); utils::write.csv(z$edges, file, row.names = FALSE, fileEncoding = "UTF-8") }
  )
  output$dl_ama_author_aac_wide <- downloadHandler(
    filename = function() paste0("ama_author_aac_first_last_wide_", Sys.Date(), ".csv"),
    content = function(file) { r <- ama_author_aac_results_val(); utils::write.csv(r$wide, file, row.names = FALSE, fileEncoding = "UTF-8") }
  )

  term_rx <- reactive({
    if (isTRUE(free_demo_run())) demo_term else input$term
  })
 print("SERVER STARTED")
 message("SERVER STARTED: ", Sys.time())
  output$ip_status <- renderPrint(list.files())
  rv <- reactiveValues(
    ip_access_type = "checking",
    ip_addr = "..."
  )

  # ---- BCa/EIS manual recompute stores ----
  # Use explicit observeEvent + reactiveVal stores. This makes the Recompute button
  # update outputs immediately and avoids silent eventReactive timing issues.
  bca_result_store <- reactiveVal(NULL)
  top1_result_store <- reactiveVal(NULL)
  bca_values_store <- reactiveVal(NULL)
  bca_source_store <- reactiveVal("Author analysis")
  bca_status_store <- reactiveVal("BCa/EIS: not manually recomputed yet. Use Author analysis, or paste values and press Recompute.")

  output$txt_bca_run_status <- renderText({ bca_status_store() })

  observeEvent(input$bca_recompute_manual, {
    x <- .parse_bca_values(input$bca_values_text)
    if (length(x) < 5) {
      bca_status_store(paste0("Manual recompute stopped: only ", length(x), " positive values found. At least 5 values are required; at least 3 followers are required for BCa/EIS."))
      showNotification("Manual BCa/EIS stopped: at least 5 positive values are required.", type = "error", duration = 6)
      return(invisible(NULL))
    }

    bca_status_store(paste0("Manual BCa/EIS running... n=", length(x), ", R=", input$bca_R, "."))

    out <- shiny::withProgress(message = "Running manual BCa/EIS", value = 0, {
      res <- scan_aac_bca(
        x = x,
        R = input$bca_R,
        cutoff = input$bca_cutoff,
        conf = input$bca_conf,
        progress = function(i, total) {
          shiny::incProgress(0.75 / max(total, 1), detail = paste0("BCa/EIS point ", i, " of ", total, " | R=", input$bca_R))
        }
      )
      shiny::incProgress(0.10, detail = "Computing Top-1 conditional BCa/EIS")
      top1 <- scan_topk_bca_eis(
        x = x,
        k = input$bca_top_k,
        R = input$bca_R,
        cutoff_or = input$bca_cutoff,
        conf = input$bca_conf
      )
      shiny::incProgress(0.15, detail = "Finalizing tables and plots")
      list(res = res, top1 = top1)
    })

    out$res$EIS_point <- suppressWarnings(as.numeric(out$res$BCa_lower)) / input$bca_cutoff
    out$res$AAC_candidate <- ifelse(is.finite(out$res$AAC) & out$res$AAC > input$bca_cutoff, "AAC_OR > cutoff", "AAC_OR <= cutoff")
    out$res$source <- "Editable values"
    out$top1$source <- "Editable values"

    bca_values_store(x)
    bca_result_store(out$res)
    top1_result_store(out$top1)
    bca_source_store("Editable values")
    bca_status_store(paste0("Manual BCa/EIS completed: n=", length(x), ", rows=", nrow(out$res), ", source=Editable values."))
    showNotification(paste0("Manual BCa/EIS completed: ", nrow(out$res), " rows."), type = "message", duration = 5)
  }, ignoreInit = TRUE)

  bca_manual_values <- eventReactive(input$bca_recompute_manual, {
    x <- .parse_bca_values(input$bca_values_text)
    shiny::validate(
      shiny::need(length(x) >= 5, "Editable BCa/EIS values require at least 5 positive numbers.")
    )
    x
  }, ignoreInit = TRUE)

  bca_eis_result <- reactive({
    stored <- bca_result_store()
    if (!is.null(stored) && is.data.frame(stored)) return(stored)

    manual_mode <- FALSE

    if (isTRUE(manual_mode)) {
      x <- bca_manual_values()
    } else {
      nodes <- rv$author_nodes %||% rv$nodes20 %||% rv$nodes_full

      shiny::validate(
        shiny::need(is.data.frame(nodes) && nrow(nodes) > 0, "Please run author analysis first, or press Recompute BCa/EIS from editable values."),
        shiny::need("value" %in% names(nodes), paste0(
          "Author publication count column `value` not found. Available columns: ",
          if (is.data.frame(nodes)) paste(names(nodes), collapse = ", ") else "none"
        ))
      )

      x <- suppressWarnings(as.numeric(nodes$value))
      x <- x[is.finite(x) & x > 0]
    }

    shiny::validate(
      shiny::need(length(x) >= 5, "BCa/EIS requires at least 5 positive values; n >= 20 is recommended for stable BCa intervals.")
    )

    src_label <- bca_source_store()
    res <- shiny::withProgress(message = "Running BCa/EIS bootstrap", value = 0, {
      scan_aac_bca(
        x = x,
        R = input$bca_R,
        cutoff = input$bca_cutoff,
        conf = input$bca_conf,
        progress = function(i, total) {
          shiny::incProgress(1 / max(total, 1), detail = paste0("Point ", i, " of ", total, " | R=", input$bca_R))
        }
      )
    })
    res$EIS_point <- suppressWarnings(as.numeric(res$BCa_lower)) / input$bca_cutoff
    res$AAC_candidate <- ifelse(is.finite(res$AAC) & res$AAC > input$bca_cutoff, "AAC_OR > cutoff", "AAC_OR <= cutoff")
    res$source <- src_label
    res
  })

  bca_eis_summary <- reactive({
    summarize_bca_eis(bca_eis_result(), cutoff = input$bca_cutoff)
  })

  output$txt_bca_eis_decision <- renderText({
    sm <- bca_eis_summary()
    if (is.null(sm) || !nrow(sm)) return("No BCa/EIS result available.")
    paste0(
      "Data source = ", bca_source_store(), "\n",
      "Auto-selected elbow point = ", sm$selected_point,
      "; publication count = ", round(sm$selected_value, 3),
      "; AAC = ", round(sm$AAC, 3),
      "; BCa CI = [", round(sm$BCa_lower, 3), ", ", round(sm$BCa_upper, 3), "]",
      "; EIS = ", round(sm$EIS, 3),
      "; decision = ", sm$classification, ".\n",
      "Criterion: EIS = BCa_lower / AAC_OR cutoff. EIS > 1 indicates a statistically supported elbow; EIS <= 1 indicates no significant elbow. ",
      sm$rule
    )
  })

  output$tbl_bca_eis_summary <- DT::renderDT({
    DT::datatable(
      bca_eis_summary(),
      options = list(dom = 't', scrollX = TRUE),
      rownames = FALSE
    )
  })

  output$tbl_bca_eis <- DT::renderDT({
    res <- bca_eis_result()
    DT::datatable(
      res,
      options = list(pageLength = 10, scrollX = TRUE),
      rownames = FALSE
    )
  })

  top1_bca_eis_result <- reactive({
    stored <- top1_result_store()
    if (!is.null(stored) && is.data.frame(stored)) {
      # Recompute stored manual result if user changes Top-k, cutoff, R, or confidence.
      stored_x <- bca_values_store()
      if (!is.null(stored_x) && length(stored_x) >= 5) {
        return(scan_topk_bca_eis(
          x = stored_x,
          k = input$bca_top_k,
          R = input$bca_R,
          cutoff_or = input$bca_cutoff,
          conf = input$bca_conf
        ))
      }
      return(stored)
    }

    nodes <- rv$author_nodes %||% rv$nodes20 %||% rv$nodes_full
    shiny::validate(
      shiny::need(is.data.frame(nodes) && nrow(nodes) > 0, "Please run author analysis first, or press Recompute BCa/EIS from editable values."),
      shiny::need("value" %in% names(nodes), paste0(
        "Author publication count column `value` not found. Available columns: ",
        if (is.data.frame(nodes)) paste(names(nodes), collapse = ", ") else "none"
      ))
    )
    x <- suppressWarnings(as.numeric(nodes$value))
    x <- sort(x[is.finite(x) & x > 0], decreasing = TRUE)
    shiny::validate(shiny::need(length(x) >= 5, "Top-k BCa/EIS requires at least 5 positive values."))

    res <- shiny::withProgress(message = "Running Top-k conditional BCa/EIS", value = 0, {
      shiny::incProgress(0.25, detail = paste0("Fixing Top ", input$bca_top_k, " and resampling followers | R=", input$bca_R))
      out <- scan_topk_bca_eis(
        x = x,
        k = input$bca_top_k,
        R = input$bca_R,
        cutoff_or = input$bca_cutoff,
        conf = input$bca_conf
      )
      shiny::incProgress(0.75, detail = "Drawing Top-k dominance plot")
      out
    })
    res$source <- bca_source_store()
    res
  })

  output$txt_top1_bca_eis <- renderText({
    z <- top1_bca_eis_result()
    if (is.null(z) || !nrow(z)) return("No Top-k BCa/EIS result available.")
    .bca_caption_text(z, source_label = bca_source_store())
  })

  output$tbl_top1_bca_eis <- DT::renderDT({
    DT::datatable(
      top1_bca_eis_result(),
      options = list(dom = 't', scrollX = TRUE),
      rownames = FALSE
    )
  })

  output$plt_top1_dominance <- renderPlot({
    stored_x <- bca_values_store()
    if (!is.null(stored_x) && length(stored_x) >= 5) {
      x <- stored_x
    } else {
      nodes <- rv$author_nodes %||% rv$nodes20 %||% rv$nodes_full
      shiny::validate(
        shiny::need(is.data.frame(nodes) && nrow(nodes) > 0, "Please run author analysis first, or press Recompute BCa/EIS from editable values."),
        shiny::need("value" %in% names(nodes), "Author publication count column `value` not found.")
      )
      x <- suppressWarnings(as.numeric(nodes$value))
      x <- x[is.finite(x) & x > 0]
    }
    x <- sort(x[is.finite(x) & x > 0], decreasing = TRUE)
    shiny::validate(shiny::need(length(x) >= 5, "Top-k dominance plot requires at least 5 positive values."))
    z <- top1_bca_eis_result()
    .plot_topk_dominance(x, z, main_title = "Top-k Dominance vs Rest: Conditional BCa/EIS")
  })

  output$plt_bca_eis <- renderPlot({
    df <- bca_eis_result()
    shiny::validate(shiny::need(is.data.frame(df) && nrow(df) > 0, "No BCa/EIS result available."))

    df$point_index <- suppressWarnings(as.numeric(df$point_index))
    df$AAC <- suppressWarnings(as.numeric(df$AAC))
    df$BCa_lower <- suppressWarnings(as.numeric(df$BCa_lower))
    df$BCa_upper <- suppressWarnings(as.numeric(df$BCa_upper))
    if (!("EIS_point" %in% names(df))) df$EIS_point <- df$BCa_lower / input$bca_cutoff
    df$EIS_point <- suppressWarnings(as.numeric(df$EIS_point))
    smry <- summarize_bca_eis(df, cutoff = input$bca_cutoff)
    sel <- which(df$point_index == smry$selected_point)
    if (!length(sel)) sel <- integer(0)

    yy <- c(df$AAC, df$BCa_lower, df$BCa_upper, input$bca_cutoff)
    yy <- yy[is.finite(yy)]
    if (!length(yy)) yy <- c(0, 1)
    y_lim <- range(yy, na.rm = TRUE)
    pad <- diff(y_lim) * 0.10
    if (!is.finite(pad) || pad == 0) pad <- 0.5
    y_lim <- y_lim + c(-pad, pad)

    plot(
      df$point_index,
      df$AAC,
      type = "b",
      pch = 19,
      lwd = 2,
      xlab = "Rank position / internal point index",
      ylab = "AAC OR ratio",
      main = "BCa/EIS: AAC_OR Curve with BCa Confidence Intervals and Smooth Trend",
      ylim = y_lim
    )

    ok_ci <- is.finite(df$point_index) & is.finite(df$BCa_lower) & is.finite(df$BCa_upper)
    if (any(ok_ci)) {
      arrows(
        x0 = df$point_index[ok_ci],
        y0 = df$BCa_lower[ok_ci],
        x1 = df$point_index[ok_ci],
        y1 = df$BCa_upper[ok_ci],
        angle = 90,
        code = 3,
        length = 0.05
      )
    }

    ok_smooth <- is.finite(df$point_index) & is.finite(df$AAC)
    if (sum(ok_smooth) >= 4) {
      sm <- stats::lowess(df$point_index[ok_smooth], df$AAC[ok_smooth], f = 2/3, iter = 1)
      lines(sm$x, sm$y, lwd = 3, lty = 1)
    }

    abline(h = input$bca_cutoff, lty = 2, lwd = 2)

    # Criterion note printed directly on the plot
    usr <- par("usr")
    text(
      x = usr[1] + 0.02 * (usr[2] - usr[1]),
      y = usr[4] - 0.06 * (usr[4] - usr[3]),
      labels = "Criterion: Significant if BCa lower > AAC_OR cutoff; EIS = BCa lower / cutoff > 1",
      adj = c(0, 1), cex = 0.8
    )

    sig <- !is.na(df$decision) & df$decision == "Significant elbow" & is.finite(df$point_index) & is.finite(df$AAC)
    if (any(sig)) {
      points(df$point_index[sig], df$AAC[sig], pch = 8, cex = 1.6, lwd = 2)
    }

    # Automatic selected elbow/candidate: always shown in red.
    if (length(sel) && is.finite(df$point_index[sel[1]]) && is.finite(df$AAC[sel[1]])) {
      points(df$point_index[sel[1]], df$AAC[sel[1]], pch = 19, cex = 2.0, col = "red")
      text(
        df$point_index[sel[1]], df$AAC[sel[1]],
        labels = paste0("Selected: ", smry$classification, "\nEIS=", round(smry$EIS, 2), " (Sig if >1)"),
        pos = 4, cex = 0.95, col = "red"
      )
    }

    legend(
      "topright",
      legend = c("AAC_OR", "BCa CI", "LOWESS smooth", paste0("AAC_OR cutoff = ", input$bca_cutoff), "Significant elbow", "Auto-selected elbow / candidate"),
      lty = c(1, 1, 1, 2, NA, NA),
      pch = c(19, NA, NA, NA, 8, 19),
      lwd = c(2, 1, 3, 2, NA, NA),
      col = c("black", "black", "black", "black", "black", "red"),
      bty = "n"
    )
  })

  output$dl_bca_eis <- downloadHandler(
    filename = function() paste0("BCa_EIS_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content = function(file) {
      res <- tryCatch(bca_eis_result(), error = function(e) data.frame(Message = conditionMessage(e)))
      utils::write.csv(res, file, row.names = FALSE, fileEncoding = "UTF-8")
    }
  )


  output$dl_bca_topk_png <- downloadHandler(
    filename = function() paste0("TopK_BCa_EIS_dominance_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
    content = function(file) {
      stored_x <- bca_values_store()
      if (!is.null(stored_x) && length(stored_x) >= 5) {
        x <- stored_x
      } else {
        nodes <- rv$author_nodes %||% rv$nodes20 %||% rv$nodes_full
        x <- if (is.data.frame(nodes) && "value" %in% names(nodes)) suppressWarnings(as.numeric(nodes$value)) else numeric(0)
      }
      x <- sort(x[is.finite(x) & x > 0], decreasing = TRUE)
      z <- top1_bca_eis_result()
      grDevices::png(file, width = 1800, height = 1200, res = 180)
      on.exit(grDevices::dev.off(), add = TRUE)
      .plot_topk_dominance(x, z, main_title = "Top-k Dominance vs Rest: Conditional BCa/EIS")
    }
  )

  output$dl_bca_topk_figure_zip <- downloadHandler(
    filename = function() paste0("TopK_BCa_EIS_figure_caption_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".zip"),
    content = function(file) {
      td <- tempfile("bca_topk_")
      dir.create(td, recursive = TRUE, showWarnings = FALSE)
      png_file <- file.path(td, "TopK_BCa_EIS_dominance.png")
      cap_file <- file.path(td, "TopK_BCa_EIS_caption.txt")
      stored_x <- bca_values_store()
      if (!is.null(stored_x) && length(stored_x) >= 5) {
        x <- stored_x
      } else {
        nodes <- rv$author_nodes %||% rv$nodes20 %||% rv$nodes_full
        x <- if (is.data.frame(nodes) && "value" %in% names(nodes)) suppressWarnings(as.numeric(nodes$value)) else numeric(0)
      }
      x <- sort(x[is.finite(x) & x > 0], decreasing = TRUE)
      z <- top1_bca_eis_result()
      grDevices::png(png_file, width = 1800, height = 1200, res = 180)
      .plot_topk_dominance(x, z, main_title = "Top-k Dominance vs Rest: Conditional BCa/EIS")
      grDevices::dev.off()
      writeLines(.bca_caption_text(z, source_label = bca_source_store()), cap_file, useBytes = TRUE)
      old <- getwd(); on.exit(setwd(old), add = TRUE)
      setwd(td)
      utils::zip(zipfile = file, files = c("TopK_BCa_EIS_dominance.png", "TopK_BCa_EIS_caption.txt"))
    }
  )

  lotka_dist_reactive <- reactive({
    if (!is.null(rv$author_lists_fl) && length(rv$author_lists_fl) > 0) {
      return(.build_lotka_input_from_author_lists(rv$author_lists_fl))
    }
    if (!is.null(rv$author_lists_all) && length(rv$author_lists_all) > 0) {
      return(.build_lotka_input_from_author_lists(rv$author_lists_all))
    }
    if (!is.null(rv$nodes_full) && is.data.frame(rv$nodes_full) && nrow(rv$nodes_full) > 0 && ("n_pubmed" %in% names(rv$nodes_full))) {
      counts <- suppressWarnings(as.numeric(rv$nodes_full$n_pubmed))
      counts <- counts[is.finite(counts) & counts > 0]
      if (!length(counts)) return(NULL)
      tb <- as.data.frame(table(counts), stringsAsFactors = FALSE)
      names(tb) <- c("papers", "authors")
      tb$papers <- as.integer(as.character(tb$papers))
      tb$authors <- as.integer(tb$authors)
      tb <- tb[order(tb$papers), , drop = FALSE]
      rownames(tb) <- NULL
      return(tb)
    }
    NULL
  })

  lotka_res_reactive <- reactive({
    df <- lotka_dist_reactive()
    if (is.null(df) || !is.data.frame(df) || nrow(df) < 3) return(NULL)
    tryCatch(test_lotka(df, make_plot = FALSE, quiet = TRUE), error = function(e) NULL)
  })

  # Always compute IP status on load
  observe({
    g <- tryCatch(
      ipm_gate_session(session,
                       cmc = input$cmc %||% "",
                       app_dir = APP_DIR,
                       inc_count_on_allow = FALSE),
      error = function(e) NULL
    )

    if (!is.null(g)) {
      rv$ip_access_type <- g$policy %||% "unknown"
      rv$ip_addr <- g$ip %||% "unknown"
    }
  })
  print("Rendering access status")
  output$ip_pass_status <- renderText({
    paste0("Access mode: ",
           rv$ip_access_type,
           " | IP: ",
           rv$ip_addr)
  })

# ---- CMC + IP trial gate (ip.txt + iplist.txt) ------------------------------------------
# Policy (as requested):
# - Demo (button): always allowed (no CMC; does not consume trial).
# - Normal Fetch / domain runs / URL autorun:
#     * If iplist.txt exists and contains IPs: only those IPs are always allowed.
#     * Otherwise: each new IP gets ONE free run (trial) when not present in ip.txt.
#       After first run, IP is recorded in ip.txt and future runs require a valid 10-digit CMC.
# - URL autorun: cmc=test is allowed only for first-time IP (not in ip.txt). Otherwise 10-digit CMC required.

`%||%` <- function(a,b){ if (!is.null(a)) a else b }

.is_cmc_10 <- function(x){
  x <- as.character(x %||% "")
  x <- trimws(x)
  isTRUE(nzchar(x) && grepl("^[0-9]{10}$", x))
}

.app_dir <- local({
  d <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NA_character_)
  if (is.null(d) || !nzchar(d) || is.na(d)) d <- getwd()
  tryCatch(normalizePath(d, winslash = "/", mustWork = FALSE), error = function(e) d)
})
.ip_file <- function(){ file.path(.app_dir, "ip.txt") }
.iplist_file <- function(){ file.path(.app_dir, "iplist.txt") }

.read_ip_txt <- function(){
  f <- .ip_file()
  if (!file.exists(f)) {
    return(data.frame(ip=character(), recent_date=character(), total_count=integer(), cmc=character(),
                      stringsAsFactors = FALSE))
  }
  df <- tryCatch(utils::read.delim(f, sep="\t", header=FALSE, stringsAsFactors=FALSE), error=function(e) NULL)
  if (is.null(df) || ncol(df) < 4) {
    return(data.frame(ip=character(), recent_date=character(), total_count=integer(), cmc=character(),
                      stringsAsFactors = FALSE))
  }
  df <- df[,1:4, drop=FALSE]
  names(df) <- c("ip","recent_date","total_count","cmc")
  df$ip <- as.character(df$ip)
  df$recent_date <- as.character(df$recent_date)
  df$total_count <- suppressWarnings(as.integer(df$total_count))
  df$total_count[!is.finite(df$total_count)] <- 0L
  df$cmc <- as.character(df$cmc)
  df
}

.write_ip_txt <- function(df){
  f <- .ip_file()
  if (!is.data.frame(df) || !nrow(df)) return(invisible(FALSE))
  df <- df[, c("ip","recent_date","total_count","cmc"), drop=FALSE]
  tryCatch(utils::write.table(df, f, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE),
           error=function(e) NULL)
  invisible(TRUE)
}

.ip_seen <- function(ip){
  df <- .read_ip_txt()
  if (!nrow(df)) return(FALSE)
  ip <- as.character(ip %||% "unknown"); if (!nzchar(ip)) ip <- "unknown"
  any(df$ip == ip)
}

.upsert_ip <- function(ip, cmc, inc_count = TRUE){
  ip <- as.character(ip %||% "unknown"); if (!nzchar(ip)) ip <- "unknown"
  cmc <- as.character(cmc %||% "")
  df <- .read_ip_txt()
  today <- as.character(Sys.Date())
  if (!nrow(df) || !(ip %in% df$ip)) {
    df <- rbind(df, data.frame(ip=ip, recent_date=today, total_count=0L, cmc=cmc, stringsAsFactors = FALSE))
  } else {
    i <- match(ip, df$ip)
    df$recent_date[i] <- today
    df$cmc[i] <- cmc
  }
  if (isTRUE(inc_count)) {
    i <- match(ip, df$ip)
    df$total_count[i] <- as.integer(df$total_count[i] %||% 0L) + 1L
  }
  .write_ip_txt(df)
  df
}

.read_iplist <- function(){
  f <- .iplist_file()
  if (!file.exists(f)) return(character())
  x <- tryCatch(readLines(f, warn = FALSE, encoding = "UTF-8"), error = function(e) character())
  x <- trimws(x)
  x <- x[nzchar(x)]
  unique(x)
}

.ip_allowlisted <- function(ip){
  allow <- .read_iplist()
  if (!length(allow)) return(TRUE)
  ip <- as.character(ip %||% "unknown"); if (!nzchar(ip)) ip <- "unknown"
  ip %in% allow
}

get_client_ip <- function(){
  hdr <- tryCatch(session$request$HTTP_X_FORWARDED_FOR, error = function(e) NULL)
  if (!is.null(hdr) && nzchar(hdr)) return(trimws(strsplit(hdr, ",")[[1]][1]))
  ip <- tryCatch(session$request$REMOTE_ADDR, error = function(e) NULL)
  if (!is.null(ip) && nzchar(ip)) return(trimws(ip))
  "unknown"
}

.gate_normal_run <- function(ip, cmc){
  if (!.ip_allowlisted(ip)) {
    showModal(modalDialog(title="Access blocked", "Your IP is not in iplist.txt.", easyClose=TRUE, footer=modalButton("OK")))
    return(FALSE)
  }
  if (.is_cmc_10(cmc)) return(TRUE)
  if (!.ip_seen(ip)) return(TRUE)
  showModal(modalDialog(title="CMC required",
                        "This IP has already used the one-time trial. Please enter a 10-digit numeric CMC.",
                        easyClose=TRUE, footer=modalButton("OK")))
  FALSE
}

.gate_url_autorun <- function(ip, cmc_param){
  if (!.ip_allowlisted(ip)) {
    showModal(modalDialog(title="Access blocked", "Your IP is not in iplist.txt.", easyClose=TRUE, footer=modalButton("OK")))
    return(FALSE)
  }
  cmc_param <- tolower(trimws(as.character(cmc_param %||% "")))
  if (.is_cmc_10(cmc_param)) return(TRUE)
  if (cmc_param == "test") {
    if (!.ip_seen(ip)) return(TRUE)
    showModal(modalDialog(title="Trial already used",
                          "cmc=test is allowed only once (first-time IP) via URL autorun. Please provide a 10-digit CMC.",
                          easyClose=TRUE, footer=modalButton("OK")))
    return(FALSE)
  }
  if (!.ip_seen(ip)) return(TRUE)
  showModal(modalDialog(title="CMC required",
                        "URL autorun requires cmc=test (first-time IP only) or a 10-digit numeric CMC.",
                        easyClose=TRUE, footer=modalButton("OK")))
  FALSE
}


  # ---- Reactive state (must be defined before any output uses it) ----
  rv <- reactiveValues(
      pubmeta = NULL,
      log = "",
      nodes_full = NULL,
      edges_full = NULL,
      # --- report stats (per-domain) ---
      author_max_items = NA_real_, author_median_items = NA_real_, author_n_ge2 = NA_integer_,
      country_max_items = NA_real_, country_median_items = NA_real_, country_n_ge2 = NA_integer_,
      stateprov_max_items = NA_real_, stateprov_median_items = NA_real_, stateprov_n_ge2 = NA_integer_,
      inst_max_items = NA_real_, inst_median_items = NA_real_, inst_n_ge2 = NA_integer_,
      dept_max_items = NA_real_, dept_median_items = NA_real_, dept_n_ge2 = NA_integer_,
      mesh_max_items = NA_real_, mesh_median_items = NA_real_, mesh_n_ge2 = NA_integer_,
      jy_max_items = NA_real_, jy_median_items = NA_real_, jy_n_ge2 = NA_integer_,
  
      AAC_value_author = NA_real_, AAC_value2_author = NA_real_, AAC_ss_author = NA_real_, AAC_a_star_author = NA_real_,
      AAC_value_country = NA_real_, AAC_value2_country = NA_real_, AAC_ss_country = NA_real_, AAC_a_star_country = NA_real_,
      AAC_value_stateprov = NA_real_, AAC_value2_stateprov = NA_real_, AAC_ss_stateprov = NA_real_, AAC_a_star_stateprov = NA_real_,
      AAC_value_inst = NA_real_, AAC_value2_inst = NA_real_, AAC_ss_inst = NA_real_, AAC_a_star_inst = NA_real_,
      AAC_value_dept = NA_real_, AAC_value2_dept = NA_real_, AAC_ss_dept = NA_real_, AAC_a_star_dept = NA_real_,
      AAC_value_mesh = NA_real_, AAC_value2_mesh = NA_real_, AAC_ss_mesh = NA_real_, AAC_a_star_mesh = NA_real_,
      AAC_value_jy = NA_real_, AAC_value2_jy = NA_real_, AAC_ss_jy = NA_real_, AAC_a_star_jy = NA_real_,
  
  
      # --- computed results (Top20) ---
      author_nodes = NULL, author_edges = NULL,
      country_nodes = NULL, country_edges = NULL,
      stateprov_nodes = NULL, stateprov_edges = NULL,
      inst_nodes = NULL, inst_edges = NULL,
      dept_nodes = NULL, dept_edges = NULL,
      mesh_nodes = NULL, mesh_edges = NULL,
  
      # --- cached parsed term lists + raw co-occ edges (for on-demand domain runs) ---
      countries_list = NULL, stateprov_list = NULL,
      inst_list = NULL, dept_list = NULL, mesh_list = NULL,
      edges_country = NULL, edges_stateprov = NULL,
      edges_inst = NULL, edges_dept = NULL, edges_mesh = NULL,
  
      report_path = NULL,
      pmids = NULL,
            icite_url = NULL,
      perf_long_auto = NULL,
      ip_access_type = NA_character_,
      ip_addr = "unknown",
      ip_total_count = NA_integer_,
      ip_recent_date = NA_character_,
      ip_cmc = "",
      taaa_df = NULL,
      taaa_freq_df = NULL,
      taaa_profile_map_df = NULL,
      taaa_kappa_df = NULL,
      taaa_conf_df = NULL
    )

  demo_mode <- reactiveVal(FALSE)  # internal: bypass CMC gate for free Demo


# ---- IP tab outputs (access status + allowlist + log) ----
output$ip_access_mode_ui <- renderUI({
  ip <- rv$ip_addr %||% "unknown"
  allow <- tryCatch(ipm_read_iplist(app_dir = APP_DIR), error = function(e) character())
  ipdf <- tryCatch(ipm_read_ip_txt(app_dir = APP_DIR), error = function(e) NULL)

  in_allow <- length(allow) && isTRUE(ip %in% allow)
  in_log <- !is.null(ipdf) && nrow(ipdf) && isTRUE(ip %in% ipdf$ip)

  # Derive a human-readable access mode
  mode <- if (in_allow) {
    "IP allowlisted (always allowed)"
  } else if (identical(rv$ip_access_type %||% "", "CMC pass")) {
    "CMC validated (allowed)"
  } else if (in_log) {
    "CMC required (trial already used)"
  } else {
    "Trial available (first upload free)"
  }

  # Color cue (green/orange/red)
  col <- if (grepl("always allowed|validated", mode, ignore.case = TRUE)) "#1b7f3a" else if (grepl("Trial", mode, ignore.case = TRUE)) "#b36b00" else "#b00020"

  tags$div(style="margin:6px 0 10px 0;",
           tags$b("Access mode: "),
           tags$span(style=paste0("color:", col, "; font-weight:700;"), mode),
           tags$div(style="color:#666; font-size:12px; margin-top:2px;",
                    "Rule of thumb: allowlist > valid CMC > trial(first time) > blocked/CMC required."))
})

output$ip_status <- renderTable({
  data.frame(
    Access = rv$ip_access_type %||% NA_character_,
    AccessMode = {
      ip <- rv$ip_addr %||% "unknown"
      allow <- tryCatch(ipm_read_iplist(app_dir = APP_DIR), error = function(e) character())
      ipdf <- tryCatch(ipm_read_ip_txt(app_dir = APP_DIR), error = function(e) NULL)
      in_allow <- length(allow) && isTRUE(ip %in% allow)
      in_log <- !is.null(ipdf) && nrow(ipdf) && isTRUE(ip %in% ipdf$ip)
      if (in_allow) "IP allowlisted" else if (identical(rv$ip_access_type %||% "", "CMC pass")) "CMC validated" else if (in_log) "CMC required" else "Trial available"
    },
    IP = rv$ip_addr %||% "unknown",
    TotalCount = rv$ip_total_count %||% NA_integer_,
    RecentDate = rv$ip_recent_date %||% NA_character_,
    CMC = rv$ip_cmc %||% "",
    check.names = FALSE
  )
})

output$ip_allowlist <- DT::renderDT({
  allow <- tryCatch(ipm_read_iplist(app_dir = APP_DIR), error = function(e) character())
  ip <- rv$ip_addr %||% "unknown"
  df <- data.frame(
    allowlisted_ip = if (length(allow)) allow else character(),
    is_current_ip = if (length(allow)) (allow == ip) else logical(),
    stringsAsFactors = FALSE
  )
  DT::datatable(df, rownames = FALSE, options = list(pageLength = 25, dom = "tip"))
})

output$ip_log <- DT::renderDT({
  df <- tryCatch(ipm_read_ip_txt(app_dir = APP_DIR), error = function(e) NULL)
  if (is.null(df) || !nrow(df)) {
    df <- data.frame(ip=character(), recent_date=character(), total_count=integer(), cmc=character(), stringsAsFactors = FALSE)
  }
  DT::datatable(df, rownames = FALSE, options = list(pageLength = 25, dom = "tip"))
})


# ==========================================================
# Performance report PNG (8 domains, Top5 + AAC) - CSV upload
# CSV columns: Domain, Element, Value (or value)
# ==========================================================
perf_data_long <- reactive({
  # Priority: uploaded CSV (manual) > auto-built from fetched PubMed
  if (!is.null(input$perf_file)) {
    df <- tryCatch({
      read.csv(input$perf_file$datapath, stringsAsFactors = FALSE, check.names = FALSE)
    }, error = function(e) NULL)
    validate(need(!is.null(df) && nrow(df) > 0, "Uploaded CSV is empty or unreadable."))

    nms <- names(df)
    nms_l <- tolower(nms)

    pick_col <- function(cands){
      cands_l <- tolower(cands)
      hit <- which(nms_l %in% cands_l)
      if (length(hit)) nms[hit[1]] else NA_character_
    }

    col_domain  <- pick_col(c("Domain","domain"))
    col_element <- pick_col(c("Element","element","Item","item","Name","name"))
    col_value   <- pick_col(c("Value","value","Count","count","N","n"))

    validate(need(!is.na(col_domain) && !is.na(col_element) && !is.na(col_value),
                  "CSV must have columns Domain, Element, and Value (or value)."))

    out <- data.frame(
      Domain  = as.character(df[[col_domain]]),
      Element = as.character(df[[col_element]]),
      Value   = suppressWarnings(as.numeric(df[[col_value]])),
      stringsAsFactors = FALSE
    )
    out <- out[!is.na(out$Domain) & nzchar(out$Domain) & !is.na(out$Element) & nzchar(out$Element), , drop=FALSE]
    return(out)
  }

  # auto
  if (!is.null(rv$perf_long_auto) && is.data.frame(rv$perf_long_auto) && nrow(rv$perf_long_auto) > 0) {
    return(rv$perf_long_auto)
  }

  validate(need(FALSE, "No performance data yet. Click 'Fetch PubMed' first (auto), or upload a long-format CSV."))
})
# AAC formula: r/(1+r)
.perf_compute_aac <- function(values) {
  values <- values[!is.na(values)]
  if (!length(values)) return(NA_real_)
  n_vals <- length(values)

  if (n_vals >= 3) {
    if (!isTRUE(values[2] != 0) || !isTRUE(values[3] != 0)) return(NA_real_)
    r <- (values[1] / values[2]) / (values[2] / values[3])
    return(r / (1 + r))
  } else if (n_vals == 2) {
    if (!isTRUE(values[2] != 0)) return(NA_real_)
    r <- values[1] / values[2]
    return(r / (1 + r))
  } else if (n_vals == 1) {
    return(0.5)
  }
  NA_real_
}

.perf_build_panel <- function(df, domain_var) {
  dom_all <- df[df$Domain == domain_var, , drop=FALSE]

  if (!nrow(dom_all)) {
    plot.new()
    title(main = paste(domain_var, "(No Data)"), font.main = 2)
    return(invisible(NULL))
  }

  # Allow a dedicated AAC row (Domain=..., Element="AAC", Value=...)
  aac_row <- dom_all[tolower(dom_all$Element) %in% c("aac","aac_value","aac(value)"), , drop=FALSE]
  aac_from_data <- if (nrow(aac_row)) suppressWarnings(as.numeric(aac_row$Value[1])) else NA_real_

  dom <- dom_all[!(tolower(dom_all$Element) %in% c("aac","aac_value","aac(value)")), , drop=FALSE]
  if (!nrow(dom)) dom <- dom_all

  # Year: show last 5 in chronological order; others: Top 5 by Value
  if (identical(domain_var, "Year")) {
    dom$._year_num <- suppressWarnings(as.numeric(as.character(dom$Element)))
    dom <- dom[order(dom$._year_num), , drop=FALSE]
    if (nrow(dom) > 5) dom <- tail(dom, 5)
  } else {
    dom <- dom[order(-dom$Value), , drop=FALSE]
    if (nrow(dom) > 5) dom <- head(dom, 5)
  }

  values <- dom$Value
  aac <- if (!is.na(aac_from_data)) aac_from_data else .perf_compute_aac(values)

  plot.new()
  title(main = domain_var, font.main = 2, cex.main = 1.9, line = 0.5, col.main = "#d62728")

  n_rows <- nrow(dom)
  text_cex <- if (n_rows <= 4) 2.4 else if (n_rows <= 7) 2.0 else 1.6
  y_positions <- seq(0.85, 0.2, length.out = n_rows)

  for (i in seq_len(n_rows)) {
    text(0.05, y_positions[i], dom$Element[i], adj = 0, cex = text_cex, font = 2)
    text(0.95, y_positions[i], dom$Value[i],   adj = 1, cex = text_cex, font = 2)
  }

  if (!is.na(aac)) {
    text(0.5, 0.08, paste("AAC =", round(aac, 4)), cex = 1.8, font = 2, col = "#d62728")
  } else {
    text(0.5, 0.08, "AAC = NA", cex = 1.8, font = 2, col = "#d62728")
  }

  invisible(NULL)
}

.perf_render_plot <- function(df) {
  domains <- c("Country","Journal","Institute","Year","Department","State/Province","Author","MeSH")
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)

  par(mfrow = c(4, 2))
  par(mar = c(1.1, 0.8, 2.2, 0.8))
  par(oma = c(0, 0, 0, 0))

  for (d in domains) .perf_build_panel(df, d)
}


# Build performance long-format from fetched PubMed results (Top5 + AAC, 8 domains)
.build_perf_long_auto <- function() {
  doms <- c("Country","Journal","Institute","Year","Department","State/Province","Author","MeSH")
  out <- data.frame(Domain=character(), Element=character(), Value=numeric(), stringsAsFactors = FALSE)

  for (d in doms) {
    x <- .get_domain_summary(d)
    df <- x$df
    aac_val <- as_scalar_or_na(x$AAC, NA)

    if (identical(d, "State/Province")) {
      # summary uses "State/Province" domain name already
    }

    if (!is.null(df) && is.data.frame(df) && nrow(df) > 0) {
      # Year: keep chronological last 5 rows; others: top 5 by value
      if (identical(d, "Year")) {
        df <- df[order(suppressWarnings(as.numeric(as.character(df$name)))), , drop=FALSE]
        if (nrow(df) > 5) df <- tail(df, 5)
      } else {
        if ("value" %in% names(df)) df <- df[order(df$value, decreasing=TRUE, na.last=TRUE), , drop=FALSE]
        df <- head(df, 5)
      }
      if (!("name" %in% names(df))) df$name <- as.character(df[[1]])
      if (!("value" %in% names(df))) df$value <- suppressWarnings(as.numeric(as.character(df[[2]])))

      out <- rbind(out, data.frame(Domain=d, Element=as.character(df$name), Value=as.numeric(df$value), stringsAsFactors=FALSE))
    }

    # add AAC row for the panel
    out <- rbind(out, data.frame(Domain=d, Element="AAC", Value=as.numeric(aac_val), stringsAsFactors=FALSE))
  }

  out
}

observeEvent(input$perf_use_auto, {
  if (is.null(rv$pmids) || !length(rv$pmids)) {
    showNotification("No fetched PubMed data yet. Click 'Fetch PubMed' first.", type="error", duration=5)
    return(invisible(NULL))
  }

  tryCatch({
    rv$perf_long_auto <- .build_perf_long_auto()

    if (is.null(rv$perf_long_auto) || !is.data.frame(rv$perf_long_auto) || nrow(rv$perf_long_auto) == 0) {
      showNotification("Fetched PubMed exists, but no performance long-format data was generated yet. Try running domains first.", type="error", duration=6)
      return(invisible(NULL))
    }

    showNotification("Performance report updated from fetched PubMed.", type="message", duration=3)
  }, error = function(e){
    # Never crash the app: show error as notification
    showNotification(paste0("Performance auto-build failed: ", conditionMessage(e)), type="error", duration=8)
    rv$perf_long_auto <- NULL
  })

}, ignoreInit = TRUE)


output$perf_plot_ui <- renderUI({
  # Keep homepage/tabs unchanged; only make the Performance report section compact when no data yet.
  # When there is no performance data, use a small plot area so the "Summary report" section stays visible.
  h <- "900px"
  df_try <- tryCatch(perf_data_long(), error = function(e) NULL)
  if (is.null(df_try) || !is.data.frame(df_try) || nrow(df_try) == 0) h <- "260px"
  plotOutput("perf_preview_plot", height = h)
})

output$perf_preview_plot <- renderPlot({
  req(perf_data_long())
  .perf_render_plot(perf_data_long())
})

output$perf_download_png <- downloadHandler(
  filename = function() "performance_report.png",
  contentType = "image/png",
  content  = function(file) {
    req(perf_data_long())
    grDevices::png(file, width = 3200, height = 2400, res = 300)
    on.exit(grDevices::dev.off(), add = TRUE)
    .perf_render_plot(perf_data_long())
  }
)



# ---- ReadMe (concise) ----
output$readme_html <- shiny::renderUI(shiny::HTML(paste0(
  "<div style='line-height:1.55'>",
  "<p><b>What this app does</b>: Fetch PubMed (by query) or parse an uploaded MEDLINE .txt, then build Top-20 networks (Author / Journal-Year / Country / Institute / Department / MeSH) and show plots (Sankey, Kano, SSplot) + export files.</p>",
  "<ol>",
  "<li><b>Choose data source</b>:<ul>",
  "<li><b>PubMed query</b>: enter a term and click <b>Fetch PubMed</b>.</li>",
  "<li><b>Uploaded MEDLINE</b>: upload .txt, keep <i>Use uploaded MEDLINE (ignore query)</i> checked, then click <b>Fetch PubMed</b>.</li>",
  "<li>If the page URL contains <code>?term=...</code>, the app will automatically <b>uncheck</b> the upload option to avoid mixing sources.</li>",
  "</ul></li>",
  "<li><b>Select a domain</b> (tabs or dropdowns) to run FLCA + MajorSampling Top-20 nodes/edges.</li>",
  "<li><b>Outputs</b>: tables (Sampled nodes/edges), plots, and SankeyMATIC text/link when available.</li>",
  "</ol>",
  "<p class='small-note'>Tip: If a tab shows 'no data yet', run Fetch PubMed first, then open the domain tab again.</p>",
  "</div>"
)))

  # =========================
  # OUTPUT BINDINGS (restore)
  # - Status (log)
  # - Sankey (Author) + SankeyMATIC
  # - World map (Country pre-FLCA distribution)
  # - USA map (State/Province pre-FLCA, USA only)
  # =========================

  # Status panel
  output$log <- shiny::renderText({
    if (!exists("rv")) return("")
    rv$log %||% ""
  })

  # Helper: safe `%||%`
  `%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# -------- Author sampled nodes / edges tables --------
output$tbl_author_nodes <- DT::renderDT({
  df <- rv$author_nodes %||% rv$nodes20
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) {
    return(DT::datatable(data.frame(Message = "No sampled nodes yet. Click 'Fetch PubMed' first."), options=list(dom='t')))
  }
  dt <- DT::datatable(df, options=list(pageLength=10, scrollX=TRUE), rownames=FALSE)
  num_cols <- which(sapply(df, is.numeric))
  if (length(num_cols) > 0) dt <- DT::formatRound(dt, columns=num_cols, digits=2)
  dt
})

output$tbl_author_edges <- DT::renderDT({
  df <- rv$author_edges %||% rv$edges20 %||% rv$data20
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) {
    return(DT::datatable(data.frame(Message = "No sampled edges yet. Click 'Fetch PubMed' first."), options=list(dom='t')))
  }
  dt <- DT::datatable(df, options=list(pageLength=10, scrollX=TRUE), rownames=FALSE)
  num_cols <- which(sapply(df, is.numeric))
  if (length(num_cols) > 0) dt <- DT::formatRound(dt, columns=num_cols, digits=2)
  dt
})


  
# -------- Metadata domain tables (Top20 sampled nodes/edges) --------
.dt_msg <- function(msg){
  DT::datatable(data.frame(Message = msg), options = list(dom = 't'), rownames = FALSE)
}

.dt_fmt <- function(df){
  dt <- DT::datatable(df, options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
  num_cols <- which(sapply(df, is.numeric))
  if (length(num_cols) > 0) dt <- DT::formatRound(dt, columns = num_cols, digits = 2)
  dt
}

# Country
output$tbl_country_nodes <- DT::renderDT({
  df <- rv$country_nodes
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) return(.dt_msg("No Country nodes yet. Click 'Run Country' after fetching/uploading."))
  .dt_fmt(df)
})
output$tbl_country_edges <- DT::renderDT({
  df <- rv$country_edges
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) return(.dt_msg("No Country edges yet. Click 'Run Country' after fetching/uploading."))
  .dt_fmt(df)
})

# State/Province
output$tbl_stateprov_nodes <- DT::renderDT({
  df <- rv$stateprov_nodes
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) return(.dt_msg("No State/Province nodes yet. Click 'Run State/Prov' after fetching/uploading."))
  .dt_fmt(df)
})
output$tbl_stateprov_edges <- DT::renderDT({
  df <- rv$stateprov_edges
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) return(.dt_msg("No State/Province edges yet. Click 'Run State/Prov' after fetching/uploading."))
  .dt_fmt(df)
})

# Institute
output$tbl_inst_nodes <- DT::renderDT({
  df <- rv$inst_nodes
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) return(.dt_msg("No Institute nodes yet. Click 'Run Institute' after fetching/uploading."))
  .dt_fmt(df)
})
output$tbl_inst_edges <- DT::renderDT({
  df <- rv$inst_edges
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) return(.dt_msg("No Institute edges yet. Click 'Run Institute' after fetching/uploading."))
  .dt_fmt(df)
})

# Department
output$tbl_dept_nodes <- DT::renderDT({
  df <- rv$dept_nodes
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) return(.dt_msg("No Department nodes yet. Click 'Run Department' after fetching/uploading."))
  .dt_fmt(df)
})
output$tbl_dept_edges <- DT::renderDT({
  df <- rv$dept_edges
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) return(.dt_msg("No Department edges yet. Click 'Run Department' after fetching/uploading."))
  .dt_fmt(df)
})

# MeSH
output$tbl_mesh_nodes <- DT::renderDT({
  df <- rv$mesh_nodes
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) return(.dt_msg("No MeSH nodes yet. Click 'Run MeSH' after fetching/uploading."))
  .dt_fmt(df)
})
output$tbl_mesh_edges <- DT::renderDT({
  df <- rv$mesh_edges
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) return(.dt_msg("No MeSH edges yet. Click 'Run MeSH' after fetching/uploading."))
  .dt_fmt(df)
})

# ---- Slope-tab extra plots ----
# The old separate analysis-tab block has been removed intentionally.
# These three tabs are rendered later by the FINAL OVERRIDE block,
# which reuses the stable Slope tab pipeline: .sg_final_build_wide(),
# .sg_final_long_from_wide(), and .sg_final_draw_domain().
# This prevents the old DT/plotly/ggplot text-label outputs from reappearing.


# Journal/Year table
output$tbl_journal_year <- DT::renderDT({
  req(rv$pubmeta)
  df <- rv$pubmeta
  if (!all(c("Year","Journal") %in% names(df))) return(.dt_msg("No Journal/Year data available."))
  tab <- df %>%
    dplyr::mutate(Year = suppressWarnings(as.integer(as.character(.data$Year)))) %>%
    dplyr::filter(!is.na(.data$Year), nzchar(.data$Journal)) %>%
    dplyr::count(Year, Journal, name = "n") %>%
    dplyr::arrange(dplyr::desc(n), dplyr::desc(Year)) %>%
    dplyr::slice_head(n = 200)
  .dt_fmt(tab)
})

# Year/Articles table
output$tbl_year_articles <- DT::renderDT({
  req(rv$pubmeta)
  df <- rv$pubmeta
  if (!("Year" %in% names(df))) return(.dt_msg("No Year data available."))
  tab <- df %>%
    dplyr::mutate(Year = suppressWarnings(as.integer(as.character(.data$Year)))) %>%
    dplyr::filter(!is.na(.data$Year)) %>%
    dplyr::count(Year, name = "n_articles") %>%
    dplyr::arrange(.data$Year)
  .dt_fmt(tab)
})


# -------- Sankey (Domain) --------
  output$ui_sankey_link <- shiny::renderUI({
    tags$a("Open SankeyMATIC builder", href = "https://sankeymatic.com/build/", target = "_blank")
  })

  output$ui_sankey_author <- shiny::renderUI({
    dom <- input$sankey_domain %||% "Author"
    dd  <- .get_domain(dom)
    if (is.null(dd$edges) || !is.data.frame(dd$edges) || nrow(dd$edges) == 0) {
      return(tags$div(class="small-note", sprintf("Sankey: no edges yet for domain '%s'. Run PubMed / upload, then click a domain button to build Top-20.", dom)))
    }
    if (!exists("render_author_sankey", mode="function")) {
      return(tags$div(class="small-note",
                      "Sankey module not loaded (missing sankey.R)."))
    }
    htmltools::tagList(render_author_sankey(dd$edges, dd$nodes))
  })

  output$txt_sankey_code <- shiny::renderText({
    dom <- input$sankey_domain %||% "Author"
    dd  <- .get_domain(dom)
    e   <- dd$edges
    n   <- dd$nodes

    if (is.null(e) || is.null(n) || !is.data.frame(e) || !is.data.frame(n) || nrow(e) == 0) {
      return(sprintf("SankeyMATIC: no data yet for domain '%s' (run PubMed or upload; then build domain Top-20).", dom))
    }

    # Normalize edges
    if (all(c("Leader","follower","WCD") %in% names(e))) {
      e2 <- dplyr::transmute(e, source = as.character(Leader), target = as.character(follower), w = suppressWarnings(as.numeric(WCD)))
    } else if (all(c("source","target","value") %in% names(e))) {
      e2 <- dplyr::transmute(e, source = as.character(source), target = as.character(target), w = suppressWarnings(as.numeric(value)))
    } else if (ncol(e) >= 3) {
      e2 <- e[,1:3,drop=FALSE]
      colnames(e2) <- c("source","target","w")
      e2$source <- as.character(e2$source); e2$target <- as.character(e2$target)
      e2$w <- suppressWarnings(as.numeric(e2$w))
    } else {
      return("SankeyMATIC: edges require Leader/follower/WCD (or source/target/value).")
    }
    e2 <- e2 %>% dplyr::filter(nzchar(source), nzchar(target), !is.na(w), w > 0)

    # SankeyMATIC link lines (force black links, as requested)
    link_lines <- paste0(e2$source, " [", signif(e2$w, 12), "] ", e2$target, " #000000")

    # Node color lines: prefer explicit node color; else cluster color; else value intensity
    node_lines <- character()

    n2 <- n
    if (!("name" %in% names(n2))) {
      # best-effort: first column as name
      n2$name <- as.character(n2[[1]])
    } else {
      n2$name <- as.character(n2$name)
    }

    # 1) explicit node color
    if ("color" %in% names(n2)) {
      colv <- as.character(n2$color)
      ok <- grepl("^#?[0-9A-Fa-f]{6}$", colv)
      colv[ok] <- ifelse(substr(colv[ok],1,1)=="#", colv[ok], paste0("#", colv[ok]))
      node_lines <- paste0(": ", n2$name[ok], " ", toupper(colv[ok]))
    }

    # 2) cluster-based color (prefer nodes$carac; else cluster/group)
    if (length(node_lines) == 0 && ("carac" %in% names(n2) || "cluster" %in% names(n2) || "group" %in% names(n2))) {
      cl <- if ("carac" %in% names(n2)) n2$carac else if ("cluster" %in% names(n2)) n2$cluster else n2$group
      cl_num <- suppressWarnings(as.integer(as.character(cl)))
      cl_key <- ifelse(!is.na(cl_num), as.character(cl_num), as.character(cl))
      lv <- unique(cl_key[!is.na(cl_key) & nzchar(cl_key)])
      lv_num <- suppressWarnings(as.integer(lv))
      if (all(!is.na(lv_num))) lv <- as.character(sort(lv_num)) else lv <- sort(lv)

      # palette: fixed series by cluster number, recycle if needed
      pal <- specified_colors
      cols <- if (all(!is.na(lv_num))) {
        setNames(vapply(as.integer(lv), function(k) pal[(k - 1L) %% length(pal) + 1L], character(1)), lv)
      } else {
        setNames(rep(pal, length.out = length(lv)), lv)
      }

      node_lines <- paste0(": ", n2$name, " ", cols[cl_key])
    }

# 3) value intensity (last resort)
    if (length(node_lines) == 0 && ("value" %in% names(n2))) {
      nv <- suppressWarnings(as.numeric(n2$value)); nv[is.na(nv)] <- 0
      s <- log1p(nv)
      mx <- max(s, na.rm = TRUE)
      idx <- if (is.finite(mx) && mx > 0) s / mx else rep(0.1, length(s))
      pal <- grDevices::colorRampPalette(c("#f7fbff","#08306b"))(100)
      pick <- pal[pmax(1, pmin(100, 1 + floor(idx * 99)))]
      hex <- toupper(substr(grDevices::rgb(t(grDevices::col2rgb(pick))/255), 2, 7))
      node_lines <- paste0(": ", n2$name, " #", hex)
    }

    paste(c(link_lines, node_lines), collapse = "\n")
  })


  # -------- World map (Country distribution, pre-FLCA) --------
  output$country_map <- shiny::renderPlot({
    # derive counts from rv$countries_list (preferred) else rv$country_nodes
    counts <- NULL
    if (!is.null(rv$countries_list)) {
      x <- rv$countries_list
      if (is.data.frame(x)) {
        # try common column names
        nm <- intersect(c("name","term","country","Country"), names(x))
        ct <- intersect(c("value","count","n","freq","Count"), names(x))
        if (length(nm) >= 1 && length(ct) >= 1) {
          counts <- dplyr::transmute(x, country = as.character(.data[[nm[1]]]), n = suppressWarnings(as.numeric(.data[[ct[1]]])))
        }
      } else if (is.list(x)) {
        v <- unlist(x, use.names = FALSE)
        v <- trimws(as.character(v))
        v <- v[nzchar(v)]
        tb <- sort(table(v), decreasing = TRUE)
        counts <- data.frame(country = names(tb), n = as.numeric(tb), stringsAsFactors = FALSE)
      }
    }
    if (is.null(counts) && !is.null(rv$country_nodes) && is.data.frame(rv$country_nodes)) {
      if (all(c("name","value") %in% names(rv$country_nodes))) {
        counts <- dplyr::transmute(rv$country_nodes, country = as.character(name), n = suppressWarnings(as.numeric(value)))
      }
    }

    if (is.null(counts) || nrow(counts) == 0) {
      plot.new(); text(0.5, 0.5, "World map: no country data yet.\nRun PubMed / upload to build Country term-list."); return()
    }

    world <- ggplot2::map_data("world")
    counts$country2 <- tolower(trimws(counts$country))
    world$region2 <- tolower(trimws(world$region))

    # aggregate
    agg <- counts %>% dplyr::group_by(country2) %>% dplyr::summarise(n = sum(n, na.rm = TRUE), .groups="drop")
    world2 <- dplyr::left_join(world, agg, by = c("region2" = "country2"))

    # make 0/NA light, non-zero use log scale
    world2$n0 <- world2$n
    world2$n0[is.na(world2$n0)] <- 0
    world2$fillv <- ifelse(world2$n0 <= 0, NA, log10(world2$n0 + 1))

    ggplot2::ggplot(world2, ggplot2::aes(long, lat, group = group)) +
      ggplot2::geom_polygon(ggplot2::aes(fill = fillv), color = "white", linewidth = 0.1) +
      ggplot2::coord_fixed(1.3) +
      ggplot2::theme_void() +
      ggplot2::labs(fill = "Count (log10)") +
      ggplot2::scale_fill_gradient(na.value = "grey95")
  })

  # -------- USA map (States, pre-FLCA) --------
  output$usa_map_placeholder <- shiny::renderUI({
    if (requireNamespace("plotly", quietly = TRUE)) return(tags$div(class="small-note", ""))
    tags$div(class="small-note",
             "USA map requires package: plotly. ",
             tags$code("install.packages('plotly')"))
  })

  output$usa_map <- plotly::renderPlotly({
    req(requireNamespace("plotly", quietly = TRUE))

    # build counts from rv$stateprov_list or rv$stateprov_nodes
    sp <- NULL
    if (!is.null(rv$stateprov_list)) {
      x <- rv$stateprov_list
      if (is.data.frame(x)) {
        nm <- intersect(c("name","term","state","State"), names(x))
        ct <- intersect(c("value","count","n","freq","Count"), names(x))
        if (length(nm) >= 1 && length(ct) >= 1) {
          sp <- dplyr::transmute(x, st = as.character(.data[[nm[1]]]), n = suppressWarnings(as.numeric(.data[[ct[1]]])))
        }
      } else if (is.list(x)) {
        v <- unlist(x, use.names = FALSE)
        v <- trimws(as.character(v))
        v <- v[nzchar(v)]
        tb <- sort(table(v), decreasing = TRUE)
        sp <- data.frame(st = names(tb), n = as.numeric(tb), stringsAsFactors = FALSE)
      }
    }
    if (is.null(sp) && !is.null(rv$stateprov_nodes) && is.data.frame(rv$stateprov_nodes)) {
      if (all(c("name","value") %in% names(rv$stateprov_nodes))) {
        sp <- dplyr::transmute(rv$stateprov_nodes, st = as.character(name), n = suppressWarnings(as.numeric(value)))
      }
    }

    if (is.null(sp) || nrow(sp) == 0) {
      sp <- data.frame(code = state.abb, state = state.name, n = 0)
    } else {
      sp$st <- trimws(sp$st)
      sp$code <- toupper(sp$st)
      # if full state name, map to abbreviation
      name_map <- setNames(state.abb, tolower(state.name))
      sp$code <- ifelse(nchar(sp$code) == 2, sp$code, name_map[tolower(sp$st)])
      sp$code[is.na(sp$code)] <- ""
      sp <- sp %>% dplyr::filter(nchar(code) == 2) %>% dplyr::group_by(code) %>% dplyr::summarise(n = sum(n, na.rm=TRUE), .groups="drop")
      sp$state <- state.name[match(sp$code, state.abb)]
    }

    sp$hover <- paste0(sp$state, "<br>Count: ", sp$n)

    g <- list(scope = "usa",
              projection = list(type = "albers usa"),
              showlakes = TRUE,
              lakecolor = "white")

    fig <- plotly::plot_geo(sp, locationmode = "USA-states") %>%
      plotly::add_trace(z = ~n, text = ~hover, locations = ~code, color = ~n, colors = "Purples") %>%
      plotly::colorbar(title = "Count") %>%
      plotly::layout(title = "USA (States) - pre-FLCA", geo = g)
    fig
  })

  output$tbl_usa_counts <- DT::renderDT({
    sp <- NULL
    if (!is.null(rv$stateprov_nodes) && is.data.frame(rv$stateprov_nodes) && all(c("name","value") %in% names(rv$stateprov_nodes))) {
      sp <- rv$stateprov_nodes %>% dplyr::transmute(state = as.character(name), count = suppressWarnings(as.numeric(value))) %>% dplyr::arrange(dplyr::desc(count))
    }
    if (is.null(sp)) sp <- data.frame(state=character(), count=numeric())
    DT::datatable(sp, options=list(pageLength=10), rownames=FALSE)
  })


  # -----------------------------
  # WebZIP: serve download/ folder
  # -----------------------------
  dl_dir0 <- .get_dl_dir(); dir.create(dl_dir0, recursive = TRUE, showWarnings = FALSE)
  shiny::addResourcePath("download", dl_dir0)
observeEvent(input$btn_refresh_webzip, {
  dl_dir <- .webzip_dl_dir()
  dir.create(dl_dir, recursive = TRUE, showWarnings = FALSE)

  # marker to prove refresh worked
  writeLines(paste("refreshed:", Sys.time()), file.path(dl_dir, "_refresh_ok.txt"), useBytes = TRUE)

  # Auto-run key domains if ready (so refresh can generate figures after PubMed fetch)
  seed <- input$seed %||% 1
  # Author is ready after PubMed fetch
  if (is.null(rv$author_nodes) || !is.data.frame(rv$author_nodes) || nrow(rv$author_nodes)==0) {
    # do nothing; Author nodes built in main run
  }
  # Run these if cached term-lists exist
  try(run_one_domain("Country", seed), silent=TRUE)
  try(run_one_domain("Institute", seed), silent=TRUE)
  try(run_one_domain("MeSH", seed), silent=TRUE)

  # Export essential figures (4 domains × 4 plots = 16 files; Sankey uses placeholder unless you upload real sankey.R)
  .webzip_export_domain_figs("Author",    rv$author_nodes,  rv$author_edges,  rv, dl_dir)
  .webzip_export_domain_figs("Country",   rv$country_nodes, rv$country_edges, rv, dl_dir)
  .webzip_export_domain_figs("Institute", rv$inst_nodes,    rv$inst_edges,    rv, dl_dir)
  .webzip_export_domain_figs("MeSH",      rv$mesh_nodes,    rv$mesh_edges,    rv, dl_dir)

  # Rebuild offline index.html (be robust to signature differences)
  if (exists(".write_offline_index", mode="function")) {
    fm <- tryCatch(names(formals(.write_offline_index)), error=function(e) character(0))
    if ("rv" %in% fm && "title" %in% fm) {
      try(.write_offline_index(download_dir = dl_dir, rv = rv, title="Author profile analysis (offline)"), silent=TRUE)
    } else if ("rv" %in% fm) {
      try(.write_offline_index(download_dir = dl_dir, rv = rv), silent=TRUE)
    } else {
      try(.write_offline_index(dl_dir), silent=TRUE)
    }
  }

  rv$webzip_stamp <- Sys.time()
}, ignoreInit = TRUE)


  output$webzip_file_list <- shiny::renderUI({
    dl_dir <- .get_dl_dir()
    if (!dir.exists(dl_dir)) return(shiny::tags$em("download/ folder not found."))
    files <- .list_download_files(dl_dir)
    if (length(files) == 0) return(shiny::tags$em("download/ is empty. Run analysis first."))
    shiny::tags$ul(lapply(files, function(f) shiny::tags$li(shiny::tags$a(href = paste0("download/", f), target="_blank", f))))
  })

  output$webzip_index_preview <- shiny::renderUI({
    dl_dir <- .get_dl_dir()
    idx <- file.path(dl_dir, "index.html")
    if (!file.exists(idx)) return(shiny::tags$em("index.html not found yet. Click Refresh."))
    shiny::tags$iframe(
      src = "download/index.html",
      style = "width:100%; height:520px; border:1px solid #ddd; border-radius:10px;"
    )
  })

  output$webzip_manifest <- shiny::renderTable({
    dl_dir <- .get_dl_dir()
    data.frame(path = .list_download_files(dl_dir), stringsAsFactors = FALSE)
  })

  # WebZIP download: zip whole download/
  output$dl_webzip <- shiny::downloadHandler(
    filename = function() paste0("webzip_", format(Sys.Date(), "%Y%m%d"), ".zip"),
    contentType = "application/zip",
    content = function(file) {
      dl_dir <- .get_dl_dir()
      if (!dir.exists(dl_dir)) stop("download/ folder not found")
      .write_offline_index(dl_dir, title="Author profile analysis (offline)")

      files_full <- list.files(dl_dir, recursive = TRUE, full.names = TRUE, all.files = TRUE, no.. = TRUE)
      if (length(files_full) == 0) stop("download/ is empty, nothing to zip")

      if (requireNamespace("zip", quietly = TRUE)) {
        zip::zipr(zipfile = file, files = files_full, root = dl_dir)
      } else {
        old <- getwd(); on.exit(setwd(old), add = TRUE)
        setwd(dl_dir)
        utils::zip(zipfile = file, files = list.files(".", recursive = TRUE, all.files = TRUE, no.. = TRUE))
      }
    }
  )



  # --- SAFE call: pass only arguments supported by current render_panel() ---
  safe_call_render_panel <- function(...) {
    if (!exists("render_panel", mode="function")) return(NULL)
    args <- list(...)
    fn <- get("render_panel", mode="function")
    keep <- intersect(names(args), names(formals(fn)))
    do.call(fn, args[keep])
  }


# ---- URL query -> input defaults (cmc/title/term) ----
observeEvent(session$clientData$url_search, {
  q <- shiny::parseQueryString(session$clientData$url_search)
  # accept both Title/title, Term/term, CMC/cmc
  qn <- setNames(as.list(q), tolower(names(q)))
  if (!is.null(qn$cmc))   updateTextInput(session, "cmc",   value = qn$cmc)
  if (!is.null(qn$title)) updateTextInput(session, "title", value = qn$title)
  if (!is.null(qn$term))  updateTextAreaInput(session, "term", value = qn$term)
}, once = TRUE)



# ---- URL autorun: ?term=...&autorun=1&cmc=test|<10digits> ----
observeEvent(session$clientData$url_search, {
  q <- shiny::parseQueryString(session$clientData$url_search)
  qn <- setNames(as.list(q), tolower(names(q)))
  auto <- as.character(qn$autorun %||% "")
  if (!nzchar(auto) || auto != "1") return()
  if (is.null(qn$term) || !nzchar(as.character(qn$term))) return()

  ip <- get_client_ip()
  cmc_param <- qn$cmc %||% ""
  if (!.gate_url_autorun(ip, cmc_param)) return()

  # record trial / cmc immediately (so repeated refresh won't keep granting)
  try({ .upsert_ip(ip, as.character(cmc_param %||% ""), inc_count = TRUE) }, silent=TRUE)

  # trigger the same run button
  session$sendCustomMessage("jsClick", "run")
}, once = TRUE)



  # ===== Summary (Top 10 + AAC) =====
  summary_domain <- reactiveVal("Author")

  .summary_domains <- c("Country","Journal","Year","Institute","Department","State/Province","Author","MeSH")

  .get_domain_objects <- function(dom){
    dom <- as.character(dom %||% "")
    if (identical(dom, "Author")) {
      list(nodes = rv$author_nodes, edges = rv$author_edges,
           AAC_value = rv$AAC_value_author, AAC_value2 = rv$AAC_value2_author,
           AAC_ss = rv$AAC_ss_author, AAC_a_star = rv$AAC_a_star_author)
    } else if (identical(dom, "Country")) {
      list(nodes = rv$country_nodes, edges = rv$country_edges,
           AAC_value = rv$AAC_value_country, AAC_value2 = rv$AAC_value2_country,
           AAC_ss = rv$AAC_ss_country, AAC_a_star = rv$AAC_a_star_country)
    } else if (identical(dom, "State/Province")) {
      list(nodes = rv$stateprov_nodes, edges = rv$stateprov_edges,
           AAC_value = rv$AAC_value_stateprov, AAC_value2 = rv$AAC_value2_stateprov,
           AAC_ss = rv$AAC_ss_stateprov, AAC_a_star = rv$AAC_a_star_stateprov)
    } else if (identical(dom, "Institute")) {
      list(nodes = rv$inst_nodes, edges = rv$inst_edges,
           AAC_value = rv$AAC_value_inst, AAC_value2 = rv$AAC_value2_inst,
           AAC_ss = rv$AAC_ss_inst, AAC_a_star = rv$AAC_a_star_inst)
    } else if (identical(dom, "Department")) {
      list(nodes = rv$dept_nodes, edges = rv$dept_edges,
           AAC_value = rv$AAC_value_dept, AAC_value2 = rv$AAC_value2_dept,
           AAC_ss = rv$AAC_ss_dept, AAC_a_star = rv$AAC_a_star_dept)
    } else if (identical(dom, "MeSH")) {
      list(nodes = rv$mesh_nodes, edges = rv$mesh_edges,
           AAC_value = rv$AAC_value_mesh, AAC_value2 = rv$AAC_value2_mesh,
           AAC_ss = rv$AAC_ss_mesh, AAC_a_star = rv$AAC_a_star_mesh)
    } else { # Journal/Year
      list(nodes = rv$jy_nodes, edges = rv$jy_edges,
           AAC_value = rv$AAC_value_jy, AAC_value2 = rv$AAC_value2_jy,
           AAC_ss = rv$AAC_ss_jy, AAC_a_star = rv$AAC_a_star_jy)
    }
  }

  output$summary_btns <- renderUI({
    doms <- .summary_domains
    btn <- function(id, lab){
      actionButton(id, lab, width="100%", class="btn btn-default")
    }
    # 2-column grid like your screenshot
    tags$div(
      style="display:grid; grid-template-columns: 1fr 1fr; gap: 12px; max-width: 720px;",
      btn("sum_country","Country"),
      btn("sum_journal","Journal"),
      btn("sum_inst","Institute"),
      btn("sum_year","Year"),
      btn("sum_dept","Department"),
      btn("sum_stateprov","State/Province"),
      btn("sum_author","Author"),
      btn("sum_mesh","MeSH Term")
    )
  })

  # Clicks just switch the view; if nodes are missing and the run buttons exist, auto-run that domain.
  observeEvent(input$sum_author,   { summary_domain("Author") }, ignoreInit = TRUE)
  observeEvent(input$sum_country,  {
    summary_domain("Country")
    if (is.null(rv$country_nodes) || is.null(rv$country_edges)) {
      if (!is.null(input$run_country)) run_one_domain("Country", input$seed)
    }
  }, ignoreInit = TRUE)
  observeEvent(input$sum_stateprov, {
    summary_domain("State/Province")
    if (is.null(rv$stateprov_nodes) || is.null(rv$stateprov_edges)) {
      if (!is.null(input$run_stateprov)) run_one_domain("State/Province", input$seed)
    }
  }, ignoreInit = TRUE)
  observeEvent(input$sum_inst, {
    summary_domain("Institute")
    if (is.null(rv$inst_nodes) || is.null(rv$inst_edges)) {
      if (!is.null(input$run_inst)) run_one_domain("Institute", input$seed)
    }
  }, ignoreInit = TRUE)
  observeEvent(input$sum_dept, {
    summary_domain("Department")
    if (is.null(rv$dept_nodes) || is.null(rv$dept_edges)) {
      if (!is.null(input$run_dept)) run_one_domain("Department", input$seed)
    }
  }, ignoreInit = TRUE)
  observeEvent(input$sum_mesh, {
    summary_domain("MeSH")
    if (is.null(rv$mesh_nodes) || is.null(rv$mesh_edges)) {
      if (!is.null(input$run_mesh)) run_one_domain("MeSH", input$seed)
    }
  }, ignoreInit = TRUE)
  observeEvent(input$sum_journal, {
    summary_domain("Journal")
    # data come from Journal/Year nodes (split by label)
  }, ignoreInit = TRUE)
  observeEvent(input$sum_year, {
    summary_domain("Year")
  }, ignoreInit = TRUE)

  output$summary_status <- renderUI({
    dom <- summary_domain()

    # Simple status: show AAC(value) based on displayed values
    aac_val <- NA
    n_nodes <- NA
    n_edges <- NA

    if (dom %in% c("Journal","Year")) {
      obj <- .get_domain_objects("Journal/Year")
      n_nodes <- if (is.null(obj$nodes)) 0 else nrow(as.data.frame(obj$nodes))
      n_edges <- if (is.null(obj$edges)) 0 else nrow(as.data.frame(obj$edges))
      x <- .get_domain_summary(dom)
      aac_val <- x$AAC
    } else if (dom == "Year") {
      x <- .get_domain_summary("Year")
      aac_val <- x$AAC
      n_nodes <- if (is.null(x$df)) 0 else nrow(x$df)
      n_edges <- 0
    } else {
      obj <- .get_domain_objects(dom)
      n_nodes <- if (is.null(obj$nodes)) 0 else nrow(as.data.frame(obj$nodes))
      n_edges <- if (is.null(obj$edges)) 0 else nrow(as.data.frame(obj$edges))
      x <- .get_domain_summary(dom)
      aac_val <- x$AAC
    }

    tags$div(
      tags$b("Selected domain: "), dom, tags$br(),
      tags$b("Top nodes: "), n_nodes, "  |  ", tags$b("Top edges: "), n_edges, tags$br(),
      tags$b("AAC: "), signif(as_scalar_or_na(aac_val, NA), 4)
    )
  })


  output$tbl_summary_top10 <- DT::renderDT({
    dom <- summary_domain()
    obj <- .get_domain_objects(dom)

    # For the table view, show only Top 5 name + value (simple)
    # Special handling: Journal and Year should be split from Journal/Year nodes
    get_simple_df <- function(dom){
      if (dom %in% c("Journal","Year")) {
        obj2 <- .get_domain_objects("Journal/Year")
        nd <- obj2$nodes
        if (!is.null(nd) && is.data.frame(nd) && nrow(nd) > 0 && all(c("name","value") %in% names(nd))) {
          df <- nd[, c("name","value"), drop=FALSE]
          df$name  <- as.character(df$name)
          df$value <- suppressWarnings(as.numeric(as.character(df$value)))
          is_year <- grepl("^[12][0-9]{3}$", trimws(df$name))
          if (dom == "Journal") df <- df[!is_year, , drop=FALSE]
          if (dom == "Year")    df <- df[ is_year, , drop=FALSE]
          df <- df[order(df$value, decreasing=TRUE, na.last=TRUE), , drop=FALSE]
          return(head(df, 5))
        }
        return(NULL)
      }

      nd <- obj$nodes
      if (is.null(nd) || !is.data.frame(nd) || nrow(nd) == 0 || !all(c("name","value") %in% names(nd))) return(NULL)
      df <- nd[, c("name","value"), drop=FALSE]
      df$name  <- as.character(df$name)
      df$value <- suppressWarnings(as.numeric(as.character(df$value)))
      df <- df[order(df$value, decreasing=TRUE, na.last=TRUE), , drop=FALSE]
      head(df, 5)
    }

    df <- get_simple_df(dom)

    if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) {
      return(DT::datatable(data.frame(Message="No data yet. Run this domain first."),
                           rownames = FALSE, options = list(dom='t')))
    }

    DT::datatable(df, rownames = FALSE, options = list(pageLength = 5, dom='t'))
  })



  
  
  # ===== Summary PNG export (Composite: 8 domains on 1 page) =====
  # Goal: simple "Top 5 + AAC" using *value only* (like the reference image)
  .summary_domains8 <- c("Country","Journal","Institute","Year","Department","State/Province","Author","MeSH")

  .is_year_label <- function(x){
    x <- trimws(as.character(x))
    grepl("^[12][0-9]{3}$", x)
  }

  .get_domain_summary <- function(dom){
    # Returns: list(df = data.frame(name,value), AAC = numeric)
    dom <- as.character(dom)

    # ---- Journal / Year split from existing "Journal/Year" nodes ----
    if (dom %in% c("Journal","Year")) {
      obj_jy <- .get_domain_objects("Journal/Year")
      nodes  <- obj_jy$nodes

      if (!is.null(nodes) && is.data.frame(nodes) && nrow(nodes) > 0 && all(c("name","value") %in% names(nodes))) {
        df <- nodes[, c("name","value"), drop=FALSE]
        df$name  <- as.character(df$name)
        df$value <- suppressWarnings(as.numeric(as.character(df$value)))

        if (dom == "Journal") df <- df[! .is_year_label(df$name), , drop=FALSE]
        if (dom == "Year")    df <- df[  .is_year_label(df$name), , drop=FALSE]

        df <- df[order(df$value, decreasing=TRUE, na.last=TRUE), , drop=FALSE]
        df <- head(df, 5)

        return(list(
          df=df,
          AAC=safe_AAC(df$value)
        ))
      }

      # If we don't have Journal/Year nodes yet, Year falls back to year/articles table
      if (dom == "Year") {
        ya <- NULL
        if (!is.null(rv$year_articles) && is.data.frame(rv$year_articles)) {
          ya <- rv$year_articles
        } else if (!is.null(rv$year_df) && is.data.frame(rv$year_df)) {
          ya <- rv$year_df
        }
        if (is.null(ya) || !nrow(ya)) return(list(df=NULL, AAC=NA))
        nms <- names(ya)
        ycol <- intersect(c("year","Year","pub_year"), nms)
        ccol <- intersect(c("n","count","freq","N","Articles","articles"), nms)
        ycol <- if (length(ycol)) ycol[1] else nms[1]
        ccol <- if (length(ccol)) ccol[1] else nms[min(2, length(nms))]
        df <- ya[, c(ycol, ccol), drop=FALSE]
        names(df) <- c("name","value")
        df$name  <- as.character(df$name)
        df$value <- suppressWarnings(as.numeric(as.character(df$value)))
        df <- df[order(df$value, decreasing=TRUE, na.last=TRUE), , drop=FALSE]
        df <- head(df, 5)
        return(list(df=df, AAC=safe_AAC(df$value)))
      }

      return(list(df=NULL, AAC=NA))
    }

    # ---- other domains (use their nodes directly) ----
    obj <- .get_domain_objects(dom)
    nodes <- obj$nodes
    if (is.null(nodes) || !is.data.frame(nodes) || nrow(nodes) == 0 || !all(c("name","value") %in% names(nodes))) {
      return(list(df=NULL, AAC=as_scalar_or_na(obj$AAC_value, NA)))
    }

    df <- nodes[, c("name","value"), drop=FALSE]
    df$name  <- as.character(df$name)
    df$value <- suppressWarnings(as.numeric(as.character(df$value)))
    df <- df[order(df$value, decreasing=TRUE, na.last=TRUE), , drop=FALSE]
    df <- head(df, 5)

    list(
      df=df,
      AAC=safe_AAC(df$value)
    )
  }

  .draw_summary_panel <- function(dom, df, aac){
    plot.new()
    title(main = dom, cex.main = 1.35, font.main = 2, col.main = "#d62728")

    if (is.null(df) || !is.data.frame(df) || !nrow(df)) {
      text(0.5, 0.55, "No data yet.\nRun this domain first.", cex = 1.05)
      mtext(paste0("AAC = ", signif(as_scalar_or_na(aac, NA), 4)), side=1, line=0.35, cex=1.0, font=2, col="#d62728")
      return(invisible(NULL))
    }

    # simple list layout
    y0 <- 0.84
    dy <- 0.13
    x_name <- 0.05
    x_val  <- 0.95

    for (i in seq_len(min(nrow(df), 5))) {
      yy <- y0 - (i-1)*dy

      nm <- as.character(df$name[i])
      if (nchar(nm) > 30) nm <- paste0(substr(nm, 1, 27), "…")
      val <- df$value[i]

      text(x_name, yy, nm, adj=0, cex=1.35, font=2)
      text(x_val,  yy, format(val, trim=TRUE, scientific=FALSE), adj=1, cex=1.35, font=2)
    }

    mtext(paste0("AAC = ", signif(as_scalar_or_na(aac, NA), 4)), side=1, line=0.35, cex=1.15, font=2, col="#d62728")
    invisible(NULL)
  }

  output$dl_summary_png <- downloadHandler(
    filename = function(){
      paste0("performance_summary_top5_aac_8domains.png")
    },
    contentType = "image/png",
    content = function(file){
      tryCatch({
        doms <- .summary_domains8

        png(file, width = 2400, height = 2700, res = 220)
        op <- par(no.readonly = TRUE)
        on.exit({par(op); dev.off()}, add=TRUE)

        par(mfrow = c(4, 2), mar = c(2.2, 1.3, 2.2, 0.9), oma = c(1.2, 1.2, 3.2, 1.2))

        for (d in doms) {
          x <- .get_domain_summary(d)
          .draw_summary_panel(d, x$df, x$AAC)
        }

        mtext("Summary report (Top 5 + AAC) - 8 domains", outer=TRUE, side=3, line=1, cex=1.55, font=2)
        invisible(NULL)

      }, error = function(e){
        # Always write a PNG (never HTML)
        png(file, width = 2000, height = 1200, res = 220)
        plot.new()
        title(main = "Summary PNG export failed", cex.main = 1.5)
        msg <- paste0(conditionMessage(e))
        text(0.5, 0.5, msg, cex = 1.0)
        dev.off()
      })
    }
  )
# ===== Report: domain summary (always show; robust names) =====
  output$tbl_report_domains <- renderTable({
  # Directly reference rv$... so Shiny tracks dependencies (no eval(parse()) here)
  make_row <- function(dom, nodes, edges,
                     max_items = NA, median_items = NA, n_ge2 = NA,
                     AAC_value = NA, AAC_value2 = NA,
                     AAC_ss = NA, AAC_a_star = NA) {

  median_items <- as_scalar_or_na(median_items, NA)
  n_ge2        <- as_scalar_or_na(n_ge2, NA)
  AAC_value    <- as_scalar_or_na(AAC_value, NA)
  AAC_value2   <- as_scalar_or_na(AAC_value2, NA)
  AAC_ss       <- as_scalar_or_na(AAC_ss, NA)
  AAC_a_star   <- as_scalar_or_na(AAC_a_star, NA)

  n_nodes <- if (is.null(nodes)) 0L else nrow(as.data.frame(nodes))
  n_edges <- if (is.null(edges)) 0L else nrow(as.data.frame(edges))

  data.frame(
    domain       = dom,
    n_nodes      = n_nodes,
    edge_number  = n_edges,
    max_items    = max_items,
    median_items = median_items,
    n_ge2        = n_ge2,
    AAC_value    = AAC_value,
    AAC_value2   = AAC_value2,
    AAC_ss       = AAC_ss,
    AAC_a_star   = AAC_a_star,
    stringsAsFactors = FALSE
  )
}

  rows <- list(
    make_row("Author", rv$author_nodes, rv$author_edges,
             rv$author_max_items, rv$author_median_items, rv$author_n_ge2,
             rv$AAC_value_author, rv$AAC_value2_author, rv$AAC_ss_author, rv$AAC_a_star_author),
    make_row("Country", rv$country_nodes, rv$country_edges,
             rv$country_max_items, rv$country_median_items, rv$country_n_ge2,
             rv$AAC_value_country, rv$AAC_value2_country, rv$AAC_ss_country, rv$AAC_a_star_country),
    make_row("State/Province", rv$stateprov_nodes, rv$stateprov_edges,
             rv$stateprov_max_items, rv$stateprov_median_items, rv$stateprov_n_ge2,
             rv$AAC_value_stateprov, rv$AAC_value2_stateprov, rv$AAC_ss_stateprov, rv$AAC_a_star_stateprov),
    make_row("Institute", rv$inst_nodes, rv$inst_edges,
             rv$inst_max_items, rv$inst_median_items, rv$inst_n_ge2,
             rv$AAC_value_inst, rv$AAC_value2_inst, rv$AAC_ss_inst, rv$AAC_a_star_inst),
    make_row("Department", rv$dept_nodes, rv$dept_edges,
             rv$dept_max_items, rv$dept_median_items, rv$dept_n_ge2,
             rv$AAC_value_dept, rv$AAC_value2_dept, rv$AAC_ss_dept, rv$AAC_a_star_dept),
    make_row("MeSH", rv$mesh_nodes, rv$mesh_edges,
             rv$mesh_max_items, rv$mesh_median_items, rv$mesh_n_ge2,
             rv$AAC_value_mesh, rv$AAC_value2_mesh, rv$AAC_ss_mesh, rv$AAC_a_star_mesh)
  )

  # Journal/Year as additional row (if objects exist)
  if (!is.null(rv$jy_nodes) || !is.null(rv$jy_edges) || !is.null(rv$journal_nodes) || !is.null(rv$journal_edges)){
    jn <- if (!is.null(rv$jy_nodes)) rv$jy_nodes else rv$journal_nodes
    je <- if (!is.null(rv$jy_edges)) rv$jy_edges else rv$journal_edges
    rows[[length(rows)+1]] <- make_row("Journal/Year", jn, je,
                                       rv$jy_max_items, rv$jy_median_items, rv$jy_n_ge2,
                                       rv$AAC_value_jy, rv$AAC_value2_jy, rv$AAC_ss_jy, rv$AAC_a_star_jy)
  }

  do.call(rbind, rows)
}, rownames = FALSE)

  # ---- Reactive wrappers (avoid accessing input$* outside reactive consumers) ----
  cmc <- reactive({ input$cmc %||% "" })
  term_rx <- reactive({ input$term %||% "" })


  # ---- Contact modal (no popup blockers) ----
  observeEvent(input$contact_btn, {
    showModal(modalDialog(
      title = "Contact authors / Request CMC",
      easyClose = TRUE,
      footer = modalButton("Close"),
      tags$p("Please choose one of the following ways to contact the authors."),
      tags$h4("1) LINE account"),
      tags$p("LINE Official Account ID:"), tags$code("@onq5657t"),
      tags$p(tags$a("Open LINE add-friend page", href="https://line.me/R/ti/p/%40onq5657t",
                    target="_blank", rel="noopener noreferrer")),
      tags$hr(),
      tags$h4("2) Email"),
      tags$p(tags$a("raschonline.service@gmail.com", href="mailto:raschonline.service@gmail.com")),
      tags$hr(),
      tags$h4("3) ChatGPT group"),
      tags$p("If you need a CMC, please send a short message with your name and purpose.")
    ))
  }, ignoreInit = TRUE)


  
  # (rv moved earlier)




  # ==== Save uploaded TXT to a physical file in project folder (added) ====
  observeEvent(input$pubmed_txt, {
    req(input$pubmed_txt)

# AMA/PubMed TXT upload is allowed without CMC.
# CMC remains required only for online PubMed Fetch/domain analysis.
    dest <- PERM_PUBMED_TXT
    ok <- tryCatch(file.copy(input$pubmed_txt$datapath, dest, overwrite = TRUE), error=function(e) { rv$log <- paste0(rv$log, "
[TXT] file.copy error: ", e$message); FALSE })
    if (ok) {
      rv$uploaded_perm_path <- dest
      rv$uploaded_txt_path <- dest  # backward compat
      rv$log <- paste0(rv$log, "\n[TXT] Saved uploaded file to: ", dest)
    } else {
      rv$log <- paste0(rv$log, "\n[TXT] Failed to save uploaded file.")
    }
  })

  # ==== Clear uploaded MEDLINE file (added) ====
  observeEvent(input$btn_clear_uploaded, {
    if (isTRUE(input$use_uploaded_txt) && file.exists(PERM_PUBMED_TXT)) {
      try(unlink(PERM_PUBMED_TXT), silent = TRUE)
    }
    rv$uploaded_perm_path <- NULL
    rv$uploaded_txt_path <- NULL
    rv$pmids_from_txt <- FALSE
    rv$log <- paste0(rv$log, "\n[TXT] Cleared uploaded file: ", basename(PERM_PUBMED_TXT))
  })

  # Optional: nodes+edges zip (if you later wire this button)
  output$dl_nodes_edges_zip <- downloadHandler(
    filename = function(){ paste0("nodes_edges_", format(Sys.Date(), "%Y%m%d"), ".zip") },
    contentType = "application/zip",
    content = function(file){
      dl_dir <- .get_dl_dir()
      dir.create(dl_dir, showWarnings = FALSE, recursive = TRUE)
      # Try best-effort export from rv
      if (!is.null(rv$nodes20)) utils::write.csv(rv$nodes20, file.path(dl_dir,"nodes.csv"), row.names = FALSE, fileEncoding="UTF-8")
      if (!is.null(rv$edges20)) utils::write.csv(rv$edges20, file.path(dl_dir,"edges.csv"), row.names = FALSE, fileEncoding="UTF-8")
      .write_offline_index(dl_dir, rv)
      .zip_download_dir(file, dl_dir)
    }
  )

  # ---- PMID query-result download and preview ----
  .pmids_current <- function(){
    pm <- rv$pmids
    if (is.null(pm) || length(pm) == 0) return(character())
    pm <- as.character(pm)
    pm <- pm[!is.na(pm) & nzchar(trimws(pm))]
    unique(pm)
  }

  .pmids_df <- function(){
    pm <- .pmids_current()
    data.frame(No = seq_along(pm), PMID = pm, stringsAsFactors = FALSE)
  }

  output$tbl_pmids <- DT::renderDT({
    pm <- .pmids_current()
    if (!length(pm)) {
      return(DT::datatable(data.frame(Message = "No PMIDs yet. Run a PubMed query or upload MEDLINE first."), options = list(dom = "t"), rownames = FALSE))
    }
    DT::datatable(.pmids_df(), options = list(pageLength = 20, scrollX = TRUE), rownames = FALSE)
  })

  output$ui_icite_link <- renderUI({
    pm <- .pmids_current()
    if (!length(pm)) return(tags$div(class = "small-note", "No PMIDs yet. Run Fetch PubMed first."))
    url <- rv$icite_url
    if (is.null(url) || !nzchar(url)) {
      url <- paste0("https://icite.od.nih.gov/analysis?pmids=", paste(pm[seq_len(min(length(pm), 900))], collapse = ","))
    }
    tags$div(
      tags$p(tags$b("Fetched PMIDs: "), length(pm)),
      tags$a(href = url, target = "_blank", "Open these PMIDs in NIH iCite"),
      tags$p(class = "small-note", "The iCite URL uses the first 900 PMIDs to keep the URL practical; the download buttons export the full fetched PMID list.")
    )
  })

  output$dl_pmids_txt <- downloadHandler(
    filename = function(){ paste0("pmids_query_result_", format(Sys.Date(), "%Y%m%d"), ".txt") },
    contentType = "text/plain",
    content = function(file){
      pm <- .pmids_current()
      shiny::validate(shiny::need(length(pm) > 0, "No PMIDs available. Run Fetch PubMed first."))
      writeLines(pm, file, useBytes = TRUE)
    }
  )

  output$dl_pmids_csv <- downloadHandler(
    filename = function(){ paste0("pmids_query_result_", format(Sys.Date(), "%Y%m%d"), ".csv") },
    contentType = "text/csv",
    content = function(file){
      pm <- .pmids_current()
      shiny::validate(shiny::need(length(pm) > 0, "No PMIDs available. Run Fetch PubMed first."))
      utils::write.csv(data.frame(PMID = pm, stringsAsFactors = FALSE), file, row.names = FALSE, fileEncoding = "UTF-8")
    }
  )


# ---- domain picker helpers (Kano/SSplot) ----
.domain_key <- function(dom){
  dom <- as.character(dom %||% "Author")
  if (dom == "Author") return(list(nd="author_nodes", ed="author_edges",
                                 aac_v="AAC_value_author", aac_v2="AAC_value2_author",
                                 aac_ss="AAC_ss_author", aac_a="AAC_a_star_author"))
  if (dom == "Journal/Year") return(list(nd="jy_nodes", ed="jy_edges",
                                       aac_v="AAC_value_jy", aac_v2="AAC_value2_jy",
                                       aac_ss="AAC_ss_jy", aac_a="AAC_a_star_jy"))
  if (dom == "Country") return(list(nd="country_nodes", ed="country_edges",
                                  aac_v="AAC_value_country", aac_v2="AAC_value2_country",
                                  aac_ss="AAC_ss_country", aac_a="AAC_a_star_country"))
  if (dom == "Institute") return(list(nd="inst_nodes", ed="inst_edges",
                                    aac_v="AAC_value_inst", aac_v2="AAC_value2_inst",
                                    aac_ss="AAC_ss_inst", aac_a="AAC_a_star_inst"))
  if (dom == "Department") return(list(nd="dept_nodes", ed="dept_edges",
                                     aac_v="AAC_value_dept", aac_v2="AAC_value2_dept",
                                     aac_ss="AAC_ss_dept", aac_a="AAC_a_star_dept"))
  if (dom == "MeSH") return(list(nd="mesh_nodes", ed="mesh_edges",
                              aac_v="AAC_value_mesh", aac_v2="AAC_value2_mesh",
                              aac_ss="AAC_ss_mesh", aac_a="AAC_a_star_mesh"))
  # fallback
  list(nd="author_nodes", ed="author_edges",
       aac_v="AAC_value_author", aac_v2="AAC_value2_author",
       aac_ss="AAC_ss_author", aac_a="AAC_a_star_author")
}

.get_domain <- function(dom){
  k <- .domain_key(dom)
  nd <- rv[[k$nd]]
  ed <- rv[[k$ed]]
  list(
    nodes = nd,
    edges = ed,
    AAC_value  = rv[[k$aac_v]]  %||% 0,
    AAC_value2 = rv[[k$aac_v2]] %||% 0,
    AAC_ss     = rv[[k$aac_ss]] %||% 0,
    AAC_a_star = rv[[k$aac_a]]  %||% 0
  )
}

  addlog <- function(...) {
    msg <- paste0(..., collapse = "")
    cat(msg, "\n")
    try(flush.console(), silent = TRUE)
    shiny::isolate({
      old <- rv$log
      if (is.null(old) || !nzchar(old)) rv$log <- msg else rv$log <- paste(old, msg, sep = "\n")
    })
  }

observeEvent(input$load_example, {
    key <- input$example_pick
    if (nzchar(key)) updateTextAreaInput(session, "term", value = examples[[key]])
  })

  # Author parsing (first-last)
  get_authors_per_article <- function(article_xml) {
    a_blocks <- regmatches(article_xml, gregexpr("<Author\\b[\\s\\S]*?</Author>", article_xml, perl=TRUE))[[1]]
    if (length(a_blocks) == 0) return(character())
    pick_tag <- function(x, tag) {
      m <- regexec(paste0("<", tag, ">([\\s\\S]*?)</", tag, ">"), x, perl=TRUE)
      r <- regmatches(x, m)[[1]]
      if (length(r) >= 2) r[2] else ""
    }
    authors <- vapply(a_blocks, function(b) {
      last <- pick_tag(b, "LastName")
      ini  <- pick_tag(b, "Initials")
      fore <- pick_tag(b, "ForeName")
      if (nzchar(last) && nzchar(ini)) return(paste(last, ini))
      if (nzchar(last) && nzchar(fore)) return(paste(last, fore))
      coll <- pick_tag(b, "CollectiveName")
      if (nzchar(coll)) return(coll)
      ""
    }, character(1))
    authors <- trimws(authors)
    authors <- authors[nzchar(authors)]
    unique(authors)
  }



  # ---- Series author AAC mode ---------------------------------------------
  .batch_split_author_names <- function(x){
    x <- as.character(x %||% "")
    x <- gsub("\r", "\n", x)
    parts <- unlist(strsplit(x, "\n|;|\\|", perl = TRUE), use.names = FALSE)
    parts <- trimws(parts)
    parts <- parts[nzchar(parts)]
    parts <- sub("\\[Author\\]$", "", parts, ignore.case = TRUE)
    parts <- gsub('^"|"$', "", trimws(parts))
    unique(parts[nzchar(parts)])
  }

  .batch_get_author_names_for_run <- function(){
    # When the AMA textarea switch is checked, use the textarea as the source.
    # This supports either one-author-per-line input or full AMA/PubMed references.
    if (isTRUE(input$use_ama_textarea)) {
      txt <- paste(input$ama_refs_text %||% "", collapse = "\n")
      txt <- enc2utf8(as.character(txt))
      txt <- gsub("\r\n|\r", "\n", txt, perl = TRUE)
      if (!nzchar(trimws(txt))) return(character())

      wide <- tryCatch(.ref_to_wide_from_text_strict(txt), error = function(e) NULL)
      if (is.data.frame(wide) && nrow(wide) > 0) {
        author_cols <- grep("^Author_", names(wide), value = TRUE)
        if (length(author_cols) > 0) {
          nm <- unlist(wide[author_cols], use.names = FALSE)
          nm <- trimws(as.character(nm))
          nm <- nm[nzchar(nm)]
          if (length(nm) > 0) return(unique(nm))
        }
      }
      return(.batch_split_author_names(txt))
    }

    .batch_split_author_names(input$term)
  }

  .batch_aac_label <- function(aac){
    if (!is.finite(aac)) return(NA_character_)
    if (aac < 0.60) return("Low dominance")
    if (aac < 0.67) return("Emerging dominance")
    if (aac < 0.70) return("Moderate dominance")
    if (aac < 0.75) return("High dominance")
    "Very high dominance"
  }

  .batch_fetch_xml_by_pmids <- function(pmids, batch_size = 200){
    if (!length(pmids)) return("")
    xml_parts <- character()
    idxs <- seq(1, length(pmids), by = batch_size)
    for (i in idxs) {
      chunk <- pmids[i:min(i + batch_size - 1, length(pmids))]
      x <- tryCatch(rentrez::entrez_fetch(db = "pubmed", id = chunk, rettype = "xml", parsed = FALSE),
                    error = function(e) "")
      if (nzchar(x)) xml_parts <- c(xml_parts, x)
      Sys.sleep(0.35)
    }
    paste(xml_parts, collapse = "\n")
  }

  .batch_compute_author_aac_one <- function(author_name, retmax = 500){
    author_name <- trimws(as.character(author_name))
    q <- if (grepl("\\[Author\\]", author_name, ignore.case = TRUE)) author_name else paste0(author_name, "[Author]")

    ss <- tryCatch(rentrez::entrez_search(db = "pubmed", term = q, retmax = retmax), error = function(e) e)
    if (inherits(ss, "error")) {
      return(data.frame(target_author = author_name, query = q, pubmed_hits = NA_integer_, fetched_pmids = 0L,
                        top1_author = NA_character_, top1_n = NA_integer_, top2_author = NA_character_, top2_n = NA_integer_,
                        top3_author = NA_character_, top3_n = NA_integer_, AAC = NA_real_, odds_ratio_r = NA_real_,
                        dominance_label = NA_character_, status = paste0("PubMed search failed: ", conditionMessage(ss)), stringsAsFactors = FALSE))
    }

    pmids <- ss$ids
    pubmed_hits <- suppressWarnings(as.numeric(ss$count))
    fetched_pmids <- length(pmids)
    if (!length(pmids)) {
      return(data.frame(target_author = author_name, query = q, pubmed_hits = as.integer(ss$count), fetched_pmids = 0L,
                        top1_author = NA_character_, top1_n = NA_integer_, top2_author = NA_character_, top2_n = NA_integer_,
                        top3_author = NA_character_, top3_n = NA_integer_, AAC = NA_real_, odds_ratio_r = NA_real_,
                        dominance_label = NA_character_, status = "No PMID fetched", stringsAsFactors = FALSE))
    }

    # AAC is valid only when all PubMed hits are actually fetched.
    # If PubMed has more hits than retmax/fetched_pmids, the author distribution is incomplete;
    # this is common for duplicated/ambiguous author names and would bias AAC.
    if (is.finite(pubmed_hits) && pubmed_hits > fetched_pmids) {
      return(data.frame(target_author = author_name, query = q, pubmed_hits = as.integer(pubmed_hits), fetched_pmids = fetched_pmids,
                        top1_author = NA_character_, top1_n = NA_integer_, top2_author = NA_character_, top2_n = NA_integer_,
                        top3_author = NA_character_, top3_n = NA_integer_, AAC = NA_real_, odds_ratio_r = NA_real_,
                        dominance_label = NA_character_,
                        status = sprintf("Skipped AAC: incomplete PubMed fetch (pubmed_hits=%s > fetched_pmids=%s). Increase retmax or refine the author query.",
                                         format(pubmed_hits, scientific = FALSE, trim = TRUE), fetched_pmids),
                        stringsAsFactors = FALSE))
    }

    xml_txt <- .batch_fetch_xml_by_pmids(pmids)
    articles <- regmatches(xml_txt, gregexpr("<PubmedArticle\\b[\\s\\S]*?</PubmedArticle>", xml_txt, perl = TRUE))[[1]]
    if (!length(articles)) {
      return(data.frame(target_author = author_name, query = q, pubmed_hits = as.integer(ss$count), fetched_pmids = length(pmids),
                        top1_author = NA_character_, top1_n = NA_integer_, top2_author = NA_character_, top2_n = NA_integer_,
                        top3_author = NA_character_, top3_n = NA_integer_, AAC = NA_real_, odds_ratio_r = NA_real_,
                        dominance_label = NA_character_, status = "No XML article blocks parsed", stringsAsFactors = FALSE))
    }

    fl <- unlist(lapply(articles, function(z){
      a <- get_authors_per_article(z)
      a <- trimws(as.character(a)); a <- a[nzchar(a)]
      if (!length(a)) return(character(0))
      if (length(a) == 1) return(a[1])
      unique(c(a[1], a[length(a)]))
    }), use.names = FALSE)
    fl <- trimws(fl); fl <- fl[nzchar(fl)]
    if (!length(fl)) {
      return(data.frame(target_author = author_name, query = q, pubmed_hits = as.integer(ss$count), fetched_pmids = length(pmids),
                        top1_author = NA_character_, top1_n = NA_integer_, top2_author = NA_character_, top2_n = NA_integer_,
                        top3_author = NA_character_, top3_n = NA_integer_, AAC = NA_real_, odds_ratio_r = NA_real_,
                        dominance_label = NA_character_, status = "No first/last authors extracted", stringsAsFactors = FALSE))
    }

    tb <- sort(table(fl), decreasing = TRUE)
    nm <- names(tb); vv <- as.numeric(tb)
    top_names <- c(nm, rep(NA_character_, 3))[1:3]
    top_vals  <- c(vv, rep(NA_real_, 3))[1:3]
    r <- NA_real_; aac <- NA_real_
    if (all(is.finite(top_vals[1:3])) && top_vals[2] > 0 && top_vals[3] > 0) {
      r <- (top_vals[1] / top_vals[2]) / (top_vals[2] / top_vals[3])
      aac <- r / (1 + r)
    }
    data.frame(target_author = author_name, query = q, pubmed_hits = as.integer(ss$count), fetched_pmids = length(pmids),
               top1_author = top_names[1], top1_n = as.integer(top_vals[1]),
               top2_author = top_names[2], top2_n = as.integer(top_vals[2]),
               top3_author = top_names[3], top3_n = as.integer(top_vals[3]),
               AAC = round(aac, 4), odds_ratio_r = round(r, 4), dominance_label = .batch_aac_label(aac),
               status = "OK", stringsAsFactors = FALSE)
  }

  output$dl_batch_author_aac <- downloadHandler(
    filename = function() paste0("series_author_AAC_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content = function(file) {
      df <- rv$batch_author_aac
      if (is.null(df) || !is.data.frame(df)) df <- data.frame(Message = "No series author AAC result yet. Check the box and click Fetch PubMed.")
      utils::write.csv(df, file, row.names = FALSE, fileEncoding = "UTF-8")
    }
  )

  # Build a domain network via: co-occurrence edges -> directed symmetric -> FLCA_run -> major sampling
  
  # Build Top20 nodes/edges for a metadata domain.
  # - If co-occurrence edges exist: run FLCA + major sampling (original behavior).
  # - If edges are empty: fall back to frequency nodes (so n_nodes won't be 0).
  build_term_freq_nodes <- function(term_list, top_n = 20) {
    if (is.null(term_list)) return(NULL)
    # term_list: list(article -> character vector of terms)
    flat <- unlist(lapply(term_list, function(x){
      x <- as.character(x)
      x <- trimws(x)
      x[nzchar(x)]
    }), use.names = FALSE)
    if (!length(flat)) return(NULL)
    tab <- sort(table(flat), decreasing = TRUE)
    top <- head(tab, top_n)
    data.frame(
      name = names(top),
      value = as.numeric(top),
      value2 = as.numeric(top),
      stringsAsFactors = FALSE
    )
  }

  build_domain_flca_top20 <- function(edges_undirected, term_list = NULL, seed = 1) {

    # ---- fallback: no edges, but still show nodes by frequency ----
    if (!is.data.frame(edges_undirected) || nrow(edges_undirected) == 0) {
      nodes <- build_term_freq_nodes(term_list, top_n = 20)
      if (is.null(nodes) || !nrow(nodes)) return(list(nodes = NULL, data = NULL))
      return(list(nodes = nodes, data = data.frame(Leader=character(), Follower=character(), WCD=numeric(), stringsAsFactors = FALSE)))
    }

    need_cols <- c("from","to","weight")
    if (!all(need_cols %in% names(edges_undirected))) {
      nodes <- build_term_freq_nodes(term_list, top_n = 20)
      if (is.null(nodes) || !nrow(nodes)) return(list(nodes = NULL, data = NULL))
      return(list(nodes = nodes, data = data.frame(Leader=character(), Follower=character(), WCD=numeric(), stringsAsFactors = FALSE)))
    }

    edges_undirected <- edges_undirected[, need_cols, drop=FALSE]
    edges_undirected$from   <- as.character(edges_undirected$from)
    edges_undirected$to     <- as.character(edges_undirected$to)
    edges_undirected$weight <- suppressWarnings(as.numeric(edges_undirected$weight))
    edges_undirected <- edges_undirected[is.finite(edges_undirected$weight) & edges_undirected$weight > 0, , drop=FALSE]
    if (!nrow(edges_undirected)) {
      nodes <- build_term_freq_nodes(term_list, top_n = 20)
      if (is.null(nodes) || !nrow(nodes)) return(list(nodes = NULL, data = NULL))
      return(list(nodes = nodes, data = data.frame(Leader=character(), Follower=character(), WCD=numeric(), stringsAsFactors = FALSE)))
    }

    # nodes
    items <- sort(unique(c(edges_undirected$from, edges_undirected$to)))
    w1 <- tapply(edges_undirected$weight, edges_undirected$from, sum, na.rm=TRUE)
    w2 <- tapply(edges_undirected$weight, edges_undirected$to,   sum, na.rm=TRUE)
    v1 <- w1[match(items, names(w1))]; v1[is.na(v1)] <- 0
    v2 <- w2[match(items, names(w2))]; v2[is.na(v2)] <- 0
    val <- v1 + v2
    nodes <- data.frame(name=items, value=as.numeric(val), value2=as.numeric(val), stringsAsFactors = FALSE)

    # ---- RUNNER: FLCA + MajorSampling + Silhouette(SS) + a* ----
    # runner expects edges: Leader/Follower/WCD (undirected OK)
    edges_lfw <- data.frame(
      Leader   = edges_undirected$from,
      Follower = edges_undirected$to,
      WCD      = as.numeric(edges_undirected$weight),
      stringsAsFactors = FALSE
    )

    cfg_runner <- list(
      top_clusters=5, base_per_cluster=4, target_n=20,
      intra_delta=2, inter_delta=5, eps=1e-9
    )

    addlog(sprintf('[DOMAIN] runner input: nodes=%d edges=%d', nrow(nodes), nrow(edges_lfw)))

    set.seed(seed)
    res <- tryCatch(run_flca_ms_sil_runner(nodes, edges_lfw, cfg_runner, verbose = TRUE), error=function(e) e)
    if (inherits(res, "error")) {
      addlog(paste0("[RUNNER ERROR] ", res$message))
      nodes_fb <- build_term_freq_nodes(term_list, top_n = 20)
      if (!is.null(nodes_fb) && nrow(nodes_fb)) {
        return(list(nodes = nodes_fb,
                    data  = data.frame(Leader=character(), Follower=character(), WCD=numeric(), stringsAsFactors = FALSE)))
      }
      return(list(nodes=NULL, data=NULL))
    }

    # Standardize runner outputs (nodes/data)
    out <- list(
      nodes = if (!is.null(res$nodes)) res$nodes else res$nodes,
      data  = res$data
    )

    if (is.null(out$nodes) || !is.data.frame(out$nodes) || !nrow(out$nodes)) {
      nodes_fb <- build_term_freq_nodes(term_list, top_n = 20)
      if (!is.null(nodes_fb) && nrow(nodes_fb)) {
        return(list(nodes = nodes_fb,
                    data  = data.frame(Leader=character(), Follower=character(), WCD=numeric(), stringsAsFactors = FALSE)))
      }
      return(list(nodes=NULL, data=NULL))
    }
    if (is.null(out$data) || !is.data.frame(out$data)) {
      out$data <- data.frame(Leader=character(), Follower=character(), WCD=numeric(), stringsAsFactors = FALSE)
    }

    addlog(sprintf('[RUNNER] done: nodes=%s edges=%s clusters=%s',
                   if(!is.null(out$nodes) && is.data.frame(out$nodes)) nrow(out$nodes) else NA,
                   if(!is.null(out$data) && is.data.frame(out$data)) nrow(out$data) else NA,
                   if(!is.null(out$nodes) && is.data.frame(out$nodes) && 'carac' %in% names(out$nodes)) length(unique(out$nodes$carac)) else NA))

    if (!("a_star1" %in% names(out$nodes)) && ("a_i" %in% names(out$nodes))) {
      out$nodes$a_star1 <- 1/(1+as.numeric(out$nodes$a_i))
    }

    try({
      ed_tmp <- out$data
      if (is.data.frame(ed_tmp) && nrow(ed_tmp) > 0) {
        if (!all(c("Leader","Follower") %in% names(ed_tmp))) {
          if (all(c("from","to") %in% names(ed_tmp))) names(ed_tmp)[match(c("from","to"), names(ed_tmp))] <- c("Leader","Follower")
        }
        if (!("WCD" %in% names(ed_tmp))) ed_tmp$WCD <- if ("value" %in% names(ed_tmp)) ed_tmp$value else 1
        sil_df <- compute_silhouette_df(out$nodes, ed_tmp)
        if (is.data.frame(sil_df) && nrow(sil_df) > 0) {
          # Robust key join: prefer full/original names if available (avoid display truncation mismatches)
          key_nodes <- if ("orig_name" %in% names(out$nodes)) out$nodes$orig_name else out$nodes$name
          key_sil   <- if ("orig_name" %in% names(sil_df))   sil_df$orig_name   else sil_df$name
          if (!is.null(key_nodes) && !is.null(key_sil)) {
            idx <- match(as.character(key_nodes), as.character(key_sil))
          } else if ("id" %in% names(out$nodes) && "id" %in% names(sil_df)) {
            idx <- match(out$nodes$id, sil_df$id)
          } else {
            idx <- rep(NA_integer_, nrow(out$nodes))
          }

          if ("ssi" %in% names(sil_df)) out$nodes$ssi <- sil_df$ssi[idx] else if ("sil_width" %in% names(sil_df)) out$nodes$ssi <- sil_df$sil_width[idx]
          if ("a_i" %in% names(sil_df))       out$nodes$a_i <- sil_df$a_i[idx]
          if ("b_i" %in% names(sil_df))       out$nodes$b_i <- sil_df$b_i[idx]
          if (!("a_star1" %in% names(out$nodes))) {
            if ("a_star1" %in% names(sil_df)) {
              out$nodes$a_star1 <- sil_df$a_star1[idx]
            } else if ("a_i" %in% names(out$nodes)) {
              out$nodes$a_star1 <- 1/(1+as.numeric(out$nodes$a_i))
            }
          }
        }
        
        if ("ssi" %in% names(out$nodes)) out$nodes$ssi[is.na(out$nodes$ssi) | !is.finite(out$nodes$ssi)] <- 0
        if ("a_star1" %in% names(out$nodes)) out$nodes$a_star1[is.na(out$nodes$a_star1) | !is.finite(out$nodes$a_star1)] <- 0
# Back-compat aliases for Report/AAC (some parts expect these names)
        if (!("SSi" %in% names(out$nodes)) && "ssi" %in% names(out$nodes)) out$nodes$SSi <- out$nodes$ssi
        if (!("a_star" %in% names(out$nodes)) && "a_star1" %in% names(out$nodes)) out$nodes$a_star <- out$nodes$a_star1
      }
    }, silent=TRUE)

    list(nodes=out$nodes, data=out$data)
  }
  


observeEvent(input$run, {
    # ---- Access gate (IP allowlist + trial + CMC) ----
ip_addr <- ipm_get_client_ip(session)
cmc_now  <- input$cmc %||% ""

# If CMC missing/invalid: force SoftwareX demo query (still counts as Trial for display)
force_demo <- !ipm_is_cmc_10(cmc_now)
if (isTRUE(force_demo)) {
  demo_term <- "SoftwareX[Journal]"
  try(updateTextAreaInput(session, "term", value = demo_term), silent = TRUE)
}

# Gate logic:
# - For forced demo: always allow (but still respect iplist allowlist if present), and record run in ip.txt
# - For valid CMC: apply iplist/trial gate via ipmodule
if (isTRUE(force_demo)) {
  allow <- tryCatch(ipm_read_iplist(app_dir = APP_DIR), error = function(e) character())
  if (length(allow) && !(ip_addr %in% allow)) {
    g <- list(ok = FALSE, reason = "ip_not_allowlisted", ip = ip_addr, policy = "iplist")
  } else {
    try(ipm_upsert_ip(ip = ip_addr, cmc = "", inc_count = TRUE, app_dir = APP_DIR), silent = TRUE)
    g <- list(ok = TRUE, reason = "demo", ip = ip_addr, policy = if (length(allow)) "iplist" else "trial")
  }
  rv$ip_access_type <- "Trial"
} else {
  g <- ipm_gate(ip = ip_addr, cmc = cmc_now, app_dir = APP_DIR, inc_count_on_allow = TRUE)
  rv$ip_access_type <- if (identical(g$policy, "iplist")) "IP pass" else if (identical(g$reason, "trial_first_time")) "Trial" else "CMC pass"
}

rv$ip_addr <- g$ip
# pull latest row from ip.txt for this IP (best effort)
ipdf <- tryCatch(ipm_read_ip_txt(app_dir = APP_DIR), error = function(e) NULL)
if (!is.null(ipdf) && nrow(ipdf)) {
  ii <- match(rv$ip_addr, ipdf$ip)
  if (is.finite(ii) && !is.na(ii)) {
    rv$ip_recent_date <- ipdf$recent_date[ii]
    rv$ip_total_count <- ipdf$total_count[ii]
    rv$ip_cmc <- ipdf$cmc[ii]
  }
}

if (!isTRUE(g$ok)) {
  if (identical(g$reason, "ip_not_allowlisted")) {
    showModal(modalDialog(
      title = "Access blocked",
      "Your IP is not in iplist.txt.",
      easyClose = TRUE,
      footer = modalButton("OK")
    ))
  } else {
    showModal(modalDialog(
      title = "CMC required",
      "This IP has already used the one-time trial.",
      "Please enter a 10-digit numeric CMC or contact the author.",
      easyClose = TRUE,
      footer = modalButton("OK")
    ))
  }
  return()
}

if (isTRUE(force_demo)) {
  showNotification('CMC missing/invalid → running demo query: "SoftwareX"[Journal].', type="message", duration=4)
}

# ---- Series author AAC mode: paste one author per line, append [Author], compute AAC CSV ----
if (isTRUE(input$batch_author_aac)) {
  author_source <- if (isTRUE(input$use_ama_textarea)) "AMA/PubMed textarea" else "PubMed query box"
  author_names <- .batch_get_author_names_for_run()
  if (!length(author_names)) {
    addlog("[BATCH AAC][STOP] No author names found. When the textarea option is checked, paste one author per line or full AMA/PubMed references in the textarea; otherwise paste one author per line in the PubMed query box.")
    return(NULL)
  }
  rv$log <- ""
  rv$batch_author_aac <- NULL
  addlog("[BATCH AAC] START | source=", author_source, " | authors=", length(author_names), " | retmax=", input$retmax, " | AAC guard=compute only when pubmed_hits <= fetched_pmids")
  res <- list()
  withProgress(message = "Series author AAC: fetching PubMed by author", value = 0, {
    for (ii in seq_along(author_names)) {
      nm <- author_names[ii]
      incProgress(1 / length(author_names), detail = paste0(ii, "/", length(author_names), ": ", nm))
      addlog("[BATCH AAC] ", ii, "/", length(author_names), " | ", nm, "[Author]")
      res[[ii]] <- .batch_compute_author_aac_one(nm, retmax = input$retmax)
      rr <- res[[ii]]
      if (is.data.frame(rr) && nrow(rr)) {
        addlog("[BATCH AAC][RESULT] ", nm,
               " | hits=", rr$pubmed_hits[1], " fetched=", rr$fetched_pmids[1],
               " | top1=", rr$top1_author[1], " n=", rr$top1_n[1],
               " | top2=", rr$top2_author[1], " n=", rr$top2_n[1],
               " | top3=", rr$top3_author[1], " n=", rr$top3_n[1],
               " | AAC=", rr$AAC[1], " | ", rr$dominance_label[1],
               " | status=", rr$status[1])
      }
    }
  })
  rv$batch_author_aac <- do.call(rbind, res)

  # Also write a physical CSV into ./download/ so it appears in the download folder and WebZIP.
  dl_dir <- .get_dl_dir()
  dir.create(dl_dir, recursive = TRUE, showWarnings = FALSE)
  batch_csv_path <- file.path(dl_dir, "series_author_AAC.csv")
  utils::write.csv(rv$batch_author_aac, batch_csv_path, row.names = FALSE, fileEncoding = "UTF-8")
  rv$batch_author_aac_csv_path <- batch_csv_path

  addlog("[BATCH AAC] DONE | CSV saved: ", normalizePath(batch_csv_path, winslash = "/", mustWork = FALSE),
         " | Download button: series author AAC CSV")
  return(NULL)
}



    rv$log <- ""
    rv$nodes_full <- rv$edges_full <- NULL
    rv$author_nodes <- rv$author_edges <- NULL
    rv$country_nodes <- rv$country_edges <- NULL
    rv$stateprov_nodes <- rv$stateprov_edges <- NULL
    rv$inst_nodes <- rv$inst_edges <- NULL
    rv$dept_nodes <- rv$dept_edges <- NULL
    rv$mesh_nodes <- rv$mesh_edges <- NULL

    # cached lists/edges for on-demand domain runs
    rv$countries_list <- rv$stateprov_list <- NULL
    rv$inst_list <- rv$dept_list <- rv$mesh_list <- NULL
    rv$taaa_df <- rv$taaa_freq_df <- NULL
    rv$taaa_profile_map_df <- rv$taaa_kappa_df <- rv$taaa_conf_df <- NULL
    

# =========================
# CO-OCCURRENCE EDGES (PAIRWISE) - ROBUST
# =========================
# term_list: list(article -> character vector or delimited string)
# returns undirected edges data.frame(from,to,weight)
normalize_terms <- function(x){
  if (is.null(x)) return(character())
  x <- as.character(x)
  x <- x[nzchar(trimws(x))]
  if (!length(x)) return(character())
  # If a single string contains delimiters, split it
  if (length(x) == 1L) {
    s <- x[[1]]
    if (grepl("[;|]", s, fixed=FALSE) || grepl("\\s{2,}", s)) {
      x <- unlist(strsplit(s, "[;|]", perl=TRUE), use.names = FALSE)
    } else if (grepl(",\\s*", s)) {
      # be conservative: only split by comma if it looks like a list (has many commas)
      if (length(gregexpr(",", s, fixed=TRUE)[[1]]) >= 3) {
        x <- unlist(strsplit(s, ",", fixed=TRUE), use.names = FALSE)
      }
    }
  }
  x <- trimws(x)
  x <- x[nzchar(x)]
  unique(x)
}

build_cooc_edges <- function(term_list){
  if (is.null(term_list) || !length(term_list)) {
    return(data.frame(from=character(), to=character(), weight=numeric(), stringsAsFactors = FALSE))
  }
  all_pairs <- list()
  k <- 0L
  n_ge2 <- 0L
  sizes <- integer(length(term_list))
  for (i in seq_along(term_list)) {
    terms <- normalize_terms(term_list[[i]])
    sizes[i] <- length(terms)
    if (length(terms) < 2L) next
    n_ge2 <- n_ge2 + 1L
    pairs <- t(combn(terms, 2))
    k <- k + 1L
    all_pairs[[k]] <- data.frame(from=pairs[,1], to=pairs[,2], weight=1, stringsAsFactors = FALSE)
  }
  if (!k) {
    return(data.frame(from=character(), to=character(), weight=numeric(), stringsAsFactors = FALSE))
  }
  ed <- do.call(rbind, all_pairs)
  # undirected: sort endpoints
  a <- pmin(ed$from, ed$to)
  b <- pmax(ed$from, ed$to)
  ed$from <- a; ed$to <- b
  agg <- aggregate(weight ~ from + to, data = ed, FUN = sum)
  agg$weight <- as.numeric(agg$weight)
  agg <- agg[order(-agg$weight), , drop=FALSE]
  attr(agg, "n_articles") <- length(term_list)
  attr(agg, "n_ge2") <- n_ge2
  attr(agg, "max_items") <- ifelse(length(sizes), max(sizes, na.rm=TRUE), 0)
  attr(agg, "median_items") <- ifelse(length(sizes), median(sizes, na.rm=TRUE), 0)
  agg
}

rv$edges_country <- rv$edges_stateprov <- NULL
    rv$edges_inst <- rv$edges_dept <- rv$edges_mesh <- NULL

    rv$report_path <- NULL
    term <- term_rx()
    retmax <- input$retmax
    seed <- input$seed

    withProgress(message = "Fetching PubMed → parsing metadata", value = 0, {
      # Decide data source: permanent uploaded MEDLINE txt (preferred) or online PubMed query
      if (isTRUE(input$use_uploaded_txt) && file.exists(PERM_PUBMED_TXT)) {
        incProgress(0.05, detail = "Reading uploaded MEDLINE…")
        addlog("[RUN] Using uploaded MEDLINE: ", basename(PERM_PUBMED_TXT))
        fi <- file.info(PERM_PUBMED_TXT)
        addlog("[RUN] Upload file size: ", fi$size, " bytes | mtime: ", format(fi$mtime))
        pmids <- parse_pmids_from_medline_txt(PERM_PUBMED_TXT)
        if (length(pmids) == 0) { addlog("[STOP] No PMID- lines found in uploaded MEDLINE file."); return(NULL) }
        rv$pmids <- pmids
        rv$pmids_from_txt <- TRUE
        # mimic entrez_search count for downstream logs
        s <- list(count = length(pmids))
      } else {
        incProgress(0.05, detail = "Searching PubMed…")
        addlog("[RUN] term = ", term)
        addlog("[RUN] retmax = ", retmax)

        s <- tryCatch(entrez_search(db="pubmed", term=term, retmax=retmax), error=function(e) e)
        if (inherits(s, "error")) { addlog("[ERR] PubMed search failed: ", conditionMessage(s)); return(NULL) }

        pmids <- s$ids
        rv$pmids <- pmids
        rv$pmids_from_txt <- FALSE
      }
      # iCite supports direct PMID list via URL: https://icite.od.nih.gov/analysis?pmids=PMID1,PMID2,... (up to ~900 is convenient for URLs)
      if (length(pmids) > 0) {
        pmids_for_url <- pmids[seq_len(min(length(pmids), 900))]
        rv$icite_url <- paste0("https://icite.od.nih.gov/analysis?pmids=", paste(pmids_for_url, collapse=",")) 
      } else {
        rv$icite_url <- NULL
      }

      addlog("[IO] PubMed hits: ", s$count, " | fetched pmids: ", length(pmids))
      if (length(pmids) == 0) { addlog("[STOP] No PMIDs fetched. Try a different query."); return(NULL) }

      # Save PMID list from the current query/upload for direct download and WebZIP.
      try({
        dl_dir <- .get_dl_dir()
        dir.create(dl_dir, recursive = TRUE, showWarnings = FALSE)
        pmid_df <- data.frame(PMID = as.character(pmids), stringsAsFactors = FALSE)
        utils::write.csv(pmid_df, file.path(dl_dir, "pmids_query_result.csv"), row.names = FALSE, fileEncoding = "UTF-8")
        writeLines(as.character(pmids), file.path(dl_dir, "pmids_query_result.txt"), useBytes = TRUE)
        addlog("[IO] PMID list saved: download/pmids_query_result.csv and download/pmids_query_result.txt")
      }, silent = TRUE)

      
incProgress(0.20, detail = "Fetching XML (batched)…")

# ---- Robust PubMed fetch: use history + batch + retry + backoff ----
fetch_pubmed_history_xml <- function(web_history, total_count, batch_size = 200,
                                     max_retries = 4, sleep_base = 0.35) {
  stopifnot(!is.null(web_history))
  # If API key is set, we can be a bit faster
  has_key <- nzchar(Sys.getenv("ENTREZ_KEY", unset=""))
  if (!has_key) sleep_base <- max(sleep_base, 0.6)

  out <- character()
  retstart <- 0
  while (retstart < total_count) {
    bs <- min(batch_size, total_count - retstart)
    ok <- FALSE
    tries <- 0
    cur_bs <- bs

    while (!ok && tries < max_retries) {
      tries <- tries + 1
      x <- tryCatch(
        rentrez::entrez_fetch(db="pubmed", web_history=web_history,
                             rettype="xml", retmode="xml",
                             retstart=retstart, retmax=cur_bs,
                             parsed=FALSE),
        error=function(e) ""
      )
      if (nzchar(x)) {
        out <- c(out, x)
        ok <- TRUE
      } else {
        # backoff + optionally reduce batch size
        Sys.sleep(sleep_base * tries * (if (has_key) 1 else 1.6))
        if (cur_bs > 25) cur_bs <- max(25, floor(cur_bs / 2))
      }
    }

    if (!ok) {
      return(list(ok=FALSE, xml=""))
    }

    retstart <- retstart + bs
    Sys.sleep(sleep_base)
  }
  list(ok=TRUE, xml=paste(out, collapse="\n"))
}

# Prefer history-based fetch (more stable on shinyapps.io)
# IMPORTANT: In uploaded-MEDLINE mode, we must NOT run any term/history search.
if (isTRUE(rv$pmids_from_txt)) {
  addlog("[RUN] Upload mode: fetching XML by uploaded PMIDs (no query/history).")
  xml_parts <- character()
  idxs <- seq(1, length(pmids), by=200)
  for (i in idxs) {
    chunk <- pmids[i:min(i+199, length(pmids))]
    x <- tryCatch(rentrez::entrez_fetch(db="pubmed", id=chunk, rettype="xml", parsed=FALSE), error=function(e) "")
    if (!nzchar(x)) { addlog("[ERR] XML chunk fetch failed or empty (PMID batch)."); return(NULL) }
    xml_parts <- c(xml_parts, x)
    Sys.sleep(0.4)
  }
  xml_txt <- paste(xml_parts, collapse="\n")
} else {
  s_hist <- tryCatch(rentrez::entrez_search(db="pubmed", term=term, retmax=0, use_history=TRUE), error=function(e) NULL)
  if (!is.null(s_hist) && !is.null(s_hist$web_history) && isTRUE(s$count > 0)) {
    fb <- fetch_pubmed_history_xml(s_hist$web_history, total_count = length(pmids), batch_size = 200)
    xml_txt <- fb$xml
    if (!isTRUE(fb$ok) || !nzchar(xml_txt)) { addlog("[ERR] XML fetch failed or empty (history batch)."); return(NULL) }
  } else {
    addlog("[WARN] History not available; falling back to PMID batching.")
    xml_parts <- character()
    idxs <- seq(1, length(pmids), by=200)
    for (i in idxs) {
      chunk <- pmids[i:min(i+199, length(pmids))]
      x <- tryCatch(rentrez::entrez_fetch(db="pubmed", id=chunk, rettype="xml", parsed=FALSE), error=function(e) "")
      if (!nzchar(x)) { addlog("[ERR] XML chunk fetch failed or empty."); return(NULL) }
      xml_parts <- c(xml_parts, x)
      Sys.sleep(0.4)
    }
    xml_txt <- paste(xml_parts, collapse="\n")
  }
}


      if (isTRUE(input$save_xml)) {
        writeLines(xml_txt, con="pubmed.xml", useBytes=TRUE)
        addlog("[OK] Wrote pubmed.xml")
      }

      incProgress(0.30, detail = "Fetching MEDLINE (optional)…")
      if (isTRUE(input$save_xml)) {
        med_txt <- tryCatch(entrez_fetch(db="pubmed", id=pmids, rettype="medline", retmode="text", parsed=FALSE), error=function(e) e)
        if (!inherits(med_txt, "error") && nzchar(med_txt)) {
          writeLines(med_txt, con="pubmed_medline.txt", useBytes=TRUE)
          addlog("[OK] Wrote pubmed_medline.txt")
        } else addlog("[WARN] MEDLINE fetch empty/unavailable (non-fatal).")
      }

      incProgress(0.45, detail = "Splitting articles…")
      articles <- regmatches(xml_txt, gregexpr("<PubmedArticle\\b[\\s\\S]*?</PubmedArticle>", xml_txt, perl=TRUE))[[1]]
      addlog("[IO] Parsed articles: ", length(articles))


# --- FA/LA affiliation helpers (stable counts for Summary) ---
.extract_author_blocks <- function(article_xml){
  # returns list of <Author>...</Author> blocks
  blks <- regmatches(article_xml, gregexpr("<Author\\b[\\s\\S]*?</Author>", article_xml, perl=TRUE))[[1]]
  if (!length(blks)) return(character())
  # prefer personal authors with LastName
  personal <- blks[grepl("<LastName>", blks, fixed=TRUE)]
  if (length(personal)) personal else blks
}

.extract_affils_from_author <- function(author_block){
  aff <- regmatches(author_block, gregexpr("<Affiliation[^>]*>([\\s\\S]*?)</Affiliation>", author_block, perl=TRUE))[[1]]
  if (!length(aff)) return(character())
  # strip tags inside affiliation if any
  aff <- gsub("^<Affiliation[^>]*>|</Affiliation>$", "", aff, perl=TRUE)
  aff <- gsub("<[^>]+>", " ", aff)
  aff <- gsub("[\\r\\n\\t]+", " ", aff)
  aff <- trimws(aff)
  aff <- aff[nzchar(aff)]
  aff
}

.biblio_from_affil_text <- function(aff_txt){
  # Use existing parse_article_biblio() mapping by wrapping as minimal PubmedArticle with Affiliation tags.
  # This keeps country/state/institute/department logic consistent with the rest of the app.
  if (is.null(aff_txt) || !nzchar(aff_txt)) {
    return(list(countries=character(), stateprov=character(), institutes=character(), departments=character(), mesh=character()))
  }
  mini <- paste0("<PubmedArticle><Affiliation>", aff_txt, "</Affiliation></PubmedArticle>")
  out <- try(parse_article_biblio(mini), silent=TRUE)
  if (inherits(out, "try-error") || is.null(out)) {
    list(countries=character(), stateprov=character(), institutes=character(), departments=character(), mesh=character())
  } else out
}

.first_term <- function(x){
  if (is.null(x) || !length(x)) return("")
  x <- trimws(as.character(x))
  x <- x[nzchar(x)]
  if (!length(x)) "" else x[1]
}

.count_fala_terms <- function(fa_vec, la_vec){
  # per-paper unique count: if FA==LA count once; else count both once
  fa <- trimws(as.character(fa_vec)); la <- trimws(as.character(la_vec))
  fa[is.na(fa)] <- ""; la[is.na(la)] <- ""
  n <- max(length(fa), length(la))
  length(fa) <- n; length(la) <- n
  keys <- vector("character", 0)
  for (i in seq_len(n)){
    a <- fa[i]; b <- la[i]
    if (!nzchar(a) && !nzchar(b)) next
    if (nzchar(a) && nzchar(b) && identical(a,b)) {
      keys <- c(keys, a)
    } else {
      if (nzchar(a)) keys <- c(keys, a)
      if (nzchar(b)) keys <- c(keys, b)
    }
  }
  tb <- sort(table(keys), decreasing=TRUE)
  data.frame(name = names(tb), value = as.numeric(tb), stringsAsFactors = FALSE)
}

# --- PMID table (Most cited tab): PMID + Title + Journal + Year ---
get_tag1 <- function(x, pat){
  m <- regexec(pat, x, perl=TRUE)
  r <- regmatches(x, m)[[1]]
  if (length(r) >= 2) r[2] else ""
}
extract_pmid <- function(article){
  pm <- get_tag1(article, "<PMID[^>]*>([0-9]+)</PMID>")
  trimws(pm)
}
extract_title <- function(article){
  tt <- get_tag1(article, "<ArticleTitle>([\\s\\S]*?)</ArticleTitle>")
  tt <- gsub("<[^>]+>", "", tt)  # strip nested tags
  tt <- gsub("[\r\n\t]+", " ", tt)
  trimws(tt)
}
extract_abstract <- function(article){
  aa <- regmatches(article, gregexpr("<AbstractText[^>]*>([\\s\\S]*?)</AbstractText>", article, perl = TRUE))[[1]]
  if (!length(aa)) return("")
  aa <- gsub("<[^>]+>", " ", aa)
  aa <- gsub("[\r\n\t]+", " ", aa)
  aa <- gsub("\\s+", " ", aa)
  trimws(paste(aa, collapse = " "))
}



      # ----- PubMeta: Journal / Year / Articles -----
      get_first <- function(x, pat){
        m <- regexec(pat, x, perl=TRUE)
        r <- regmatches(x, m)[[1]]
        if (length(r) >= 2) r[2] else ""
      }
      extract_year <- function(article){
        y <- get_first(article, "<PubDate>[\\s\\S]*?<Year>([0-9]{4})</Year>")
        if (nzchar(y)) return(y)
        y <- get_first(article, "<ArticleDate[^>]*>[\\s\\S]*?<Year>([0-9]{4})</Year>")
        if (nzchar(y)) return(y)
        md <- get_first(article, "<MedlineDate>([0-9]{4})[^<]*</MedlineDate>")
        if (nzchar(md)) return(md)
        ""
      }
      extract_journal <- function(article){
        j <- get_first(article, "<ISOAbbreviation>([\\s\\S]*?)</ISOAbbreviation>")
        if (nzchar(j)) return(gsub("[\\r\\n\\t]+", " ", j))
        j <- get_first(article, "<Title>([\\s\\S]*?)</Title>")
        j <- gsub("[\\r\\n\\t]+", " ", j)
        trimws(j)
      }

# Build pmid table for Most cited tab
pmid_vec  <- vapply(articles, extract_pmid, character(1))
title_vec <- vapply(articles, extract_title, character(1))
abstract_vec <- vapply(articles, extract_abstract, character(1))
year_vec  <- vapply(articles, extract_year, character(1))
jour_vec  <- vapply(articles, extract_journal, character(1))
keep <- nzchar(pmid_vec)
rv$pmid_tbl <- data.frame(
  PMID = pmid_vec[keep],
  Year = year_vec[keep],
  Journal = jour_vec[keep],
  Title = title_vec[keep],
  Abstract = abstract_vec[keep],
  stringsAsFactors = FALSE
)

# Store per-article Year vector aligned with articles (for slopegraphs)
rv$article_years <- year_vec

      rv$pubmeta <- data.frame(
        PMID    = pmid_vec,
        Title   = title_vec,
        Abstract = abstract_vec,
        Journal = vapply(articles, extract_journal, character(1)),
        Year    = vapply(articles, extract_year, character(1)),
        stringsAsFactors = FALSE
      )
      rv$pubmeta <- rv$pubmeta[nzchar(rv$pubmeta$Year) & nzchar(rv$pubmeta$Journal), , drop=FALSE]
      addlog("[IO] pubmeta rows: ", nrow(rv$pubmeta))
# ---- Journal/Year network objects for Report tab (Top20 via FLCA) ----
rv$jy_nodes <- rv$jy_edges <- NULL
rv$jy_max_items <- rv$jy_median_items <- rv$jy_n_ge2 <- NA
try({
  if (is.data.frame(rv$pubmeta) && nrow(rv$pubmeta) > 0) {
    jy_list <- lapply(seq_len(nrow(rv$pubmeta)), function(i) c(rv$pubmeta$Journal[i], rv$pubmeta$Year[i]))
    rv$jy_list <- jy_list
    rv$edges_jy <- build_cooc_edges(jy_list)
    ts_jy <- .term_stats(jy_list)
    rv$jy_max_items <- ts_jy$max_items
    rv$jy_median_items <- ts_jy$median_items
    rv$jy_n_ge2 <- ts_jy$n_ge2

    # build Top20 + edges (keeps bipartite structure if present)
    jy_out <- tryCatch(build_domain_flca_top20(rv$edges_jy, term_list=jy_list, seed=seed), error=function(e) NULL)
    if (!is.null(jy_out)) {
      rv$jy_nodes <- jy_out$nodes
      rv$jy_edges <- jy_out$data
    }

    # AAC metrics for Journal/Year (computed when JY network exists)
    .update_domain_aac(rv, 'jy', rv$jy_nodes, rv$jy_edges)
    addlog("[JY] nodes=", if(is.data.frame(rv$jy_nodes)) nrow(rv$jy_nodes) else 0,
           " edges=", if(is.data.frame(rv$jy_edges)) nrow(rv$jy_edges) else 0,
           " max_items=", rv$jy_max_items)
  }
}, silent=TRUE)




      if (length(articles) == 0) { addlog("[STOP] No <PubmedArticle> blocks found (unexpected)."); return(NULL) }

      # ----- Author network -----
      incProgress(0.55, detail = "Building author edges/nodes…")
      author_lists <- lapply(articles, get_authors_per_article)
      author_lists_all <- lapply(articles, get_authors_per_article)
      rv$author_lists_all <- author_lists_all  # keep singles for year-frequency slope
      # 1st+Last author list per article (for Term/Year)
      rv$author_lists_fl <- lapply(author_lists_all, function(a){
        a <- unique(trimws(as.character(a)))
        a <- a[nzchar(a)]
        if (length(a) == 0) return(character(0))
        if (length(a) == 1) return(a)
        c(a[1], a[length(a)])
      })


      # n(PubMed) and byline counts must use ALL authors in the byline.
      # First/last authors are used only for value(first/last) and author edges.
      # Example: Chien TW can have n(any byline)=111 but n(first/last)=18.
      author_lists_for_counts <- author_lists_all
      n_pubmed_tab <- table(unlist(lapply(author_lists_for_counts, function(a){
        a <- as.character(a); a <- trimws(a); a <- a[nzchar(a)]
        unique(a)  # one PubMed record counts an author once
      })))

      n_byline_tab <- table(unlist(lapply(author_lists_for_counts, function(a){
        a <- as.character(a); a <- trimws(a); a <- a[nzchar(a)]
        unique(a)  # any byline appearance by record, not just first/last
      })))
      # single-author articles (first==last, length==1)
      single_tab <- table(unlist(lapply(author_lists_all, function(a){
        a <- as.character(a); a <- trimws(a); a <- a[nzchar(a)]
        if (length(a)==1) return(a[1]) else character(0)
      })))
      author_lists <- author_lists[vapply(author_lists, length, integer(1)) >= 2]
      rv$author_lists <- author_lists
      addlog("[IO] Articles with >=2 authors: ", length(author_lists))
      if (length(author_lists) == 0) { addlog("[STOP] All fetched articles have <2 authors; cannot form first-last edges."); return(NULL) }

      edges_raw <- do.call(rbind, lapply(author_lists, function(a) data.frame(Leader=a[1], Follower=a[length(a)], stringsAsFactors=FALSE)))
      edges_raw <- edges_raw[edges_raw$Leader != edges_raw$Follower, , drop=FALSE]
      addlog("[EDGE] Raw edges rows: ", nrow(edges_raw))
      if (nrow(edges_raw) == 0) { addlog("[STOP] Edge table is empty after self-loop removal."); return(NULL) }

      key  <- paste(edges_raw$Leader, edges_raw$Follower, sep="|||")
      tab  <- table(key)
      keys <- names(tab)
      data_edges <- data.frame(
        Leader   = sub("\\|\\|\\|.*$", "", keys),
        Follower = sub("^.*\\|\\|\\|", "", keys),
        WCD      = as.integer(tab),
        stringsAsFactors = FALSE
      )
      data_edges <- data_edges[order(-data_edges$WCD), , drop=FALSE]
      rv$author_edges_full <- data_edges  # full first-last edge table (global counts)

      addlog("[EDGE] Aggregated edges rows: ", nrow(data_edges))

      node_names <- sort(unique(c(data_edges$Leader, data_edges$Follower)))
      out_w <- tapply(data_edges$WCD, data_edges$Leader,   sum)
      in_w  <- tapply(data_edges$WCD, data_edges$Follower, sum)
      out_val <- out_w[match(node_names, names(out_w))]; out_val[is.na(out_val)] <- 0
      in_val  <- in_w[match(node_names, names(in_w))];   in_val[is.na(in_val)] <- 0

      nodes <- data.frame(name=node_names, stringsAsFactors=FALSE)

      # --- Author metrics (FA/LA + single) ---
      # FA (first-author coword strength): sum(WCD) where node is Leader
      # LA (last-author coword strength):  sum(WCD) where node is Follower
      str_out <- tapply(as.numeric(data_edges$WCD), data_edges$Leader, sum)
      str_in  <- tapply(as.numeric(data_edges$WCD), data_edges$Follower, sum)

      nodes$FA <- as.integer(str_out[node_names]); nodes$FA[is.na(nodes$FA)] <- 0L
      nodes$LA <- as.integer(str_in[node_names]);  nodes$LA[is.na(nodes$LA)] <- 0L

      # value2 = count(edge) (unique partners on FA/LA edges; self-loop excluded)
      out_deg <- tapply(data_edges$Follower, data_edges$Leader, function(v) length(unique(v)))
      in_deg  <- tapply(data_edges$Leader,   data_edges$Follower, function(v) length(unique(v)))
      nodes$out_deg <- as.integer(out_deg[node_names]); nodes$out_deg[is.na(nodes$out_deg)] <- 0L
      nodes$in_deg  <- as.integer(in_deg[node_names]);  nodes$in_deg[is.na(nodes$in_deg)]  <- 0L
      nodes$value2  <- as.integer(nodes$out_deg + nodes$in_deg)
      
      # value_strength = FA + LA (coword strength; self not included)
      nodes$value_strength <- as.integer(nodes$FA + nodes$LA)

      # value2_strength = AAC (computed later from top3 of value2 within the chosen domain)
      nodes$value2_strength <- NA_real_
# n(PubMed): byline appearances (each article counts author once)
      nodes$n_pubmed <- as.integer(n_pubmed_tab[node_names]); nodes$n_pubmed[is.na(nodes$n_pubmed)] <- 0L

      nodes$n_byline <- as.integer(n_byline_tab[node_names]); nodes$n_byline[is.na(nodes$n_byline)] <- 0L

      # single-author articles count
      nodes$single_author <- as.integer(single_tab[node_names]); nodes$single_author[is.na(nodes$single_author)] <- 0L

      nodes$self_coword <- nodes$single_author

      # value: (FA+LA) + single_author (self coword)
      nodes$value <- as.integer(nodes$value_strength + nodes$single_author)
      nodes$value[is.na(nodes$value)] <- 0L

      # AAC from top-3 value2 (domain-level, attached for convenience)
      nodes$value2_strength <- (.AAC_INLINE(nodes$value2))

      nodes <- nodes[order(-nodes$value, nodes$name), , drop=FALSE]

      rv$nodes_full <- nodes
      rv$edges_full <- data_edges
      write.csv(nodes, "nodes.csv", row.names=FALSE, fileEncoding="UTF-8")
      write.csv(data_edges, "data_edges.csv", row.names=FALSE, fileEncoding="UTF-8")
      addlog("[OK] Wrote nodes.csv and data_edges.csv")

      # normalize node/value, edge schema (Leader/Follower/WCD)
      nodes <- normalize_nodes(nodes)
      data_edges <- fix_edge_cols(normalize_edges(data_edges))

      # Use the internal wrapper: FLCA + MajorSampling Top20 + silhouette/a*.
      # IMPORTANT: do not call run_flca_ms_sil() directly here; the module runner
      # can expect a list(input) and may fail with "not all is.data.frame(nodes)".
      cfg <- list(top_clusters=5, base_per_cluster=4, target_n=20, intra_delta=2, inter_delta=5, eps=1e-9)
      res <- tryCatch({
        # Prefer the uploaded flca_ms_sil_module.R runner. It requires two
        # data.frames (nodes, edges0), not a list. It returns modes + data.
        if (exists("run_flca_ms_sil_runner", mode = "function")) {
          z <- run_flca_ms_sil_runner(nodes, data_edges, cfg = cfg, verbose = TRUE)
          if (!is.null(z$modes) && is.data.frame(z$modes)) {
            list(nodes = z$modes, data = z$data, raw = z, engine = "run_flca_ms_sil_runner")
          } else if (!is.null(z$nodes) && is.data.frame(z$nodes)) {
            list(nodes = z$nodes, data = z$data, raw = z, engine = "run_flca_ms_sil_runner")
          } else {
            stop("run_flca_ms_sil_runner returned no modes/nodes data.frame")
          }
        } else {
          run_flca_ms_sil_internal(nodes, data_edges, cfg = cfg, verbose = TRUE)
        }
      }, error = function(e) e)
      if (inherits(res, 'error') || is.null(res) || !is.data.frame(res$nodes) || !nrow(res$nodes)) {
        .res_msg <- if (inherits(res, "error")) conditionMessage(res) else "no valid Top20 result"
        warning("[WARN] Author FLCA-MA-SIL Top20 failed; using deterministic Top20 fallback: ", .res_msg, call. = FALSE)
        nodes20 <- nodes[order(-suppressWarnings(as.numeric(nodes$value)), as.character(nodes$name)), , drop = FALSE]
        nodes20 <- utils::head(nodes20, 20)
        if (!('carac' %in% names(nodes20))) nodes20$carac <- 1L
        if (!('ssi' %in% names(nodes20))) nodes20$ssi <- 0
        if (!('a_star1' %in% names(nodes20))) nodes20$a_star1 <- 0
        if (!('a_i' %in% names(nodes20))) nodes20$a_i <- 0
        if (!('b_i' %in% names(nodes20))) nodes20$b_i <- 0
        top_names <- as.character(nodes20$name)
        edges20 <- data_edges[data_edges$Leader %in% top_names & data_edges$Follower %in% top_names, , drop = FALSE]
      } else {
        nodes20 <- res$nodes
        edges20 <- res$data
        nodes20 <- nodes20[order(-suppressWarnings(as.numeric(nodes20$value)), as.character(nodes20$name)), , drop = FALSE]
        nodes20 <- utils::head(nodes20, 20)
        top_names <- as.character(nodes20$name)
        if (is.data.frame(edges20) && nrow(edges20)) {
          if ("follower" %in% names(edges20) && !("Follower" %in% names(edges20))) names(edges20)[names(edges20) == "follower"] <- "Follower"
          edges20 <- edges20[edges20$Leader %in% top_names & edges20$Follower %in% top_names, , drop = FALSE]
        } else {
          edges20 <- data.frame(Leader=character(), Follower=character(), WCD=numeric(), stringsAsFactors = FALSE)
        }
      }
      # Enforce FLCA one-link structure before SS: each follower has only one leader.
      edges20 <- .author_enforce_single_edge_per_follower(nodes20, edges20)
      addlog("[AUTHOR one-link] followers=", length(unique(edges20$Follower)),
             " edges=", nrow(edges20),
             "; duplicated followers=", sum(duplicated(edges20$Follower)))

      # Make sure the displayed Top20 has at least two cluster labels for SS.
      fixed_author <- .force_author_min2_clusters_and_ss(nodes20, edges20, data_edges = data_edges, cfg = cfg)
      nodes20 <- fixed_author$nodes
      edges20 <- .author_enforce_single_edge_per_follower(nodes20, fixed_author$edges)
      addlog("[AUTHOR Top20] nodes=", if(is.data.frame(nodes20)) nrow(nodes20) else 0,
             " clusters=", if(is.data.frame(nodes20) && "carac" %in% names(nodes20)) length(unique(nodes20$carac)) else 0,
             " edges=", if(is.data.frame(edges20)) nrow(edges20) else 0,
             if (isTRUE(fixed_author$changed)) paste0("; cluster repair: ", fixed_author$reason) else "")
      # keep in rv for downstream tables/plots
      rv$author_nodes <- nodes20
      rv$author_edges <- edges20

      

      # Augment Author metrics using the FULL byline and FULL first/last edges.
      # n_pubmed / n_byline = any byline appearance by PubMed record.
      # value = first/last author appearance including single-author papers once.
      author_lists_metric <- rv$author_lists_all %||% author_lists_all
      if (!is.null(author_lists_metric) && length(author_lists_metric) > 0) {
        alist_any <- lapply(author_lists_metric, function(a){
          a <- unique(trimws(as.character(a)))
          a <- a[nzchar(a)]
          a
        })
        flat_any <- unlist(alist_any, use.names = FALSE)
        if (length(flat_any) > 0) {
          tab_any <- sort(table(flat_any), decreasing = TRUE)

          # single-author papers: length==1; count once in value(first/last incl. single)
          singles <- unlist(lapply(alist_any, function(a) if (length(a)==1) a else character(0)), use.names = FALSE)
          tab_s <- if (length(singles) > 0) table(singles) else table(character(0))

          # Full first/last non-single edge strengths from the full author edge table.
          ed_full <- rv$author_edges_full %||% data_edges
          if (!is.null(ed_full) && is.data.frame(ed_full) && nrow(ed_full) > 0) {
            if (!("Leader" %in% names(ed_full)) && "leader" %in% names(ed_full)) ed_full$Leader <- ed_full$leader
            if (!("Follower" %in% names(ed_full)) && "follower" %in% names(ed_full)) ed_full$Follower <- ed_full$follower
            if (!("WCD" %in% names(ed_full))) ed_full$WCD <- 1
            ed_full$Leader <- trimws(as.character(ed_full$Leader))
            ed_full$Follower <- trimws(as.character(ed_full$Follower))
            ed_full$WCD <- suppressWarnings(as.numeric(ed_full$WCD))
            ed_full <- ed_full[nzchar(ed_full$Leader) & nzchar(ed_full$Follower) & is.finite(ed_full$WCD) & ed_full$WCD > 0, , drop=FALSE]
          } else {
            ed_full <- data.frame(Leader=character(), Follower=character(), WCD=numeric(), stringsAsFactors=FALSE)
          }
          fa_full <- if (nrow(ed_full)) tapply(ed_full$WCD, ed_full$Leader, sum, na.rm=TRUE) else numeric(0)
          la_full <- if (nrow(ed_full)) tapply(ed_full$WCD, ed_full$Follower, sum, na.rm=TRUE) else numeric(0)

          all_names_metric <- sort(unique(c(names(tab_any), names(tab_s), names(fa_full), names(la_full))))
          df_pub <- data.frame(name = all_names_metric, stringsAsFactors = FALSE)
          df_pub$n_pubmed <- as.integer(tab_any[df_pub$name]); df_pub$n_pubmed[is.na(df_pub$n_pubmed)] <- 0L
          df_pub$n_byline <- df_pub$n_pubmed
          df_pub$single_author <- as.integer(tab_s[df_pub$name]); df_pub$single_author[is.na(df_pub$single_author)] <- 0L
          df_pub$FA <- as.numeric(fa_full[df_pub$name]); df_pub$FA[!is.finite(df_pub$FA)] <- 0
          df_pub$LA <- as.numeric(la_full[df_pub$name]); df_pub$LA[!is.finite(df_pub$LA)] <- 0
          df_pub$value2 <- as.numeric(df_pub$FA + df_pub$LA)
          df_pub$value <- as.numeric(df_pub$value2 + df_pub$single_author)
          df_pub$value_strength <- df_pub$value2
          df_pub$self_coword <- df_pub$single_author

          # merge into sampled nodes and replace stale Top20/one-link metrics.
          if (!("name" %in% names(rv$author_nodes))) {
            if ("label" %in% names(rv$author_nodes)) rv$author_nodes$name <- as.character(rv$author_nodes$label)
            else if ("id" %in% names(rv$author_nodes)) rv$author_nodes$name <- as.character(rv$author_nodes$id)
          }
          rv$author_nodes$name <- as.character(rv$author_nodes$name)
          for (.cc in c("n_pubmed","n_byline","single_author","FA","LA","value2","value","value_strength","self_coword")) {
            if (.cc %in% names(rv$author_nodes)) rv$author_nodes[[.cc]] <- NULL
          }
          rv$author_nodes <- merge(rv$author_nodes, df_pub, by="name", all.x=TRUE)
          for (.cc in c("n_pubmed","n_byline","single_author","FA","LA","value2","value","value_strength","self_coword")) {
            if (.cc %in% names(rv$author_nodes)) {
              rv$author_nodes[[.cc]] <- suppressWarnings(as.numeric(rv$author_nodes[[.cc]]))
              rv$author_nodes[[.cc]][!is.finite(rv$author_nodes[[.cc]])] <- 0
            }
          }
        }
      }

      # Recompute displayed Author metrics from the same one-link Top20 edges.
      # Full first/last edges are kept separately as rv$author_edges_full for downloads only.
      rv$author_edges <- .author_enforce_single_edge_per_follower(rv$author_nodes, rv$author_edges)
      rv$author_nodes <- .augment_author_nodes(rv$author_nodes, rv$author_edges)
# debug top-3 (console)
      cat('\n[CHECK] Top-3 author nodes20 (name,value,value2,carac,ssi,a_star1):\n')
      print(utils::head(rv$author_nodes[, intersect(c('name','value','value2','carac','ssi','a_star1'), names(rv$author_nodes)), drop=FALSE], 3))
      # ----- AAC metrics for Report (author) -----
      .update_domain_aac(rv, 'author', rv$author_nodes, rv$author_edges)

      # ----- Sankeymatic code + URL -----
      rv$sankey_code <- ""
      rv$sankey_url  <- ""
      try({
        ed_tmp <- rv$author_edges
        if (is.data.frame(ed_tmp) && nrow(ed_tmp)>0) {
          if (!all(c("Leader","Follower") %in% names(ed_tmp))) {
            if (all(c("from","to") %in% names(ed_tmp))) names(ed_tmp)[match(c("from","to"), names(ed_tmp))] <- c("Leader","Follower")
          }
          if (!("WCD" %in% names(ed_tmp))) ed_tmp$WCD <- if ("value" %in% names(ed_tmp)) ed_tmp$value else 1
          sank_lines <- paste0(ed_tmp$Leader, " [", ed_tmp$WCD, "] ", ed_tmp$Follower)
          rv$sankey_code <- paste(sank_lines, collapse = "\n")
          # Sankeymatic expects encoded text after '#'
          rv$sankey_url <- paste0("https://sankeymatic.com/build/?i=", URLencode(rv$sankey_code, reserved=TRUE))
        }
      }, silent=TRUE)


      # ----- Parse biblio -----
      incProgress(0.78, detail = "Parsing Country/Institute/MeSH…")
      biblio <- lapply(articles, parse_article_biblio)

      countries_list <- lapply(biblio, `[[`, "countries")     # already dictionary filtered (no cities)
      stateprov_list <- lapply(biblio, `[[`, "stateprov")  # US states / China provinces / else country
      inst_list      <- lapply(biblio, `[[`, "institutes")
      dept_list      <- lapply(biblio, `[[`, "departments")
      mesh_list_raw  <- lapply(biblio, `[[`, "mesh")
      mesh_list      <- lapply(mesh_list_raw, filter_mesh_terms)

      # store term lists for report diagnostics (max_items/median_items/n_ge2)
      rv$countries_list <- countries_list
      rv$stateprov_list <- stateprov_list
      rv$inst_list      <- inst_list
      rv$dept_list      <- dept_list
      rv$mesh_list      <- mesh_list

      # Keep per-article MeSH headings in rv$pubmeta for MeSH/year slopegraphs.
      # The MeSH string is semicolon-separated so headings with commas remain intact.
      try({
        mesh_chr <- vapply(mesh_list, function(z){
          z <- unique(trimws(as.character(z)))
          z <- z[!is.na(z) & nzchar(z)]
          paste(z, collapse = "; ")
        }, character(1))
        if (is.data.frame(rv$pubmeta) && nrow(rv$pubmeta) > 0) {
          n_pm <- nrow(rv$pubmeta)
          if (length(mesh_chr) >= n_pm) {
            rv$pubmeta$MeSH <- mesh_chr[seq_len(n_pm)]
          } else {
            rv$pubmeta$MeSH <- c(mesh_chr, rep("", n_pm - length(mesh_chr)))
          }
          addlog("[IO] pubmeta MeSH column added: ", sum(nzchar(rv$pubmeta$MeSH)), " article rows with MeSH")
        }
      }, silent = TRUE)
  # remove demographics/general terms

      
      # ----- FA/LA (First/Last author) stable counts for Summary (NOT Top20) -----
      # We derive a single best-mapped term for FA and LA from their first available affiliation.
      fa_country <- character(length(articles)); la_country <- character(length(articles))
      fa_stateprov <- character(length(articles)); la_stateprov <- character(length(articles))
      fa_inst <- character(length(articles)); la_inst <- character(length(articles))
      fa_dept <- character(length(articles)); la_dept <- character(length(articles))

      for (ii in seq_along(articles)) {
        blks <- .extract_author_blocks(articles[[ii]])
        if (!length(blks)) next
        fa_blk <- blks[1]
        la_blk <- blks[length(blks)]

        fa_aff <- .extract_affils_from_author(fa_blk)
        la_aff <- .extract_affils_from_author(la_blk)

        fa_b <- .biblio_from_affil_text(if (length(fa_aff)) fa_aff[1] else "")
        la_b <- .biblio_from_affil_text(if (length(la_aff)) la_aff[1] else "")

        fa_country[ii]   <- .first_term(fa_b$countries)
        la_country[ii]   <- .first_term(la_b$countries)
        fa_stateprov[ii] <- .first_term(fa_b$stateprov)
        la_stateprov[ii] <- .first_term(la_b$stateprov)
        fa_inst[ii]      <- .first_term(fa_b$institutes)
        la_inst[ii]      <- .first_term(la_b$institutes)
        fa_dept[ii]      <- .first_term(fa_b$departments)
        la_dept[ii]      <- .first_term(la_b$departments)
      }

      rv$sum_country_fala   <- .count_fala_terms(fa_country, la_country)
      rv$sum_stateprov_fala <- .count_fala_terms(fa_stateprov, la_stateprov)
      rv$sum_inst_fala      <- .count_fala_terms(fa_inst, la_inst)
      rv$sum_dept_fala      <- .count_fala_terms(fa_dept, la_dept)

rv$edges_country   <- build_cooc_edges(countries_list)
      rv$edges_stateprov <- build_cooc_edges(stateprov_list)
      rv$edges_inst      <- build_cooc_edges(inst_list)
      rv$edges_dept      <- build_cooc_edges(dept_list)
      rv$edges_mesh      <- build_cooc_edges(mesh_list)
# ---- Fill Report stats from edge attributes + term lists ----
# term-list stats (items per article)
ts_country <- .term_stats(countries_list)
ts_state  <- .term_stats(stateprov_list)
ts_inst   <- .term_stats(inst_list)
ts_dept   <- .term_stats(dept_list)
ts_mesh   <- .term_stats(mesh_list)

rv$country_max_items   <- ts_country$max_items
rv$country_median_items<- ts_country$median_items
rv$country_n_ge2       <- ts_country$n_ge2

rv$stateprov_max_items   <- ts_state$max_items
rv$stateprov_median_items<- ts_state$median_items
rv$stateprov_n_ge2       <- ts_state$n_ge2

rv$inst_max_items   <- ts_inst$max_items
rv$inst_median_items<- ts_inst$median_items
rv$inst_n_ge2       <- ts_inst$n_ge2

rv$dept_max_items   <- ts_dept$max_items
rv$dept_median_items<- ts_dept$median_items
rv$dept_n_ge2       <- ts_dept$n_ge2

rv$mesh_max_items   <- ts_mesh$max_items
rv$mesh_median_items<- ts_mesh$median_items
rv$mesh_n_ge2       <- ts_mesh$n_ge2

# ---- Author stats (items per article = #authors) ----
if (!is.null(author_lists) && length(author_lists)) {
  lens_auth <- vapply(author_lists, function(x) length(unique(normalize_terms(x))), integer(1))
  rv$author_max_items    <- if (length(lens_auth)) max(lens_auth) else NA
  rv$author_median_items <- if (length(lens_auth)) as.numeric(stats::median(lens_auth)) else NA
  rv$author_n_ge2        <- sum(lens_auth >= 2, na.rm=TRUE)
}


# ---- Diagnostics (why edges=0?) ----
try({
  addlog("[DIAG] Country items per article: max=", attr(rv$edges_country, "max_items"),
         " median=", attr(rv$edges_country, "median_items"),
         " n_ge2=", attr(rv$edges_country, "n_ge2"),
         " edges=", nrow(rv$edges_country))
  addlog("[DIAG] State/Province items per article: max=", attr(rv$edges_stateprov, "max_items"),
         " median=", attr(rv$edges_stateprov, "median_items"),
         " n_ge2=", attr(rv$edges_stateprov, "n_ge2"),
         " edges=", nrow(rv$edges_stateprov))
  addlog("[DIAG] Institute items per article: max=", attr(rv$edges_inst, "max_items"),
         " median=", attr(rv$edges_inst, "median_items"),
         " n_ge2=", attr(rv$edges_inst, "n_ge2"),
         " edges=", nrow(rv$edges_inst))
  addlog("[DIAG] Department items per article: max=", attr(rv$edges_dept, "max_items"),
         " median=", attr(rv$edges_dept, "median_items"),
         " n_ge2=", attr(rv$edges_dept, "n_ge2"),
         " edges=", nrow(rv$edges_dept))
  addlog("[DIAG] MeSH items per article: max=", attr(rv$edges_mesh, "max_items"),
         " median=", attr(rv$edges_mesh, "median_items"),
         " n_ge2=", attr(rv$edges_mesh, "n_ge2"),
         " edges=", nrow(rv$edges_mesh))
}, silent=TRUE)


            # --- On-demand domain runs (faster) ---
      incProgress(0.85, detail = "Ready. Run a single domain button (Country/State/Institute/Department/MeSH)…")
      addlog("[READY] PubMed fetched + parsed. Author network ready. Domain term-lists cached; click a domain button to run FLCA Top20 only for that domain.")
      incProgress(1.0, detail="Ready.")

    })
      # ---- Auto-run 8 domains for Performance PNG (Top5 + AAC) ----
      # Author + Journal/Year are usually built during fetch; others may be on-demand.
      try({
        if (is.null(rv$country_nodes) && !is.null(rv$edges_country) && !is.null(rv$countries_list)) run_one_domain("Country", input$seed)
        if (is.null(rv$stateprov_nodes) && !is.null(rv$edges_stateprov) && !is.null(rv$stateprov_list)) run_one_domain("State/Province", input$seed)
        if (is.null(rv$inst_nodes) && !is.null(rv$edges_inst) && !is.null(rv$inst_list)) run_one_domain("Institute", input$seed)
        if (is.null(rv$dept_nodes) && !is.null(rv$edges_dept) && !is.null(rv$dept_list)) run_one_domain("Department", input$seed)
        if (is.null(rv$mesh_nodes) && !is.null(rv$edges_mesh) && !is.null(rv$mesh_list)) run_one_domain("MeSH", input$seed)
      }, silent = TRUE)

      # Build long-format for performance report panels
      try({
        if (!is.null(rv$mesh_list) && !is.null(rv$mesh_nodes) && is.data.frame(rv$mesh_nodes) && nrow(rv$mesh_nodes) > 0) {
          taaa_res <- .build_taaa_from_term_list(rv$mesh_list, rv$mesh_nodes)
          if (is.list(taaa_res)) {
            rv$taaa_df <- taaa_res$row_table %||% NULL
            rv$taaa_freq_df <- taaa_res$freq_table %||% NULL
            rv$taaa_profile_map_df <- taaa_res$profile_map_table %||% NULL
            rv$taaa_kappa_df <- taaa_res$kappa_table %||% NULL
            rv$taaa_conf_df <- taaa_res$confusion_table %||% NULL
            addlog("[TAAA] rows=", ifelse(is.null(rv$taaa_df), 0, nrow(rv$taaa_df)))
          }
        }
      }, silent = TRUE)

      try({
        rv$perf_long_auto <- .build_perf_long_auto()
      }, silent = TRUE)

  })

  # ---- Run a single domain on demand (faster, shows progress bar) ----
  run_one_domain <- function(domain, seed){
    if (is.null(domain) || !nzchar(domain)) return(invisible(NULL))
    withProgress(message = paste0("Running domain: ", domain), value = 0, {
      incProgress(0.10, detail = "Preparing inputs…")
      if (identical(domain, "Country")) {
        req(rv$edges_country, rv$countries_list)
        incProgress(0.35, detail = "FLCA + major sampling (Top20)…")
        out <- build_domain_flca_top20(rv$edges_country, term_list=rv$countries_list, seed=seed)
        rv$country_nodes <- out$nodes; rv$country_edges <- out$data

        # AAC metrics (compute immediately after run)
        .update_domain_aac(rv, 'country', rv$country_nodes, rv$country_edges)

        # Abbreviate US/UK in Country dashboard labels
        try({
          if (!is.null(rv$country_nodes) && is.data.frame(rv$country_nodes) && "name" %in% names(rv$country_nodes)) {
            rv$country_nodes$name <- ifelse(rv$country_nodes$name == "United States", "US", rv$country_nodes$name)
            rv$country_nodes$name <- ifelse(rv$country_nodes$name == "United Kingdom", "UK", rv$country_nodes$name)
          }
          if (!is.null(rv$country_edges) && is.data.frame(rv$country_edges)) {
            if ("Leader" %in% names(rv$country_edges)) rv$country_edges$Leader <- ifelse(rv$country_edges$Leader == "United States", "US", rv$country_edges$Leader)
            if ("Follower" %in% names(rv$country_edges)) rv$country_edges$Follower <- ifelse(rv$country_edges$Follower == "United States", "US", rv$country_edges$Follower)
            if ("Leader" %in% names(rv$country_edges)) rv$country_edges$Leader <- ifelse(rv$country_edges$Leader == "United Kingdom", "UK", rv$country_edges$Leader)
            if ("Follower" %in% names(rv$country_edges)) rv$country_edges$Follower <- ifelse(rv$country_edges$Follower == "United Kingdom", "UK", rv$country_edges$Follower)
          }
        }, silent=TRUE)

        addlog("[OK] Country Top20: ", ifelse(is.null(rv$country_nodes), 0, nrow(rv$country_nodes)))
      } else if (identical(domain, "State/Province")) {
        req(rv$edges_stateprov, rv$stateprov_list)
        incProgress(0.35, detail = "FLCA + major sampling (Top20)…")
        out <- build_domain_flca_top20(rv$edges_stateprov, term_list=rv$stateprov_list, seed=seed)
        rv$stateprov_nodes <- out$nodes; rv$stateprov_edges <- out$data
        .update_domain_aac(rv, 'stateprov', rv$stateprov_nodes, rv$stateprov_edges)
        addlog("[OK] State/Province Top20: ", ifelse(is.null(rv$stateprov_nodes), 0, nrow(rv$stateprov_nodes)))
      } else if (identical(domain, "Institute")) {
        req(rv$edges_inst, rv$inst_list)
        incProgress(0.35, detail = "FLCA + major sampling (Top20)…")
        out <- build_domain_flca_top20(rv$edges_inst, term_list=rv$inst_list, seed=seed)
        rv$inst_nodes <- out$nodes; rv$inst_edges <- out$data
        .update_domain_aac(rv, 'inst', rv$inst_nodes, rv$inst_edges)
        addlog("[OK] Institute Top20: ", ifelse(is.null(rv$inst_nodes), 0, nrow(rv$inst_nodes)))
      } else if (identical(domain, "Department")) {
        req(rv$edges_dept, rv$dept_list)
        incProgress(0.35, detail = "FLCA + major sampling (Top20)…")
        out <- build_domain_flca_top20(rv$edges_dept, term_list=rv$dept_list, seed=seed)
        rv$dept_nodes <- out$nodes; rv$dept_edges <- out$data
        .update_domain_aac(rv, 'dept', rv$dept_nodes, rv$dept_edges)
        addlog("[OK] Department Top20: ", ifelse(is.null(rv$dept_nodes), 0, nrow(rv$dept_nodes)))
      } else if (identical(domain, "MeSH")) {
        req(rv$edges_mesh, rv$mesh_list)
        incProgress(0.35, detail = "FLCA + major sampling (Top20)…")
        out <- build_domain_flca_top20(rv$edges_mesh, term_list=rv$mesh_list, seed=seed)
        rv$mesh_nodes <- out$nodes; rv$mesh_edges <- out$data
        .update_domain_aac(rv, 'mesh', rv$mesh_nodes, rv$mesh_edges)
        try({
          taaa_res <- .build_taaa_from_term_list(rv$mesh_list, rv$mesh_nodes)
          if (is.list(taaa_res)) {
            rv$taaa_df <- taaa_res$row_table %||% NULL
            rv$taaa_freq_df <- taaa_res$freq_table %||% NULL
            rv$taaa_profile_map_df <- taaa_res$profile_map_table %||% NULL
            rv$taaa_kappa_df <- taaa_res$kappa_table %||% NULL
            rv$taaa_conf_df <- taaa_res$confusion_table %||% NULL
          }
        }, silent = TRUE)
        addlog("[OK] MeSH Top20: ", ifelse(is.null(rv$mesh_nodes), 0, nrow(rv$mesh_nodes)))
      } else {
        addlog("[WARN] Unknown domain: ", domain)
      }
      incProgress(1.0, detail = "Done.")
    })
  }
    # Record this run for IP-trial / CMC tracking (skip for FREE demo runs)
    if (!isTRUE(bypass_record)) {
      try({ .upsert_ip(ip, as.character(cmc_now %||% ""), inc_count = TRUE) }, silent=TRUE)
    }

  observeEvent(input$run_country, {
    seed <- input$seed
    run_one_domain("Country", seed)
    updateTabsetPanel(session, "main_tabs", selected = "Country")
  }, ignoreInit = TRUE)

  observeEvent(input$run_stateprov, {
    seed <- input$seed
    run_one_domain("State/Province", seed)
    updateTabsetPanel(session, "main_tabs", selected = "State/Province")
  }, ignoreInit = TRUE)

  observeEvent(input$run_inst, {
    seed <- input$seed
    run_one_domain("Institute", seed)
    updateTabsetPanel(session, "main_tabs", selected = "Institute")
  }, ignoreInit = TRUE)

  observeEvent(input$run_dept, {
    seed <- input$seed
    run_one_domain("Department", seed)
    updateTabsetPanel(session, "main_tabs", selected = "Department")
  }, ignoreInit = TRUE)

  observeEvent(input$run_mesh, {
    seed <- input$seed
    run_one_domain("MeSH", seed)
    updateTabsetPanel(session, "main_tabs", selected = "MeSH")
  }, ignoreInit = TRUE)



  # ---- visNetwork helper ----
  render_vn <- function(nodes, data, directed = FALSE, label_size = 20, bold = TRUE, n_pubmed_total = NA_real_) {
    # Robust to column-name variations (Leader/Follower/WCD vs leader/follower/wcd)
    if (is.null(nodes) || is.null(data) || !is.data.frame(nodes) || !nrow(nodes)) {
      return(visNetwork(data.frame(id=1, label="No data", stringsAsFactors = FALSE), data.frame()))
    }

    # nodes
    if (!("name" %in% names(nodes))) {
      nm <- intersect(c("Name","ID","id","label"), names(nodes))
      if (length(nm)) nodes$name <- as.character(nodes[[nm[1]]]) else nodes$name <- as.character(seq_len(nrow(nodes)))
    }
    if (!("value" %in% names(nodes))) nodes$value <- 1
    if (!("carac" %in% names(nodes))) nodes$carac <- "1"

    
# ---- ensure common node columns across all metadata ----
if (!("value2" %in% names(nodes))) nodes$value2 <- NA_real_
# SS + a*
if (!("ss" %in% names(nodes))) {
  if ("SSi" %in% names(nodes)) nodes$ss <- suppressWarnings(as.numeric(nodes$SSi))
  else if ("ssi" %in% names(nodes)) nodes$ss <- suppressWarnings(as.numeric(nodes$ssi))
  else nodes$ss <- NA_real_
}
if (!("a_star" %in% names(nodes))) {
  if ("a_star1" %in% names(nodes)) nodes$a_star <- suppressWarnings(as.numeric(nodes$a_star1))
  else if ("a_star" %in% names(nodes)) nodes$a_star <- suppressWarnings(as.numeric(nodes$a_star))
  else nodes$a_star <- NA_real_
}

# ---- nearest vertex/cluster (by strongest WCD, undirected) ----
nodes$nearest_vertex <- NA_character_
nodes$nearest_cluster <- NA_character_
if (is.data.frame(data) && nrow(data)) {
  leader_col   <- intersect(c("Leader","leader","from","From"), names(data))
  follower_col <- intersect(c("Follower","follower","to","To"), names(data))
  wcd_col      <- intersect(c("WCD","wcd","weight","Weight","value"), names(data))
  if (length(leader_col) && length(follower_col) && length(wcd_col)) {
    ed2 <- data.frame(
      a = as.character(data[[leader_col[1]]]),
      b = as.character(data[[follower_col[1]]]),
      w = suppressWarnings(as.numeric(data[[wcd_col[1]]])),
      stringsAsFactors = FALSE
    )
    ed2 <- ed2[is.finite(ed2$w) & nzchar(ed2$a) & nzchar(ed2$b), , drop=FALSE]
    if (nrow(ed2)) {
      # make undirected duplicates
      edU <- rbind(ed2, data.frame(a=ed2$b, b=ed2$a, w=ed2$w, stringsAsFactors=FALSE))
      # for each a, pick b with max w
      ord <- order(edU$a, -edU$w)
      edU <- edU[ord, , drop=FALSE]
      first_idx <- !duplicated(edU$a)
      best <- edU[first_idx, , drop=FALSE]
      idx <- match(nodes$name, best$a)
      nodes$nearest_vertex <- ifelse(is.na(idx), NA_character_, best$b[idx])
      # nearest cluster from nodes table
      if ("carac" %in% names(nodes)) {
        map_carac <- setNames(as.character(nodes$carac), as.character(nodes$name))
        nodes$nearest_cluster <- unname(map_carac[nodes$nearest_vertex])
      }
    }
  }
}


    # ---- global n(PubMed) total for ratio (scalar only) ----
    n_total <- suppressWarnings(as.numeric(n_pubmed_total))
    if (!is.finite(n_total) || length(n_total) != 1) n_total <- NA_real_
    ratio_main <- rep(NA_real_, nrow(nodes))
    warn_main  <- rep(NA_character_, nrow(nodes))
    if (is.finite(n_total) && n_total > 0) {
      vv <- suppressWarnings(as.numeric(nodes$value))
      ratio_main <- vv / n_total
      warn_main  <- ifelse(is.finite(ratio_main) & ratio_main < 0.5, "UNHEALTHFUL (primary-author ratio <0.5)", "OK")
    } else {
      warn_main <- "NA"
    }
vn_nodes <- data.frame(
      id    = as.character(nodes$name),
      label = as.character(nodes$name),
      value = pmax(1, suppressWarnings(as.numeric(nodes$value))),
      group = as.character(nodes$carac),
      title = paste0(
        "<b>", nodes$name, "</b>",
        "<br>n(PubMed total)=", ifelse(is.finite(n_total), as.integer(n_total), "NA"),
        "<br>n(any byline appearance)=", nodes$n_byline,
        "<br>n(first/last author appearance)=", nodes$value,
        "<br>value2(first/last excl. single)=", nodes$value2,
        "<br>single=value-value2=", nodes$single_author,
        "<br>ratio=value/nPubMedTotal=", ifelse(is.finite(ratio_main), sprintf("%.2f", ratio_main), "NA"),
        "<br>middle-author warning=", warn_main,
        "<br>value2_strength=", sprintf("%.0f", nodes$value2_strength),
        "<br>ss=", ifelse(is.finite(nodes$ss), sprintf("%.2f", nodes$ss), "NA"),
        "<br>a*=", ifelse(is.finite(nodes$a_star), sprintf("%.2f", nodes$a_star), "NA"),
        "<br>nearest_cluster=", nodes$nearest_cluster,
        "<br>nearest_vertex=", nodes$nearest_vertex,
        "<br>cluster=", nodes$carac
      ),
      stringsAsFactors = FALSE
    )

    # groups/colors: make cluster colors explicit and stable
    grp_levels <- sort(unique(vn_nodes$group))
    grp_cols <- NULL
    if (length(grp_levels) > 0) {
      grp_cols <- grDevices::hcl(
        h = seq(15, 375, length.out = length(grp_levels) + 1)[-(length(grp_levels) + 1)],
        c = 100, l = 55
      )
      names(grp_cols) <- grp_levels
    }

    # edges
    if (!is.data.frame(data) || !nrow(data)) {
      return(visNetwork(vn_nodes, data.frame()))
    }
    leader_col   <- intersect(c("Leader","leader","from","From"), names(data))
    follower_col <- intersect(c("Follower","follower","to","To"), names(data))
    wcd_col      <- intersect(c("WCD","wcd","weight","Weight","value"), names(data))

    from_vec <- if (length(leader_col)) as.character(data[[leader_col[1]]]) else character(0)
    to_vec   <- if (length(follower_col)) as.character(data[[follower_col[1]]]) else character(0)
    w_vec    <- if (length(wcd_col)) suppressWarnings(as.numeric(data[[wcd_col[1]]])) else rep(1, length(from_vec))

    # If we still have mismatched lengths, drop to a safe empty-edge view instead of crashing
    if (length(from_vec) != length(to_vec) || length(from_vec) != length(w_vec)) {
      warning("render_vn(): edge columns have different lengths; rendering nodes only.")
      return(visNetwork(vn_nodes, data.frame()))
    }

    vn_edges <- data.frame(
      from   = from_vec,
      to     = to_vec,
      value  = pmax(1, w_vec),
      width  = if (length(unique(w_vec)) > 1) {
        1 + 6 * (w_vec - min(w_vec, na.rm=TRUE)) / (max(w_vec, na.rm=TRUE) - min(w_vec, na.rm=TRUE))
      } else {
        rep(2, length(w_vec))
      },
      title  = paste0("WCD=", w_vec),
      arrows = if (directed) "to" else "",
      stringsAsFactors = FALSE
    )

    p <- visNetwork(vn_nodes, vn_edges)
    if (!is.null(grp_cols) && length(grp_cols) > 0) {
      for (g in names(grp_cols)) {
        p <- visGroups(p, groupname = g, color = list(background = grp_cols[[g]], border = grp_cols[[g]]))
      }
    }
    p <- visOptions(p, highlightNearest = TRUE, nodesIdSelection = TRUE)
    p <- visInteraction(p, navigationButtons = TRUE)
    p <- visPhysics(p, stabilization = TRUE)
    p <- visEdges(p, smooth = list(enabled = TRUE, type = "dynamic"))

    if (bold) {
      p <- visNodes(p, font = list(size = label_size, face = "bold", color = "black"))
    } else {
      p <- visNodes(p, font = list(size = label_size, color = "black"))
    }
    p
  }

  output$vn_author <- renderVisNetwork({
    req(rv$author_nodes, rv$author_edges)
    render_vn(rv$author_nodes, rv$author_edges, directed = TRUE, label_size = input$label_size, bold = TRUE, n_pubmed_total = length(rv$pmids))
  })
  output$vn_country <- renderVisNetwork({
    req(rv$country_nodes, rv$country_edges)
    render_vn(rv$country_nodes, rv$country_edges, directed = FALSE, label_size = 18, bold = TRUE)
  })
  output$vn_stateprov <- renderVisNetwork({
    req(rv$stateprov_nodes, rv$stateprov_edges)
    render_vn(rv$stateprov_nodes, rv$stateprov_edges, directed = FALSE, label_size = 16, bold = TRUE)
  })

  # ---- China map (Provinces & Municipalities) ----
  .china_lookup <- data.frame(
    name_cn = c("北京","天津","河北","山西","内蒙古","辽宁","吉林","黑龙江","上海","江苏","浙江","安徽","福建","江西","山东","河南","湖北","湖南","广东","广西","海南","重庆","四川","贵州","云南","西藏","陕西","甘肃","青海","宁夏","新疆","台湾","香港","澳门","南海诸岛"),
    name_en = c("Beijing","Tianjin","Hebei","Shanxi","Inner Mongolia","Liaoning","Jilin","Heilongjiang","Shanghai","Jiangsu","Zhejiang","Anhui","Fujian","Jiangxi","Shandong","Henan","Hubei","Hunan","Guangdong","Guangxi","Hainan","Chongqing","Sichuan","Guizhou","Yunnan","Tibet","Shaanxi","Gansu","Qinghai","Ningxia","Xinjiang","Taiwan","Hong Kong","Macau","South China Sea Islands"),
    stringsAsFactors = FALSE
  )

  .match_china_province <- function(x) {
    x0 <- trimws(as.character(x))
    if (!nzchar(x0)) return(NA_character_)
    # exact CN match
    w <- match(x0, .china_lookup$name_cn)
    if (!is.na(w)) return(.china_lookup$name_cn[w])
    # EN match (case-insensitive, allow a few common variants)
    x1 <- tolower(x0)
    en <- tolower(.china_lookup$name_en)
    # normalize a couple of variants
    x1 <- gsub("\\s+", " ", x1)
    x1 <- gsub("south_?china_?sea_?islands|south china see ilands", "south china sea islands", x1)
    x1 <- gsub("^hong$", "hong kong", x1)
    x1 <- gsub("^inner$", "inner mongolia", x1)
    w2 <- match(x1, en)
    if (!is.na(w2)) return(.china_lookup$name_cn[w2])
    # partial fallback for inner mongolia
    if (grepl("inner mongolia", x1, fixed=TRUE)) return("内蒙古")
    if (grepl("hong kong", x1, fixed=TRUE)) return("香港")
    if (grepl("macao|macau", x1)) return("澳门")
    if (grepl("taiwan", x1, fixed=TRUE)) return("台湾")
    NA_character_
  }

  china_counts <- shiny::reactive({
    # Prefer nodes$value from the State/Province nodes table (pre-FLCA distribution),
    # fallback to raw term-list when nodes are not available.
    sp <- NULL

    if (!is.null(rv$stateprov_nodes) && is.data.frame(rv$stateprov_nodes) &&
        all(c("name","value") %in% names(rv$stateprov_nodes))) {
      sp <- dplyr::transmute(
        rv$stateprov_nodes,
        st = trimws(as.character(name)),
        n  = suppressWarnings(as.numeric(value))
      )
    } else if (!is.null(rv$stateprov_list)) {
      x <- rv$stateprov_list
      if (is.data.frame(x)) {
        nm <- intersect(c("name","term","state","State"), names(x))
        ct <- intersect(c("value","count","n","freq","Count"), names(x))
        if (length(nm) >= 1 && length(ct) >= 1) {
          sp <- dplyr::transmute(x,
                                 st = trimws(as.character(.data[[nm[1]]])),
                                 n  = suppressWarnings(as.numeric(.data[[ct[1]]])))
        }
      } else if (is.list(x)) {
        v <- unlist(x, use.names = FALSE)
        v <- trimws(as.character(v))
        v <- v[nzchar(v)]
        if (length(v)) {
          tb <- sort(table(v), decreasing = TRUE)
          sp <- data.frame(st = names(tb), n = as.numeric(tb), stringsAsFactors = FALSE)
        }
      }
    }

    if (is.null(sp) || !is.data.frame(sp) || nrow(sp) == 0) {
      return(data.frame(name_cn=character(), name_en=character(), value=numeric(), stringsAsFactors=FALSE))
    }

    sp$st <- trimws(as.character(sp$st))
    sp$n  <- suppressWarnings(as.numeric(sp$n))
    sp <- sp[is.finite(sp$n) & nzchar(sp$st), , drop=FALSE]
    if (nrow(sp) == 0) {
      return(data.frame(name_cn=character(), name_en=character(), value=numeric(), stringsAsFactors=FALSE))
    }

    # Match to China provinces/SARs (accept CN/EN/with-country strings)
    hit <- vapply(sp$st, .match_china_province, character(1))
    keep <- !is.na(hit)
    if (!any(keep)) {
      return(data.frame(name_cn=character(), name_en=character(), value=numeric(), stringsAsFactors=FALSE))
    }

    sp2 <- data.frame(name_cn = hit[keep], value = sp$n[keep], stringsAsFactors = FALSE)
    sp2 <- sp2 %>%
      dplyr::group_by(name_cn) %>%
      dplyr::summarise(value = sum(value, na.rm = TRUE), .groups = "drop") %>%
      dplyr::arrange(dplyr::desc(value))

    out <- merge(sp2, .china_lookup, by.x="name_cn", by.y="name_cn", all.x=TRUE, sort=FALSE)
    out <- out[order(out$value, decreasing=TRUE), c("name_cn","name_en","value")]
    out
  })

  output$tbl_china_counts <- DT::renderDT({
    df <- china_counts()
    if (!is.data.frame(df) || nrow(df)==0) {
      return(DT::datatable(data.frame(Message="No China province/SAR terms found in State/Province list.", stringsAsFactors=FALSE),
                           options=list(dom='t'), rownames=FALSE))
    }
    DT::datatable(df, options=list(pageLength=20), rownames=FALSE)
  })

  output$china_map <- shiny::renderUI({
    df <- china_counts()

    if (!requireNamespace("hchinamap", quietly = TRUE)) {
      return(shiny::tags$div(class = "alert alert-warning",
                            "China map requires package: hchinamap. Please install.packages('hchinamap')."))
    }

    if (!is.data.frame(df) || nrow(df) == 0) {
      return(shiny::tags$div(class = "alert alert-info",
                            "China map: no China province/SAR terms found."))
    }

    # hchinamap expects Chinese province names in `name`
    hchinamap::hchinamap(
      name   = df$name_cn,
      value  = df$value,
      region = "China",
      width  = "100%",
      height = "650px",
      title  = "Map of China",
      minColor = "#f1eef6",
      maxColor = "#980043"
    )
  })

output$vn_inst <- renderVisNetwork({
    req(rv$inst_nodes, rv$inst_edges)
    render_vn(rv$inst_nodes, rv$inst_edges, directed = FALSE, label_size = 16, bold = TRUE)
  })
  output$vn_dept <- renderVisNetwork({
    req(rv$dept_nodes, rv$dept_edges)
    render_vn(rv$dept_nodes, rv$dept_edges, directed = FALSE, label_size = 16, bold = TRUE)
  })
  output$vn_mesh <- renderVisNetwork({
    req(rv$mesh_nodes, rv$mesh_edges)
    render_vn(rv$mesh_nodes, rv$mesh_edges, directed = FALSE, label_size = 16, bold = TRUE)
  })

  # ==========================================================
  # Combo (single selector across metadata domains)
  # ==========================================================
  .combo_get <- reactive({
    dom <- input$combo_domain %||% "Author"
    nodes <- switch(dom,
      "Author"          = rv$author_nodes,
      "Country"         = rv$country_nodes,
      "State/Province"  = rv$stateprov_nodes,
      "Institute"       = rv$inst_nodes,
      "Department"      = rv$dept_nodes,
      "MeSH"            = rv$mesh_nodes,
      "Journal/Year"    = (rv$jy_nodes %||% rv$journal_nodes),
      rv$author_nodes
    )
    edges <- switch(dom,
      "Author"          = rv$author_edges,
      "Country"         = rv$country_edges,
      "State/Province"  = rv$stateprov_edges,
      "Institute"       = rv$inst_edges,
      "Department"      = rv$dept_edges,
      "MeSH"            = rv$mesh_edges,
      "Journal/Year"    = (rv$jy_edges %||% rv$journal_edges),
      rv$author_edges
    )
    list(dom=dom, nodes=nodes, edges=edges)
  })

  output$vn_combo <- renderVisNetwork({
    x <- .combo_get()
    req(x$nodes, x$edges)
    directed <- identical(x$dom, "Author")
    render_vn(x$nodes, x$edges, directed = directed, label_size = input$label_size %||% 16, bold = TRUE)
  })

  output$combo_nodes <- renderDT({
    x <- .combo_get()
    if (is.null(x$nodes) || !is.data.frame(x$nodes)) return(DT::datatable(data.frame()))
    DT::datatable(x$nodes, options = list(pageLength = 20, scrollX = TRUE))
  })
  output$combo_edges <- renderDT({
    x <- .combo_get()
    if (is.null(x$edges) || !is.data.frame(x$edges)) return(DT::datatable(data.frame()))
    DT::datatable(x$edges, options = list(pageLength = 20, scrollX = TRUE))
  })


  # ==========================================================
  # Chord diagram (from selected domain Top20 nodes/edges)
  # ==========================================================
  .chord_cluster_base_palette <- c(
    "#96ab03", "#6495ED", "#FF7F50", "#CCCCFF", "#9FE2BF",
    "#40E0D0", "#FFBF00", "#f00000", "#7ABF1B", "#C6CCBD",
    "#8A2BE2", "#00CED1", "#FFD700", "#DDA0DD", "#00FF7F",
    "#1E90FF", "#DC143C", "#00BFFF", "#FF69B4", "#A52A2A"
  )
  .chord_cluster_palette_map <- function(cluster_levels){
    clv <- unique(as.character(cluster_levels))
    clv <- clv[!is.na(clv) & nzchar(clv)]
    if (!length(clv)) clv <- "C1"
    cols <- .chord_cluster_base_palette[((seq_along(clv)-1) %% length(.chord_cluster_base_palette)) + 1]
    setNames(cols, clv)
  }

    # ---- Chord helpers (ported from appwos.R, adapted to this app's node/edge objects) ----

  .build_chord_matrix <- function(nodes, edges20) {
    if (is.null(nodes) || !is.data.frame(nodes) || nrow(nodes) < 2) return(NULL)
    if (is.null(edges20) || !is.data.frame(edges20) || nrow(edges20) < 1) return(NULL)

    # clean labels to avoid mismatch (e.g., "name<br>cluster=...")
    if ("name" %in% names(nodes)) nodes$name <- gsub("<br>.*$", "", as.character(nodes$name))
    else if ("label" %in% names(nodes)) nodes$name <- gsub("<br>.*$", "", as.character(nodes$label))
    else if ("id" %in% names(nodes)) nodes$name <- as.character(nodes$id)
    else nodes$name <- as.character(nodes[[1]])

    nodes <- nodes[!is.na(nodes$name) & nzchar(nodes$name), , drop = FALSE]
    if (nrow(nodes) < 2) return(NULL)

    # ensure unique node names for sectors
    nodes <- nodes[!duplicated(nodes$name), , drop = FALSE]
    if (nrow(nodes) < 2) return(NULL)

    # normalize edges columns
    nms <- names(edges20)
    from_col <- intersect(nms, c("Leader","from","From","Source","source","node1","Node1"))[1]
    to_col   <- intersect(nms, c("Follower","to","To","Target","target","node2","Node2"))[1]
    w_col    <- intersect(nms, c("WCD","weight","Weight","value","Value","n_pubmed","n","freq"))[1]
    if (is.na(from_col) || is.na(to_col)) return(NULL)
    if (is.na(w_col)) { edges20$..w <- 1; w_col <- "..w" }

    from <- as.character(edges20[[from_col]])
    to   <- as.character(edges20[[to_col]])
    w    <- suppressWarnings(as.numeric(edges20[[w_col]]))
    w[!is.finite(w)] <- 0

    keep <- from %in% nodes$name & to %in% nodes$name & w > 0
    from <- from[keep]; to <- to[keep]; w <- w[keep]
    if (!length(w)) return(NULL)

    # build symmetric matrix (undirected)
    mat <- matrix(0, nrow = nrow(nodes), ncol = nrow(nodes),
                  dimnames = list(nodes$name, nodes$name))
    for (i in seq_along(w)) {
      a <- match(from[i], nodes$name)
      b <- match(to[i], nodes$name)
      if (!is.na(a) && !is.na(b) && a != b) {
        mat[a, b] <- mat[a, b] + w[i]
        mat[b, a] <- mat[b, a] + w[i]
      }
    }
    mat
  }

  # Safe wrapper for chorddiag across versions
  .safe_chorddiag_widget <- function(mat, group, groupColors, groupnamePadding = 20) {
    if (!requireNamespace("chorddiag", quietly = TRUE)) return(NULL)
    fn <- chorddiag::chorddiag
    fml <- names(formals(fn))

    args <- list()

    if ("x" %in% fml) {
      args$x <- mat
    } else if ("mat" %in% fml) {
      args$mat <- mat
    } else {
      args[[1]] <- mat
    }

    if ("group" %in% fml) args$group <- group

    if ("groupColors" %in% fml) {
      args$groupColors <- groupColors
    } else if ("groupColours" %in% fml) {
      args$groupColours <- groupColors
    } else if ("col" %in% fml) {
      args$col <- groupColors
    }

    if ("groupnamePadding" %in% fml) {
      args$groupnamePadding <- groupnamePadding
    } else if ("groupNamePadding" %in% fml) {
      args$groupNamePadding <- groupnamePadding
    }

    do.call(fn, args)
  }

  .cluster_colors_from_carac <- function(nodes) {
    if (is.null(nodes) || !is.data.frame(nodes) || nrow(nodes) < 1) return(NULL)

    # accept either carac or cluster
    if (!("carac" %in% names(nodes)) && ("cluster" %in% names(nodes))) nodes$carac <- nodes$cluster
    if (!("carac" %in% names(nodes))) nodes$carac <- "C1"

    # node names
    if (!("name" %in% names(nodes))) {
      if ("label" %in% names(nodes)) nodes$name <- as.character(nodes$label)
      else if ("id" %in% names(nodes)) nodes$name <- as.character(nodes$id)
      else nodes$name <- as.character(nodes[[1]])
    }

    nodes$name <- gsub("<br>.*$", "", as.character(nodes$name))
    nodes$carac <- as.character(nodes$carac)
    nodes$carac[is.na(nodes$carac) | !nzchar(nodes$carac)] <- "C1"

    # unique node names; keep first occurrence
    keep <- !is.na(nodes$name) & nzchar(nodes$name) & !duplicated(nodes$name)
    nodes <- nodes[keep, , drop = FALSE]
    if (nrow(nodes) < 1) return(NULL)

    pal <- .chord_cluster_palette_map(nodes$carac)
    cols <- unname(pal[nodes$carac])
    names(cols) <- nodes$name
    cols
  }

  # Sync chord domain selector with combo network selector (like appwos.R)
  observeEvent(input$combo_domain, {
    req(input$combo_domain)
    if (isTRUE(input$chord_follow_combo)) {
      updateSelectInput(session, "chord_domain", selected = input$combo_domain)
    }
  }, ignoreInit = TRUE)

  observeEvent(input$chord_domain, {
    req(input$chord_domain)
    if (isTRUE(input$chord_follow_combo)) {
      updateSelectInput(session, "combo_domain", selected = input$chord_domain)
    }
  }, ignoreInit = TRUE)

  .chord_get <- reactive({
    dom <- if (isTRUE(input$chord_follow_combo)) (input$combo_domain %||% "Author") else (input$chord_domain %||% "Author")

    nodes <- switch(dom,
      "Author"          = rv$author_nodes,
      "Country"         = rv$country_nodes,
      "State/Province"  = rv$stateprov_nodes,
      "Institute"       = rv$inst_nodes,
      "Department"      = rv$dept_nodes,
      "MeSH"            = rv$mesh_nodes,
      "Journal/Year"    = (rv$jy_nodes %||% rv$journal_nodes),
      rv$author_nodes
    )

    edges20 <- switch(dom,
      "Author"          = rv$author_edges,
      "Country"         = rv$country_edges,
      "State/Province"  = rv$stateprov_edges,
      "Institute"       = rv$inst_edges,
      "Department"      = rv$dept_edges,
      "MeSH"            = rv$mesh_edges,
      "Journal/Year"    = (rv$jy_edges %||% rv$journal_edges),
      rv$author_edges
    )

    list(dom=dom, nodes=nodes, edges20=edges20)
  })

  # ---- Chord payload (matrix + carac colors aligned to sectors) ----
  chord_payload <- reactive({
    x <- .chord_get()
    if (is.null(x$nodes) || !is.data.frame(x$nodes) || is.null(x$edges20) || !is.data.frame(x$edges20)) {
      return(list(ok=FALSE, reason="no data"))
    }

    mat <- .build_chord_matrix(x$nodes, x$edges20)
    if (is.null(mat) || nrow(mat) < 2 || sum(mat) <= 0) return(list(ok=FALSE, reason="not enough"))

    # colors from nodes$carac aligned to rownames(mat)
    cols0 <- .cluster_colors_from_carac(x$nodes)
    rn <- rownames(mat)
    cols <- NULL
    if (!is.null(cols0) && length(cols0)) {
      cols <- cols0[match(rn, names(cols0))]
      names(cols) <- rn
      cols[is.na(cols) | !nzchar(cols)] <- "#B0B0B0"
    }

    # build map for debug
    nodes2 <- x$nodes
    if (!("name" %in% names(nodes2))) {
      if ("label" %in% names(nodes2)) nodes2$name <- as.character(nodes2$label)
      else if ("id" %in% names(nodes2)) nodes2$name <- as.character(nodes2$id)
      else nodes2$name <- as.character(nodes2[[1]])
    }
    nodes2$name <- gsub("<br>.*$", "", as.character(nodes2$name))
    if (!("carac" %in% names(nodes2)) && ("cluster" %in% names(nodes2))) nodes2$carac <- nodes2$cluster
    if (!("carac" %in% names(nodes2))) nodes2$carac <- "C1"
    carac_by_name <- setNames(as.character(nodes2$carac), as.character(nodes2$name))

    df_map <- data.frame(
      node = rn,
      carac = unname(carac_by_name[rn]),
      color = if (is.null(cols)) NA_character_ else unname(cols),
      stringsAsFactors = FALSE
    )
    df_map$carac[is.na(df_map$carac) | !nzchar(df_map$carac)] <- "C1"

    list(
      ok = TRUE,
      dom = x$dom,
      mat = mat,
      cols = cols,
      map = df_map,
      edges_cols = names(x$edges20),
      n_nodes = nrow(x$nodes),
      n_edges = nrow(x$edges20)
    )
  })

  output$chord_debug <- renderText({
    payload <- chord_payload()
    if (is.null(payload) || isFALSE(payload$ok)) {
      x <- .chord_get()
      return(paste0(
        "Chord: not enough data. domain=", x$dom,
        " | nodes=", if (is.null(x$nodes)) "NULL" else nrow(x$nodes),
        " | edges20=", if (is.null(x$edges20)) "NULL" else nrow(x$edges20),
        " | edges cols=", if (is.null(x$edges20)) "NULL" else paste(names(x$edges20), collapse=",")
      ))
    }
    paste0(
      "Chord OK | domain=", payload$dom,
      " | sectors=", nrow(payload$mat),
      " | sumW=", sum(payload$mat),
      " | colored=", if (is.null(payload$cols)) 0 else sum(!is.na(payload$cols)),
      " | edges cols=", paste(payload$edges_cols, collapse=",")
    )
  })

  output$chord_ui <- renderUI({
    payload <- chord_payload()
    if (is.null(payload) || isFALSE(payload$ok)) {
      return(tags$div(class="small-note","Chord: not enough data yet (need nodes + edges)."))
    }

    # choose rendering path
    if (requireNamespace("chorddiag", quietly = TRUE)) {
      rn <- rownames(payload$mat)

      # group per sector = nodes$carac (cluster) aligned to rn
      grp <- payload$map$carac
      names(grp) <- payload$map$node
      grp <- as.character(grp[rn])
      grp[is.na(grp) | !nzchar(grp)] <- "C1"
      grp_levels <- unique(grp)

      # sector colors aligned to rn (fallback gray)
      sector_cols <- payload$cols
      if (is.null(sector_cols) || !length(sector_cols)) {
        sector_cols <- rep("#B0B0B0", length(rn))
        names(sector_cols) <- rn
      } else {
        sector_cols <- sector_cols[rn]
        sector_cols[is.na(sector_cols) | !nzchar(sector_cols)] <- "#B0B0B0"
        names(sector_cols) <- rn
      }

      # groupColors: one color per group (take first sector color in that group)
      group_cols <- vapply(grp_levels, function(g) {
        idx <- which(grp == g)[1]
        unname(sector_cols[idx])
      }, FUN.VALUE = character(1))

      w <- tryCatch({
        .safe_chorddiag_widget(
          mat = payload$mat,
          group = factor(grp, levels = grp_levels),
          groupColors = unname(group_cols),
          groupnamePadding = 20
        )
      }, error=function(e) NULL)

      if (!is.null(w)) return(w)
    }

    if (requireNamespace("circlize", quietly = TRUE)) {
      plotOutput("chord_plot", height = "720px")
    } else {
      tags$div(
        tags$b("Chord needs a package:"),
        tags$br(),
        "Install either ", tags$code("install.packages('chorddiag')"),
        " (interactive) or ", tags$code("install.packages('circlize')"), " (static)."
      )
    }
  })

  output$chord_plot <- renderPlot({
    payload <- chord_payload()
    validate(need(!is.null(payload) && isTRUE(payload$ok), "Chord: not enough data."))
    validate(need(requireNamespace("circlize", quietly = TRUE), "Install 'circlize' for static chord."))
    circlize::circos.clear()
    on.exit(circlize::circos.clear(), add = TRUE)
    circlize::circos.par(start.degree = 90, gap.after = 2)

    if (!is.null(payload$cols)) {
      circlize::chordDiagram(payload$mat, grid.col = payload$cols, transparency = 0.25, annotationTrack = "grid")
    } else {
      circlize::chordDiagram(payload$mat, transparency = 0.25, annotationTrack = "grid")
    }
  })

# ---- AAC + Kano plot (Author only) ----

  # ---- Safe Top20 and plot helpers for SSplot/Kano tabs ----
  .safe_top20_nodes <- function(df, n = 20) {
    if (is.null(df) || !is.data.frame(df) || !nrow(df)) return(df)
    df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
    if (!("name" %in% names(df))) {
      if ("id" %in% names(df)) df$name <- as.character(df$id) else df$name <- as.character(df[[1]])
    }
    if (!("rank" %in% names(df))) {
      if ("value" %in% names(df)) {
        vv <- suppressWarnings(as.numeric(df$value))
        ord <- order(ifelse(is.finite(vv), vv, -Inf), decreasing = TRUE, na.last = TRUE)
      } else {
        ord <- seq_len(nrow(df))
      }
      df <- df[ord, , drop = FALSE]
      df$rank <- seq_len(nrow(df))
    } else {
      rr <- suppressWarnings(as.numeric(df$rank))
      vv <- if ("value" %in% names(df)) suppressWarnings(as.numeric(df$value)) else rep(NA_real_, nrow(df))
      ord <- order(ifelse(is.finite(rr), rr, Inf),
                   ifelse(is.finite(vv), -vv, Inf),
                   na.last = TRUE)
      df <- df[ord, , drop = FALSE]
    }
    utils::head(df, n)
  }

  .safe_draw_kano_xy <- function(nd, xcol, ycol, sizecol = "value",
                                 title_txt = "Kano plot",
                                 xlab = xcol, ylab = ycol,
                                 label_size = 4,
                                 vertical = FALSE) {
    # SAFE BASE-R KANO DRAWER
    # Purpose: avoid HTML/list/ggplot objects being passed to Shiny text/HTML handlers,
    # which caused: "不是所有的 is.character(txt) 都是 TRUE".
    nd <- .safe_top20_nodes(nd, 20)
    if (!is.data.frame(nd) || nrow(nd) < 2) {
      plot.new(); text(0.5, 0.5, "Need at least 2 Top20 nodes for Kano plot.", cex = 1.3, font = 2)
      return(invisible(NULL))
    }
    if (!("name" %in% names(nd))) nd$name <- as.character(seq_len(nrow(nd)))
    if (!("carac" %in% names(nd))) nd$carac <- 1
    if (!(xcol %in% names(nd))) nd[[xcol]] <- NA_real_
    if (!(ycol %in% names(nd))) nd[[ycol]] <- NA_real_
    if (!(sizecol %in% names(nd))) nd[[sizecol]] <- 1

    nd$name <- as.character(nd$name)
    nd$carac <- as.character(nd$carac)
    nd[[xcol]] <- suppressWarnings(as.numeric(nd[[xcol]]))
    nd[[ycol]] <- suppressWarnings(as.numeric(nd[[ycol]]))
    nd[[sizecol]] <- suppressWarnings(as.numeric(nd[[sizecol]]))
    nd[[sizecol]][!is.finite(nd[[sizecol]]) | is.na(nd[[sizecol]])] <- 1

    ok <- is.finite(nd[[xcol]]) & is.finite(nd[[ycol]])
    nd <- nd[ok, , drop = FALSE]
    if (nrow(nd) < 2) {
      plot.new(); text(0.5, 0.5, "Need >=2 finite Top20 points for Kano plot.", cex = 1.3, font = 2)
      return(invisible(NULL))
    }

    x <- nd[[xcol]]; y <- nd[[ycol]]
    sx <- nd[[sizecol]]
    sx[!is.finite(sx)] <- 1
    sx0 <- sqrt(pmax(sx, 0))
    mx <- max(sx0, na.rm = TRUE)
    if (!is.finite(mx) || mx <= 0) mx <- 1
    cex_pt <- 1.5 + 4.0 * sx0 / mx

    medx <- stats::median(x, na.rm = TRUE)
    medy <- stats::median(y, na.rm = TRUE)
    ux <- range(x, medx, finite = TRUE)
    uy <- range(y, medy, finite = TRUE)
    dx <- diff(ux); dy <- diff(uy)
    if (!is.finite(dx) || dx == 0) dx <- max(abs(ux), 1)
    if (!is.finite(dy) || dy == 0) dy <- max(abs(uy), 1)
    xlim <- ux + c(-0.16, 0.22) * dx
    ylim <- uy + c(-0.20, 0.22) * dy

    cls <- sort(unique(nd$carac))
    pal <- grDevices::hcl.colors(max(3, length(cls)), "Dark 3")
    col_map <- setNames(pal[seq_along(cls)], cls)
    bg <- col_map[nd$carac]

    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar), add = TRUE)
    par(mar = c(6.0, 6.5, 5.5, 3.0), xpd = NA, family = "sans")
    plot(x, y, type = "n", xlim = xlim, ylim = ylim,
         xlab = xlab, ylab = ylab, main = title_txt,
         cex.main = 1.8, cex.lab = 1.55, cex.axis = 1.25,
         font.main = 2, font.lab = 2)
    grid(col = "grey88", lty = 1)
    abline(v = medx, h = medy, lty = 2, col = "red", lwd = 2.4)
    if (length(unique(x)) >= 2 && length(unique(y)) >= 2) {
      try(abline(stats::lm(y ~ x), col = "blue", lwd = 3.6), silent = TRUE)
    }

    # Kano-style two wings through the median cross, scaled to current data range.
    slope1 <- dy / dx
    abline(a = medy - slope1 * medx, b = slope1, lty = 3, col = "grey45", lwd = 1.6)
    abline(a = medy + slope1 * medx, b = -slope1, lty = 3, col = "grey45", lwd = 1.6)

    # Two reference circles/ellipses around the center. Base R circles distort if x/y scales differ;
    # use ellipse coordinates so the guides stay visible in a tall plot.
    th <- seq(0, 2*pi, length.out = 241)
    for (rr in c(0.28, 0.50)) {
      lines(medx + rr * dx * cos(th), medy + rr * dy * sin(th), lty = 2, col = "grey50", lwd = 1.3)
    }

    points(x, y, pch = 21, bg = bg, col = "black", cex = cex_pt, lwd = 1.8)
    text(x, y, labels = nd$name, pos = 3,
         cex = max(0.75, min(1.3, as.numeric(label_size) / 4.0)),
         font = 2, col = "black")
    legend("topright", legend = paste("Cluster", cls), pt.bg = col_map[cls],
           pch = 21, bty = "n", cex = 1.15, text.font = 2, title = "Cluster")
    invisible(TRUE)
  }

  output$kano_aac_header <- renderText({
    d <- .get_domain(input$kano_domain)
    sprintf("AAC(value)=%.2f | AAC(value2)=%.2f | AAC(SS)=%.2f | AAC(a*)=%.2f",
            as.numeric(d$AAC_value), as.numeric(d$AAC_value2), as.numeric(d$AAC_ss), as.numeric(d$AAC_a_star))
  })
  
  output$kano_vv2_plot <- renderPlot({
    d <- .get_domain(input$kano_domain)
    nd <- .safe_top20_nodes(d$nodes, 20)
    ed <- d$edges
    validate(need(is.data.frame(nd) && nrow(nd) > 1, "Need at least 2 Top20 nodes for Kano plot."))
    if (!("value2" %in% names(nd))) {
      if (exists("compute_value2_strength", mode="function")) {
        nd <- tryCatch(compute_value2_strength(nd, ed), error=function(e) nd)
      }
    }
    if (!("value2" %in% names(nd))) nd$value2 <- nd$value
    if (!("value" %in% names(nd))) nd$value <- 1
    title_txt <- sprintf("Kano: value vs value2 | AAC=%.2f | AAC2=%.2f",
                         as.numeric(d$AAC_value), as.numeric(d$AAC_value2))
    .safe_draw_kano_xy(nd, xcol = "value2", ycol = "value", sizecol = "value",
                       title_txt = title_txt, xlab = "value2", ylab = "value",
                       label_size = input$kano_label_size %||% 4,
                       vertical = TRUE)
  }, height = 1050)

output$kano_ss_astar_plot <- renderPlot({
    d <- .get_domain(input$kano_domain)
    nd <- .safe_top20_nodes(d$nodes, 20)
    validate(need(is.data.frame(nd) && nrow(nd) > 1, "Need at least 2 Top20 nodes for Kano plot."))

    if (!("ssi" %in% names(nd)) && ("SSi" %in% names(nd))) nd$ssi <- nd$SSi
    if (!("a_star1" %in% names(nd)) && ("a_star" %in% names(nd))) nd$a_star1 <- nd$a_star
    if (!("a_star1" %in% names(nd))) {
      if ("a_i" %in% names(nd)) nd$a_star1 <- 1/(1+as.numeric(nd$a_i)) else nd$a_star1 <- NA_real_
    }
    nd$ssi <- suppressWarnings(as.numeric(nd$ssi))
    nd$a_star1 <- suppressWarnings(as.numeric(nd$a_star1))
    if (!("value" %in% names(nd))) nd$value <- 1
    aac_ss  <- if (!is.null(d$AAC_ss)) as.numeric(d$AAC_ss) else safe_AAC(nd$ssi)
    aac_ast <- if (!is.null(d$AAC_a_star)) as.numeric(d$AAC_a_star) else safe_AAC(nd$a_star1)
    title_txt <- sprintf("Kano: SS vs a* | AAC(SS)=%.2f | AAC(a*)=%.2f", aac_ss, aac_ast)

    .safe_draw_kano_xy(nd, xcol = "a_star1", ycol = "ssi", sizecol = "value",
                       title_txt = title_txt, xlab = "a*", ylab = "SS",
                       label_size = input$kano_label_size %||% 4,
                       vertical = TRUE)
  }, height = 1050)


  # ---- SS Kano tab  # ---- SS Kano tab (real Kano: 2 wings + 2 circles; bubble=value; color=cluster) ----
  output$ss_kano_plot <- renderPlot({
    dom <- input$ss_kano_domain %||% "Author"
    d <- .get_domain(dom)
    nd <- .safe_top20_nodes(d$nodes, 20)
    validate(need(is.data.frame(nd) && nrow(nd) > 1, "Need at least 2 Top20 nodes."))

    if (!("ssi" %in% names(nd)) && ("SSi" %in% names(nd))) nd$ssi <- nd$SSi
    if (!("a_star1" %in% names(nd)) && ("a_star" %in% names(nd))) nd$a_star1 <- nd$a_star
    if (!("a_star1" %in% names(nd)) && ("a_i" %in% names(nd))) nd$a_star1 <- 1/(1+as.numeric(nd$a_i))
    nd$ssi     <- suppressWarnings(as.numeric(nd$ssi))
    nd$a_star1 <- suppressWarnings(as.numeric(nd$a_star1))
    if (!("value" %in% names(nd))) nd$value <- 1
    aac_ss  <- safe_AAC(nd$ssi)
    aac_ast <- safe_AAC(nd$a_star1)
    title_txt <- sprintf("SS Kano: SS vs a* | AAC(SS)=%.2f | AAC(a*)=%.2f", aac_ss, aac_ast)

    .safe_draw_kano_xy(nd, xcol = "a_star1", ycol = "ssi", sizecol = "value",
                       title_txt = title_txt, xlab = "a*", ylab = "SS",
                       label_size = input$ss_kano_label_size %||% 4,
                       vertical = TRUE)
  }, height = 1050)


output$taaa_table <- renderTable({
  validate(need(!is.null(rv$taaa_df) && is.data.frame(rv$taaa_df) && nrow(rv$taaa_df) > 0,
                "Run PubMed / upload first, then build the MeSH domain to generate the TAAA table."))
  rv$taaa_df
}, striped = TRUE, bordered = TRUE, hover = TRUE, spacing = "s")

output$taaa_theme_freq_table <- renderTable({
  validate(need(!is.null(rv$taaa_freq_df) && is.data.frame(rv$taaa_freq_df) && nrow(rv$taaa_freq_df) > 0,
                "Run PubMed / upload first, then build the MeSH domain to generate the theme-frequency table."))
  rv$taaa_freq_df
}, striped = TRUE, bordered = TRUE, hover = TRUE, spacing = "s")

output$taaa_profile_map_table <- renderTable({
  validate(need(!is.null(rv$taaa_profile_map_df) && is.data.frame(rv$taaa_profile_map_df) && nrow(rv$taaa_profile_map_df) > 0,
                "Shown when a Profile column exists and theme count matches the profile count."))
  rv$taaa_profile_map_df
}, striped = TRUE, bordered = TRUE, hover = TRUE, spacing = "s")

output$taaa_kappa_table <- renderTable({
  validate(need(!is.null(rv$taaa_kappa_df) && is.data.frame(rv$taaa_kappa_df) && nrow(rv$taaa_kappa_df) > 0,
                "Shown when a Profile column exists and theme count matches the profile count."))
  rv$taaa_kappa_df
}, striped = TRUE, bordered = TRUE, hover = TRUE, spacing = "s")

output$taaa_confusion_table <- renderTable({
  validate(need(!is.null(rv$taaa_conf_df) && is.data.frame(rv$taaa_conf_df) && nrow(rv$taaa_conf_df) > 0,
                "Shown when a Profile column exists and theme count matches the profile count."))
  rv$taaa_conf_df
}, striped = TRUE, bordered = TRUE, hover = TRUE, spacing = "s")


# ---- TAAA2: PubMedBERT/SPECTER2 semantic clustering demo ----
.taaa2_clean_text <- function(x) {
  x <- tolower(as.character(x %||% ""))
  x <- gsub("[’']s\\b", "", x, perl = TRUE)
  x <- gsub("[^a-z0-9\\- ]+", " ", x, perl = TRUE)
  x <- gsub("\\s+", " ", x, perl = TRUE)
  trimws(x)
}

# TAAA2 term purification: aligned with en_tail_keyword_engine.R style.
# Removes verbs, weak adjectives/adverbs, metadata words, numbers, and generic fragments
# before TF-IDF cluster labels are generated.
.taaa2_stop_tokens <- unique(c(
  "the", "and", "for", "with", "using", "from", "into", "can", "could", "may", "might",
  "in", "of", "to", "a", "an", "on", "by", "as", "or", "if", "via", "per", "than",
  "is", "are", "was", "were", "be", "been", "being", "am", "do", "does", "did",
  "this", "that", "these", "those", "which", "whose", "while", "whereas", "because",
  "before", "after", "between", "among", "within", "without", "under", "over", "through",
  "also", "moreover", "furthermore", "however", "therefore", "although", "another", "many", "all",
  "used", "use", "uses", "based", "shown", "found", "performed", "included", "including", "according",
  "study", "studies", "article", "articles", "paper", "papers", "research", "result", "results",
  "method", "methods", "material", "materials", "discussion", "conclusion", "abstract", "introduction",
  "keyword", "keywords", "reference", "references", "bibliography", "metadata", "author", "journal",
  "doi", "pmid", "pmcid", "http", "https", "www", "html", "xml", "json", "csv", "tsv", "xlsx", "xls",
  "pdf", "docx", "org", "com", "io", "gov", "edu", "raw", "github", "open", "access",
  "analysis", "model", "models", "data", "dataset", "datasets", "online", "assessment", "evaluate", "evaluation",
  "number", "numbers", "view", "views", "hospital", "hospitals", "cat", "0", "1", "2", "3", "4", "5",
  "significant", "different", "important", "possible", "potential", "overall", "individual", "given", "short"
))

.taaa2_bad_starts <- c(
  "accompanied", "automatically", "collapsible", "correspondence", "facilitating", "further",
  "generating", "hierarchical", "inputs", "institutional", "intuitive", "optional", "despite",
  "including", "maintains", "demonstrating", "analyzing", "combining", "providing", "enabling",
  "supporting", "suitability", "predicting", "evaluating", "assessing"
)

.taaa2_good_single_suffix <- "(ability|ibility|bility|ization|isation|ology|omics|itis|osis|therapy|diagnosis|dermatology|bibliometrics|scientometrics)$"

.taaa2_is_bad_term <- function(term) {
  x <- .taaa2_clean_text(term)
  if (!nzchar(x)) return(TRUE)
  if (grepl("[0-9]", x, perl = TRUE)) return(TRUE)

  toks <- unlist(strsplit(x, "\\s+", perl = TRUE), use.names = FALSE)
  toks <- toks[nzchar(toks)]
  if (!length(toks)) return(TRUE)

  # Remove phrases containing generic stop/verb/meta tokens.
  if (any(toks %in% .taaa2_stop_tokens)) return(TRUE)

  # Remove very short non-abbreviation tokens.
  keep_short <- c("ai", "nlp", "r", "us", "uk")
  if (any(nchar(toks) <= 2L & !(toks %in% keep_short))) return(TRUE)

  # Remove broken fragments and non-word fragments.
  bad_frag <- c("tion", "zation", "ization", "sion", "ment", "ity", "ware", "analy", "port", "im", "sum", "marization")
  if (any(toks %in% bad_frag)) return(TRUE)
  if (any(!grepl("^[a-z][a-z-]*$", toks, perl = TRUE))) return(TRUE)

  # Remove phrases starting with weak verbs/adjectives/adverbs.
  if (grepl(paste0("^(" , paste(.taaa2_bad_starts, collapse = "|"), ")\\b"), x, perl = TRUE)) return(TRUE)

  # Single-word terms must be domain-like or meaningful suffix terms.
  if (length(toks) == 1L) {
    if (toks %in% c("ai", "nlp", "dermatology", "teledermatology", "eczema", "psoriasis", "melanoma", "biologics", "keratinocytes", "bibliometrics", "citations", "rasch", "kidmap")) return(FALSE)
    if (!grepl(.taaa2_good_single_suffix, x, perl = TRUE)) return(TRUE)
  }

  FALSE
}

.taaa2_make_ngrams <- function(text, n = 1:4) {
  words <- unlist(strsplit(.taaa2_clean_text(text), "\\s+", perl = TRUE), use.names = FALSE)
  words <- words[nzchar(words)]
  words <- words[!words %in% .taaa2_stop_tokens]

  out <- character(0)
  for (k in n) {
    if (length(words) >= k) {
      out <- c(out, sapply(seq_len(length(words) - k + 1), function(i) {
        paste(words[i:(i + k - 1)], collapse = " ")
      }))
    }
  }
  out <- unique(out)
  out <- out[!vapply(out, .taaa2_is_bad_term, logical(1))]
  out
}

.taaa2_cosine_similarity <- function(x) {
  x_norm <- x / sqrt(rowSums(x^2))
  x_norm %*% t(x_norm)
}

.taaa2_demo_data <- function() {
  data.frame(
    id = paste0("P", 1:8),
    title = c(
      "Artificial intelligence in dermatology diagnosis",
      "Deep learning for melanoma image classification",
      "Atopic dermatitis and immune dysregulation",
      "Biologic therapy for psoriasis patients",
      "Teledermatology for remote skin disease consultation",
      "Digital health platforms in dermatology care",
      "Keratinocyte inflammation in eczema",
      "Machine learning prediction of skin cancer"
    ),
    abstract = c(
      "AI models can support clinical diagnosis using dermoscopic images.",
      "Convolutional neural networks improve melanoma detection from skin images.",
      "Atopic dermatitis involves immune pathways and chronic inflammation.",
      "Biologics targeting cytokines improve outcomes in psoriasis treatment.",
      "Teledermatology enables remote consultation and access to dermatology care.",
      "Digital health tools improve communication and dermatology service delivery.",
      "Keratinocytes contribute to inflammatory responses in eczema lesions.",
      "Machine learning algorithms predict skin cancer risk using clinical features."
    ),
    stringsAsFactors = FALSE
  )
}


.taaa2_find_col <- function(df, candidates) {
  if (!is.data.frame(df)) return(NA_character_)
  nm <- names(df)
  low <- tolower(nm)
  hit <- match(tolower(candidates), low)
  hit <- hit[!is.na(hit)]
  if (length(hit)) nm[hit[1]] else NA_character_
}

.taaa2_prepare_pubmed_data <- function(rv) {
  candidates <- list(pubmeta = rv$pubmeta, pmid_tbl = rv$pmid_tbl)

  for (nm in names(candidates)) {
    df <- candidates[[nm]]
    if (!is.data.frame(df) || nrow(df) == 0) next

    title_col <- .taaa2_find_col(df, c("Title", "ArticleTitle", "article_title", "TI", "title"))
    abs_col   <- .taaa2_find_col(df, c("Abstract", "AbstractText", "AB", "abstract"))
    pmid_col  <- .taaa2_find_col(df, c("PMID", "pmid", "id", "ID"))
    if (is.na(title_col) && is.na(abs_col)) next

    title <- if (!is.na(title_col)) as.character(df[[title_col]]) else rep("", nrow(df))
    abstract <- if (!is.na(abs_col)) as.character(df[[abs_col]]) else rep("", nrow(df))
    id <- if (!is.na(pmid_col)) as.character(df[[pmid_col]]) else paste0("P", seq_len(nrow(df)))

    out <- data.frame(id = id, title = title, abstract = abstract, source = nm, stringsAsFactors = FALSE)
    out$title[is.na(out$title)] <- ""
    out$abstract[is.na(out$abstract)] <- ""
    out$text <- trimws(paste(out$title, out$abstract))
    out <- out[nzchar(out$text), , drop = FALSE]
    out <- out[!duplicated(out$id), , drop = FALSE]
    if (nrow(out) >= 3) return(out)
  }

  out <- .taaa2_demo_data()
  out$source <- "demo_fallback"
  out$text <- paste(out$title, out$abstract)
  out
}

taaa2_result <- eventReactive(input$run_taaa2, {
  req(input$taaa2_model)

  shiny::withProgress(message = "Running TAAA2 semantic clustering", value = 0, {
    shiny::incProgress(0.03, detail = "Checking required R/Python packages")

    validate(
      need(requireNamespace("reticulate", quietly = TRUE), "Package reticulate is required. Please run install.packages('reticulate')."),
      need(requireNamespace("uwot", quietly = TRUE), "Package uwot is required."),
      need(requireNamespace("igraph", quietly = TRUE), "Package igraph is required.")
    )

    shiny::incProgress(0.07, detail = "Preparing PubMed title and abstract data")
    papers2 <- .taaa2_prepare_pubmed_data(rv)
    validate(need(nrow(papers2) >= 3, "TAAA2 requires at least 3 PubMed records with title and/or abstract. Please run Fetch PubMed or upload PubMed data first."))

    model_name <- input$taaa2_model
    k <- max(2, min(as.integer(input$taaa2_k %||% 3), nrow(papers2) - 1))
    top_n_edges <- max(1, as.integer(input$taaa2_top_edges %||% 2))

    shiny::incProgress(0.10, detail = "Loading reticulate Python environment")
    reticulate::use_condaenv("r-transformers", required = TRUE)
    transformers <- reticulate::import("transformers")
    torch <- reticulate::import("torch")

    shiny::incProgress(0.10, detail = "Loading PubMedBERT / SPECTER2 model")
    tokenizer <- transformers$AutoTokenizer$from_pretrained(model_name, use_fast = FALSE)
    model <- transformers$AutoModel$from_pretrained(model_name, use_safetensors = TRUE)
    model$eval()

    py$tokenizer <- tokenizer
    py$model <- model
    py$torch <- torch
    py_run_string("\ndef get_taaa2_cls_embedding(text):\n    encoded = tokenizer(text, padding=True, truncation=True, max_length=512, return_tensors='pt')\n    with torch.no_grad():\n        output = model(input_ids=encoded['input_ids'], attention_mask=encoded['attention_mask'])\n    cls_emb = output.last_hidden_state[:, 0, :]\n    return cls_emb.squeeze().detach().cpu().numpy()\n")

    shiny::incProgress(0.25, detail = paste0("Generating semantic embeddings for ", nrow(papers2), " articles"))
    get_embedding <- function(text) as.numeric(py$get_taaa2_cls_embedding(text))
    embeddings <- do.call(rbind, lapply(seq_along(papers2$text), function(i) {
      if (i %% 5 == 0 || i == nrow(papers2)) {
        shiny::incProgress(0.001, detail = paste0("Embedding article ", i, " of ", nrow(papers2)))
      }
      get_embedding(papers2$text[[i]])
    }))
    rownames(embeddings) <- papers2$id

    shiny::incProgress(0.10, detail = "Computing cosine similarity and clusters")
    sim_mat <- .taaa2_cosine_similarity(embeddings)
    dist_mat <- as.dist(1 - sim_mat)

    set.seed(123)
    hc <- hclust(dist_mat, method = "ward.D2")
    papers2$cluster <- cutree(hc, k = k)

    shiny::incProgress(0.10, detail = "Extracting TF-IDF element terms for cluster labels")
    cluster_docs <- papers2 %>%
      dplyr::group_by(cluster) %>%
      dplyr::summarise(cluster_text = paste(text, collapse = " "), .groups = "drop")

    all_terms <- lapply(cluster_docs$cluster_text, .taaa2_make_ngrams)
    names(all_terms) <- cluster_docs$cluster

    tfidf_table <- do.call(rbind, lapply(seq_along(all_terms), function(i) {
      terms <- all_terms[[i]]
      tab <- table(terms)
      df <- sapply(names(tab), function(term) sum(sapply(all_terms, function(x) term %in% x)))
      tf <- as.numeric(tab)
      idf <- log((length(all_terms) + 1) / (df + 1)) + 1
      data.frame(cluster = as.integer(names(all_terms)[i]), term = names(tab), tfidf = tf * idf, stringsAsFactors = FALSE)
    }))

    tfidf_table <- tfidf_table %>%
      dplyr::filter(!vapply(term, .taaa2_is_bad_term, logical(1)))

    # Fallback: if a cluster becomes empty after strict filtering, use cleaned title words
    # from its centroid article rather than allowing stopwords such as was/were/article.
    if (!nrow(tfidf_table)) {
      tfidf_table <- data.frame(cluster = unique(papers2$cluster), term = paste0("Cluster ", unique(papers2$cluster)), tfidf = 1, stringsAsFactors = FALSE)
    }

    top_terms <- tfidf_table %>%
      dplyr::group_by(cluster) %>%
      dplyr::arrange(desc(tfidf), .by_group = TRUE) %>%
      dplyr::slice_head(n = 5) %>%
      dplyr::summarise(top_terms = paste(term, collapse = "; "), .groups = "drop")

    nearest_to_centroid <- function(cluster_id) {
      idx <- which(papers2$cluster == cluster_id)
      emb <- embeddings[idx, , drop = FALSE]
      centroid <- colMeans(emb)
      sims <- as.numeric(emb %*% centroid / (sqrt(rowSums(emb^2)) * sqrt(sum(centroid^2))))
      idx[which.max(sims)]
    }

    centroid_titles <- data.frame(
      cluster = sort(unique(papers2$cluster)),
      centroid_title = sapply(sort(unique(papers2$cluster)), function(cl) papers2$title[nearest_to_centroid(cl)]),
      stringsAsFactors = FALSE
    )

    cluster_labels <- top_terms %>%
      dplyr::left_join(centroid_titles, by = "cluster") %>%
      dplyr::mutate(
        cluster_name = dplyr::case_when(
          grepl("artificial|deep learning|machine learning|teledermatology|digital|melanoma", top_terms) ~ "AI and Digital Dermatology",
          grepl("atopic|eczema|keratinocyte|inflammation|immune", top_terms) ~ "Inflammatory Skin Disease",
          grepl("psoriasis|biologic|therapy|cytokines", top_terms) ~ "Psoriasis Biologic Therapy",
          TRUE ~ paste("Cluster", cluster)
        )
      )

    papers2 <- papers2 %>% dplyr::left_join(cluster_labels[, c("cluster", "cluster_name", "top_terms")], by = "cluster")

    shiny::incProgress(0.10, detail = "Generating UMAP semantic map")
    set.seed(123)
    umap_xy <- uwot::umap(embeddings, n_neighbors = 3, min_dist = 0.1, metric = "cosine")
    papers2$x <- umap_xy[, 1]
    papers2$y <- umap_xy[, 2]

    shiny::incProgress(0.10, detail = "Building semantic network edges and labels")
    edge_df <- as.data.frame(as.table(sim_mat))
    colnames(edge_df) <- c("from", "to", "weight")
    edge_df <- edge_df %>%
      dplyr::mutate(from = as.character(from), to = as.character(to), weight = as.numeric(weight)) %>%
      dplyr::filter(from != to) %>%
      dplyr::group_by(from) %>%
      dplyr::arrange(desc(weight), .by_group = TRUE) %>%
      dplyr::slice_head(n = top_n_edges) %>%
      dplyr::ungroup() %>%
      dplyr::rowwise() %>%
      dplyr::mutate(edge_id = paste(sort(c(from, to)), collapse = "_")) %>%
      dplyr::ungroup() %>%
      dplyr::distinct(edge_id, .keep_all = TRUE) %>%
      dplyr::select(from, to, weight)

    node_df <- papers2 %>%
      dplyr::select(id, title, cluster, cluster_name, top_terms) %>%
      dplyr::mutate(
        short_terms = sapply(strsplit(top_terms, ";\\s*"), function(x) paste(head(x, 3), collapse = "; ")),
        label = paste0(id, "\n", cluster_name, "\n", short_terms)
      ) %>%
      dplyr::rename(name = id)

    shiny::incProgress(0.05, detail = "Finalizing TAAA2 outputs")
    list(
      papers = papers2,
      source = unique(papers2$source)[1],
      labels = cluster_labels,
      embeddings = embeddings,
      sim_mat = sim_mat,
      edge_df = edge_df,
      node_df = node_df,
      model_name = model_name,
      top_n_edges = top_n_edges
    )
  })
})

output$taaa2_umap_plot <- renderPlot({
  res <- taaa2_result()
  validate(need(!is.null(res$papers), "Click Run TAAA2 to generate the semantic clustering map."))
  ggplot2::ggplot(res$papers, ggplot2::aes(x = x, y = y, label = id, color = cluster_name)) +
    ggplot2::geom_point(size = 4) +
    ggplot2::geom_text(vjust = -0.8, fontface = "bold") +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = paste("TAAA2 Document Clustering Using", res$model_name), color = "Cluster label")
}, height = 560)

output$taaa2_network_plot <- renderPlot({
  res <- taaa2_result()
  validate(need(!is.null(res$edge_df) && nrow(res$edge_df) > 0, "No network edges were generated."))
  g <- igraph::graph_from_data_frame(res$edge_df, vertices = res$node_df, directed = FALSE)
  set.seed(123)
  lay <- igraph::layout_with_fr(g)
  cls <- as.factor(igraph::V(g)$cluster_name)
  pal <- grDevices::rainbow(length(levels(cls)))
  cols <- pal[as.integer(cls)]
  plot(g,
       layout = lay,
       vertex.color = cols,
       vertex.size = 18,
       vertex.label = igraph::V(g)$label,
       vertex.label.cex = 0.75,
       vertex.label.font = 2,
       edge.width = pmax(1, igraph::E(g)$weight * 2.5),
       edge.color = grDevices::adjustcolor("gray40", alpha.f = 0.55),
       main = paste0("TAAA2 Semantic Similarity Network: top ", res$top_n_edges, " links per article"))
  legend("topright", legend = levels(cls), col = pal, pch = 19, bty = "n", cex = 0.9)
}, height = 660)

output$taaa2_article_table <- renderTable({
  res <- taaa2_result()
  validate(need(!is.null(res$papers), "Click Run TAAA2 to generate the article table."))
  res$papers[, c("id", "title", "cluster", "cluster_name")]
}, striped = TRUE, bordered = TRUE, hover = TRUE, spacing = "s")

output$taaa2_cluster_label_table <- renderTable({
  res <- taaa2_result()
  validate(need(!is.null(res$labels), "Click Run TAAA2 to generate cluster labels."))
  res$labels
}, striped = TRUE, bordered = TRUE, hover = TRUE, spacing = "s")

output$lotka_plot <- renderPlot({
  res <- lotka_res_reactive()
  validate(need(!is.null(res), "Run routines first to generate the Lotka tab."))
  .plot_lotka_result(res)
}, height = 460)

output$lotka_test_table <- renderTable({
  res <- lotka_res_reactive()
  validate(need(!is.null(res), "Run routines first to generate the Lotka tab."))
  .lotka_summary_table(res)
}, striped = TRUE, bordered = TRUE, spacing = "s", width = "100%")

output$lotka_observed_expected_table <- renderTable({
  res <- lotka_res_reactive()
  validate(need(!is.null(res), "Run routines first to generate the Lotka tab."))
  df <- as.data.frame(res$observed_expected_table, stringsAsFactors = FALSE)
  if ("expected_raw" %in% names(df)) df$expected_raw <- round(as.numeric(df$expected_raw), 4)
  if ("expected" %in% names(df)) df$expected <- round(as.numeric(df$expected), 4)
  if ("residual" %in% names(df)) df$residual <- round(as.numeric(df$residual), 4)
  df
}, striped = TRUE, bordered = TRUE, spacing = "s", width = "100%")

output$lotka_chisq_table <- renderTable({
  res <- lotka_res_reactive()
  validate(need(!is.null(res), "Run routines first to generate the Lotka tab."))
  df <- as.data.frame(res$test_table, stringsAsFactors = FALSE)
  if ("expected" %in% names(df)) df$expected <- round(as.numeric(df$expected), 4)
  df
}, striped = TRUE, bordered = TRUE, spacing = "s", width = "100%")

# --- Year/Articles & Journal/Year (requires rv$pubmeta) ---
  output$plt_year_articles <- renderPlot({
      # ---- harden against small plotting devices ----
      op <- par(no.readonly = TRUE)
      on.exit(par(op), add = TRUE)
      par(mfrow = c(1,1))
      par(mar = c(3,3,2,1))
      par(xpd = NA)
      w <- session$clientData$output_plt_year_articles_width
      h <- session$clientData$output_plt_year_articles_height
      req(is.numeric(w), is.numeric(h), w > 400, h > 300)

    req(rv$pubmeta)
    df <- rv$pubmeta
    validate(need("Year" %in% names(df), "Year not available in this run."))
    tb <- as.data.frame(table(df$Year), stringsAsFactors=FALSE)
    colnames(tb) <- c("Year","Articles")
    tb$Year <- suppressWarnings(as.integer(as.character(tb$Year)))
    tb <- tb[order(tb$Year), , drop=FALSE]
    ggplot(tb, aes(x=Year, y=Articles)) + geom_col() + theme_minimal()
  })
  output$tbl_year_articles <- renderDT({
    req(rv$pubmeta)
    df <- rv$pubmeta
    if (!("Year" %in% names(df))) return(datatable(data.frame()))
    tb <- as.data.frame(table(df$Year), stringsAsFactors=FALSE)
    colnames(tb) <- c("Year","Articles")
    datatable(tb, options=list(pageLength=10, scrollX=TRUE))
  })
  
  # --- Year-frequency slopegraph (combo domains Top20) ---
  output$yearcount_ui <- renderUI({
    tagList(
      fluidRow(
        column(6,
               selectizeInput(
                 "yc_domains", "Domain (Top10 each; multi-select)",
                 choices = c("Author","Journal","Country","State/Province","Institute","Department","MeSH"),
                 selected = c("Author","Country"),
                 multiple = TRUE,
                 options = list(plugins = list("remove_button"))
               )
        ),
        column(6,
               sliderInput(
                 "yc_recent_years", "Recent years window",
                 min = 5, max = 20, value = 10, step = 1
               )
        )
      ),
      fluidRow(
        column(12,
               selectizeInput(
                 "yc_items", "Optional: select items (Domain::Item). Leave blank = show all Top10",
                 choices = character(0),
                 multiple = TRUE,
                 options = list(placeholder = "e.g., Author::Smith J | Country::Taiwan", maxOptions = 5000)
               )
        )
      )
    )
  })

  
  # --- Term/Year (1st+Last authors) : identical slopegraph pipeline ---
  output$yearcount_ui_fl <- renderUI({
    tagList(
      fluidRow(
        column(6,
               selectizeInput(
                 "yc_domains_fl", "Domain (Term/Year; multi-select)",
                 choices = c("Author (1st+Last)","Country","Journal","State/Province","Institute","Department","MeSH"),
                 selected = c("Author (1st+Last)","Country"),
                 multiple = TRUE,
                 options = list(plugins = list("remove_button"))
               )
        ),
        column(6,
               sliderInput(
                 "yc_recent_years_fl", "Recent years window",
                 min = 5, max = 20, value = 10, step = 1
               )
        )
      ),
      fluidRow(
        column(12,
               selectizeInput(
                 "yc_items_fl", "Optional: select items (Domain::Item). Leave blank = show all Top10",
                 choices = character(0),
                 multiple = TRUE,
                 options = list(placeholder = "e.g., Author (1st+Last)::Smith J | Country::Taiwan", maxOptions = 5000)
               )
        )
      )
    )
  })

  observeEvent(list(input$yc_domains_fl, rv$author_nodes, rv$country_nodes, rv$stateprov_nodes, rv$inst_nodes, rv$dept_nodes, rv$mesh_nodes, rv$jy_nodes), {
    doms <- input$yc_domains_fl %||% character(0)
    if (!length(doms)) doms <- c("Author (1st+Last)","Country")
    ch <- character(0)
    for (d in doms) {
      top <- if (tolower(d) %in% c("mesh", "mesh term", "meshterm")) {
        .get_mesh_top_by_pubmeta_or_list(rv, n = 10)
      } else {
        .get_domain_top20(rv, d, n = 10)
      }
      if (length(top)) ch <- c(ch, paste0(d, "::", top))
    }
    ch <- unique(ch)
    updateSelectizeInput(session, "yc_items_fl", choices = ch, selected = character(0), server = TRUE)
  }, ignoreInit = TRUE)

  output$plt_slope_top2_author_fl <- renderPlot({
    req(rv$pubmeta)
    doms <- input$yc_domains_fl %||% c("Author (1st+Last)","Country")
    yc <- compute_combo_year_counts(rv, domains = doms, recent_n_years = input$yc_recent_years_fl %||% 10)
    validate(need(!is.null(yc) && nrow(yc) > 0, "No data for slopegraph (need Year + term lists)."))

    sel <- input$yc_items_fl %||% character(0)
    if (length(sel)) {
      dom_sel  <- sub("::.*$", "", sel)
      item_sel <- sub("^.*?::", "", sel)
      yc <- yc[yc$domain %in% dom_sel & yc$item %in% item_sel, , drop=FALSE]
      validate(need(nrow(yc) > 0, "Selected items have no counts in recent years window."))
    }

    st <- summarize_prepost_ttest2(yc)
    yc$group_label <- paste0(yc$domain, "::", yc$item)

    df_s <- yc[, c("group_label","Year","Count")]
    names(df_s) <- c("itemlab","Year","Count")
    df_t <- tufte_sort2(df_s, x="Year", y="Count", group="itemlab", min.space=0.05)

    yrs <- sort(unique(as.character(df_s$Year)))
    if (all(grepl("^\\d{4}$", yrs))) yrs <- as.character(sort(as.integer(yrs)))
    df_t$x <- factor(as.character(df_t$x), levels = yrs)

    plot_slopegraph2(
      df_t,
      st = st,
      title = paste0("Term/Year slopegraph (recent ", input$yc_recent_years_fl %||% 10, " years)"),
      label_digits = 0
    )
  })

  output$tbl_yearcount_top20_fl <- DT::renderDT({
    req(rv$pubmeta)
    doms <- input$yc_domains_fl %||% c("Author (1st+Last)","Country")
    yc <- compute_combo_year_counts(rv, domains = doms, recent_n_years = input$yc_recent_years_fl %||% 10)
    validate(need(!is.null(yc) && nrow(yc) > 0, "No year-count data."))

    sel <- input$yc_items_fl %||% character(0)
    if (length(sel)) {
      dom_sel <- sub("::.*$", "", sel)
      item_sel <- sub("^.*?::", "", sel)
      yc <- yc[yc$domain %in% dom_sel & yc$item %in% item_sel, , drop=FALSE]
      validate(need(nrow(yc) > 0, "Selected items have no counts in recent years window."))
    }
    yc$group_label <- paste0(yc$domain, "::", yc$item)

    st <- summarize_prepost_ttest2(yc)
    if (is.null(st) || !nrow(st)) {
      st <- unique(yc[, c("domain","item")])
      st$mean_pre <- st$mean_post <- st$pval <- NA_real_
      st$pval_fmt <- NA_character_
      st$trend <- "flat"
    }
    st$group_label <- paste0(st$domain, "::", st$item)

    wide <- reshape2::dcast(yc, group_label + domain + item ~ Year, value.var="Count", fill=0)
    out <- merge(st[, c("group_label","mean_pre","mean_post","pval","pval_fmt","trend")], wide, by="group_label", all.x=TRUE)
    out$pval_num <- out$pval
    out$pval <- out$pval_fmt
    out <- out[order(is.na(out$pval_num), out$pval_num, -out$mean_post), , drop=FALSE]

    DT::datatable(out, options=list(pageLength=10, scrollX=TRUE))
  })

observeEvent(list(input$yc_domains, rv$author_nodes, rv$country_nodes, rv$stateprov_nodes, rv$inst_nodes, rv$dept_nodes, rv$mesh_nodes, rv$jy_nodes), {
    doms <- input$yc_domains %||% character(0)
    if (!length(doms)) doms <- c("Author","Country")
    # build choices = Domain::Item (Top20 each)
    ch <- character(0)
    for (d in doms) {
      top <- if (tolower(d) %in% c("mesh", "mesh term", "meshterm")) {
        .get_mesh_top_by_pubmeta_or_list(rv, n = 10)
      } else {
        .get_domain_top20(rv, d, n = 10)
      }
      if (length(top)) {
        ch <- c(ch, paste0(d, "::", top))
      }
    }
    ch <- unique(ch)
    updateSelectizeInput(session, "yc_items", choices = ch, selected = character(0), server = TRUE)
  }, ignoreInit = TRUE)

  output$plt_slope_top2_country <- renderPlot({
    req(rv$pubmeta)
    doms <- input$yc_domains %||% c("Author","Country")
    yc <- compute_combo_year_counts(rv, domains = doms, recent_n_years = input$yc_recent_years %||% 10)

    validate(need(!is.null(yc) && nrow(yc) > 0, "No data for slopegraph (need Year + term lists)."))

    # Optional filtering by selected items (Domain::Item)
    sel <- input$yc_items %||% character(0)
    if (length(sel)) {
      dom_sel  <- sub("::.*$", "", sel)
      item_sel <- sub("^.*?::", "", sel)
      yc <- yc[yc$domain %in% dom_sel & yc$item %in% item_sel, , drop=FALSE]
      validate(need(nrow(yc) > 0, "Selected items have no counts in recent years window."))
    }

    # stats (p-value two decimals) + trend
    st <- summarize_prepost_ttest2(yc)

    # Build group label for plotting
    yc$group_label <- paste0(yc$domain, "::", yc$item)

    df_s <- yc[, c("group_label","Year","Count")]
    names(df_s) <- c("itemlab","Year","Count")

    # Tufte spacing
    df_t <- tufte_sort2(df_s, x="Year", y="Count", group="itemlab", min.space=0.05)

    # keep Year order
    yrs <- sort(unique(as.character(df_s$Year)))
    if (all(grepl("^\\d{4}$", yrs))) yrs <- as.character(sort(as.integer(yrs)))
    df_t$x <- factor(as.character(df_t$x), levels = yrs)

    plot_slopegraph2(
      df_t,
      st = st,
      title = paste0("Top10 combo slopegraph (recent ", input$yc_recent_years %||% 10, " years)"),
      label_digits = 0
    )
  })
  output$plt_slope_top2_author <- renderPlot({
    req(rv$pubmeta)
    yc <- compute_combo_year_counts(rv, domains = c("Author"), recent_n_years = input$yc_recent_years %||% 10)
    validate(need(!is.null(yc) && nrow(yc) > 0, "No Author slopegraph data (need Year + Author metadata)."))

    sel <- input$yc_items %||% character(0)
    sel_author <- sel[grepl("^Author::", sel)]
    if (length(sel_author)) {
      item_sel <- sub("^.*?::", "", sel_author)
      yc <- yc[yc$item %in% item_sel, , drop=FALSE]
      validate(need(nrow(yc) > 0, "Selected Author items have no counts in recent years window."))
    }

    st <- summarize_prepost_ttest2(yc)
    yc$group_label <- paste0(yc$domain, "::", yc$item)
    df_s <- yc[, c("group_label","Year","Count")]
    names(df_s) <- c("itemlab","Year","Count")
    df_t <- tufte_sort2(df_s, x="Year", y="Count", group="itemlab", min.space=0.05)

    yrs <- sort(unique(as.character(df_s$Year)))
    if (all(grepl("^[0-9]{4}$", yrs))) yrs <- as.character(sort(as.integer(yrs)))
    df_t$x <- factor(as.character(df_t$x), levels = yrs)

    plot_slopegraph2(
      df_t,
      st = st,
      title = paste0("Author metadata slopegraph (recent ", input$yc_recent_years %||% 10, " years)"),
      label_digits = 0
    )
  })

  output$tbl_yearcount_top20_author <- DT::renderDT({
    req(rv$pubmeta)
    yc <- compute_combo_year_counts(rv, domains = c("Author"), recent_n_years = input$yc_recent_years %||% 10)
    validate(need(!is.null(yc) && nrow(yc) > 0, "No Author year-count data."))

    sel <- input$yc_items %||% character(0)
    sel_author <- sel[grepl("^Author::", sel)]
    if (length(sel_author)) {
      item_sel <- sub("^.*?::", "", sel_author)
      yc <- yc[yc$item %in% item_sel, , drop=FALSE]
      validate(need(nrow(yc) > 0, "Selected Author items have no counts in recent years window."))
    }
    yc$group_label <- paste0(yc$domain, "::", yc$item)

    st <- summarize_prepost_ttest2(yc)
    if (is.null(st) || !nrow(st)) {
      st <- unique(yc[, c("domain","item")])
      st$mean_pre <- st$mean_post <- st$pval <- NA_real_
      st$pval_fmt <- NA_character_
      st$trend <- "flat"
    }
    st$group_label <- paste0(st$domain, "::", st$item)

    wide <- reshape2::dcast(yc, group_label + domain + item ~ Year, value.var="Count", fill=0)
    out <- merge(st[, c("group_label","mean_pre","mean_post","pval","pval_fmt","trend")], wide, by="group_label", all.x=TRUE)
    out$pval_num <- out$pval
    out$pval <- out$pval_fmt
    out <- out[order(is.na(out$pval_num), out$pval_num, -out$mean_post), , drop=FALSE]

    DT::datatable(out, options=list(pageLength=10, scrollX=TRUE))
  })


  output$plt_slope_top2_mesh <- renderPlot({
    req(rv$pubmeta)
    yc <- compute_combo_year_counts(rv, domains = c("MeSH"), recent_n_years = input$yc_recent_years %||% 10)
    validate(need(!is.null(yc) && nrow(yc) > 0, "No MeSH slopegraph data (need Year + MeSH metadata)."))

    sel <- input$yc_items %||% character(0)
    sel_mesh <- sel[grepl("^MeSH::", sel)]
    if (length(sel_mesh)) {
      item_sel <- sub("^.*?::", "", sel_mesh)
      yc <- yc[yc$item %in% item_sel, , drop=FALSE]
      validate(need(nrow(yc) > 0, "Selected MeSH items have no counts in recent years window."))
    }

    st <- summarize_prepost_ttest2(yc)
    yc$group_label <- paste0(yc$domain, "::", yc$item)
    df_s <- yc[, c("group_label","Year","Count")]
    names(df_s) <- c("itemlab","Year","Count")
    df_t <- tufte_sort2(df_s, x="Year", y="Count", group="itemlab", min.space=0.05)

    yrs <- sort(unique(as.character(df_s$Year)))
    if (all(grepl("^[0-9]{4}$", yrs))) yrs <- as.character(sort(as.integer(yrs)))
    df_t$x <- factor(as.character(df_t$x), levels = yrs)

    plot_slopegraph2(
      df_t,
      st = st,
      title = paste0("MeSH metadata slopegraph (recent ", input$yc_recent_years %||% 10, " years)"),
      label_digits = 0
    )
  })

  output$tbl_yearcount_top20_mesh <- DT::renderDT({
    req(rv$pubmeta)
    yc <- compute_combo_year_counts(rv, domains = c("MeSH"), recent_n_years = input$yc_recent_years %||% 10)
    validate(need(!is.null(yc) && nrow(yc) > 0, "No MeSH year-count data."))

    sel <- input$yc_items %||% character(0)
    sel_mesh <- sel[grepl("^MeSH::", sel)]
    if (length(sel_mesh)) {
      item_sel <- sub("^.*?::", "", sel_mesh)
      yc <- yc[yc$item %in% item_sel, , drop=FALSE]
      validate(need(nrow(yc) > 0, "Selected MeSH items have no counts in recent years window."))
    }
    yc$group_label <- paste0(yc$domain, "::", yc$item)

    st <- summarize_prepost_ttest2(yc)
    if (is.null(st) || !nrow(st)) {
      st <- unique(yc[, c("domain","item")])
      st$mean_pre <- st$mean_post <- st$pval <- NA_real_
      st$pval_fmt <- NA_character_
      st$trend <- "flat"
    }
    st$group_label <- paste0(st$domain, "::", st$item)

    wide <- reshape2::dcast(yc, group_label + domain + item ~ Year, value.var="Count", fill=0)
    out <- merge(st[, c("group_label","mean_pre","mean_post","pval","pval_fmt","trend")], wide, by="group_label", all.x=TRUE)
    out$pval_num <- out$pval
    out$pval <- out$pval_fmt
    out <- out[order(is.na(out$pval_num), out$pval_num, -out$mean_post), , drop=FALSE]

    DT::datatable(out, options=list(pageLength=10, scrollX=TRUE))
  })

  output$tbl_yearcount_top20 <- DT::renderDT({
    req(rv$pubmeta)
    doms <- input$yc_domains %||% c("Author","Country")
    yc <- compute_combo_year_counts(rv, domains = doms, recent_n_years = input$yc_recent_years %||% 10)
    validate(need(!is.null(yc) && nrow(yc) > 0, "No year-count data."))

    sel <- input$yc_items %||% character(0)
    if (length(sel)) {
      dom_sel <- sub("::.*$", "", sel)
      item_sel <- sub("^.*?::", "", sel)
      yc <- yc[yc$domain %in% dom_sel & yc$item %in% item_sel, , drop=FALSE]
      validate(need(nrow(yc) > 0, "Selected items have no counts in recent years window."))
    }
    yc$group_label <- paste0(yc$domain, "::", yc$item)

    st <- summarize_prepost_ttest2(yc)
    if (is.null(st) || !nrow(st)) {
      st <- unique(yc[, c("domain","item")])
      st$mean_pre <- st$mean_post <- st$pval <- NA_real_
      st$trend <- "flat"
    }
    st$group_label <- paste0(st$domain, "::", st$item)

    wide <- reshape2::dcast(yc, group_label + domain + item ~ Year, value.var="Count", fill=0)
    out <- merge(st[, c("group_label","mean_pre","mean_post","pval","pval_fmt","trend")], wide, by="group_label", all.x=TRUE)
    out$pval_num <- out$pval
    out$pval <- out$pval_fmt
    out <- out[order(is.na(out$pval_num), out$pval_num, -out$mean_post), , drop=FALSE]
    DT::datatable(out, options=list(pageLength=10, scrollX=TRUE))
  })


  # ---- Slope tab: app(708).R-style combo entity Top10 over years ----
  .slope_combo_selected_domain <- reactive({
    dom <- input$slope_combo_entity_domain %||% "Author"
    if (isTRUE(input$slope_combo_sync_combo)) {
      cd <- input$combo_domain %||% dom
      cd <- as.character(cd)
      # Combo tab may use Journal/Year; slope year-count needs Journal.
      if (identical(cd, "Journal/Year")) cd <- "Journal"
      if (identical(cd, "Year/Articles")) cd <- "Year"
      if (identical(cd, "Term/Year")) cd <- "MeSH"
      allowed <- c("Author","Journal","Country","State/Province","Institute","Department","MeSH")
      if (cd %in% allowed) dom <- cd
    }
    dom
  })

  .build_slope_combo_entity_counts <- reactive({
    req(rv$pubmeta)
    dom <- .slope_combo_selected_domain()
    recent_n <- input$slope_combo_recent_years %||% 10
    yc <- compute_combo_year_counts(
      rv,
      domains = c(dom),
      recent_n_years = recent_n
    )
    if (is.null(yc) || !is.data.frame(yc) || !nrow(yc)) return(data.frame())

    # Guarantee Top10 by total count in selected entity domain.
    top_items <- yc |>
      dplyr::group_by(domain, item) |>
      dplyr::summarise(total = sum(Count, na.rm = TRUE), .groups = "drop") |>
      dplyr::arrange(dplyr::desc(total), item) |>
      dplyr::slice_head(n = 10)

    yc <- yc[yc$domain %in% top_items$domain & yc$item %in% top_items$item, , drop = FALSE]
    yc
  })

  output$plot_slope_combo_entity_top10 <- renderPlot({
    yc <- .build_slope_combo_entity_counts()
    validate(need(is.data.frame(yc) && nrow(yc) > 0,
                  "No slopegraph data yet. Run analysis first and select an entity domain with Year data."))

    st <- summarize_prepost_ttest2(yc)
    yc$group_label <- paste0(yc$domain, "::", yc$item)

    df_s <- yc[, c("group_label", "Year", "Count")]
    names(df_s) <- c("itemlab", "Year", "Count")

    df_t <- tufte_sort2(df_s, x = "Year", y = "Count", group = "itemlab", min.space = 0.05)

    yrs <- sort(unique(as.character(df_s$Year)))
    if (all(grepl("^[0-9]{4}$", yrs))) yrs <- as.character(sort(as.integer(yrs)))
    df_t$x <- factor(as.character(df_t$x), levels = yrs)

    print(plot_slopegraph2(
      df_t,
      st = st,
      title = paste0(.slope_combo_selected_domain(), " Top10 slopegraph over recent ",
                     input$slope_combo_recent_years %||% 10, " years"),
      label_digits = 0
    ))
  })

  output$tbl_slope_combo_entity_top10 <- DT::renderDT({
    yc <- .build_slope_combo_entity_counts()
    validate(need(is.data.frame(yc) && nrow(yc) > 0, "No year-count data."))

    yc$group_label <- paste0(yc$domain, "::", yc$item)
    st <- summarize_prepost_ttest2(yc)
    if (is.null(st) || !nrow(st)) {
      st <- unique(yc[, c("domain", "item")])
      st$mean_pre <- st$mean_post <- st$pval <- NA_real_
      st$pval_fmt <- NA_character_
      st$trend <- "flat"
    }
    st$group_label <- paste0(st$domain, "::", st$item)

    wide <- reshape2::dcast(yc, group_label + domain + item ~ Year, value.var = "Count", fill = 0)
    out <- merge(
      st[, c("group_label", "mean_pre", "mean_post", "pval", "pval_fmt", "trend")],
      wide,
      by = "group_label",
      all.x = TRUE
    )
    out$pval_num <- out$pval
    out$pval <- out$pval_fmt
    out <- out[order(is.na(out$pval_num), out$pval_num, -out$mean_post), , drop = FALSE]
    DT::datatable(out, options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
  })



output$vn_journal_year <- renderVisNetwork({
      req(rv$pubmeta)
      df <- rv$pubmeta
      if (!all(c("Year","Journal") %in% names(df))) {
        # nothing to show
        return(visNetwork::visNetwork(data.frame(id=1,label="No Journal/Year data"), data.frame()))
      }
      # Top 20 by publication count
      agg <- df %>%
        dplyr::filter(!is.na(Year), !is.na(Journal), nzchar(Year), nzchar(Journal)) %>%
        dplyr::mutate(Year = as.character(Year), Journal = as.character(Journal)) %>%
        dplyr::count(Journal, Year, name="n")
      if (nrow(agg) == 0) return(visNetwork::visNetwork(data.frame(id=1,label="No Journal/Year data"), data.frame()))
      topJ <- agg %>% dplyr::group_by(Journal) %>% dplyr::summarise(n=sum(n), .groups="drop") %>% dplyr::arrange(dplyr::desc(n)) %>% dplyr::slice_head(n=20)
      agg2 <- agg %>% dplyr::semi_join(topJ, by="Journal")
      # bipartite graph: J:xxx <-> Y:yyyy
      nodesJ <- data.frame(id=paste0("J:", topJ$Journal), label=topJ$Journal, group="Journal", value=topJ$n, stringsAsFactors = FALSE)
      nodesY <- data.frame(id=paste0("Y:", sort(unique(agg2$Year))), label=sort(unique(agg2$Year)), group="Year", value=1, stringsAsFactors = FALSE)
      nodes <- rbind(nodesJ, nodesY)
      edges <- data.frame(from=paste0("J:", agg2$Journal), to=paste0("Y:", agg2$Year), value=agg2$n, stringsAsFactors = FALSE)

      visNetwork::visNetwork(nodes, edges, height="520px") %>%
        visNetwork::visEdges(smooth=FALSE) %>%
        visNetwork::visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
        visNetwork::visLegend()
    })


  # --- SSplot (Author) --- (Author) ---
  
  output$kano_ss_astar <- renderPlot({
    req(rv$author_nodes)
    nd <- rv$author_nodes
    if (!is.data.frame(nd) || nrow(nd)==0) { plot.new(); text(0.5,0.5,"No nodes"); return() }
    # compute if missing
    if (!("a_star1" %in% names(nd)) && ("a_i" %in% names(nd))) nd$a_star1 <- 1/(1+as.numeric(nd$a_i))
    if (!("ssi" %in% names(nd))) {
      plot.new(); text(0.5,0.5,"ssi not available (run FLCA)"); return()
    }
    ok <- is.finite(nd$ssi) & is.finite(nd$a_star1)
    if (!any(ok)) { plot.new(); text(0.5,0.5,"Need ssi and a*"); return() }
    plot(nd$a_star1[ok], nd$ssi[ok], xlab="a* (1/(1+a_i))", ylab="SS (silhouette width)", main="Kano: SS vs a*")
    lab <- if ("id" %in% names(nd)) nd$id else if ("name" %in% names(nd)) nd$name else rep("", nrow(nd))
    text(nd$a_star1[ok], nd$ssi[ok], labels=lab[ok], pos=4, cex=0.7)
  }, height = 560)

output$ssplot_panel <- renderPlot({

  dom <- input$ss_domain %||% "Author"

  get_nodes <- function(dom){
    switch(dom,
      "Author"          = rv$author_nodes,
      "Country"         = rv$country_nodes,
      "State/Province"  = rv$stateprov_nodes,
      "Institute"       = rv$inst_nodes,
      "Department"      = rv$dept_nodes,
      "MeSH"            = rv$mesh_nodes,
      "Journal/Year"    = (rv$jy_nodes %||% rv$journal_nodes),
      NULL
    )
  }
  nd <- get_nodes(dom)

  if (is.null(nd) || !is.data.frame(nd) || nrow(nd) == 0) {
    plot.new()
    text(0.5, 0.5, paste0("SSplot: domain '", dom, "' not ready"), cex = 1.05)
    return(invisible(NULL))
  }

  sil_df <- nd
  if (!("name" %in% names(sil_df))) {
    if ("id" %in% names(sil_df)) sil_df$name <- as.character(sil_df$id) else sil_df$name <- as.character(sil_df[[1]])
  }
  if (!("carac" %in% names(sil_df))) sil_df$carac <- "C1"

  # RenderSSPlus expects sil_width
  if (!("sil_width" %in% names(sil_df))) {
    if ("ssi" %in% names(sil_df)) sil_df$sil_width <- as.numeric(sil_df$ssi)
    else if ("SSi" %in% names(sil_df)) sil_df$sil_width <- as.numeric(sil_df$SSi)
    else sil_df$sil_width <- 0
  }

  # ---- NEVER-NA SS width (fix SSplot NA) ----
  sil_df$sil_width <- suppressWarnings(as.numeric(sil_df$sil_width))
  sil_df$sil_width[!is.finite(sil_df$sil_width)] <- 0

  # keep compatibility columns for downstream renderers
  if ("ssi" %in% names(sil_df)) {
    sil_df$ssi <- suppressWarnings(as.numeric(sil_df$ssi))
    sil_df$ssi[!is.finite(sil_df$ssi)] <- 0
  }
  if ("SSi" %in% names(sil_df)) {
    sil_df$SSi <- suppressWarnings(as.numeric(sil_df$SSi))
    sil_df$SSi[!is.finite(sil_df$SSi)] <- 0
  }
  if ("a_star1" %in% names(sil_df)) {
    sil_df$a_star1 <- suppressWarnings(as.numeric(sil_df$a_star1))
    sil_df$a_star1[!is.finite(sil_df$a_star1)] <- 0
  }
  if ("a_star" %in% names(sil_df)) {
    sil_df$a_star <- suppressWarnings(as.numeric(sil_df$a_star))
    sil_df$a_star[!is.finite(sil_df$a_star)] <- 0
  }

# --- normalize cluster labels (carac) to numeric cluster ids; keep Cluster labels as C# in results ---
.normClusterNum <- function(x){
  if (is.null(x)) return(rep(NA_integer_, nrow(sil_df)))
  x <- as.character(x)
  suppressWarnings(as.integer(gsub("[^0-9]", "", toupper(x))))
}
sil_df$carac <- .normClusterNum(sil_df$carac)
sil_df$carac[is.na(sil_df$carac)] <- 1L

# --- ensure neighborC exists (needed by SSplot panel tooltip/table) ---
if (!("neighborC" %in% names(sil_df)) || all(is.na(sil_df$neighborC))) {
  if ("nn_cluster" %in% names(sil_df)) sil_df$neighborC <- sil_df$nn_cluster else sil_df$neighborC <- sil_df$carac
} else {
  if (!("neighborC" %in% names(sil_df))) sil_df$neighborC <- sil_df$carac
}
sil_df$neighborC <- .normClusterNum(sil_df$neighborC)
sil_df$neighborC[is.na(sil_df$neighborC)] <- sil_df$carac[is.na(sil_df$neighborC)]

# Force SSplot inputs to app-defined values: value=document count, value2=sum(edge)
if (is.data.frame(nd) && nrow(nd)) {
  nd_name_clean  <- trimws(sub("#[0-9]+$", "", sub("#C[0-9]+$", "", as.character(nd$name))))
  sil_name_clean <- trimws(sub("#[0-9]+$", "", sub("#C[0-9]+$", "", as.character(sil_df$name))))
  ixv <- match(sil_name_clean, nd_name_clean)
  if ("value" %in% names(nd))  sil_df$value  <- suppressWarnings(as.numeric(nd$value[ixv]))
  if ("value2" %in% names(nd)) sil_df$value2 <- suppressWarnings(as.numeric(nd$value2[ixv]))
  if ("carac" %in% names(nd))  sil_df$carac  <- suppressWarnings(as.integer(gsub("[^0-9]", "", as.character(nd$carac[ixv]))))
}
sil_df$value[!is.finite(sil_df$value) | sil_df$value < 1] <- 1
sil_df$value2[!is.finite(sil_df$value2)] <- 0
if ("value" %in% names(nd)) {
  nd$value <- suppressWarnings(as.numeric(nd$value))
  nd$value[!is.finite(nd$value) | nd$value < 1] <- 1
}
if ("value2" %in% names(nd)) {
  nd$value2 <- suppressWarnings(as.numeric(nd$value2))
  nd$value2[!is.finite(nd$value2)] <- 0
}

# Force SSplot to show only FLCA-MA-SIL Top20 nodes.
# nodes0 keeps the full domain for footer totals; sil_df/nodes used by render_panel are Top20.
nd_full <- nd
sil_df <- .safe_top20_nodes(sil_df, 20)
sil_df$value[!is.finite(sil_df$value) | sil_df$value < 1] <- 1
sil_df$value2[!is.finite(sil_df$value2)] <- 0
nd <- sil_df

# clusterwise summary for SSplot; add Q summaries from appwos.R
clv <- unique(stats::na.omit(sil_df$carac))
results <- do.call(rbind, lapply(clv, function(cc){
  sub <- sil_df[sil_df$carac == cc, , drop = FALSE]
  qwc <- suppressWarnings(as.numeric(sub$Qw_cluster %||% sub$Qw %||% NA_real_))
  quc <- suppressWarnings(as.numeric(sub$Qu_cluster %||% sub$Qu %||% NA_real_))
  qc  <- suppressWarnings(as.numeric(sub$Q_cluster  %||% sub$Q   %||% NA_real_))
  data.frame(
    Cluster = paste0("C", cc),
    SS = mean(sub$sil_width, na.rm = TRUE),
    Qw = if (all(!is.finite(qwc))) 0 else mean(qwc[is.finite(qwc)], na.rm = TRUE),
    Qu = if (all(!is.finite(quc))) 0 else mean(quc[is.finite(quc)], na.rm = TRUE),
    Q  = if (all(!is.finite(qc)))  0 else mean(qc[is.finite(qc)],  na.rm = TRUE),
    n = nrow(sub),
    stringsAsFactors = FALSE
  )
}))
results$SS[!is.finite(results$SS)] <- 0
results$Qw[!is.finite(results$Qw)] <- 0
results$Qu[!is.finite(results$Qu)] <- 0
results$Q[!is.finite(results$Q)] <- 0
results$SS_total <- mean(sil_df$sil_width, na.rm = TRUE)
results$Qw_total <- sum(results$Qw, na.rm = TRUE)
results$Qu_total <- sum(results$Qu, na.rm = TRUE)
results$Q_total  <- suppressWarnings(mean(as.numeric(sil_df$Q_total %||% sil_df$Q_cluster %||% sil_df$Q %||% results$Q), na.rm = TRUE))
results$SS_total[!is.finite(results$SS_total)] <- 0
results$Qw_total[!is.finite(results$Qw_total)] <- 0
results$Qu_total[!is.finite(results$Qu_total)] <- 0
results$Q_total[!is.finite(results$Q_total)] <- 0
results <- rbind(
  results,
  data.frame(
    Cluster = "OVERALL",
    SS = results$SS_total[1],
    Qw = results$Qw_total[1],
    Qu = results$Qu_total[1],
    Q = results$Q_total[1],
    n = nrow(sil_df),
    SS_total = results$SS_total[1],
    Qw_total = results$Qw_total[1],
    Qu_total = results$Qu_total[1],
    Q_total = results$Q_total[1],
    stringsAsFactors = FALSE
  )
)

  # prefer Render SS Plus if available
  if (exists("render_panel")) {
    tryCatch(
      render_panel(
        sil_df = sil_df,
        nodes0 = nd_full,
        results = results,
        nodes = nd,
        top_n = 20,
        font_scale = (input$ss_font_scale %||% 1.3)
      ),
      error = function(e){
        plot.new()
        text(0.5, 0.5, paste("SSplot render error:", e$message), cex = 1.0)
      }
    )
  } else if (exists(".render_panel_ss")) {
    tryCatch(
      .render_panel_ss(
        sil_df = sil_df,
        results = results,
        nodes0 = nd_full,
        nodes = nd,
        top_n = 20,
        font_scale = (input$ss_font_scale %||% 1.3)
      ),
      error = function(e){
        plot.new()
        text(0.5, 0.5, paste("SSplot render error:", e$message), cex = 1.0)
      }
    )
  } else {
    plot.new()
    text(0.5, 0.5, "SSplot: renderSSplot.R not loaded", cex = 1.05)
  }

}, height = 820)

 .aac_from_values <- function(v){
  v <- as.numeric(v)
  v <- v[is.finite(v) & !is.na(v)]
  if (!length(v)) return(NA_real_)
  v <- sort(v, decreasing = TRUE)
  if (length(v) == 1) return(0.5)
  if (length(v) == 2) {
    if (v[2] == 0) return(NA_real_)
    r12 <- v[1] / v[2]
    g <- r12
    return(g / (1 + g))
  }
  # length >= 3
  if (v[2] == 0 || v[3] == 0) return(NA_real_)
  r12 <- v[1] / v[2]
  r23 <- v[2] / v[3]
  if (!is.finite(r12) || !is.finite(r23) || r23 == 0) return(NA_real_)
  g <- r12 / r23
  g / (1 + g)
}
 .top10_name_value <- function(df, name_col = NULL, value_col = NULL){
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) return(data.frame(name=character(), value=numeric()))
  nms <- names(df)
  if (is.null(name_col)) {
    name_col <- intersect(c("name","label","term","Entity","entity","Country","Journal","Year"), nms)[1]
  }
  if (is.null(value_col)) {
    value_col <- intersect(c("value","count","n","Freq","freq","Value"), nms)[1]
  }
  if (is.na(name_col) || is.na(value_col)) return(data.frame(name=character(), value=numeric()))
  out <- data.frame(
    name  = as.character(df[[name_col]]),
    value = suppressWarnings(as.numeric(df[[value_col]])),
    stringsAsFactors = FALSE
  )
  out <- out[is.finite(out$value) & !is.na(out$value) & nzchar(out$name), , drop=FALSE]
  out <- out[order(out$value, decreasing = TRUE), , drop=FALSE]
  utils::head(out, 10)
}
 .render_sum_dt <- function(df){
  DT::datatable(
    df,
    rownames = FALSE,
    options = list(dom='t', pageLength = 10, ordering = FALSE),
    selection = "none"
  )
}
 .year_table <- reactive({
  y <- rv$article_years
  if (is.null(y) || !length(y)) return(data.frame(name=character(), value=numeric()))
  y <- y[!is.na(y) & nzchar(as.character(y))]
  tb <- sort(table(y), decreasing = TRUE)
  data.frame(name = names(tb), value = as.numeric(tb), stringsAsFactors = FALSE)[1:min(10, length(tb)), , drop=FALSE]
})
 # reactive summaries# ---- Presence counts helper (per-article unique, never > n_pubmed) ----
.top10_presence_from_list <- function(lst){
  if (is.null(lst) || length(lst)==0) return(data.frame(name=character(), value=numeric(), stringsAsFactors=FALSE))
  # lst: list of character vectors, one per article
  per_article <- lapply(lst, function(v){
    v <- as.character(v)
    v <- trimws(v)
    v <- v[!is.na(v) & nzchar(v)]
    unique(v)
  })
  all_terms <- unlist(per_article, use.names=FALSE)
  if (!length(all_terms)) return(data.frame(name=character(), value=numeric(), stringsAsFactors=FALSE))
  tb <- sort(table(all_terms), decreasing=TRUE)
  data.frame(name=names(tb), value=as.numeric(tb), stringsAsFactors=FALSE)[1:min(10, length(tb)), , drop=FALSE]
}


sum_author    <- reactive(.top10_name_value(rv$author_nodes))
sum_country   <- reactive(.top10_presence_from_list(rv$countries_list))
sum_stateprov <- reactive(.top10_presence_from_list(rv$stateprov_list))
sum_inst      <- reactive(.top10_presence_from_list(rv$inst_list))
sum_dept      <- reactive(.top10_presence_from_list(rv$dept_list))
sum_mesh      <- reactive(.top10_presence_from_list(rv$mesh_list))
sum_journal   <- reactive({
  pm <- rv$pubmeta
  if (!is.data.frame(pm) || !nrow(pm) || !("Journal" %in% names(pm))) {
    return(data.frame(name=character(), value=numeric(), stringsAsFactors = FALSE))
  }
  j <- trimws(as.character(pm$Journal))
  j <- j[!is.na(j) & nzchar(j)]
  tb <- sort(table(j), decreasing = TRUE)
  if (!length(tb)) return(data.frame(name=character(), value=numeric(), stringsAsFactors = FALSE))
  utils::head(data.frame(name=names(tb), value=as.numeric(tb), stringsAsFactors = FALSE), 10)
})
sum_year      <- reactive(.year_table())
 # AAC text outputs
output$aac_author    <- renderText({ sprintf("AAC=%.2f", .aac_from_values(sum_author()$value)) })
output$aac_country   <- renderText({ sprintf("AAC=%.2f", .aac_from_values(sum_country()$value)) })
output$aac_stateprov <- renderText({ sprintf("AAC=%.2f", .aac_from_values(sum_stateprov()$value)) })
output$aac_inst      <- renderText({ sprintf("AAC=%.2f", .aac_from_values(sum_inst()$value)) })
output$aac_dept      <- renderText({ sprintf("AAC=%.2f", .aac_from_values(sum_dept()$value)) })
output$aac_mesh      <- renderText({ sprintf("AAC=%.2f", .aac_from_values(sum_mesh()$value)) })
output$aac_journal   <- renderText({ sprintf("AAC=%.2f", .aac_from_values(sum_journal()$value)) })
output$aac_year      <- renderText({ sprintf("AAC=%.2f", .aac_from_values(sum_year()$value)) })
 # DT outputs
output$sum_author    <- DT::renderDT(.render_sum_dt(sum_author()))
output$sum_country   <- DT::renderDT(.render_sum_dt(sum_country()))
output$sum_stateprov <- DT::renderDT(.render_sum_dt(sum_stateprov()))
output$sum_inst      <- DT::renderDT(.render_sum_dt(sum_inst()))
output$sum_dept      <- DT::renderDT(.render_sum_dt(sum_dept()))
output$sum_mesh      <- DT::renderDT(.render_sum_dt(sum_mesh()))
output$sum_journal   <- DT::renderDT(.render_sum_dt(sum_journal()))
output$sum_year      <- DT::renderDT(.render_sum_dt(sum_year()))

# ---- Summary HTML report (preview + download) ----
make_summary_html <- reactive({
# build a WoS-style summary HTML using the current Top10 + AAC blocks
blocks <- list(
  list(dom="Country",       aac=.aac_from_values(sum_country()$value),    df=sum_country()),
  list(dom="Institute",     aac=.aac_from_values(sum_inst()$value),       df=sum_inst()),
  list(dom="Department",    aac=.aac_from_values(sum_dept()$value),       df=sum_dept()),
  list(dom="Author",        aac=.aac_from_values(sum_author()$value),     df=sum_author()),
  list(dom="Journal",       aac=.aac_from_values(sum_journal()$value),    df=sum_journal()),
  list(dom="Year",          aac=.aac_from_values(sum_year()$value),       df=sum_year()),
  list(dom="State/Province",aac=.aac_from_values(sum_stateprov()$value),  df=sum_stateprov()),
  list(dom="MeSH Term",     aac=.aac_from_values(sum_mesh()$value),       df=sum_mesh())
)
 domain_block <- function(dom, aac, df){
  if (is.null(df) || !is.data.frame(df) || !nrow(df)) return(NULL)
  dd <- df
  dd$name <- as.character(dd$name)
  dd$value <- suppressWarnings(as.numeric(dd$value))
  dd <- dd[!is.na(dd$name) & nzchar(dd$name), , drop=FALSE]
  dd <- dd[order(dd$value, decreasing = TRUE), , drop=FALSE]
  dd <- utils::head(dd, 5)
   htmltools::tags$div(
    style="border:1px solid #e6e6e6; border-radius:12px; padding:10px; background:#fff;",
    htmltools::tags$div(
      style="display:flex; align-items:flex-end; justify-content:space-between; margin-bottom:6px;",
      htmltools::tags$div(dom, style="color:#d32f2f; font-weight:700; font-size:15px;"),
      htmltools::tags$div(
        style="font-size:12px; color:#111; font-weight:700;",
        sprintf("AAC = %s", ifelse(is.finite(aac), format(round(aac, 2), nsmall=2), "NA"))
      )
    ),
    htmltools::tags$table(
      style="border-collapse:collapse; width:100%; font-size:12px;",
      htmltools::tags$tbody(
        lapply(seq_len(nrow(dd)), function(i){
          htmltools::tags$tr(
            htmltools::tags$td(htmltools::htmlEscape(dd$name[i]),
                               style="padding:2px 4px; text-align:left; white-space:normal; overflow-wrap:break-word; word-break:break-word;"),
            htmltools::tags$td(ifelse(is.finite(dd$value[i]), sprintf("%d", as.integer(dd$value[i])), ""),
                               style="padding:2px 4px; text-align:right; width:60px;")
          )
        })
      )
    )
  )
}
 grid <- htmltools::tags$div(
  style="display:grid; grid-template-columns: 1fr 1fr; gap:12px;",
  lapply(blocks, function(b) domain_block(b$dom, b$aac, b$df))
)
 page <- htmltools::tags$html(
  htmltools::tags$head(
    tags$script(HTML("(function(){var tries=0;var t=setInterval(function(){tries++;try{if(window.Shiny&&Shiny.addCustomMessageHandler){Shiny.addCustomMessageHandler('jsClick',function(id){try{var el=document.getElementById(id); if(el) el.click();}catch(e){}});clearInterval(t);}}catch(e){} if(tries>50){clearInterval(t);} },100);})();")),
    tags$script(HTML('Shiny.addCustomMessageHandler("jsClick", function(id){var el=document.getElementById(id); if(el){el.click();}});')),
    htmltools::tags$meta(charset="utf-8"),
    htmltools::tags$title("Summary Report")
  ),
  htmltools::tags$body(
    style="font-family: Arial, sans-serif; margin: 12px;",
    htmltools::tags$h3("Summary Report (Top5 + AAC)"),
    htmltools::tags$p(sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))),
    grid
  )
)
htmltools::renderTags(page)$html
})

output$summary_html_preview <- renderUI({
# Always show an HTML report preview (same "full page" behavior as appwos.R)
# Use iframe(srcdoc=...) to reliably render full HTML (<html>/<head>/<style> etc.)
html <- NULL
msg <- NULL
 # detect whether any summary domain has data
has_any <- FALSE
for (df in list(sum_country(), sum_inst(), sum_dept(), sum_author(), sum_journal(), sum_year(), sum_stateprov(), sum_mesh())) {
  if (is.data.frame(df) && nrow(df) > 0) { has_any <- TRUE; break }
}
 if (!isTRUE(has_any)) {
  # still render a valid HTML page so the iframe isn't blank
  html <- paste0(
    "<!doctype html><html><head><meta charset='utf-8'><title>Summary Report</title></head><body style='font-family:Arial,sans-serif;margin:12px;'>",
    "<h3>Summary Report (Top5 + AAC)</h3>",
    "<p style='color:#666'>No summary yet. Please run at least one domain/network first.</p>",
    "</body></html>"
  )
} else {
  html <- make_summary_html()
}
 tags$iframe(
  style="width:100%; height:820px; border:1px solid #ddd; border-radius:12px; background:#fff;",
  srcdoc = html
)
})
output$dl_summary_html <- downloadHandler(
filename = function(){ sprintf("summary_%s.html", format(Sys.time(), "%Y%m%d_%H%M%S")) },
contentType = "text/html",
content  = function(file){
  writeLines(make_summary_html(), con = file, useBytes = TRUE)
}
)
 # PNG dashboard download
output$download_summary_png <- downloadHandler(
  filename = function(){ sprintf("summary_%s.png", format(Sys.time(), "%Y%m%d_%H%M%S")) },
  contentType = "image/png",
  content = function(file){
    if (!requireNamespace("grid", quietly = TRUE)) stop("grid package missing")
    if (!requireNamespace("gridExtra", quietly = TRUE)) stop("gridExtra package missing; install.packages('gridExtra')")
     make_panel <- function(title, aac_txt, df){
      grob_tbl <- gridExtra::tableGrob(df, rows = NULL,
                                      theme = gridExtra::ttheme_minimal(
                                        base_size = 9,
                                        padding = grid::unit(c(2, 2), "mm")
                                      ))
      header <- grid::textGrob(
        sprintf("%s  (%s)", title, aac_txt),
        x = 0, just = "left",
        gp = grid::gpar(fontface = "bold", fontsize = 12)
      )
      gridExtra::arrangeGrob(header, grob_tbl, ncol = 1, heights = grid::unit.c(grid::unit(6, "mm"), grid::unit(1, "null")))
    }
     # collect panels in desired order
    panels <- list(
      make_panel("Country",    sprintf("AAC=%.2f", .aac_from_values(sum_country()$value)),   sum_country()),
      make_panel("Institute",  sprintf("AAC=%.2f", .aac_from_values(sum_inst()$value)),      sum_inst()),
      make_panel("Department", sprintf("AAC=%.2f", .aac_from_values(sum_dept()$value)),      sum_dept()),
      make_panel("Author",     sprintf("AAC=%.2f", .aac_from_values(sum_author()$value)),    sum_author()),
      make_panel("Journal",    sprintf("AAC=%.2f", .aac_from_values(sum_journal()$value)),   sum_journal()),
      make_panel("Year",       sprintf("AAC=%.2f", .aac_from_values(sum_year()$value)),      sum_year()),
      make_panel("State/Province", sprintf("AAC=%.2f", .aac_from_values(sum_stateprov()$value)), sum_stateprov()),
      make_panel("MeSH Term",  sprintf("AAC=%.2f", .aac_from_values(sum_mesh()$value)),      sum_mesh())
    )
     grDevices::png(file, width = 1400, height = 1800, res = 150)
    on.exit(grDevices::dev.off(), add = TRUE)
    grid::grid.newpage()
    # 4 rows x 2 cols dashboard
    gridExtra::grid.arrange(grobs = panels, ncol = 2)
  }
)

# =========================
# YEAR-FREQUENCY SLOPEGRAPH (COMBO DOMAINS)
# =========================
.get_domain_term_list <- function(rv, domain){
  d <- tolower(domain)
  # Author domain is restricted to first/last authors for consistency with Author network.
  if (d %in% c("author","authors","author_byline","author(byline)","author (byline)")) return(rv$author_lists_fl %||% rv$author_lists)
  # Author 1st+Last
  if (grepl("1st", d) || grepl("first", d) || grepl("last", d) || d %in% c("author_fl","author_first_last","author(1st+last)","author (1st+last)")) {
    return(rv$author_lists_fl %||% rv$author_lists_all %||% rv$author_lists)
  }
  if (d %in% c("journal")) {
    if (is.data.frame(rv$pubmeta) && "Journal" %in% names(rv$pubmeta)) {
      return(as.list(as.character(rv$pubmeta$Journal)))
    }
    return(NULL)
  }
  if (d %in% c("country")) return(rv$countries_list)
  if (d %in% c("state/province","stateprov","state")) return(rv$stateprov_list)
  if (d %in% c("institute","inst")) return(rv$inst_list)
  if (d %in% c("department","dept")) return(rv$dept_list)
  if (d %in% c("mesh","mesh term","meshterm")) return(rv$mesh_list)
  return(NULL)
}

.get_domain_top20 <- function(rv, domain, n=20){
  d <- tolower(domain)
  take_names <- function(df){
    if (is.data.frame(df) && "name" %in% names(df)) {
      x <- as.character(df$name)
      x <- x[!is.na(x) & nzchar(x)]
      x <- x[.is_single_metadata_item(x)]
      return(utils::head(unique(x), n))
    }
    character(0)
  }
  if (d=="author" || grepl("author", d)) return(take_names(rv$author_nodes))
  if (d=="journal") {
    # prefer JY nodes but keep only journals (exclude pure years)
    if (is.data.frame(rv$jy_nodes) && "name" %in% names(rv$jy_nodes)) {
      x <- as.character(rv$jy_nodes$name)
      x <- x[nzchar(x)]
      # remove year-like tokens
      x <- x[!grepl("^\\d{4}$", x)]
      return(utils::head(unique(x), n))
    }
    if (is.data.frame(rv$pubmeta) && "Journal" %in% names(rv$pubmeta)) {
      x <- sort(table(rv$pubmeta$Journal), decreasing=TRUE)
      xn <- names(x)
      xn <- xn[.is_single_metadata_item(xn)]
      return(utils::head(xn, n))
    }
    return(character(0))
  }
  if (d=="country") return(take_names(rv$country_nodes))
  if (d %in% c("state/province","stateprov","state")) return(take_names(rv$stateprov_nodes))
  if (d %in% c("institute","inst")) return(take_names(rv$inst_nodes))
  if (d %in% c("department","dept")) return(take_names(rv$dept_nodes))
  if (d %in% c("mesh","mesh term","meshterm")) return(take_names(rv$mesh_nodes))
  character(0)
}

# split terms helper (vectorized-friendly)
.split_terms2 <- function(x){
  if (is.null(x)) return(character(0))
  x <- as.character(x)
  x <- x[nzchar(trimws(x))]
  if (!length(x)) return(character(0))
  if (length(x)==1L) {
    s <- x[[1]]
    if (grepl("[;|]", s, perl=TRUE)) x <- unlist(strsplit(s, "[;|]", perl=TRUE), use.names=FALSE)
  }
  x <- trimws(x)
  x <- x[nzchar(x)]
  unique(x)
}

.is_single_metadata_item <- function(x){
  x <- as.character(x)
  x <- trimws(x)
  nzchar(x) & !grepl("[;|]", x, perl=TRUE)
}

.norm_slope_term <- function(x, domain=NULL){
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("\\s+", " ", x)
  x <- gsub("^\\*+|\\*+$", "", x)
  x <- tolower(x)
  d <- tolower(domain %||% "")
  if (d %in% c("mesh", "mesh term", "meshterm")) {
    x <- gsub("\\s*/\\s*.*$", "", x)
  }
  x
}


# ---- MeSH-safe helpers for year-count / slopegraph ----
# MeSH headings may contain commas (e.g., "Ethics, Professional"),
# so MeSH must be split only by semicolon or pipe, never by comma.
.split_mesh_terms_safe <- function(x){
  if (is.null(x)) return(character(0))
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(trimws(x))]
  if (!length(x)) return(character(0))
  unlist(strsplit(x, "\\s*[;|]\\s*", perl = TRUE), use.names = FALSE) |>
    trimws() |>
    unique()
}

.get_mesh_top_by_pubmeta_or_list <- function(rv, n = 10){
  terms <- character(0)
  if (is.data.frame(rv$pubmeta)) {
    cn <- intersect(c("MeSH", "Mesh", "MESH", "mesh"), names(rv$pubmeta))
    if (length(cn)) terms <- .split_mesh_terms_safe(rv$pubmeta[[cn[1]]])
  }
  if (!length(terms) && !is.null(rv$mesh_list)) {
    terms <- unique(unlist(lapply(rv$mesh_list, .split_mesh_terms_safe), use.names = FALSE))
  }
  terms <- trimws(terms)
  terms <- terms[!is.na(terms) & nzchar(terms)]
  if (!length(terms)) return(character(0))
  tb <- sort(table(terms), decreasing = TRUE)
  utils::head(names(tb), min(n, length(tb)))
}

compute_combo_year_counts <- function(rv, domains, items_keep=NULL, years_keep=NULL, recent_n_years=10){
  # Returns long table: domain,item,Year,Count
  # Robust: derives term-lists from rv$*_list first, then falls back to rv$pubmeta columns.
  # Completes missing years with Count=0 for each (domain,item) so slopegraph labels are stable.

  # years vector
  yrs <- rv$article_years
  if (is.null(yrs) && is.data.frame(rv$pubmeta) && "Year" %in% names(rv$pubmeta)) yrs <- rv$pubmeta$Year
  yrs <- suppressWarnings(as.integer(as.character(yrs)))
  ok_year <- is.finite(yrs)
  if (!any(ok_year)) return(NULL)

  years_all <- sort(unique(yrs[ok_year]))
  years_all <- years_all[is.finite(years_all)]
  if (!length(years_all)) return(NULL)

  if (is.null(years_keep)) {
    years_keep <- utils::tail(years_all, min(recent_n_years, length(years_all)))
  } else {
    years_keep <- suppressWarnings(as.integer(years_keep))
    years_keep <- years_keep[is.finite(years_keep)]
  }
  if (length(years_keep) < 2) return(NULL)

  # helper to get term list (vector/list length ~ n articles)
  get_terms_fallback <- function(dom){
    tl <- .get_domain_term_list(rv, dom)
    if (!is.null(tl)) return(tl)

    # fallback to pubmeta columns
    pm <- rv$pubmeta
    if (!is.data.frame(pm)) return(NULL)

    d <- tolower(dom)
    # possible column name candidates
    cand <- switch(
      d,
      "author" = c("Author","Authors","AU","author","authors"),
      "country" = c("Country","Countries","country"),
      "state/province" = c("State/Province","State","Province","state","province","state_province","stateprov"),
      "stateprov" = c("State/Province","State","Province","state","province","state_province","stateprov"),
      "institute" = c("Institute","Inst","institute","inst"),
      "department" = c("Department","Dept","department","dept"),
      "mesh" = c("MeSH","Mesh","MESH","mesh"),
      "journal" = c("Journal","journal"),
      c()
    )
    cn <- intersect(cand, names(pm))
    if (!length(cn)) return(NULL)
    as.list(as.character(pm[[cn[1]]]))
  }

  out <- list()
  for (dom in domains) {
    tl <- get_terms_fallback(dom)

    # For MeSH, prefer article-level rv$pubmeta$MeSH when available.
    # This keeps each MeSH string aligned with the article year and prevents
    # zero-count slopegraphs caused by mesh_nodes/mesh_list mismatch.
    if (tolower(dom) %in% c("mesh", "mesh term", "meshterm") && is.data.frame(rv$pubmeta)) {
      cn_mesh <- intersect(c("MeSH", "Mesh", "MESH", "mesh"), names(rv$pubmeta))
      if (length(cn_mesh)) tl <- as.list(as.character(rv$pubmeta[[cn_mesh[1]]]))
    }
    if (is.null(tl)) next

    # align lengths with years
    n <- min(length(tl), length(yrs))
    tl <- tl[seq_len(n)]
    yv <- yrs[seq_len(n)]

    # choose top20 items. For MeSH, select from article-level MeSH terms,
    # not from rv$mesh_nodes alone, because slope counts must match Year by article.
    top20 <- if (tolower(dom) %in% c("mesh", "mesh term", "meshterm")) {
      .get_mesh_top_by_pubmeta_or_list(rv, n = 10)
    } else {
      .get_domain_top20(rv, dom, n = 10)
    }
    if (!length(top20)) {
      splitter <- if (tolower(dom) %in% c("mesh", "mesh term", "meshterm")) .split_mesh_terms_safe else .split_terms2
      top20 <- utils::head(unique(unlist(lapply(tl, splitter), use.names = FALSE)), 20)
    }
    top20 <- as.character(top20)
    top20 <- top20[!is.na(top20) & nzchar(top20)]
    top20 <- top20[.is_single_metadata_item(top20)]
    top20 <- utils::head(unique(top20), 10)
    if (!length(top20)) next

    if (!is.null(items_keep)) {
      top20 <- intersect(top20, items_keep)
      if (!length(top20)) next
    }

    # count per item per year
    for (it in top20) {
      it_key <- .norm_slope_term(it, dom)
      ct <- tapply(seq_len(n), yv, function(ix){
        sum(vapply(ix, function(i){
          terms <- if (tolower(dom) %in% c("mesh", "mesh term", "meshterm")) {
            .split_mesh_terms_safe(tl[[i]])
          } else {
            .split_terms2(tl[[i]])
          }
          term_keys <- .norm_slope_term(terms, dom)
          any(term_keys == it_key)
        }, logical(1)))
      })

      df <- data.frame(
        domain = dom,
        item   = it,
        Year   = as.integer(names(ct)),
        Count  = as.numeric(ct),
        stringsAsFactors = FALSE
      )
      df <- df[df$Year %in% years_keep, , drop=FALSE]

      # complete missing years with 0
      if (nrow(df) > 0) {
        miss <- setdiff(years_keep, df$Year)
        if (length(miss)) {
          df <- rbind(
            df,
            data.frame(domain=dom, item=it, Year=as.integer(miss), Count=0, stringsAsFactors=FALSE)
          )
        }
        df <- df[order(df$Year), , drop=FALSE]
      } else {
        # if none in window, still keep completed zeros so labels appear
        df <- data.frame(domain=dom, item=it, Year=as.integer(years_keep), Count=0, stringsAsFactors=FALSE)
      }

      out[[length(out)+1]] <- df
    }
  }

  if (!length(out)) return(NULL)
  yc <- do.call(rbind, out)
  yc <- yc[is.finite(yc$Year) & is.finite(yc$Count), , drop=FALSE]
  if (!nrow(yc)) return(NULL)
  yc
}

summarize_prepost_ttest2 <- function(yc){
  # yc: data.frame with columns domain, item, Year, Count (Year can be numeric or character)
  if (is.null(yc) || !is.data.frame(yc) || nrow(yc) == 0) return(NULL)
  if (!all(c("domain","item","Year","Count") %in% names(yc))) return(NULL)

  yrs <- unique(yc$Year)
  # Sort years smartly: if they look like 4-digit years, sort numerically; else lexicographically.
  if (all(grepl("^\\d{4}$", as.character(yrs)))) {
    yrs <- as.character(sort(as.integer(yrs)))
  } else {
    yrs <- sort(as.character(yrs))
  }

  if (length(yrs) < 6) {
    tmp <- unique(yc[, c("domain","item")])
    tmp$mean_pre <- tmp$mean_post <- NA_real_
    tmp$pval <- NA_real_
    tmp$pval_fmt <- NA_character_
    tmp$trend <- "flat"
    rownames(tmp) <- NULL
    return(tmp)
  }

  k <- floor(length(yrs)/2)
  pre_years  <- yrs[seq_len(k)]
  post_years <- yrs[(k+1):length(yrs)]

  out <- yc %>%
    dplyr::mutate(
      domain = as.character(domain),
      item   = as.character(item),
      Year   = as.character(Year),
      Count  = suppressWarnings(as.numeric(Count))
    ) %>%
    dplyr::group_by(domain, item) %>%
    dplyr::summarise(
      mean_pre  = mean(Count[Year %in% pre_years],  na.rm = TRUE),
      mean_post = mean(Count[Year %in% post_years], na.rm = TRUE),
      pval = {
        pre  <- Count[Year %in% pre_years]
        post <- Count[Year %in% post_years]
        pre  <- pre[is.finite(pre)]
        post <- post[is.finite(post)]
        if (length(pre) < 2 || length(post) < 2) NA_real_
        else tryCatch(stats::t.test(pre, post)$p.value, error=function(e) NA_real_)
      },
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      pval_fmt = dplyr::if_else(is.finite(pval), sprintf("%.2f", pval), NA_character_),
      trend = dplyr::case_when(
        is.finite(pval) & pval < 0.05 & mean_post > mean_pre ~ "up",
        is.finite(pval) & pval < 0.05 & mean_post < mean_pre ~ "down",
        TRUE ~ "flat"
      )
    )

  as.data.frame(out)
}

# Tufte-style slope spacing (adapted from jkeirstead slopegraph)
tufte_sort2 <- function(df, x="Year", y="Count", group="item", min.space=0.05) {
  ids <- match(c(x, y, group), names(df))
  df <- df[, ids]
  names(df) <- c("x", "y", "group")
  tmp <- expand.grid(x = unique(df$x), group = unique(df$group))
  tmp <- merge(df, tmp, all.y = TRUE)
  tmp$y <- ifelse(is.na(tmp$y), 0, tmp$y)
  # wide
  wide <- reshape2::dcast(tmp, group ~ x, value.var = "y")
  ord <- order(wide[, 2])
  wide <- wide[ord, ]
  min.space <- min.space * diff(range(wide[, -1]))
  yshift <- numeric(nrow(wide))
  for (i in 2:nrow(wide)) {
    mat <- as.matrix(wide[(i-1):i, -1])
    d.min <- min(diff(mat))
    yshift[i] <- ifelse(d.min < min.space, min.space - d.min, 0)
  }
  wide$yshift <- cumsum(yshift)
  long <- reshape2::melt(wide, id=c("group","yshift"), variable.name="x", value.name="y")
  long$ypos <- long$y + long$yshift
  long
}

plot_slopegraph2 <- function(df, st=NULL, title="Slopegraph (recent years)", label_digits=0){
  # df: output of tufte_sort2 with columns x,y,ypos,group
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) {
    return(ggplot2::ggplot() + ggplot2::theme_void())
  }

  # Build per-group trend from st if available (expects st has domain,item,trend)
  trend_map <- NULL
  if (!is.null(st) && is.data.frame(st) && nrow(st) > 0 && all(c("domain","item","trend") %in% names(st))) {
    trend_map <- setNames(as.character(st$trend), paste0(st$domain, "::", st$item))
  }

  df$group_label <- as.character(df$group)
  df$trend <- if (!is.null(trend_map)) as.character(trend_map[df$group_label]) else NA_character_
  df$trend[is.na(df$trend)] <- "flat"

  df$line_color <- ifelse(df$trend=="up","inc", ifelse(df$trend=="down","dec","flat"))

  # --- Hide zeros: treat Count==0 as missing for plotting (break line, no point/label) ---
  df$has_value <- is.finite(df$y) & (as.numeric(df$y) != 0)
  df$ypos_plot <- ifelse(df$has_value, df$ypos, NA_real_)
  df$label <- ifelse(df$has_value, sprintf(paste0("%.", label_digits, "f"), df$y), "")

  # y-axis labels from the first x positions (only if not too many lines)
  x_levels <- if (is.factor(df$x)) levels(df$x) else sort(unique(as.character(df$x)))
  first_x <- x_levels[1]
  left <- df[df$x == first_x & is.finite(df$ypos), c("group_label","ypos"), drop=FALSE]
  left <- left[order(left$ypos), , drop=FALSE]
  show_names <- nrow(left) > 0 && nrow(left) <= 35

  p <- ggplot2::ggplot(df, ggplot2::aes(x=x, y=ypos_plot, group=group_label)) +
    ggplot2::geom_line(ggplot2::aes(color=line_color), linewidth=1, na.rm=TRUE) +
    ggplot2::geom_point(data=df[df$has_value, , drop=FALSE], color="white", size=5, na.rm=TRUE) +
    ggplot2::geom_text(data=df[df$has_value, , drop=FALSE], ggplot2::aes(label=label), size=3, na.rm=TRUE) +
    ggplot2::scale_color_manual(values=c(inc="red", dec="blue", flat="black"), guide="none") +
    ggplot2::labs(title=title, x=NULL, y=NULL) +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5, face="bold"))

  if (show_names) {
    p <- p + ggplot2::scale_y_continuous(breaks = left$ypos, labels = left$group_label)
  } else {
    p <- p + ggplot2::theme(axis.text.y=ggplot2::element_blank(),
                            axis.ticks.y=ggplot2::element_blank())
  }

  p
}

# ----------------------------
# Helpers: AAC + Author metric augmentation
# ----------------------------
.calc_aac_from_top3 <- function(x){
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  x <- sort(x, decreasing = TRUE)
  if (length(x) < 3) return(NA_real_)
  v1 <- x[1]; v2 <- x[2]; v3 <- x[3]
  if (v2 == 0 || v3 == 0) return(NA_real_)
  r <- (v1 * v3) / (v2 * v2)  # (v1/v2)/(v2/v3)
  if (!is.finite(r)) return(NA_real_)
  r / (1 + r)
}




# ----------------------------
# Top1 (metadomain combo)
# ----------------------------
.top1_nodes_df <- reactive({
  req(input$top1_metadomain)
  dom <- input$top1_metadomain

  # pick nodes by domain (stable rv names)
  df <- switch(dom,
    "Author" = rv$author_nodes,
    "Country" = rv$country_nodes,
    "State/Province" = rv$stateprov_nodes,
    "Institute" = rv$inst_nodes,
    "Department" = rv$dept_nodes,
    "MeSH" = rv$mesh_nodes,
    # Some builds store Journal/Year as a combined domain (rv$jy_nodes)
    "Journal" = rv$jy_nodes %||% rv$journal_nodes %||% rv$jour_nodes,
    "Year" = rv$jy_nodes %||% rv$year_nodes,
    NULL
  )
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) return(NULL)

  # normalize node name
  if (!("name" %in% names(df))) {
    if ("label" %in% names(df)) df$name <- as.character(df$label)
    else if ("id" %in% names(df)) df$name <- as.character(df$id)
  }
  df$name <- as.character(df$name)

  # ensure core columns exist
  if (!("n_pubmed" %in% names(df))) df$n_pubmed <- suppressWarnings(as.numeric(df$`n(Pubmed total)` %||% df$nPubmed %||% NA))
  if (!("single_author" %in% names(df))) {
    # try common variants; if absent, derive from value - value2 when possible
    if ("single" %in% names(df)) df$single_author <- suppressWarnings(as.numeric(df$single))
    else if (all(c("value","value2") %in% names(df))) df$single_author <- pmax(0, suppressWarnings(as.numeric(df$value) - as.numeric(df$value2)))
    else df$single_author <- NA_real_
  }
  if (!("n_byline" %in% names(df))) {
    df$n_byline <- ifelse(is.na(df$n_pubmed), NA_real_, as.numeric(df$n_pubmed))
  }

  # For Author domain: enforce your definition
  # value2 = FA count (out-strength), LA count shown separately
  if (identical(dom, "Author") && !is.null(rv$author_edges_full %||% rv$author_edges) && is.data.frame(rv$author_edges_full %||% rv$author_edges) && nrow(rv$author_edges_full %||% rv$author_edges) > 0) {
    ed <- (rv$author_edges_full %||% rv$author_edges)
    # normalize column names
    if (!("Leader" %in% names(ed)) && "leader" %in% names(ed)) ed$Leader <- ed$leader
    if (!("Follower" %in% names(ed)) && "follower" %in% names(ed)) ed$Follower <- ed$follower
    wcol <- if ("WCD" %in% names(ed)) "WCD" else if ("wcd" %in% names(ed)) "wcd" else NULL
    if (is.null(wcol)) {
      ed$WCD <- 1
      wcol <- "WCD"
    }
    ed$Leader <- as.character(ed$Leader)
    ed$Follower <- as.character(ed$Follower)
    ed[[wcol]] <- suppressWarnings(as.numeric(ed[[wcol]]))
    ed[[wcol]][is.na(ed[[wcol]])] <- 0

    # out/in strength
    fa <- aggregate(ed[[wcol]], by=list(name=ed$Leader), FUN=sum, na.rm=TRUE)
    la <- aggregate(ed[[wcol]], by=list(name=ed$Follower), FUN=sum, na.rm=TRUE)
    names(fa)[2] <- "FA"
    names(la)[2] <- "LA"
    df <- merge(df, fa, by="name", all.x=TRUE)
    df <- merge(df, la, by="name", all.x=TRUE)
    df$FA[is.na(df$FA)] <- 0
    df$LA[is.na(df$LA)] <- 0

    # enforce columns per your definitions:
    # value2 = edge cowork (FA/LA) = FA + LA
    df$value2 <- as.numeric(df$FA + df$LA)

    # value = value2 + self_coword (single-author count)
    df$self_coword <- df$single_author
    df$value <- ifelse(is.na(df$single_author), df$value2, df$value2 + as.numeric(df$single_author))

    # value2_strength = AAC computed from top-3 of value2 within this domain
    df$value2_strength <- (.AAC_INLINE(df$value2))
  } else {
    # Non-author: keep existing if present; derive value where possible
    if (!("value2" %in% names(df))) df$value2 <- NA_real_
    if (!("value" %in% names(df))) {
      if (all(c("value2","single_author") %in% names(df))) df$value <- df$value2 + df$single_author
      else df$value <- NA_real_
    }
    if (!("value2_strength" %in% names(df))) df$value2_strength <- NA_real_
    if (!("FA" %in% names(df))) df$FA <- NA_real_
    if (!("LA" %in% names(df))) df$LA <- NA_real_
  }

  df
})

output$top1_title <- renderText({
  df <- .top1_nodes_df()
  if (is.null(df) || nrow(df) == 0) return("Top1: (no data yet)")
  dom <- input$top1_metadomain %||% ""
  metric <- input$top1_metric %||% ""
  paste0("Top1: ", dom, "  |  Metric: ", metric)
})

output$top1_table <- renderTable({
  df <- .top1_nodes_df()
  if (is.null(df) || nrow(df) == 0) return(NULL)

  metric_col <- input$top1_metric %||% "value"
  if (!(metric_col %in% names(df))) return(NULL)

  # filter NA metric unless requested
  if (!isTRUE(input$top1_include_na)) df <- df[!is.na(df[[metric_col]]), , drop=FALSE]
  if (nrow(df) == 0) return(NULL)

  df <- df[order(df[[metric_col]], decreasing = TRUE), , drop=FALSE]
  top <- df[1, , drop=FALSE]

  aac <- (.AAC_INLINE(df[[metric_col]]))

  .getnum <- function(x, col) {
    if (is.null(x) || !(col %in% names(x))) return(NA_real_)
    as.numeric(x[[col]][1])
  }

  out <- data.frame(
    Element = top$name,
    MetricValue = .getnum(top, metric_col),
    N_pubmed = .getnum(top, "n_pubmed"),
    `n(any byline appearance)` = .getnum(top, "n_byline"),
    value2 = .getnum(top, "value2"),
    self_coword = .getnum(top, "single_author"),
    `n(first/last author appearance)` = .getnum(top, "value"),
    FA = .getnum(top, "FA"),
    LA = .getnum(top, "LA"),
    value2_strength = .getnum(top, "value2_strength"),
    AAC = as.numeric(aac),
    check.names = FALSE
  )
  out
}, striped = TRUE, bordered = TRUE, hover = TRUE, spacing = "s")


  # ============================================================
  # FINAL REAL PLOT OVERRIDE (patched: uses renderSSplot(80).R + kano(63).R when present)
  # ============================================================
  # FINAL REAL PLOT OVERRIDE: use renderSSplot.R and kano.R when available
  # ------------------------------------------------------------
  # Purpose:
  # - Do NOT force fake/base SSplot when render_panel() exists.
  # - Do NOT force fake/base Kano when kano.R functions exist.
  # - Keep slopegraph safe so HTML/text objects do not trigger is.character(txt).
  # ============================================================

  .real_source_real_modules <- function() {
    for (.f in c("renderSSplot(80).R", "renderSSplot(79).R", "renderSSplot.R", "renderSSplot(78).R")) {
      if (file.exists(.f)) {
        try(source(.f, local = FALSE, encoding = "UTF-8"), silent = TRUE)
        break
      }
    }
    for (.f in c("kano(63).R", "kano.R", "kano(62).R")) {
      if (file.exists(.f)) {
        try(source(.f, local = FALSE, encoding = "UTF-8"), silent = TRUE)
        break
      }
    }
    invisible(TRUE)
  }
  .real_source_real_modules()

  .real_chr <- function(x) {
    x <- as.character(x)
    x[is.na(x)] <- ""
    x
  }

  .real_num <- function(x, default = 0) {
    z <- suppressWarnings(as.numeric(x))
    z[!is.finite(z) | is.na(z)] <- default
    z
  }

  .real_get_domain <- function(dom) {
    dom <- .real_chr(dom)[1]
    if (!nzchar(dom)) dom <- "Author"
    out <- tryCatch(.get_domain(dom), error = function(e) NULL)
    if (is.null(out)) out <- list(nodes = data.frame(), edges = data.frame())
    out
  }

  .real_top20_nodes <- function(nd, n = 20) {
    if (is.null(nd) || !is.data.frame(nd) || !nrow(nd)) return(data.frame())
    nd <- as.data.frame(nd, stringsAsFactors = FALSE, check.names = FALSE)
    if (!"name" %in% names(nd)) nd$name <- if ("id" %in% names(nd)) .real_chr(nd$id) else .real_chr(nd[[1]])
    nd$name <- trimws(.real_chr(nd$name))
    nd <- nd[nzchar(nd$name), , drop = FALSE]
    nd <- nd[!duplicated(nd$name), , drop = FALSE]
    if (!"value" %in% names(nd)) nd$value <- 1
    if (!"value2" %in% names(nd)) nd$value2 <- nd$value
    if (!"carac" %in% names(nd)) nd$carac <- 1
    if (!"ssi" %in% names(nd)) {
      if ("SSi" %in% names(nd)) nd$ssi <- nd$SSi else if ("sil_width" %in% names(nd)) nd$ssi <- nd$sil_width else nd$ssi <- 0
    }
    if (!"sil_width" %in% names(nd)) nd$sil_width <- nd$ssi
    if (!"a_star1" %in% names(nd)) {
      if ("a_star" %in% names(nd)) nd$a_star1 <- nd$a_star else if ("a_i" %in% names(nd)) nd$a_star1 <- 1/(1 + .real_num(nd$a_i, 0)) else nd$a_star1 <- 0
    }
    nd$value <- .real_num(nd$value, 0)
    nd$value2 <- .real_num(nd$value2, 0)
    nd$ssi <- .real_num(nd$ssi, 0)
    nd$sil_width <- .real_num(nd$sil_width, nd$ssi)
    nd$a_star1 <- .real_num(nd$a_star1, 0)
    nd$carac <- suppressWarnings(as.integer(gsub("[^0-9-]", "", .real_chr(nd$carac))))
    nd$carac[!is.finite(nd$carac) | is.na(nd$carac)] <- 1L
    if (!"rank" %in% names(nd)) nd$rank <- seq_len(nrow(nd))
    rr <- suppressWarnings(as.numeric(nd$rank))
    ord <- order(ifelse(is.finite(rr), rr, Inf), -nd$value, na.last = TRUE)
    nd <- nd[ord, , drop = FALSE]
    nd <- utils::head(nd, n)
    nd$rank <- seq_len(nrow(nd))
    nd
  }

  .real_edges_for_nodes <- function(ed, nd) {
    if (is.null(ed) || !is.data.frame(ed) || !nrow(ed) || is.null(nd) || !nrow(nd)) {
      return(data.frame(Leader=character(), Follower=character(), WCD=numeric(), stringsAsFactors = FALSE))
    }
    ed <- as.data.frame(ed, stringsAsFactors = FALSE, check.names = FALSE)
    if (all(c("from","to") %in% names(ed)) && !all(c("Leader","Follower") %in% names(ed))) {
      names(ed)[match(c("from","to"), names(ed))] <- c("Leader","Follower")
    }
    if ("follower" %in% names(ed) && !"Follower" %in% names(ed)) names(ed)[names(ed)=="follower"] <- "Follower"
    if (!"WCD" %in% names(ed)) {
      if ("weight" %in% names(ed)) ed$WCD <- ed$weight else if ("value" %in% names(ed)) ed$WCD <- ed$value else ed$WCD <- 1
    }
    if (!all(c("Leader","Follower","WCD") %in% names(ed))) {
      return(data.frame(Leader=character(), Follower=character(), WCD=numeric(), stringsAsFactors = FALSE))
    }
    ed <- ed[, c("Leader","Follower","WCD"), drop = FALSE]
    ed$Leader <- trimws(.real_chr(ed$Leader)); ed$Follower <- trimws(.real_chr(ed$Follower)); ed$WCD <- .real_num(ed$WCD, 1)
    keep <- nd$name
    ed <- ed[nzchar(ed$Leader) & nzchar(ed$Follower) & ed$Leader %in% keep & ed$Follower %in% keep & ed$WCD > 0, , drop = FALSE]
    ed
  }

  .real_modularity_results <- function(nd, ed) {
    nd <- .real_top20_nodes(nd, 20)
    ed <- .real_edges_for_nodes(ed, nd)
    clv <- sort(unique(stats::na.omit(nd$carac)))
    if (!length(clv)) clv <- 1L
    out <- do.call(rbind, lapply(clv, function(cc) {
      sub <- nd[nd$carac == cc, , drop = FALSE]
      data.frame(Cluster = paste0("C", cc),
                 SS = mean(sub$sil_width, na.rm = TRUE),
                 Qw = 0, Qu = 0,
                 D_GiniSimpson = NA_real_, Q_over_D = NA_real_,
                 OneMinus_1_over_k = NA_real_, Q_over_Dmax_eff = NA_real_,
                 n = nrow(sub), stringsAsFactors = FALSE)
    }))
    if (!nrow(ed) || nrow(nd) < 2) {
      overall <- data.frame(Cluster="OVERALL", SS=mean(nd$sil_width, na.rm=TRUE), Qw=0, Qu=0,
                            D_GiniSimpson=NA_real_, Q_over_D=NA_real_, OneMinus_1_over_k=NA_real_, Q_over_Dmax_eff=NA_real_,
                            n=nrow(nd), stringsAsFactors=FALSE)
      return(rbind(overall, out))
    }
    calcQ <- function(weighted = TRUE) {
      W <- matrix(0, nrow(nd), nrow(nd), dimnames = list(nd$name, nd$name))
      for (i in seq_len(nrow(ed))) {
        a <- ed$Leader[i]; b <- ed$Follower[i]; w <- if (weighted) ed$WCD[i] else 1
        if (a %in% nd$name && b %in% nd$name && a != b) { W[a,b] <- W[a,b] + w; W[b,a] <- W[b,a] + w }
      }
      m2 <- sum(W)
      if (!is.finite(m2) || m2 <= 0) return(list(Q=0, by=setNames(rep(0,length(clv)), paste0("C",clv))))
      k <- rowSums(W)
      cl <- setNames(nd$carac, nd$name)
      q_by <- setNames(numeric(length(clv)), paste0("C", clv))
      for (cc in clv) {
        ids <- names(cl)[cl == cc]
        if (length(ids)) q_by[paste0("C",cc)] <- sum(W[ids, ids, drop=FALSE] - outer(k[ids], k[ids]) / m2) / m2
      }
      list(Q=sum(q_by), by=q_by)
    }
    qw <- calcQ(TRUE); qu <- calcQ(FALSE)
    out$Qw <- as.numeric(qw$by[out$Cluster]); out$Qu <- as.numeric(qu$by[out$Cluster])
    out$Qw[!is.finite(out$Qw)] <- 0; out$Qu[!is.finite(out$Qu)] <- 0
    cc_all <- nd$carac; pk <- as.numeric(table(cc_all)) / length(cc_all); D <- 1 - sum(pk^2); Dmax <- 1 - 1 / max(1, length(unique(cc_all)))
    overall <- data.frame(Cluster="OVERALL", SS=mean(nd$sil_width, na.rm=TRUE), Qw=qw$Q, Qu=qu$Q,
                          D_GiniSimpson=D, Q_over_D=ifelse(is.finite(D) && D>0, qu$Q/D, NA_real_),
                          OneMinus_1_over_k=Dmax, Q_over_Dmax_eff=ifelse(is.finite(Dmax) && Dmax>0, qu$Q/Dmax, NA_real_),
                          n=nrow(nd), stringsAsFactors=FALSE)
    rbind(overall, out)
  }

  .real_print_plot_object <- function(obj) {
    if (is.null(obj)) return(invisible(TRUE))
    if (inherits(obj, "recordedplot")) { replayPlot(obj); return(invisible(TRUE)) }
    if (inherits(obj, "ggplot")) { print(obj); return(invisible(TRUE)) }
    if (inherits(obj, "grob") || inherits(obj, "gTree") || inherits(obj, "gtable")) {
      if (requireNamespace("grid", quietly = TRUE)) grid::grid.draw(obj)
      return(invisible(TRUE))
    }
    if (is.function(obj)) { obj(); return(invisible(TRUE)) }
    invisible(TRUE)
  }

  .real_draw_kano_value_value2 <- function(nd, ed, title_txt = "Kano: value vs value2") {
    .real_source_real_modules()
    nd <- .real_top20_nodes(nd, 20); ed <- .real_edges_for_nodes(ed, nd)
    if (!nrow(nd)) { plot.new(); text(0.5,0.5,"No Top20 nodes available"); return(invisible()) }
    # Real Kano from kano.R: two blue wings plus hub/outer circles.
    if (exists("plot_kano_real", mode = "function")) {
      obj <- plot_kano_real(nodes = nd, edges = ed, title_txt = title_txt, visual_ratio = 1/1.5, label_size = 4)
      .real_print_plot_object(obj); return(invisible(TRUE))
    }
    if (exists("kano_plot", mode = "function")) {
      obj <- kano_plot(nd, ed, title_txt = title_txt, xlab = "value2", ylab = "value", visual_ratio = 1/1.5)
      .real_print_plot_object(obj); return(invisible(TRUE))
    }
    plot.new(); text(0.5, 0.5, "kano.R not found: place kano.R/kano(63).R in the app folder to draw the real two-wing Kano plot.", cex = 1.15, font = 2, col = "red")
  }

  .real_draw_kano_ss_astar <- function(nd, ed, title_txt = "Kano: SS vs a*") {
    .real_source_real_modules()
    nd <- .real_top20_nodes(nd, 20); ed <- .real_edges_for_nodes(ed, nd)
    if (!nrow(nd)) { plot.new(); text(0.5,0.5,"No Top20 nodes available"); return(invisible()) }
    nd$ss <- if ("sil_width" %in% names(nd)) nd$sil_width else nd$ssi
    nd$ss <- .real_num(nd$ss, 0)
    nd$a_star <- if ("a_star1" %in% names(nd)) nd$a_star1 else if ("a_star" %in% names(nd)) nd$a_star else 0
    nd$a_star <- .real_num(nd$a_star, 0)
    # Real SS-Kano from kano.R: map x=SS, y=a* into plot_kano_real() so the two wings/circles remain.
    if (exists("kano_plot_ss_astar", mode = "function")) {
      obj <- kano_plot_ss_astar(nd, ed, xlab = "SS", ylab = "a*", title_txt = title_txt, visual_ratio = 1/1.5)
      .real_print_plot_object(obj); return(invisible(TRUE))
    }
    if (exists("plot_kano_ss_astar", mode = "function")) {
      obj <- plot_kano_ss_astar(nd, title_txt = title_txt)
      .real_print_plot_object(obj); return(invisible(TRUE))
    }
    if (exists("plot_kano_real", mode = "function")) {
      nd2 <- nd
      nd2$value <- nd$a_star
      nd2$value2 <- nd$ss
      obj <- plot_kano_real(nodes = nd2, edges = ed, title_txt = title_txt, visual_ratio = 1/1.5, label_size = 4)
      .real_print_plot_object(obj); return(invisible(TRUE))
    }
    plot.new(); text(0.5, 0.5, "kano.R not found: place kano.R/kano(63).R in the app folder to draw the real SS-Kano plot.", cex = 1.15, font = 2, col = "red")
  }

  .real_draw_ssplot <- function(nd, ed, title_txt = "SSplot") {
    nd20 <- .real_top20_nodes(nd, 20)
    ed20 <- .real_edges_for_nodes(ed, nd20)
    if (!nrow(nd20)) { plot.new(); text(0.5, 0.5, "No Top20 nodes for SSplot"); return(invisible()) }
    nd20$sil_width <- nd20$ssi
    if (!"wsel" %in% names(nd20)) nd20$wsel <- NA_real_
    if (!"role" %in% names(nd20)) nd20$role <- NA_character_
    if (!"neighbor_name" %in% names(nd20)) nd20$neighbor_name <- nd20$name
    if (!"neighborC" %in% names(nd20)) nd20$neighborC <- nd20$carac
    if (nrow(ed20)) {
      for (i in seq_len(nrow(ed20))) {
        j <- match(ed20$Follower[i], nd20$name)
        if (!is.na(j)) {
          nd20$wsel[j] <- ed20$WCD[i]
          nd20$role[j] <- "follower"
          nd20$neighbor_name[j] <- ed20$Leader[i]
          nd20$neighborC[j] <- nd20$carac[match(ed20$Leader[i], nd20$name)]
        }
      }
    }
    leader_idx <- tapply(seq_len(nrow(nd20)), nd20$carac, function(ii) ii[order(-nd20$value2[ii], -nd20$value[ii], nd20$rank[ii])][1])
    nd20$role[as.integer(leader_idx)] <- "leader"
    res <- .real_modularity_results(nd20, ed20)
    if (!exists("render_panel", mode = "function")) stop("render_panel() not found. Put renderSSplot.R/renderSSplot(79).R in the app folder.")
    render_panel(sil_df = nd20, nodes0 = nd20, results = res, nodes = nd20,
                 top_n = nrow(nd20), font_scale = input$ss_font_scale %||% 1.3,
                 aac_side = "left", neighbor_side = "right", neighbor_on_bar = TRUE,
                 footer_label = title_txt)
  }

  .real_base_slope <- function(yc, title = "Slopegraph") {
    if (is.null(yc) || !is.data.frame(yc) || !nrow(yc)) { plot.new(); text(0.5, 0.5, "No slopegraph data", cex = 1.4, font = 2); return(invisible(NULL)) }
    yc <- as.data.frame(yc, stringsAsFactors = FALSE)
    if (!all(c("Year", "Count", "item") %in% names(yc))) { plot.new(); text(0.5, 0.5, "Slopegraph needs Year, item, Count", cex = 1.2, font = 2); return(invisible(NULL)) }
    yc$Year <- suppressWarnings(as.integer(as.character(yc$Year)))
    yc$Count <- .real_num(yc$Count, 0)
    yc$item <- .real_chr(yc$item)
    yc <- yc[is.finite(yc$Year) & nzchar(yc$item), , drop = FALSE]
    if (!nrow(yc)) { plot.new(); text(0.5, 0.5, "No finite slopegraph data", cex = 1.4, font = 2); return(invisible(NULL)) }
    tot <- aggregate(Count ~ item, data = yc, sum, na.rm = TRUE)
    keep <- head(tot$item[order(-tot$Count)], 10)
    yc <- yc[yc$item %in% keep, , drop = FALSE]
    yrs <- sort(unique(yc$Year))
    old <- par(no.readonly = TRUE); on.exit(par(old), add = TRUE)
    par(mar = c(5.5, 9.0, 4.5, 9.0), family = "sans", xpd = NA)
    ylim <- range(yc$Count, finite = TRUE); if (diff(ylim) == 0) ylim <- ylim + c(-1,1)
    plot(range(yrs), ylim, type = "n", xlab = "Year", ylab = "Count", main = title, cex.main = 1.45, font.main = 2, cex.lab = 1.25, font.lab = 2)
    grid(col = "grey88")
    pal <- grDevices::hcl.colors(max(3, length(keep)), "Dark 3")
    for (i in seq_along(keep)) {
      d <- yc[yc$item == keep[i], , drop = FALSE]
      d <- d[order(d$Year), , drop = FALSE]
      lines(d$Year, d$Count, type = "b", lwd = 2.2, pch = 19, col = pal[i])
      if (nrow(d)) {
        text(d$Year[1], d$Count[1], labels = keep[i], pos = 2, cex = 0.85, font = 2, col = pal[i])
        text(d$Year[nrow(d)], d$Count[nrow(d)], labels = keep[i], pos = 4, cex = 0.85, font = 2, col = pal[i])
      }
    }
    invisible(NULL)
  }

  output$ssplot_panel <- renderPlot({
    tryCatch({
      dom <- .real_chr(input$ss_domain)[1]; if (!nzchar(dom)) dom <- "Author"
      d <- .real_get_domain(dom)
      .real_draw_ssplot(d$nodes, d$edges, title_txt = paste0("SSplot: FLCA-MA-SIL Top20 - ", dom))
    }, error = function(e) { plot.new(); text(0.5, 0.5, paste("SSplot error:", conditionMessage(e)), col = "red", cex = 1.1, font = 2) })
  }, height = 950)

  output$kano_vv2_plot <- renderPlot({
    tryCatch({
      dom <- .real_chr(input$kano_domain)[1]; if (!nzchar(dom)) dom <- "Author"
      d <- .real_get_domain(dom)
      .real_draw_kano_value_value2(d$nodes, d$edges, title_txt = paste0("Kano: value vs value2 - ", dom))
    }, error = function(e) { plot.new(); text(0.5, 0.5, paste("Kano error:", conditionMessage(e)), col = "red", cex = 1.1, font = 2) })
  }, height = 1250)

  output$kano_ss_astar_plot <- renderPlot({
    tryCatch({
      dom <- .real_chr(input$kano_domain)[1]; if (!nzchar(dom)) dom <- "Author"
      d <- .real_get_domain(dom)
      .real_draw_kano_ss_astar(d$nodes, d$edges, title_txt = paste0("Kano: SS vs a* - ", dom))
    }, error = function(e) { plot.new(); text(0.5, 0.5, paste("Kano SS/a* error:", conditionMessage(e)), col = "red", cex = 1.1, font = 2) })
  }, height = 1250)

  output$ss_kano_plot <- renderPlot({
    tryCatch({
      dom <- .real_chr(input$ss_kano_domain)[1]; if (!nzchar(dom)) dom <- "Author"
      d <- .real_get_domain(dom)
      .real_draw_kano_ss_astar(d$nodes, d$edges, title_txt = paste0("SS Kano: SS vs a* - ", dom))
    }, error = function(e) { plot.new(); text(0.5, 0.5, paste("SSKano error:", conditionMessage(e)), col = "red", cex = 1.1, font = 2) })
  }, height = 1250)

  output$plt_slope_top2_country <- renderPlot({
    tryCatch({
      req(rv$pubmeta)
      doms <- input$yc_domains; if (is.null(doms) || !length(doms)) doms <- c("Author", "Country")
      yc <- compute_combo_year_counts(rv, domains = doms, recent_n_years = if (is.null(input$yc_recent_years)) 10 else input$yc_recent_years)
      .real_base_slope(yc, title = "Top10 combo slopegraph")
    }, error = function(e) { plot.new(); text(0.5, 0.5, paste("Slopegraph error:", conditionMessage(e)), col = "red", cex = 1.1, font = 2) })
  }, height = 620)

  output$plt_slope_top2_author <- renderPlot({
    tryCatch({
      req(rv$pubmeta)
      yc <- compute_combo_year_counts(rv, domains = c("Author"), recent_n_years = if (is.null(input$yc_recent_years)) 10 else input$yc_recent_years)
      .real_base_slope(yc, title = "Author metadata slopegraph")
    }, error = function(e) { plot.new(); text(0.5, 0.5, paste("Author slopegraph error:", conditionMessage(e)), col = "red", cex = 1.1, font = 2) })
  }, height = 620)

  output$plt_slope_top2_mesh <- renderPlot({
    tryCatch({
      req(rv$pubmeta)
      yc <- compute_combo_year_counts(rv, domains = c("MeSH"), recent_n_years = if (is.null(input$yc_recent_years)) 10 else input$yc_recent_years)
      .real_base_slope(yc, title = "MeSH metadata slopegraph")
    }, error = function(e) { plot.new(); text(0.5, 0.5, paste("MeSH slopegraph error:", conditionMessage(e)), col = "red", cex = 1.1, font = 2) })
  }, height = 620)

  output$plt_slope_top2_author_fl <- renderPlot({
    tryCatch({
      req(rv$pubmeta)
      doms <- input$yc_domains_fl; if (is.null(doms) || !length(doms)) doms <- c("Author (1st+Last)", "Country")
      yc <- compute_combo_year_counts(rv, domains = doms, recent_n_years = if (is.null(input$yc_recent_years_fl)) 10 else input$yc_recent_years_fl)
      .real_base_slope(yc, title = "Term/Year slopegraph")
    }, error = function(e) { plot.new(); text(0.5, 0.5, paste("Term/Year slopegraph error:", conditionMessage(e)), col = "red", cex = 1.1, font = 2) })
  }, height = 620)



  # ============================================================
  # FINAL FIX: real Tufte-form slopegraph, Top10 for each selected entity/domain
  # - Top10 is recalculated within the selected entity/domain (Author, Country, etc.).
  # - Data path follows the requested form: wide data -> long data -> tufte_sort -> slope plot.
  # - Uses base plotting for the final rendering to avoid htmltools/grid errors such as:
  #   "不是所有的 is.character(txt) 都是 TRUE".
  # ============================================================
  .sg_final_chr <- function(x) {
    x <- as.character(x)
    x[is.na(x)] <- ""
    x <- iconv(x, from = "", to = "UTF-8", sub = "")
    x[is.na(x)] <- ""
    trimws(x)
  }

  .sg_final_num <- function(x, default = 0) {
    x <- suppressWarnings(as.numeric(as.character(x)))
    x[!is.finite(x) | is.na(x)] <- default
    x
  }

  .sg_final_one <- function(x, default = NULL) {
    if (is.null(x) || length(x) == 0 || all(is.na(x))) return(default)
    x[[1]]
  }

  .sg_final_domain <- function(x) {
    z <- .sg_final_chr(.sg_final_one(x, "Author"))[1]
    allowed <- c("Author", "Journal", "Country", "State/Province", "Institute", "Department", "MeSH")
    if (!nzchar(z) || !z %in% allowed) z <- "Author"
    z
  }

  .sg_final_recent <- function(x) {
    z <- suppressWarnings(as.integer(.sg_final_one(x, 10)))
    if (!is.finite(z) || is.na(z) || z < 2) z <- 10L
    z
  }

  .sg_final_blank_plot <- function(msg = "No slopegraph data available.") {
    graphics::plot.new()
    graphics::text(0.5, 0.5, .sg_final_chr(msg)[1], col = "red", cex = 1.15, font = 2)
    invisible(NULL)
  }

  .sg_final_build_wide <- function(domain = "Author", recent_n_years = 10, top_n = 10) {
    domain <- .sg_final_domain(domain)
    recent_n_years <- .sg_final_recent(recent_n_years)

    # Build Top10 directly from article-level metadata lists.
    # This avoids compute_combo_year_counts(), which can inherit node order or
    # FLCA/network filtering and make Author look alphabetically sorted.
    yrs <- rv$article_years
    if (is.null(yrs) && is.data.frame(rv$pubmeta) && "Year" %in% names(rv$pubmeta)) yrs <- rv$pubmeta$Year
    yrs <- suppressWarnings(as.integer(as.character(yrs)))
    if (!any(is.finite(yrs))) return(data.frame())

    years_all <- sort(unique(yrs[is.finite(yrs)]))
    years_keep <- utils::tail(years_all, min(recent_n_years, length(years_all)))
    if (length(years_keep) < 2) return(data.frame())

    tl <- tryCatch(.get_domain_term_list(rv, domain), error = function(e) NULL)

    # MeSH can be stored in pubmeta and should remain article-year aligned.
    if (tolower(domain) %in% c("mesh", "mesh term", "meshterm") && is.data.frame(rv$pubmeta)) {
      cn_mesh <- intersect(c("MeSH", "Mesh", "MESH", "mesh"), names(rv$pubmeta))
      if (length(cn_mesh)) tl <- as.list(as.character(rv$pubmeta[[cn_mesh[1]]]))
    }

    # Fallback to pubmeta columns only if the article-level list does not exist.
    if (is.null(tl) && is.data.frame(rv$pubmeta)) {
      d <- tolower(domain)
      cand <- switch(
        d,
        "author" = c("Author", "Authors", "AU", "author", "authors"),
        "journal" = c("Journal", "Source", "Source Title", "journal"),
        "country" = c("Country", "Countries", "country"),
        "state/province" = c("State/Province", "State", "Province", "state", "province", "state_province", "stateprov"),
        "institute" = c("Institute", "Inst", "Affiliation", "institute", "inst"),
        "department" = c("Department", "Dept", "department", "dept"),
        "mesh" = c("MeSH", "Mesh", "MESH", "mesh"),
        c()
      )
      cn <- intersect(cand, names(rv$pubmeta))
      if (length(cn)) tl <- as.list(as.character(rv$pubmeta[[cn[1]]]))
    }
    if (is.null(tl) || !length(tl)) return(data.frame())

    n <- min(length(tl), length(yrs))
    tl <- tl[seq_len(n)]
    yv <- yrs[seq_len(n)]
    keep_i <- which(is.finite(yv) & yv %in% years_keep)
    if (!length(keep_i)) return(data.frame())

    splitter <- if (tolower(domain) %in% c("mesh", "mesh term", "meshterm")) .split_mesh_terms_safe else .split_terms2

    rec <- vector("list", length(keep_i))
    jj <- 0L
    for (i in keep_i) {
      terms <- tryCatch(splitter(tl[[i]]), error = function(e) character(0))
      terms <- .sg_final_chr(terms)
      terms <- terms[nzchar(terms)]
      terms <- terms[.is_single_metadata_item(terms)]
      terms <- unique(terms)  # one publication contributes once per entity per year
      if (!length(terms)) next

      jj <- jj + 1L
      rec[[jj]] <- data.frame(
        item = terms,
        Year = rep(as.integer(yv[i]), length(terms)),
        Count = rep(1, length(terms)),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
    if (jj == 0L) return(data.frame())
    rec <- rec[seq_len(jj)]

    yc_raw <- do.call(rbind, rec)
    yc_raw$item <- .sg_final_chr(yc_raw$item)
    yc_raw$Year <- suppressWarnings(as.integer(as.character(yc_raw$Year)))
    yc_raw$Count <- .sg_final_num(yc_raw$Count, 0)
    yc_raw <- yc_raw[nzchar(yc_raw$item) & is.finite(yc_raw$Year) & yc_raw$Year %in% years_keep, , drop = FALSE]
    if (!nrow(yc_raw)) return(data.frame())

    # Select Top10 within the selected domain by TOTAL count in the chosen year window.
    totals <- stats::aggregate(Count ~ item, data = yc_raw, FUN = sum, na.rm = TRUE)
    names(totals)[names(totals) == "Count"] <- "Total"
    totals$item <- .sg_final_chr(totals$item)
    totals$Total <- .sg_final_num(totals$Total, 0)
    totals <- totals[totals$Total > 0 & nzchar(totals$item), , drop = FALSE]
    totals <- totals[order(-totals$Total, totals$item), , drop = FALSE]
    totals <- utils::head(totals, min(top_n, nrow(totals)))
    if (!nrow(totals)) return(data.frame())
    totals$OverallRank <- seq_len(nrow(totals))

    # Now compute yearly counts only for those Top10 items.
    yc <- yc_raw[yc_raw$item %in% totals$item, , drop = FALSE]
    yc <- stats::aggregate(Count ~ item + Year, data = yc, FUN = sum, na.rm = TRUE)

    full <- expand.grid(item = totals$item, Year = years_keep, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    full <- merge(full, yc, by = c("item", "Year"), all.x = TRUE, sort = FALSE)
    full$Count <- .sg_final_num(full$Count, 0)
    full <- merge(full, totals[, c("item", "OverallRank", "Total"), drop = FALSE], by = "item", all.x = TRUE, sort = FALSE)
    full <- full[order(full$OverallRank, full$Year), , drop = FALSE]

    mat <- xtabs(Count ~ item + Year, data = full)
    mat <- mat[as.character(totals$item), as.character(years_keep), drop = FALSE]
    wide <- data.frame(group = rownames(mat), as.data.frame.matrix(mat),
                       stringsAsFactors = FALSE, check.names = FALSE)
    names(wide)[seq_along(years_keep) + 1L] <- as.character(years_keep)
    wide <- merge(wide, totals[, c("item", "OverallRank", "Total"), drop = FALSE],
                  by.x = "group", by.y = "item", all.x = TRUE, sort = FALSE)
    wide <- wide[order(wide$OverallRank), , drop = FALSE]
    rownames(wide) <- NULL
    wide
  }


  .sg_final_long_from_wide <- function(data_wide) {
    if (is.null(data_wide) || !is.data.frame(data_wide) || !nrow(data_wide)) return(data.frame())
    data_wide <- as.data.frame(data_wide, stringsAsFactors = FALSE, check.names = FALSE)
    if (!"group" %in% names(data_wide)) names(data_wide)[1] <- "group"
    year_cols <- setdiff(names(data_wide), c("group", "OverallRank", "Total"))
    year_cols <- year_cols[grepl("^[0-9]{4}$", year_cols)]
    if (!length(year_cols)) return(data.frame())

    out <- do.call(rbind, lapply(year_cols, function(cc) {
      data.frame(group = .sg_final_chr(data_wide$group),
                 year  = .sg_final_chr(cc),
                 value = .sg_final_num(data_wide[[cc]], 0),
                 stringsAsFactors = FALSE, check.names = FALSE)
    }))
    out$group <- .sg_final_chr(out$group)
    out$year <- .sg_final_chr(out$year)
    out$value <- .sg_final_num(out$value, 0)
    attr(out, "year_cols") <- year_cols
    out
  }

  .sg_final_tufte_sort <- function(df, x = "year", y = "value", group = "group", method = "tufte", min.space = 0.05) {
    ids <- match(c(x, y, group), names(df))
    if (any(is.na(ids))) return(data.frame())
    df <- df[, ids, drop = FALSE]
    names(df) <- c("x", "y", "group")
    df$x <- .sg_final_chr(df$x)
    df$group <- .sg_final_chr(df$group)
    df$y <- .sg_final_num(df$y, 0)

    x_levels <- unique(df$x)
    if (all(grepl("^[0-9]{4}$", x_levels))) x_levels <- as.character(sort(as.integer(x_levels)))

    # group order arrives from .sg_final_build_wide(): Total Count descending.
    # To display the highest-count Top10 item at the top in a base plot,
    # the internal Tufte stacking starts from the lowest-count item.
    group_levels_desc <- unique(df$group)
    group_levels_bottom_to_top <- rev(group_levels_desc)

    tmp_grid <- expand.grid(x = x_levels, group = group_levels_bottom_to_top, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    tmp <- merge(df, tmp_grid, by = c("x", "group"), all.y = TRUE, sort = FALSE)
    tmp$y <- .sg_final_num(tmp$y, 0)

    mat <- xtabs(y ~ group + x, data = tmp)
    mat <- mat[group_levels_bottom_to_top, x_levels, drop = FALSE]
    tmp_wide <- data.frame(group = rownames(mat), as.data.frame.matrix(mat),
                           stringsAsFactors = FALSE, check.names = FALSE)
    names(tmp_wide)[seq_along(x_levels) + 1L] <- x_levels

    value_mat <- as.matrix(tmp_wide[, x_levels, drop = FALSE])
    storage.mode(value_mat) <- "numeric"
    rng <- range(value_mat, na.rm = TRUE, finite = TRUE)
    span <- diff(rng)
    if (!is.finite(span) || span <= 0) span <- 1
    gap <- max(min.space * span, 0.85)

    yshift <- numeric(nrow(tmp_wide))
    if (nrow(tmp_wide) >= 2) {
      for (i in 2:nrow(tmp_wide)) {
        m2 <- as.matrix(tmp_wide[(i - 1):i, x_levels, drop = FALSE])
        storage.mode(m2) <- "numeric"
        d_min <- suppressWarnings(min(diff(m2), na.rm = TRUE))
        if (!is.finite(d_min)) d_min <- 0
        yshift[i] <- ifelse(d_min < gap, gap - d_min, 0)
      }
    }
    tmp_wide$yshift <- cumsum(yshift)

    out <- do.call(rbind, lapply(x_levels, function(xx) {
      data.frame(group = .sg_final_chr(tmp_wide$group),
                 yshift = .sg_final_num(tmp_wide$yshift, 0),
                 x = .sg_final_chr(xx),
                 y = .sg_final_num(tmp_wide[[xx]], 0),
                 stringsAsFactors = FALSE, check.names = FALSE)
    }))
    out$ypos <- out$y + out$yshift
    out$x <- factor(out$x, levels = x_levels, labels = x_levels)
    out
  }


  .sg_final_base_tufte_plot <- function(df, title = "Estimates of Values over Years", label_digits = 0) {
    if (is.null(df) || !is.data.frame(df) || !nrow(df)) return(.sg_final_blank_plot())
    df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
    df$group <- .sg_final_chr(df$group)
    df$x_chr <- .sg_final_chr(as.character(df$x))
    df$y <- .sg_final_num(df$y, 0)
    df$ypos <- .sg_final_num(df$ypos, 0)

    x_levels <- unique(df$x_chr)
    x_pos <- seq_along(x_levels)
    names(x_pos) <- x_levels

    first_x <- x_levels[1]
    left <- df[df$x_chr == first_x & is.finite(df$ypos), c("group", "ypos"), drop = FALSE]

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(oldpar), add = TRUE)
    graphics::par(mar = c(4.5, 20, 4.5, 2.5), xpd = NA)

    y_range <- range(df$ypos, na.rm = TRUE, finite = TRUE)
    if (!all(is.finite(y_range)) || diff(y_range) <= 0) y_range <- y_range + c(-1, 1)

    graphics::plot(x_pos[df$x_chr], df$ypos, type = "n", axes = FALSE,
                   xlab = "", ylab = "", main = .sg_final_chr(title)[1],
                   xlim = c(min(x_pos) - 0.15, max(x_pos) + 0.15),
                   ylim = y_range)
    graphics::axis(1, at = x_pos, labels = .sg_final_chr(x_levels), font = 2, cex.axis = 1.0)

    # Larger bold entity labels on the left. Labels are not drawn through axis()
    # because axis() can shrink them and makes long Author labels hard to read.
    graphics::axis(2, at = left$ypos, labels = FALSE, las = 1, tick = FALSE)
    graphics::text(rep(min(x_pos) - 0.10, nrow(left)), left$ypos,
                   labels = .sg_final_chr(left$group),
                   adj = 1, font = 2, cex = 1.08)

    for (g in unique(df$group)) {
      d <- df[df$group == g, , drop = FALSE]
      d <- d[match(x_levels, d$x_chr), , drop = FALSE]
      xp <- x_pos[d$x_chr]
      graphics::lines(xp, d$ypos, col = "red", lwd = 2)
      graphics::points(xp, d$ypos, pch = 21, bg = "white", col = "black", cex = 2.6, lwd = 1)

      # Blank labels for zero counts: keep the white circle, remove the "0".
      nz <- which(is.finite(d$y) & d$y != 0)
      if (length(nz)) {
        lab <- sprintf(paste0("%.", as.integer(label_digits), "f"), d$y[nz])
        graphics::text(xp[nz], d$ypos[nz], labels = .sg_final_chr(lab), cex = 0.78)
      }
    }
    graphics::box(bty = "l")
    invisible(NULL)
  }


  .sg_final_draw_domain <- function(domain = "Author", recent_n_years = 10) {
    domain <- .sg_final_domain(domain)
    recent_n_years <- .sg_final_recent(recent_n_years)
    wide <- .sg_final_build_wide(domain, recent_n_years, top_n = 10)
    if (is.null(wide) || !is.data.frame(wide) || !nrow(wide)) {
      return(.sg_final_blank_plot(paste0("No Top10 year-count data for ", domain, ". Run analysis first or choose another entity.")))
    }
    data_long <- .sg_final_long_from_wide(wide)
    df <- .sg_final_tufte_sort(data_long, x = "year", y = "value", group = "group", method = "tufte", min.space = 0.05)
    .sg_final_base_tufte_plot(df, title = paste0(domain, " Top10 real slopegraph over recent ", recent_n_years, " years"), label_digits = 0)
  }

  .sg_final_table <- function(domain = "Author", recent_n_years = 10) {
    wide <- .sg_final_build_wide(domain, recent_n_years, top_n = 10)
    if (is.null(wide) || !is.data.frame(wide) || !nrow(wide)) {
      return(data.frame(Message = "No Top10 year-count data.", stringsAsFactors = FALSE))
    }
    wide
  }

  output$slope_combo_entity_note <- renderText({
    "Real Tufte-style slopegraph: Top10 is recalculated within the selected entity/domain; y positions use yearly Count values plus Tufte spacing."
  })

  output$plot_slope_combo_entity_top10 <- renderPlot({
    req(rv$pubmeta)
    .sg_final_draw_domain(.sg_final_domain(input$slope_combo_entity_domain), .sg_final_recent(input$slope_combo_recent_years))
  }, height = 760)

  output$tbl_slope_combo_entity_top10 <- DT::renderDT({
    req(rv$pubmeta)
    DT::datatable(
      .sg_final_table(.sg_final_domain(input$slope_combo_entity_domain), .sg_final_recent(input$slope_combo_recent_years)),
      rownames = FALSE,
      options = list(pageLength = 10, scrollX = TRUE)
    )
  })

  output$slope_combo_entity_plot <- renderPlot({
    req(rv$pubmeta)
    .sg_final_draw_domain(.sg_final_domain(input$slope_simple_entity_domain), .sg_final_recent(input$slope_simple_recent_years))
  }, height = 760)

  output$slope_combo_entity_table <- DT::renderDT({
    req(rv$pubmeta)
    DT::datatable(
      .sg_final_table(.sg_final_domain(input$slope_simple_entity_domain), .sg_final_recent(input$slope_simple_recent_years)),
      rownames = FALSE,
      options = list(pageLength = 10, scrollX = TRUE)
    )
  })



  # The former standalone extra-analysis output block has been deleted.
  # These graphics are now rendered only inside the Slope tab below.

  # ============================================================
  # Slope-tab integrated extra plots
  # User request: draw these plots inside Slope only.
  # All three plots use the stable Slope Top10 long table as the single source.
  # ============================================================
  .sg_slope_extra_domain <- function() {
    z <- input$slope_combo_entity_domain
    if (is.null(z) || length(z) == 0 || !nzchar(.sg_final_chr(z)[1])) z <- input$slope_simple_entity_domain
    if (is.null(z) || length(z) == 0 || !nzchar(.sg_final_chr(z)[1])) z <- input$combo_domain
    if (is.null(z) || length(z) == 0 || !nzchar(.sg_final_chr(z)[1])) z <- "Author"
    z <- .sg_final_chr(z)[1]
    if (identical(z, "Journal/Year")) z <- "Journal"
    .sg_final_domain(z)
  }

  .sg_slope_extra_recent <- function() {
    z <- input$slope_combo_recent_years
    if (is.null(z) || length(z) == 0) z <- input$slope_simple_recent_years
    if (is.null(z) || length(z) == 0) z <- 10
    .sg_final_recent(z)
  }

  .sg_slope_extra_long <- function() {
    dom <- .sg_slope_extra_domain()
    wide <- .sg_final_build_wide(dom, .sg_slope_extra_recent(), top_n = 10)
    if (is.null(wide) || !is.data.frame(wide) || !nrow(wide)) return(data.frame())
    long <- .sg_final_long_from_wide(wide)
    if (is.null(long) || !is.data.frame(long) || !nrow(long)) return(data.frame())

    out <- data.frame(
      Domain = rep(dom, nrow(long)),
      term   = .sg_final_chr(long$group),
      Year   = suppressWarnings(as.integer(.sg_final_chr(long$year))),
      Count  = .sg_final_num(long$value, 0),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    out <- out[nzchar(out$term) & is.finite(out$Year), , drop = FALSE]
    if (!nrow(out)) return(data.frame())

    term_order <- .sg_final_chr(wide$group)
    term_order <- term_order[nzchar(term_order)]
    out$term_order <- match(out$term, term_order)
    out <- out[order(out$term_order, out$Year), , drop = FALSE]

    mean_df <- stats::aggregate(Count ~ term, data = out, FUN = mean, na.rm = TRUE)
    names(mean_df)[names(mean_df) == "Count"] <- "mean_count"
    out <- merge(out, mean_df, by = "term", all.x = TRUE, sort = FALSE)
    out$mean_count <- .sg_final_num(out$mean_count, 0)
    out$is_burst <- out$Count > out$mean_count
    out <- out[order(out$term_order, out$Year), c("Domain","term","Year","Count","mean_count","is_burst"), drop = FALSE]
    rownames(out) <- NULL
    out
  }

  .sg_slope_extra_blank <- function(msg) {
    graphics::plot.new()
    graphics::text(0.5, 0.5, .sg_final_chr(msg)[1], col = "red", cex = 1.15, font = 2)
    invisible(NULL)
  }

  .sg_slope_extra_draw_heatmap <- function() {
    df <- .sg_slope_extra_long()
    dom <- .sg_slope_extra_domain()
    if (is.null(df) || !is.data.frame(df) || !nrow(df)) {
      return(.sg_slope_extra_blank(paste0("No Top10 long table data for ", dom, ". Run analysis first or choose another Slope entity.")))
    }
    terms <- unique(.sg_final_chr(df$term))
    yrs <- sort(unique(suppressWarnings(as.integer(df$Year))))
    if (!length(terms) || !length(yrs)) return(.sg_slope_extra_blank("No finite heatmap data."))
    df$term <- factor(.sg_final_chr(df$term), levels = terms)
    df$Year <- factor(as.character(as.integer(df$Year)), levels = as.character(yrs))
    mat <- stats::xtabs(Count ~ term + Year, data = df)
    mat <- mat[terms, as.character(yrs), drop = FALSE]

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(oldpar), add = TRUE)
    graphics::par(mar = c(5.2, 18, 4.2, 2), xpd = NA)
    z <- t(as.matrix(mat)); storage.mode(z) <- "numeric"
    graphics::image(seq_along(yrs), seq_along(terms), z, axes = FALSE,
                    xlab = "Year", ylab = "",
                    main = paste0(dom, " Top10 Heatmap"),
                    col = grDevices::colorRampPalette(c("blue", "white", "red"))(64))
    graphics::axis(1, at = seq_along(yrs), labels = .sg_final_chr(as.character(yrs)), las = 2, font = 2)
    graphics::axis(2, at = seq_along(terms), labels = .sg_final_chr(terms), las = 1, font = 2, cex.axis = 0.85)
    graphics::box()
    burst <- df[df$is_burst %in% TRUE, , drop = FALSE]
    if (nrow(burst)) {
      bx <- match(.sg_final_chr(as.character(burst$Year)), .sg_final_chr(as.character(yrs)))
      by <- match(.sg_final_chr(as.character(burst$term)), .sg_final_chr(terms))
      ok <- is.finite(bx) & is.finite(by)
      if (any(ok)) graphics::points(bx[ok], by[ok], pch = 21, bg = "red", col = "black", cex = 1.25, lwd = 0.8)
    }
    graphics::legend("topright", legend = "Count > term mean", pch = 21, pt.bg = "red", bty = "n", cex = 0.9)
    invisible(NULL)
  }

  .sg_slope_extra_draw_timeline <- function() {
    df <- .sg_slope_extra_long()
    dom <- .sg_slope_extra_domain()
    if (is.null(df) || !is.data.frame(df) || !nrow(df)) {
      return(.sg_slope_extra_blank(paste0("No Top10 timeline data for ", dom, ".")))
    }
    terms <- unique(.sg_final_chr(df$term))
    yrs <- sort(unique(suppressWarnings(as.integer(df$Year))))
    if (!length(terms) || !length(yrs)) return(.sg_slope_extra_blank("No finite timeline data."))
    df$term <- .sg_final_chr(df$term)
    df$Count <- .sg_final_num(df$Count, 0)
    max_count <- max(df$Count, na.rm = TRUE); if (!is.finite(max_count) || max_count <= 0) max_count <- 1

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(oldpar), add = TRUE)
    graphics::par(mar = c(5.2, 18, 4.5, 2), xpd = NA)
    graphics::plot(range(yrs), c(1, length(terms)), type = "n", axes = FALSE,
                   xlab = "Year", ylab = "",
                   main = paste0(dom, " Top10 year-count spot plot"))
    graphics::axis(1, at = yrs, labels = .sg_final_chr(as.character(yrs)), las = 2, font = 2)
    graphics::axis(2, at = seq_along(terms), labels = .sg_final_chr(terms), las = 1, font = 2, cex.axis = 0.85)
    graphics::grid(nx = NA, ny = NULL, col = "grey88")
    for (i in seq_along(terms)) {
      d <- df[df$term == terms[i], , drop = FALSE]
      d <- d[order(d$Year), , drop = FALSE]
      if (!nrow(d)) next
      graphics::lines(d$Year, rep(i, nrow(d)), col = "grey70", lwd = 1)
      cex <- 0.55 + 2.3 * sqrt(pmax(0, d$Count) / max_count)
      bg <- ifelse(d$is_burst %in% TRUE, "red", "white")
      graphics::points(d$Year, rep(i, nrow(d)), pch = 21, bg = bg, col = "black", cex = cex, lwd = 0.8)
    }
    graphics::legend("topright", legend = c("High", "Other"), pch = 21, pt.bg = c("red", "white"), bty = "n")
    invisible(NULL)
  }

  .sg_slope_extra_draw_3d <- function() {
    df <- .sg_slope_extra_long()
    dom <- .sg_slope_extra_domain()
    if (!requireNamespace("plotly", quietly = TRUE)) {
      return(plotly::plot_ly() |> plotly::layout(
        title = list(text = "Package plotly is required. Please run install.packages('plotly').")
      ))
    }
    if (is.null(df) || !is.data.frame(df) || !nrow(df)) {
      return(plotly::plot_ly() |> plotly::layout(
        title = list(text = paste0("No Top10 3D dashboard data for ", dom, ". Run analysis first or choose another Slope entity."))
      ))
    }

    terms <- unique(.sg_final_chr(df$term))

    # Use current domain FLCA nodes for maturity/influence when available.
    nodes <- tryCatch(.real_get_domain(if (dom == "Journal") "Journal/Year" else dom)$nodes, error = function(e) NULL)
    if (!is.null(nodes) && is.data.frame(nodes) && nrow(nodes) && "name" %in% names(nodes)) {
      nodes$name <- .sg_final_chr(nodes$name)
      node_df <- nodes[nodes$name %in% terms, , drop = FALSE]
      if (!"value" %in% names(node_df)) node_df$value <- NA_real_
      if (!"value2" %in% names(node_df)) node_df$value2 <- NA_real_
      node_df <- node_df[, intersect(c("name", "value", "value2"), names(node_df)), drop = FALSE]
      names(node_df)[names(node_df) == "name"] <- "term"
    } else {
      node_df <- data.frame(term = terms, value = NA_real_, value2 = NA_real_, stringsAsFactors = FALSE)
    }

    # Fallback maturity/influence from Top10 long table if FLCA node values are unavailable.
    agg_total <- stats::aggregate(Count ~ term, data = df, FUN = sum, na.rm = TRUE)
    agg_peak  <- stats::aggregate(Count ~ term, data = df, FUN = max, na.rm = TRUE)
    names(agg_total)[2] <- "total_count"
    names(agg_peak)[2]  <- "peak_count"

    met <- merge(data.frame(term = terms, stringsAsFactors = FALSE), node_df, by = "term", all.x = TRUE, sort = FALSE)
    met <- merge(met, agg_total, by = "term", all.x = TRUE, sort = FALSE)
    met <- merge(met, agg_peak,  by = "term", all.x = TRUE, sort = FALSE)
    met$value  <- .sg_final_num(met$value, NA_real_)
    met$value2 <- .sg_final_num(met$value2, NA_real_)
    met$maturity  <- ifelse(is.finite(met$value),  met$value,  .sg_final_num(met$total_count, 0))
    met$influence <- ifelse(is.finite(met$value2), met$value2, .sg_final_num(met$peak_count, 0))

    # Trend / recency = correlation coefficient between rank number 1:n and recent counts.
    # Use up to the most recent 4 timepoints for each Top10 element.
    met$trend <- vapply(met$term, function(tt) {
      d <- df[df$term == tt, , drop = FALSE]
      d <- d[order(d$Year), , drop = FALSE]
      d <- utils::tail(d, min(4, nrow(d)))
      y <- .sg_final_num(d$Count, 0)
      if (length(y) < 2 || !is.finite(stats::sd(y)) || stats::sd(y) == 0) return(0)
      r <- suppressWarnings(stats::cor(seq_along(y), y, use = "complete.obs"))
      if (is.finite(r)) r else 0
    }, numeric(1))

    # Keep exact Top10 order from the slope/long-table source.
    met$term <- .sg_final_chr(met$term)
    met <- met[nzchar(met$term), , drop = FALSE]
    if (!nrow(met)) {
      return(plotly::plot_ly() |> plotly::layout(title = list(text = "No finite Top10 3D dashboard values.")))
    }

    x <- .sg_final_num(met$maturity, 0)
    y <- .sg_final_num(met$influence, 0)
    z <- .sg_final_num(met$trend, 0)

    xmin <- min(x, na.rm = TRUE); xmax <- max(x, na.rm = TRUE)
    ymin <- min(y, na.rm = TRUE); ymax <- max(y, na.rm = TRUE)
    zmin <- min(-1, z, na.rm = TRUE); zmax <- max(1, z, na.rm = TRUE)
    if (!is.finite(xmin) || !is.finite(xmax) || xmin == xmax) { xmin <- xmin - 1; xmax <- xmax + 1 }
    if (!is.finite(ymin) || !is.finite(ymax) || ymin == ymax) { ymin <- ymin - 1; ymax <- ymax + 1 }
    if (!is.finite(zmin) || !is.finite(zmax) || zmin == zmax) { zmin <- -1; zmax <- 1 }

    # Put the vertical Trend/Z guide at the left-front corner of the maturity–influence plane.
    x_axis <- xmin
    y_axis <- ymin

    size_val <- .sg_final_num(met$maturity, 0)
    if (length(size_val) && max(size_val, na.rm = TRUE) > min(size_val, na.rm = TRUE)) {
      marker_size <- 10 + 24 * (size_val - min(size_val, na.rm = TRUE)) / (max(size_val, na.rm = TRUE) - min(size_val, na.rm = TRUE) + 1e-9)
    } else {
      marker_size <- rep(16, nrow(met))
    }

    hover_txt <- paste0(
      "<b>", htmltools::htmlEscape(met$term), "</b>",
      "<br>Maturity(value): ", round(x, 3),
      "<br>Influence(value2): ", round(y, 3),
      "<br>Trend/recency r: ", round(z, 3)
    )

    p <- plotly::plot_ly()

    # Blue maturity-influence plane at z=0.
    p <- plotly::add_trace(
      p,
      x = c(xmin, xmax, xmax, xmin),
      y = c(ymin, ymin, ymax, ymax),
      z = c(0, 0, 0, 0),
      type = "mesh3d",
      i = c(0, 0), j = c(1, 2), k = c(2, 3),
      opacity = 0.13,
      color = "lightblue",
      name = "Maturity–Influence plane",
      hoverinfo = "skip",
      showlegend = FALSE
    )

    # Red vertical Z/trend axis at left side.
    p <- plotly::add_trace(
      p,
      x = c(x_axis, x_axis),
      y = c(y_axis, y_axis),
      z = c(zmin, zmax),
      type = "scatter3d",
      mode = "lines",
      line = list(color = "red", width = 10),
      name = "Trend / recency axis",
      hoverinfo = "skip",
      showlegend = TRUE
    )

    p <- plotly::add_trace(
      p,
      x = x,
      y = y,
      z = z,
      type = "scatter3d",
      mode = "markers+text",
      text = as.character(seq_len(nrow(met))),
      textposition = "top center",
      hovertext = hover_txt,
      hoverinfo = "text",
      marker = list(
        size = marker_size,
        color = z,
        colorscale = list(c(0, "blue"), c(0.5, "white"), c(1, "red")),
        cmin = -1,
        cmax = 1,
        opacity = 0.88,
        line = list(color = "black", width = 1)
      ),
      name = "Top10 elements",
      showlegend = FALSE
    )

    # Add a compact numeric legend in the right margin.
    legend_txt <- paste0(seq_len(nrow(met)), ". ", htmltools::htmlEscape(met$term), " (r=", sprintf("%.2f", z), ")", collapse = "<br>")

    plotly::layout(
      p,
      title = list(text = paste0(dom, " Top10 interactive 3D dashboard")),
      scene = list(
        xaxis = list(title = "Maturity = value", range = c(xmin, xmax), backgroundcolor = "rgba(245,245,255,0.75)"),
        yaxis = list(title = "Influence = value2", range = c(ymin, ymax), backgroundcolor = "rgba(245,255,245,0.75)"),
        zaxis = list(
          title = list(text = "Trend / recency = cor(1:n, recent counts)", font = list(color = "red")),
          range = c(zmin, zmax),
          tickfont = list(color = "red"),
          gridcolor = "red",
          zerolinecolor = "red",
          backgroundcolor = "rgba(255,245,245,0.80)"
        ),
        aspectmode = "cube",
        camera = list(eye = list(x = 1.85, y = -2.35, z = 1.25))
      ),
      annotations = list(list(
        x = 1.02, y = 0.98, xref = "paper", yref = "paper",
        text = legend_txt,
        showarrow = FALSE,
        align = "left",
        xanchor = "left",
        yanchor = "top",
        font = list(size = 11)
      )),
      margin = list(l = 0, r = 260, b = 0, t = 55)
    )
  }

  output$slope_burst_heatmap_plot <- renderPlot({
    req(rv$pubmeta)
    .sg_slope_extra_draw_heatmap()
  }, height = 650)

  output$slope_burst_timeline_plot <- renderPlot({
    req(rv$pubmeta)
    .sg_slope_extra_draw_timeline()
  }, height = 650)

  output$slope_mesh_3d_plot <- plotly::renderPlotly({
    req(rv$pubmeta)
    .sg_slope_extra_draw_3d()
  })

  output$slope_burst_long_table <- renderPrint({
    req(rv$pubmeta)
    df <- .sg_slope_extra_long()
    if (is.null(df) || !is.data.frame(df) || !nrow(df)) {
      cat("No Top10 long table data. Run analysis first or choose another Slope entity.\n")
    } else {
      print(df, row.names = FALSE)
    }
  })

}


shinyApp(ui = ui, server = server)

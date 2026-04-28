if (file.exists(".Renviron")) {
  readRenviron(".Renviron")
}

# ---- IP module (shared gate) ----
APP_DIR <- tryCatch(normalizePath(getwd(), winslash = "/", mustWork = FALSE), error = function(e) getwd())
if (file.exists(file.path(APP_DIR, "ipmodule.R"))) {
  source(file.path(APP_DIR, "ipmodule.R"), local = FALSE)
}
PERM_PUBMED_XML <- file.path(tempdir(), "pubmed.xml") 
PERM_PUBMED_MEDLINE <- file.path(tempdir(), "pubmed_medline.txt") 
PERM_PUBMED_TXT <- file.path(tempdir(), "uploaded_pubmed_medline.txt")
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
PERM_PUBMED_TXT <- file.path(getwd(), "uploaded_pubmed_medline.txt")

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

PERM_PUBMED_TXT <- if (is_cloud()) file.path(tempdir(), "uploaded_pubmed_medline.txt") else PERM_PUBMED_TXT
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
tryCatch({ source("kano.R", local = FALSE); cat("[BOOT] sourced: kano.R
") }, error=function(e){ cat("[BOOT][ERR] kano.R:", e$message, "
") })
try({ source("sankey.R"); cat("[BOOT] sourced: sankey.R\n") }, silent=TRUE)

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
if (file.exists("renderSSplot.R")) {
  source("renderSSplot.R", local = FALSE, encoding = "UTF-8")
  cat("[BOOT] sourced: renderSSplot.R\\n")
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
          tags$div(class="small-note","After upload, the file is saved to ./uploaded_pubmed.txt automatically. Then click Fetch PubMed to run using this file."),
          tags$hr()
      ),
      div(class="card",
          h4("1) PubMed query"),
          textAreaInput("term", "Search term", value = default_term, rows = 3, width = "100%"),
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
         h4("Top10 各名稱年度次數 / Slopegraph (10年門檻 + 前5/後5 t-test；single-term only)"),
         uiOutput("yearcount_ui"),
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
         DT::DTOutput("tbl_yearcount_top20_mesh")
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
                 plotOutput("ss_kano_plot", height = "700px")
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
                 plotOutput("kano_vv2_plot", height = "650px"),
                 tags$hr(),
                 h5("Kano: SS vs a*"),
                 plotOutput("kano_ss_astar_plot", height = "650px")
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
  if (is.null(nodes) || !is.data.frame(nodes) || nrow(nodes)==0) return(nodes)
  if (is.null(edges) || !is.data.frame(edges) || nrow(edges)==0) {
    # still ensure n_byline exists if possible
    if (!("n_byline" %in% names(nodes)) && all(c("n_pubmed","single_author") %in% names(nodes))){
      nodes$n_byline <- pmax(0, as.numeric(nodes$n_pubmed) - as.numeric(nodes$single_author))
    }
    return(nodes)
  }
  nd <- nodes
  if (!("name" %in% names(nd))) {
    if ("label" %in% names(nd)) nd$name <- as.character(nd$label)
    else if ("id" %in% names(nd)) nd$name <- as.character(nd$id)
  }
  nd$name <- as.character(nd$name)

  ed <- edges
  if (!("Leader" %in% names(ed)) && "leader" %in% names(ed)) ed$Leader <- ed$leader
  if (!("Follower" %in% names(ed)) && "follower" %in% names(ed)) ed$Follower <- ed$follower
  wcol <- if ("WCD" %in% names(ed)) "WCD" else if ("wcd" %in% names(ed)) "wcd" else NULL
  if (is.null(wcol)) { ed$WCD <- 1; wcol <- "WCD" }
  ed$Leader <- as.character(ed$Leader)
  ed$Follower <- as.character(ed$Follower)
  ed[[wcol]] <- suppressWarnings(as.numeric(ed[[wcol]]))
  ed[[wcol]][is.na(ed[[wcol]])] <- 0

  fa <- aggregate(ed[[wcol]], by=list(name=ed$Leader), FUN=sum, na.rm=TRUE)
  la <- aggregate(ed[[wcol]], by=list(name=ed$Follower), FUN=sum, na.rm=TRUE)
  names(fa)[2] <- "FA"
  names(la)[2] <- "LA"
  nd <- merge(nd, fa, by="name", all.x=TRUE)
  nd <- merge(nd, la, by="name", all.x=TRUE)
  nd$FA[is.na(nd$FA)] <- 0
  nd$LA[is.na(nd$LA)] <- 0

  # Ensure single_author (self coword) exists; if absent, derive from value - value2 when possible
  if (!("single_author" %in% names(nd))) {
    if (all(c("value","value2") %in% names(nd))) nd$single_author <- pmax(0, suppressWarnings(as.numeric(nd$value) - as.numeric(nd$value2)))
    else nd$single_author <- 0
  }
  nd$single_author[is.na(nd$single_author)] <- 0

  # value2 = first+last appearances (exclude single-author self-loops)
  nd$value2 <- as.numeric(nd$FA) + as.numeric(nd$LA)

  # value = first+last appearances (include single-author papers)
  # Single-author paper counts as BOTH first and last -> +2
  nd$value  <- nd$value2 + 2 * as.numeric(nd$single_author)

  # n_pubmed = byline appearances across all authors (any position).
  # If it's not already present (or is NA), keep as NA; we will try to merge it in earlier from rv$author_lists_all.
  if (!("n_pubmed" %in% names(nd))) nd$n_pubmed <- NA_real_

  # n_byline = byline appearances (including single-author papers)
  if (!("n_byline" %in% names(nd))) {
    nd$n_byline <- nd$n_pubmed
  } else {
    nd$n_byline <- ifelse(is.na(nd$n_byline) & !is.na(nd$n_pubmed), nd$n_pubmed, nd$n_byline)
  }

  # optional: byline appearances excluding single-author papers
  if (!("n_byline_non_single" %in% names(nd))) {
    nd$n_byline_non_single <- ifelse(is.na(nd$n_pubmed), NA_real_, pmax(0, as.numeric(nd$n_pubmed) - as.numeric(nd$single_author)))
  }

  # AAC over top-3 FA values (domain-level, attached to all nodes for hover convenience)
  aac <- (.AAC_INLINE(nd$value2))
  nd$value2_strength <- aac

  nd
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

server <- function(input, output, session) {
  aac_server("aac")
  free_demo_run <- reactiveVal(TRUE)
  demo_term <- "asthma[Title/Abstract]"   # Demo 預設查詢（可改）
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

  lotka_dist_reactive <- reactive({
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
  title(main = domain_var, font.main = 2, cex.main = 1.9, line = 0.5)

  n_rows <- nrow(dom)
  text_cex <- if (n_rows <= 4) 2.4 else if (n_rows <= 7) 2.0 else 1.6
  y_positions <- seq(0.85, 0.2, length.out = n_rows)

  for (i in seq_len(n_rows)) {
    text(0.05, y_positions[i], dom$Element[i], adj = 0, cex = text_cex, font = 2)
    text(0.95, y_positions[i], dom$Value[i],   adj = 1, cex = text_cex, font = 2)
  }

  if (!is.na(aac)) {
    text(0.5, 0.08, paste("AAC =", round(aac, 4)), cex = 1.8, font = 2)
  } else {
    text(0.5, 0.08, "AAC = NA", cex = 1.8, font = 2)
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
    aac_val <- as_scalar_or_na(x$AAC$value, NA)

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
    title(main = dom, cex.main = 1.35, font.main = 2)

    if (is.null(df) || !is.data.frame(df) || !nrow(df)) {
      text(0.5, 0.55, "No data yet.\nRun this domain first.", cex = 1.05)
      mtext(paste0("AAC = ", signif(as_scalar_or_na(aac, NA), 4)), side=1, line=0.35, cex=1.0, font=2)
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

    mtext(paste0("AAC = ", signif(as_scalar_or_na(aac, NA), 4)), side=1, line=0.35, cex=1.15, font=2)
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

# ---- Upload access gate (iplist/trial/CMC) ----
g_up <- ipm_gate_session(session, cmc = input$cmc %||% "", app_dir = APP_DIR, inc_count_on_allow = TRUE)

rv$ip_access_type <- if (identical(g_up$policy, "iplist")) "IP pass" else if (identical(g_up$reason, "trial_first_time")) "Trial" else "CMC pass"
rv$ip_addr <- g_up$ip

if (!isTRUE(g_up$ok)) {
  if (identical(g_up$reason, "ip_not_allowlisted")) {
    showModal(modalDialog(
      title = "Access blocked",
      "Your IP is not in iplist.txt.",
      easyClose = TRUE,
      footer = modalButton("OK")
    ))
  } else {
    showModal(modalDialog(
      title = "CMC required",
      "Upload requires a valid 10-digit numeric CMC after the one-time trial.",
      "Please enter CMC or contact the author.",
      easyClose = TRUE,
      footer = modalButton("OK")
    ))
  }
  return()
}

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
    content = function(file){
      dir.create(file.path(RUNTIME_DIR,"download"), showWarnings = FALSE, recursive = TRUE)
      # Try best-effort export from rv
      if (!is.null(rv$nodes20)) utils::write.csv(rv$nodes20, file.path(RUNTIME_DIR,"download","nodes.csv"), row.names = FALSE, fileEncoding="UTF-8")
      if (!is.null(rv$edges20)) utils::write.csv(rv$edges20, file.path(RUNTIME_DIR,"download","edges.csv"), row.names = FALSE, fileEncoding="UTF-8")
      .write_offline_index("download", rv)
      .zip_download_dir(file, "download")
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
      nodes = if (!is.null(res$modes)) res$modes else res$nodes,
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
  demo_term <- "\"SoftwareX\"[Journal]"
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
year_vec  <- vapply(articles, extract_year, character(1))
jour_vec  <- vapply(articles, extract_journal, character(1))
keep <- nzchar(pmid_vec)
rv$pmid_tbl <- data.frame(
  PMID = pmid_vec[keep],
  Year = year_vec[keep],
  Journal = jour_vec[keep],
  Title = title_vec[keep],
  stringsAsFactors = FALSE
)

# Store per-article Year vector aligned with articles (for slopegraphs)
rv$article_years <- year_vec

      rv$pubmeta <- data.frame(
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


            # n(PubMed): number of publications per author (count of records containing the author)
      n_pubmed_tab <- table(unlist(lapply(author_lists_all, function(a){
        a <- as.character(a); a <- trimws(a); a <- a[nzchar(a)]
        unique(a)
      })))

      # n(byline appearances): raw count in byline lists (kept separate in case of duplicates)
      n_byline_tab <- table(unlist(lapply(author_lists_all, function(a){
        a <- as.character(a); a <- trimws(a); a <- a[nzchar(a)]
        a
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

      # Use the module runner: FLCA + major sampling + silhouette in ONE step
      cfg <- list(top_clusters=5, base_per_cluster=4, target_n=20, intra_delta=2, inter_delta=5, eps=1e-9)
      res <- tryCatch(run_flca_ms_sil(nodes, data_edges, cfg, verbose = TRUE), error=function(e) e)
      if (inherits(res, 'error')) {
        warning(res$message)
        nodes20 <- nodes
        if (!('carac' %in% names(nodes20))) nodes20$carac <- 1L
        nodes20$ssi <- 0
        nodes20$a_star1 <- 0
        edges20 <- data_edges
      } else {
        nodes20 <- res$modes
        edges20 <- res$data
      }

      # keep in rv for downstream tables/plots
      rv$author_nodes <- nodes20
      rv$author_edges <- edges20

      

      # Augment Author metrics (FA/LA, self, n_byline, AAC)
      # ---- byline appearances (all-author positions) + single-author count ----
      # Use rv$author_lists_all (all authors per article) created during PubMed fetch.
      if (!is.null(rv$author_lists_all) && length(rv$author_lists_all) > 0) {
        alist <- lapply(rv$author_lists_all, function(a){
          a <- unique(trimws(as.character(a)))
          a <- a[nzchar(a)]
          a
        })
        flat <- unlist(alist, use.names = FALSE)
        if (length(flat) > 0) {
          tab_all <- sort(table(flat), decreasing = TRUE)
          df_pub <- data.frame(
            name = names(tab_all),
            n_pubmed = as.integer(tab_all),
            stringsAsFactors = FALSE
          )
          # single-author papers: length==1
          singles <- unlist(lapply(alist, function(a) if (length(a)==1) a else character(0)), use.names = FALSE)
          if (length(singles) > 0) {
            tab_s <- table(singles)
            df_pub$single_author <- as.integer(tab_s[df_pub$name])
            df_pub$single_author[is.na(df_pub$single_author)] <- 0L
          } else {
            df_pub$single_author <- 0L
          }

          # merge into sampled nodes
          if (!("name" %in% names(rv$author_nodes))) {
            if ("label" %in% names(rv$author_nodes)) rv$author_nodes$name <- as.character(rv$author_nodes$label)
            else if ("id" %in% names(rv$author_nodes)) rv$author_nodes$name <- as.character(rv$author_nodes$id)
          }
          rv$author_nodes$name <- as.character(rv$author_nodes$name)
          # avoid duplicate columns from prior runs
          if ("n_pubmed" %in% names(rv$author_nodes)) rv$author_nodes$n_pubmed <- NULL
          if ("single_author" %in% names(rv$author_nodes)) rv$author_nodes$single_author <- NULL
          rv$author_nodes <- merge(rv$author_nodes, df_pub, by="name", all.x=TRUE)
          rv$author_nodes$n_pubmed[is.na(rv$author_nodes$n_pubmed)] <- 0L
          rv$author_nodes$single_author[is.na(rv$author_nodes$single_author)] <- 0L

          # defensive: if merge created suffix columns, coalesce
          if ("n_pubmed.x" %in% names(rv$author_nodes) || "n_pubmed.y" %in% names(rv$author_nodes)) {
            nx <- if ("n_pubmed.x" %in% names(rv$author_nodes)) rv$author_nodes[["n_pubmed.x"]] else NA
            ny <- if ("n_pubmed.y" %in% names(rv$author_nodes)) rv$author_nodes[["n_pubmed.y"]] else NA
            rv$author_nodes$n_pubmed <- ifelse(!is.na(ny), ny, nx)
            rv$author_nodes$n_pubmed[is.na(rv$author_nodes$n_pubmed)] <- 0L
            rv$author_nodes[["n_pubmed.x"]] <- NULL
            rv$author_nodes[["n_pubmed.y"]] <- NULL
          }
          if ("single_author.x" %in% names(rv$author_nodes) || "single_author.y" %in% names(rv$author_nodes)) {
            sx <- if ("single_author.x" %in% names(rv$author_nodes)) rv$author_nodes[["single_author.x"]] else 0L
            sy <- if ("single_author.y" %in% names(rv$author_nodes)) rv$author_nodes[["single_author.y"]] else 0L
            rv$author_nodes$single_author <- ifelse(!is.na(sy), sy, sx)
            rv$author_nodes$single_author[is.na(rv$author_nodes$single_author)] <- 0L
            rv$author_nodes[["single_author.x"]] <- NULL
            rv$author_nodes[["single_author.y"]] <- NULL
          }
        }
      }

      rv$author_nodes <- .augment_author_nodes(rv$author_nodes, rv$author_edges_full %||% rv$author_edges)
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
        "<br>n(byline appearances)=", nodes$n_byline,
        "<br>value(first/last incl. single)=", nodes$value,
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
  output$kano_aac_header <- renderText({
    d <- .get_domain(input$kano_domain)
    sprintf("AAC(value)=%.2f | AAC(value2)=%.2f | AAC(SS)=%.2f | AAC(a*)=%.2f",
            as.numeric(d$AAC_value), as.numeric(d$AAC_value2), as.numeric(d$AAC_ss), as.numeric(d$AAC_a_star))
  })
  
  output$kano_vv2_plot <- renderPlot({
    d <- .get_domain(input$kano_domain)
    nd <- d$nodes
    ed <- d$edges

    validate(need(is.data.frame(nd) && nrow(nd) > 1, "Need at least 2 nodes for Kano plot."))
    if (!("value2" %in% names(nd))) {
      if (exists("compute_value2_strength", mode="function")) {
        nd <- tryCatch(compute_value2_strength(nd, ed), error=function(e) nd)
      } else {
        plot.new(); text(0.5,0.5,"value2 not available"); return()
      }
    }

    # ensure carac exists
    if (!("carac" %in% names(nd))) nd$carac <- 1

    title_txt <- sprintf("Kano: value vs value2 | AAC=%.2f | AAC2=%.2f", as.numeric(d$AAC_value), as.numeric(d$AAC_value2))
    if (exists("plot_kano_real_xy", mode="function")) {
      plot_kano_real_xy(nd, edges = ed, xcol = "value2", ycol = "value", sizecol = "value",
                        title_txt = title_txt, label_size = input$kano_label_size %||% 4)
    } else if (exists("plot_kano_real", mode="function")) {
      print(plot_kano_real(nd, edges = ed, title_txt = title_txt))
    } else {
      plot(nd$value2, nd$value, main=title_txt, xlab="value2", ylab="value")
    }
  }, height = 650)
output$kano_ss_astar_plot <- renderPlot({
    d <- .get_domain(input$kano_domain)
    nd <- d$nodes
    validate(need(is.data.frame(nd) && nrow(nd) > 1, "Need at least 2 nodes for Kano plot."))

    # normalize columns
    if (!("ssi" %in% names(nd)) && ("SSi" %in% names(nd))) nd$ssi <- nd$SSi
    if (!("a_star1" %in% names(nd)) && ("a_star" %in% names(nd))) nd$a_star1 <- nd$a_star

    # a* fallback
    if (!("a_star1" %in% names(nd))) {
      if ("a_i" %in% names(nd)) nd$a_star1 <- 1/(1+as.numeric(nd$a_i)) else nd$a_star1 <- NA_real_
    }

    # carac safeguard (Kano core requires it)
    if (!("carac" %in% names(nd))) nd$carac <- 1L
    nd$carac <- suppressWarnings(as.integer(as.character(nd$carac)))
    nd$carac[is.na(nd$carac)] <- 1L

    nd$ssi     <- suppressWarnings(as.numeric(nd$ssi))
    nd$a_star1 <- suppressWarnings(as.numeric(nd$a_star1))

    ok <- is.finite(nd$ssi) & is.finite(nd$a_star1)
    validate(need(sum(ok) >= 2, "Need >=2 finite points for SS vs a* Kano."))

    aac_ss  <- if (!is.null(d$AAC_ss)) as.numeric(d$AAC_ss) else safe_AAC(nd$ssi)
    aac_ast <- if (!is.null(d$AAC_a_star)) as.numeric(d$AAC_a_star) else if (!is.null(d$AAC_a_star1)) as.numeric(d$AAC_a_star1) else safe_AAC(nd$a_star1)
    title_txt <- sprintf("Kano: SS vs a* | AAC(SS)=%.2f | AAC(a*)=%.2f", aac_ss, aac_ast)

    tryCatch({
      if (exists("plot_kano_real_xy", mode="function")) {
        # true Kano: x=a*, y=SS
        plot_kano_real_xy(nd[ok, , drop=FALSE], edges = NULL,
                          xcol = "a_star1", ycol = "ssi", sizecol = "ssi",
                          title_txt = title_txt, label_size = input$kano_label_size %||% 4,
                          xlab="a*", ylab="SS")
      } else {
        plot(nd$a_star1[ok], nd$ssi[ok], xlab="a*", ylab="SS", main=title_txt)
        text(nd$a_star1[ok], nd$ssi[ok], labels=nd$name[ok], pos=4, cex=0.8)
      }
    }, error = function(e){
      plot(nd$a_star1[ok], nd$ssi[ok], xlab="a*", ylab="SS", main=title_txt)
      text(nd$a_star1[ok], nd$ssi[ok], labels=nd$name[ok], pos=4, cex=0.8)
    })
  }, height = 650)

  # ---- SS Kano tab (real Kano: 2 wings + 2 circles; bubble=value; color=cluster) ----
  output$ss_kano_plot <- renderPlot({
    dom <- input$ss_kano_domain %||% "Author"
    d <- .get_domain(dom)
    nd <- d$nodes
    validate(need(is.data.frame(nd) && nrow(nd) > 1, "Need at least 2 nodes."))

    # normalize columns
    if (!("ssi" %in% names(nd)) && ("SSi" %in% names(nd))) nd$ssi <- nd$SSi
    if (!("a_star1" %in% names(nd)) && ("a_star" %in% names(nd))) nd$a_star1 <- nd$a_star
    if (!("a_star1" %in% names(nd)) && ("a_i" %in% names(nd))) nd$a_star1 <- 1/(1+as.numeric(nd$a_i))

    nd$ssi     <- suppressWarnings(as.numeric(nd$ssi))
    nd$a_star1 <- suppressWarnings(as.numeric(nd$a_star1))
    nd$value   <- suppressWarnings(as.numeric(nd$value))
    if (!("carac" %in% names(nd))) nd$carac <- "C1"
    nd$carac <- as.character(nd$carac)  # keep label; do NOT coerce to integer

    ok <- is.finite(nd$ssi) & is.finite(nd$a_star1)
    validate(need(sum(ok) >= 2, "Need >=2 finite points."))

    aac_ss  <- safe_AAC(nd$ssi[ok])
    aac_ast <- safe_AAC(nd$a_star1[ok])
    title_txt <- sprintf("SS Kano: SS vs a* | AAC(SS)=%.2f | AAC(a*)=%.2f", aac_ss, aac_ast)

    if (exists("plot_kano_real_xy", mode="function")) {
      plot_kano_real_xy(nd[ok, , drop=FALSE], edges = NULL,
                        xcol="a_star1", ycol="ssi", sizecol="value",
                        title_txt = title_txt, label_size = input$ss_kano_label_size %||% 4,
                        xlab="a*", ylab="SS")
    } else if (exists("plot_kano_real", mode="function")) {
      # fallback: rename to canonical value/value2 so plot_kano_real can draw wings/circles
      tmp <- nd[ok, , drop=FALSE]
      tmp$value2 <- tmp$a_star1
      tmp$value  <- tmp$ssi
      tmp$size_plot <- pmax(1, tmp$value)  # bubble by value
      print(plot_kano_real(tmp, edges = NULL, title_txt = title_txt))
    } else {
      plot(nd$a_star1[ok], nd$ssi[ok], xlab="a*", ylab="SS", main=title_txt)
      text(nd$a_star1[ok], nd$ssi[ok], labels=nd$name[ok], pos=4, cex=0.8)
    }
  }, height = 700)


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
      top <- .get_domain_top20(rv, d, n = 10)
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
      top <- .get_domain_top20(rv, d, n = 10)
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
        nodes0 = nd,
        results = results,
        nodes = nd,
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
        nodes0 = nd,
        nodes = nd,
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
content  = function(file){
  writeLines(make_summary_html(), con = file, useBytes = TRUE)
}
)
 # PNG dashboard download
output$download_summary_png <- downloadHandler(
  filename = function(){ sprintf("summary_%s.png", format(Sys.time(), "%Y%m%d_%H%M%S")) },
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
  # Author byline (all authors)
  if (d %in% c("author","authors","author_byline","author(byline)","author (byline)")) return(rv$author_lists_all %||% rv$author_lists)
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
    if (is.null(tl)) next

    # align lengths with years
    n <- min(length(tl), length(yrs))
    tl <- tl[seq_len(n)]
    yv <- yrs[seq_len(n)]

    # choose top20 items
    top20 <- .get_domain_top20(rv, dom, n = 10)
    if (!length(top20)) {
      top20 <- utils::head(unique(unlist(lapply(tl, .split_terms2))), 20)
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
          terms <- .split_terms2(tl[[i]])
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
    n_byline = .getnum(top, "n_byline"),
    value2 = .getnum(top, "value2"),
    self_coword = .getnum(top, "single_author"),
    value = .getnum(top, "value"),
    FA = .getnum(top, "FA"),
    LA = .getnum(top, "LA"),
    value2_strength = .getnum(top, "value2_strength"),
    AAC = as.numeric(aac),
    check.names = FALSE
  )
  out
}, striped = TRUE, bordered = TRUE, hover = TRUE, spacing = "s")


}


`%||%` <- function(a, b) if (is.null(a)) b else a

shinyApp(ui = ui, server = server)

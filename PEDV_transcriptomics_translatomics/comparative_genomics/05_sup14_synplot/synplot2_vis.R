#!/usr/bin/env Rscript
# synplot2_vis.R
#
# Usage:
#   Rscript synplot2_vis.R [prefix] [window] [alignment.fas]
#
# {prefix}_main.pdf   — synplot + per-sequence codon heatmap (3 frames)
# {prefix}_pairs.pdf  — per-pair diagnostics

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

args     <- commandArgs(trailingOnly = TRUE)
prefix   <- if (length(args) >= 1) args[1] else "frasynplot"
window   <- if (length(args) >= 2) as.integer(args[2]) else 25L
fas_path <- if (length(args) >= 3) args[3] else NULL

# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

read_plot <- function(prefix) {
  f <- sprintf("%s_plot.txt", prefix)
  if (!file.exists(f)) stop("Cannot find: ", f)
  read.table(f, comment.char="#", header=FALSE,
             col.names=c("codon","obs","exp","stddev"))
}

read_data <- function(prefix) {
  f <- sprintf("%s_data.txt", prefix)
  if (!file.exists(f)) stop("Cannot find: ", f)
  read.table(f, comment.char="#", header=FALSE,
             col.names=c("seq_a","seq_b","codlen","synlen",
                         "nmuts","ntrans","t","expnmuts","expntrans","flag"))
}

roll_mean <- function(x, k) {
  h <- floor(k/2); n <- length(x)
  vapply(seq_len(n), function(i) mean(x[max(1,i-h):min(n,i+h)], na.rm=TRUE), numeric(1))
}

sliding_pval <- function(obs, exp, stddev, w) {
  n <- length(obs); pv <- rep(NA_real_, n); var <- stddev^2; h <- floor(w/2)
  for (i in seq_len(n)) {
    lo <- max(1,i-h); hi <- min(n,i+h)
    mv <- mean(var[lo:hi], na.rm=TRUE)
    if (is.na(mv) || mv <= 0) next
    pv[i] <- pnorm(sqrt(w) * mean(obs[lo:hi]-exp[lo:hi], na.rm=TRUE) / sqrt(mv))
  }
  pv
}

# ─────────────────────────────────────────────────────────────────────────────
# FASTA reader + per-sequence codon scanner
# ─────────────────────────────────────────────────────────────────────────────

read_fasta <- function(path) {
  raw   <- readLines(path)
  raw   <- gsub("\r","", raw)
  heads <- which(startsWith(raw, ">"))
  ends  <- c(heads[-1]-1L, length(raw))
  names <- sub("^>\\s*","", raw[heads]); names <- sub("\\s.*","", names)
  seqs  <- vapply(seq_along(heads), function(i) {
    bl <- raw[(heads[i]+1L):ends[i]]
    toupper(paste(bl[nchar(bl)>0], collapse=""))
  }, character(1))
  list(names=names, seqs=seqs)
}

# Scan one sequence in one frame.
# Returns integer vector of 1-based nt positions for each codon type.
scan_frame <- function(seq, frame) {
  seq <- chartr("U","T", seq)
  L   <- nchar(seq)
  i   <- frame + 1L
  start_pos <- integer(); stop_pos <- integer()
  while (i + 2L <= L) {
    cdn <- substr(seq, i, i+2L)
    if (!grepl("[-.]", cdn)) {
      if (cdn == "ATG")                      start_pos <- c(start_pos, i)
      else if (cdn %in% c("TAA","TAG","TGA")) stop_pos  <- c(stop_pos,  i)
    }
    i <- i + 3L
  }
  list(start=start_pos, stop=stop_pos)
}

# Build long data.frame for the heatmap:
#   seq_idx  row index (1 = first seq in file)
#   frame    0 / 1 / 2
#   codon_x  position in frame-0 codon units  (= (nt_pos-1)/3 + 1)
#   type     "start" | "stop"
build_codon_df <- function(fas_path) {
  message("  Reading FASTA: ", fas_path)
  fa <- read_fasta(fas_path)
  n  <- length(fa$seqs)
  L  <- nchar(fa$seqs[1])
  message("  ", n, " sequences, ", L, " nt")

  rows <- list()
  for (si in seq_len(n)) {
    for (fr in 0:2) {
      hits <- scan_frame(fa$seqs[si], fr)
      for (tp in c("start","stop")) {
        pos <- hits[[tp]]
        if (length(pos) == 0) next
        rows[[length(rows)+1]] <- data.frame(
          seq_idx = si,
          frame   = fr,
          codon_x = (pos - 1) / 3 + 1,   # align to frame-0 codon axis
          type    = tp,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  df        <- do.call(rbind, rows)
  df$n_seqs <- n
  df$seq_names <- fa$names[df$seq_idx]
  df
}

# ─────────────────────────────────────────────────────────────────────────────
# Codon heatmap — one panel per reading frame
# ─────────────────────────────────────────────────────────────────────────────

make_heatmap <- function(codon_df, fr, n_codons, show_x = FALSE) {
  n_seqs <- codon_df$n_seqs[1]
  df     <- codon_df[codon_df$frame == fr, ]

  # Tile width in codon units: show each codon as 1 unit wide
  p <- ggplot(df, aes(x=codon_x, y=seq_idx, fill=type)) +
    geom_tile(width=1, height=1, colour=NA) +
    scale_fill_manual(values = c("start"="#27ae60", "stop"="#7f8c8d"),
                      labels = c("start"="ATG (start)", "stop"="Stop codon"),
                      name   = NULL,
                      drop   = FALSE) +
    scale_x_continuous(limits  = c(0.5, n_codons+0.5),
                       expand  = expansion(0),
                       breaks  = seq(0, n_codons, by=500)) +
    scale_y_continuous(limits  = c(0.5, n_seqs+0.5),
                       expand  = expansion(0),
                       breaks  = NULL) +
    labs(y = sprintf("Frame +%d", fr+1),
         x = if (show_x) sprintf("Codon position (frame-0 units)") else NULL) +
    theme_classic(base_size=9) +
    theme(
      axis.line      = element_blank(),
      axis.ticks.y   = element_blank(),
      axis.text.y    = element_blank(),
      axis.title.y   = element_text(angle=0, vjust=0.5, size=8, colour="grey30",
                                    margin=margin(r=4)),
      axis.text.x    = if (show_x) element_text(size=7) else element_blank(),
      axis.ticks.x   = if (show_x) element_line() else element_blank(),
      axis.title.x   = element_text(size=8),
      panel.background = element_rect(fill="grey97", colour="grey75"),
      panel.border   = element_rect(fill=NA, colour="grey75"),
      legend.position = "none",
      plot.margin    = margin(0,8,0,8)
    )
  p
}

# ─────────────────────────────────────────────────────────────────────────────
# Main synplot panels
# ─────────────────────────────────────────────────────────────────────────────

make_synplot_panels <- function(prefix, window) {
  d    <- read_plot(prefix)
  n    <- nrow(d)
  pval <- sliding_pval(d$obs, d$exp, d$stddev, window)
  bonf <- 0.05 / max(sum(!is.na(pval)), 1)

  ratio <- ifelse(roll_mean(d$exp,window) > 0,
                  roll_mean(d$obs,window) / roll_mean(d$exp,window), NA_real_)

  pval[pval<=0] <- NA_real_
  nlp <- -log10(pval)
  xbr <- seq(0, n, by=500)

  p_pval <- ggplot(data.frame(x=d$codon, y=nlp), aes(x,y)) +
    geom_line(colour="#c0392b", linewidth=0.55, na.rm=TRUE) +
    geom_hline(yintercept=-log10(bonf), linetype="dashed",
               colour="grey50", linewidth=0.45) +
    annotate("text", x=n*0.99, y=-log10(bonf),
             label=sprintf("Bonferroni (%.1e)",bonf),
             hjust=1, vjust=-0.4, size=2.5, colour="grey40") +
    scale_x_continuous(breaks=xbr, expand=expansion(0)) +
    scale_y_continuous(expand=expansion(mult=c(0,0.05))) +
    labs(y=expression(-log[10](p))) +
    theme_classic(base_size=11) +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
          axis.ticks.x=element_line(),
          panel.grid.major.x=element_line(colour="grey93",linewidth=0.3),
          plot.margin=margin(4,8,1,8))

  ratio_max <- max(ceiling(max(ratio,na.rm=TRUE)*2)/2, 2)
  p_ratio <- ggplot(data.frame(x=d$codon, y=ratio), aes(x,y)) +
    geom_line(colour="#6c3483", linewidth=0.55, na.rm=TRUE) +
    geom_hline(yintercept=1, linetype="dashed", colour="grey40", linewidth=0.45) +
    scale_y_continuous(limits=c(0,ratio_max), breaks=seq(0,ratio_max,0.5),
                       expand=expansion(mult=c(0,0.02))) +
    scale_x_continuous(breaks=xbr, expand=expansion(0)) +
    labs(x=NULL, y="obs / exp") +
    theme_classic(base_size=11) +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
          axis.ticks.x=element_line(),
          panel.grid.major.x=element_line(colour="grey93",linewidth=0.3),
          plot.margin=margin(1,8,1,8))

  list(p_pval=p_pval, p_ratio=p_ratio, n=n, bonf=bonf)
}

# ─────────────────────────────────────────────────────────────────────────────
# Assemble + save main figure
# ─────────────────────────────────────────────────────────────────────────────

plot_main <- function(prefix, window, fas_path=NULL) {
  sp <- make_synplot_panels(prefix, window)
  n  <- sp$n

  if (!is.null(fas_path) && file.exists(fas_path)) {
    message("Building codon heatmaps ...")
    cdf    <- build_codon_df(fas_path)
    n_seqs <- cdf$n_seqs[1]

    hm0 <- make_heatmap(cdf, 0, n, show_x=FALSE)
    hm1 <- make_heatmap(cdf, 1, n, show_x=FALSE)
    hm2 <- make_heatmap(cdf, 2, n, show_x=TRUE)

    # Height of each heatmap row scales with number of sequences
    hm_h <- max(0.05 * n_seqs, 0.5)   # inches per heatmap

    fig <- sp$p_pval / sp$p_ratio / hm0 / hm1 / hm2 +
      plot_layout(heights = c(1.8, 1, hm_h, hm_h, hm_h),
                  guides  = "collect") +
      plot_annotation(
        title    = sprintf("Synonymous site conservation  |  %s", prefix),
        subtitle = sprintf(
          "window = %d codons  |  Bonferroni = %.2e  |  heatmap: %d seqs × 3 frames",
          window, sp$bonf, n_seqs),
        theme = theme(
          plot.title    = element_text(size=12, face="bold"),
          plot.subtitle = element_text(size=8, colour="grey40"),
          legend.position = "bottom",
          legend.key.size = unit(0.35,"cm"),
          legend.text     = element_text(size=9)
        )
      )

    fig_h <- 5.5 + 3 * hm_h + 0.8   # synplot + 3 heatmaps + legend

  } else {
    if (!is.null(fas_path)) warning("FASTA not found: ", fas_path)
    fig <- sp$p_pval / sp$p_ratio +
      plot_layout(heights=c(1.6,1)) +
      plot_annotation(
        title    = sprintf("Synonymous site conservation  |  %s", prefix),
        subtitle = sprintf("window = %d  |  Bonferroni = %.2e", window, sp$bonf),
        theme = theme(plot.title=element_text(size=12,face="bold"),
                      plot.subtitle=element_text(size=9,colour="grey40"))
      )
    fig_h <- 5.5
  }

  out   <- sprintf("%s_main.pdf", prefix)
  fig_w <- max(8, 0.5 + n/100)
  ggsave(out, fig, width=fig_w, height=fig_h, device="pdf", limitsize=FALSE)
  message(sprintf("Wrote: %s  (%.1f × %.1f in)", out, fig_w, fig_h))
}

# ─────────────────────────────────────────────────────────────────────────────
# Per-pair diagnostics
# ─────────────────────────────────────────────────────────────────────────────

plot_pairs <- function(prefix) {
  d <- read_data(prefix)
  d$pair_idx  <- seq_len(nrow(d))
  d$converged <- factor(ifelse(d$flag==0,"OK","flag=1"), levels=c("OK","flag=1"))

  pA <- ggplot(d, aes(pair_idx, t, colour=converged)) +
    geom_point(size=1.8, alpha=0.85) +
    geom_segment(aes(xend=pair_idx, yend=0), linewidth=0.35, alpha=0.5) +
    scale_colour_manual(values=c("OK"="#2980b9","flag=1"="#e74c3c"), name=NULL) +
    labs(x="Pair index", y="t", title="Fitted t per sequence pair") +
    theme_classic(base_size=10) +
    theme(legend.position="top",
          panel.grid.major.y=element_line(colour="grey93",linewidth=0.3))

  dd  <- subset(d, synlen>0)
  lim <- max(dd$nmuts, dd$expnmuts, na.rm=TRUE) * 1.05
  pB  <- ggplot(dd, aes(expnmuts, nmuts, colour=converged)) +
    geom_abline(slope=1, intercept=0, linetype="dashed", colour="grey55") +
    geom_point(size=1.8, alpha=0.85) +
    scale_colour_manual(values=c("OK"="#27ae60","flag=1"="#e74c3c"), name=NULL) +
    coord_fixed(xlim=c(0,lim), ylim=c(0,lim)) +
    labs(x="Expected mutations", y="Observed mutations",
         title="Obs vs Exp synonymous mutations") +
    theme_classic(base_size=10) +
    theme(legend.position="top",
          panel.grid.major=element_line(colour="grey93",linewidth=0.3))

  pC <- ggplot(d, aes(pair_idx, synlen)) +
    geom_col(fill="#1a5276", width=0.7, alpha=0.8) +
    geom_hline(yintercept=median(d$synlen,na.rm=TRUE),
               linetype="dashed", colour="#e67e22", linewidth=0.5) +
    annotate("text", x=nrow(d)*0.99, y=median(d$synlen,na.rm=TRUE),
             label=sprintf("median=%d", as.integer(median(d$synlen))),
             hjust=1, vjust=-0.4, size=2.8, colour="#e67e22") +
    labs(x="Pair index", y="Synonymous sites",
         title="Synonymous sites per pair") +
    theme_classic(base_size=10) +
    theme(panel.grid.major.y=element_line(colour="grey93",linewidth=0.3))

  out <- sprintf("%s_pairs.pdf", prefix)
  ggsave(out,
    (pA | pB) / pC +
      plot_layout(heights=c(1,1)) +
      plot_annotation(
        title    = sprintf("Per-pair diagnostics  |  %s", prefix),
        subtitle = sprintf("%d pairs  |  %d flagged", nrow(d), sum(d$flag!=0)),
        theme = theme(plot.title=element_text(size=12,face="bold"),
                      plot.subtitle=element_text(size=9,colour="grey40"))
      ),
    width=12, height=8, device="pdf")
  message("Wrote: ", out)
}

# ─────────────────────────────────────────────────────────────────────────────
# Run
# ─────────────────────────────────────────────────────────────────────────────

message(sprintf("prefix='%s'  window=%d  fasta=%s",
                prefix, window, if(is.null(fas_path)) "(none)" else fas_path))
plot_main(prefix, window, fas_path)
plot_pairs(prefix)
message("Done.")

# UC CosMx 6K: mixed-effects logistic regression of SPP1, MMP9 and receptor detection in SELENOP+ macrophages on distance
# to the nearest (ligand-expressing) OGN+RSPO3+ fibroblast (Fig 3H; plotted by 12).
# Input:  Fig2H_cosmx_fib_mac_distance/intermediate/ligand_sender_cells.csv.gz (07)
# Output: intermediate/glmm_final.csv (read by 12)
# Run:    Rscript <file>, or RUN_SECTIONS=<id>[,<id>] Rscript <file>
.RUN_SECTIONS <- Filter(nzchar, trimws(strsplit(Sys.getenv("RUN_SECTIONS", ""), ",")[[1]]))
.run <- function(section) {
  if (length(.RUN_SECTIONS) > 0 && !(section %in% .RUN_SECTIONS)) return(FALSE)
  message("===== section ", section, " ====="); TRUE
}


# Section F3H_distance_glmm: binomial GLMMs per outcome and model (UC, HC, ALL, INT); BH within family per model
if (.run("F3H_distance_glmm")) {
# CosMx distance models. Sender definition 'union-both': OGN+RSPO3+ fibroblast =
# fine label OGN+RSPO3+ OR any fibroblast detecting BOTH OGN and RSPO3; ligand sender additionally expresses the ligand.
# Columns dub_* from 07_mac_fib_spatial_input.py. Outcomes: SPP1, MMP9; IL33 IL1RL1, IL1RAP; CXCL12 CXCR4, ACKR3 (ACKR3 is a
# decoy); FGF7 FGFR1, FGFR2, NRP1; IGF1 IGF1R, IGF2R - each receptor modelled on its own; response program for all 4 ligands.
# Models (binomial glmer; x = distance to nearest sender / 50 um, cap 100; o = distance to nearest other fibroblast / 50):
#   UC : y ~ x + o + lt + (1 + x || patient) + (1 | patient:sample)          (8 UC patients; main)
#   HC : y ~ x + o + lt + (1 | patient) + (1 | patient:sample)               (3 healthy controls; descriptive)
#   ALL: y ~ x + condition + o + lt + (1 + x || patient) + (1 | patient:sample)  (all 11; UC + healthy combined)
#   INT: y ~ x * condition + o + lt + (1 + x || patient) + (1 | patient:sample)  (all 11; x:condition = UC vs healthy difference)
# BH-FDR within family (SMS markers; each ligand's displayed outcomes) per model. Output intermediate/glmm_final.csv
suppressPackageStartupMessages(library(lme4))
O <- '/share/fsmresfiles/UC/revision_jenny/Fig2H_cosmx_fib_mac_distance'
D <- read.csv(gzfile(file.path(O,'intermediate/ligand_sender_cells.csv.gz'))); D$lt <- as.numeric(scale(D$log_total))
D$condition <- factor(D$condition, levels=c('Control','UC'))
FAM <- list('SMS markers'=list(x='dub_ogn', o='dub_otherfib', y=c('SPP1','MMP9')),
            IL33=list(x='dub_IL33', o='dub_other_IL33', y=c('IL1RL1','IL1RAP','prog_IL33')),
            CXCL12=list(x='dub_CXCL12', o='dub_other_CXCL12', y=c('CXCR4','ACKR3','prog_CXCL12')),
            FGF7=list(x='dub_FGF7', o='dub_other_FGF7', y=c('FGFR1','FGFR2','NRP1','prog_FGF7')),
            IGF1=list(x='dub_IGF1', o='dub_other_IGF1', y=c('IGF1R','IGF2R','prog_IGF1')))
ctl <- glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=2e5))
fit <- function(d, y, xc, oc, model){
  d$x <- d[[xc]]/50; d$o <- d[[oc]]/50; prog <- startsWith(y,'prog_')
  lhs <- if (prog) sprintf('cbind(%s_k, %s_n - %s_k)', y, y, y) else y
  rhs <- switch(model, UC='x + o + lt + (1 + x || patient) + (1 | patient:sample)',
                       HC='x + o + lt + (1 | patient) + (1 | patient:sample)',
                       ALL='x + condition + o + lt + (1 + x || patient) + (1 | patient:sample)',
                       INT='x * condition + o + lt + (1 + x || patient) + (1 | patient:sample)')
  m <- tryCatch(glmer(as.formula(paste(lhs,'~',rhs)), data=d, family=binomial, control=ctl), error=function(e) NULL)
  if (is.null(m) && model!='HC') m <- tryCatch(glmer(as.formula(paste(lhs,'~',sub('\\(1 \\+ x \\|\\| patient\\)','(1 | patient)',rhs))), data=d, family=binomial, control=ctl), error=function(e) NULL)
  if (is.null(m)) return(c(NA,NA,NA))
  cf <- summary(m)$coefficients; term <- if (model=='INT') 'x:conditionUC' else 'x'
  c(cf[term,'Estimate'], cf[term,'Std. Error'], cf[term,'Pr(>|z|)'])
}
out <- list()
for (model in c('UC','HC','ALL','INT')){
  d <- switch(model, UC=D[D$condition=='UC',], HC=D[D$condition=='Control',], D)
  for (fam in names(FAM)) for (y in FAM[[fam]]$y){
    r <- fit(d, y, FAM[[fam]]$x, FAM[[fam]]$o, model)
    out[[length(out)+1]] <- data.frame(model=model, family=fam, outcome=y, beta=r[1], se=r[2], p=r[3],
      OR=exp(-r[1]), OR_lo=exp(-(r[1]+1.96*r[2])), OR_hi=exp(-(r[1]-1.96*r[2])), patients=length(unique(d$patient)), cells=nrow(d))
  }
}
R <- do.call(rbind,out); R$q <- ave(R$p, R$model, R$family, FUN=function(p) p.adjust(p,'BH'))
write.csv(R, file.path(O,'intermediate/glmm_final.csv'), row.names=FALSE)
}  # end section F3H_distance_glmm

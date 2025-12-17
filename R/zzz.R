utils::globalVariables(c(
  "CV", "Index", "count", "dnorm", "group", "highlight_label", "is",
  "label", "label_y", "logFC", "mean_log2_intensity", "movingSD", "sd",
  "setNames", "sig.FDR.fct", "sig.pval.fct", "status", "terms", "tot.int",
  "tot.num", "value", "values"
))

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("A", "M", "deltaM", "method", "residual"))
}

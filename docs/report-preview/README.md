# Synthetic report preview

Open `mock-report.html` in a browser to inspect the workflow report layout. All data are invented, and the preview is prominently labelled. The source copies `assets/templates/report_template.qmd` and adds representative sections using the workflow's table, figure and tab patterns. It preserves the existing report title, author, theme, CSS and table of contents. Plot content and section selection are illustrative.

The HTML embeds its assets and can be opened without a server. This preview was rendered with local Quarto 1.7.31 and R, rather than the pinned production container.

To regenerate with Quarto and the R packages `knitr`, `rmarkdown` and `kableExtra` installed:

```bash
quarto render docs/report-preview/mock-report.qmd --to html --no-cache
```

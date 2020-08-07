### Downstream analysis

I used methylKit for the downstream analysis. After bismark and sorting, bam files were fed into ```processBismarkAln``` and ```getCoverageStats```, ```getMethylationStats``` functions for statistics, which generates histograms of % of methylation and coverage.

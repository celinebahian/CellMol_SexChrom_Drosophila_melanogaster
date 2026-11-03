# CellMol_SexChrom_Drosophila_melanogaster

**Name:** Bahian, Celine R.

**Section:** B

**Selected species:** **Drosophila melanogaster*

**Study question:** How can sequencing read depth across the autosomes and sex chromosomes be used to infer the possible chromosomal sex of a Drosophila melanogaster WGS sample?

**Assembly accession:** GCF_000001215.4 https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001215.4/

**Assembly BioSample accession:** SAMN02803731 https://www.ncbi.nlm.nih.gov/biosample/SAMN02803731

**WGS accession (Run):** SRR6702604 https://www.ncbi.nlm.nih.gov/sra/SRR6702604

**WGS BioSample accession:** SAMN08511563 https://www.ncbi.nlm.nih.gov/biosample/SAMN08511563

## Galaxy Steps

<img width="1362" height="640" alt="image" src="https://github.com/user-attachments/assets/09377405-3c64-478c-b8d2-30f52f9070df" />

**Figure 1. Uploading the reference genome assembly files to Galaxy.** Two files associated with the **Drosophila melanogaster* reference genome assembly (GCF_000001215.4) were added to the upload queue in Galaxy before proceeding with the genomic analysis.

<img width="1073" height="630" alt="image" src="https://github.com/user-attachments/assets/711cb65a-59aa-4f70-be0a-ba064404085d" />

**Figure 2. Downloading WGS reads from NCBI SRA using Faster Download and Extract Reads in FASTQ format.** The SRA accession SRR6702604 was entered into Galaxy to retrieve the whole-genome sequencing reads in FASTQ format

<img width="1078" height="637" alt="image" src="https://github.com/user-attachments/assets/6a1fe52c-f82a-45d5-8ecf-cafe0d3f322a" />

**Figure 3. Performing read quality assessment using FastQC.** The downloaded paired-end sequencing reads were selected as input for FastQC to evaluate their quality before proceeding with the mapping process.


<img width="1075" height="633" alt="image" src="https://github.com/user-attachments/assets/8985f28c-a1d2-4741-8479-b6f598e28fbd" />

**Figure 4. Separating paired-end sequencing reads using Unzip Collection.** The paired-end FASTQ collection was separated into forward and reverse read datasets for individual processing during genome mapping.


<img width="1079" height="634" alt="image" src="https://github.com/user-attachments/assets/558e4cb7-2351-4261-bfd7-318818358ed9" />

**Figure 5. Building the reference genome index using BWA-MEM2.** The **Drosophila melanogaster* reference genome assembly (GCF_000001215.4) was selected to build an index required for read alignment.


<img width="1142" height="637" alt="image" src="https://github.com/user-attachments/assets/33d55fa3-9a3e-4617-9b0c-67294dad0576" />

**Figure 6. Mapping paired-end reads to the reference genome using BWA-MEM2.** The forward and reverse sequencing reads were aligned against the **Drosophila melanogaster* reference genome to determine their corresponding genomic positions.


<img width="1135" height="632" alt="image" src="https://github.com/user-attachments/assets/2830e24b-d133-4e6f-b6b5-e5c385459b3f" />

**Figure 7. Calculating chromosome coverage using Samtools Coverage.** The mapped BAM file was used as input to calculate sequencing coverage across the reference chromosomes. The minimum mapping quality and minimum base quality were both set to 20 to filter low-quality alignments and bases.

https://usegalaxy.org/u/celinebahian/h/sex-dosage-drosophila-melanogaster

### Chromosome Groups

| Chromosome Group | Accession Number | Reason for Selection |
|---|---|---|
| Autosomes (A) | NC_004353.4 | Used as a reference for calculating relative sequencing depth. |
| Autosomes (A) | NT_033777.3 | Included in the autosomal baseline for comparison. |
| Autosomes (A) | NT_033778.4 | Included in the autosomal baseline for comparison. |
| Autosomes (A) | NT_033779.5 | Included in the autosomal baseline for comparison. |
| Autosomes (A) | NT_037436.4 | Included in the autosomal baseline for comparison. |
| X Chromosome | NC_004354.4 | Used to calculate the X:A ratio and assess X-chromosome representation. |
| Y Chromosome | NC_024512.1 | Used to calculate the Y:A ratio and assess Y-chromosome representation. |

Five autosomal reference sequences were selected for the autosomal depth calculation. These include chromosome 4 and the left and right arms of chromosomes 2 and 3. The X chromosome was used to calculate the X:A ratio, while the Y chromosome was included to examine Y-chromosome read representation. 

## Complete Rcode

#Sex-chromosome dosage from whole-genome DNA read depth
#Input: Samtools coverage TSV + student chromosome labels
cov <- read.delim("chromosome_coverage.tsv", comment.char = "",
                  check.names = FALSE, stringsAsFactors = FALSE)
names(cov)[1] <- "rname" # Accepts the '#rname' header
if (!"meandepth" %in% names(cov)) stop("Missing meandepth column")
key <- read.csv("chromosome_groups.csv", stringsAsFactors = FALSE)
if (!all(c("rname", "group") %in% names(key))) stop("Bad key columns")
key$rname <- trimws(key$rname)
key$group <- toupper(trimws(key$group))
if (anyDuplicated(key$rname)) stop("Repeated rname in key file")
if (!all(key$group %in% c("A", "X", "Y", "Z", "W")))
  stop("Use only A, X, Y, Z or W in the group column")

x <- merge(key, cov[, c("rname", "meandepth")],
           by = "rname", all.x = TRUE, sort = FALSE)
x$meandepth <- suppressWarnings(as.numeric(x$meandepth))
if (anyNA(x$meandepth)) stop("Missing depths: check rname matches")
if (any(x$meandepth < 0)) stop("Negative depth is not valid")
if (sum(x$group == "A") < 5) stop("Label at least 5 autosomes")
if ("X" %in% x$group && "Z" %in% x$group)
  stop("Choose either an XY or a ZW species, not both")
if (sum(x$group == "X") == 1 && !any(x$group == "Z")) {
  target <- "X"
  one_copy <- "XY-like"
  two_copies <- "XX-like"
} else if (sum(x$group == "Z") == 1 && !any(x$group == "X")) {
  target <- "Z"
  one_copy <- "ZW-like"
  two_copies <- "ZZ-like"
} else stop("Label exactly one X or exactly one Z chromosome")

A <- median(x$meandepth[x$group == "A"])
if (!is.finite(A) || A <= 0) stop("No usable autosomal coverage")
x$relative_depth <- x$meandepth / A
ratio <- x$relative_depth[x$group == target]
call <- if (ratio >= 0.80) two_copies else if (ratio <= 0.65)
  one_copy else "Ambiguous"
if (A < 1) call <- "Inconclusive: very low autosomal depth"

#Additional Y or W evidence can support, but not prove, a call
partner <- if (target == "X") "Y" else "W"
partner_ratio <- if (partner %in% x$group)
  x$relative_depth[x$group == partner][1] else NA_real_

#Export results and a simple, reproducible figure
dir.create("results", showWarnings = FALSE)
write.csv(x, "results/chromosome_relative_depth.csv", row.names = FALSE)
png("results/chromosome_depth.png", width = 1000, height = 580)
par(mar = c(8, 4, 3, 1))
barplot(x$relative_depth, names.arg = x$rname, las = 2,
        ylab = "Depth relative to median autosome",
        main = "Autosomal and sex-chromosome depth", cex.names = 0.75)
abline(h = c(0.5, 1.0), lty = 2)
dev.off()

summary <- c(paste("Median autosome depth:", round(A, 3)),
             paste(target, "/ A ratio:", round(ratio, 3)),
             paste(partner, "/ A ratio:", round(partner_ratio, 3)),
             paste("Dosage interpretation:", call),
             "Compare with BioSample sex and discuss limitations.")
writeLines(summary, "results/sex_inference.txt")
cat(paste(summary, collapse = "\n"), "\n")

### Coverage Table

| Accession Number | Chromosome Group | Mean Depth | Relative Depth |
|---|---|---:|---:|
| NC_004353.4 | Autosome (A) | 32.3133 | 1.0007 |
| NC_004354.4 | X Chromosome | 25.7999 | 0.7990 |
| NT_033777.3 | Autosome (A) | 32.2893 | 1.0000 |
| NT_033778.4 | Autosome (A) | 31.5393 | 0.9768 |
| NT_033779.5 | Autosome (A) | 33.6238 | 1.0413 |
| NT_037436.4 | Autosome (A) | 31.3983 | 0.9724 |
| NC_024512.1 | Y Chromosome | 1.4915 | 0.0462 |

## Questions for the Report

**1. What species did you choose, and what evidence supports a diploid XX/XY or ZZ/ZW system?**

The selected species is **Drosophila melanogaster* (fruit fly). It follows an XX/XY sex-chromosome system, where females generally have XX chromosomes and males have XY chromosomes. The X and Y chromosome sequences were included in the analysis to examine their sequencing depth and determine the possible chromosomal sex of the WGS sample.

**2. What are your assembly accession, WGS run, and both BioSample accessions? Are the donors the same?**

The reference genome assembly used was GCF_000001215.4, while the whole-genome sequencing run was SRR6702604. Their corresponding BioSample accessions are SAMN02803731 and SAMN08511563, respectively. These BioSamples have different accession numbers and describe different biological samples. Therefore, it cannot be confirmed that the reference assembly and WGS reads came from the same donor.

**3. Which sequences were used as autosomes and as X or Z? How did you verify their identities?**

| Chromosome Group | Accession Number | Purpose |
|---|---|---|
| Autosomes (A) | NC_004353.4 | Used as part of the autosomal reference for sequencing-depth comparison. |
| Autosomes (A) | NT_033777.3 | Used as part of the autosomal reference for sequencing-depth comparison. |
| Autosomes (A) | NT_033778.4 | Used as part of the autosomal reference for sequencing-depth comparison. |
| Autosomes (A) | NT_033779.5 | Used as part of the autosomal reference for sequencing-depth comparison. |
| Autosomes (A) | NT_037436.4 | Used as part of the autosomal reference for sequencing-depth comparison. |
| X Chromosome | NC_004354.4 | Used to calculate the X:A ratio and assess X-chromosome coverage. |
| Y Chromosome | NC_024512.1 | Used to calculate the Y:A ratio and assess Y-chromosome coverage. |

The chromosome sequences were classified using their accession numbers and chromosome annotations in the NCBI reference genome assembly report. Five autosomal reference sequences (2L, 2R, 3L, 3R, and chromosome 4) were grouped as A. The X chromosome (NC_004354.4) was classified as X, while the Y chromosome (NC_024512.1) was classified as Y.

**4. What are the median autosomal depth, X:A or Z:A ratio, and inferred sex-chromosome complement?**

| Measurement | Result | Interpretation |
|---|---:|---|
| Median Autosomal Depth | 32.289 | Used as the baseline for comparing chromosome coverage. |
| X:A Ratio | 0.799 | Falls between the expected XX (1.0) and XY (0.5) values. |
| Y:A Ratio | 0.046 | Indicates very low Y-chromosome sequencing depth. |
| Dosage Interpretation | Ambiguous | The result does not confidently establish an XX or XY complement. |

Under the ordinary diploid XX/XY system, an X:A ratio near 1.0 generally indicates XX, while a ratio near 0.5 generally indicates XY. Since the observed X:A ratio of 0.799 falls between these expected values, the result is considered ambiguous. The very low Y:A ratio indicates limited Y-chromosome read representation, but it does not independently confirm the absence of the Y chromosome. Therefore, the available sequencing-depth evidence is insufficient to confidently assign an XX or XY complement. The ambiguous result should also be considered in relation to the pooled nature of the WGS sample, which contained DNA from both male and female flies.

**5. Does your result agree with the recorded sex of the WGS donor, if reported? Explain any discrepancy or uncertainty.**

The NCBI BioSample record for the WGS sample (SAMN08511563) identifies the biological material as pooled male and female **Drosophila melanogaster*. Therefore, the sequencing data do not represent a single individual with one identifiable XX or XY chromosome complement. The observed X:A ratio of 0.799 may reflect the mixture of male and female DNA, along with possible technical factors affecting sequencing coverage and mapping. The low Y:A ratio does not prove that the Y chromosome is absent. Thus, the results cannot be used to confidently determine the chromosomal sex of an individual fly.

**6. Why is a chromosome-level assembly alone insufficient for a reliable read-depth sex inference?**

A chromosome-level assembly alone is insufficient because it only provides the reference genome sequences and does not show the actual sequencing depth of an individual sample. To determine the possible chromosomal sex, WGS reads must be mapped to the reference genome so that the coverage of the X and Y chromosomes can be compared with the autosomes. Therefore, both the genome assembly and WGS sequencing data are needed to calculate the ratios and make a more reliable interpretation.

## Limitations

The analysis has several limitations that may affect the results. First, the reference genome assembly (GCF_000001215.4) and WGS run (SRR6702604) have different BioSample accessions, so they cannot be confirmed as originating from the same biological sample. Another important limitation is that the WGS sample consisted of pooled male and female flies rather than one individual. Consequently, the sequencing-depth measurements represent a mixture of DNA from both sexes, making it inappropriate to infer the chromosomal complement of a single individual from the observed X:A ratio. The X:A ratio of 0.799 falls between the expected values for XX and XY individuals, resulting in an ambiguous interpretation. The Y:A ratio of 0.046 also indicates very low Y-chromosome read representation. This may be influenced by repetitive DNA sequences, mapping difficulties, or other sequencing limitations. Therefore, the low Y:A ratio does not automatically indicate the absence of the Y chromosome. Overall, the results should be interpreted with caution because of the pooled sample, differences in BioSample accessions, and possible sequencing and mapping limitations.

## Interpretation of the WGS Chromosome-Depth Analysis

The WGS analysis showed an X:A ratio of 0.799, a Y:A ratio of 0.046, and a median autosomal sequencing depth of 32.289.

In a normal diploid system, the expected X:A ratio is around 1.0 for XX (female) and 0.5 for XY (male). The observed X:A ratio of 0.799 falls between these values. Therefore, the result is considered ambiguous, and the chromosomal complement of the WGS donor cannot be confidently identified as either XX or XY.

The Y:A ratio of 0.046 shows that the sequencing depth of the Y chromosome is very low compared with the autosomes. However, this does not necessarily mean that the Y chromosome is absent. Several factors may have affected the results, such as low sequencing coverage, mapping bias, repetitive DNA sequences, pseudoautosomal regions, and missing or incomplete Y-chromosome sequences in the reference genome. These factors may have influenced the accuracy of chromosome-depth measurements.

**Comparison with BioSample Information:** The reported sex of the WGS BioSample was not available in the results. Therefore, the inferred chromosomal complement cannot yet be compared with the donor's reported sex.

The results show an intermediate X:A ratio and very low Y:A ratio. However, these values are not enough to determine whether the WGS donor has an XX or XY chromosomal complement. Additional BioSample information and consideration of possible sequencing and mapping limitations are needed to make a more reliable interpretation.

## References

1. Chang, C. H., & Larracuente, A. M. (2019). Heterochromatin-enriched assemblies reveal the sequence and organization of the *Drosophila melanogaster* Y chromosome. *Genetics, 211*(1), 333–348. https://doi.org/10.1534/genetics.118.301765

2. Cabrera, A. A., et al. (2022). How low can you go? Introducing SeXY: Sex identification from low-quantity sequencing data despite lacking assembled sex chromosomes. *Ecology and Evolution*. https://doi.org/10.1002/ece3.9185

3. National Center for Biotechnology Information. (2014). *Drosophila melanogaster: Release 6 plus ISO1 MT (GCF_000001215.4)* [Genome assembly]. NCBI Assembly. https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001215.4/

4. National Center for Biotechnology Information. (2018). *Illumina sequencing data for ISO-1 (SRR6702604)* [Sequence read archive]. NCBI Sequence Read Archive. https://www.ncbi.nlm.nih.gov/sra/SRR6702604

5. National Center for Biotechnology Information. (n.d.). Drosophila melanogaster (SAMN08511563) [BioSample]. NCBI BioSample. https://www.ncbi.nlm.nih.gov/biosample/SAMN08511563

7. National Center for Biotechnology Information. (n.d.). *Drosophila melanogaster* (SAMN08511563) [BioSample]. NCBI BioSample. https://www.ncbi.nlm.nih.gov/biosample/SAMN08511563

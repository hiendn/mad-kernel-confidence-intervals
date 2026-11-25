# Common-grid MAD manuscript revision

The revised joint-limit diagnostic uses sample sizes 500, 1,000, 2,000 and 5,000 with 5,000 Monte Carlo datasets per design for each of the normal and lognormal populations. Table 3 reports all eight designs. The reproduction archive also includes the principal and extension coverage results.

- [Manuscript PDF](MAD_Common_Grid_Manuscript.pdf)
- [Overleaf source](MAD_Common_Grid_Overleaf_Source.zip). Extract the ZIP and compile the main LaTeX file with pdfLaTeX. The numerical sections are editable and the accepted theory and proofs are included as PDF pages.
- Reproduction archive [part 1](MAD_Common_Grid_Reproduction.zip.001) and [part 2](MAD_Common_Grid_Reproduction.zip.002)

Download both archive parts and the checksum file to the same directory. Join the parts before checking and extracting the archive.

```sh
cat MAD_Common_Grid_Reproduction.zip.001 MAD_Common_Grid_Reproduction.zip.002 > MAD_Common_Grid_Reproduction.zip
unzip -t MAD_Common_Grid_Reproduction.zip
unzip MAD_Common_Grid_Reproduction.zip
```

[SHA256SUMS](SHA256SUMS) gives the checksums of the downloads and the reconstructed ZIP. The README inside the ZIP gives the replication commands and recorded design settings.

The original 15 November computational release is pinned to commit [0818d7624630](https://github.com/hiendn/mad-kernel-confidence-intervals/tree/0818d7624630ea7a9c765545127429921b275a63).

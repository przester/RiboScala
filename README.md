# RiboScala

A Ribo-Seq utility for calculating P-site offset metaplots and extracting P-site coverages.

### Usage
On Linux: `chmod +x gradlew && ./gradlew run --args="<arguments>"`

On Windows `.\gradlew.bat run --args="<arguments>"`

The utility takes three required file path arguments:
- genome feature annotation GTF file (`--genes <filename.gtf>`)
- BAM file containing reads for analysis (`--bam <filename.bam>`)
- target location for a CSV output file containing read coverages (`--out <filename.csv>`)

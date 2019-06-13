package ribo

import java.nio.file.{Files, Path}

import htsjdk.samtools.{QueryInterval, SAMSequenceDictionary, SamReaderFactory}

import scala.collection.JavaConverters._

class PSiteCoverage(refs: SAMSequenceDictionary, refCoverage: Array[Array[Array[Int]]]) {
  /**
    * Fetch coverage array, indexed by biocoordinates.
    * @param refIdx         reference sequence index (as provided in SAM header)
    * @param negativeStrand specifies what strand to fetch coverage for
    * @return
    */
  def apply(refIdx: Int, negativeStrand: Boolean): Array[Int] =
    refCoverage(refIdx)(if (!negativeStrand) 0 else 1)

  /**
    * Create a CSV file containing P-site coverage for each nucleotide in each reference sequence strand.
    * If the sequences have different lengths, the shorter ones are end-padded with zeros to yield a rectangular shape.
    * @param outFile     target file path
    * @param writeHeader insert a header row with seqName[+|-] column names
    */
  def saveCSV(outFile: Path, writeHeader: Boolean = true): Unit = {
    val maxRefLen = refs.getSequences.asScala.map(_.getSequenceLength).max
    val refNames = (0 until refs.size()).map(refs.getSequence(_).getSequenceName)
    val csv = Files.newBufferedWriter(outFile)
    if (writeHeader)
      csv.write((refNames.map(_ ++ "+") ++ refNames.map(_ ++ "-")).mkString(",") ++ "\n")
    (1 to maxRefLen).foreach { refPos =>
      val line =
        Vector(0, 1).flatMap { strandIdx =>
          (0 until refs.size()).map { refIdx =>
            if (refCoverage(refIdx)(strandIdx).isDefinedAt(refPos))
              refCoverage(refIdx)(strandIdx)(refPos)
            else
              0
          }}
          .mkString(",")
      csv.write(line ++ "\n")
    }
    csv.close()
  }
}

object PSiteCoverage {
  /**
    * Create a PSiteCoverage object using all mapped reads from samFile and offsets provided by metagene.
    * @param metagene
    * @param samFile
    * @return
    */
  def apply(metagene: MetageneProfile, samFile: Path): PSiteCoverage = {
    val refs = SamReaderFactory.makeDefault().open(samFile).getFileHeader.getSequenceDictionary
    val refCoverage = Array.ofDim[Array[Array[Int]]](refs.size())
    (0 until refs.size()).par
      .foreach { refIdx =>
        val samReader = SamReaderFactory.makeDefault().open(samFile)
        val refLength = refs.getSequence(refIdx).getSequenceLength
        refCoverage(refIdx) = Array.ofDim[Int](2, refLength + 1)
        samReader.query(Array(new QueryInterval(refIdx, 1, refLength)), false)
          .stream().iterator().asScala
          .filter(!_.getReadUnmappedFlag)
          .foreach { entry =>
            val isNegative = entry.getReadNegativeStrandFlag
            val pSitePos =
              if (!isNegative)
                entry.getAlignmentStart + metagene.offset(entry.getReadLength)
              else
                entry.getAlignmentEnd - metagene.offset(entry.getReadLength)
            refCoverage(refIdx)(if (!isNegative) 0 else 1)(pSitePos) += 1
          }
      }
    new PSiteCoverage(refs, refCoverage)
  }
}


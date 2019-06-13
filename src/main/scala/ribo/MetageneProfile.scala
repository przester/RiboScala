package ribo

import java.nio.file.Path

import htsjdk.samtools.{QueryInterval, SamReaderFactory}
import org.biojava.nbio.genome.parsers.gff.{FeatureI, FeatureList}

import scala.collection.JavaConverters._

/**
  * @param profile internal representation for easy merging - read count for each (read length, offset) pair
  */
class MetageneProfile(val profile: Map[(Int, Int), Int]) {
  /**
    * Merge two metagene profiles, summing their distributions.
    */
  def ++(that: MetageneProfile) =
    new MetageneProfile(this.profile ++ that.profile.map { case (k, v) => k -> (v + this.profile.getOrElse(k, 0)) })

  /**
    * A vector of all encountered read lengths (sorted in ascending order).
    */
  lazy val readLengths: Vector[Int] =
    profile.toVector.map(_._1._1).distinct.sorted

  /**
    * An (offset -> read count) map for each read length.
    */
  lazy val offsetPlot: Map[Int, Map[Int, Int]] =
    profile.toVector
      .groupBy(_._1._1)
      .mapValues(_.map { case ((_, off), count) => (off, count) }.toMap)

  /**
    * Total read count for each read length.
    */
  lazy val readCount: Map[Int, Int] =
    offsetPlot.mapValues(_.values.sum)

  /**
    * Most abundant offset for each read length.
    */
  lazy val offsetMap: Map[Int, Int] =
    offsetPlot.mapValues(_.maxBy(_._2)._1)

  /**
    * Same as offsetMap, but with considerably faster access.
    */
  lazy val offset: Vector[Int] =
    (0 to readLengths.max).map(offsetMap.withDefaultValue(0)(_)).toVector

  /**
    * Periodicity score for each read length.
    * <p>
    * The score is defined here as fraction of reads that are in-frame with the most abundant offset.
    */
  lazy val periodicity: Map[Int, Float] =
    offsetMap.map { case (readLen, bestOffset) =>
      val framePlot = offsetPlot(readLen).groupBy(_._1 % 3)
        .map { case (frame, frameReads) => frame -> frameReads.values.sum }
      readLen -> framePlot(bestOffset % 3).toFloat / readCount(readLen).toFloat
    }
}

object MetageneProfile {
  /**
    * Calculate a metagene profile containing read offset distribution upstream of every gene's start codon.
    * @param samFile SAM file containing the reads
    * @param genes   genome feature list
    * @return
    */
  def apply(samFile: Path, genes: FeatureList): MetageneProfile = {
    genes.asScala.par
      .filter(_.`type`() == "CDS")
      .filter(_.getAttribute("exon_number") == "1")
      .map(MetageneProfile.fromGene(samFile, _))
      .reduce((p1, p2) => p1 ++ p2)
  }

  def fromGene(samFile: Path, gene: FeatureI): MetageneProfile = {
    val samReader = SamReaderFactory.makeDefault().open(samFile)
    val isNegativeStrand = gene.location().isNegative
    val refIndex = samReader.getFileHeader.getSequenceDictionary.getSequenceIndex(gene.seqname())
    val refPos = if (!isNegativeStrand) gene.location().bioStart() else gene.location().bioEnd()
    val entries = samReader.query(Array(new QueryInterval(refIndex, refPos, refPos)), false)
    val profile =
      entries.stream().iterator().asScala
        .filter(_.getReadNegativeStrandFlag == isNegativeStrand)
        .map(e => (e.getReadLength,
          if (!isNegativeStrand) refPos - e.getAlignmentStart else e.getAlignmentEnd - refPos))
        .filter(_._2 > 0)
        .toVector
        .groupBy(identity)
        .mapValues(_.size)
    new MetageneProfile(profile)
  }
}

import java.io.File
import java.nio.file.{FileSystems, Path}

import scopt.OParser
import org.apache.logging.log4j.scala.Logging
import org.biojava.nbio.genome.parsers.gff.GFF3Reader
import ribo.{MetageneProfile, PSiteCoverage}

case class Config(geneFile: File = new File(""),
                  samFile: File = new File(""),
                  coverageFile: File = new File(""))

object Config {
  def parse(args: Array[String]): Option[Config] = {
    val builder = OParser.builder[Config]
    val parser = {
      import builder._
      OParser.sequence(
        programName("RiboScala"),
        opt[File]("genes")
          .required()
          .action((x, c) => c.copy(geneFile = x)),
        opt[File]("bam")
          .required()
          .action((x, c) => c.copy(samFile = x)),
        opt[File]("out")
          .required()
          .action((x, c) => c.copy(coverageFile = x))
      )
    }
    OParser.parse(parser, args, Config())
  }
}

object RiboScala extends App with Logging {
  Config.parse(args) match {
    case Some(config) =>
      val fs = FileSystems.getDefault

      println(config.geneFile)

      val genes = GFF3Reader.read(config.geneFile.getAbsolutePath)
      val samFile = fs.getPath(config.samFile.getAbsolutePath)
      val coverageFile = fs.getPath(config.coverageFile.getAbsolutePath)

      logger.info("Calculating metagene offsets...")
      val metagene = MetageneProfile(samFile, genes)

      logger.info(Vector("Read length", "Read count", "Read %\t", "Periodicity", "Offset").mkString("\t"))
      metagene.readLengths.foreach { readLen =>
        logger.info(
          Vector(
            s"$readLen",
            s"${metagene.readCount(readLen)}",
            s"${100*metagene.readCount(readLen)/metagene.readCount.values.sum}%",
            f"${metagene.periodicity(readLen)}%.2f",
            s"${metagene.offset(readLen)}"
          ).mkString("\t\t"))
      }

      logger.info("Calculating P-site coverage...")
      val coverage = PSiteCoverage(metagene, samFile)

      logger.info("Saving P-site coverage...")
      coverage.saveCSV(coverageFile)

      logger.info("Finished.")
    case _ => ()
  }
}

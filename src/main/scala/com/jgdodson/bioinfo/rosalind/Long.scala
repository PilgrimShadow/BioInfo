package com.jgdodson.bioinfo.rosalind

import com.jgdodson.bioinfo.utils.Utils.readFastaFile

object Long {

  def main(args: Array[String]): Unit = {

    val reads = readFastaFile(args(0)).map(_._2).toSet
    println(assemble(reads))

  }

  // TODO: cut this in half
  def assemble(reads: Set[String]): String = {

    def attemptPrepend(read: String, genome: String): Option[String] = {

      val len = read.length
      val startInd = math.ceil(len / 2.0).toInt - 1

      (startInd to 0 by -1).find { i =>
        read.substring(i) == genome.substring(0, len - i)
      } map { j =>
        read.substring(0, j) ++ genome
      }
    }

    def attemptAppend(genome: String, read: String): Option[String] = {

      val readLen = read.length
      val genLen = genome.length
      val startInd = math.floor(readLen / 2.0).toInt + 1

      (startInd to readLen).find { i =>
        genome.substring(genLen - i) == read.substring(0, i)
      } map { j =>
        genome.substring(0, genLen - j) ++ read
      }
    }

    def firstStage(reads: Set[String]): (String, Set[String]) = {

      def loop(genome: String, remaining: Set[String]): (String, Set[String]) = {

        val res = remaining.find { s =>
          attemptAppend(genome, s).isDefined
        }

        res match {
          case Some(s) => loop(attemptAppend(genome, s).get, remaining - s)
          case None => (genome, remaining)
        }
      }

      val longest = reads.maxBy(_.length)

      loop(longest, reads - longest)
    }

    def secondStage(genome: String, reads: Set[String]): String = {

      def loop(genome: String, remaining: Set[String]): String = {

        val res = remaining.find { s =>
          attemptPrepend(s, genome).isDefined
        }

        res match {
          case Some(s) => loop(attemptPrepend(s, genome).get, remaining - s)
          case None => genome
        }
      }

      loop(genome, reads)
    }

    val (partialGenome, remaining) = firstStage(reads)
    secondStage(partialGenome, remaining)
  }

}

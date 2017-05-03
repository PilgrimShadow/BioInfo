package com.jgdodson.rosalind

import utils.Utils

object Prob {

  def main(args: Array[String]): Unit = {
    val lines = Utils.readLines(args(0))
    val pattern = lines.head
    val gc = lines(1).split(" ").map(_.toDouble)
    val logProbs = gc.map(item => prob(item, pattern))
    println(logProbs.mkString(" "))
  }

  def prob(gcContent: Double, dna: String): Double = {

    assert(0 <= gcContent && gcContent <= 1.0)

    val gcProbLog = math.log10(gcContent / 2)
    val atProbLog = math.log10((1 - gcContent) / 2)

    dna.foldLeft(0.0) { (acc, next) =>
      acc + (if (next == 'G' || next == 'C') gcProbLog else atProbLog)
    }
  }
}

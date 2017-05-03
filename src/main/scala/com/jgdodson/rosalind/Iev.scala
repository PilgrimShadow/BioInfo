package com.jgdodson.rosalind

import utils.Utils

object Iev {

  def main(args: Array[String]): Unit = {
    val counts = Utils.readLines(args(0)).head.split(" ").map(Integer.parseInt)
    println(expectedOffspring(counts))
  }

  def expectedOffspring(counts: Seq[Int]): Double = {

    assert(counts.length == 6)

    val condProbs = Seq(1.0, 1.0, 1.0, .75, .5, 0.0)

    2 * (0 until 6).foldLeft(0.0) { (sum, i) => sum + condProbs(i) * counts(i) }
  }
}

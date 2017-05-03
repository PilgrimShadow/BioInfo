package com.jgdodson.rosalind

import utils.Utils

object Grph {

  def main(args: Array[String]): Unit = {
    val strings = Utils.readFastaFile(args(0))
    val overlaps = findOverlaps(3, strings)
    val formatted = overlaps.map(entry => entry._1 + " " + entry._2).mkString("\n")
    println(formatted)
  }

  def findOverlaps(n: Int, strings: Vector[(String, String)]): Seq[(String, String)] = {

    def overlap(k: Int, s1: String, s2: String): Boolean = {
      s1.substring(s1.length - 3, s1.length) == s2.substring(0, 3)
    }

    for {
      i <- 0 until strings.length
      j <- 0 until strings.length if i != j && overlap(3, strings(i)._2, strings(j)._2)
    } yield (strings(i)._1, strings(j)._1)
  }

}

package com.jgdodson.rosalind

import utils.Utils.readFastaFile

object Tran {

  def main(args: Array[String]): Unit = {

    val seqs = readFastaFile(args(0)).map(_._2)
    println(ttRatio(seqs(0), seqs(1)))
  }

  def ttRatio(seq1: String, seq2: String): Double = {

    assert(seq1.length == seq2.length)

    def isTransition(base1: Char, base2: Char): Boolean = (base1, base2) match {
      case ('A', 'G') => true
      case ('G', 'A') => true
      case ('C', 'T') => true
      case ('T', 'C') => true
      case _ => false
    }

    def isTransversion(base1: Char, base2: Char): Boolean = (base1, base2) match {
      case ('A', 'C') => true
      case ('C', 'A') => true
      case ('G', 'T') => true
      case ('T', 'G') => true
      case ('A', 'T') => true
      case ('T', 'A') => true
      case ('C', 'G') => true
      case ('G', 'C') => true
      case _ => false
    }

    val ttCount = (0 until seq1.length).foldLeft(0, 0) { (acc, i) =>
      if (isTransition(seq1(i), seq2(i))) (acc._1 + 1, acc._2)
      else if (isTransversion(seq1(i), seq2(i))) (acc._1, acc._2 + 1)
      else acc
    }

    ttCount._1.toDouble / ttCount._2
  }

}

package com.jgdodson.rosalind

import utils.Utils

object Cons {

  def main(args: Array[String]): Unit = {
    val pairs = Utils.readFastaFile(args(0))
    val strings = pairs.map(_._2)
    val profile = makeProfile(strings)
    println(makeConsensus(profile))
    println(formatProfile(profile))
  }

  def makeProfile(strings: Vector[String]): Vector[Vector[Int]] = {

    val len = strings.head.length

    val counts = Array.fill(4) {
      Array.fill[Int](len)(0)
    }

    for (string <- strings; i <- 0 until string.length) {
      if (string(i) == 'A') counts(0)(i) += 1
      else if (string(i) == 'C') counts(1)(i) += 1
      else if (string(i) == 'G') counts(2)(i) += 1
      else counts(3)(i) += 1
    }

    counts.map(_.toVector).toVector
  }

  def makeConsensus(profile: Vector[Vector[Int]]): String = {

    def decide(A: Int, C: Int, G: Int, T: Int): Char = {

      val max = Vector(A, C, G, T).max
      if (max == A) 'A'
      else if (max == C) 'C'
      else if (max == G) 'G'
      else 'T'
    }

    (for (i <- 0 until profile.head.length) yield {
      decide(profile(0)(i), profile(1)(i), profile(2)(i), profile(3)(i))
    }).mkString("")

  }

  def formatProfile(profile: Vector[Vector[Int]]): String = {
    val labels = List("A: ", "C: ", "G: ", "T: ").iterator
    profile.map(line => labels.next + line.mkString(" ")).mkString("\n")
  }


}

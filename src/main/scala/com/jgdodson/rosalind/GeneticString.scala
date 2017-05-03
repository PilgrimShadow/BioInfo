package com.jgdodson.rosalind

// It's weired how the types work out just right... (Curiously Recurring Template Pattern)
abstract class GeneticString[T <: GeneticString[T]] {

  val seq: String

  // The order of the characters is significant as it is used as a default
  // ordering when manipulating sequences.
  def alphabet: Seq[Char]

  def length: Int = seq.length

  def masses: Map[Char, Double]

  def mass: Double = seq.foldLeft(0.0)(_ + masses(_))

  def reverse: T

  def substring(start: Int, end: Int): T

  // TODO: convert this into an explicit motif finding function
  def failureArray: Vector[Int] = {

    (1 until seq.length).foldLeft((Vector[Int](0), Set[Int](0))) { (acc, next) =>
      val updated = acc._2.filter(i => seq(next) == seq(i)).map(i => i + 1) + 0
      (acc._1 :+ updated.max, updated)
    }._1
  }

  // Based on the Knuth-Morris-Pratt algorithm.
  def findMotif(other: T): Vector[Int] = {

    (0 until seq.length).foldLeft((Vector[Int](), Set[Int](0))) { (acc, next) =>
      val updated = acc._2.filter(i => seq(next) == other.seq(i)).map(_ + 1) + 0

      if (updated.max == other.length) (acc._1 :+ (next - updated.max + 1), updated - updated.max)
      else (acc._1, updated)
    }._1
  }

  def kmerComposition(k: Int): Seq[Int] = {

    val indices = Lexf.enumerateFixed(k, alphabet).zipWithIndex.toMap

    (0 to length - k).foldLeft(for (_ <- Seq.range(0, indices.size)) yield 0) { (acc, next) =>
      val index = indices(seq.substring(next, next + k))
      acc.updated(index, acc(index) + 1)
    }
  }

  // Allows for strings of different length
  def hammingDistance(other: T): Int = {
    val min = Math.min(length, other.length)
    val max = Math.max(length, other.length)

    (0 until min).count(i => seq(i) != other.seq(i)) + (max - min)
  }

  def findSplicedMotif(motif: T): Option[Vector[Int]] = {

    val res = (0 until seq.length).foldLeft(motif.seq, Vector[Int]()) { (acc, i) =>
      if (acc._1.isEmpty) acc
      else if (seq(i) == acc._1.head) (acc._1.tail, acc._2 :+ i)
      else acc
    }

    if (res._1.isEmpty) Some(res._2)
    else None
  }

  def transitionTable(): Vector[Vector[Int]] = {
    val inds = alphabet.zipWithIndex.toMap
    val res = Array.fill(alphabet.length, alphabet.length)(0)
    for (pair <- seq.sliding(2)) {
      res(inds(pair(1)))(inds(pair(0))) += 1
    }
    res.map(_.toVector).toVector
  }

  // TODO: Finish this
  def longestSharedSplicedMotif(other: T): T = {

    val minLen = math.min(seq.length, other.seq.length)

    (0 until minLen).foldLeft(((Set[Char](), Set[Char]()), "")) { (acc, next) =>
      ???
    }
    ???
  }
}

package com.jgdodson.rosalind


/**
  *
  * It's weird how the types work out just right... (Curiously Recurring Template Pattern)
  *
  * @tparam T
  */
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


  /**
    * Find all occurences of the given motif within this genetic string
    *
    * Based on the Knuth-Morris-Pratt algorithm.
    *
    * @param motif
    * @return
    */
  def findMotif(motif: T): Vector[Int] = {

    (0 until seq.length).foldLeft((Vector[Int](), Set[Int](0))) { (acc, next) =>
      val updated = acc._2.filter(i => seq(next) == motif.seq(i)).map(_ + 1) + 0

      if (updated.max == motif.length) (acc._1 :+ (next - updated.max + 1), updated - updated.max)
      else (acc._1, updated)
    }._1
  }

  /**
    * Alternate implementation for testing purposes
    *
    * @param motif
    * @return
    */
  def findMotif2(motif: T): Vector[Int] = {

    val candidates = collection.mutable.Set[Int](0)
    val matches = collection.mutable.Seq[Int]()

    for (i <- 0 until seq.length) {
      candidates.filter(j => seq(i) == motif.seq(j)).map(_ + 1) + 0

      if (candidates.max == motif.length) {
        matches :+ (i - motif.length + 1)
      }
    }

    // Return the vector of matching indices
    matches.toVector
  }

  // TODO: Write a Regex-based motif finder. Compare speeds.

  /**
    *
    * @param k
    * @return
    */
  def kmerComposition(k: Int): Seq[Int] = {

    val indices = Lexf.enumerateFixed(k, alphabet).zipWithIndex.toMap

    (0 to length - k).foldLeft(for (_ <- Seq.range(0, indices.size)) yield 0) { (acc, next) =>
      val index = indices(seq.substring(next, next + k))
      acc.updated(index, acc(index) + 1)
    }
  }


  /**
    * Compute the Hamming distance between this genetic string and another
    *
    * Allows for strings of different length
    *
    * @param other The other genetic string
    * @return
    */
  def hammingDistance(other: T): Int = {
    val min = Math.min(length, other.length)
    val max = Math.max(length, other.length)

    (0 until min).count(i => seq(i) != other.seq(i)) + (max - min)
  }


  /**
    *
    * @param motif
    * @return
    */
  def findSplicedMotif(motif: T): Option[Vector[Int]] = {

    val res = (0 until seq.length).foldLeft(motif.seq, Vector[Int]()) { (acc, i) =>
      if (acc._1.isEmpty) acc
      else if (seq(i) == acc._1.head) (acc._1.tail, acc._2 :+ i)
      else acc
    }

    if (res._1.isEmpty) Some(res._2)
    else None
  }

  /**
    *
    * @return
    */
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

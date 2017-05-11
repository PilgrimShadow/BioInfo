package com.jgdodson.bioinfo

// TODO: Remove this dependency
import com.jgdodson.bioinfo.rosalind.Lexf

import scala.collection.mutable

/**
  *
  * It's weird how the types work out just right... (Curiously Recurring Template Pattern)
  *
  * @tparam T
  */
abstract class GeneticString[T <: GeneticString[T]] {

  val seq: String


  /**
    *
    * @param i The index of the character to return
    * @return
    */
  def apply(i: Int): Char = seq(i)


  /**
    * The number of base pairs in this genetic string
    *
    * @return
    */
  def length: Int = seq.length


  // The order of the characters is significant as it is used as a default
  // ordering when manipulating sequences.
  protected def alphabet: Seq[Char]

  protected def masses: Map[Char, Double]

  def mass: Double = seq.foldLeft(0.0)(_ + masses(_))

  def reverse: T

  def substring(start: Int, end: Int): T


  /**
    * Indicate whether this GeneticString contains a spliced version of the given motif
    *
    * @param motif The motif to search for
    * @return
    */
  def containsSplicedMotif(motif: T): Boolean = {

    // Motif is not longer than this sequence
    // Should we return false in this case or throw?
    assert(motif.length <= this.length)

    var rest = 0

    for (char <- motif.seq) {
      rest = this.seq.indexOf(char, rest) + 1

      if (rest == 0) {
        return false
      }
    }

    true
  }


  // TODO: convert this into an explicit motif finding function
  def failureArray: Vector[Int] = {

    (1 until seq.length).foldLeft((Vector[Int](0), Set[Int](0))) { (acc, next) =>
      val updated = acc._2.filter(i => seq(next) == seq(i)).map(i => i + 1) + 0
      (acc._1 :+ updated.max, updated)
    }._1
  }


  /**
    * Find all occurrences of the given motif within this genetic string
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

    val candidates = mutable.Set[Int](0)
    val matches = mutable.Seq[Int]()

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
    * Finds the starting index of all spliced occurrences of the given motif
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
    * Find one of the longest shared spliced motifs
    *
    * TODO: A function that returns ALL of the longest SSMs
    *
    * @param other The sequence to search against this sequence
    * @return
    */
  def longestSharedSplicedMotif(other: T): String = {

    // Determine which sequence is shorter/longer
    val minSeq = if (this.length == other.length) this else List(this, other).minBy(_.length)
    val maxSeq = if (this.length == other.length) other else List(this, other).maxBy(_.length)

    // Store partial solutions as we go
    val t = mutable.Map[Int, String](-1 -> "")

    for (char <- minSeq.seq) {

      // Extend all of our partial SSMs with the next char.
      // Ensure we don't modify the same t that we iterate over
      // TODO: Likely quicker to create an update map and then modify t after this for loop (well, maybe...)
      for ((k, v) <- t.toMap) {

        // The next match for this SSM
        val x = maxSeq.seq.indexOf(char, k + 1)

        // If a match was found
        if (x != -1) {

          // Update the SSM ending at index x if necessary
          if (!t.contains(x) || (v.length >= t(x).length)) {
            // TODO: This will modify the update map instead
            t(x) = v + char
          }

        }

      }

      // TODO: Update t with the update map here

    }

    // Return one of the longest SSMs
    t.maxBy(_._2.length)._2
  }


  /**
    *
    * @param other
    * @return
    */
  def shortestCommonSuperSequence(other: T): String = {

    val lssm = this.longestSharedSplicedMotif(other)

    var i = 0
    var j = 0

    val res = new mutable.StringBuilder()

    for (c <- lssm) {

      val s = this.seq.substring(i, this.seq.indexOf(c, i))
      val t = other.seq.substring(j, other.seq.indexOf(c, j) + 1)

      // Must add t after s
      res ++= s
      res ++= t

      i += s.length + 1
      j += t.length
    }

    // Add any leftovers
    res ++= this.seq.substring(i)
    res ++= other.seq.substring(j)

    res.toString
  }


  def longestSharedSplicedMotif3(other: T): String = {

    val minSeq = if (this.length == other.length) this else List(this, other).minBy(_.length)
    val maxSeq = if (this.length == other.length) other else List(this, other).maxBy(_.length)

    val t = Vector.fill(minSeq.length)(mutable.Set[Set[Int]]())

    val r = Range(0, minSeq.length).dropWhile(i => maxSeq.seq.indexOf(minSeq(i)) == -1)

    // The first shared substring
    t(r.min).add(Set(maxSeq.seq.indexOf(minSeq(r.min))))

    for (i <- r.drop(1)) {

      val updated = t.slice(0, i).map(_.map(ssm => {
        val x = maxSeq.seq.indexOf(minSeq(i), ssm.max + 1)

        if (x == -1) ssm else ssm + x
      }))

      val maxSize = updated.map(_.map(_.size).max).max

      val n: Vector[Set[Int]] = updated.flatMap(_.filter(ssm => ssm.size == maxSize))

      for (ssm <- n) {
        t(i).add(ssm)
      }
    }

    t.last.head.toVector.sorted.map(i => maxSeq(i)).mkString
  }


  /**
    * Compute the edit distance between this genetic string and another
    *
    * @param other The other genetic string
    * @return
    */
  def editDistance(other: T): Int = {

    val lssm = this.longestSharedSplicedMotif(other)

    var i = 0
    var j = 0

    var res = 0

    for (c <- lssm) {

      val s = this.seq.substring(i, this.seq.indexOf(c, i))
      val t = other.seq.substring(j, other.seq.indexOf(c, j))

      res += math.max(s.length, t.length)

      i += s.length + 1
      j += t.length + 1
    }

    // Add any leftovers
    res += math.max(this.seq.substring(i).length, other.seq.substring(j).length)

    res
  }


  /**
    * Find one of the longest shared spliced motifs
    *
    * TODO: This should return a T
    *
    * @param other
    * @return
    */
  def longestSharedSplicedMotif2(other: T): String = {

    val minSeq = if (this.length == other.length) this else List(this, other).minBy(_.length)
    val maxSeq = if (this.length == other.length) other else List(this, other).maxBy(_.length)

    val t = Vector.fill(minSeq.length)(mutable.Set[Set[Int]]())

    for ((char, i) <- minSeq.seq.zipWithIndex) {

      val x = maxSeq.seq.indexOf(char)

      if (x != -1) {
        t(i).add(Set(x))
      }

      // Update all previous SSM candidates
      for (j <- 0 until i) {
        for (ssm <- t(j)) {

          val x = maxSeq.seq.indexOf(char, ssm.max + 1)
          if (x != -1) {
            t(j).add(ssm + x)
          }
        }
      }
    }

    val q = t.map(_.maxBy(_.size)).maxBy(_.size)
    q.toVector.sorted.map(i => maxSeq(i)).mkString
  }
}
package com.jgdodson.bioinfo

case class DNAString(seq: String) extends GeneticString[DNAString] {

  // TODO: What about Set(seq) >: alphabet. Is one faster?

  // Initialization error checking
  if (seq.exists(ch => !alphabet.contains(ch))) {
    throw new Error("DNA string contains invalid character.")
  } else if (seq.length == 0) {
    throw new Error("DNAString cannot be empty")
  }

  protected def alphabet: Seq[Char] = DNAString.alphabet

  protected def masses = DNAString.masses


  /**
    * Return an RNAString where all T's in this DNAString are replaced with U's
    *
    * @return
    */
  def toRNAString: RNAString = RNAString(seq.replace('T', 'U'))


  def substring(start: Int, end: Int): DNAString = {
    DNAString(seq.substring(start, end))
  }

  /**
    * Return the complement of this DNAString
    *
    * @return
    */
  def complement: DNAString = DNAString(seq.map(DNAString.complements))


  /**
    * The fraction of C-G base pairs
    *
    * @return
    */
  lazy val GCcontent: Double = seq.foldLeft(0D) { (acc, next) =>
    if (Set('C', 'G').contains(next)) acc + 1
    else acc
  } / length


  /**
    * The fraction of A-T base pairs
    *
    * @return
    */
  lazy val ATcontent: Double = 1 - GCcontent


  /**
    * Return the reverse of this DNAString
    *
    * @return
    */
  def reverse: DNAString = DNAString(seq.reverse)


  /**
    * Return the reverse complement of this DNAString
    *
    * @return
    */
  def reverseComplement: DNAString = DNAString(seq.reverseMap(DNAString.complements))


  /**
    * Compute the number of transitions between this DNAString and another
    *
    * @param other
    * @return
    */
  def transitionCount(other: DNAString): Int = {

    assert(this.length == other.length)

    def isTransition(p: (Char, Char)): Boolean = p match {
      case ('A', 'G') => true
      case ('G', 'A') => true
      case ('C', 'T') => true
      case ('T', 'C') => true
      case _ => false
    }

    this.seq.zip(other.seq).count(isTransition)
  }


  /**
    * Compute the number of transversions between this DNAString and another
    *
    * @param other
    * @return
    */
  def transversionCount(other: DNAString): Int = {

    assert(this.length == other.length)

    def isTransversion(p: (Char, Char)): Boolean = p match {
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

    this.seq.zip(other.seq).count(isTransversion)
  }


  /**
    * Compute the transition/transversion ratio between this DNAString and another
    *
    * @param other
    * @return
    */
  def ttRatio(other: DNAString): Double = {
    transitionCount(other).toDouble / transversionCount(other).toDouble
  }


  /**
    *
    * @param introns
    * @return
    */
  def removeIntrons(introns: Seq[DNAString]): DNAString = {
    DNAString(seq.split(introns.mkString("|")).mkString(""))
  }

  // Review this function
  def restrictionSites: Set[(Int, Int)] = {

    (0 until seq.length).foldLeft(Set[(Int, Int)]()) { (acc, i) =>

      val res = (4 to 12 by 2).filter { j =>

        val end = i + j

        if (end > seq.length) false
        else {
          val test = seq.substring(i, end)
          test.reverseMap(DNAString.complements) == test
        }
      }

      // Convert to one-based indexing
      acc ++ res.map(entry => (i + 1, entry)).toSet
    }

  }

}

/**
  * Companion object for the DNAString class
  */
object DNAString {

  val alphabet: Seq[Char] = Seq('A', 'C', 'G', 'T')

  // Might need to change these
  val masses = Map(
    'A' -> 135.054489,
    'C' -> 111.043259,
    'G' -> 151.049408,
    'T' -> 126.042931)

  val complements = Map(
    'A' -> 'T',
    'C' -> 'G',
    'G' -> 'C',
    'T' -> 'A'
  )


  /**
    * Compute the profile matrix for a group of DNAStrings
    *
    * @param group
    * @return
    */
  def profileMatrix(group: Seq[DNAString]): Vector[Vector[Int]] = {

    assert(Set(group.map(_.length)).size == 1)

    val counts = Vector.fill(4) {
      Array.fill[Int](group.head.length)(0)
    }

    for (dna <- group; i <- 0 until dna.length) {
      if (dna(i) == 'A') counts(0)(i) += 1
      else if (dna(i) == 'C') counts(1)(i) += 1
      else if (dna(i) == 'G') counts(2)(i) += 1
      else counts(3)(i) += 1
    }

    counts.map(_.toVector)
  }


  /**
    *
    * @param group
    * @return
    */
  def profileMap(group: Seq[DNAString]): Map[Char, Vector[Int]] = {

    assert(Set(group.map(_.length)).size == 1)

    val counts = List('A', 'C', 'G', 'T').map(c => (c, Array.fill[Int](group.head.length)(0))).toMap

    for (dna <- group; i <- 0 until dna.length) {
      counts(dna(i))(i) += 1
    }

    counts.mapValues(_.toVector)
  }

  def consensus(group: Seq[DNAString]): DNAString = {
    ???
  }


  // TODO: Refactor this
  def consensusFromProfile(profile: Vector[Vector[Int]]): DNAString = {

    def decide(A: Int, C: Int, G: Int, T: Int): Char = {

      val max = Vector(A, C, G, T).max
      if (max == A) 'A'
      else if (max == C) 'C'
      else if (max == G) 'G'
      else 'T'
    }

    DNAString(profile.head.indices.map(i => decide(profile(0)(i), profile(1)(i), profile(2)(i), profile(3)(i))).mkString(""))
  }


  /**
    * Compute a matrix of Hamming distances for the group of DNAStrings
    *
    * TODO: VectorBuilder might be useful here
    *
    * @param group
    * @return
    */
  def distanceMatrix(group: Seq[DNAString]): Vector[Vector[Int]] = {

    // All sequences have same length
    assert(Set(group.map(_.length)).size == 1)

    val distMatrix = Vector.fill(group.length) {
      Array.fill[Int](group.length)(0)
    }

    for (i <- group.indices; j <- i + 1 until group.length) {

      val dist = group(i).hammingDistance(group(j))

      distMatrix(i)(j) = dist
      distMatrix(j)(i) = dist
    }

    distMatrix.map(_.toVector)
  }

  def normalizedDistanceMatrix(group: Seq[DNAString]): Vector[Vector[Double]] = {

    distanceMatrix(group).map(_.map(dist => dist.toDouble / group.head.length))
  }

}

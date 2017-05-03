package com.jgdodson.rosalind

case class DNAString(seq: String) extends GeneticString[DNAString] {

  // Error checking during initialization.
  if (seq.exists(ch => !alphabet.contains(ch))) {
    throw new Error("DNA string contains invalid character.")
  } else if (seq.length == 0) {
    throw new Error("DNAString cannot be empty")
  }

  def alphabet: Seq[Char] = DNAString.alphabet

  def masses = DNAString.masses

  def toRNAString: RNAString = RNAString(seq.replace('T', 'U'))

  def substring(start: Int, end: Int): DNAString = {
    DNAString(seq.substring(start, end))
  }

  def complement: DNAString = DNAString(seq.map(DNAString.complements))

  def gcContent: Double =
    seq.foldLeft(0.0) { (acc, next) =>
      if (Set('C', 'G').contains(next)) acc + 1
      else acc
    } / length

  def atContent: Double = 1 - gcContent

  def reverse: DNAString = DNAString(seq.reverse)

  def reverseComplement: DNAString = DNAString(seq.reverseMap(DNAString.complements))


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

}

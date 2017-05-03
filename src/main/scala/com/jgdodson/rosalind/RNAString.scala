package com.jgdodson.rosalind

case class RNAString(seq: String) extends GeneticString[RNAString] {

  // Error checking during initialization.
  if (seq.exists(ch => !alphabet.contains(ch))) {
    throw new Error("RNA string contains invalid character.")
  } else if (seq.length == 0) {
    throw new Error("RNAString cannot be empty")
  }

  def alphabet: Seq[Char] = RNAString.alphabet

  def masses = RNAString.masses

  def toDNAString: DNAString = DNAString(seq.replace('U', 'T'))

  // TODO: improve this
  def toProteinString: ProteinString = {

    val raw = seq.grouped(3).foldLeft("") { (acc, next) =>
      if (next.length == 3) acc + RNAString.rnaCodonTable(next)
      else acc
    }

    ProteinString(raw)
  }

  def reverse: RNAString = RNAString(seq.reverse)

  def substring(start: Int, end: Int): RNAString = {
    RNAString(seq.substring(start, end))
  }

  def openFrames: Set[(Int, Int)] = {

    def isStartCodon(codon: String): Boolean = {
      RNAString.rnaCodonTable(codon) == 'M'
    }

    def isStopCodon(codon: String): Boolean = {
      RNAString.rnaCodonTable(codon) == 'X'
    }

    // TODO: This is WAY too complicated. Clean it up, son.
    def onepass(seq: String) = {
      (0 to seq.length - 3).foldLeft((for (_ <- 0 to 2) yield Set[Int](), Set[(Int, Int)]())) { (acc, i) =>
        val binNum = i % 3
        val codon = seq.substring(i, i + 3)
        if (isStartCodon(codon)) {
          (acc._1.updated(binNum, acc._1(binNum) + i), acc._2)
        } else if (isStopCodon(codon)) {
          (acc._1.updated(binNum, Set[Int]()), acc._2 ++ acc._1(binNum).map((_, i)))
        } else acc
      }._2
    }

    onepass(seq) ++ onepass(seq.reverse)
  }


}

object RNAString {

  val alphabet = Seq('A', 'C', 'G', 'U')

  val masses = Map(
    'A' -> 135.054489,
    'C' -> 111.043259,
    'G' -> 151.049408,
    'U' -> 112.027275
  )

  // X is used for stop codons
  val rnaCodonTable: Map[String, Char] =
    Map(
      "UUU" -> 'F', "CUU" -> 'L', "AUU" -> 'I', "GUU" -> 'V',
      "UUC" -> 'F', "CUC" -> 'L', "AUC" -> 'I', "GUC" -> 'V',
      "UUA" -> 'L', "CUA" -> 'L', "AUA" -> 'I', "GUA" -> 'V',
      "UUG" -> 'L', "CUG" -> 'L', "AUG" -> 'M', "GUG" -> 'V',
      "UCU" -> 'S', "CCU" -> 'P', "ACU" -> 'T', "GCU" -> 'A',
      "UCC" -> 'S', "CCC" -> 'P', "ACC" -> 'T', "GCC" -> 'A',
      "UCA" -> 'S', "CCA" -> 'P', "ACA" -> 'T', "GCA" -> 'A',
      "UCG" -> 'S', "CCG" -> 'P', "ACG" -> 'T', "GCG" -> 'A',
      "UAU" -> 'Y', "CAU" -> 'H', "AAU" -> 'N', "GAU" -> 'D',
      "UAC" -> 'Y', "CAC" -> 'H', "AAC" -> 'N', "GAC" -> 'D',
      "UAA" -> 'X', "CAA" -> 'Q', "AAA" -> 'K', "GAA" -> 'E',
      "UAG" -> 'X', "CAG" -> 'Q', "AAG" -> 'K', "GAG" -> 'E',
      "UGU" -> 'C', "CGU" -> 'R', "AGU" -> 'S', "GGU" -> 'G',
      "UGC" -> 'C', "CGC" -> 'R', "AGC" -> 'S', "GGC" -> 'G',
      "UGA" -> 'X', "CGA" -> 'R', "AGA" -> 'R', "GGA" -> 'G',
      "UGG" -> 'W', "CGG" -> 'R', "AGG" -> 'R', "GGG" -> 'G'
    )

}
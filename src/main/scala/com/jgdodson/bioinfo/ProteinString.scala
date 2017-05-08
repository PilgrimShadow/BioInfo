package com.jgdodson.bioinfo

import utils.Utils.nearestMatch

case class ProteinString(seq: String) extends GeneticString[ProteinString] {

  // Error checking during initialization.
  if (seq.exists(ch => !alphabet.contains(ch))) {
    throw new Error("Peptide string contains invalid character.")
  } else if (seq.length == 0) {
    throw new Error("ProteinString cannot be empty")
  }

  protected def alphabet: Seq[Char] = ProteinString.alphabet

  protected def masses: Map[Char, Double] = ProteinString.masses

  def substring(start: Int, end: Int): ProteinString = {
    ProteinString(seq.substring(start, end))
  }

  def reverse: ProteinString = ProteinString(seq.reverse)

}

object ProteinString {

  def nearestAmino(searchMass: Double): Char = {

    masses.mapValues(actualMass => math.abs(actualMass - searchMass)).minBy(_._2)._1
  }


  def fromSpectrum(prefixMasses: Seq[Double]): ProteinString = {

    ProteinString(prefixMasses.indices.drop(1).map(i => nearestAmino(prefixMasses(i) - prefixMasses(i - 1))).mkString)
  }


  def fromFullSpectrum(ionMasses: Seq[Double], peptideMass: Double): ProteinString = {

    assert(ionMasses.length % 2 == 0)

    ???
  }


  // The character 'X' is used for the stop codon.
  val alphabet = Seq('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
    'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y') :+ 'X'

  // These are monoisotopic residue masses.
  val masses = Map(
    'A' -> 71.03711,
    'C' -> 103.00919,
    'D' -> 115.02694,
    'E' -> 129.04259,
    'F' -> 147.06841,
    'G' -> 57.02146,
    'H' -> 137.05891,
    'I' -> 113.08406,
    'K' -> 128.09496,
    'L' -> 113.08406,
    'M' -> 131.04049,
    'N' -> 114.04293,
    'P' -> 97.05276,
    'Q' -> 128.05858,
    'R' -> 156.10111,
    'S' -> 87.03203,
    'T' -> 101.04768,
    'V' -> 99.06841,
    'W' -> 186.07931,
    'Y' -> 163.06333
  )

}

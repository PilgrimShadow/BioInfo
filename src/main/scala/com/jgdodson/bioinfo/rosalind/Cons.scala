package com.jgdodson.bioinfo.rosalind

import com.jgdodson.bioinfo.DNAString
import com.jgdodson.bioinfo.utils.Utils

object Cons {

  def main(args: Array[String]): Unit = {
    val pairs = Utils.readFastaFile(args(0))

    val strings = pairs.map(_._2).map(DNAString(_))
    val profile = DNAString.profileMatrix(strings)

    println(DNAString.consensusFromProfile(profile))
    println(formatProfile(profile))
  }

  def formatProfile(profile: Vector[Vector[Int]]): String = {
    val labels = List("A: ", "C: ", "G: ", "T: ").iterator
    profile.map(line => labels.next + line.mkString(" ")).mkString("\n")
  }


}

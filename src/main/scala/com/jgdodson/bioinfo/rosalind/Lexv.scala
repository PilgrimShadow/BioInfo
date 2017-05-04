package com.jgdodson.bioinfo.rosalind

object Lexv {

  def main(args: Array[String]): Unit = {

    val lines = io.Source.fromFile(args(0)).mkString.split("\n")
    val alpha = lines(0).split(" ").map(a => a(0))
    val order = Lexf.makeOrdering(alpha)
    val len = Integer.parseInt(lines(1).trim())

    // Loosen up types to prevent these conversions
    val words = enumerateVar(len, alpha)
    val sorted = words.toVector.sorted(order)
    println(sorted.mkString("\n"))
  }

  def enumerateVar[T](maxLen: Int, alphabet: Seq[T]): Set[String] = {
    (1 to maxLen).foldLeft(Set[String]())((acc, next) => acc ++ Lexf.enumerateFixed(next, alphabet))
  }


}

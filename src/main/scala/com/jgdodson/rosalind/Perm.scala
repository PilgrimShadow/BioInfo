package com.jgdodson.rosalind

object Perm {

  def main(args: Array[String]): Unit = {

    val n = Integer.parseInt(args(0))

    val startSet = (1 to n).toSet

    val ans = formatResult(permutations(startSet))

    println(ans)
  }

  /**
   * Returns all permutation strings for the given set of Ints.
   *
   * @param ss The set of Ints to be permuted.
   * @return
   */
  def permutations[T](ss: Set[T]): Set[String] = {

    // Return a set with the given integer prepended to all strings.
    def prependAll(n: T, ss: Set[String]): Set[String] = {
      ss.map(item => n.toString + item)
    }

    if (ss.size > 1) {
      ss.foldLeft(Set[String]())((acc, next) => acc ++ prependAll(next, permutations(ss - next)))
    } else ss.map(_.toString)
  }

  def permutations2[T](ss: Set[T]): Set[Vector[T]] = {

    def prependAll(n: T, ss: Set[Vector[T]]): Set[Vector[T]] = {
      ss.map(item => n +: item)
    }

    if (ss.size > 1) {
      ss.foldLeft(Set[Vector[T]]()) { (acc, next) =>
        acc ++ prependAll(next, permutations2(ss - next))
      }
    } else ss.map(item => Vector(item))

  }

  /**
   * Format output for Rosalind problem: Enumerating Gene Orders.
   * @param ss
   * @return
   */
  def formatResult(ss: Set[String]): String = {

    val withSpaces = ss.map(str => str.mkString(" "))

    ss.size.toString + withSpaces.mkString("\n", "\n", "")
  }

}

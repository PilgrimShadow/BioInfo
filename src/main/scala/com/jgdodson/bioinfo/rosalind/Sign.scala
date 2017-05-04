package com.jgdodson.bioinfo.rosalind

object Sign {

  def main(args: Array[String]): Unit = {
    val n = args(0).toInt
    val count = numSignedPermutations(n)
    val perms = Perm.permutations2((1 to n).toSet).flatMap(signed).map(_.mkString(" "))
    println(count)
    println(perms.mkString("\n"))
  }

  def numSignedPermutations(n: Int): Int = {
    math.pow(2, n).toInt * (2 to n).product
  }

  def signed(perm: Vector[Int]): Set[Vector[Int]] = {

    def prependAll(n: Int, ss: Set[Vector[Int]]): Set[Vector[Int]] = {
      ss.map(item => n +: item)
    }

    if (perm.length == 1) Set(perm, perm.map(i => -i))
    else {
      prependAll(perm.head, signed(perm.tail)) ++ prependAll(-perm.head, signed(perm.tail))
    }
  }

}

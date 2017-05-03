package com.jgdodson.rosalind

import utils.Utils

object Pper {

  def main(args: Array[String]): Unit = {
    val line = Utils.readLines(args(0)).head.split(" ").map(_.toInt)
    val ans = partialPerm(line(0), line(1), 1000000)
    println(ans)
  }

  def partialPerm(n: Int, k: Int, mod: Int): Int = {
    ((n - k + 1) to n).foldLeft(1) { (acc, next) =>
      (acc * next) % mod
    }
  }

}

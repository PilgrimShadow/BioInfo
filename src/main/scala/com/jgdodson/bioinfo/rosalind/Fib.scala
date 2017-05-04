package com.jgdodson.bioinfo.rosalind

object Fib {

  def main(args: Array[String]): Unit = {
    val n = Integer.parseInt(args(0))
    val k = Integer.parseInt(args(1))
    println(genFib(k, n))
  }

  // genFib(k,n) = genFib(k, n-1) + k * genFib(k, n-2)
  def genFib(k: Long, n: Long): Long = {

    def loop(n: Long, n1: Long, n2: Long): Long = {
      if (n == 1) n2
      else loop(n - 1, n1 + k * n2, n1)
    }

    loop(n, 1, 1)
  }
}

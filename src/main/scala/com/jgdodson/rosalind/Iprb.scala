package com.jgdodson.rosalind

object Iprb {

  def main(args: Array[String]): Unit = {
    val k = Integer.parseInt(args(0))
    val m = Integer.parseInt(args(1))
    val n = Integer.parseInt(args(2))

    val cases = List(HomoDom, Hetero, HomoRec)

    val condProb = condMap(k, m, n)
    val domProb = (for (a <- cases; b <- cases) yield condProb(a, b) * parentProb(a, b)).sum
    println(domProb)
  }

  def parentProb(a: Genotype, b: Genotype): Double = (a, b) match {
    case (HomoDom, _) => 1.0
    case (_, HomoDom) => 1.0
    case (Hetero, Hetero) => 0.75
    case (Hetero, _) => 0.5
    case (_, Hetero) => 0.5
    case _ => 0.0
  }

  def condMap(k: Int, m: Int, n: Int): Map[(Genotype, Genotype), Double] = {
    val total: Double = k + m + n
    val norm = total * (total - 1)
    Map(
      (HomoDom, HomoDom) -> k * (k - 1) / norm,
      (HomoDom, Hetero) -> k * m / norm,
      (HomoDom, HomoRec) -> k * n / norm,
      (Hetero, HomoDom) -> m * k / norm,
      (Hetero, Hetero) -> m * (m - 1) / norm,
      (Hetero, HomoRec) -> m * n / norm,
      (HomoRec, HomoDom) -> n * k / norm,
      (HomoRec, Hetero) -> n * m / norm,
      (HomoRec, HomoRec) -> n * (n - 1) / norm
    )
  }

}

sealed abstract class Genotype

case object HomoDom extends Genotype

case object Hetero extends Genotype

case object HomoRec extends Genotype

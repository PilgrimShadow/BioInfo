package com.jgdodson.bioinfo.rosalind

object Lexf {

  def main(args: Array[String]): Unit = {

    val lines = io.Source.fromFile(args(0)).mkString.split("\n")
    val alpha = lines(0).split(" ").map(a => a(0))
    val order = makeOrdering(alpha)
    val len = Integer.parseInt(lines(1).trim())

    // Loosen up types to prevent these conversions (all we need is an iterable, right?)
    val words = enumerateFixed(len, alpha)
    val sorted = words.toVector.sorted(order)
    println(sorted.mkString("\n"))
  }

  def enumerateFixed[T](len: Int, alphabet: Seq[T]): Seq[String] = {

    assert(len >= 0)

    if (len == 0) Seq[String]()
    else if (len == 1) alphabet.map(_.toString)
    else alphabet.flatMap((next: T) => enumerateFixed(len - 1, alphabet).map(item => next.toString + item))
  }

  def makeOrdering(alphabet: Seq[Char]): Ordering[String] = new Ordering[String] {

    // Throw error if strings are not from the alphabet?

    def compare(xs: String, ys: String): Int = {

      // Could throw an error if strings are not upper/lower case instead of casefolding.

      val xLen = xs.length
      val yLen = ys.length
      val lim = math.min(xLen, yLen)

      for (i <- 0 until lim) {

        if (alphabet.indexOf(xs(i)) < alphabet.indexOf(ys(i))) return -1
        else if (alphabet.indexOf(xs(i)) > alphabet.indexOf(ys(i))) return 1
      }

      // If the intersection of the strings was equal, we compare the lengths.
      if (xLen < yLen) -1
      else if (yLen < xLen) 1
      else 0
    }
  }

}

package com.jgdodson.bioinfo.rosalind

// http://rosalind.info/problems/mrna/
object Mrna {

  def main(args: Array[String]): Unit = {

    val counts = makeCountMap(Prot.rnaCodonTable)

    val aminos = io.Source.fromFile(args(0)).mkString.replace("\n", "")

    val ans = aminos.foldLeft(1)((acc, next) => acc * counts(next.toString) % 1000000) * 3

    println(counts)

    println(ans)

  }

  def makeCountMap(table: Map[String, String]): Map[String, Int] = {

    def loop(table: Map[String, String], counts: Map[String, Int]): Map[String, Int] = {
      if (table.isEmpty) counts
      else {
        val key = table.head._2
        if (counts.isDefinedAt(key)) loop(table.tail, counts.updated(key, counts(key) + 1))
        else loop(table.tail, counts.updated(key, 1))
      }
    }

    loop(table, Map())
  }

}

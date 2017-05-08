package com.jgdodson.bioinfo.utils

import scalaj.http._

object Utils {

  // add a readFastas method for reading multiple entries and have readFasta
  // handle the singular case.

  def readFasta(raw: String): String = {
    val pat = "(?s)>([^\\n]*)\\n([^>]*)".r
    raw match {
      case pat(_, seq) => seq.trim
    }
  }

  def readFastas(raw: String): Vector[(String, String)] = {

    def processEntry(raw: String): (String, String) = {

      val lines = raw.split("\n")

      val name = lines(0).substring(1)

      val data = lines.slice(1, lines.length).mkString("")

      (name, data)
    }

    val pat = "(?s)>([^\\n]*)\\n([^>]*)".r

    val entries = pat.findAllIn(raw)

    entries.foldLeft(Vector[(String, String)]())((acc, next) => acc :+ processEntry(next))
  }

  def readFastaFile(fileName: String): Vector[(String, String)] = {
    val raw = io.Source.fromFile(fileName).mkString
    readFastas(raw)
  }

  def readLines(fileName: String): Vector[String] = {
    io.Source.fromFile(fileName).mkString.split("\n").toVector
  }

  def uniProtSequence(proteinID: String): String = {
    val response = Http("http://www.uniprot.org/uniprot/" + proteinID + ".fasta").asString.body
    readFasta(response)
  }


  def nearestMatch(target: Double, pool: Iterable[Double]): Double = {

    pool.minBy(candidate => math.abs(candidate - target))
  }

}

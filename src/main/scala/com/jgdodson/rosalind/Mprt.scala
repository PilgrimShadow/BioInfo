package com.jgdodson.rosalind

import scalaj.http.{Http, HttpOptions}
import utils.Utils


object Mprt {

  def main(args: Array[String]): Unit = {
    val lines = io.Source.fromFile(args(0)).mkString.split("\n")
    val aminos = lines.map(search)
    val indices = aminos.map(entry => (entry._1, findGlyMotif(entry._2)))
    val nonempty = indices.filter(entry => entry._2.length > 0)
    val formatted = nonempty.map(entry => entry._1 + "\n" + entry._2.mkString(" ")).mkString("\n")
    println(formatted)
  }


  def search(searchString: String): (String, String) = {

    val response = Http("http://www.uniprot.org/uniprot/" + searchString + ".fasta").
      option(HttpOptions.followRedirects(shouldFollow = true)).asString

    (searchString, Utils.readFastas(response.body).head._2)
  }

  def findGlyMotif(aminos: String): Vector[Int] = {
    val motif = "N[^P][ST][^P]".r
    Vector.range(0, aminos.length).filter(i => motif.findPrefixOf(aminos.substring(i)).isDefined).map(_ + 1)
  }
}

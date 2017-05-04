package com.jgdodson.bioinfo.rosalind

import com.jgdodson.bioinfo.utils.Utils
import com.jgdodson.bioinfo.DNAString

object Pdst {

  def main(args: Array[String]): Unit = {

    val group = Utils.readFastaFile(args(0)).map(fastaEntry => DNAString(fastaEntry._2))

    println(DNAString.normalizedDistanceMatrix(group).map(_.mkString(" ")).mkString("\n"))
  }

}

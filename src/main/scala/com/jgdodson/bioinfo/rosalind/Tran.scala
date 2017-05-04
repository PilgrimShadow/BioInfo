package com.jgdodson.bioinfo.rosalind

import com.jgdodson.bioinfo.DNAString
import com.jgdodson.bioinfo.utils.Utils.readFastaFile

object Tran {

  def main(args: Array[String]): Unit = {

    val seqs = readFastaFile(args(0)).map(_._2)

    val d1 = DNAString(seqs(0))
    val d2 = DNAString(seqs(1))

    println(d1.ttRatio(d2))
  }

}
